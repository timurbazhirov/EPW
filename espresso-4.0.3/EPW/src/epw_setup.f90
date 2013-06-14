  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_setup - Quantum-ESPRESSO group                
  !-----------------------------------------------------------------------
  SUBROUTINE epw_setup
  !-----------------------------------------------------------------------
  !
  ! 10/2009 looks clearly like phq_setup.  a lot of the 'fat'
  !         could be trimmed
  !
#include "f_defs.h"
  !
  USE kinds,         ONLY : DP
  USE ions_base,     ONLY : tau, nat, ntyp => nsp, ityp
  USE cell_base,     ONLY : at, bg  
  USE io_global,     ONLY : stdout, ionode, ionode_id
  USE ener,          ONLY : Ef
  USE klist,         ONLY : xk, lgauss, degauss, ngauss, nks, nelec, nkstot
  USE ktetra,        ONLY : ltetra, tetra
  USE lsda_mod,      ONLY : nspin, lsda, starting_magnetization
  USE scf,           ONLY : v, vrs, vltot, rho, rho_core, kedtau
  USE gvect,         ONLY : nrxx, ngm
  USE gsmooth,       ONLY : doublegrid
  USE symme,         ONLY : nsym, s, ftau, irt, t_rev, time_reversal
  USE uspp_param,    ONLY : upf
  USE spin_orb,      ONLY : domag
  USE constants,     ONLY : degspin, pi
  USE noncollin_module,     ONLY : noncolin, m_loc, angle1, angle2, ux
  USE wvfct,         ONLY : nbnd, et
  USE rap_point_group,      ONLY : code_group, nclass, nelem, elem, which_irr,&
                                  char_mat, name_rap, gname, name_class, ir_ram
  USE rap_point_group_is,   ONLY : code_group_is, gname_is
  use phcom
  USE control_flags, ONLY : iverbosity, modenum
  USE funct,         ONLY : dmxc, dmxc_spin, dmxc_nc, dft_is_gradient
  USE mp,            ONLY : mp_max, mp_min, mp_bcast
  USE mp_global,     ONLY : inter_pool_comm
  USE epwcom,        ONLY : xk_cryst
  implicit none

  real(DP) :: rhotot, rhoup, rhodw, target, small, fac, xmax, emin, emax
  ! total charge
  ! total up charge
  ! total down charge
  ! auxiliary variables used
  ! to set nbnd_occ in the metallic case
  ! minimum band energy
  ! maximum band energy

  real(DP) :: sr(3,3,48), sr_is(3,3,48)

  integer :: ir, table (48, 48), isym, jsym, irot, ik, jk, ibnd, ipol, &
       mu, nu, imode0, irr, ipert, na, it, is, js, nsym_is 
  ! counter on mesh points
  ! the multiplication table of the point g
  ! counter on symmetries
  ! counter on symmetries
  ! counter on rotations
  ! counter on k points
  ! counter on k points
  ! counter on bands
  ! counter on polarizations
  ! counter on modes
  ! the starting mode
  ! counter on representation and perturbat
  ! counter on atoms
  ! counter on iterations
  ! counter on atomic type

  real(DP) :: auxdmuxc(4,4)

  logical :: sym (48), is_symmorphic
  ! the symmetry operations

  CALL start_clock ('epw_setup')
  !
  ! 0) Set up list of kpoints in crystal coordinates
  !
  DO jk = 1, nkstot
     xk_cryst(:,jk) = xk(:,jk)
  END DO
  CALL cryst_to_cart (nkstot, xk_cryst, at, -1)
  CALL mp_bcast(xk_cryst,ionode_id)
  !
  ! 1) Computes the total local potential (external+scf) on the smooth grid
  !
  CALL set_vrs (vrs, vltot, v%of_r, kedtau, v%kin_r, nrxx, nspin, doublegrid)
  !
  ! 2) Set non linear core correction stuff
  !
  nlcc_any = ANY ( upf(1:ntyp)%nlcc )
  IF (nlcc_any) allocate (drc( ngm, ntyp))    
  !
  ! 3) Computes the derivative of the xc potential
  !
  dmuxc(:,:,:) = 0.d0
  IF (lsda) THEN
     DO ir = 1, nrxx
        rhoup = rho%of_r (ir, 1) + 0.5d0 * rho_core (ir)
        rhodw = rho%of_r (ir, 2) + 0.5d0 * rho_core (ir)
        CALL dmxc_spin (rhoup, rhodw, dmuxc(ir,1,1), dmuxc(ir,2,1), &
                                      dmuxc(ir,1,2), dmuxc(ir,2,2) )
     ENDDO
  ELSE
     IF (noncolin.and.domag) THEN
        DO ir = 1, nrxx
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           CALL dmxc_nc (rhotot, rho%of_r(ir,2), rho%of_r(ir,3), rho%of_r(ir,4), auxdmuxc)
           DO is=1,nspin
              DO js=1,nspin
                 dmuxc(ir,is,js)=auxdmuxc(is,js)
              END DO
           END DO
        ENDDO
     ELSE
        DO ir = 1, nrxx
           rhotot = rho%of_r (ir, 1) + rho_core (ir)
           IF (rhotot.gt.1.d-30) dmuxc (ir, 1, 1) = dmxc (rhotot)
           IF (rhotot.lt. - 1.d-30) dmuxc (ir, 1, 1) = - dmxc ( - rhotot)
        ENDDO
     END IF
  ENDIF
  !
  ! 3.1) Setup all gradient correction stuff
  !
  CALL setup_dgc
  !
  ! 4) Computes the inverse of each matrix
  !
  CALL multable (nsym, s, table)
  DO isym = 1, nsym
     DO jsym = 1, nsym
        IF (table (isym, jsym) == 1) invs (isym) = jsym
     ENDDO
  ENDDO
  !
  ! 5) Computes the number of occupied bands for each k point
  !
  IF (lgauss) THEN
     !
     ! discard conduction bands such that w0gauss(x,n) < small
     !
     ! hint:
     !   small = 1.0333492677046d-2  ! corresponds to 2 gaussian sigma
     !   small = 6.9626525973374d-5  ! corresponds to 3 gaussian sigma
     !   small = 6.3491173359333d-8  ! corresponds to 4 gaussian sigma
     !
     small = 6.9626525973374d-5
     !
     ! - appropriate limit for gaussian broadening (used for all ngauss)
     !
     xmax = sqrt ( - log (sqrt (pi) * small) )
     !
     ! - appropriate limit for Fermi-Dirac
     !
     IF (ngauss.eq. - 99) THEN
        fac = 1.d0 / sqrt (small)
        xmax = 2.d0 * log (0.5d0 * (fac + sqrt (fac * fac - 4.d0) ) )
     ENDIF
     target = ef + xmax * degauss
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           IF (et (ibnd, ik) .lt.target) nbnd_occ (ik) = ibnd
        ENDDO
        IF (nbnd_occ (ik) .eq.nbnd) WRITE( stdout, '(5x,/,&
             &"Possibly too few bands at point ", i4,3f10.5)') &
             ik,  (xk (ipol, ik) , ipol = 1, 3)
     ENDDO
  ELSEIF (ltetra) THEN
     CALL errore('epw_setup','phonon + tetrahedra not implemented', 1)
  ELSE
     IF (lsda) call infomsg('epw_setup','occupation numbers probably wrong')
     IF (noncolin) THEN
        nbnd_occ = nint (nelec) 
     ELSE
        DO ik = 1, nks
           nbnd_occ (ik) = nint (nelec) / degspin
        ENDDO
     ENDIF
  ENDIF
  !
  ! 6) Computes alpha_pv
  !
  emin = et (1, 1)
  DO ik = 1, nks
     DO ibnd = 1, nbnd
        emin = min (emin, et (ibnd, ik) )
     ENDDO
  ENDDO
#ifdef __PARA
  ! find the minimum across pools
  CALL mp_min( emin, inter_pool_comm )
#endif
  IF (lgauss) THEN
     emax = target
     alpha_pv = emax - emin
  ELSE
     emax = et (1, 1)
     DO ik = 1, nks
        DO ibnd = 1, nbnd_occ(ik)
           emax = max (emax, et (ibnd, ik) )
        ENDDO
     ENDDO
#ifdef __PARA
     ! find the maximum across pools
     CALL mp_max( emax, inter_pool_comm )
#endif
     alpha_pv = 2.d0 * (emax - emin)
  ENDIF
  ! avoid zero value for alpha_pv
  alpha_pv = max (alpha_pv, 1.0d-2)
  !
  ! 7) set all the variables needed to use the pattern representation
  !
  time_reversal = .NOT. (noncolin .AND. domag)
  ! 
  ! allocate and calculate rtau, the rotated position of each atom
  !
  DO isym = 1, nsym
     sym (isym) = .true.
  ENDDO

  CALL sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  nmodes = 3 * nat
  ! if minus_q=.t. set_irr will search for
  minus_q = (modenum .eq. 0)
  ! Sq=-q+G symmetry. On output minus_q=.t.
  ! if such a symmetry has been found
  IF (modenum .ne. 0) THEN
     CALL set_irr_mode (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
          nirr, gi, gimq, iverbosity, modenum)
  ELSE
     IF (nsym > 1.and..not.lgamma_gamma) THEN
        CALL set_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
             nirr, gi, gimq, iverbosity)
     ELSE
        CALL set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
             irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
             nirr, gi, gimq, iverbosity)
     ENDIF
  ENDIF
  is_symmorphic=.true.
  DO isym=1,nsymq
     is_symmorphic=( is_symmorphic.and.(ftau(1,irgq(isym))==0).and.  &
                                       (ftau(2,irgq(isym))==0).and.  &
                                       (ftau(3,irgq(isym))==0) )
  
  END DO
  search_sym=.true.
  IF (.not.is_symmorphic) THEN
     DO isym=1,nsymq
        search_sym=( search_sym.and.(abs(gi(1,irgq(isym)))<1.d-8).and.  &
                                    (abs(gi(2,irgq(isym)))<1.d-8).and.  &
                                    (abs(gi(3,irgq(isym)))<1.d-8) )
     END DO
  END IF
  IF (search_sym) THEN
     DO isym=1,nsym
        CALL s_axis_to_cart (s(1,1,isym), sr(1,1,isym), at, bg)
     END DO
     CALL find_group(nsym,sr,gname,code_group)
     CALL set_irr_rap(code_group,nclass,char_mat,name_rap,name_class,ir_ram)
     CALL divide_class(code_group,nsym,sr,nclass,nelem,elem,which_irr)
     IF (noncolin .and. domag) THEN
        nsym_is=0.d0
        DO isym=1,nsym
           IF (t_rev(isym)==0) THEN
              nsym_is=nsym_is+1
              CALL s_axis_to_cart (s(1,1,isym), sr_is(1,1,nsym_is), at, bg)
           ENDIF
        END DO
        CALL find_group(nsym_is,sr_is,gname_is,code_group_is)
     ENDIF
  ENDIF

  IF (lgamma_gamma) THEN
     ALLOCATE(has_equivalent(nat))
     ALLOCATE(n_equiv_atoms(nat))
     ALLOCATE(equiv_atoms(nat,nat))
     CALL find_equiv_sites (nat,nat,nsym,irt,has_equivalent,n_diff_sites, &
                       n_equiv_atoms,equiv_atoms)

     IF (n_diff_sites .LE. 0 .OR. n_diff_sites .GT. nat)            &
          &      CALL errore('epw_setup','problem with n_diff_sites',1)
     !
     ! look if ASR can be exploited to reduce the number of calculations
     ! we need to locate an independent atom with no equivalent atoms
     nasr=0
     IF (asr.AND.n_diff_sites.GT.1) THEN
        DO na = 1, n_diff_sites
           IF (n_equiv_atoms(na).EQ.1 ) THEN
              nasr = equiv_atoms(na, 1)
              GO TO 1
           END IF
        END DO
 1      CONTINUE
     END IF
  END IF


  IF (fildrho.ne.' ') call io_pattern (fildrho,nirr,npert,u,+1)

  !
  !  set maxirr if not already set
  !
  IF (maxirr.le.0.or.maxirr.gt.nirr) maxirr = nirr + 1
  IF (niter_ph.lt.maxter) maxirr = 1
  !
  !  set the alpha_mix parameter
  !
  DO it = 2, niter_ph
     IF (alpha_mix (it) .eq.0.d0) alpha_mix (it) = alpha_mix (it - 1)
  ENDDO
  !
  ! 8) Set the ubar
  !

  ubar(:) =( 0.d0,0.d0)
  !
  !   NB: the following instructions are for testing purposes of delta rho
  !       the user must know how many atoms there are in the system
  !
  !      ubar(1)=(1.d-3,0.d0)
  !      ubar(5)=(1.d0,0.d0)
  !      ubar(6)=(1.d0,0.d0)
  !
  !  9) set the variables needed for the partial computation
  !
     IF (nat_todo.eq.0) THEN
        !
        !    The partial computation option is not used, compute all atoms
        !
        DO na = 1, nat
           atomo (na) = na
        ENDDO
        nat_todo = nat
     ENDIF
     !
     !   Sets the atoms which must be computed: the requested atoms and all
     !   the symmetry related atoms
     !
     DO na = 1, nat
        IFat (na) = 0
     ENDDO
     DO na = 1, nat_todo
        IFat (atomo (na) ) = 1
        DO isym = 1, nsymq
           irot = irgq (isym)
           IFat (irt (irot, atomo (na) ) ) = 1
        ENDDO
     ENDDO
     !
     !    Computes again nat_todo, prepare the list atomo and sets all_comp
     !
     nat_todo = 0
     DO na = 1, nat
        IF (ifat (na) .eq.1) THEN
           nat_todo = nat_todo + 1
           atomo (nat_todo) = na
        ENDIF
     ENDDO
     !
     !     Find the irreducible representations to be computed
     !
     imode0 = 0
     DO irr = 1, nirr
        comp_irr (irr) = 0
        DO ipert = 1, npert (irr)
           mu = imode0 + ipert
           DO na = 1, nat
              IF (ifat (na) == 1 .and. comp_irr (irr) == 0) THEN
                 DO ipol = 1, 3
                    nu = 3 * (na - 1) + ipol
                    IF (abs (u (nu, mu) ) > 1.d-6)  comp_irr (irr) = 1
                 ENDDO
              ENDIF
           ENDDO
        ENDDO
        imode0 = imode0 + npert (irr)
     ENDDO
  !
  !   Initialize done_irr, find max dimension of the irreps
  !
  IF (lgamma_gamma) THEN
     comp_irr=0
     DO na=1,nat
        IF (has_equivalent(na)==0) THEN
            DO ipol=1,3
               comp_irr(3*(na-1)+ipol)=1
            ENDDO
        ENDIF
     ENDDO
     IF (nasr>0) THEN
        DO ipol=1,3
           comp_irr(3*(nasr-1)+ipol)=0
        ENDDO
     ENDIF     
  ENDIF
  all_comp = nat_todo.eq.nat
  npertx = 0
  DO irr = 1, nirr
     done_irr (irr) = 0
     npertx = max (npertx, npert (irr) )
  ENDDO

  CALL stop_clock ('epw_setup')
  !
  END SUBROUTINE epw_setup
