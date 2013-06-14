  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE elphon_shuffle_wrap
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation with Wannier functions: load all phonon q's
  !
  ! 09/2009 This SUBROUTINE is the main driver of the electron-phonon 
  ! calculation.  It first calculates the electron-phonon matrix elements
  ! on the coarse mesh and then passes the data off to ephwann_shuffle
  ! to perform the interpolation.
  !
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  !
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : my_pool_id, inter_pool_comm, root_pool, &
                            mpime, intra_pool_comm,npool
  USE mp,            ONLY : mp_barrier, mp_bcast
#endif
  USE us,            ONLY : nqxq, dq, qrad
  USE gvect,         ONLY : gcutm
  USE cellmd,        ONLY : cell_factor
  USE klist,         ONLY : xqq
  USE uspp_param,    ONLY : lmaxq, nbetam
  USE ions_base,     ONLY : nsp
  USE io_files,      ONLY : prefix, tmp_dir
  USE wavefunctions_module, ONLY: evc
  USE ions_base,     ONLY : nat, ntyp => nsp, tau, ityp
  USE control_flags, ONLY : iverbosity, modenum, nosym, noinv
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : nr1, nr2, nr3, et, xk, irt, s, nsym, ef,    &
                            ibrav, ftau, symm_type, igk, at, bg, sname, &
                            nks, nbnd, nkstot, ngm, nk1, nk2, nk3
  USE phcom,         ONLY : npert, u, gi, nirr, t, irotmq, max_irr_dim, &
                            tmq, dpsi, dvpsi, gimq, igkq, evq, nq1, nq3,&
                            nq2, nmodes,  minus_q, rtau, nsymq, irgq, &
                            lgamma, invs, xq
  USE symme,         ONLY : time_reversal
  USE epwcom,        ONLY : epbread, epbwrite, epwread,iuetf,phinterp, iunepmatf, &
                            phonselfen, elecselfen, nbndsub, elinterp,    &
                            iswitch, parallel_k, parallel_q, wannierize,  &
                            kmaps, nest_fn, fly, indabs, eig_read, fileig, &
                            band_plot, a2f
  USE el_phon,       ONLY : epmatq, dynq, et_all, xk_all, et_mb, et_ks
  implicit none
  real(kind=DP), allocatable :: xqc_irr(:,:), wqlist_irr(:), xqc(:,:), wqlist(:)
  ! the qpoints in the irr wedge
  ! the corresponding weigths
  ! the qpoints in the uniform mesh
  ! the corresponding weigths
  real(kind=DP), PARAMETER :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058,         &
       two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  integer :: nqc_irr, nqc, iuepb, max, nqxq_tmp, iueig, ibnd, ik, ios, &
             dummy1, dummy2, ik_start, ik_stop
  ! number of qpoints in the irreducible wedge
  ! number of qpoints on the uniform grid
  ! unit for .epb files
  ! 
  ! symmetry-related variables
  !
  integer :: gmapsym(ngm,48)
  !  correspondence G -> S(G)
  complex(kind=DP), PARAMETER :: czero = (0.d0,0.d0)
  complex(kind=DP) :: eigv (ngm, 48)
  ! e^{ iGv} for 1...nsym (v the fractional translation)
  complex(kind=DP) :: cz1( nmodes, nmodes), cz2(nmodes, nmodes)
  !  the eigenvectors for the first q in the star
  !  the rotated eigenvectors, for the current q in the star
  !
  integer :: nq, isq (48), imq 
  ! degeneracy of the star of q
  ! index of q in the star of a given sym.op.
  ! index of -q in the star of q (0 if not present)
  integer :: sym_sgq(48)
  ! the symmetries giving the q point iq in the star
  real(kind=DP) :: sxq (3, 48), et_tmp(nbnd, nkstot)
  ! list of vectors in the star of q
  integer :: i, j, iq, iq_irr, isym, &
     iq_first, jsym, ism1, nsq, ipol
  real(kind=DP) xq0(3), aq(3), saq(3), raq(3), sr(3,3), ft1, ft2, ft3
  logical :: sym(48),  eqvect_strict, nog, symmo, exst
  character (len=256) :: tempfile
  character (len=3) :: filelab
  !
  IF ( elinterp .and. (.not.phinterp ) ) CALL errore &
        ('elphon_shuffle_wrap','elinterp requires phinterp' ,1)
  !
  ! READ qpoint list from stdin
  !
#ifdef __PARA
  IF (mpime.eq.ionode_id) &
#endif
  READ(5,*) nqc_irr
#ifdef __PARA
  CALL mp_bcast (nqc_irr, ionode_id, inter_pool_comm)
  CALL mp_bcast (nqc_irr, root_pool, intra_pool_comm)
#endif
  allocate ( xqc_irr(3,nqc_irr), wqlist_irr(nqc_irr) )
  allocate ( xqc(3,nq1*nq2*nq3), wqlist(nq1*nq2*nq3) )
  !  
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
    DO iq = 1, nqc_irr
      READ (5,*) xqc_irr (:,iq), wqlist_irr (iq)
    ENDDO
#ifdef __PARA
  ENDIF
  CALL mp_bcast (xqc_irr, ionode_id, inter_pool_comm)
  CALL mp_bcast (xqc_irr, root_pool, intra_pool_comm)
#endif
  !
  ! fix for uspp
  max = nqxq
  DO iq = 1, nqc_irr
     nqxq_tmp = INT( ( (sqrt(gcutm) + sqrt(xqc_irr(1,iq)**2 + &
          xqc_irr(2,iq)**2 + xqc_irr(3,iq)**2) ) &
          / dq + 4) * cell_factor )
     IF (nqxq_tmp .gt. max)  max = nqxq_tmp
  ENDDO
  IF (max .gt. nqxq) then
     IF (allocated(qrad)) deallocate(qrad)
     allocate (qrad (max, nbetam*(nbetam+1)/2,lmaxq, nsp))
  ENDIF
  IF (nkstot .ne. nk1*nk2*nk3 ) &
       CALL errore('nscf run inconsistent with epw input',1)  
  !
  ! READ in external electronic eigenvalues. e.g. GW 
  !
  ALLOCATE(et_ks(nbnd,nks), et_mb(nbnd,nks))
  IF (eig_read) then
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
   WRITE (6,'(5x,a,i5,a,i5,a)') "Reading external electronic eigenvalues (", &
        nbnd, ",", nkstot,")"
   tempfile=trim(prefix)//'.eig'
   OPEN(1, file=tempfile, form='formatted', action='read', iostat=ios)
   IF (ios /= 0) CALL errore ('run_wannier','error OPENing' // tempfile, 1)
    DO ik = 1, nkstot
       DO ibnd = 1, nbnd
          READ (1,*) dummy1, dummy2, et_tmp (ibnd,ik)
          IF (dummy1.ne.ibnd) CALL errore('elphon_shuf_wrap', "Incorrect eigenvalue file", 1)
          IF (dummy2.ne.ik)   CALL errore('elphon_shuf_wrap', "Incorrect eigenvalue file", 1)
       ENDDO
    ENDDO
    CLOSE(1)
    et_tmp = et_tmp / ryd2ev
#ifdef __PARA
    ENDIF
    CALL mp_bcast (et_tmp, ionode_id, inter_pool_comm)
    CALL mp_bcast (et_tmp, root_pool, intra_pool_comm)
#endif
    !
    CALL ckbounds(ik_start, ik_stop)
    et_ks(:,:) = et(:,1:nks)
    et(:,1:nks) = et_tmp(:,ik_start:ik_stop)
    et_mb(:,:)  = et(:,1:nks)
 ENDIF
  !
  ! compute coarse grid dipole matrix elements.  Very fast 
  CALL compute_pmn_para
  !
  ! unit for the .epf and .etf file[s], temporary storage for interpolation 
  ! and for the .epb files (e-ph matrix in bloch representation)
  ! 
  iunepmatf = 76
  iuetf = 78
  iuepb = 85
  iueig = 86
  !
  !  gather electronic eigenvalues for subsequent shuffle
  !  
  allocate ( xk_all( 3, nkstot) , et_all( nbnd, nkstot) )
  CALL poolgather (    3, nkstot, nks, xk(:,1:nks),      xk_all)
  CALL poolgather ( nbnd, nkstot, nks, et(1:nbnd,1:nks), et_all)
  !
  !
  IF (.not.kmaps) then
     CALL start_clock('kmaps')
     CALL createkmap_pw2(xk_all,nkstot, xq0)
     CALL stop_clock('kmaps')
     CALL print_clock('kmaps')
  ELSE
     WRITE (stdout, '(/5x,a)')  'Using kmap and kgmap from disk'
  ENDIF
#ifdef __PARA
  CALL mp_barrier()
#endif
  !  if we start with a gamma point calculation, ../PW/set_kplusq.f90
  !  is not active and the gmap has not been produced...
  !
  IF (lgamma) CALL errore &
    ('elphon_shuffle_wrap','tshuffle2 requires q!=0 starting nscf calculation',1)
  !
  !  allocate dynamical matrix and ep matrix for all q's
  !
  allocate ( dynq (nmodes, nmodes, nq1*nq2*nq3), &
       epmatq (nbnd, nbnd, nks, nmodes, nq1*nq2*nq3))
  !
  !  NB: the symmetries are q-dependent and we should reset them
  !  for every q in the star. However, since we use the rotation of the
  !  wavefunctions (see notes), the dvscf_q_nu must be exactly the same
  !  as for the originating qpoints, and so the patterns etc.
  !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! ftau comes from the modules. Problem: when I call star_q to reset
! symmetries, in star_q.f90 ftau is a local variable. Therefore, in the
! call to sgama.f90 within star_q.f90 (which calls sgama_at.f90 where
! teh fractional translations are defined) the ftau are reset but this
! SUBROUTINE does not know that. And here I am getting the ftau corresponding
! to the weird q I USE in the nscf calculation. This means that we have to
! get the ftau's from star_q.
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !
  ! if we do not have epmatq already on file, we calculate the sandwiches here
  !
  IF (.not. epbread) THEN
    !
    nqc = 0
    iq_first = 1
    DO iq_irr = 1, nqc_irr
      !  
      WRITE(6,'(//5x,a)') repeat('=',67) 
      WRITE(6,'(5x,"irreducible q point # ",i4)') iq_irr
      WRITE(6,'(5x,a/)') repeat('=',67) 
      CALL flush(6)
      !
      xq = xqc_irr(:,iq_irr)
      !
      ! store the first q in the star
      !
      ! In the ph.x, xq0 would be the input q.  It is now the current q on the input list
      !
      xq0 = xq
      !
      !  reset symmetries of the crystal (not the small group of q)
      !  and determine the star of q
      !
      CALL star_q2 (xq, at, bg, ibrav, symm_type, nat, tau, ityp, nr1, nr2, nr3, &
        nsym, s, invs, irt, rtau, nq, sxq, isq, imq,  modenum, time_reversal, ftau)
      !
      !  determine the G vector map S(G) -> G 
      !  NB: after star_q we have the global symms, not those of the small group of q: 
      !  do not move this CALL.
      !
      IF (iq_irr.eq.1) then
        CALL gmap_sym ( nsym, s, ftau, gmapsym, eigv, invs)
        !  I checked that gmapsym(gmapsym(ig,isym),invs(isym)) = ig
      ENDIF
      !
      CALL readmat_shuffle2 ( iq_irr, nqc_irr, nq, iq_first, sxq, imq)
      ! now dynq is the cartesian dyn mat (NOT divided by the masses)
      !
      !  re-set the variables needed for the pattern representation
      !  and the symmetries of the small group of irr-q
      !  (from phq_setup.f90)
      !
      DO isym = 1, nsym
        sym (isym) = .true.
      ENDDO
      CALL sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)

      minus_q = (iswitch .gt. -3)  
      CALL set_irr       &
        (nat, at, bg, xq, s, invs, nsym, rtau, irt, irgq, nsymq, minus_q, &
         irotmq, t, tmq, max_irr_dim, u, npert, nirr, gi, gimq, iverbosity)
      !
      !
      !  if smallgq.f90 determined nsymq.eq.1, then the original phonon
      !  calculation has been done with the patterns set to the identity
      !  matrix, and we need to USE the same convention here 
      !
      IF ( nsymq .eq. 1 ) CALL set_irr_nosym  &
        (nat, at, bg, xq, s, invs, nsym, rtau, irt, irgq, nsymq, minus_q, &
         irotmq, t, tmq, max_irr_dim, u, npert, nirr, gi, gimq, iverbosity)
      !
      !  re-initialize the q-dependent quantities for the local
      !  and for the nonlocal part of the pseudopotentials
      !
      ! this initialization is done now on the q in the list 
      CALL epw_init(.false.)
      !
      !  loop over the q points of the star 
      !
      DO iq = 1, nq
        !
        ! retrieve the q in the star
        xq = sxq(:,iq)                               
        !
        ! and popolate the uniform grid
        nqc = nqc + 1
        xqc(:,nqc) = xq
        !
        IF (iq.eq.1) write(6,*)
        WRITE( stdout, 5) nqc, xq
        !
        !  prepare the gmap for the refolding
        !
        CALL createkmap ( xq )                      
        !
        IF (iverbosity.eq.1) then
          !
          !   description of symmetries
          !
          WRITE( stdout, '(36x,"s",24x,"frac. trans.")')
          DO isym = 1, nsym
            WRITE( stdout, '(/6x,"isym = ",i2,5x,a45/)') isym, sname(isym)
            CALL s_axis_to_cart (s(1,1,isym), sr, at, bg)
            IF (ftau(1,isym).ne.0.or.ftau(2,isym).ne.0.or.ftau(3,isym).ne.0) then
                ft1 = at(1,1)*ftau(1,isym)/nr1 + at(1,2)*ftau(2,isym)/nr2 + &
                      at(1,3)*ftau(3,isym)/nr3
                ft2 = at(2,1)*ftau(1,isym)/nr1 + at(2,2)*ftau(2,isym)/nr2 + &
                      at(2,3)*ftau(3,isym)/nr3
                ft3 = at(3,1)*ftau(1,isym)/nr1 + at(3,2)*ftau(2,isym)/nr2 + &
                      at(3,3)*ftau(3,isym)/nr3
                WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (s(1,ipol,isym),ipol=1,3), dble(ftau(1,isym))/dble(nr1)
                WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )")') &
                            (s(2,ipol,isym),ipol=1,3), dble(ftau(2,isym))/dble(nr2)
                WRITE( stdout, '(17x," (",3(i6,5x), " )       ( ",f10.7," )"/)') &
                            (s(3,ipol,isym),ipol=1,3), dble(ftau(3,isym))/dble(nr3)
                WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7, &
                      &        " )    f =( ",f10.7," )")') &
                      isym, (sr(1,ipol),ipol=1,3), ft1
                WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )")') &
                            (sr(2,ipol),ipol=1,3), ft2
                WRITE( stdout, '(17x," (",3f11.7, " )       ( ",f10.7," )"/)') &
                            (sr(3,ipol),ipol=1,3), ft3
             ELSE
                WRITE( stdout, '(1x,"cryst.",3x,"s(",i2,") = (",3(i6,5x), " )")') &
                                       isym,  (s (1, ipol, isym) , ipol = 1,3)
                WRITE( stdout, '(17x," (",3(i6,5x)," )")')  (s(2,ipol,isym), ipol=1,3)
                WRITE( stdout, '(17x," (",3(i6,5x)," )"/)') (s(3,ipol,isym), ipol=1,3)
                WRITE( stdout, '(1x,"cart. ",3x,"s(",i2,") = (",3f11.7," )")') &
                                              isym,  (sr (1, ipol) , ipol = 1, 3)
                WRITE( stdout, '(17x," (",3f11.7," )")')  (sr (2, ipol) , ipol = 1, 3)
                WRITE( stdout, '(17x," (",3f11.7," )"/)') (sr (3, ipol) , ipol = 1, 3)
             ENDIF
          ENDDO
          !
        ENDIF
        !
        ! isq(isym)=iq means: when we apply symmetry isym to the originating q 
        ! of the star, we get the iq-th member of the star. There are as many 
        ! matches as the degeneracy of the star.
        !
        ! pick up the first q in the small group of q*
        ! [it is important to select the first q in the small group, since
        ! it really corresponds to Sq. If we choose another element in the small group
        ! the actual q-point may be Sq+G and we screw up the q-vector below to generate
        ! k+q from k and for the KB projectors
        !
        nsq = 0 ! nsq is the degeneracy of the small group for this iq in the star
        !
        DO jsym = 1, nsym
          IF ( isq(jsym) .eq. iq ) then
             nsq = nsq + 1
             sym_sgq(nsq) = jsym
          ENDIF
        ENDDO
        IF ( nsq*nq .ne. nsym ) CALL errore ('elphon_shuffle_wrap', 'wrong degeneracy', iq)
        ! 
        IF (iverbosity.eq.1) then
          !
          WRITE(6,*) 'iq, i, isym, nog, symmo'
          DO i = 1, nsq
            !
            isym = sym_sgq(i)
            ism1 = invs (isym)
            !
            !  check for G such that Sq = q* + G 
            ! 
            aq  = xq0
            saq = xq
            CALL cryst_to_cart (1, aq, at, -1)
            DO j = 1, 3
              raq (j) = s (j, 1, ism1) * aq (1) &
                      + s (j, 2, ism1) * aq (2) &
                      + s (j, 3, ism1) * aq (3)
            ENDDO
            CALL cryst_to_cart (1, saq, at, -1)
            nog = eqvect_strict( raq, saq) 
            !
            !  check whether the symmetry belongs to a symmorphic group
            !
            symmo = (ftau(1,isym).eq.0 .and. ftau(2,isym).eq.0 .and. ftau(3,isym).eq.0)
            !
            WRITE(6,'(3i5,2a)') iq, i, isym, nog, symmo
            !
          ENDDO  
          !
        ENDIF
        ! 
        isym = sym_sgq (1)  ! primo q del piccolo gruppo, assicura che G=0
        !
        !   calculate the sandwiches
        !
        ! a more accurate way of doing is to symmetrize the matrix element w.r.t.
        ! the small group of the given q in the star. I'm not doint this here.
        ! (but I checked that even without symm the result of full zone and irr zone
        ! are equal to 5+ digits)
        ! For any volunteers, please write to giustino@civet.berkeley.edu
        !
        CALL elphon_shuffle ( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, invs, xq0, .false. )
        !
        !  bring epmatq in the mode representation of iq_first, 
        !  and then in the cartesian representation of iq
        !
        CALL rotate_eigenm ( iq_first, iq, nqc, isym, nsym, s, invs, irt, &
           rtau, xq, isq, cz1, cz2, .false. )
        !
        CALL rotate_epmat ( cz1, cz2, xq, nqc )
        !
! To be removed when the imq.eq.0 section is tested and corrected
goto 123
        IF (imq .eq. 0) then
          ! NOT CURRENTLY FUNCTIONAL 
          CALL errore('elphon_shuffle_wrap', 'imq.eq.0 not correct!',1)
          ! 
          ! -q is not in the star: we introduce it here using time-reversal symmetry
          ! [Eq. (4.4) Maradudin&Vsoko RMP]
          !
          ! retrieve the -q in the star
          nqc = nqc + 1
          xq = - xq
          xqc(:,nqc) = xq  
          !
          WRITE( stdout, 5) nqc, xq
          !
          !
          !  prepare the gmap for the refolding
          !
          CALL createkmap ( xq )
          !
          CALL elphon_shuffle ( iq_irr, nqc_irr, nqc, gmapsym, eigv, isym, invs, xq0, .true. )
          !
          !  bring epmatq in the mode representation of iq_first,
          !  and then in the cartesian representation of iq
          !  taking into account time reversal
          !
          CALL rotate_eigenm ( iq_first, iq, nqc, isym, nsym, s, invs, irt, &
             rtau, xqc(:,nqc-1), isq, cz1, cz2, .true. )
          !
          !
          CALL rotate_epmat ( cz1, cz2, xq, nqc )
          !
        ENDIF
! To be removed when the imq.eq.0 section is tested and corrected
123 continue
        !
        CALL flush(stdout) ! opteron ONLY
        !
      ENDDO  
      !
      iq_first = iq_first + nq
! To be added back in when the imq.eq.0 section is tested and corrected
!      if (imq .eq. 0) iq_first = iq_first + nq
      !
    ENDDO
    !
  ENDIF
  !
  IF ( epbread .or. epbwrite ) then
    !
    ! write the e-ph matrix elements and other info in the Bloch representation
    ! in .epb files (one for each pool)
    !
    tempfile = trim(tmp_dir) // trim(prefix) // '.epb' 
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
     tempfile = trim(tmp_dir) // trim(prefix) // '.epb' // filelab
#endif
     !
     IF (epbread)  THEN
        inquire(file = tempfile, exist=exst)
        IF (.not. exst ) CALL errore( 'elphon_shuffle_wrap', 'epb files not found ', 1)
        OPEN  (iuepb, file = tempfile, form = 'unformatted')
        WRITE(6,'(/5x,"Reading epmatq from .epb files"/)') 
        READ  (iuepb) nqc, xqc, et, dynq, epmatq
        CLOSE (iuepb)
     ENDIF
     !
     IF (epbwrite) THEN
        OPEN  (iuepb, file = tempfile, form = 'unformatted')
        WRITE(6,'(/5x,"Writing epmatq on .epb files"/)') 
        WRITE (iuepb) nqc, xqc, et, dynq, epmatq
        CLOSE (iuepb)
     ENDIF
  ENDIF
  !
  !
#ifdef __PARA
  CALL mp_barrier()
#endif
  !
  IF (nqc.ne.nq1*nq2*nq3) &
    CALL errore ('elphon_shuffle_wrap','nqc .ne. nq1*nq2*nq3',nqc)
  wqlist = float(1)/float(nqc)
  !
  !   now dynq is the cartesian dyn mat ( NOT divided by the masses)
  !   and epmatq is the epmat in cartesian representation (rotation in elphon_shuffle)
  !
  !
  ! free up some memory
  !
  IF ( ASSOCIATED (evq)  ) NULLIFY    (evq)
  IF ( ALLOCATED  (evc)  ) DEALLOCATE (evc)
  IF ( ASSOCIATED (igkq) ) NULLIFY    (igkq)
  IF ( ALLOCATED  (igk)  ) DEALLOCATE (igk)
  IF ( ALLOCATED  (dvpsi)) DEALLOCATE (dvpsi)
  IF ( ALLOCATED  (dpsi) ) DEALLOCATE (dpsi)
  !
  ! the electron-phonon wannier interpolation
  !
  CALL ephwann_shuffle ( nqc, xqc )
  !
  ! the calculation of the self energy
  ! 
  !
  IF (elecselfen  .and. .not.fly) CALL selfen_elec
  IF (phonselfen  .and. .not.fly) CALL selfen_phon
  IF (nest_fn     .and. .not.fly) CALL nesting_fn
  IF (band_plot   .and. .not.fly) CALL plot_band
  !
5 format (8x,"q(",i5," ) = (",3f12.7," )") 
  !
  RETURN
  END SUBROUTINE elphon_shuffle_wrap
  !
  !---------------------------------------------------------------------------
  SUBROUTINE irotate ( x, s, sx)
  !---------------------------------------------------------------------------
  !
  ! a simple symmetry operation in crystal coordinates ( s is integer!)
  !
  USE kinds, ONLY : DP
  implicit none
  real(kind=DP) x(3), sx(3)
  integer :: s(3,3),i
  !
  DO i = 1, 3
     sx (i) = float(s (i, 1)) * x (1) &
            + float(s (i, 2)) * x (2) &
            + float(s (i, 3)) * x (3)
  ENDDO
  !
  RETURN
  end
  !---------------------------------------------------------------------------
  logical function eqvect_strict (x, y)
  !-----------------------------------------------------------------------
  !
  !   This function test if two tridimensional vectors are equal
  !
  USE kinds
  implicit none
  real(kind=DP) :: x (3), y (3)
  ! input: input vector
  ! input: second input vector
  real(kind=DP) :: accep
  ! acceptance PARAMETER
  PARAMETER (accep = 1.0d-5)
  !
  eqvect_strict = abs( x(1)-y(1) ).lt.accep .and. &
                  abs( x(2)-y(2) ).lt.accep .and. &
                  abs( x(3)-y(3) ).lt.accep
  END FUNCTION eqvect_strict
