  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !----------------------------------------------------------------------
  SUBROUTINE ephwann_shuffle( nqc, xqc )
  !---------------------------------------------------------------------
  !
  !  Wannier interpolation of electron-phonon vertex
  !
  !  Scalar implementation   Feb 2006
  !  Parallel version        May 2006
  !  Disentenglement         Oct 2006
  !  Compact formalism       Dec 2006
  !  Phonon irreducible zone Mar 2007
  !
  !  No ultrasoft now
  !  No spin polarization
  !
  !-----------------------------------------------------------------------
  !
  !
  !
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : nbnd, nks, nkstot, nk1, nk2, nk3, &
                            et, xk, at, bg, ef,  celldm, nelec,amconv
  USE ions_base,     ONLY : amass 
  USE phcom,         ONLY : dyn, lgamma, nq1, nq2, nq3, nmodes, w2
  USE epwcom,        ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, & 
                            iuetf, lretf, epwread, epwwrite,       &
                            iunepmatwp, iunepmatwe, ngaussw, degaussw, &
                            nbndskip, parallel_k, parallel_q, epf_mem, etf_mem, &
                            phonselfen, nest_fn, fly, a2f, indabs, &
                            epexst, wepexst, vme, eig_read, emaxabs, &
                            eminabs, nsmear, delta_smear, twophoton, ephmatwrite
  USE control_flags, ONLY : iverbosity
  USE io_files,      ONLY : prefix, tmp_dir
  USE ions_base,     ONLY : ityp
  USE io_global,     ONLY : stdout
  USE el_phon,       ONLY : nrr_k, nrr_q, cu, cuq, irvec, ndegen_k, ndegen_q,   &
                            wslen, chw, chw_ks, cvmew, cdmew, rdw, epmatwp, epmatq, &
                            wf, etf, etfq, etf_ks, xqf, xkf, wkf, wqf, epmatw17, &
                            dynq, nxqf, nksf, epf17, nksqf, et_ks, &
                            ibndmin, ibndmax, lambda_all, lambda_v_all, dmec, dmef, vmef
#ifdef __PARA
  USE mp,            ONLY : mp_barrier, mp_bcast, mp_sum
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : mpime, nproc, my_pool_id, nproc_pool, intra_image_comm, &
                            inter_pool_comm, me_pool, root_pool, mpime, intra_pool_comm, &
                            my_pool_id
#endif
  !
  implicit none
  !
  complex(kind=DP), ALLOCATABLE :: &
    epmatwe  (:,:,:,:),         &! e-p matrix  in wannier basis - electrons
    epmatwef (:,:,:,:)           ! e-p matrix  in el wannier - fine Bloch phonon grid
  complex(kind=DP), ALLOCATABLE :: &
    epmatf( :, :, :),           &! e-p matrix  in smooth Bloch basis, fine mesh
    cufkk ( :, :),              &! Rotation matrix, fine mesh, points k
    cufkq ( :, :),              &! the same, for points k+q
    uf    ( :, :)                ! Rotation matrix for phonons
  integer :: &
    nqc,                        &! number of qpoints in the coarse grid
    iunepwdata                   ! number of qpoints in the coarse grid
  real(kind=DP) :: &
    xqc (3, nqc)                 ! qpoint list, coarse mesh
  !
  integer :: iq, ik, ikk, ikq, ibnd, jbnd, imode, ir, na, nu, mu, &
    fermicount, nrec, indnew, indold, lrepmatw, ios, unf_recl, &
    statb(13)
  integer, ALLOCATABLE ::  do_q(:)
  integer :: fstat
  logical :: exst, opnd
  character (len=256) :: filint, filename, tempfile
  real(kind=DP) :: xxq(3), xxk(3), xkk(3), xkq(3), size_m
  real(kind=DP), PARAMETER :: rydcm1 = 13.60580 * 8065.5d0, &
     ryd2mev = 13605.8, ryd2ev = 13.6058, two = 2.d0
  !
  real(kind=DP), PARAMETER :: bohr2ang = 0.5291772108
  ! 
  iunepwdata = 79
  !
  IF (nbndsub.ne.nbnd) &
       WRITE(stdout, '(/,14x,a,i4)' ) 'band disentanglement is used:  nbndsub = ', nbndsub
  !
  ALLOCATE ( cu ( nbnd, nbndsub, nks), & 
             cuq ( nbnd, nbndsub, nks), & 
             irvec (3, 20*nk1*nk2*nk3), &
             ndegen_k (20*nk1*nk2*nk3), &
             ndegen_q (20*nq1*nq2*nq3), &
             wslen(20*nk1*nk2*nk3)      )
  !
  CALL start_clock ( 'ephwann' )
 !
  IF (parallel_k) THEN
     CALL loadqmesh_serial
     CALL loadkmesh_para
  ELSE
     CALL errore('ephwann_shuffle', "parallel q not (yet) implemented",1)
  ENDIF
  !
  !
  ! determine Wigner-Seitz points
  !
  CALL wigner_seitz2 &
       ( nk1, nk2, nk3, nq1, nq2, nq3, nrr_k, nrr_q, irvec, wslen, ndegen_k, ndegen_q )
  !
  !
  ! allocate dipole matrix elements after getting grid size
  !
  ALLOCATE ( dmef(3, nbndsub, nbndsub, 2 * nksqf) )
  IF (vme) ALLOCATE ( vmef(3, nbndsub, nbndsub, 2 * nksqf) )
  !
  ! open the .epmatwe and .epmatwp file[s] with the proper record length
  !
  iunepmatwp = 111
  iunepmatwe = 112
  lrepmatw   = 2 * nbndsub * nbndsub * nrr_k * nmodes 
  !
  filint    = trim(prefix)//'.epmatwe'
  IF (.not.epwread) CALL diropn (iunepmatwe, 'epmatwe', lrepmatw, exst)
  !
  filint    = trim(prefix)//'.epmatwp'
  CALL diropn (iunepmatwp, 'epmatwp', lrepmatw, exst)
  !
  ! at this point, we will interpolate the Wannier rep to the Bloch rep 
  !
  IF ( epwread ) THEN
     !
     !  read all quantities in Wannier representation from file
     !  in parallel case all pools read the same file
     !
     CALL epw_read(iunepwdata)
     !
  ELSE !if not epwread (i.e. need to calculate fmt file)
     !
     xxq = 0.d0 
     CALL loadumat &
          ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq )  
     !
     ! ------------------------------------------------------
     !   Bloch to Wannier transform
     ! ------------------------------------------------------
     !
     ALLOCATE ( epmatwe ( nbndsub, nbndsub, nrr_k, nmodes), &
          epmatwp ( nbndsub, nbndsub, nrr_k, nmodes), &
          chw     ( nbndsub, nbndsub, nrr_k ),        &
          chw_ks     ( nbndsub, nbndsub, nrr_k ),        &
          cdmew   ( 3, nbndsub, nbndsub, nrr_k ),        &
          rdw     ( nmodes,  nmodes,  nrr_q ) )
     IF (vme) ALLOCATE(cvmew   ( 3, nbndsub, nbndsub, nrr_k ) )
     !
     ! Hamiltonian
     !
     CALL hambloch2wan &
          ( nbnd, nbndsub, nks, nkstot, lgamma, et,    xk, cu, nrr_k, irvec, wslen, chw )
     !
     ! Kohn-Sham eigenvalues
     !
     IF (eig_read) THEN
        WRITE (6,'(5x,a)') "Interpolating MB and KS eigenvalues"
        CALL hambloch2wan &
             ( nbnd, nbndsub, nks, nkstot, lgamma, et_ks, xk, cu, nrr_k, irvec, wslen, chw_ks )
     ENDIF
     !
     ! Dipole
     !
     CALL dmebloch2wan &
          ( nbnd, nbndsub, nks, nkstot, nkstot, lgamma, dmec, xk, cu, nrr_k, irvec, wslen )
     !
     ! Dynamical Matrix 
     !
     CALL dynbloch2wan &
          ( nmodes, nqc, xqc, dynq, nrr_q, irvec, wslen )
     !
     ! Transform of position matrix elements
     ! PRB 74 195118  (2006)
     IF (vme) CALL vmebloch2wan &
         ( nbnd, nbndsub, nks, nks, nkstot, lgamma, xk, cu, nrr_k, irvec, wslen )
     !
     ! Electron-Phonon vertex (Bloch el and Bloch ph -> Wannier el and Bloch ph)
     !
     DO iq = 1, nqc
        !
        xxq = xqc (:, iq)
        !
        ! we need the cu again for the k+q points, we generate the map here
        !
        CALL loadumat ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq )
        !
        DO imode = 1, nmodes
           !
           CALL ephbloch2wane &
                ( nbnd, nbndsub, nks, nkstot, lgamma, xk, cu, cuq, &
                epmatq (:,:,:,imode,iq), nrr_k, irvec, wslen, epmatwe(:,:,:,imode) )
           !
        ENDDO
        !
        ! direct write of epmatwe for this iq 
        IF (.not.epexst) CALL rwepmatw ( epmatwe, nbndsub, nrr_k, nmodes, iq, iunepmatwe, +1)
        !
     ENDDO
     !
     ! Electron-Phonon vertex (Wannier el and Bloch ph -> Wannier el and Wannier ph)
     !
     CALL ephbloch2wanp &
          ( nbndsub, nmodes, xqc, nqc, dynq, irvec, wslen, nrr_k, nrr_q, epmatwe )
     !
     IF ( epwwrite ) THEN
        CALL epw_write(iunepwdata) 
        CALL epw_read (iunepwdata)
     ENDIF
     !
  ENDIF
  !
  !
  IF ( ALLOCATED (epmatwe) ) DEALLOCATE (epmatwe)
  IF ( ALLOCATED (epmatq) )  DEALLOCATE (epmatq)
  IF ( ALLOCATED (cu) )      DEALLOCATE (cu)
  IF ( ALLOCATED (cuq) )     DEALLOCATE (cuq)
  !
  ! at this point, we will interpolate the Wannier rep to the Bloch rep 
  ! for electrons, phonons and the ep-matrix
  !
  !  need to add some sort of parallelization (on g-vectors?)  what
  !  else can be done when we don't ever see the wfcs??
  !
  ALLOCATE ( epmatwef( nbndsub, nbndsub, nrr_k, nmodes),             &
       wf ( nmodes,  nxqf ), etf ( nbndsub, nksf),                   &
       etf_ks ( nbndsub, nksf),                   &
       epmatf( nbndsub, nbndsub, nmodes), cufkk ( nbndsub, nbndsub), &
       cufkq ( nbndsub, nbndsub), uf ( nmodes, nmodes))
  !
  IF (fly) THEN
     ALLOCATE ( etfq(nbndsub, nksf,    1) )
  ELSE
     ALLOCATE ( etfq(nbndsub, nksf, nxqf) )
  ENDIF
  !
  ! this loops over the fine mesh of q points.
  ! if parallel_k then this is the entire q-list (nxqftot)
  ! if parallel_q then this is nxqftot/npool
  !
  DO iq = 1, nxqf
     !
     xxq = xqf (:, iq)
     !
     ! ------------------------------------------------------
     ! dynamical matrix : Wannier -> Bloch
     ! ------------------------------------------------------
     !
     CALL dynwan2bloch &
          ( nmodes, nrr_q, irvec, ndegen_q, xxq, uf, w2 )
     !
     DO nu = 1, nmodes
        IF ( w2 (nu) .gt. 0.d0 ) THEN
           wf(nu,iq) =  sqrt(abs( w2 (nu) ))
        ELSE
           wf(nu,iq) = -sqrt(abs( w2 (nu) ))
        ENDIF
     ENDDO
     !
  ENDDO
  !
  ! ------------------------------------------------------
  ! hamiltonian : Wannier -> Bloch (preliminary)
  ! ------------------------------------------------------
  !
  ! we here perform a preliminary interpolation of the hamiltonian
  ! in order to determine the fermi window ibndmin:ibndmax for later use.
  ! We will interpolate again afterwards, for each k and k+q separately
  !
  xxq = 0.d0
  !
  ! nksf is the number of kpoints in the pool
  ! parallel_k case = nkstotf/npool
  ! parallel_q case = nkstotf
  DO ik = 1, nksf
     !
     xxk = xkf (:, ik)
     !
     IF ( 2*(ik/2).eq.ik ) THEN
        !
        !  this is a k+q point : redefine as xkf (:, ik-1) + xxq
        !
        CALL cryst_to_cart ( 1, xxq, at,-1 )
        xxk = xkf (:, ik-1) + xxq
        CALL cryst_to_cart ( 1, xxq, bg, 1 )
        !
     ENDIF
     !
     CALL hamwan2bloch &
          ( nbndsub, nrr_k, irvec, ndegen_k, xxk, cufkk, etf (:, ik), chw)
     !
     !
  ENDDO
  !
  ! identify the bands within fsthick from the Fermi level
  ! (in shuffle mode this actually does not depend on q)
  !
  CALL fermiwindow
  !
  ! get the size of the matrix elements stored in each pool
  ! for informational purposes.  Not necessary
  !
  CALL mem_size(ibndmin, ibndmax, nmodes, nksqf, nxqf, fly) 
  !
  IF (etf_mem) THEN
     ! Fine mesh set of g-matrices.  It is large for memory storage
     IF (fly) THEN
        ALLOCATE ( epf17 (nksqf, 1,    ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes) )
     ELSE
        ALLOCATE ( epf17 (nksqf, nxqf, ibndmax-ibndmin+1, ibndmax-ibndmin+1, nmodes) )
     ENDIF
     !
  ELSE
     !
     !  open epf and etf files with the correct record length
     !
     lrepmatf  = 2 * (ibndmax-ibndmin+1) * (ibndmax-ibndmin+1)
     CALL diropn (iunepmatf, 'epf', lrepmatf, exst)
     !
     lretf     = (ibndmax-ibndmin+1) 
     CALL diropn (iuetf, 'etf', lretf, exst)
     !
  ENDIF
  !
  !
  IF (epf_mem) THEN
     !
     ALLOCATE ( epmatw17 ( nbndsub, nbndsub, nrr_k, nrr_q, nmodes) )
     !  
     !  direct read of epmatw17 - wannier matrices on disk from epwwrite
     ios = fstat ( iunepmatwp, statb)
     size_m = float(statb(8))/(1024.d0**2)
     WRITE (6,'(5x,a,f10.2,a)') "Loading Wannier rep into memory: ", size_m, " MB"
     DO ir = 1, nrr_q
        CALL rwepmatw ( epmatw17(:,:,:,ir,:), nbndsub, nrr_k, nmodes, ir, iunepmatwp, -1)
     ENDDO
     !
  ENDIF
  !
  !  xqf must be in crystal coordinates
  !
  ALLOCATE(do_q(nxqf))
  CALL collect_q_list(do_q, nxqf)
  !
  DO iq = 1, nxqf
     !   
     IF (do_q(iq) .eq. 0) CYCLE
     !
     CALL start_clock ( 'ep-interp' )
     !
     IF (.not.fly) THEN
        IF (iverbosity .ge. 1) THEN
           WRITE(6,'(/5x,"Interpolating ",i6," out of ",i6)') iq, nxqf
           CALL flush(6)
        ELSE
           !
           IF (iq.eq.1) THEN
              WRITE(6,'(/5x,"Interpolation progress: ")',advance='no')
              indold = 0
           ENDIF
           indnew = nint(float(iq)/float(nxqf)*40)
           IF (indnew.ne.indold) WRITE(6,'(a)',advance='no') '#'
           indold = indnew
           CALL flush(6) ! ONLY on opteron
       ENDIF
     ENDIF
     !
     xxq = xqf (:, iq)
     !
     ! ------------------------------------------------------
     ! dynamical matrix : Wannier -> Bloch
     ! ------------------------------------------------------
     !
     CALL dynwan2bloch &
          ( nmodes, nrr_q, irvec, ndegen_q, xxq, uf, w2 )
     !
     ! ...then take into account the mass factors and square-root the frequencies...
     !
     DO nu = 1, nmodes
        !
        ! wf are the interpolated eigenfrequencies
        ! (omega on fine grid)
        !
        IF ( w2 (nu) .gt. 0.d0 ) THEN
           wf(nu,iq) =  sqrt(abs( w2 (nu) ))
        ELSE
           wf(nu,iq) = -sqrt(abs( w2 (nu) ))
        ENDIF
        !
        DO mu = 1, nmodes
           na = (mu - 1) / 3 + 1
           uf (mu, nu) = uf (mu, nu) / sqrt(amass(ityp(na)))
        ENDDO
     ENDDO
     !
     ! --------------------------------------------------------------
     ! epmat : Wannier el and Wannier ph -> Wannier el and Bloch ph
     ! --------------------------------------------------------------
     !
     CALL ephwan2blochp &
          ( nmodes, xxq, irvec, ndegen_q, nrr_q, uf, epmatwef, nbndsub, nrr_k )
     !
     !
     !  number of k points with a band on the Fermi surface
     fermicount = 0
     !
     ! this is a loop over k blocks in the pool
     ! (size of the local k-set)
     DO ik = 1, nksqf
        !
        ! xkf is assumed to be in crys coord
        !
        ikk = 2 * ik - 1
        ikq = ikk + 1
        !
        xkk = xkf(:, ikk)
        xkq = xkk + xxq
        !
        ! ------------------------------------------------------        
        ! hamiltonian : Wannier -> Bloch 
        ! ------------------------------------------------------
        !
        ! Kohn-Sham first, then get the rotation matricies for following interp.
        IF (eig_read) THEN
           CALL hamwan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf_ks (:, ikk), chw_ks)
           CALL hamwan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf_ks (:, ikq), chw_ks)
        ENDIF
        !
        CALL hamwan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etf (:, ikk), chw)
        CALL hamwan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, etf (:, ikq), chw)
        !
        ! ------------------------------------------------------        
        !  dipole: Wannier -> Bloch
        ! ------------------------------------------------------        
        !
        CALL dmewan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, dmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk))
        CALL dmewan2bloch &
             ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, dmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq))
        !
        ! ------------------------------------------------------        
        !  velocity: Wannier -> Bloch
        ! ------------------------------------------------------        
        !
        IF (vme) THEN
           IF (eig_read) THEN
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw_ks)
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw_ks)
           ELSE
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, vmef (:,:,:, ikk), etf(:,ikk), etf_ks(:,ikk), chw)
              CALL vmewan2bloch &
                   ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkq, vmef (:,:,:, ikq), etf(:,ikq), etf_ks(:,ikq), chw)
           ENDIF
        ENDIF
        !
        !
        IF (etf_mem) THEN
           ! store in mem, otherwise the parall in q is going to read a mess on file
           ! this is an array size (ibndmax-ibndmin+1)*(k blocks in the pool)*(total number of qs on fine mesh)
           ! parallel in K!
           IF (fly) THEN
              etfq(:,ikk, 1) = etf (:, ikk)
              etfq(:,ikq, 1) = etf (:, ikq)
           ELSE
              etfq(:,ikk,iq) = etf (:, ikk)
              etfq(:,ikq,iq) = etf (:, ikq)
           ENDIF
           !
        ELSE
           nrec = (iq-1) * nksf + ikk
           IF (fly) nrec  = ikk
           CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, + 1)
           nrec = (iq-1) * nksf + ikq
           IF (fly) nrec  = ikq
           CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, + 1)
           !
        ENDIF
        !
        ! interpolate ONLY when (k,k+q) both have at least one band 
        ! within a Fermi shell of size fsthick 
        !
        IF ( (( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
             ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick )) .or. indabs ) THEN
           !
           !  fermicount = fermicount + 1
           !
           ! --------------------------------------------------------------
           ! epmat : Wannier el and Bloch ph -> Bloch el and Bloch ph
           ! --------------------------------------------------------------
           !
           !
           CALL ephwan2bloch &
                ( nbndsub, nrr_k, irvec, ndegen_k, epmatwef, xkk, cufkk, cufkq, epmatf, nmodes )
           !
           ! 
           ! write epmatf to file / store in memory
           !
           !
           DO imode = 1, nmodes
              !
              IF (etf_mem) THEN
                 !
                 DO jbnd = ibndmin, ibndmax
                    DO ibnd = ibndmin, ibndmax
                       !
                       IF (fly) THEN
                          epf17(ik, 1,jbnd-ibndmin+1,ibnd-ibndmin+1,imode) = epmatf(jbnd,ibnd,imode)
                       ELSE
                          epf17(ik,iq,jbnd-ibndmin+1,ibnd-ibndmin+1,imode) = epmatf(jbnd,ibnd,imode)
                       ENDIF
                       !
                    ENDDO
                 ENDDO
                 !
              ELSE
                 !
                 nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
                 IF (fly)  nrec = (imode-1) * nksqf + ik
                 CALL dasmio ( epmatf(ibndmin:ibndmax,ibndmin:ibndmax,imode), &
                      ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, +1)
                 !
              ENDIF
              !
           ENDDO
           !
        ENDIF
        !
     ENDDO  ! end loop over k points
     !
     IF (phonselfen  .and. fly) CALL selfen_phon_fly( iq )
     IF (nest_fn     .and. fly) CALL nesting_fn_fly( iq )
!     IF (indabs      .and. fly) CALL indabs_fly (iq)
!     IF (twophoton   .and. fly) CALL twophoton_fly (iq)
     IF (ephmatwrite .and. fly) CALL write_ephmat_fly( iq ) 
     !
     CALL stop_clock ( 'ep-interp' )
     !
  ENDDO  ! end loop over q points
  !
  IF (a2f .and. fly) CALL eliashberg_a2f( lambda_all(:,:,1), lambda_v_all(:,:,1))
  !
  CALL stop_clock ( 'ephwann' )
  WRITE (stdout,*) ''
  CALL print_clock ( 'ep-interp' )
  CALL print_clock ( 'ephwann' )
  !
  END SUBROUTINE ephwann_shuffle



!-------------------------------------------
SUBROUTINE epw_write (iunepwdata)
!-------------------------------------------
  !
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nbndsub, iunepmatwp, vme, eig_read
  USE pwcom,     ONLY : ef
  USE el_phon,   ONLY : nrr_k, nrr_q, chw, rdw, epmatwp, cdmew, cvmew, chw_ks
  USE phcom,     ONLY : nmodes  
  USE mp_global, ONLY : mpime, nproc, my_pool_id
  USE io_global, ONLY : ionode_id
  USE io_files,  ONLY : find_free_unit
  USE mp, ONLY : mp_barrier
  !
  implicit none
  integer   ::  iunepwdata, iunvmedata, iundmedata, iunksdata, &
       ibnd, jbnd, jmode, imode, irk, irq, ipol
  !
  WRITE(6,'(/5x,"Writing Hamiltonian, Dynamical matrix and EP vertex in Wann rep to file"/)')
  !
#ifdef __PARA
     IF (mpime.eq.ionode_id) THEN
#endif     
       !
       OPEN(unit=iunepwdata,file='epwdata.fmt')
       iundmedata = find_free_unit()
       OPEN(unit=iundmedata,file='dmedata.fmt')
       IF (vme) iunvmedata = find_free_unit()
       IF (vme) OPEN(unit=iunvmedata,file='vmedata.fmt')
       IF (eig_read) iunksdata = find_free_unit()
       IF (eig_read) OPEN(unit=iunksdata,file='ksdata.fmt')
       WRITE (iunepwdata,*) ef
       WRITE (iunepwdata,*) nbndsub, nrr_k, nmodes, nrr_q
       DO ibnd = 1, nbndsub
        DO jbnd = 1, nbndsub
         DO irk = 1, nrr_k
          WRITE (iunepwdata,*) chw(ibnd,jbnd,irk)
          IF (eig_read) WRITE (iunksdata,*) chw_ks(ibnd,jbnd,irk)
         ENDDO
        ENDDO
       ENDDO
       !
       DO imode = 1, nmodes
        DO jmode = 1, nmodes
         DO irq = 1, nrr_q
          WRITE (iunepwdata,*) rdw(imode,jmode,irq) 
         ENDDO
        ENDDO
       ENDDO
       !
       DO ibnd = 1, nbndsub
          DO jbnd = 1, nbndsub
             DO irk = 1, nrr_k
                DO ipol = 1,3
                   WRITE (iundmedata,*) cdmew(ipol, ibnd,jbnd,irk)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO ibnd = 1, nbndsub
          DO jbnd = 1, nbndsub
             DO irk = 1, nrr_k
                DO ipol = 1,3
                   IF (vme) WRITE (iunvmedata,*) cvmew(ipol, ibnd,jbnd,irk)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !
       DO irq = 1, nrr_q
          !
          ! direct read epmatwp for this irq 
          CALL rwepmatwp ( nbndsub, nrr_k, nmodes, irq, iunepmatwp, -1)
          !
          DO ibnd = 1, nbndsub
             DO jbnd = 1, nbndsub
                DO irk = 1, nrr_k
                   DO imode = 1, nmodes
                      WRITE (iunepwdata,*) epmatwp(ibnd,jbnd,irk,imode) 
                   ENDDO
                ENDDO
             ENDDO
          ENDDO
          !
       ENDDO
       CLOSE(iunepwdata)
       CLOSE(iundmedata)
       IF (vme) CLOSE(iunvmedata)
       IF (eig_read) CLOSE(iunksdata)
       !
#ifdef __PARA
    ENDIF
    CALL mp_barrier()
#endif     
     !

!---------------------------------
END SUBROUTINE epw_write
!---------------------------------


!---------------------------------
SUBROUTINE epw_read(iunepwdata)
!---------------------------------
  USE kinds,     ONLY : DP
  USE epwcom,    ONLY : nbndsub, iunepmatwp, wepexst, vme, eig_read
  USE pwcom,     ONLY : ef
  USE el_phon,   ONLY : nrr_k, nrr_q, chw, rdw, epmatwp, &
       cdmew, cvmew, chw_ks
  USE phcom,     ONLY : nmodes  
  USE mp_global, ONLY : mpime, nproc, my_pool_id, &
     intra_pool_comm, inter_pool_comm, root_pool
  USE io_global, ONLY : ionode_id
  USE mp,        ONLY : mp_barrier, mp_bcast
  USE io_files,  ONLY : find_free_unit
  !
  implicit none
  !
  integer   ::  iunepwdata, iunvmedata, iundmedata, iunksdata, &
       ibnd, jbnd, jmode, imode, irk, irq, indold, indnew, ipol
  integer :: ios
     !
     WRITE(6,'(/5x,"Reading Hamiltonian, Dynamical matrix and EP vertex in Wann rep from file"/)')
     call flush(6)
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      !
      OPEN(unit=iunepwdata,file='epwdata.fmt',status='old',iostat=ios)
      iundmedata = find_free_unit()
      OPEN(unit=iundmedata,file='dmedata.fmt',status='old',iostat=ios)
      IF (eig_read) iunksdata = find_free_unit()
      IF (eig_read) OPEN(unit=iunksdata,file='ksdata.fmt',status='old',iostat=ios)
      IF (vme) iunvmedata = find_free_unit()
      IF (vme) OPEN(unit=iunvmedata,file='vmedata.fmt',status='old',iostat=ios)
      IF (ios /= 0) call errore ('ephwann_shuffle', 'error opening epwdata.fmt',iunepwdata)
      READ (iunepwdata,*) ef
      READ (iunepwdata,*) nbndsub, nrr_k, nmodes, nrr_q
      ! 
#ifdef __PARA
    ENDIF
    CALL mp_bcast (ef, ionode_id, inter_pool_comm)
    CALL mp_bcast (ef, root_pool, intra_pool_comm)
    !
    CALL mp_bcast (nbndsub, ionode_id, inter_pool_comm)
    CALL mp_bcast (nbndsub, root_pool, intra_pool_comm)
    !
    CALL mp_bcast (nrr_k, ionode_id, inter_pool_comm)
    CALL mp_bcast (nrr_k, root_pool, intra_pool_comm)
    !
    CALL mp_bcast (nmodes, ionode_id, inter_pool_comm)
    CALL mp_bcast (nmodes, root_pool, intra_pool_comm)
    !
    CALL mp_bcast (nrr_q, ionode_id, inter_pool_comm)
    CALL mp_bcast (nrr_q, root_pool, intra_pool_comm)
#endif
    !
    IF (.not. ALLOCATED(epmatwp)) ALLOCATE ( epmatwp ( nbndsub, nbndsub, nrr_k, nmodes) )
    IF (.not. ALLOCATED(chw)    ) ALLOCATE ( chw ( nbndsub, nbndsub, nrr_k )            )
    IF (.not. ALLOCATED(chw_ks) ) ALLOCATE ( chw_ks ( nbndsub, nbndsub, nrr_k )            )
    IF (.not. ALLOCATED(rdw)    ) ALLOCATE ( rdw( nmodes,  nmodes,  nrr_q )             )
    IF (.not. ALLOCATED(cdmew)  ) ALLOCATE ( cdmew   ( 3, nbndsub, nbndsub, nrr_k )     )
    IF (vme .and. (.not.ALLOCATED(cvmew))  ) ALLOCATE ( cvmew   ( 3, nbndsub, nbndsub, nrr_k )     )

    !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      !
      DO ibnd = 1, nbndsub
       DO jbnd = 1, nbndsub
        DO irk = 1, nrr_k
         READ (iunepwdata,*) chw(ibnd,jbnd,irk)
         IF (eig_read) READ (iunksdata,*) chw_ks(ibnd,jbnd,irk)
        ENDDO
       ENDDO
      ENDDO
      !     
      DO imode = 1, nmodes
       DO jmode = 1, nmodes
        DO irq = 1, nrr_q
         READ (iunepwdata,*) rdw(imode,jmode,irq) 
        ENDDO
       ENDDO
      ENDDO
      !
       DO ibnd = 1, nbndsub
          DO jbnd = 1, nbndsub
             DO irk = 1, nrr_k
                DO ipol = 1,3
                   READ (iundmedata,*) cdmew(ipol, ibnd,jbnd,irk)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       DO ibnd = 1, nbndsub
          DO jbnd = 1, nbndsub
             DO irk = 1, nrr_k
                DO ipol = 1,3
                   IF (vme) READ (iunvmedata,*) cvmew(ipol, ibnd,jbnd,irk)
                ENDDO
             ENDDO
          ENDDO
       ENDDO
       !
#ifdef __PARA
    ENDIF
    !
    CALL mp_bcast (chw, ionode_id, inter_pool_comm)
    CALL mp_bcast (chw, root_pool, intra_pool_comm)
    !
    IF (eig_read) CALL mp_bcast (chw_ks, ionode_id, inter_pool_comm)
    IF (eig_read) CALL mp_bcast (chw_ks, root_pool, intra_pool_comm)
    !
    CALL mp_bcast (rdw, ionode_id, inter_pool_comm)
    CALL mp_bcast (rdw, root_pool, intra_pool_comm)
    !
    CALL mp_bcast (cdmew, ionode_id, inter_pool_comm)
    CALL mp_bcast (cdmew, root_pool, intra_pool_comm)
    !
    IF (vme) CALL mp_bcast (cvmew, ionode_id, inter_pool_comm)
    IF (vme) CALL mp_bcast (cvmew, root_pool, intra_pool_comm)
    !
#endif
    !
    IF (.not.wepexst) THEN
    DO irq = 1, nrr_q
      !
#ifdef __PARA
      IF (mpime.eq.ionode_id) THEN
#endif
         !
         DO ibnd = 1, nbndsub
          DO jbnd = 1, nbndsub
           DO irk = 1, nrr_k
            DO imode = 1, nmodes
             READ (iunepwdata,*) epmatwp(ibnd,jbnd,irk,imode)
            ENDDO
           ENDDO
          ENDDO
         ENDDO
         !
#ifdef __PARA
      ENDIF
      !
      CALL mp_bcast (epmatwp, ionode_id, inter_pool_comm)
      CALL mp_bcast (epmatwp, root_pool, intra_pool_comm)
      !
#endif
      !
      ! direct write of epmatwp irq 
      CALL rwepmatwp ( nbndsub, nrr_k, nmodes, irq, iunepmatwp, +1)
      !
      IF (irq.eq.1) THEN
         WRITE(6,'(/5x,"Loading epwdata.fmt: ")',advance='no')
        indold = 0
      ENDIF
      indnew = nint(float(irq)/float(nrr_q)*43)
      IF (indnew.ne.indold) WRITE(6,'(a)',advance='no') '#'
      indold = indnew
      CALL flush(6) ! ONLY on opteron
      !
    ENDDO
    ENDIF
#ifdef __PARA
    CALL mp_barrier()
    IF (mpime.eq.ionode_id) &
#endif
    CLOSE(iunepwdata)
    CLOSE(iundmedata)
    IF (vme) CLOSE(iunvmedata)
    !
    WRITE(6,'(/5x,"Finished reading Wann rep data from file"/)')
    !
!---------------------------------
END SUBROUTINE epw_read
!---------------------------------

!---------------------------------
SUBROUTINE collect_q_list (do_q, nxqf)
  !
  USE kinds,         ONLY : DP
  USE epwcom, ONLY : nbndsub, indabs, nsmear, delta_smear, emaxabs, eminabs, &
       degaussw, twophoton
  USE el_phon,       ONLY : nrr_k, nrr_q, cu, cuq, irvec, ndegen_k, ndegen_q,   &
                            wslen, chw,  wf, &
                            etf, etfq, xqf, xkf, wkf, wqf, &
                            nksf, epf17, nksqf, et_ks
#ifdef __PARA
  USE mp,            ONLY : mp_sum
  USE mp_global,     ONLY : inter_pool_comm, my_pool_id
#endif
  !
  implicit none
  !
  integer :: do_q(nxqf), nxqf
  integer :: total_pairs, ik, iq, ikk, ikq, ibnd, jbnd
  real(kind=DP) :: vbm, cbm, dvc, xkk(3), xkq(3), &
       etk (nbndsub), etkq (nbndsub)
  complex(kind=DP) :: cufkk ( nbndsub, nbndsub)

  real(kind=DP), PARAMETER :: ryd2ev = 13.6058
  !
!  IF (indabs.or.twophoton) THEN
!     do_q = 0
!     total_pairs = 0
!     vbm = 5.36 / ryd2ev ! hard coded for Si
!     !
!     DO ik = 1, nksqf
!        !
!        ikk = 2 * ik - 1
!        ikq = ikk + 1
!        !
!        xkk = xkf(:, ikk)
!        !
!        CALL hamwan2bloch &
!             ( nbndsub, nrr_k, irvec, ndegen_k, xkk, cufkk, etk, chw)
!        !
!        DO iq = 1, nxqf
!           !
!           ! conduction band maximum (bad name, I know), is the top of the conduction band +
!           ! the largest value in the indabs spectra + a few times the multiple of the
!           ! smearing used on KS eigenvalues to make sure not to miss anything at the edge
!           ! of the window
!           dvc = (emaxabs)/ryd2ev + 6 * ((nsmear-1) * delta_smear + degaussw) + MAXVAL(wf(:,iq))
!           cbm = vbm + dvc
!           !
!           xkq = xkk + xqf (:, iq)
!           !
!           CALL hamwan2bloch &
!                ( nbndsub, nrr_k, irvec, ndegen_k, xkq, cufkk, etkq, chw)
!           !
!           IF (iq.eq.1) total_pairs = total_pairs + nbndsub*nbndsub
!           DO ibnd = 1, nbndsub
!              DO jbnd = 1, nbndsub
!                 IF (etk(ibnd).lt.vbm .and. etkq(jbnd).gt.vbm .and. etkq(jbnd).lt.cbm .and. abs(etk(ibnd)-etkq(jbnd)).le.dvc) THEN
!                    do_q(iq) = do_q(iq) + 1
!                 ENDIF
!                 IF (etkq(jbnd).lt.vbm .and. etk(ibnd).gt.vbm .and. etk(ibnd).lt.cbm .and.  abs(etk(ibnd)-etkq(jbnd)).le.dvc) THEN
!                    do_q(iq) = do_q(iq) + 1
!                 ENDIF
!                 !
!              ENDDO
!           ENDDO
!        ENDDO ! iq loop
!     ENDDO ! ik loop
!#ifdef __PARA
!     CALL mp_sum(do_q, inter_pool_comm)
!     CALL mp_sum(total_pairs, inter_pool_comm)
!#endif
!     do_q(1) = 1
!     do_q(nxqf) = 1
!     total_pairs = nxqf
!     DO iq = 1, nxqf
!        IF ( do_q(iq) .eq. 0) total_pairs = total_pairs - 1
!     ENDDO
!     WRITE(6,'(5x,a,i6,a,i6)') "Of ", nxqf, " phonon wavevectors, we are using ", total_pairs
!  ELSE
     do_q = 1
!  ENDIF
!---------------------------------
END SUBROUTINE collect_q_list
!---------------------------------


!---------------------------------
SUBROUTINE mem_size(ibndmin, ibndmax, nmodes, nksqf, nxqf, fly) 
!---------------------------------
!
!  SUBROUTINE estimates the amount of memory taken up by 
!  the <k+q| dV_q,nu |k> on the fine meshes and prints 
!  out a useful(?) message   
!
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  !
  implicit none
  !
  integer :: imelt, ibndmin, ibndmax, nmodes, nksqf, nxqf
  real(kind=DP)    :: rmelt
  character (len=256) :: chunit
  logical :: fly
  !
  imelt = (ibndmax-ibndmin+1)**2 * nmodes * nksqf
  IF (.not.fly) imelt = imelt * nxqf
  rmelt = imelt * 8 / 1048576.d0 ! 8 bytes per number, value in Mb
  IF (rmelt .lt. 1000.0 ) THEN
     chunit =  ' Mb '
     IF (rmelt .lt. 1.0 ) THEN
        chunit = ' Kb '
        rmelt  = rmelt * 1024.d0
     ENDIF
  ELSE
     rmelt = rmelt / 1024.d0
     chunit = ' Gb '
  ENDIF
  WRITE(stdout,'(/,5x,a, i13, a,f7.2,a,a)') "Number of ep-matrix elements per pool :", &
       imelt, " ~= ", rmelt, trim(chunit), " (@ 8 bytes/ DP)"
  !

!---------------------------------
END SUBROUTINE mem_size
!---------------------------------
