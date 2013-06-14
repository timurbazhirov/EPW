  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !---------------------------------------------------------------------
  SUBROUTINE elphel2_shuffle (npe, imode0, dvscfins, gmapsym, eigv, isym, invs, xq0, timerev)
  !---------------------------------------------------------------------
  !
  !      Calculation of the electron-phonon matrix elements el_ph_mat
  !      <\psi(k+q)|dV_{SCF}/du^q_{i a}|\psi(k)>
  !
  !      Written by Feliciano Giustino based on the routine elphel. 
  !      Main difference w.r.t. to original routine is gauge fixing, 
  !      shuffle (umklapp) mode and all-q implementation.
  !
  !      Shuffle mode implemented on may 7 2006
  !
  !      No ultrasoft now
  !
  !      Nota Bene: this SUBROUTINE is intended only for one proc per pool, 
  !      i.e. with no G-vector parallelization (some work on the igkq is 
  !      required for that in the g-mapping)
  !
  !      In order to allow a pool reading the wfc file of another
  !      pool, I had to modify the bound npwx in PW/n_plane_waves.f90
  !      which is now the max across all pools. In this way lrwfc is
  !      the same for all pools.
  !
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  !
#ifdef __PARA
  USE mp_global,   ONLY : nproc, my_pool_id, nproc_pool,    & 
                          intra_image_comm, intra_pool_comm, &
                          me_pool, root_pool, mpime, inter_pool_comm
  USE mp,          ONLY : mp_barrier, mp_bcast, mp_put,mp_sum
#endif
  USE kinds,       ONLY : DP
  USE io_global,   ONLY : stdout, ionode_id, ionode
  USE wavefunctions_module,  ONLY: evc
  USE io_files,    ONLY: prefix, iunigk, tmp_dir, nd_nmbr
  USE wvfct,       ONLY : npwx
  USE pwcom,       ONLY : current_spin, isk, nls, tpiba, g, &
                          lsda, nrx1s, nrx2s, nrx3s, &
                          nr1s, nr2s, nr3s, nbnd, npw, xk, &
                          nrxxs, nspin, ngm, s, nkb,vkb, et, igk, nks
  USE phcom,       ONLY : alphap, becp1,  u, lrwfc, dvpsi, nlcc_any, iuwfc, &
                          lgamma, npwq, igkq, evq, xq
  USE epwcom,      ONLY : tphases
  USE el_phon,     ONLY : shift, gmap, el_ph_mat, umat, umatq, &
                          umat_all, xk_all, et_all, xkq, etq
  USE control_flags, ONLY : iverbosity
  USE becmod,      ONLY : calbec
  USE klist,       ONLY : nkstot
  ! 
  implicit none
  !
  integer :: npe, imode0, i
  !  pert for this irrep
  !  mode number
  !  unit for e-ph matrix elements
  complex(kind=DP) :: dvscfins (nrxxs, nspin, npe)
  !  delta scf potential

  logical :: exst
  !
  ! work variables
  !
  integer :: ik, ipert, mode, ibnd, jbnd, ir, ig, nkq, ipool, &
       ik0, igkq_tmp (npwx), imap, &
       ipooltmp, nkq_abs, ism1, ipol
  complex(kind=DP), ALLOCATABLE :: aux1 (:), elphmat (:,:,:), eptmp (:,:), aux2(:,:)
  complex(kind=DP) :: ZDOTC
  real(kind=DP) :: xktmp(3), invsxk(3)
  complex(kind=DP), PARAMETER :: czero = (0.d0, 0.d0), cone = (1.d0, 0.d0), ci = (0.d0, 1.d0)
  !
  ! variables for folding of k+q grid
  !
!  REAL(kind=DP) :: g0vec_all_r(3,27) 
  REAL(kind=DP) :: g0vec_all_r(3,125) 
 
  !   G-vectors needed to fold the k+q grid into the k grid, cartesian coord.
  INTEGER :: ng0vec, ngxx           
  !   number of inequivalent such translations
  !   bound for the allocation of the array gmap
  !
  ! variables for rotating the wavefunctions (in order to USE q in the irr wedge)
  !
  INTEGER :: gmapsym ( ngm, 48 ), isym, invs(48)
  ! the map G->S(G)
  ! the symmetry which generates the current q in the star
  ! index of the inverse operation
  complex(kind=DP) :: eigv (ngm, 48)
  real(kind=DP) :: xq0(3)
  real(kind=DP) :: zero_vect(3) 
  integer :: nkk, nkk_abs
  !  the fractional traslation
  !  work variable
  !  the first q in the star (cartesian)
  real(kind=DP), PARAMETER :: twopi = 6.28318530717959
  logical :: timerev
  !  true if we are using time reversal
  !
  ALLOCATE (elphmat( nbnd, nbnd, npe), eptmp (nbnd, nbnd), aux1(nrxxs),&
            aux2( npwx, nbnd) )
  zero_vect = 0.d0
  IF (ALLOCATED(xkq) ) DEALLOCATE (xkq)                  
  IF (.not. ALLOCATED(xkq) ) ALLOCATE (xkq (3, nkstot) ) 
#ifdef __PARA
  IF (nproc_pool>1) call errore &
    ('elphel2_shuffle', 'ONLY one proc per pool in shuffle mode', 1)
#endif

  IF (.not.lgamma) THEN
    !
    ! setup for k+q folding
    !
    CALL kpointdivision ( ik0 )
    CALL readgmap ( nkstot, ngxx, ng0vec, g0vec_all_r )
    !
    IF (imode0.eq.0 .and. iverbosity.eq.1) WRITE(stdout, 5) ngxx
5   FORMAT (5x,'Estimated size of gmap: ngxx =',i5)
    !
  ENDIF
  !
  ! close all sequential files in order to re-open them as direct access
  ! close all .wfc files in order to prepare shuffled read
  !
  CLOSE (unit = iunigk, status = 'keep')
  CLOSE (unit = iuwfc,  status = 'keep')
#ifdef __PARA
  ! never remove this barrier
  CALL mp_barrier()
#endif
  !
  ism1 = invs (isym)
  !
!  ! fractional traslation, cartesian coord
!  !
!  v(1) = - float (ftau (1, isym) )  / float (nr1)
!  v(2) = - float (ftau (2, isym) )  / float (nr2)
!  v(3) = - float (ftau (3, isym) )  / float (nr3)
!  call cryst_to_cart (1, v, bg, +1)
  !
  !  Start the loop over the k-points
  !
  DO ik = 1, nks
     !
     !
     ! find index, and possibly pool, of k+q 
     ! the index nkq (nkq_abs) takes into account the even/odd ordering 
     ! of the nscf calc
     ! we also redefine the ikq points and the corresponding energies
     ! (we need to make sure that xk(:,ikq) is really k+q for the KB projectors
     ! below and that also that the eigenvalues are taken correctly in ephwann)
     !
#ifdef __PARA
     ipooltmp= my_pool_id+1
#endif
     !
     !
     CALL ktokpmq ( xk (:, ik), xq, +1, ipool, nkq, nkq_abs)
     !
     !   we define xkq(:,ik) and etq(:,ik) for the current xq
     !
     IF (allocated(etq)) DEALLOCATE(etq)
     ALLOCATE (etq (nbnd, nks) )
     xkq(:, ik) = xk_all(:, nkq_abs)
     etq(:, ik) = et_all(:, nkq_abs) 
     !
     !
     ! in serial execution ipool is not USEd in the called SUBROUTINEs, 
     ! in parallel mypool is for k and ipool is for k+q
     !
     CALL readwfc (  ipooltmp, ik, evc)
     CALL readigk (  ipooltmp, ik, npw, igk)
     !
     CALL readwfc ( ipool, nkq, evq)
     CALL readigk ( ipool, nkq, npwq, igkq)
     !
#ifdef __PARA
     IF (.not.lgamma.and.nks.gt.1.and.maxval(igkq(1:npwq)).gt.ngxx) &
          CALL errore ('elphel2_shuffle', 'ngxx too small', 1 )
#endif
     !
     ! ----------------------------------------------------------------
     ! Set the gauge for the eigenstates: unitary transform and phases
     ! ----------------------------------------------------------------
     !
     ! With this option, different compilers and different machines
     ! should always give the same wavefunctions.
     !
     CALL ktokpmq ( xk(:,ik),  zero_vect, +1, ipool, nkk, nkk_abs)
     CALL ktokpmq ( xkq(:,ik), zero_vect, +1, ipool, nkk, nkq_abs)
     IF ( .not. ALLOCATED (umat) ) ALLOCATE ( umat(nbnd,nbnd,nks) )
     IF ( .not. ALLOCATED (umatq) ) ALLOCATE ( umatq(nbnd,nbnd,nks) )
     umat(:,:,ik)  = umat_all(:,:, nkk_abs)
     umatq(:,:,ik) = umat_all(:,:, nkq_abs)
    !
    ! the k-vector needed for the KB projectors
    xktmp = xkq (:, ik)
    !
    ! --------------------------------------------------
    !   Fourier translation of the G-sphere igkq
    ! --------------------------------------------------
    !
    !
    !  Translate by G_0 the G-sphere where evq is defined, 
    !  none of the G-points are lost.
    !
    DO ig = 1, npwq
       imap = ng0vec * ( igkq (ig) - 1 ) + shift( ik+ik0 )
       igkq_tmp (ig) = gmap( imap )
       !  the old matrix version... 
       !  igkq_tmp (ig) = gmap( igkq (ig),  shift( ik+ik0 ) )
    ENDDO
    igkq = igkq_tmp
    !
    !  find k+q from k+q+G_0
    !  (this is needed in the calculation of the KB terms
    !  for nonlocal pseudos)
    !
!      xktmp = xk (:, ikq) - g0vec_all_r ( :, shift( ik+ik0 ) )
    xktmp = xkq (:, ik) - g0vec_all_r ( :, shift( ik+ik0 ) )
    !
    !
    ! ---------------------------------------------------------------------
    ! phase factor arising from fractional traslations
    ! ---------------------------------------------------------------------
    !
    !  u_{k+q+G_0} carries and additional factor e^{i G_0 v}
    !
!   arg = twopi * dot_product ( g0vec_all_r ( :, shift( ik+ik0 ) ), v ) 
!   eig0v = dcmplx ( cos(arg), sin(arg) )
!   call fractrasl ( nbnd, npwq, npwx, ngm, igkq, evq, eigv (:,isym), eig0v)
    !
    CALL fractrasl ( nbnd, npw,  npwx, ngm, igk,  evc, eigv (:,isym), cone)
    CALL fractrasl ( nbnd, npwq, npwx, ngm, igkq, evq, eigv (:,isym), cone)
    !
    ! ---------------------------------------------------------------------
    ! wave function rotation to generate matrix elements for the star of q
    ! ---------------------------------------------------------------------
    !
    ! ps. don't USE npwx instead of npw, npwq since the unUSEd elements
    ! may be large and blow up gmapsym (personal experience)
    !
    igk (1:npw ) = gmapsym ( igk (1:npw ), isym )
    igkq(1:npwq) = gmapsym ( igkq(1:npwq), isym )
    !
    ! ---------------------------------------------------------------------
    ! In case of time-reversal symmetry we take the c.c. of the wfs
    ! ---------------------------------------------------------------------
    ! 
    ! we need the c.c. of dV_hartree + dV_loc + dV_NL, but it is easier to do c.c. of 
    ! the (phase-fixed) wfs and then take the c.c. of the final matrix element
    !
    IF (timerev) THEN
       evc = CONJG ( evc )
       evq = CONJG ( evq )
       umat(:,:,ik) = CONJG ( umat(:,:,ik) )
       umatq(:,:,ik) = CONJG ( umatq(:,:,ik) )
    ENDIF
    !
    !
    ! In dvqpsi_us_ONLY3 we need becp1 and alphap for the rotated wfs. 
    ! The other quantities (deeq and qq) do not depend on the wfs, in
    ! particular in the KB case (not ultrasoft), the deeq's are the
    ! unscreened coefficients, and the qq's are zero.
    !
    ! For the KB part, remember dV_NL[q_0] ~ |S^-1(k)+q_0> <S^-1(k)|
    ! the total momentum transfer must be q_0 and the rotation 
    ! tranforms k+Sq_0 into S^-1(k)+q_0, k into S^-1(k)
    ! [see Eqs. (A9),(A14) Baroni et al. RMP]
    !
    CALL rotate_cart( xk(:,ik), s (1, 1, isym), invsxk)
    !
    ! here we generate vkb on the igk() set and for k ...
    CALL init_us_2 (npw, igk, invsxk, vkb)
    !
    ! ... and we recompute the becp terms with the wfs (rotated through igk)
    !
    CALL calbec(npw,vkb,evc,becp1(:,:,ik)) 
    !
    ! we also recompute the derivative of the becp terms with the (rotated) wfs
    !
    DO ipol = 1, 3
      DO ibnd = 1, nbnd
        DO ig = 1, npw
           aux2(ig,ibnd) = evc(ig,ibnd) * tpiba * ci * ( invsxk(ipol) + g(ipol,igk(ig)) )
        ENDDO
      ENDDO
      CALL calbec(npw,vkb,aux2,alphap(:,:,ipol,ik)) 
    ENDDO
    !
    ! now we generate vkb on the igkq() set becaUSE dvpsi is needed on that set
    ! we need S^-1(k)+q_0 in the KB projector: total momentum transfer must be q_0
    !
    xktmp = invsxk + xq0
    CALL init_us_2 (npwq, igkq, xktmp, vkb)
    !
    ! --------------------------------------------------
    !   Calculation of the matrix element
    ! --------------------------------------------------
    !
    IF (lsda) current_spin = isk (ik)
    !
    DO ipert = 1, npe
       !
       !  recalculate dvbare_q*psi_k 
       !  the call to dvqpsi_us3 differs from the old one to dvqpsi_us 
       !  ony by the xktmp passed. 
       !
       !  we have to USE the first q in the star in the dvqpsi_us3 call below (xq0)
       !  
       mode = imode0 + ipert
       CALL dvqpsi_us3 (ik, mode, u (1, mode), nlcc_any, xktmp, xq0 ) 
       !
       !  calculate dvscf_q*psi_k
       !
       DO ibnd = 1, nbnd
          aux1(:) = (0.d0, 0.d0)
          DO ig = 1, npw
             aux1 (nls (igk (ig) ) ) = evc ( ig, ibnd)
          ENDDO
          CALL cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
          DO ir = 1, nrxxs
             aux1 (ir) = aux1 (ir) * dvscfins (ir, current_spin, ipert)
          ENDDO
          CALL cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
          DO ig = 1, npwq
            dvpsi (ig, ibnd) = dvpsi (ig, ibnd) + aux1 (nls (igkq (ig) ) )
          ENDDO
       END DO
       !
       ! calculate elphmat(j,i)=<psi_{k+q,j}|dvscf_q*psi_{k,i}> for this pertur
       !
       DO ibnd =1, nbnd
         DO jbnd = 1, nbnd
           elphmat (jbnd, ibnd, ipert) = &
              ZDOTC (npwq, evq(1, jbnd), 1, dvpsi(1, ibnd), 1)
         ENDDO
       ENDDO
       !
    ENDDO
#ifdef __PARA
    CALL mp_sum(elphmat, intra_pool_comm)
#endif
    !
    !  Rotate elphmat with the gauge matrices (this should be equivalent 
    !  to calculate elphmat with the truely rotated eigenstates)

    DO ipert = 1, npe
       !
       ! the two zgemm call perform the following ops:
       !  elphmat = umat(k+q)^\dagger * [ elphmat * umat(k) ]
       !
       CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, elphmat(:,:,ipert), & 
            nbnd, umat(:,:,ik), nbnd, czero, eptmp, nbnd)
       CALL zgemm ('c', 'n', nbnd, nbnd, nbnd, cone, umatq(:,:,ik), & 
            nbnd, eptmp, nbnd, czero, elphmat(:,:,ipert), nbnd)
       !
    ENDDO
    !
    ! remember we have to take the c.c. of the final matrix element 
    ! when we USE time-reversal
    !
    IF (timerev) elphmat = CONJG ( elphmat ) 
    !
    !  save eph matrix elements into el_ph_mat
    !
    DO ipert = 1, npe
      DO jbnd = 1, nbnd
      DO ibnd = 1, nbnd
        el_ph_mat (ibnd, jbnd, ik, ipert + imode0) = & 
          elphmat (ibnd, jbnd, ipert)
      ENDDO
      ENDDO
    ENDDO
    !
    !
  ENDDO
  !
  !  restore original configuration of files
  !
  CALL seqopn (iunigk, 'igk', 'unformatted', exst)
  CALL diropn (iuwfc, 'wfc', lrwfc, exst) 
#ifdef __PARA
  ! never remove this barrier - > insures that wfcs are restored to each pool before moving on
  CALL mp_barrier()
#endif
  !
  DEALLOCATE (elphmat, eptmp, aux1, aux2)
  DEALLOCATE (gmap, shift)
  !
  END SUBROUTINE elphel2_shuffle
  !
  !------------------------------------------------------------
  SUBROUTINE fractrasl ( nbnd, npw, npwx, ngm, igk, evc, eigv1, eig0v)
  !------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  implicit none
  integer :: nbnd, npw, npwx, ngm, igk (npwx)
  complex(kind=DP) :: evc (npwx, nbnd), eigv1 (ngm), eig0v
  integer :: ig, ibnd
  ! 
  DO ibnd = 1, nbnd
    DO ig = 1, npw
      evc (ig, ibnd) = evc (ig, ibnd) * eigv1 ( igk(ig) ) * eig0v
    ENDDO
  ENDDO

  END SUBROUTINE fractrasl
  !
  !------------------------------------------------------------
  SUBROUTINE rotate_cart( x, s, sx)
  !------------------------------------------------------------
  !
  ! a simple symmetry operation in cartesian coordinates 
  ! ( s is integer and in crystal coord!)
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  implicit none
  real(kind=DP) :: x(3), sx(3), xcrys(3)
  integer :: s(3,3),i
  !
  xcrys = x
  CALL cryst_to_cart (1, xcrys, at, -1)
  DO i = 1, 3
     sx (i) = float(s (i, 1)) * xcrys (1) &
            + float(s (i, 2)) * xcrys (2) &
            + float(s (i, 3)) * xcrys (3)
  ENDDO
  CALL cryst_to_cart (1, sx, bg, +1)
  !

  END SUBROUTINE rotate_cart
  !------------------------------------------------------------
  ! 
  ! I USEd the following to check the correctness of the symmetry
  ! operation on the wavefunctions: symmetrization in real-space 
  !
  !     integer :: i,j,k,ri,rj,rk
  !     complex(kind=DP), allocatable :: aux2(:)
  !     allocate (aux2(nrxxs))
  !
  !     do ibnd = 1, nbnd
  !       !
  !       aux1 = czero
  !       aux2 = czero
  !       do ig = 1, npw
  !          aux1 (nls (igk (ig) ) ) = evc ( ig, ibnd)
  !       enddo
  !       call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
  !       do k = 1, nr3
  !        do j = 1, nr2
  !         do i = 1, nr1
  !            call ruotaijk (s (1, 1, isym), ftau (1, isym), &
  !                i, j, k, nr1, nr2, nr3, ri, rj, rk )
  !            aux2 (1 + ( i-1) + nr1*( j-1) + nr1*nr2*( k-1)) = &
  !            aux1 (1 + (ri-1) + nr1*(rj-1) + nr1*nr2*(rk-1)) 
  !         enddo
  !        enddo
  !       enddo
  !       call cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
  !       do ig = 1, npw
  !         evc ( ig, ibnd) = aux2 (nls (igk (ig) ) )
  !       enddo
  !       !
  !       aux1 = czero
  !       aux2 = czero
  !       do ig = 1, npwq
  !          aux1 (nls (igkq (ig) ) ) = evq ( ig, ibnd)
  !       enddo
  !       call cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
  !       do k = 1, nr3
  !        do j = 1, nr2
  !         do i = 1, nr1
  !            call ruotaijk (s (1, 1, isym), ftau (1, isym), &
  !                i, j, k, nr1, nr2, nr3, ri, rj, rk )
  !            aux2 (1 + ( i-1) + nr1*( j-1) + nr1*nr2*( k-1)) = &
  !            aux1 (1 + (ri-1) + nr1*(rj-1) + nr1*nr2*(rk-1)) 
  !         enddo
  !        enddo
  !       enddo
  !       call cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
  !       do ig = 1, npwq
  !         evq ( ig, ibnd) = aux2 (nls (igkq (ig) ) )
  !       enddo
  !       !
  !     enddo
  !
  !  Instead of rotating the wfs, we can also perform the operation
  !  directly on the dvscf: here is how it goes. Note that in order
  !  to use this we have to switch off the local and the nonlocal
  !  pseudo contributions (i.e. we keep only the hartree term)
  !  [set dvpsi = czero after the dvqpsi_us3 call]
  !       do k = 1, nr3
  !        do j = 1, nr2
  !         do i = 1, nr1
  !            call ruotaijk (s (1, 1, ism1), ftau (1, ism1), &
  !                i, j, k, nr1, nr2, nr3, ri, rj, rk )
  !            aux2 (1 + ( i-1) + nr1*( j-1) + nr1*nr2*( k-1)) = &
  !            dvscfins (1 + (ri-1) + nr1*(rj-1) + nr1*nr2*(rk-1), current_spin, ipert)
  !         enddo
  !        enddo
  !       enddo
  !       do ir = 1, nrxxs
  !         aux1 (ir) = aux1 (ir) * aux2 (ir)
  !       enddo
  !
