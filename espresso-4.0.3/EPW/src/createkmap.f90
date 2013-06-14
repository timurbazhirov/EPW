  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE createkmap ( xq )
  !-----------------------------------------------------------------------
  !  
  !  This subroutine is called from elphon_shuffle_wrap for each
  !  nq1*nq2*nq3 phonon on the coarse mesh.    
  !
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  USE klist, ONLY : nkstot, nks, xk
  USE pwcom, ONLY : nk1, nk2, nk3
  USE epwcom, ONLY : xk_cryst
  USE io_global, ONLY : stdout
  USE io_files, ONLY : prefix
  USE klist_epw, ONLY : kmap
  USE kfold
  USE control_flags, ONLY : iverbosity
  USE control_ph, ONLY : lgamma
#ifdef __PARA
  USE mp_global, ONLY : my_pool_id, mpime
  USE mp, ONLY : mp_barrier
#endif
  USE el_phon, ONLY : xkq
  implicit none
  !
  real(kind=DP) :: eps, xq (3), xk_q(3)
  !
  integer :: ik, jk, j
  !
  !  obsolete: 
  !
  !  variables for folding of the k+q mesh into the kmesh @FG
  !  there may be at most 9 different G_0 to refold k+q into k
  !  of the first BZ (including G_0=0)
  !  [actually there are 3^3=27 translations but for a fixed q ONLY 
  !   2^3=8 out of these 27 are possible since q has a definite direction]
  !
  !  new:
  !
  !  I need all q-vectors int the same run, so I consider all the
  !  27 possibilities by setting ng0vec = 27 in ../PW/set_kplusq.f90 
  !
  !  newest:
  !
  !  apparently we need 125 for the cuprates.  
  ! 
  !  future:
  !
  !  maybe the size can be dynamical?  
  !
  !
  real(kind=DP) :: xx, yy, zz, xx_n, yy_n, zz_n, xx_c, yy_c, zz_c
  integer :: i, k, n, g0vec(3), ig0, iukmap, ig1, ig2, ig3
  logical :: in_the_list, bingo, found
  !
  IF (.not. ALLOCATED(xkq) ) ALLOCATE(xkq (3, nkstot) )
  !
#ifdef __PARA
  IF (mpime.eq.0) THEN
#endif
    !
    !  the first proc keeps a copy of all kpoints !
    !
    IF ( .not. ALLOCATED (shift) ) ALLOCATE ( shift(nkstot) )
    !
    !  Now fold k+q back into the k-grid for wannier interpolation.
    !  Since this is done before divide&impera, every pool has all the kpoints. @FG
    !
    ! bring q in crystal coordinates and check commensuration
    ! 
    !  loosy tolerance: not important since k+q is defined trhough nint() 
    eps = 1.d-5 
    ! 
    CALL cryst_to_cart (1,xq,at,-1)
    !
    xx = xq(1)*nk1 
    yy = xq(2)*nk2 
    zz = xq(3)*nk3 
    in_the_list = abs(xx-nint(xx)).le.eps .and. &
                  abs(yy-nint(yy)).le.eps .and. &
                  abs(zz-nint(zz)).le.eps
    IF (iverbosity.eq.1) &
      WRITE(stdout,'(a,3i3)') '  q in integer coord:',nint(xx),nint(yy),nint(zz)
    IF (.not.in_the_list) CALL errore('create_kmap','q-vec not commensurate',1)
    !
    !  bring all the k's in crystal coordinates 
    !
    CALL cryst_to_cart ( nkstot, xk, at, -1)
    !
    ! previously the 27 possible translations
    ! now we have 125 possibilities for LSCO, becaUSE the unit cell is far from cubic... sob!
    ! ng0vec must be redefined 
    !
    ng0vec = 0
!    do ig1 = -1,1
!     do ig2 = -1,1
!      do ig3 = -1,1
    DO ig1 = -2,2
     DO ig2 = -2,2
      DO ig3 = -2,2
         ng0vec = ng0vec + 1
         g0vec_all(1,ng0vec) = ig1
         g0vec_all(2,ng0vec) = ig2
         g0vec_all(3,ng0vec) = ig3
      ENDDO
     ENDDO
    ENDDO
    !
    DO ik = 1, nkstot
      !
      !  check that the k's are actually on a uniform mesh centered at gamma
      !
      xx = xk(1, ik)*nk1
      yy = xk(2, ik)*nk2
      zz = xk(3, ik)*nk3
      IF (iverbosity.eq.1) &
        WRITE(stdout,'(a,i3,a,3i3)') 'ik = ',ik,',   k   in integer coord:',nint(xx),nint(yy),nint(zz)
      in_the_list = abs(xx-nint(xx)).le.eps .and. &
                    abs(yy-nint(yy)).le.eps .and. &
                    abs(zz-nint(zz)).le.eps
      IF (.not.in_the_list) CALL errore('createkmap','is this a uniform k-mesh?',1)
      !
      !  now add the phonon wavevector and check that k+q falls again on the k grid
      !
      xk_q (:) = xk (:, ik) + xq (:)
      !
      xx = xk_q(1)*nk1
      yy = xk_q(2)*nk2
      zz = xk_q(3)*nk3
      IF (iverbosity.eq.1) &
        WRITE(stdout,'(a,i3,a,3i3)') 'ik = ',ik,',   k+q in integer coord:',nint(xx),nint(yy),nint(zz)
      in_the_list = abs(xx-nint(xx)).le.eps .and. &
                    abs(yy-nint(yy)).le.eps .and. &
                    abs(zz-nint(zz)).le.eps
      IF (.not.in_the_list) CALL errore('createkmap','k+q does not fall on k-grid',1)
      !
      !  find the index of this k+q in the k-grid
      !
      i = mod ( nint ( xx + 2*nk1), nk1 ) 
      j = mod ( nint ( yy + 2*nk2), nk2 ) 
      k = mod ( nint ( zz + 2*nk3), nk3 ) 
      IF (iverbosity.eq.1) &
        WRITE(stdout,'(a,i3,a,3i3)') 'ik = ',ik,', f-k+q in integer coord:',i,j,k
      !
      xx_n = xx
      yy_n = yy
      zz_n = zz
      !
      !  make sure xx_n, yy_n and zz_n are in 1st BZ
      !
      CALL backtoBZ( xx_n, yy_n, zz_n, nk1, nk2, nk3 )
      !
      n = 0
      found = .false. 
      DO jk = 1, nkstot
         IF ( .not. found ) THEN
            xx_c = xk_cryst(1,jk)*nk1
            yy_c = xk_cryst(2,jk)*nk2
            zz_c = xk_cryst(3,jk)*nk3
            !
            ! check that the k-mesh was defined in the positive region of 1st BZ
            !
            IF ( xx_c .lt. -eps .or. yy_c .lt. -eps .or. zz_c .lt. -eps ) &
               call errore('ktokpmq','coarse k-mesh needs to be strictly positive in 1st BZ',1)
            !
            bingo = nint(xx_c) .eq. nint(xx_n) .and. &
                    nint(yy_c) .eq. nint(yy_n) .and. &
                    nint(zz_c) .eq. nint(zz_n) 
            IF (bingo) THEN 
               n = jk
               found = .true. 
            END IF
          END IF
      END DO
      !
      IF (n .eq. 0) call errore('createkmap','problem indexing k+q',1)
      !
      kmap( ik ) = n
      !
      !  determine the G_0 such that k+q+G_0 belongs to the first BZ
      !
      g0vec(1) = ( i - nint(xx) )/ nk1
      g0vec(2) = ( j - nint(yy) )/ nk2
      g0vec(3) = ( k - nint(zz) )/ nk3
!      IF ( ik .eq. 3 ) THEN 
!         WRITE(*,'(a,8i6,3f14.7)') 'here ', ik, ng0vec, i, j, k, nint(xx), nint(yy), nint(zz), xk_q(:)
!      ENDIF

      IF (iverbosity.eq.1) THEN
        WRITE(stdout,'(a,i3,a,3i3)') 'ik = ',ik,',   G_0 in integer coord:',g0vec(:)
        WRITE(stdout,'(2i5)') ik, kmap(ik)
      ENDIF
      !
      !  now store the shift for this k+q point
      !
      in_the_list = .false.
      ig0 = 0
      DO WHILE ((ig0.le.ng0vec).and.(.not.in_the_list))
        ig0 = ig0 + 1
        in_the_list = (g0vec(1) .eq. g0vec_all(1,ig0)) .and. &
                      (g0vec(2) .eq. g0vec_all(2,ig0)) .and. &
                      (g0vec(3) .eq. g0vec_all(3,ig0))
      ENDDO
      shift( ik ) = ig0
      !
      IF (.not.in_the_list) CALL errore &
         ('createkmap','cannot find the folding vector in the list',1)
      IF (iverbosity.eq.1) THEN
        WRITE(stdout,'(a,i3,10x,a,i3)') 'ik = ',ik,'shift code:', shift(ik)
        WRITE(stdout,*) 
      ENDIF
      !
      !  obsolete:
      !
      !  very important: now redefine k+q through the corresponding kpoint on the k mesh
      !  Note that this will require using the periodic gauge in the calculation of the
      !  electron-phonon matrix elements (factor e^iG_0*r if G_0 is the vector USEd for 
      !  refolding)
      !
      xkq (:, ik) = xk (:, n) 
      !
    ENDDO
    !
    ! bring everybody back to cartesian coordinates 
    !
    CALL cryst_to_cart ( 1, xq, bg, 1) 
    CALL cryst_to_cart ( nkstot, xk, bg, 1)
    !
    g0vec_all_r = float ( g0vec_all )
    CALL cryst_to_cart ( ng0vec, g0vec_all_r, bg, 1)
    !
    IF (iverbosity.eq.1) THEN
      WRITE(stdout,*) 
      WRITE(stdout,*) '  ik, kmap(ik), shift(ik)'
      DO ik = 1, nkstot
        WRITE(stdout,'(3i4)') ik, kmap(ik), shift(ik)
      ENDDO
      WRITE(stdout,*) '  There are ',ng0vec,'inequivalent folding G_0 vectors'
      WRITE(stdout,*) '  crystal coordinates: '
      DO ig0 = 1,ng0vec
        WRITE(stdout,'(a,i3,a,3i3)') 'g0vec( ',ig0,') = ',g0vec_all(:,ig0)
      ENDDO
      WRITE(stdout,*) '  cartesian coordinates: '
      DO ig0 = 1,ng0vec
        WRITE(stdout,'(3f9.5)') g0vec_all_r(:,ig0)
      ENDDO
    ENDIF
    !
    !  the unit with kmap(ik) and shift(ik)
    ! 
    iukmap  = 97
    !
    OPEN ( unit = iukmap, file = TRIM(prefix)//'.kmap', form = 'formatted')
    DO ik = 1, nkstot
      WRITE (iukmap, '(3i6)') ik, kmap(ik), shift(ik) 
    ENDDO
    CLOSE ( unit = iukmap )
    !
#ifdef __PARA
  ENDIF
  CALL mp_barrier()
#endif
  !
  END SUBROUTINE createkmap


  !-----------------------------------------------------------------------
  SUBROUTINE createkmap2 ( xxq )
  !-----------------------------------------------------------------------
  !
  !  generate the map k+q --> k for folding the rotation matrix U(k+q) 
  !  
  !  in parallel case, this subroutine must be called ONLY by first proc 
  !  (which has all the kpoints)
  !
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  USE klist, ONLY : nkstot, nks, xk
  USE klist_epw, ONLY : kmap  
  USE pwcom, ONLY : nk1, nk2, nk3
  USE epwcom, ONLY : xk_cryst
  USE io_global, ONLY : stdout
  USE control_flags, ONLY : iverbosity
  USE control_ph, ONLY : lgamma
  USE el_phon,   ONLY : xkq
  implicit none
  !
  real(kind=DP) :: eps, xxq (3), xx, yy, zz, xx_c, yy_c, zz_c
  integer :: i, k, n, ik, jk, j
  logical :: in_the_list, bingo, found
  !
  !  the first proc keeps a copy of all kpoints !
  !
  ! bring q in crystal coordinates and check commensuration
  ! 
  !  loosy tolerance: not important since k+q is defined trhough nint() 
  eps = 1.d-4 
  ! 
  CALL cryst_to_cart (1,xxq,at,-1)
  !
  xx = xxq(1)*nk1 
  yy = xxq(2)*nk2 
  zz = xxq(3)*nk3 
  in_the_list = abs(xx-nint(xx)).le.eps .and. &
                abs(yy-nint(yy)).le.eps .and. &
                abs(zz-nint(zz)).le.eps
  IF (iverbosity.eq.1) &
    WRITE(stdout,'(a,3i3)') '  q in integer coord:',nint(xx),nint(yy),nint(zz)
  IF (.not.in_the_list) CALL errore('create_kmap2','q-vec not commensurate',1)
  IF (.not. allocated(xkq) ) ALLOCATE(xkq (3, nkstot) )
  !
  !  bring all the k's in crystal coordinates 
  !
  CALL cryst_to_cart ( nkstot, xk, at, -1)
  !
  DO ik = 1, nkstot
    !
    !  check that the k's are actually on a uniform mesh centered at gamma
    !
    xx = xk(1, ik)*nk1
    yy = xk(2, ik)*nk2
    zz = xk(3, ik)*nk3
    IF (iverbosity.eq.1) &
      WRITE(stdout,'(a,i3,a,3i3)') 'ik = ',ik,',   k   in integer coord:',nint(xx),nint(yy),nint(zz)
    in_the_list = abs(xx-nint(xx)).le.eps .and. &
                  abs(yy-nint(yy)).le.eps .and. &
                  abs(zz-nint(zz)).le.eps
    IF (.not.in_the_list) CALL errore('createkmap2','is this a uniform k-mesh?',1)
    !
    !  now add the phonon wavevector and check that k+q falls again on the k grid
    !
    xkq (:, ik) = xk (:, ik) + xxq (:)
    !
    xx = xkq(1, ik)*nk1
    yy = xkq(2, ik)*nk2
    zz = xkq(3, ik)*nk3
    IF (iverbosity.eq.1) &
      WRITE(stdout,'(a,i3,a,3i3)') 'ik = ',ik,',   k+q in integer coord:',nint(xx),nint(yy),nint(zz)
    in_the_list = abs(xx-nint(xx)).le.eps .and. &
                  abs(yy-nint(yy)).le.eps .and. &
                  abs(zz-nint(zz)).le.eps
    IF (.not.in_the_list) CALL errore('createkmap2','k+q does not fall on k-grid',1)
    !
    !  find the index of this k+q in the k-grid
    !
    IF (iverbosity.eq.1) &
      WRITE(stdout,'(a,i3,a,3i3)') 'ik = ',ik,', f-k+q in integer coord:',i,j,k
    !
    ! make sure xx, yy and zz are in 1st BZ
    !
    CALL backtoBZ( xx, yy, zz, nk1, nk2, nk3 )
    !
    n = 0
    found = .false.
    DO jk = 1, nkstot
       IF ( .not. found ) THEN
          xx_c = xk_cryst(1,jk)*nk1
          yy_c = xk_cryst(2,jk)*nk2
          zz_c = xk_cryst(3,jk)*nk3
          !
          ! check that the k-mesh was defined in the positive region of 1st BZ
          !
          IF ( xx_c .lt. -eps .or. yy_c .lt. -eps .or. zz_c .lt. -eps ) &
             call errore('ktokpmq','coarse k-mesh needs to be strictly positive in 1st BZ',1)
          !
          bingo = nint(xx_c) .eq. nint(xx) .and. &
                  nint(yy_c) .eq. nint(yy) .and. &
                  nint(zz_c) .eq. nint(zz)
          IF (bingo) THEN  
             n = jk
             found = .true.
          END IF
        END IF
    END DO
    !
    IF (n .eq. 0) call errore('createkmap2','problem indexing k+q',1)
    !
    kmap( ik ) = n
    !
  ENDDO
  !
  ! bring everybody back to cartesian coordinates 
  !
  CALL cryst_to_cart ( 1, xxq, bg, 1) 
  CALL cryst_to_cart ( nkstot, xk, bg, 1)
  !
  RETURN
  END SUBROUTINE createkmap2

  !-------------------------------------------------------------------------
  SUBROUTINE createkmap_pw2(xk_all,nkstot,xq0)
  !-------------------------------------------------------------------------
  !
  ! Creates the first instances of [prefix].kmap and [prefix].kgmap. Previously
  ! this was done in PW2 (or set_kplusq, refold, etc. even earlier).
  !
  !-------------------------------------------------------------------------
USE kinds, ONLY        : DP
USE pwcom,    ONLY     : nk1,nk2,nk3,at,bg
USE epwcom,   ONLY     : xk_cryst
USE io_global, ONLY    : stdout
USE io_files, ONLY     : prefix
USE gsmooth, ONLY      : ngms, gcutms, ngms_g,nr1s,nr2s,nr3s,&
                         nrx1s,nrx3s,nls
USE gvect, ONLY        : g,gg,ngm,ngm_g,ngm_l,nr1,nr2,nr3,&
                         gcutm, nrx1,nrx2,nrx3,ig1,ig2,ig3,&
                         nl, gstart, gl, ngl, igtongl
USE reciprocal_vectors, ONLY : ig_l2g
USE control_flags, ONLY: iverbosity,gamma_only
USE constants, ONLY    : eps8
USE fft_base, ONLY     : dfftp,dffts
USE kfold
#ifdef __PARA
  USE mp_global, ONLY : my_pool_id, mpime
  USE mp, ONLY        : mp_barrier
#endif

IMPLICIT NONE

INTEGER, INTENT(IN) :: nkstot
REAL(kind=DP), INTENT(IN), DIMENSION(3,nkstot) :: xk_all
real(kind=DP) :: xkq_all(3,nkstot) 
! Local variables
REAL(kind=DP) :: xq0(3)
INTEGER:: ik, jk, j
REAL(KIND=DP) :: xx, yy, zz, xx_n, yy_n, zz_n, xx_c, yy_c, zz_c, eps
INTEGER :: i,k,n,g0vec(3),ig0,iukmap,iukgmap,ig_1,ig_2,ig_3
LOGICAL :: in_the_list, bingo, found

! variables for modified ggen.f90 section of this routine
INTEGER, ALLOCATABLE :: igsrt(:)
INTEGER, ALLOCATABLE :: mill_g(:,:)
REAL(KIND=DP),ALLOCATABLE :: g2sort_g(:),esort(:)
REAL(KIND=DP) :: t(3), tt, swap
INTEGER :: ngmx,n1,n2,n3
INTEGER :: ipol,ng,iswap,indsw, indnew, indold
! New local variables  
REAL(KIND=DP) :: gg_2(ngm), g_2(3,ngm)
INTEGER :: ngm_2
INTEGER :: nl_2(ngm)

#ifdef __PARA
  INTEGER :: m1,m2,mc
  IF (mpime==0) THEN
#endif 
!
eps = 1.d-5
     !
     iukmap=97   ! unit for the file prefix.kmap
     iukgmap=96  ! unit for the file prefix.kgmap
     !
     !
     WRITE(stdout, '(/5x,a)') 'Calculating kmap and kgmap'
     CALL flush(stdout)
     !
     OPEN(UNIT=iukmap, file = trim(prefix)//'.kmap', form='formatted')
     OPEN(UNIT=iukgmap,file = trim(prefix)//'.kgmap',form='formatted')
     ! 
     IF ( .NOT. ALLOCATED(shift)) ALLOCATE ( shift(nkstot) )

! Fold the k+q back into the k grid for Wannier interpolation. This is 
! meaningful ONLY when q is commensurate.

! bring the kpoints into crystal coordinates 
DO ik = 1, nkstot
   xkq_all(:,ik) = xk_all(:,ik) + xq0(:)
ENDDO
!
CALL cryst_to_cart(nkstot,xk_all,at,-1)
CALL cryst_to_cart(nkstot,xkq_all,at,-1)
! the 5^3 possible G translations
ng0vec=0
DO ig_1=-2,2
   DO ig_2=-2,2
      DO ig_3=-2,2
         ng0vec=ng0vec+1
         g0vec_all(1,ng0vec)=ig_1
         g0vec_all(2,ng0vec)=ig_2
         g0vec_all(3,ng0vec)=ig_3
      END DO
   END DO
END DO
!
! go through all the k+q points and find their index in the k grid
! ONLY the even points are from the result of a k+q operation
DO ik=1,nkstot       
!
           IF (ik.eq.1) THEN
              WRITE(6,'(5x,"Progress  kmap: ")',advance='no')
              indold = 0
           ENDIF
           indnew = nint(float(ik)/float(nkstot)*40)
           IF (indnew.ne.indold) WRITE(6,'(a)',advance='no') '#'
           indold = indnew
           CALL flush(6) ! ONLY on opteron
!
   xx=xkq_all(1,ik)*nk1
   yy=xkq_all(2,ik)*nk2
   zz=xkq_all(3,ik)*nk3
   IF (iverbosity==1) &
        WRITE(stdout,'(a,i3,a,3i3)') 'ik = ', ik, ',   k in integer coord: ', nint(xx), nint(yy),nint(zz)
   i=mod ( nint (xx+ 2*nk1), nk1 )
   j=mod ( nint (yy+ 2*nk2), nk2 )
   k=mod ( nint (zz+ 2*nk3), nk3 )
   xx_n = xx
   yy_n = yy
   zz_n = zz
   !
   !  make sure xx_n, yy_n and zz_n are in 1st BZ
   !
   CALL backtoBZ( xx_n, yy_n, zz_n, nk1, nk2, nk3)
   !
   n = 0
   found = .false.
   DO jk = 1, nkstot
      IF ( .not. found ) THEN
         xx_c = xk_cryst(1,jk)*nk1
         yy_c = xk_cryst(2,jk)*nk2
         zz_c = xk_cryst(3,jk)*nk3
         !
         ! check that the k-mesh was defined in the positive region of 1st BZ
         !
         IF ( xx_c .lt. -eps .or. yy_c .lt. -eps .or. zz_c .lt. -eps ) &
            call errore('ktokpmq','coarse k-mesh needs to be strictly positive in 1st BZ',1)
         !
         bingo = nint(xx_c) .eq. nint(xx_n) .and. &
                 nint(yy_c) .eq. nint(yy_n) .and. &
                 nint(zz_c) .eq. nint(zz_n) 
         IF (bingo) THEN  
            n = jk
            found = .true. 
         END IF
       END IF
   END DO
   WRITE(iukmap,*) ik,n
   ! determine the G_0 such that k+q+G_0 belongs to the first BZ
   g0vec(1) = ( i - nint(xx) )/ nk1
   g0vec(2) = ( j - nint(yy) )/ nk2
   g0vec(3) = ( k - nint(zz) )/ nk3
   in_the_list = .false.
   ig0=0
   DO WHILE ((ig0<ng0vec) .and. (.not. in_the_list))
      ig0=ig0+1
      in_the_list = (g0vec(1) == g0vec_all(1,ig0)) .and. &
                    (g0vec(2) == g0vec_all(2,ig0)) .and. &
                    (g0vec(3) == g0vec_all(3,ig0)) 
   END DO
   shift( ik ) = ig0
END DO

! go back to cartesian coordinates
CALL cryst_to_cart(nkstot,xk_all,bg,1)
CALL cryst_to_cart(nkstot,xkq_all,bg,1)

g0vec_all_r = float ( g0vec_all )
CALL cryst_to_cart ( ng0vec, g0vec_all_r,bg,1)

   

DO ik=1,nkstot
   WRITE(iukgmap, '(3i6)') ik, shift(ik)
END DO
WRITE(iukgmap, '(i5)') ng0vec
DO ig0=1,ng0vec
   WRITE(iukgmap,'(3f20.15)') g0vec_all_r (:,ig0)
END DO

IF (allocated(shift) ) DEALLOCATE(shift)
CLOSE(iukmap)

#ifdef __PARA
END IF  
CALL mp_barrier()
#endif

! below are the routines previously found in the modified ggen.f90 code 



!gg(:) =gcutm +1.d0
nl_2(:)=nl(:)
gg_2(:)=gg(:)
gg_2(:) = gcutm + 1.d0

IF (ALLOCATED(ig_l2g)) DEALLOCATE(ig_l2g)

ALLOCATE( ig_l2g( ngm_l ) )
ALLOCATE( mill_g( 3, ngm_g ) )
ALLOCATE( igsrt( ngm_g ) )
ALLOCATE( g2sort_g( ngm_g ) )
g2sort_g(:) = 1.0d20
n1=nr1+1
n2=nr2+1
n3=nr3+1
ngmx=ngm

ngm_2 = 0
ngms = 0


  DO i = - n1, n1
     !
     ! Gamma-ONLY: exclude space with x < 0
     !
     IF ( gamma_only .AND. i < 0) go to 10
     DO j = - n2, n2
        !
        ! exclude plane with x = 0, y < 0
        !
        IF ( gamma_only .AND. i == 0 .AND. j < 0) go to 11
        DO k = - n3, n3
           !
           ! exclude line with x = 0, y = 0, z < 0
           !
           IF ( gamma_only .AND. i == 0 .AND. j == 0 .AND. k < 0) go to 12
           tt = 0.d0
           DO ipol = 1, 3
              t (ipol) = i * bg (ipol, 1) + j * bg (ipol, 2) + k * bg (ipol, 3)
              tt = tt + t (ipol) * t (ipol)
           ENDDO
           IF (tt <= gcutm) THEN
              ngm_2 = ngm_2 + 1 
              IF (tt <= gcutms) ngms = ngms + 1
              IF (ngm_2 > ngm_g) CALL errore ('ggen', 'too many g-vectors', ngm_2)
              mill_g( 1, ngm_2 ) = i
              mill_g( 2, ngm_2 ) = j
              mill_g( 3, ngm_2 ) = k
             IF ( tt > eps8 ) THEN
                 g2sort_g(ngm_2) = tt
              ELSE
                 g2sort_g(ngm_2) = 0.d0
              ENDIF
           END IF
12         CONTINUE
        ENDDO
11      CONTINUE
     ENDDO
10   CONTINUE
  ENDDO
  IF (ngm_2  /= ngm_g ) &
       CALL errore ('ggen', 'g-vectors missing !', ABS(ngm_2 - ngm_g))
  IF (ngms /= ngms_g) &
       CALL errore ('ggen', 'smooth g-vectors missing !', ABS(ngms - ngms_g))
  igsrt(1) = 0


  CALL hpsort_eps( ngm_g, g2sort_g, igsrt, eps8 )
  DO ng = 1, ngm_g-1
    indsw = ng
7   IF(igsrt(indsw) /= ng) THEN
! ..  swap indices
      DO i = 1, 3
        iswap = mill_g(i,indsw)
        mill_g(i,indsw) = mill_g(i,igsrt(indsw))
        mill_g(i,igsrt(indsw)) = iswap
      END DO
! ..  swap indices
      iswap = indsw; indsw = igsrt(indsw); igsrt(iswap) = iswap
      IF(igsrt(indsw) == ng) THEN
        igsrt(indsw)=indsw
      ELSE
        GOTO 7
      END IF
    END IF
  END DO
  ALLOCATE(esort(ngm_2) )
  esort(:) = 1.0d20
  ngm_2 = 0
  ngms = 0
  DO ng = 1, ngm_g
    i = mill_g(1, ng)
    j = mill_g(2, ng)
    k = mill_g(3, ng)

#ifdef __PARA
    m1 = MOD (i, nr1) + 1
    IF (m1.LT.1) m1 = m1 + nr1
    m2 = MOD (j, nr2) + 1
    IF (m2.LT.1) m2 = m2 + nr2
    mc = m1 + (m2 - 1) * nrx1
    IF ( dfftp%isind ( mc ) .EQ.0) GOTO 1
#endif

    tt = 0.d0
    DO ipol = 1, 3
      t (ipol) = i * bg (ipol, 1) + j * bg (ipol, 2) + k * bg (ipol, 3)
      tt = tt + t (ipol) * t (ipol)
    ENDDO

    ngm_2 = ngm_2 + 1
    IF (tt <= gcutms) ngms = ngms + 1
    IF (ngm_2 > ngmx) CALL errore ('ggen', 'too many g-vectors', ngm_2)
    !
    !  Here map local and global g index !!!
    !
    ig_l2g( ngm_2 ) = ng
    !

    g_2 (1:3, ngm_2) = t (1:3)
    gg_2 (ngm_2) = tt

    IF (tt > eps8) THEN
      esort (ngm_2) = tt 
    ELSE
      esort (ngm_2) = 0.d0
    ENDIF

1   CONTINUE
  ENDDO

     IF (ngm_2.NE.ngmx) &
          CALL errore ('ggen', 'g-vectors missing !', ABS(ngm_2 - ngmx))
     !
     !   reorder the g's in order of increasing magnitude. On exit
     !   from hpsort esort is ordered, and nl contains the new order.
     !
     !   initialize the index inside sorting routine

     nl_2 (1) = 0
     CALL hpsort_eps ( ngm_2, esort, nl_2, eps8 )
     !
     DEALLOCATE( esort  )
     !
     !   reorder also the g vectors, and nl
     !
     DO ng = 1, ngm_2 - 1
20      indsw = nl_2 (ng)
        IF (indsw.NE.ng) THEN
           DO ipol = 1, 3
              swap = g_2 (ipol, indsw)
              g_2 (ipol, indsw) = g_2 (ipol, nl_2 (indsw) )
              g_2 (ipol, nl_2 (indsw) ) = swap
              !
              i=mill_g(ipol,indsw)
              mill_g(ipol,indsw)=mill_g( ipol , nl_2(indsw) )
              mill_g(ipol , nl_2(indsw) )=i
              !
           ENDDO
           swap = gg_2 (indsw)
           gg_2 (indsw) = gg_2 (nl_2 (indsw) )
           gg_2 (nl_2 (indsw) ) = swap
          !
          !  Remember: ig_l2g is the index of a given G vectors in the
          !  sorted global array containing all G vectors, it is used to
          !  collect all wave function components
          !
          iswap = ig_l2g( indsw )
          ig_l2g( indsw ) = ig_l2g( nl_2(indsw) )
          ig_l2g( nl_2(indsw) ) = iswap

           iswap = nl_2 (ng)
           nl_2 (ng) = nl_2 (indsw)
           nl_2 (indsw) = iswap

           GOTO 20
        ENDIF
     ENDDO
     !
     CALL refold( ngm_g, mill_g)
     !

#ifdef __PARA
  CALL mp_barrier
#endif
  DEALLOCATE(ig_l2g,mill_g,igsrt,g2sort_g)
  END SUBROUTINE createkmap_pw2
