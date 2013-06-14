  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
! The mesh loading needs to be cleaner.  The above code is a nightmare
! Each load should be succinct and have the option of
! 1) read from file
! 2) generate on a grid (either uniform or in the wedge)
! 3) generate a random number of points
!
! The total end result of the mesh loading should be
!  xkf, wkf, xqf, wqf, nxqf, nkstotf

!-----------------------------------------------------------------------
SUBROUTINE loadkmesh_para
!-----------------------------------------------------------------------
!
!  load fine k mesh and distribute among pools
!
!-----------------------------------------------------------------------
#include "f_defs.h"
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : mpime, inter_pool_comm, &
       my_pool_id, npool
  USE mp,        ONLY : mp_bcast, mp_sum
#endif
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : filkf, nkf1, nkf2, nkf3, skf1, skf2, skf3, &
       rand_k, rand_nk, mp_mesh_k
  USE el_phon,   ONLY : nkstotf, nksf, xkf, wkf, &
       nksqf
  USE pwcom,     ONLY : at, bg, ibrav, s, sname, symm_type
  USE symme,     ONLY : t_rev, time_reversal
  implicit none
  !
  logical :: exst
  integer, PARAMETER :: iunkkf  = 70         
  ! unit with fine k point mesh (crys coord)
  real(kind=DP), ALLOCATABLE :: xkf_(:,:), wkf_(:), xkf_tmp(:,:), &
       wkf_tmp(:)
  integer :: ik, ikk, ikq, lower_bnd, upper_bnd, rest, &
       ipol, i, j, k, nrot, tipo

  

#ifdef __PARA
 IF (mpime .eq. ionode_id) THEN
#endif
    INQUIRE(file = filkf, exist=exst)
    IF (exst) THEN ! load from file (crystal coordinates)
       !
       WRITE (stdout, *) '     Using k-mesh file: ', trim(filkf)
       OPEN ( unit = iunkkf, file = filkf, status = 'old', form = 'formatted')
       READ(iunkkf, *) nkstotf 
       !
       ALLOCATE (xkf_ (3, 2*nkstotf), wkf_(2*nkstotf))
       !
       DO ik = 1, nkstotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          READ (iunkkf, *) xkf_ (:, ikk ), wkf_ (ikk)
          !
          !  bring the k point to crystal coordinates
          ! CALL cryst_to_cart ( 1, xkf_ (:,ikk), at, -1)
          !
          xkf_ (:, ikq) = xkf_ (:, ikk) 
          wkf_ ( ikq ) = 0.d0
          !
       ENDDO
       CLOSE(iunkkf)
       !
       ! redefine nkstotf to include the k+q points
       !
       nkstotf = 2 * nkstotf
       !
    ELSEIF ( (nkf1.ne.0) .and. (nkf2.ne.0) .and. (nkf3.ne.0) ) THEN ! generate grid
       IF (mp_mesh_k) THEN
           ! get size of the mp_mesh in the irr wedge 
          WRITE (stdout, '(a,3i4)') '     Using uniform MP k-mesh: ', nkf1, nkf2, nkf3
          IF (ibrav == 4 .or. ibrav == 5) THEN
             CALL hexsym (at, s, sname, nrot)
             tipo = 2
          ELSEIF (ibrav >= 1 .and. ibrav <= 14) THEN
             CALL cubicsym (at, s, sname, nrot)
             tipo = 1
          ELSEIF (ibrav == 0) THEN
             IF (symm_type == 'cubic') THEN
                CALL cubicsym (at, s, sname, nrot)
                tipo = 1
             ENDIF
             IF (symm_type == 'hexagonal') THEN
                CALL hexsym (at, s, sname, nrot)
                tipo = 2
             ENDIF
          ELSE
             CALL errore ('loadmesh', 'wrong ibrav', 1)
          ENDIF
          !                                         
          !
          ALLOCATE ( xkf_ (3, 2*nkf1*nkf2*nkf3), wkf_(2*nkf1*nkf2*nkf3) )
          ! the result of this call is just nkstotf
          CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nkf1*nkf2*nkf3, &
               0,0,0, nkf1,nkf2,nkf3, nkstotf, xkf_, wkf_)
          DEALLOCATE (xkf_, wkf_)
          ALLOCATE ( xkf_ (3, 2*nkstotf), wkf_(2*nkstotf)) 
          ALLOCATE (xkf_tmp (3,nkstotf), wkf_tmp(nkstotf))
          CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nkf1*nkf2*nkf3, &
               0,0,0, nkf1,nkf2,nkf3, nkstotf, xkf_tmp, wkf_tmp)
          !  
          ! assign to k and k+q for xkf and wkf 
          ! 
          DO ik = 1, nkstotf
             ikk = 2 * ik - 1
             ikq = ikk + 1
             xkf_(:,ikk) = xkf_tmp(:,ik)
             xkf_(:,ikq) = xkf_tmp(:,ik)
             wkf_(ikk)   = 2.d0 * wkf_tmp(ik)
             wkf_(ikq)   = 0.d0
          ENDDO
          DEALLOCATE (xkf_tmp, wkf_tmp)
          !       
          ! bring the k point to crystal coordinates       
          CALL cryst_to_cart (2*nkstotf, xkf_, at, -1)
          !
          nkstotf = 2 * nkstotf
          !
       ELSE
          !
          IF ( (skf1.gt.1) .and. (skf2.gt.1) .and. (skf3.gt.1) ) THEN
             WRITE (stdout, '(a,3i4)') '     Using uniform shifted k-mesh: ', nkf1, nkf2, nkf3
          ELSE
             WRITE (stdout, '(a,3i4)') '     Using uniform k-mesh: ', nkf1, nkf2, nkf3
          ENDIF
          !
          nkstotf = 2 * nkf1 * nkf2 * nkf3
          ALLOCATE ( xkf_ (3, nkstotf), wkf_(nkstotf) )
          wkf_(:) = 0.d0
          DO ik = 1, nkf1 * nkf2 * nkf3
             wkf_(2*ik-1) = 2.d0/(float(nkstotf/2))
          ENDDO
          DO i = 1, nkf1
             DO j = 1, nkf2
                DO k = 1, nkf3
                   ik = (i-1)*nkf2*nkf3 + (j-1)*nkf3 + k
                   ikk = 2 * ik - 1
                   ikq = ikk + 1
                   xkf_(1, ikk) = float(i-1)/float(nkf1)
                   xkf_(2, ikk) = float(j-1)/float(nkf2)
                   xkf_(3, ikk) = float(k-1)/float(nkf3)
                   IF ( skf1 .gt. 1) xkf_(1, ikk) = xkf_(1, ikk) + 1.d0/float(nkf1)/float(skf1)
                   IF ( skf2 .gt. 1) xkf_(2, ikk) = xkf_(2, ikk) + 1.d0/float(nkf2)/float(skf2)
                   IF ( skf3 .gt. 1) xkf_(3, ikk) = xkf_(3, ikk) + 1.d0/float(nkf3)/float(skf3)
                   xkf_(1, ikq) = xkf_(1, ikk)
                   xkf_(2, ikq) = xkf_(2, ikk)
                   xkf_(3, ikq) = xkf_(3, ikk) 
                ENDDO
             ENDDO
          ENDDO
          !WRITE(stdout,'(a)') '  '
          DO ik = 1, nkf1 * nkf2 * nkf3
             ikk = 2 * ik - 1
             !WRITE(stdout,'(4f12.7)') xkf_(:,ikk), wkf_(ikk)
          ENDDO
          !WRITE(stdout,'(a)') '  '
          !
       ENDIF
       !
    ELSEIF (rand_k) THEN  ! random points
       ! random grid
       WRITE (stdout, *) '    Using random k-mesh: ', rand_nk
       !
       nkstotf = rand_nk
       ALLOCATE (xkf_ (3, 2*nkstotf), wkf_(2*nkstotf))
       !WRITE(stdout,'(a)') '  '
       DO ik = 1, nkstotf
          !
          ikk = 2 * ik - 1
          ikq = ikk + 1
          !
          wkf_(ikk) = 2.d0/ float(nkstotf)
          wkf_(ikq) = 0.d0
          CALL random_number(xkf_(:,ikk))
          xkf_(:,ikq) = xkf_(:,ikk)
          !WRITE(stdout,'(4f12.7)') xkf_(:,ikk), wkf_(ikk)
       ENDDO
       !WRITE(stdout,'(a)') '  '
       !
       ! redefine nkstotf to include the k+q points
       !
       nkstotf = 2 * nkstotf
       !
    ELSE ! don't know how to get grid
       CALL errore('loadkmesh_para', "Cannot load fine k points", 1)
    ENDIF
#ifdef __PARA
 ENDIF
 CALL mp_bcast (nkstotf, ionode_id, inter_pool_comm)
 !
 !  scatter the k points of the fine mesh across the pools
 !
 nksf = 2 * ( nkstotf / 2 / npool )
 rest = ( nkstotf - nksf * npool ) / 2
 IF (my_pool_id < rest ) THEN
    nksf=nksf+2
    lower_bnd = my_pool_id*nksf + 1
    upper_bnd = lower_bnd + nksf - 1
 ELSE
    lower_bnd = rest*(nksf+2)+(my_pool_id-rest)*nksf + 1
    upper_bnd = lower_bnd + nksf - 1
 ENDIF
 !
 nksqf = nksf / 2 
 IF (.not.ALLOCATEd(xkf_)) ALLOCATE (xkf_(3,nkstotf))
 IF (.not.ALLOCATEd(wkf_)) ALLOCATE (wkf_(  nkstotf))
 CALL mp_bcast(xkf_, ionode_id)
 CALL mp_bcast(wkf_, ionode_id)
 !
#else
 !
 ! In serial the definitions are much easier 
 !
 nksf = nkstotf
 nksqf = nksf / 2 
 lower_bnd = 1
 upper_bnd = nksf
 !
#endif
 !
 !  Assign the weights and vectors to the correct bounds
 !
 ALLOCATE(xkf(3,nksf))
 ALLOCATE(wkf(  nksf))
 xkf(:,:) = xkf_ (:, lower_bnd:upper_bnd)
 wkf(  :) = wkf_ (   lower_bnd:upper_bnd)
 !
 WRITE( stdout, '(5x,"Size of k point mesh for interpolation: ",i10)' ) nkstotf 
 WRITE( stdout, '(5x,"Max number of k points per pool:",7x,i10)' ) nksf 
 !
 IF (ALLOCATED(xkf_)) DEALLOCATE(xkf_)
 IF (ALLOCATED(wkf_)) DEALLOCATE(wkf_)
 !
END SUBROUTINE loadkmesh_para
!-----------------------------------------------------------------------

!SUBROUTINE loadkmesh_serial
!end SUBROUTINE loadkmesh_serial



!-----------------------------------------------------------------------
SUBROUTINE loadqmesh_para
!-----------------------------------------------------------------------
!
!  load fine q mesh and distribute among pools
!
!-----------------------------------------------------------------------
#include "f_defs.h"
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : mpime
#endif
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : filqf, nqf1, nqf2, nqf3, rand_q
  USE el_phon,   ONLY : nxqf
  implicit none
  !
  logical :: exst
  integer, PARAMETER :: iunqf  = 77         
  ! unit with fine q point mesh (crys coord)

#ifdef __PARA
 IF (mpime .eq. ionode_id) THEN
#endif

 IF (filqf .ne. '') THEN ! load from file (crystal coordinates)
    inquire(file = filqf, exist=exst)
    IF (.not. exst) CALL errore('loadmesh', 'filqf does not exist',1)
    !
    WRITE (stdout, *) '     Using q-mesh file: ', trim(filqf)
    OPEN ( unit = iunqf, file = filqf, status = 'old', form = 'formatted')
    READ(iunqf, *) nxqf 
    !
 ELSEIF ( (nqf1.ne.0) .and. (nqf2.ne.0) .and. (nqf3.ne.0) ) THEN ! generate grid

 ELSEIF (rand_q) THEN  ! random points
    
 ELSE ! don't know how to get grid
    CALL errore('loadqmesh_para', "Cannot load fine q points", 1)
 ENDIF

#ifdef __PARA
 ENDIF
#endif
END SUBROUTINE loadqmesh_para
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
SUBROUTINE loadqmesh_serial
!-----------------------------------------------------------------------
!
!  load fine q mesh on each pool
!
!-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  USE epwcom,    ONLY : filqf, nqf1, nqf2, nqf3, sqf1, sqf2, sqf3, &
       rand_q, rand_nq, mp_mesh_q
  USE el_phon,   ONLY : xqf, wqf, nksqf, nxqf
  USE pwcom,     ONLY : at, bg, ibrav, s, sname, symm_type
  USE symme,     ONLY : t_rev, time_reversal
  implicit none
  !
  logical :: exst
  integer, PARAMETER :: iunqf  = 77         
  ! unit with fine k point mesh (crys coord)
  integer :: iq , i, j, k, nrot, tipo


  inquire(file = filqf, exist=exst)
  IF (exst) THEN ! load from file
     !
     ! Each pool gets its own copy from the action=read statement
     !
     WRITE (stdout, *) '     Using q-mesh file: ', trim(filqf)
     OPEN ( unit = iunqf, file = filqf, status = 'old', form = 'formatted', action='read')
     READ(iunqf, *) nxqf
     ALLOCATE (xqf(3, nxqf), wqf(nxqf))
     DO iq = 1, nxqf
        READ (iunqf, *) xqf (:, iq), wqf(iq)
     ENDDO
     CLOSE(iunqf)
     !
     ! bring xqf in crystal coordinates
     ! CALL cryst_to_cart (nxqf, xqf, at, -1)
     !
 ELSEIF ( (nqf1.ne.0) .and. (nqf2.ne.0) .and. (nqf3.ne.0) ) THEN ! generate grid
    IF (mp_mesh_q) THEN
       ! get size of the mp_mesh in the irr wedge 
       WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
       IF (ibrav == 4 .or. ibrav == 5) THEN
          CALL hexsym (at, s, sname, nrot)
          tipo = 2
       ELSEIF (ibrav >= 1 .and. ibrav <= 14) THEN
          CALL cubicsym (at, s, sname, nrot)
          tipo = 1
       ELSEIF (ibrav == 0) THEN
          IF (symm_type == 'cubic') THEN
             CALL cubicsym (at, s, sname, nrot)
             tipo = 1
          ENDIF
          IF (symm_type == 'hexagonal') THEN
             CALL hexsym (at, s, sname, nrot)
             tipo = 2
          ENDIF
       ELSE
          CALL errore ('loadmesh', 'wrong ibrav', 1)
       ENDIF
       !                                         
       !
       ALLOCATE ( xqf (3, nqf1*nqf2*nqf3), wqf(nqf1*nqf2*nqf3) )
       ! the result of this call is just nkstotf
       CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nqf1*nqf2*nqf3, &
            0,0,0, nqf1,nqf2,nqf3, nxqf, xqf, wqf)
       DEALLOCATE ( xqf, wqf) 
       ALLOCATE ( xqf(3, nxqf), wqf(nxqf)) 
       CALL kpoint_grid ( nrot, time_reversal, s, t_rev, bg, nqf1*nqf2*nqf3, &
            0,0,0, nqf1,nqf2,nqf3, nxqf, xqf, wqf)
       !
       ! bring xqf in crystal coordinates       
       CALL cryst_to_cart (nxqf, xqf, at, -1)
       !
    ELSE
       ! currently no offset.  
       ! q's are in crystal coordinates in xqf
       IF ( (sqf1.gt.1) .and. (sqf2.gt.1) .and. (sqf3.gt.1) ) THEN
          WRITE (stdout, '(a,3i4)') '     Using uniform shifted q-mesh: ', nqf1, nqf2, nqf3
       ELSE
          WRITE (stdout, '(a,3i4)') '     Using uniform q-mesh: ', nqf1, nqf2, nqf3
       ENDIF
       !
       nxqf = nqf1 * nqf2 * nqf3
       ALLOCATE ( xqf (3, nxqf), wqf(nxqf) )
       wqf(:) = 1.d0/(float(nxqf))
       DO i = 1, nqf1
          DO j = 1, nqf2
             DO k = 1, nqf3
                iq = (i-1)*nqf2*nqf3 + (j-1)*nqf3 + k
                xqf(1, iq) = float(i-1)/float(nqf1)
                xqf(2, iq) = float(j-1)/float(nqf2)
                xqf(3, iq) = float(k-1)/float(nqf3)
                IF ( sqf1 .gt. 1 ) xqf(1, iq) = xqf(1, iq) + 1.d0/float(nqf1)/float(sqf1)
                IF ( sqf2 .gt. 1 ) xqf(2, iq) = xqf(2, iq) + 1.d0/float(nqf2)/float(sqf2)
                IF ( sqf3 .gt. 1 ) xqf(3, iq) = xqf(3, iq) + 1.d0/float(nqf3)/float(sqf3)
             ENDDO
          ENDDO
       ENDDO
       !
       !WRITE(stdout,'(a)') '  '
       DO iq = 1, nxqf
          !WRITE(stdout,'(4f12.7)') xqf(:,iq), wqf(iq)
       ENDDO
       !WRITE(stdout,'(a)') '  '
       !
    ENDIF
 ELSEIF (rand_q) THEN  ! random points
    ! random grid
    WRITE (stdout, *) '    Using random q-mesh: ', rand_nq
    !
    nxqf = rand_nq
    ALLOCATE (xqf(3, nxqf), wqf(nxqf))
    !WRITE(stdout,'(a)') '  '
    wqf(:) = 1.d0/(float(nxqf))
    DO iq = 1, nxqf
       !
       CALL random_number(xqf(:,iq))
       !WRITE(stdout,'(4f12.7)') xqf(:,iq), wqf(iq)
    ENDDO
    !WRITE(stdout,'(a)') '  '
    !
 ELSE ! don't know how to get grid
    CALL errore('loadqmesh_serial', "Cannot load fine q points", 1)
 ENDIF
 !
 IF (abs(sum (wqf) - 1.d0) .gt. 1.d-4 ) &
      CALL errore('loadqmesh_serial',"q-point wieghts do not add up to 1", 1)
 !
 WRITE( stdout, '(5x,"Size of q point mesh for interpolation: ",i10)' ) nxqf
 !
end SUBROUTINE loadqmesh_serial
!-----------------------------------------------------------------------
