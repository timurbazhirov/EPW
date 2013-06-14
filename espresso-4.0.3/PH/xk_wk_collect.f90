!
! Copyright (C) 2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE xk_wk_collect( xk_start, wk_start, xk, wk, nkstot, nks )
  !----------------------------------------------------------------------------
  !
  ! ... This routine collects the k points (with granularity kunit) among nodes
  ! ... and sets the variable xk_start and wk_start with the total number of 
  !     k-points
  !
  USE io_global, only : stdout
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : my_pool_id, npool, kunit
  USE mp_global, ONLY : inter_pool_comm, intra_pool_comm
  USE mp,        ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: nkstot, nks
    ! total number of k-points
    ! number of k-points per pool
  REAL (DP) :: xk(3,nks), wk(nks)
  REAL (DP) :: xk_start(3,nkstot), wk_start(nkstot)
    ! k-points
    ! k-point weights
  !
#if defined (__PARA)
  !
  INTEGER :: ik, nbase, rest, nks1
  !
  xk_start=0.d0
  !
  wk_start=0.d0
  !
  nks1    = kunit * ( nkstot / kunit / npool )
  !
  rest = ( nkstot - nks1 * npool ) / kunit
  !
  IF ( ( my_pool_id + 1 ) <= rest ) nks1 = nks1 + kunit
  !
  IF (nks1.ne.nks) &
     call errore('xk_wk_collect','problems with nks1',1)
  !
  ! ... calculates nbase = the position in the list of the first point that
  ! ...                    belong to this npool - 1
  !
  nbase = nks * my_pool_id
  !
  IF ( ( my_pool_id + 1 ) > rest ) nbase = nbase + rest * kunit
  !
  ! copy the original points in the correct position of the list
  !
  xk_start(:,nbase+1:nbase+nks) = xk(:,1:nks)
  !
  wk_start(nbase+1:nbase+nks)=wk(1:nks)
  !
  CALL mp_sum( xk_start, inter_pool_comm )
  !
  CALL mp_sum( wk_start, inter_pool_comm )
  !
#endif
  !
  RETURN
  !
END SUBROUTINE xk_wk_collect
