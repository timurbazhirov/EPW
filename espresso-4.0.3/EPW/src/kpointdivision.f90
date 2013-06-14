  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !---------------------------------------------------------------------
  subroutine kpointdivision ( ik0 )
  !---------------------------------------------------------------------
  !
  ! This is just to find the first kpoint block in the pool
  !
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  !
#ifdef __PARA
!  use para
  USE mp_global,   ONLY : my_pool_id,npool
#endif
  use pwcom, only :  nkstot
  use phcom, only : lgamma
  ! 
  implicit none
  integer :: nkl, nkr, ik0, iks
  !
#ifdef __PARA
  !
  !   number of kpoint blocks, kpoints per pool and reminder
  !
  nkl   = 1 * ( nkstot / npool )
  nkr   = ( nkstot - nkl * npool ) / 1
  !
  !  the reminder goes to the first nkr pools  
  !   (0...nkr-1)
  !
  IF ( my_pool_id < nkr ) nkl = nkl + 1 !kunit
  !
  !  the index of the first k point in this pool
  !
  iks = nkl * my_pool_id + 1
  IF ( my_pool_id >= nkr ) iks = iks + nkr * 1 !kunit
  !
  !  the index of the first k point block in this pool - 1
  !  (I will need the index of ik, not ikk)
  !
  ik0 = ( iks - 1 )
  !
#else
  ik0 = 0
#endif
  !
  end subroutine kpointdivision

