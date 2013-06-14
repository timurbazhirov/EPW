!
! modified from stop_ph
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine stop_epw (flag)
  !-----------------------------------------------------------------------
  !
  ! Close all files and synchronize processes before stopping.
  ! Called at the end of the run with flag=.true. (removes 'recover')
  ! or during execution with flag=.false. (does not remove 'recover')
  !
  use pwcom
  USE kinds, only : DP
  use phcom
  use epwcom
#ifdef __PARA
  use mp, only: mp_end, mp_barrier
  USE parallel_include
#endif
  implicit none
  logical :: flag
  !
  CALL print_clock_epw
#ifdef __PARA
  CALL mp_barrier()
  CALL mp_end()
#endif


#ifdef __T3E
  !
  ! set streambuffers off
  !
  CALL set_d_stream (0)
#endif

  CALL deallocate_part

  IF (flag) then
     stop
  ELSE
     stop 1
  ENDIF

end subroutine stop_epw
