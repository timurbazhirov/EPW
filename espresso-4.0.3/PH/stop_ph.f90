!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE stop_ph( flag )
  !----------------------------------------------------------------------------
  !
  ! ... Synchronize processes before stopping.
  !
  USE kinds, ONLY : DP
  USE mp,    ONLY : mp_end, mp_barrier
  !
  USE parallel_include
  !
  IMPLICIT NONE
  !
  INTEGER :: info
  LOGICAL :: flag, exst
  !
  !
  CALL print_clock_ph()
  !
  CALL mp_barrier()
  !
  CALL mp_end()
  !
#if defined (__T3E)
  !
  ! ... set streambuffers off
  !
  CALL set_d_stream( 0 )
  !
#endif
  !
  CALL deallocate_part()
  !
  IF ( flag ) THEN
     !
     STOP
     !
  ELSE
     !
     STOP 1
     !
  ENDIF
  !
END SUBROUTINE stop_ph
