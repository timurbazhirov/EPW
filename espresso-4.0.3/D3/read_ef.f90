!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------
SUBROUTINE read_ef()
  !-----------------------------------------------------------------------
  !
  ! Reads the shift of the Fermi Energy
  !
  USE pwcom
  USE d3com
  USE io_global, ONLY : ionode, ionode_id
  USE mp,        ONLY : mp_bcast
  !
  IMPLICIT NONE
  !
  INTEGER :: ios
  !
  IF (degauss == 0.d0 ) RETURN
  !
  IF ( ionode ) THEN
     !
     REWIND (unit = iuef)
     READ (iuef, err = 100, iostat = ios) ef_sh
     !
100  CALL errore ('d3_valence', 'reading iuef', ABS (ios) )
     !
  END IF

  CALL mp_bcast( ef_sh, ionode_id )

  RETURN
END SUBROUTINE read_ef
