!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
SUBROUTINE d3_readin()
  !-----------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program phononq. T
  !    input is read from unit 5. A namelist is used on the machine which
  !    allows it. A second routine readfile reads the variables saved
  !    on the data file by the self-consistent program.
  !
  USE ions_base,     ONLY : nat, ntyp => nsp, amass
  USE pwcom
  USE control_flags, ONLY : iverbosity
  USE phcom
  USE d3com
  USE noncollin_module, ONLY : noncolin
  USE io_files,      ONLY : tmp_dir, prefix, trimcheck
  USE io_global,     ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER :: ios, ipol, iter, na, it, ii
  ! counters
  CHARACTER(len=256) :: outdir

  NAMELIST / inputph / ethr_ph, amass, iverbosity, outdir, prefix, &
       fildyn, fildrho, fild0rho, q0mode_todo, wraux, recv, istop, &
       testflag, testint, testreal
  ! convergence threshold
  ! atomic masses
  ! write control
  ! directory for temporary files
  ! the punch file produced by pwscf
  ! the file with the dynamical matrix
  ! the file with the deltarho
  ! the file with q=0 deltarho
  ! list of the q=0 modes to be computed
  ! .true.==> writes some auxiliary
  ! .true.==> this is a recover run
  ! to stop the program at a given point
  ! variables used for testing purposes

  IF ( ionode ) THEN
     !
     !    Read the first line of the input file
     !
     READ (5, '(a)', err = 100, iostat = ios) title
100  CALL errore ('d3_readin', 'reading title ', ABS (ios) )
     !
     !   set default values for variables in namelist
     !
     ethr_ph = 1.d-5
     iverbosity = 0
     CALL get_env( 'ESPRESSO_TMPDIR', outdir )
     IF ( TRIM( outdir ) == ' ' ) outdir = './'
     prefix = 'pwscf'
     fildyn = 'd3dyn'
     fildrho = ' '
     fild0rho = ' '
     DO ii = 1, 300
        q0mode_todo (ii) = 0
     ENDDO
     wraux = .FALSE.
     recv = .FALSE.
     istop = 0
     DO ii = 1, 50
        testflag (ii) = .FALSE.
     ENDDO
     !
     !     reading the namelist inputph
     !
     !   Note: for AIX machine (xlf compiler version 3.0 or higher):
     !   The variable XLFRTEOPTS must be set to "namelist=old"
     !   in order to have "&end" to end the namelist
     !
     READ (5, inputph, err = 200, iostat = ios)

200  CALL errore ('d3_readin', 'reading inputph namelist', ABS (ios) )
     !
     !     Check all namelist variables
     !
     IF (ethr_ph.LE.0.d0) CALL errore (' d3_readin', ' Wrong ethr_ph ', 1)
     IF (iverbosity.NE.0.AND.iverbosity.NE.1) &
          CALL errore ('d3_readin', ' Wrong iverbosity ', 1)
     IF (fildrho.EQ.' ') CALL errore ('d3_readin', ' Wrong fildrho ', 1)
     IF (fild0rho.EQ.' ') CALL errore ('d3_readin', ' Wrong fild0rho ', 1)
     !
     !    reads the q point
     !
     READ (5, *, err = 300, iostat = ios) (xq (ipol), ipol = 1, 3)
300  CALL errore ('d3_readin', 'reading xq', ABS (ios) )

     lgamma = xq (1) .EQ.0.d0.AND.xq (2) .EQ.0.d0.AND.xq (3) .EQ.0.d0
     tmp_dir = trimcheck (outdir)
     !
  END IF
  !
  CALL bcast_d3_input()
  !
  xqq (:) = xq(:)
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  !
  IF (lgamma) THEN
     nksq = nks
  ELSE
     nksq = nks / 2
  ENDIF
  !
  IF (lsda) CALL errore ('d3_readin', 'lsda not implemented', 1)
  IF (okvan) CALL errore ('d3_readin', 'US not implemented', 1)
  IF (noncolin) call errore('d3_readin', &
     'd3 is not working in the noncolinear case', 1)

  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  CALL allocate_part

  !IF (.NOT.(lphonon.OR.lgamma)) &
  !   CALL errore ('d3_readin', ' Wrong calculation ', 1 )
  DO it = 1, ntyp
     IF (amass (it) .LE.0.d0) CALL errore ('d3_readin', 'Wrong masses', &
          it)
  ENDDO
  IF (MOD (nks, 2) .NE.0.AND..NOT.lgamma) CALL errore ('d3_readin', &
       'k-points are odd', nks)
  !
  ! q0mode, and q0mode_todo are not allocated dynamically. Their
  ! dimension is fixed to 300
  !

  IF (3 * nat.GT.300) CALL errore ('d3_readin', 'wrong dimension of &
       &q0mode variable', 1)
  DO ii = 1, 3 * nat
     IF (q0mode_todo (ii) .GT.3 * nat) CALL errore ('d3_readin', ' wrong &
          & q0mode_todo ', 1)
  ENDDO
  RETURN
END SUBROUTINE d3_readin
