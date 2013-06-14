  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from PH/ph.f90  
  !-----------------------------------------------------------------------
  PROGRAM epw
  !-----------------------------------------------------------------------
  !
  ! This is the main EPW driver which sets the phases on the wavefunctions,
  ! calls wannierize and elphon_shuffle_wrap
  !
  ! 8/13/08 removed epsil variables, fildyn, ldisp
  !
  ! 8/14/08 lnscf is unnecessary, as is nqs,iq_start
  ! 8/15/08 recover has been cut
  ! 02/2009 in state of v0.2  Barely resembles phonon.f90
  !         
  USE kinds,           ONLY : DP
  USE io_global,       ONLY : stdout
  USE io_files,        ONLY : nd_nmbr
#ifdef __PARA
  USE mp,              ONLY : mp_bcast, mp_barrier
  USE mp_global,       ONLY : my_pool_id
#endif
  USE uspp_param,      ONLY : upf
  USE control_flags,   ONLY : gamma_only
  USE control_ph,      ONLY : elph, maxirr
  USE control_epw,     ONLY : phinterp, epstrict, tshuffle2,wannierize
  USE global_version,  ONLY : version_number
  USE epwcom,          ONLY : filukk
  !
  !
  implicit none
  !
  CHARACTER (LEN=12)   :: code = 'EPW'
  !
  version_number = '2.3.5'
  !
  !
  CALL init_clocks( .TRUE. )
  !
  CALL start_clock( 'EPW' )
  !
  gamma_only = .FALSE.
  !
  CALL startup( nd_nmbr, code, version_number )
  !
  WRITE( stdout, '(/5x,"Ultrasoft (Vanderbilt) Pseudopotentials")' )
  !
  ! read in the input file
  !
  CALL epw_readin
  !
  maxirr = 0
  !
  CALL allocate_epwq
  CALL epw_setup
  !
  !  Print run info to stdout
  !
  CALL epw_summary
  !
  CALL openfilepw
  CALL print_clock( 'EPW' )
  !
  CALL epw_init(.true.)
  CALL print_clock( 'EPW' )
  !
  CALL print_clock( 'EPW' )
  !
  !  Generates the perturbation matrix which fixes the gauge of 
  !  the calculated wavefunctions
  !
  CALL setphases_wrap
  !
  !  Create U(k, k') localization matrix 
  !
  IF (wannierize) THEN
     CALL wann_run
  ELSE
     !
     ! Read Wannier matrix from a previous run
     !
     WRITE(stdout,'(/,5x,a,/,3a,/,5x,a,/)') repeat('-',67), '     Using ', &
          trim(filukk) , ' from disk', repeat('-',67) 
  ENDIF
  !
  IF ( elph ) THEN
     !
     CALL dvanqq2()
     !
     CALL elphon_shuffle_wrap()
     !
     IF ( .NOT. tshuffle2 .AND. phinterp ) THEN
        !
        IF (epstrict) THEN
           !
           CALL elphsum3_strict()
           !
        ELSE
           !
           CALL elphsum3()
           !
        ENDIF
        !
     END IF
     !
  END IF
  !
  ! ... cleanup of the variables
  !
  CALL clean_pw( .FALSE. )
  CALL deallocate_epw
  !
  ! ... Close the files
  !
  CALL close_epw()
  !
  ! ... Print statistics and exit gracefully    
  !
  CALL stop_epw( .TRUE. )
  !
  STOP
  !
  END PROGRAM epw
