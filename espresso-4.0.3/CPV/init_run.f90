!
! Copyright (C) 2002-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE init_run()
  !----------------------------------------------------------------------------
  !
  ! ... this routine initialise the cp code and allocates (calling the
  ! ... appropriate routines) the memory
  !
  USE kinds,                    ONLY : DP
  USE control_flags,            ONLY : nbeg, nomore, lwf, iprsta, &
                                       ndr, tfor, tprnfor, tpre, program_name, &
                                       force_pairing, newnfi, tnewnfi, ndw
  USE cp_electronic_mass,       ONLY : emass, emass_cutoff
  USE ions_base,                ONLY : na, nax, nat, nsp, iforce, pmass, ityp, cdms
  USE ions_positions,           ONLY : tau0, taum, taup, taus, tausm, tausp, &
                                       vels, velsm, velsp, fion, fionm,      &
                                       atoms0, atomsm, atomsp
  USE gvecw,                    ONLY : ngw, ecutw, ngwt, ggp
  USE gvecb,                    ONLY : ngb
  USE gvecs,                    ONLY : ngs
  USE gvecp,                    ONLY : ngm
  USE reciprocal_vectors,       ONLY : gzero
  USE grid_dimensions,          ONLY : nnrx, nr1, nr2, nr3
  USE fft_base,                 ONLY : dfftp
  USE electrons_base,           ONLY : nspin, nbsp, nbspx, nupdwn, f
  USE electrons_module,         ONLY : n_emp
  USE uspp,                     ONLY : nkb, vkb, deeq, becsum,nkbus
  USE core,                     ONLY : rhoc
  USE smooth_grid_dimensions,   ONLY : nnrsx
  USE wavefunctions_module,     ONLY : c0, cm, cp
  USE cdvan,                    ONLY : dbec, drhovan
  USE derho,                    ONLY : drhor, drhog
  USE ensemble_dft,             ONLY : tens, z0t
  USE cg_module,                ONLY : tcg
  USE electrons_base,           ONLY : nudx, nbnd
  USE parameters,               ONLY : nspinx
  USE efield_module,            ONLY : tefield, tefield2
  USE uspp_param,               ONLY : nhm
  USE ions_nose,                ONLY : xnhp0, xnhpm, vnhp, nhpcl, nhpdim
  USE cell_base,                ONLY : h, hold, hnew, velh, tpiba2, ibrav, &
                                       alat, celldm, a1, a2, a3, b1, b2, b3
  USE cp_main_variables,        ONLY : lambda, lambdam, lambdap, ema0bg, bec,  &
                                       becdr, sfac, eigr, ei1, ei2, ei3, taub, &
                                       irb, eigrb, rhog, rhos, rhor, bephi,    &
                                       becp, acc, acc_this_run, wfill, &
                                       edft, nfi, vpot, ht0, htm
  USE cp_main_variables,        ONLY : allocate_mainvar, nlax, descla, nrlx
  USE energies,                 ONLY : eself, enl, ekin, etot, enthal, ekincm
  USE dener,                    ONLY : detot
  USE time_step,                ONLY : dt2, delt, tps
  USE electrons_nose,           ONLY : xnhe0, xnhem, vnhe
  USE cell_nose,                ONLY : xnhh0, xnhhm, vnhh
  USE gvecp,                    ONLY : ecutp
  USE funct,                    ONLY : dft_is_meta
  USE metagga,                  ONLY : crosstaus, dkedtaus, gradwfc
  USE pseudopotential,          ONLY : nsanl
  !
  USE efcalc,                   ONLY : clear_nbeg
  USE local_pseudo,             ONLY : allocate_local_pseudo
  USE cp_electronic_mass,       ONLY : emass_precond
  USE wannier_subroutines,      ONLY : wannier_startup
  USE from_scratch_module,      ONLY : from_scratch
  USE from_restart_module,      ONLY : from_restart  
  USE cp_interfaces,            ONLY : readfile
  USE ions_base,                ONLY : ions_cofmass
  USE ensemble_dft,             ONLY : id_matrix_init, allocate_ensemble_dft, h_matrix_init
  USE efield_module,            ONLY : allocate_efield, allocate_efield2
  USE cg_module,                ONLY : allocate_cg
  USE wannier_module,           ONLY : allocate_wannier  
  USE io_files,                 ONLY : outdir, prefix
  USE io_global,                ONLY : ionode, stdout
  USE printout_base,            ONLY : printout_base_init
  USE wave_types,               ONLY : wave_descriptor_info
  USE xml_io_base,              ONLY : restart_dir, create_directory
  USE orthogonalize_base,       ONLY : mesure_diag_perf, mesure_mmul_perf
  USE step_constraint,          ONLY : step_con
  USE ions_base,                ONLY : ions_reference_positions, cdmi, taui
  USE ldau
  !
  IMPLICIT NONE
  !
  INTEGER            :: neupdwn( nspinx )
  INTEGER            :: lds_wfc, i
  CHARACTER(LEN=256) :: dirname
  !
  !
  CALL start_clock( 'initialize' )
  !
  ! ... initialize directories
  !
  CALL printout_base_init( outdir, prefix )
  !
  dirname = restart_dir( outdir, ndw )
  !
  ! ... Create main restart directory
  !
  CALL create_directory( dirname )
  !
  ! ... initialize g-vectors, fft grids 
  ! ... The number of g-vectors are based on the input celldm!
  !
  CALL init_dimensions()
  !
  ! ... initialize atomic positions and cell
  !
  CALL init_geometry()
  !
  ! ... mesure performances of parallel routines
  !
  CALL mesure_mmul_perf( nudx )
  !
  CALL mesure_diag_perf( nudx )
  !
  IF ( lwf ) CALL clear_nbeg( nbeg )
  !
  !=======================================================================
  !     allocate and initialize nonlocal potentials
  !=======================================================================
  !
  CALL nlinit()
  !
  !=======================================================================
  !     allocation of all arrays not already allocated in init and nlinit
  !=======================================================================
  !
  CALL allocate_mainvar( ngw, ngwt, ngb, ngs, ngm, nr1, nr2, nr3, dfftp%nr1x, &
                         dfftp%nr2x, dfftp%npl, nnrx, nnrsx, nat, nax, nsp,   &
                         nspin, nbsp, nbspx, n_emp, nupdwn, nkb, gzero, nudx )
  !
  CALL allocate_local_pseudo( ngs, nsp )
  !
  !  initialize wave functions descriptors and allocate wf
  !
  ALLOCATE( c0( ngw, nbspx ) )
  ALLOCATE( cm( ngw, nbspx ) )
  ALLOCATE( cp( ngw, nbspx ) )
  !
  IF ( iprsta > 2 ) THEN
     !
     CALL wave_descriptor_info( wfill, 'wfill', stdout )
     !
  END IF
  !
  acc          = 0.D0
  acc_this_run = 0.D0
  !
  edft%ent  = 0.D0
  edft%esr  = 0.D0
  edft%evdw = 0.D0
  edft%ekin = 0.D0
  edft%enl  = 0.D0
  edft%etot = 0.D0
  !
  ALLOCATE( becsum(  nhm*(nhm+1)/2, nat, nspin ) )
  ALLOCATE( deeq( nhm, nhm, nat, nspin ) )
  IF ( tpre ) THEN
     ALLOCATE( dbec( nkb, nbsp, 3, 3 ) )
     ALLOCATE( drhovan( nhm*(nhm+1)/2, nat, nspin, 3, 3 ) )
  END IF
  !
  ALLOCATE( vkb( ngw, nkb ) )
  !
  IF( program_name == 'CP90' ) THEN
     !
     IF ( tpre ) THEN
        ALLOCATE( drhog( ngm,  nspin, 3, 3 ) )
        ALLOCATE( drhor( nnrx, nspin, 3, 3 ) )
     END IF
     !
  END IF
  !
  IF ( dft_is_meta() .AND. tens ) &
     CALL errore( 'cprmain ', 'ensemble_dft not implimented for metaGGA', 1 )
  !
  IF ( dft_is_meta() .AND. tpre ) THEN
     !
     ALLOCATE( crosstaus( nnrsx, 6, nspin ) )
     ALLOCATE( dkedtaus(  nnrsx, 3, 3, nspin ) )
     ALLOCATE( gradwfc(   nnrsx, 3 ) )
     !
  END IF
  !
  IF ( lwf ) CALL allocate_wannier( nbsp, nnrsx, nspin, ngm )
  !
  IF ( tens .OR. tcg ) &
     CALL allocate_ensemble_dft( nkb, nbsp, ngw, nudx, nspin, nbspx, nnrsx, nat, nlax, nrlx )
  !
  IF ( tcg ) CALL allocate_cg( ngw, nbspx,nkbus )
  !
  IF ( tefield ) CALL allocate_efield( ngw, nbspx, nhm, nax, nsp )
  IF ( tefield2 ) CALL allocate_efield2( ngw, nbspx, nhm, nax, nsp )
  !
  IF ( ALLOCATED( deeq ) ) deeq(:,:,:,:) = 0.D0
  !
  IF ( ALLOCATED( lambda  ) ) lambda  = 0.D0
  IF ( ALLOCATED( lambdam ) ) lambdam = 0.D0
  !
  taum  = tau0
  taup  = 0.D0
  tausm = taus
  tausp = 0.D0
  vels  = 0.D0
  velsm = 0.D0
  velsp = 0.D0
  !
  hnew = h
  !
  cm = ( 0.D0, 0.D0 )
  c0 = ( 0.D0, 0.D0 )
  cp = ( 0.D0, 0.D0 )
  !
  IF ( tens ) then
     CALL id_matrix_init( descla, nspin )
     CALL h_matrix_init( descla, nspin )
  ENDIF
  !
  IF ( lwf ) CALL wannier_startup( ibrav, alat, a1, a2, a3, b1, b2, b3 )
  !
  ! ... Calculate: ema0bg = ecutmass /  MAX( 1.0d0, (2pi/alat)^2 * |G|^2 )
  !
  CALL emass_precond( ema0bg, ggp, ngw, tpiba2, emass_cutoff )
  !
  CALL print_legend( )

  step_con = .FALSE.

  CALL ldau_init()

  IF ( nbeg < 0 ) THEN
     !
     !======================================================================
     !     Initialize from scratch nbeg = -1
     !======================================================================
     !
     nfi = 0
     !
     CALL from_scratch( sfac, eigr, ei1, ei2, ei3, bec, becdr,    &
                        taub, irb, eigrb, b1, b2, b3, nfi,   &
                        rhog, rhor, rhos, rhoc, lambda, lambdam, lambdap, ema0bg,   &
                        dbec, bephi, becp, atoms0, edft, ht0, vpot )
     !
  ELSE
     !
     !======================================================================
     !     nbeg = 0, nbeg = 1
     !======================================================================
     !
     IF( program_name == 'CP90' ) THEN
        !
        i = 1  
        CALL readfile( i, ndr, h, hold, nfi, c0, cm, taus,   &
                       tausm, vels, velsm, acc, lambda, lambdam, xnhe0, xnhem, &
                       vnhe, xnhp0, xnhpm, vnhp,nhpcl,nhpdim,ekincm, xnhh0, xnhhm,&
                       vnhh, velh, ecutp, ecutw, delt, pmass, ibrav, celldm,   &
                       fion, tps, z0t, f )
        !
     ELSE IF( program_name == 'FPMD' ) THEN
        !
        CALL readfile( nfi, tps, c0, cm, f, atoms0, atomsm, acc,     &
                       taui, cdmi, htm, ht0, rhor, vpot, lambda )
        !
     END IF
     !
     CALL from_restart( sfac, eigr, ei1, ei2, ei3, bec, becdr, &
                        taub, irb, eigrb, b1, b2, b3, nfi, rhog, rhor, rhos, &
                        rhoc, lambda, lambdam, lambdap, ema0bg, dbec, bephi, becp, &
                        htm, ht0, vpot, atoms0, edft )
     !
  END IF
  !
  !=======================================================================
  !     restart with new averages and nfi=0
  !=======================================================================
  !
  ! ... reset some variables if nbeg < 0 
  ! ... ( new simulation or step counter reset to 0 )
  !
  IF ( nbeg <= 0 ) THEN
     !
     acc = 0.D0
     nfi = 0
     !
  END IF
  !
  IF ( .NOT. tfor .AND. .NOT. tprnfor ) fion(:,:) = 0.D0
  !
  IF ( tnewnfi ) nfi = newnfi 
  !
  nomore = nomore + nfi
  !
  !  Set center of mass for scaled coordinates
  !
  CALL ions_cofmass( taus, pmass, na, nsp, cdms )
  !
  IF ( nbeg <= 0 .OR. lwf ) THEN
     !
     CALL ions_reference_positions( tau0 )
     !
  END IF
  !
  CALL stop_clock( 'initialize' )
  !
  RETURN
  !
END SUBROUTINE init_run
