!
! Copyright (C) 2001-2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!----------------------------------------------------------------------------
SUBROUTINE phq_readin()
  !----------------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program phononq.
  !    from standard input (unit 5).
  !    A second routine readfile reads the variables saved on a file
  !    by the self-consistent program.
  !
  !
  USE kinds,         ONLY : DP
  USE parameters,    ONLY : nsx
  USE constants,     ONLY : amconv
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE io_global,     ONLY : ionode_id
  USE mp,            ONLY : mp_bcast
  USE input_parameters, ONLY : max_seconds
  USE ions_base,     ONLY : amass, pmass, atm
  USE klist,         ONLY : xqq, xk, nks, nkstot, lgauss, two_fermi_energies
  USE control_flags, ONLY : gamma_only, tqr
  USE uspp,          ONLY : okvan
  USE fixed_occ,     ONLY : tfixed_occ
  USE lsda_mod,      ONLY : lsda, nspin
  USE printout_base, ONLY : title
  USE control_ph,    ONLY : maxter, alpha_mix, lgamma, lgamma_gamma, epsil, &
                            zue, trans, reduce_io, &
                            elph, tr2_ph, niter_ph, nmix_ph, maxirr, lnscf, &
                            ldisp, recover, lrpa, lnoloc
  USE gamma_gamma,   ONLY : asr
  USE qpoint,        ONLY : nksq, xq
  USE partial,       ONLY : atomo, list, nat_todo, nrapp
  USE output,        ONLY : fildyn, fildvscf, fildrho
  USE disp,          ONLY : nq1, nq2, nq3, iq1, iq2, iq3
  USE io_files,      ONLY : tmp_dir, prefix, trimcheck
  USE noncollin_module, ONLY : i_cons
  USE ldaU,          ONLY : lda_plus_u
  USE control_flags, ONLY : iverbosity, modenum
  USE io_global,     ONLY : ionode
  USE mp_global,     ONLY : nproc, nproc_pool, nproc_file, nproc_pool_file
  USE control_flags, ONLY : twfcollect
  USE paw_variables, ONLY : okpaw
  USE ramanm,        ONLY : eth_rps, eth_ns, lraman, elop, dek
  USE freq_ph
  !
  IMPLICIT NONE
  !
  INTEGER :: ios, ipol, iter, na, it
    ! integer variable for I/O control
    ! counter on polarizations
    ! counter on iterations
    ! counter on atoms
    ! counter on types
  REAL(DP) :: amass_input(nsx)
    ! save masses read from input here
  CHARACTER (LEN=256) :: outdir
  !
  CHARACTER(LEN=256)         :: input_line
  CHARACTER(LEN=80)          :: card
  CHARACTER(LEN=1), EXTERNAL :: capital
  LOGICAL                    :: tend
  LOGICAL                    :: end_of_file
  INTEGER                    :: i
  INTEGER, EXTERNAL :: atomic_number
  REAL(DP), EXTERNAL :: atom_weight

  !
  NAMELIST / INPUTPH / tr2_ph, amass, alpha_mix, niter_ph, nmix_ph,  &
                       maxirr, nat_todo, iverbosity, outdir, epsil,  &
                       trans, elph, zue, nrapp, max_seconds, reduce_io, &
                       prefix, fildyn, fildvscf, fildrho,            &
                       lnscf, ldisp, nq1, nq2, nq3, iq1, iq2, iq3,   &
                       eth_rps, eth_ns, lraman, elop, dek, recover,  &
                       fpol, asr, lrpa, lnoloc
  ! tr2_ph       : convergence threshold
  ! amass        : atomic masses
  ! alpha_mix    : the mixing parameter
  ! niter_ph     : maximum number of iterations
  ! nmix_ph      : number of previous iterations used in mixing
  ! maxirr       : the number of irreducible representations
  ! nat_todo     : number of atom to be displaced
  ! iverbosity   : verbosity control
  ! outdir       : directory where input, output, temporary files reside
  ! epsil        : if true calculate dielectric constant
  ! trans        : if true calculate phonon
  ! elph         : if true calculate electron-phonon coefficients
  ! zue          : if true calculate effective charges (alternate way)
  ! lraman       : if true calculate raman tensor
  ! elop         : if true calculate electro-optic tensor
  ! nrapp        : the representations to do
  ! max_seconds  : maximum cputime for this run
  ! reduce_io    : reduce I/O to the strict minimum
  ! prefix       : the prefix of files produced by pwscf
  ! fildyn       : output file for the dynamical matrix
  ! fildvscf     : output file containing deltavsc
  ! fildrho      : output file containing deltarho
  ! eth_rps      : threshold for calculation of  Pc R |psi> (Raman)
  ! eth_ns       : threshold for non-scf wavefunction calculation (Raman)
  ! dek          : delta_xk used for wavefunctions derivation (Raman)
  ! recover      : recover=.true. to restart from an interrupted run
  ! asr          : in the gamma_gamma case apply acoustic sum rule
  !
  IF ( .NOT. ionode ) GOTO 400
  !
  ! ... Input from file ?
  !
  CALL input_from_file ( )
  !
  ! ... Read the first line of the input file
  !
  READ( 5, '(A)', ERR = 100, IOSTAT = ios ) title
  !
100 CALL errore( 'phq_readin', 'reading title ', ABS( ios ) )
  !
  ! ... set default values for variables in namelist
  !
  tr2_ph       = 1.D-10
  eth_rps      = 1.D-9
  eth_ns       = 1.D-12
  amass(:)     = 0.D0
  alpha_mix(:) = 0.D0
  alpha_mix(1) = 0.7D0
  niter_ph     = maxter
  nmix_ph      = 4
  maxirr       = 0
  nat_todo     = 0
  nrapp        = 0
  iverbosity   = 0
  trans        = .TRUE.
  lrpa         = .FALSE.
  lnoloc       = .FALSE.
  epsil        = .FALSE.
  zue          = .FALSE.
  fpol         = .FALSE.
  elph         = .FALSE.
  lraman       = .FALSE.
  elop         = .FALSE.
  max_seconds  =  1.E+7_DP
  reduce_io    = .FALSE.
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  prefix       = 'pwscf'
  fildyn       = 'matdyn'
  fildrho      = ' '
  fildvscf     = ' '
  lnscf        = .FALSE.
  ldisp        = .FALSE.
  nq1          = 0
  nq2          = 0
  nq3          = 0
  iq1          = 0
  iq2          = 0
  iq3          = 0
  dek          = 1.0d-3
  recover      = .FALSE.
  asr          = .FALSE.
  !
  ! ...  reading the namelist inputph
  !
  READ( 5, INPUTPH, ERR = 200, IOSTAT = ios )
  !
200 CALL errore( 'phq_readin', 'reading inputph namelist', ABS( ios ) )
  !
  ! ... Check all namelist variables
  !
  IF (tr2_ph <= 0.D0) CALL errore (' phq_readin', ' Wrong tr2_ph ', 1)
  IF (eth_rps<= 0.D0) CALL errore ( 'phq_readin', ' Wrong eth_rps', 1)
  IF (eth_ns <= 0.D0) CALL errore ( 'phq_readin', ' Wrong eth_ns ', 1)

  DO iter = 1, maxter
     IF (alpha_mix (iter) .LT.0.D0.OR.alpha_mix (iter) .GT.1.D0) CALL &
          errore ('phq_readin', ' Wrong alpha_mix ', iter)
  ENDDO
  IF (niter_ph.LT.1.OR.niter_ph.GT.maxter) CALL errore ('phq_readin', &
       ' Wrong niter_ph ', 1)
  IF (nmix_ph.LT.1.OR.nmix_ph.GT.5) CALL errore ('phq_readin', ' Wrong &
       &nmix_ph ', 1)
  IF (iverbosity.NE.0.AND.iverbosity.NE.1) CALL errore ('phq_readin', &
       &' Wrong  iverbosity ', 1)
  IF (fildyn.EQ.' ') CALL errore ('phq_readin', ' Wrong fildyn ', 1)
  IF (max_seconds.LT.1.D0) CALL errore ('phq_readin', ' Wrong max_seconds', 1)

  IF (nat_todo.NE.0.AND.nrapp.NE.0) CALL errore ('phq_readin', &
       &' incompatible flags', 1)
  IF (dek <= 0.d0) CALL errore ( 'phq_readin', ' Wrong dek ', 1)
  epsil = epsil .OR. lraman .OR. elop
  IF ( (lraman.OR.elop) .AND. fildrho == ' ') fildrho = 'drho'
  !
  !    reads the q point (just if ldisp = .false.)
  !
  IF (.NOT. ldisp) THEN
     READ (5, *, err = 300, iostat = ios) (xq (ipol), ipol = 1, 3)
  END IF
300 CALL errore ('phq_readin', 'reading xq', ABS (ios) )
  lgamma = xq (1) .EQ.0.D0.AND.xq (2) .EQ.0.D0.AND.xq (3) .EQ.0.D0
  IF ( (epsil.OR.zue) .AND..NOT.lgamma) CALL errore ('phq_readin', &
       'gamma is needed for elec.field', 1)
  IF (zue.AND..NOT.trans) CALL errore ('phq_readin', 'trans must be &
       &.t. for Zue calc.', 1)

  IF (trans.AND.(lrpa.OR.lnoloc)) CALL errore('phq_readin', &
                    'only dielectric constant with lrpa or lnoloc',1)
  !
  ! reads the frequencies ( just if fpol = .true. )
  !
  IF ( fpol ) THEN
     READ (5, *, err = 10, iostat = ios) card
     IF ( TRIM(card) == 'FREQUENCIES' ) THEN
        READ (5, *, err = 10, iostat = ios) nfs
        DO i = 1, nfs
           READ (5, *, err = 10, iostat = ios) fiu(i)
        END DO
     END IF
  END IF
10 CALL errore ('phq_readin', 'reading FREQUENCIES card', ABS(ios) )
  !
  tmp_dir = trimcheck (outdir)
  !
400 CONTINUE
  CALL bcast_ph_input ( ) 
  call mp_bcast ( fpol, ionode_id )
  xqq(:) = xq(:) 
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  !   amass will also be read from file:
  !   save its content in auxiliary variables
  !
  amass_input(:)= amass(:)
  !
  CALL read_file ( )
  !
  IF (gamma_only) CALL errore('phq_readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)

  IF (lda_plus_u) CALL errore('phq_readin',&
     'The phonon code with LDA+U is not yet available',1)

  IF (okpaw)  CALL errore('phq_reading','phonon with paw not yet available',1)

  IF (nproc /= nproc_file .and. .not. twfcollect)  &
     CALL errore('phq_readin',&
     'pw.x run with a different number of processors. Use wf_collect=.true.',1)

  IF (nproc_pool /= nproc_pool_file .and. .not. twfcollect)  &
     CALL errore('phq_readin',&
     'pw.x run with a different number of pools. Use wf_collect=.true.',1)

  IF (two_fermi_energies.or.i_cons /= 0) &
     CALL errore('phq_readin',&
     'The phonon code with constrained magnetization is not yet available',1)

  IF (tqr) CALL errore('phq_readin',&
     'The phonon code with Q in real space not available',1)
  !
  !  set masses to values read from input, if available;
  !  leave values read from file otherwise
  !
  DO it = 1, ntyp
     IF (amass_input(it) < 0.0_DP) amass_input(it)= &
              atom_weight(atomic_number(TRIM(atm(it))))
     IF (amass_input(it) > 0.D0) amass(it) = amass_input(it)
     IF (amass(it) <= 0.D0) CALL errore ('phq_readin', 'Wrong masses', it)
     !
     !  convert masses to a.u.
     !
     pmass(it) = amconv * amass(it)
  ENDDO
  lgamma_gamma=.FALSE.
  IF (nkstot==1.OR.(nkstot==2.AND.nspin==2)) THEN
     lgamma_gamma=(lgamma.AND.(ABS(xk(1,1))<1.D-12) &
                         .AND.(ABS(xk(2,1))<1.D-12) &
                         .AND.(ABS(xk(3,1))<1.D-12) )
  ENDIF
  !
  IF (lgamma) THEN
     nksq = nks
  ELSE
     nksq = nks / 2
  ENDIF
  !
  IF (tfixed_occ) &
     CALL errore('phq_readin','phonon with arbitrary occupations not tested',1)
  !
  IF (elph.AND..NOT.lgauss) CALL errore ('phq_readin', 'Electron-&
       &phonon only for metals', 1)
  IF (elph.AND.lsda) CALL errore ('phq_readin', 'El-ph and spin not &
       &implemented', 1)
  IF (elph.AND.fildvscf.EQ.' ') CALL errore ('phq_readin', 'El-ph needs &
       &a DeltaVscf file', 1)
  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  CALL allocate_part ( )
  !
  IF ( .NOT. ionode ) GOTO 800

  IF (nat_todo.LT.0.OR.nat_todo.GT.nat) CALL errore ('phq_readin', &
       'nat_todo is wrong', 1)
  IF (nat_todo.NE.0) THEN
     READ (5, *, err = 600, iostat = ios) (atomo (na), na = 1, &
          nat_todo)
600  CALL errore ('phq_readin', 'reading atoms', ABS (ios) )
  ENDIF
  IF (nrapp.LT.0.OR.nrapp.GT.3 * nat) CALL errore ('phq_readin', &
       'nrapp is wrong', 1)
  IF (nrapp.NE.0) THEN
     READ (5, *, err = 700, iostat = ios) (list (na), na = 1, nrapp)
700  CALL errore ('phq_readin', 'reading list', ABS (ios) )

  ENDIF

800 CONTINUE

  CALL bcast_ph_input1 ( ) 

  IF (epsil.AND.lgauss) &
        CALL errore ('phq_readin', 'no elec. field with metals', 1)
  IF (maxirr.LT.0.OR.maxirr.GT.3 * nat) CALL errore ('phq_readin', ' &
       &Wrong maxirr ', ABS (maxirr) )
  IF (MOD (nkstot, 2) .NE.0.AND..NOT.lgamma.and..not.lnscf) &
           CALL errore ('phq_readin', 'k-points are odd', nkstot)
  IF (modenum > 0) THEN
     IF ( ldisp ) &
          CALL errore('phq_readin','Dispersion calculation and &
          & single mode calculation not possibile !',1)
     nrapp = 1
     nat_todo = 0
     list (1) = modenum
  ENDIF
  
  IF (modenum > 0 .OR. ldisp .OR. lraman ) lgamma_gamma=.FALSE.
  IF (.not.lgamma_gamma) asr=.FALSE.
  !
  !  broadcast the values of nq1, nq2, nq3, iq1, iq2, iq3
  !
  CALL mp_bcast( nq1, ionode_id )
  CALL mp_bcast( nq2, ionode_id )
  CALL mp_bcast( nq3, ionode_id )
  CALL mp_bcast( iq1, ionode_id )
  CALL mp_bcast( iq2, ionode_id )
  CALL mp_bcast( iq3, ionode_id )
  !
  IF (ldisp .AND. (nq1 .LE. 0 .OR. nq2 .LE. 0 .OR. nq3 .LE. 0)) &
       CALL errore('phq_readin','nq1, nq2, and nq3 must be greater than 0',1)
  !
  RETURN
  !
END SUBROUTINE phq_readin
