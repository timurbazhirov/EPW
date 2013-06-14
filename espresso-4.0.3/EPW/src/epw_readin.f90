  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_readin - Quantum-ESPRESSO group               
  !-----------------------------------------------------------------------
  SUBROUTINE epw_readin
  !-----------------------------------------------------------------------
  !
  !    This routine reads the control variables for the program phononq.
  !    from standard input (unit 5).
  !    A second routine readfile reads the variables saved on a file
  !    by the self-consistent program.
  !
#include "f_defs.h"
  !
  USE ions_base,     ONLY : nat, ntyp => nsp
  USE io_global,     ONLY: ionode_id
  USE mp,            ONLY: mp_bcast, mp_barrier
  USE pwcom
  USE check_stop,    ONLY : time_max => max_seconds
  USE kinds,         ONLY : DP
  USE phcom
  USE epwcom
  USE io_files,      ONLY : tmp_dir, prefix
  USE control_flags, ONLY : iverbosity, modenum, gamma_only
  USE printout_base, ONLY : title
  USE ions_base,     ONLY : amass
  USE mp_global,     ONLY : my_pool_id, me_pool 


  implicit none
  integer :: ios, na, it
  ! integer variable for I/O control
  ! counter on polarizations
  ! counter on iterations
  ! counter on atoms
  ! counter on types
  real(kind=DP), PARAMETER :: ryd2ev = 13.6058
  integer :: modenum_aux, i 
  ! auxilary variable for saving the modenum
  integer :: nk1tmp,nk2tmp,nk3tmp   ! temp vars for saving kgrid info
  character(len=256) :: outdir
  namelist / inputepw / &
       amass, outdir, prefix, iverbosity, time_max, fildvscf,                  &
       tshuffle, tshuffle2, phinterp,elinterp, epstrict,                       &
       elph,  nq1, nq2, nq3,nk1,nk2,nk3,                                       &
       nbndskip,  nbndsub,                                                     &
       tphases, fildvscf0, filukk, filukq,                                     &
       epbread, epbwrite,epwread, epwwrite, epf_mem, etf_mem, kmaps,           &
       eig_read, wepexst, epexst, vme,                                         &
       degaussw, fsthick, eptemp,  nsmear, delta_smear, eminabs, emaxabs,      &
       deltaeabs, dvscf_dir, ngaussw,                                          &
       wannierize, dis_win_max,  dis_win_min, dis_froz_min, dis_froz_max,      &
       num_iter, proj, wdata,iprint,                                           &
       write_wfn, wmin, wmax, nw,                                              &
       eps_acustic, a2f, nest_fn, indabs, twophoton, fly, specfun,             & 
       selfen_type, elecselfen, phonselfen, parallel_k, parallel_q,            &
       rand_q, rand_nq, rand_k, rand_nk,                                       &
       nqf1, nqf2, nqf3, nkf1, nkf2, nkf3,                                     &
       sqf1, sqf2, sqf3, skf1, skf2, skf3,                                     &
       mp_mesh_k, mp_mesh_q, filqf, filkf, filelph, ephmatwrite, band_plot,    &     
       degaussq, delta_qsmear, nqsmear, nqstep
  !
  ! amass    : atomic masses
  ! iverbosity   : verbosity control
  ! outdir   : directory where input, output, temporary files reside
  ! elph     : if true calculate electron-phonon coefficients
  ! time_max : maximum cputime for this run
  ! prefix   : the prefix of files produced by pwscf
  ! filelph  : output file for electron-phonon coefficients
  ! fildvscf : output file containing deltavsc
  ! fildrho  : output file containing deltarho
  !
  ! added by @ FG
  !
  ! tshuffle : elphel2 calculates e-ph matrix elements by using kpoints
  !            of the 1st BZ only (k+q folded back into k+q+G_0)
  ! phinterp : calls elphsum3 instead of elphsum, to calculate e-p matrix
  !            for the small group of q
  ! elinterp : if true perform electron interpolation of e-p matrix
  ! ngaussw  : smearing type for FS average after wann interp
  ! degaussw : corresponding width
  ! filqf    : file with fine q kmesh for interpolation
  ! filkf    : file with fine kmesh for interpolation
  ! filukk   : file with rotation matrix U(k) for interpolation
  ! filukq   : file with rotation matrix U(k+q) for interpolation
  ! tphases  : if true set absolute unitary gauge for eigenvectors
  ! epstrict : if true use strict selection rule for phonon linewidht calculation
  ! fsthick  : the thickness of the Fermi shell for averaging the e-ph matrix elements
  ! eptemp   : temperature for the electronic Fermi occupations in the e-p calculation (units of Kelvin)
  ! fildvscf0: file containing deltavscf to be used as fake perturbation to set phases
  ! nw       : nu. of bins for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  ! wmin     : min frequency for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  ! wmax     : max    "  "  "                         
  ! selfen_type : choice of real/imag part of phonon selfenergy calcuation when epstrict = .true.
  ! nbndsub  : number of bands in the optimal subspace (when disentanglement is used)
  ! tshuffle2: shuffle mode for electrons + load all phonons at once
  ! elecselfen: if .TRUE. calculate imaginary part of electron selfenergy due to e-p interaction
  ! phonselfen: if .TRUE. calculate imaginary part of phonon selfenergy due to e-p interaction
  ! dvscf_dir: the dir containing all the .dvscf and .dyn files
  ! epbread  : read epmatq array from .epb files
  ! epbwrite : write epmatq array to .epb files
  ! nbndskip : number of bands to be skipped from the original Hamitonian (nfirstwin-1 in Marzari's notation)
  ! epwread  : read all quantities in Wannier representation from file epwdata.fmt
  ! epwwrite : write all quantities in Wannier representation to file epwdata.fmt
  ! specfun  : if .TRUE. calculate electron spectral function due to e-p interaction
  !
  !  added by @jn
  !
  ! wannierize : if .TRUE. run the wannier90 code to maximally localize the WFs
  ! dis_win_min : lower bound on wannier90 disentanglement window
  ! dis_win_max : upper bound on wannier90 disentanglement window
  ! dis_froz_min : lower bound on frozen wannier90 disentanglement window
  ! dis_froz_max : upper bound on frozen wannier90 disentanglement window
  ! num_iter     : number of iterations used in the wannier90 minimisation
  ! proj         : initial projections (states) of the wannier functions before minimization
  ! wdata        : Empty array that can be used to pass extra info to prefix.win file, for things not explicitly declared here 
  ! iprint       : verbosity of the wannier90 code
  ! write_wfn    : writes out UNK files from pwscf run for plotting of XSF files
  ! kmaps        : if true, read kmap and kgmap from disk (prior run)
  ! nest_fn      : if true, calculate the nesting function for a given set of q's
  ! fly          : if true, calculate one selfenergy/nestingfn/etc at a time
  ! nsmear       : number of smearing values to use for the selfen_phon call
  ! delta_smear  : change in energy for each additional nsmear
  !
  ! added by RM
  !
  ! ephmatwrite : if true write el-phonon matrix elements on the fine mesh to file
  ! eps_acustic : min phonon frequency for e-p and a2f calculations
  ! band_plot : if true write files to plot band structure and phonon dispersion
  ! degaussq : smearing for sum over q in e-ph coupling (units of meV)
  ! delta_qsmear : change in energy for each additional smearing in the a2f (units of meV)
  ! nqsmear : number of smearings used to calculate a2f
  ! nqstep : number of steps used to calculate a2f
  !
  CHARACTER (LEN=80)  :: input_file
  INTEGER             :: nargs, iiarg, ierr
  !
  INTEGER(4) :: iargc
  !
  nk1tmp = 0
  nk2tmp = 0
  nk3tmp = 0
#ifdef __PARA
  IF (me_pool /=0 .or. my_pool_id /=0) goto 400 
#endif
  !
  !
  ! ... Input from file ?
  !
  nargs = iargc() 
  !
  DO iiarg = 1, ( nargs - 1 )
     !
     CALL getarg( iiarg, input_file )  
     IF ( TRIM( input_file ) == '-input' .OR. &
          TRIM( input_file ) == '-inp'   .OR. &
          TRIM( input_file ) == '-in' ) THEN
        !
        CALL getarg( ( iiarg + 1 ) , input_file )  
        OPEN ( UNIT = 5, FILE = input_file, FORM = 'FORMATTED', &
               STATUS = 'OLD', IOSTAT = ierr )
        CALL errore( 'iosys', 'input file ' // TRIM( input_file ) // &
                   & ' not found' , ierr )
        !
     END IF
     !
  END DO
  !
  !
  !    Read the first line of the input file
  !
  READ (5, '(a)', err = 100, iostat = ios) title
100 CALL errore ('epw_readin', 'reading title ', abs (ios) )
  !
  !   set default values for variables in namelist
  !
  amass(:)     = 0.d0
  iverbosity   = 0
  elph         = .false.
  tshuffle     = .false.
  tshuffle2    = .true.
  elecselfen   = .false.
  phonselfen   = .false.
  specfun      = .false.
  epbread      = .false.
  epbwrite     = .false.
  epwread      = .false.
  epwwrite     = .false.
  phinterp     = .true.
  wannierize   = .false.
  write_wfn    = .false.
  kmaps        = .false.
  nest_fn      = .false.
  wepexst      = .false.
  epexst       = .false.
  indabs       = .false.
  twophoton    = .false.
  eig_read     = .false.
  fly          = .true.
  dis_win_max  = 1d3
  dis_win_min  = -1d3
  dis_froz_max = -0.9d3
  dis_froz_min = -1d3
  num_iter     = 200
  proj(:)      = ''
  wdata(:)     = ''
  iprint       = 2
  wmin         = 0.d0
  wmax         = 0.02
  eps_acustic  = 20.d0 
  nw           = 10
  fsthick      = 1.d10
  eptemp(:)    = 0.000
  eptemp(1)    = 300.
  degaussw     = 0.025
  epstrict     = .false.
  selfen_type  = 2
  elinterp     = .true.
  tphases      = .false.
  parallel_k   = .true.
  parallel_q   = .false.
  a2f          = .false.
  epf_mem      = .false.
  etf_mem      = .true.
  fildvscf0    = ' '
  ngaussw      = 1
  time_max     = 10000000.d0
  outdir       = '.'
  dvscf_dir    = '.'
  prefix       = 'pwscf'
  filelph      = ' '
  fildrho      = ' '
  fildvscf     = ' '
  filukk       = ' '
  rand_q       = .false.
  rand_nq       = 1
  rand_k       = .false.
  rand_nk       = 1
  nq1          = 0
  nq2          = 0
  nq3          = 0
  nk1          = 0
  nk2          = 0
  nk3          = 0
  nqf1         = 0
  nqf2         = 0
  nqf3         = 0
  nkf1         = 0
  nkf2         = 0
  nkf3         = 0
  mp_mesh_k    = .false.
  mp_mesh_q    = .false.
  nbndsub      = 0
  nbndskip     = 0
  nsmear       = 1
  delta_smear  = 0.001
  modenum = -1
  vme = .false.
  eminabs = 0.d0
  emaxabs = 3.d0
  deltaeabs = 0.05
  ephmatwrite = .false.
  band_plot = .false.
  nqsmear = 10
  nqstep = 500
  delta_qsmear = 0.05 ! meV 
  degaussq = 0.05 ! meV
  !
  !     reading the namelist inputepw
  !
#ifdef CRAYY
  !   The Cray does not accept "err" and "iostat" together with a namelist
  READ (5, inputepw)
  ios = 0
#else
  !
  READ (5, inputepw, err = 200, iostat = ios)
#endif

200 CALL errore ('epw_readin', 'reading input_epw namelist', abs (ios) )
  !
  nk1tmp = nk1
  nk2tmp = nk2
  nk3tmp = nk3
  !
  !
  !     Check all namelist variables
  !
  IF (filukk.eq.' ') filukk=trim(prefix)//'.ukk'
  IF (nsmear .lt. 1) CALL errore ('epw_readin', &
       & 'Wrong number of nsmears',1)
  IF (iverbosity.ne.0.and.iverbosity.ne.1) CALL errore ('epw_readin', &
       &' Wrong  iverbosity ', 1)
  IF (time_max.lt.1.d0) CALL errore ('epw_readin', ' Wrong time_max', 1)
  IF (tphases.and.fildvscf0.eq.' ') CALL errore ('epw_readin', &
       &' tphases requires fildvscf0', 1)
  IF (epbread.and.epbwrite) CALL errore ('epw_readin', &
       &' epbread cannot be used with epbwrite', 1)
  IF (degaussw*4.d0.gt.fsthick) CALL errore ('epw_readin', &
       &' degaussw too close to fsthick', 1)
  IF ( nbndskip .lt. 0) CALL errore('epw_readin', &
       &' nbndskip must not be less than 0', 1)
  IF ((nw.lt.1).or.(nw.gt.1000)) CALL errore ('epw_readin', &
       &' unreasonable nw', 1)
  IF (parallel_k .and. parallel_q) CALL errore('epw_readin', &
       &'can only parallelize over k OR q',1)
  IF (.not.(parallel_k .or. parallel_q)) CALL errore('epw_readin', &
       &'must parallelize over k OR q',1)
  IF (parallel_k .and. elecselfen) CALL errore('epw_readin', &
       &'Electron selfenergy is more efficient with k_para',-1)
  IF (fly .and. elecselfen) CALL errore('epw_readin', &
       &'Electron selfenergy requires fly =.false.',-1)
  IF (parallel_q .and. .not.epf_mem) CALL errore('epw_readin', &
       &'Must store ep-matrix in memory for parallel q',1)
  IF (a2f .and. .not.phonselfen) CALL errore('epw_readin', &
       &'a2f requires phonoselfen',1)
!  if (rand_q .and. .not.fly) CALL errore ('epw_readin', &
!       'rand_q meant to be run with fly=true', 1)
  IF (parallel_q .and. phonselfen) CALL errore('epw_readin', &
       &'Phonon selfenergy is more efficient with q_para',-1)
  IF (parallel_q) CALL errore('epw_readin', &
       &'WARNING: Parallel q not tested!', -1)
  IF ( (elph.and.wannierize) .and. (epwread) ) CALL errore('epw_readin', &
       & 'must use same w90 rotation matrix for entire run', 1)
  IF ((wmin.gt.wmax)) &
       CALL errore ('epw_readin', ' check wmin, wmax ', 1)
  IF ((nbndsub.gt.200)) CALL errore ('epw_readin', & 
       ' too many wannier functions increase size of projx', 1)
  IF (indabs .and. twophoton)  CALL errore ('epw_readin', 'twophoton and indabs not used together',1 )
  IF ( elph .and. ( mp_mesh_k .or. mp_mesh_q )) CALL errore('epw_readin', &
       &'can only work with full uniform mesh',1)
  IF (ephmatwrite .and. .not.fly ) CALL errore('epw_readin', &
      &'ephmatwrite requires fly=.true.',1)
  IF (fly .and. band_plot) CALL errore('epw_readin', &
      &'plot band structure and phonon dispersion requires fly =.false.',1)
  !
  !
  DO i = 1, ntempxx
     IF (eptemp(i) .gt. 0.d0) THEN
        neptemp = i
        ! 1 K in eV = 8.6173423e-5
        ! from K to Ryd
        eptemp(i) = eptemp(i) * 0.000086173423 / ryd2eV
     ENDIF
  ENDDO
  ! from cm-1 to Ryd
  eps_acustic = eps_acustic / 13.6058d0 / 8065.5d0
  !
  !    reads the q point (just if ldisp = .false.)
  !
xq(:) = 0.d0
  !
300 CALL errore ('epw_readin', 'reading xq', abs (ios) )
  lgamma = .false.
  tmp_dir = trim(outdir)
  dvscf_dir = trim(dvscf_dir)//'/'
  !
#ifdef __PARA
400 continue
  CALL bcast_ph_input
#endif
  xqq(:) = xq(:) 
  !
  !   Here we finished the reading of the input file.
  !   Now allocate space for pwscf variables, read and check them.
  !
  modenum_aux = modenum
  !
  CALL read_file
  !
  ! nbnd comes out of readfile
  IF (nbndsub.eq.0) nbndsub = nbnd
  !
#ifdef __PARA
  IF (.not.(me_pool /=0 .or. my_pool_id /=0)) THEN
     nk1 = nk1tmp
     nk2 = nk2tmp
     nk3 = nk3tmp
  ENDIF
#else
     nk1 = nk1tmp
     nk2 = nk2tmp
     nk3 = nk3tmp
#endif

  !
  IF (gamma_only) CALL errore('epw_readin',&
     'cannot start from pw.x data file using Gamma-point tricks',1)
  !
  IF (modenum_aux .ne. -1) THEN
     modenum = modenum_aux
     iswitch = -4
  ELSE IF (modenum .eq. 0) THEN
     iswitch = -2
  ELSE
     iswitch = -4
  END IF
  !
  CALL mp_bcast ( iswitch, ionode_id )
  CALL mp_bcast ( modenum, ionode_id )
  !
  IF (tfixed_occ) &
     CALL errore('epw_readin','phonon with arbitrary occupations not tested',1)
  !
  IF (elph.and.lsda) CALL errore ('epw_readin', 'El-ph and spin not &
       &implemented', 1)
  !
  !   There might be other variables in the input file which describe
  !   partial computation of the dynamical matrix. Read them here
  !
  CALL allocate_part
#ifdef __PARA

   IF (me_pool /= 0 .or. my_pool_id /=0) goto 800 
#endif
  IF (nrapp.lt.0.or.nrapp.gt.3 * nat) CALL errore ('epw_readin', &
       'nrapp is wrong', 1)
  IF (nrapp.ne.0) THEN
     read (5, *, err = 700, iostat = ios) (list (na), na = 1, nrapp)
700  CALL errore ('epw_readin', 'reading list', abs (ios) )

  ENDIF
#ifdef __PARA
800 continue
  CALL bcast_ph_input1
#endif


  DO it = 1, ntyp
     IF (amass (it) <= 0.d0) CALL errore ('epw_readin', 'Wrong masses', it)
  ENDDO

  IF (maxirr.lt.0.or.maxirr.gt.3 * nat) CALL errore ('epw_readin', ' &
       &Wrong maxirr ', abs (maxirr) )
! No k+q anymore, can we cut?
  IF (mod (nks, 2) .ne.0) CALL errore ('epw_readin', &
      'k-points are odd', 0)
  IF (iswitch .eq. -4 ) THEN
     nrapp = 1
     list(1) = modenum
  ENDIF
  !
#ifdef __PARA
  !
  !  broadcast the values of nq1, nq2, nq3
  !
  CALL mp_bcast( nq1, ionode_id )
  CALL mp_bcast( nq2, ionode_id )
  CALL mp_bcast( nq3, ionode_id )
  CALL mp_bcast( nk1, ionode_id )
  CALL mp_bcast( nk2, ionode_id )
  CALL mp_bcast( nk3, ionode_id )
#endif
  !
  amass = amconv * amass
  !
  END SUBROUTINE epw_readin
