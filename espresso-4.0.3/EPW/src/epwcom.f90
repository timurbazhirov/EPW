  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phcom - Quantum-ESPRESSO group                
  !-----------------------------------------------------------------------
  MODULE control_epw
  !-----------------------------------------------------------------------
  !
  ! Common variables for the epw program
  !  
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY :  DP
  USE parameters, ONLY: npk
  !
  ! ... the variable controlling the phonon run
  !
  SAVE
  ! 
  INTEGER :: ngaussw, nw, selfen_type, nbndsub, nbndskip, num_iter, iprint, &
       nsmear, rand_nq, rand_nk, nqf1, nqf2, nqf3, nkf1, nkf2, nkf3, &
       sqf1, sqf2, sqf3, skf1, skf2, skf3, nqsmear, nqstep, neptemp
  ! smearing type for Fermi surface average in e-ph coupling after wann. interp.
  ! nu. of bins for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  ! choice of real/imag part of phonon selfenergy when epstrict = .true.
  ! number of bands in the optimal subspace (when disentanglement is used)
  ! the number of bands to be skipped when we use only a subspace (this is nfirstwin-1 in Marzari's notation)
  ! @JN
  ! the number of steps used in creating the maximally localized Wannier functions
  ! the verbosity of the wannier90 code
  ! number of temperature-like smearings calculated in selfen_phon, indabs
  ! use random points for the fine q-mesh
  ! use random points for the fine k-mesh
  ! kx,ky,kz sizes of the uniform electron fine mesh to be used
  ! qx,qy,qz sizes of the uniform phonon fine mesh to be used
  ! sqf1, sqf2, sqf3 shift of the uniform phonon fine mesh
  ! skf1, skf2, skf3 shift of the uniform electron fine mesh
  ! number of smearings used to calculate a2f 
  ! number of steps used to calculate a2f 
  !
  INTEGER :: nwanxx = 200, ntempxx = 25
  ! parameter used in writing prefix.win file.  Maximum number of wannier functions
  REAL (KIND=DP) :: degaussw, fsthick, wmin, wmax, ecutse, dis_win_min, dis_win_max, &
       dis_froz_min, dis_froz_max, delta_smear, eminabs, emaxabs, deltaeabs, eps_acustic, &
       degaussq, delta_qsmear
  ! smearing width for Fermi surface average in e-ph coupling after wann interp
  ! thickness of the Fermi shell for averaging the e-ph matrix element
  ! temperature for the electronic Fermi occupations in the e-p calculation 
  ! min/max frequency for frequency scan in \delta( e_k - e_k+q - w ) when strict sel. rule is applied
  ! minimum energy of the Wannier disentanglement window
  ! maximum energy of the Wannier disentanglement window
  ! minimum energy of the frozen Wannier disentanglement window
  ! maximum energy of the frozen Wannier disentanglement window
  ! change in energy for each additional smearing in the selfen_phon and indabs
  ! min. phonon frequency for e-p and a2f calculations
  ! smearing for sum over q in e-ph coupling  
  ! change in energy for each additional smearing in the a2f
  !
  logical :: phinterp , elinterp, tphases, epstrict, tshuffle2, elecselfen, phonselfen, &
       epbread, epbwrite, epwread, epwwrite, specfun, wannierize, parallel_k,      &
       parallel_q, a2f, epf_mem, etf_mem, write_wfn, kmaps, nest_fn, fly, rand_q, rand_k, &
       mp_mesh_q, mp_mesh_k, indabs, eig_read, wepexst, epexst, vme, twophoton, ephmatwrite, & 
       band_plot
  ! if .TRUE. call elphsum3 instead of elphsum, to calculate the e-p matrix
  !           for the small group of q
  ! if .TRUE. perform phonon interpolation of e-p matrix
  ! if .TRUE. perform electron interpolation of e-p matrix
  ! if .TRUE. set absolute reference for unitary gauge of the eigenvectors
  ! if .TRUE. use strict selection rule in phonon linewidth calculation
  ! if .TRUE. calculate electron selfenergy due to e-p interaction
  ! if .TRUE. calculate phonon selfenergy due to e-p interaction
  ! if .TRUE. read epmatq from files .epb
  ! if .TRUE. write epmatq to files .epb
  ! if .TRUE. read all quantities in Wannier representation from file epwdata.fmt
  ! if .TRUE. write all quantities in Wannier representation to file epwdata.fmt
  ! if .TRUE. calculate spectral electron function due to e-p interaction
  ! if .TRUE. run the wannier90 code
  ! if .TRUE. scatter the electron k-points on the fine mesh among pools (not q)
  ! if .TRUE. scatter the phonon q-points on the fine mesh among pools (not k)
  ! if .TRUE. calculate Eliashberg spectral electron function from selfen_phon
  ! if .TRUE. store fine mesh eigenvalues and matrix elements in memory
  ! if .TRUE. write out UNK files in wannier90
  ! if .TRUE. read kmap and kgmap from disk.  Do not calculate
  ! if .TRUE. calculate the electronic nesting function (metals only)
  ! if .TRUE. do not interpolate all k and q's prior to post-proc calls.  saves disk space
  ! if .TRUE. use random points for the fine q-mesh
  ! if .TRUE. use random points for the fine k-mesh
  ! if .TRUE. use points in the irreducible wedge for the uniform fine mesh
  ! if .TRUE. calculated indirect absorption spectrum
  ! if .true. (eig_read) then readin a set of electronic eigenvalues in eV to replace the calcualted ones
  ! if .TRUE. prefix.epmatwe files are already on disk.  don't recalculate. debugging param
  ! if .TRUE. prefix.epmatwp files are already on disk.  don't recalculate  debugging param
  ! if .TRUE. write el-ph matrix elements on the fine mesh to file
  ! if .TRUE. write filrs to plot band structure and phonon dispersion
  CHARACTER(len=80) :: dvscf_dir ='./' ! directory for .dvscf and .dyn files (wannier interpolation)
  CHARACTER(len=80) :: filelph, fileig ! output file for the electron-phonon coefficients
  character(len=256), dimension(200) :: proj, wdata ! projections and any extra info for W90 
  REAL (kind=DP), dimension(25) :: eptemp 
  integer :: iswitch

END MODULE control_epw
!
!
MODULE klist_epw
! The variable for the k-points 
USE kinds, ONLY: DP
USE parameters, ONLY :npk
!
SAVE
!
INTEGER :: kmap(npk)  ! map of k+q grid into k grid (only when tshuffle=.true.)
LOGICAL :: tshuffle   ! if .TRUE. refold k+q grid into k grid
REAL(DP) :: xk_cryst(3,npk) ! List of all kpoints in crystal coordinates

END MODULE klist_epw

!
MODULE units_epw
  !
  ! ... the units of the files and the record lengths
  !
  SAVE
  !
  INTEGER :: iudvscf0, iuncuf, iunepmatf, &
       lrcuf, lrepmatf, lretf, iuetf, iunepmatwe, iunepmatwp
  !
  !
  ! the unit where the delta Vscf is read to generate the fake perturbation 
  ! unit with rotation matrix on fine mesh
  ! unit with rotation matrix on fine mesh
  ! the length of the record for the rotation matrix
  ! the length of the record for the electron-phonon matrix
  ! the length of the record for the interpolated hamiltonian eigenvalues
  ! the unit with the interpolated hamiltonian eigenvalues
  ! the unit with the e-ph matrix in Wannier-Bloch representation
  ! the unit with the e-ph matrix in Wannier-Wannier representation
  !
  logical, ALLOCATABLE :: this_dvkb3_is_on_file(:), &
                          this_pcxpsi_is_on_file(:,:)
  !
END MODULE units_epw
!
!
MODULE output_epw
  !
  ! ... the name of the files
  !
  SAVE
  !
  CHARACTER (LEN=80) :: filqf, filkf, filukk, filukq, fildvscf0
  !
  ! input  file for the fine q mesh
  ! input  file for the fine k mesh
  ! input  file for the rotation matrix U(k)
  ! input  file for the rotation matrix U(k+q)
  ! input  file for wigner-seitx cell for wannier
  ! output file for the deltavscf used as a fake perturbation to set phases
  !
END MODULE output_epw
!
MODULE epwcom
  USE control_epw
  USE units_epw
  USE output_epw
  USE klist_epw
  END MODULE epwcom
