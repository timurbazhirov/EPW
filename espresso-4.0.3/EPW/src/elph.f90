  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------------
  MODULE el_phon
  !-------------------------------------------------------------------------- 
  USE kinds, ONLY :  DP
  !
  SAVE
  !
  COMPLEX(KIND=DP), ALLOCATABLE :: &
       el_ph_mat  (:,:,:,:),   &!  e-p matrix  (nbnd, nbnd, nks, 3*nat)
       cu(:,:,:),              &!  rot matrix for wannier interpolation of k point, coarse mesh (nbnd*nbnd*nkstot)
       cuq(:,:,:),             &!  rot matrix for wannier interpolation of k+q point, coarse mesh (nbnd*nbnd*nkstot)
       chw (:,:,:),            &!  Hamiltonian in wannier basis 
       chw_ks (:,:,:),         &!  Hamiltonian in wannier basis  (Kohn-Sham eigenvalues if many-body eigenvalues are read in)
       cdmew (:,:,:,:),        &!  Dipole matrix in wannier basis 
       cvmew (:,:,:,:),        &!  Velocity matrix in wannier basis 
       rdw (:,:,:),            &!  dynamical matrix in wannier basis (real) 
       epmatwp (:,:,:,:),      &!  e-p matrix  in wannier basis - electrons and phonons
       umat(:,:,:),            &!  the rotation matrix for the unique setting of the wfs gauge -- on the local pool
       umatq(:,:,:),           &!  the rotation matrix for the unique setting of the wfs gauge for the k + q-- on the local pool
       umat_all(:,:,:),        &!  the rotation matrix for the unique setting of the wfs gauge -- for all k points
       umatq_all(:,:,:),       &!  the rotation matrix for the unique setting of the wfs gauge -- for all k+q points
       dynq  (:,:,:),          &!  dynamical matrix for every q (nmode, nmodes, nqtot)
       epmatq (:,:,:,:,:),     &!  e-p matrix for every q (nbnd, nbnd, nks, nmodes, nqtot)
       epmatw17 (:,:,:,:,:),   &!  the ep matrix in the wannier representation (nwann, nwann, nrr_k, nrr_q, nmodes) 
       epf17 (:, :, :, :, :),  &!  full ep matrix in bloch rep stored in mem (nkstotf, nqstotf, nbnd, nbnd, nmodes)-nbnd inside wndw 
       dmec(:,:,:,:),          &!  dipole matrix elements on the course mesh (ipol, nbnd, nbnd, nks)
       dmef(:,:,:,:),          &!  dipole matrix elements on the fine   mesh (ipol, nbnd, nbnd, nks)
       vmef(:,:,:,:)            !  velcity matrix elements on the fine   mesh (ipol, nbnd, nbnd, nks)
  REAL(KIND=DP), ALLOCATABLE ::&
       xk_all(:,:),            &!  full k point grid, coarse (3, nkstot)
       et_all(:,:),            &!  full eigenvalue list, coarse (nbnd, nkstot)
       et_ks(:,:),            &!  lda eigenvalues
       et_mb(:,:),             &!  gw eigenvalues
       xkq(:,:),               &!  local k+q grid, coarse (3, nks)
       etq(:,:),               &!  eigenvalues of k+q wavefunctions
       xkf(:,:),               &!  fine k point grid (3, nksf)
       wkf(:),                 &!  weights on the fine grid (nksf)
       xqf(:,:),               &!  fine q point grid 
       wqf(:),                 &!  weights on the fine q grid 
       etf(:,:),               &!  interpolated eigenvalues (nbnd, nksf)
       etf_ks(:,:),            &!  interpolated eigenvalues (nbnd, nksf) KS eigenvalues in the case of eig_read
       wf(:,:),                &!  interpolated eigenfrequencies 
       wslen(:),               &!  length of the wigner seitz points in units of alat
       etfq(:,:,:),            &!  energies
       lambda_all(:,:,:),      &
       lambda_v_all(:,:,:),    &
       jdos(:),                &
       spectra(:,:,:,:,:,:)     !  dipole absorption spectra, polarizations, nomega, nsmear, dme/vme, absorption/emission
  INTEGER :: &
       nksf,                   &!  number of k points in the pool (fine grid)
       nksqf,                  &!  number of k blocks in the pool (fine grid)
       nkstotf,                &!  total number of k points (fine grid)
       nxqf,                   &!  total number of q points (fine grid)
       nrr,                    &!  number of wigner-seitz points (elec interp only)
       nrr_k,                  &!  number of wigner-seitz points for electrons
       nrr_q,                  &!  number of wigner-seitz points for phonons
       ibndmin,                &!  band bounds for slimming down electron-phonon matrix 
       ibndmax                  !
  INTEGER, ALLOCATABLE ::      & 
       irvec(:,:),             &!  crys coordinates of wigner-seitz vectors (both elec and phon)
       ndegen(:),              &!  corresponding degeneragy, electrons (old version)
       ndegen_k(:),            &!  corresponding degeneragy, electrons
       ndegen_q(:)              !  corresponding degeneragy, phonons
  INTEGER, allocatable ::       &
       shift (:),              &!  for every k+q, index of the G0 which folds k+q into k+q+G0 of the first BZ
       gmap(:)                  !  the map G -> G-G_0 in the large (density) G vectors set, for every G_0
  END MODULE el_phon
