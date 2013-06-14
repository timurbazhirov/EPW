!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
! ... Common variables for conductance calculation  
!
!
!
MODULE geomcell_cond
  USE kinds, only : DP
  !
  SAVE
  !
  INTEGER  ::  &
     nrx,      &   !  number of mesh points in the x direction
     nry,      &   !              -||-                 y direction
     nrzl,     &   !  number of slabsh in the z direction for the left lead 
     nrzs,     &   !              -||-                    for the scatt. region
     nrzr,     &   !              -||-                    for the right lead 
     nrzpl,    &   !  number of slabs per CPU for the left lead  
     nrzps,    &   !              -||-        for the scatt. region 
     nrzpr,    &   !              -||-        for the right lead 
     ngper,    &   !  number of perpendicular G vectors
     ngpsh,    &   !  number of shells for G
     nkpts,    &   !  number of kpts in the perpendicular direction
     n2d,      &   !  dimension of reduced vector space in XY
     nz1           !  number of subslabs in the slab
  INTEGER   :: &
            nk1ts, nk2ts, &    ! k-point mesh dimensions
            k1ts, k2ts         ! k-point mesh shift
  INTEGER, ALLOCATABLE :: &
     ninsh(:),   & !  number of G in shell
     nl_2d(:)      !  correspondence G list 2D fft_mesh
  REAL(DP)    :: &
     bdl,      &   !  right boundary of the left lead
     bds,      &   !        -||-     of the scatt. region 
     bdr,      &   !        -||-     of the right lead
     sarea         !  the cross section
  REAL(DP), ALLOCATABLE ::  &
     zl(:),    &   !  the division in the z direction of the left lead
     zs(:),    &   !             -||-                 of the scatt. reg.
     zr(:),    &   !             -||-                 of the right lead 
     xyk(:,:), &   !  coordinates of perpendicular k points
     wkpt(:),  &   !  the weight of k point
     gper(:,:),&   !  coordinates of perpendicular G
     gnsh(:)       !  the norm of the G shell
END MODULE geomcell_cond
!
!
MODULE orbcell_cond
  USE parameters, only : npsx 
  use radial_grids, only: ndmx
  USE kinds, only : DP
  !
  !  description of nonlocal orbitals
  SAVE
  !
  INTEGER, PARAMETER :: nbrx = 14
  !
  INTEGER     ::     &
     norbl,          &  ! number of orbitals for the left lead
     norbs,          &  !           -||-      for the scatt. region
     norbr,          &  !           -||-      for the right lead
     nocrosl,        &  ! number of crossing orbitals for left lead
     nocrosr,        &  !           -||-      for the right lead
     noinsl,         &  ! number of interior orbitals for the left lead
     noinss,         &  !           -||-      for the scatt. region
     noinsr,         &  !           -||-      for the right lead
     norbf              ! max number of needed orbitals
  INTEGER, ALLOCATABLE :: &
     tblml(:,:),       &  ! the type/beta/l/m of each orbital for the left lead
     tblms(:,:),       &  !         -||-                      for the scatt. reg.
     tblmr(:,:),       &  !         -||-                      for the right lead
     crosl(:,:),       &  ! 1 if the orbital crosses the slab - for the left lead
     cross(:,:),       &  !                -||-                 for the scatt. reg. 
     crosr(:,:)           !                -||-                 for the right lead
  REAL(DP) ::     &
     rl(ndmx,npsx),          &  ! radial mesh for the left lead
     rs(ndmx,npsx),          &  !      -||-   for the scatt. reg.
     rr(ndmx,npsx),          &  !      -||-   for the right lead
     rabl(ndmx,npsx),        &  ! log. mesh for the left lead
     rabs(ndmx,npsx),        &  !      -||- for the scatt. reg.
     rabr(ndmx,npsx),        &  !      -||- for the right lead
     betarl(ndmx,nbrx,npsx), &  ! beta functions for the left lead
     betars(ndmx,nbrx,npsx), &  !      -||-      for the scatt. reg.
     betarr(ndmx,nbrx,npsx)     !      -||-      for the right lead
  REAL(DP), ALLOCATABLE :: &
     taunewl(:,:),  &  ! center of each orbital and its radius - left lead
     taunews(:,:),  &  !                 -||-                  - scatt. reg. 
     taunewr(:,:),  &  !                 -||-                  - right lead
     zpseul(:,:,:), &  ! coefficients of nonlocal pseudopotential - left lead
     zpseus(:,:,:), &  !                 -||-                     - scatt. reg.
     zpseur(:,:,:)     !                 -||-                     - right lead
  COMPLEX(DP), ALLOCATABLE :: &
     zpseul_nc(:,:,:,:), &! coefficients of nonlocal PP (nc case) - left lead
     zpseus_nc(:,:,:,:), &!              -||-                     - scatt. reg. 
     zpseur_nc(:,:,:,:)   !              -||-                     - right lead 
END MODULE orbcell_cond
!
!
MODULE eigen_cond
  USE kinds, only : DP
  !
  !   Eigenvalue equation for local potential
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE  :: &
     vppotl(:,:,:,:),  & ! Fourier comp. of local potential in each slab - left lead
     vppots(:,:,:,:),  & !             -||-                              - scatt. reg. 
     vppotr(:,:,:,:),  & !             -||-                              - right lead
     psiperl(:,:,:),   & ! eigenvectors in each slab - left lead
     psipers(:,:,:),   & !             -||-          - scatt. reg. 
     psiperr(:,:,:),   & !             -||-          - right lead
     zkl(:,:),         & ! the k for each eigenvalue (computed through zkr) - left lead
     zks(:,:),         & !             -||-                                 - scatt. reg.
     zkr(:,:),         & !             -||-                                 - right lead
     newbg(:,:)          ! reduced basis set --> exp(G)
  REAL(DP), ALLOCATABLE ::  &
     zkrl(:,:),        & ! 2d eigenvalues - left lead
     zkrs(:,:),        & !     -||-       - scatt. reg. 
     zkrr(:,:)           !     -||-       - right lead 
END MODULE eigen_cond
!
!
MODULE control_cond
  USE kinds, only : DP
  !
  !   control of the run
  SAVE
  !
  INTEGER     :: &
     orbj_in, orbj_fin, &
     ikind,       &    ! the kind of calculation
     nenergy,     &    ! number of energies computed
     iofspin           ! spin index for calculation
  REAL(DP)  :: &
     efl,           &  ! the Ef of the left lead
     efs,           &  ! the Ef of the scatt. reg.
     efr,           &  ! the Ef of the right lead 
     energy0,       &  ! initial energy
     eryd,          &  ! the current energy in Ry 
     denergy,       &  ! delta of energy
     ecut2d,        &  ! 2D cutoff
     ewind,         &  ! the window above energy for 2D computation
     delgep,        &  ! infinitesimal for GEP
     epsproj,       &  ! accuracy of n2d reduction
     cutplot           ! cutoff of Im(k) for CB plotting
  REAL(DP), ALLOCATABLE  :: &
     earr(:),       &  ! energy array
     tran_tot(:)       ! transmission array
  LOGICAL        :: &
     lorb,          &  ! if .t. calculate the states in all slabs
     lwrite_loc,    &  ! if .t. save eigenproblem result on fil_loc
     lread_loc,     &  ! if .t. read eigenproblem result from fil_loc
     lwrite_cond,   &  ! if .t. save variables needed for pwcond
     lread_cond,    &  ! if .t. read variables needed for pwcond
     llocal            ! if .t. the local implementation
END MODULE control_cond
!
!
MODULE scattnl_cond
  USE kinds
  !
  ! ... The variables computed by scatter_forw
  !
  SAVE
  !
  COMPLEX(DP), ALLOCATABLE :: &
       fun0(:,:),      &!  local fun. on         left boundary
       fun1(:,:),      &!       --              right boundary
       fund0(:,:),     &!  local fun.' on        left boundary
       fund1(:,:),     &!       --              right boundary
       funl0(:,:),     &!  nonloc. fun. on       left boundary
       funl1(:,:),     &!       --              right boundary
       fundl0(:,:),    &!  nonlocal fun.' on     left boundary
       fundl1(:,:),    &!       --              right boundary
       funz0(:,:,:),   &!  local+nonlocal fun. on all slabs
       intw1(:,:),     &!  integrals with beta-fun. of loc. fun.
       intw2(:,:)       !       --                   nonloc fun. 
  !
END MODULE scattnl_cond                        
!
!
MODULE cb_cond
  USE kinds
  !
  ! ... Some variables of CBS for the leads needed for matching
  !
  SAVE
  !
  INTEGER :: &
       nchanl,         &!  number of prop. channels in the left lead
       nchanr           !            -- || --              right lead     
  COMPLEX(DP), ALLOCATABLE :: &
       kvall(:),       &!  k           for the left lead 
       kfunl(:,:),     &!  phi_k(z=d)  for the left lead
       kfundl(:,:),    &!  phi_k'(z=d) for the left lead
       kintl(:,:),     &!  integral of phi_k with beta-fun.
       kcoefl(:,:),    &!  coeff. of phi_k over nonloc. fun.
       kvalr(:),       &!  k           for the right lead
       kfunr(:,:),     &!  phi_k(z=0)  for the right lead
       kfundr(:,:),    &!  phi_k'(z=0) for the right lead
       kintr(:,:),     &!  integral of phi_k with beta-fun.
       kcoefr(:,:)      !  coeff. of phi_k over nonloc. fun.
  !
END MODULE cb_cond
!
!
MODULE cond
  use geomcell_cond
  USE orbcell_cond
  USE eigen_cond
  USE control_cond
  USE scattnl_cond
  USE cb_cond
END MODULE cond               

