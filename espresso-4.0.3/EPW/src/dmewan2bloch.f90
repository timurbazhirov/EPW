
  ! 
  !--------------------------------------------------------------------------
  SUBROUTINE dmewan2bloch ( nbnd, nrr, irvec, ndegen, xk, cuf, dmef, etf, etf_ks)
  !--------------------------------------------------------------------------
  !
  !  From the Hamiltonian in Wannier representation, find the corresponding
  !  Hamiltonian in Bloch representation for a given k point
  !  
  !  input  : number of bands nbnd
  !           number of WS vectors, coordinates and degeneracy 
  !           Hamiltonian in Wannier representation chw(nbnd, nbnd, nrr)
  !           kpoint coordinate xk(3)
  !
  !  output : interpolated dipole matrix elements (dmef)
  !
  !  JN, EK 09/2010
  !
  !--------------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE mp_global, ONLY : my_pool_id
  USE el_phon,   ONLY : cdmew
  USE epwcom,    ONLY : eig_read
  !
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nrr, irvec (3, nrr), ndegen (nrr), ipol, ibnd, jbnd
  ! number of bands (possibly of the optimal subspace)
  ! kpoint number for the interpolation
  ! record length and unit for direct write of rotation matrix
  ! number of WS points, crystal coordinates, degeneracy
  !
  ! Hamiltonian in wannier basis
  !
  real(kind=DP) :: xk (3), etf(nbnd), etf_ks(nbnd)
  ! kpoint coordinates for the interpolation
  !
  ! output variables
  !  
  complex(kind=DP) :: dmef (3,nbnd,nbnd)
  ! interpolated hamiltonian eigenvalues for this kpoint 
  complex(kind=DP) :: cuf(nbnd, nbnd)
  ! Rotation matrix, fine mesh 
  !
  real(kind=DP), PARAMETER :: twopi = 6.28318530717959
  complex(kind=DP), PARAMETER :: ci = (0.d0,1.d0), czero = (0.d0, 0.d0), &
       cone = (1.d0, 0.d0)
  !
  complex(kind=DP) :: cdmef(3,nbnd, nbnd), cdmef_tmp(nbnd, nbnd)
  ! dipole matrix elements in Bloch basis, fine mesh
  integer :: ir
  real(kind=DP) :: rdotk
  complex(kind=DP) :: cfac
  !
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform to fine k and k+q meshes
  !----------------------------------------------------------
  !
  !  H~_k'   = sum_R 1/ndegen(R) e^{-ik'R    } H_k(R)
  !  H~_k'+q = sum_R 1/ndegen(R) e^{-i(k'+q)R} H_k+q(R)
  !
  !  H~_k   is chf ( nbnd, nbnd, 2*ik-1 )
  !  H~_k+q is chf ( nbnd, nbnd, 2*ik   )
  !
  cdmef (:,:,:) = czero
  !  
  DO ir = 1, nrr
     !
     rdotk = twopi * dot_product( xk, float(irvec( :, ir) ))
     cfac = exp( ci*rdotk ) / float( ndegen(ir) )
     DO ipol = 1,3
        cdmef (ipol,:,:) = cdmef (ipol,:,:) + cfac * cdmew (ipol, :,:, ir )
     ENDDO
     !
  ENDDO
  !
  ! pmn(k) = U p~ U^dagger
  ! cuf,  passed from hamwan2bloch.
  !
  DO ipol = 1,3
     cdmef_tmp(:,:) = matmul( cdmef(ipol,:,:) ,  conjg(transpose(cuf(:,:)))  )
     dmef(ipol, :,:) = matmul(cuf(:,:) , cdmef_tmp(:,:) )
  ENDDO
  !
  ! Satisfy
  ! Phys. Rev. B 62, 4927–4944 (2000) , Eq. (30)
  !
  IF (eig_read) THEN
     DO ibnd = 1, nbnd
        DO jbnd = 1, nbnd
           IF (abs(etf_ks(ibnd) - etf_ks(jbnd)) .gt. 1.d-4) THEN
              dmef(:, ibnd, jbnd) = dmef(:,ibnd, jbnd) * &
                   ( etf(ibnd)    - etf(jbnd) )/ &
                   ( etf_ks(ibnd) - etf_ks(jbnd) )
           ENDIF
        ENDDO
     ENDDO
  ENDIF
  !
  END SUBROUTINE dmewan2bloch
  !
