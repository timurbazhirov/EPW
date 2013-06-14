  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !---------------------------------------------------------------------------
  subroutine rotate_epmat ( cz1, cz2, xq, iq)
  !---------------------------------------------------------------------------
  !
  ! 1). rotate the electron-phonon matrix from the cartesian representation
  !    of the first qpoint of the star to the eigenmode representation 
  !    (using cz1).
  ! 
  ! 2). rotate the electron-phonon matrix from the eigenmode representation
  !     to the cartesian representation of the qpoint iq (with cz2).
  !
  !
  !--------------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, only : DP
  use io_global, only : stdout
  use el_phon, only : epmatq
  USE pwcom, ONLY : amconv
  use phcom, only : nmodes
  use control_flags, only : iverbosity
  use pwcom, only : nbnd, nks, nr1, nr2, nr3, at, bg
  USE ions_base, ONLY : amass, tau, nat, ntyp => nsp, ityp
  USE mp_global, ONLY : my_pool_id
  implicit none
  !
  integer :: iq
  !  the current qpoint
  real(kind=DP) :: xq(3)
  !  the rotated q vector
  complex(kind=DP) :: cz1( nmodes, nmodes), cz2(nmodes, nmodes)
  !  the eigenvectors for the first q in the star
  !  the rotated eigenvectors, for the current q in the star
  !
  complex(kind=DP), parameter :: cone = (1.d0, 0.d0), czero = (0.d0, 0.d0)
  real(kind=DP), parameter :: zero = 0.d0, twopi = 6.28318530717959, &
      rydcm1 = 13.6058d0 * 8065.5d0
  !
  ! work variables 
  !
  complex(kind=DP) :: eptmp( nmodes)
  integer :: mu, na, ik, ibnd, jbnd
  real(kind=DP) :: massfac
  complex(kind=DP) :: cz_tmp(nmodes,nmodes)
  !
  ! the mass factors: 
  !  1/sqrt(M) for the  direct transform
  !  sqrt(M)   for the inverse transform 
  !
  ! if we set cz1 = cz2 here and we calculate below
  ! cz1 * cz2 we find the identity
  !
  DO mu = 1, nmodes
    na = (mu - 1) / 3 + 1
    massfac = sqrt(amass(ityp(na)))
    cz1 (mu, :) = cz1 (mu, :) / massfac
    cz2 (mu, :) = cz2 (mu, :) * massfac
  ENDDO
  !
  ! the inverse transform also requires the hermitian conjugate
  !
  cz_tmp = conjg ( transpose ( cz2 ) )
  cz2 = cz_tmp
  ! 
  !  ep_mode (j) = cfac * sum_i ep_cart(i) * u(i,j)
  !
  DO ibnd = 1, nbnd
   DO jbnd = 1, nbnd
    DO ik = 1, nks
       !
       !  bring e-p matrix from the cartesian representation of the
       !  first q in the star to the corresponding eigenmode representation
       !
       CALL zgemv ('t', nmodes, nmodes, cone, cz1, nmodes,  &
          epmatq (ibnd, jbnd, ik, :, iq), 1, czero, eptmp, 1 )
       !
       ! rotate epmat in the cartesian representation for this q in the star
       !
       CALL zgemv ('t', nmodes, nmodes, cone, cz2, nmodes, &
          eptmp, 1, czero, epmatq (ibnd, jbnd, ik, :, iq), 1 )
       !
    ENDDO
   ENDDO
  ENDDO
  !
  end subroutine rotate_epmat

