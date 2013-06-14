  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dynwan2bloch ( nbnd, nrr, irvec, ndegen, xq, cuf, eig)
  !--------------------------------------------------------------------------
  !
  !
  !  WARNING: this SUBROUTINE is identical to hamwan2bloch.f90, except
  !           that here rdw is a real array, not a complex one. This is
  !           required to obtain proper phonon dispersion interpolation
  !           and corresponds to the reality of the interatomic force
  !           constants
  !
  ! -------------------------------------------------------------------------
  !
  !  From the Hamiltonian in Wannier representation, find the corresponding
  !  Hamiltonian in Bloch representation for a given k point
  !
  !  input  : number of bands nbnd
  !           number of WS vectors, coordinates and degeneracy
  !           Hamiltonian in Wannier representation chw(nbnd, nbnd, nrr)
  !           qpoint coordinate xq(3)
  !
  !  output : rotation matrix cuf (nbnd, nbnd)
  !           interpolated hamiltonian eigenvalues eig(nbnd)
  !
  !
  !--------------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE ions_base, ONLY : amass, nat, ntyp => nsp, ityp
  USE el_phon,   ONLY : rdw
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nrr, irvec (3, nrr), ndegen (nrr)
  ! number of bands (possibly of the optimal subspace)
  ! kpoint number for the interpolation
  ! record length and unit for direct write of rotation matrix
  ! number of WS points, crystal coordinates, degeneracy
  !
  ! Hamiltonian in wannier basis
  !
  real(kind=DP) :: xq (3)
  ! kpoint coordinates for the interpolation
  !
  ! output variables
  !
  real(kind=DP) :: eig (nbnd)
  ! interpolated hamiltonian eigenvalues for this kpoint
  complex(kind=DP) :: cuf(nbnd, nbnd)
  ! Rotation matrix, fine mesh
  !
  Real(kind=DP), PARAMETER :: twopi = 6.28318530717959
  complex(kind=DP), PARAMETER :: ci = (0.d0,1.d0), czero = (0.d0, 0.d0)
  !
  ! variables for lapack ZHPEVX
  !
  integer :: neig, info, ifail( nbnd ), iwork( 5*nbnd )
  real(kind=DP) :: w( nbnd ), rwork( 7*nbnd )
  complex(kind=DP) :: champ( nbnd*(nbnd+1)/2 ), &
    cwork( 2*nbnd ), cz( nbnd, nbnd)
  !
  ! work variables
  !
  complex(kind=DP) :: chf(nbnd, nbnd)
  ! Hamiltonian in Bloch basis, fine mesh
  integer :: ibnd, jbnd, ir, na, nb
  real(kind=DP) :: rdotk, massfac
  complex(kind=DP) :: cfac
  !
  CALL start_clock ( 'DynW2B' )
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
  chf (:,:) = czero
  !
  DO ir = 1, nrr
     !
     rdotk = twopi * dot_product( xq, float(irvec( :, ir) ))
     cfac = exp( ci*rdotk ) / float( ndegen(ir) )
     chf (:,:) = chf (:,:) + cfac * rdw (:,:, ir )
     !
  ENDDO
  !  divide by the square root of masses 
  !
  DO na = 1, nat
  DO nb = 1, nat
    massfac = 1.d0 / sqrt ( amass(ityp(na)) * amass(ityp(nb)) )

    chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) = &
        chf(3*(na-1)+1:3*na, 3*(nb-1)+1:3*nb) * massfac
    
  END DO
  END DO
  !
  !---------------------------------------------------------------------
  !  STEP 4: diagonalize smooth Hamiltonian on k points of the fine grid
  !---------------------------------------------------------------------
  !
  ! champ: complex hamiltonian packed (upper triangular part for zhpevx)
  ! after hermitian-ization
  !
  DO jbnd = 1, nbnd
   DO ibnd = 1, jbnd
      champ (ibnd + (jbnd - 1) * jbnd/2 ) = &
      ( chf ( ibnd, jbnd) + conjg ( chf ( jbnd, ibnd) ) ) / 2.d0
   ENDDO
  ENDDO
  !
  CALL zhpevx ('V', 'A', 'U', nbnd, champ , 0.0, 0.0, &
               0, 0,-1.0, neig, w, cz, nbnd, cwork, &
               rwork, iwork, ifail, info)
  !
  ! rotation matrix and Ham eigenvalues
  ! [in Ry, mind when comparing with wannier code]
  !
  cuf = cz
  eig = w
  !
  CALL stop_clock ( 'DynW2B' )
  !   
  END SUBROUTINE dynwan2bloch

