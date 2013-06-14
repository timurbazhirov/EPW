  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !------------------------------------------------------------------------
  SUBROUTINE dynbloch2wan ( nmodes, nq, xk, dynq, nrr, irvec, wslen )
  !------------------------------------------------------------------------
  !
  !  From the Dynamical Matrix in Bloch representation (coarse mesh), 
  !  find the corresponding matrix in Wannier representation 
  !
  !  NOTA BENE: it seems to be very important that the matrix is kept real
  !  physically these are truely the interatomic force constants.
  !  If you use a complex matrix instead, you may get some spurious 
  !  oscillations when you interpolate the phonon dispersions.
  !
  !  Note also that the acoustic sum rule for the q=0 case has been imposed
  !  already in readmat_shuffle.f90
  !
  !---------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds, ONLY : DP
  USE pwcom, ONLY : at, bg, celldm
  USE control_flags, ONLY : iverbosity
  USE el_phon, ONLY : rdw
#ifdef __PARA
  USE io_global, ONLY : ionode_id
  USE mp_global, ONLY : mpime
  USE mp, ONLY : mp_barrier
#endif
  implicit none
  !
  !  input variables
  !
  integer :: nmodes, nq, nrr, irvec (3, nrr)
  ! number of branches
  ! number of qpoints
  ! number of WS points and coordinates
  complex(kind=DP) :: dynq (nmodes, nmodes, nq)
  ! dynamical matrix in bloch representation (Cartesian coordinates)
  real(kind=DP) :: xk (3, nq), wslen (nrr) 
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  !
  !  output variables
  !
  ! work variables 
  !
  real(kind=DP), PARAMETER :: bohr2ang = 0.5291772108, twopi = 6.28318530717959
  complex(kind=DP), PARAMETER :: ci = (0.d0,1.d0), czero = (0.d0, 0.d0)
  integer :: ik, ir
  real(kind=DP) :: rdotk, tmp
  complex(kind=DP) :: cfac
  !
  !----------------------------------------------------------
  !  Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  !  D (R) = (1/nk) sum_k e^{-ikR} D (k)
  !  rdw (nmodes, nmodes, ir) is D (R)
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nq, xk, at, -1)
  !
  rdw ( :, :, :) = 0.d0
  ! 
  DO ir = 1, nrr
    !
    DO ik = 1, nq
       !
       rdotk = twopi * dot_product( xk ( :, ik), float(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / float(nq)
 rdw ( :, :, ir ) = rdw ( :, :, ir ) +  cfac * dynq ( :, :, ik ) 
!     rdw ( :, :, ir ) = rdw ( :, :, ir ) + real ( cfac * dynq ( :, :, ik ) ) 
       !                                     ^^^^
       !                                   note this
       !
    ENDDO
    !
  ENDDO
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nq, xk, bg, 1)
  !
  !
  !  check spatial decay of dynamical matrix in Wannier basis
  !  the unit in r-space is angstrom, and I am plotting
  !  the matrix for the first mode only
  !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      OPEN(unit=302,file='decay.dynmat')
      WRITE(302, '(/3x,a/)') '#Spatial decay of Dynamical matrix in Wannier basis'
      DO ir = 1, nrr
        !
        tmp =  maxval ( abs( rdw (:,:,ir)) )
        WRITE(302, *) wslen(ir) * celldm (1) * bohr2ang, tmp
        !
      ENDDO
      CLOSE(302)
#ifdef __PARA
    ENDIF
    CALL mp_barrier()
#endif
  !
  END SUBROUTINE dynbloch2wan
  !-----------------------------------------------------
