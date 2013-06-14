  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !--------------------------------------------------------------------------
  SUBROUTINE dmebloch2wan ( nbnd, nbndsub, nks, nkstot, nkbl, lgamma, dme, xk, cu, &
     nrr, irvec, wslen )
  !--------------------------------------------------------------------------
  !
  !  From the Hamiltonian in Bloch representationi (coarse mesh), 
  !  find the corresponding Hamiltonian in Wannier representation 
  !
  !--------------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, celldm
  USE control_flags, ONLY : iverbosity
  USE el_phon,       ONLY : cdmew
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : mpime,inter_pool_comm
  USE mp,            ONLY : mp_barrier,mp_sum
#endif
  implicit none
  !
  !  input variables
  !
  integer :: nbnd, nbndsub, nks, nkstot, nkbl, nrr, irvec (3, nrr), ipol
  ! number of bands 
  ! number of bands in the optimal subspace 
  ! number of kpoints
  ! number of kpoint blocks, in the pool
  ! number of kpoint blocks, total 
  ! number of WS points and coordinates
  logical :: lgamma
  ! true if we are working with q=0
  real(kind=DP) ::  xk (3, nks), wslen (nrr) 
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  complex(kind=DP) :: dme (3,nbnd, nbnd,nks)  
  ! Dipole matrix elements on coarse mesh
  complex(kind=DP) :: cu (nbnd, nbndsub, nks)
  ! rotation matrix from wannier code
  !
  !  output variables
  !
  ! work variables 
  !
  complex(kind=DP) :: cps(3, nbndsub, nbndsub, nks)
  ! Hamiltonian in smooth Bloch basis, coarse mesh 
  real(kind=DP), PARAMETER :: bohr2ang = 0.5291772108, twopi = 6.28318530717959
  complex(kind=DP), PARAMETER :: ci = (0.d0,1.d0), czero = (0.d0, 0.d0), &
       cone = (1.d0, 0.d0)
  integer :: ik, ir
  real(kind=DP) :: rdotk, tmp
  complex(kind=DP) :: cfac
  complex(kind=DP) :: dme_utmp(nbnd,nbndsub)
  !
  !----------------------------------------------------------
  !    STEP 1: rotation to optimally smooth Bloch states
  !----------------------------------------------------------
  !
  CALL start_clock ( 'Dipole: step 1' )
  !
  !  p~ (k) = U(k)^\dagger * p(k) * U(k)
  !  p~ (k) is cps( ipol, nbndsub, nbndsub, ik )
  !
  DO ik = 1, nks
     DO ipol = 1, 3
        ! copied from ephbloch2wane.  produce equivalent results
        !
        dme_utmp(:,:) = matmul( dme(ipol,:,:,ik), cu(:,:,ik) )
        cps (ipol, :,:, ik) = matmul ( conjg(transpose( cu(:,:,ik))), dme_utmp (:,:) )
        !
    ENDDO
  ENDDO
  !
  !
  !----------------------------------------------------------
  !  STEP 2: Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  !  p (R) = (1/nk) sum_k e^{-ikR} p~ (k)
  !  cdmew (ipol, nbndsub, nbndsub, ir) is p (R)
  !
  CALL start_clock ( 'Dipole: step 2' )
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nks, xk, at, -1)
  !
  cdmew ( :, :, :, :) = czero 
  ! 
  DO ir = 1, nrr
    !
    DO ik = 1, nks
       !
       rdotk = twopi * dot_product( xk ( :, ik), float(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / float(nkbl)
       DO ipol = 1, 3
          cdmew ( ipol, :, :, ir ) = cdmew ( ipol, :, :, ir ) + cfac * cps ( ipol, :, :, ik )
       ENDDO
       !
    ENDDO
    !
  ENDDO
#ifdef __PARA
  CALL mp_sum(cdmew,inter_pool_comm) 
#endif
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nks, xk, bg, 1)
  !
  !
  !  check spatial decay of Dipole in Wannier basis
  !  the unit in r-space is angstrom
  !
#ifdef __PARA
    IF (mpime.eq.ionode_id) then
#endif
       OPEN(unit=300,file='decay.P')
       WRITE(300, '(/3x,a/)') '#Spatial decay of Dipole in Wannier basis'
       DO ir = 1, nrr
          !
          tmp =  maxval ( abs( cdmew (:, :,:,ir)) )
          WRITE(300, *) wslen(ir) * celldm (1) * bohr2ang, tmp
          !
       ENDDO
       close(300)
#ifdef __PARA
    ENDIF
    CALL mp_barrier()
#endif
  !
  CALL stop_clock ( 'Dipole: step 2' )
  !
  END SUBROUTINE dmebloch2wan

