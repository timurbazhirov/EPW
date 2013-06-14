  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !--------------------------------------------------------------------------
  SUBROUTINE ephbloch2wanp ( nbnd, nmodes, xk, nq, dynq, irvec, wslen, &
    nrk, nrr, epmatwe )
  !--------------------------------------------------------------------------
  !
  !  From the EP Matrix in Electron Bloch representation (coarse mesh), 
  !  find the corresponding matrix in Phonon Wannier representation 
  !
  !--------------------------------------------------------------------------------
  !
#include "f_defs.h"
  USE io_global,     ONLY : stdout
  USE kinds,         ONLY : DP
  USE pwcom,         ONLY : at, bg, celldm
  USE epwcom,        ONLY : iunepmatwe, iunepmatwp
  USE control_flags, ONLY : iverbosity
  USE el_phon,       ONLY : epmatwp
#ifdef __PARA
  USE io_global,     ONLY : ionode_id
  USE mp_global,     ONLY : mpime
  USE mp, ONLY : mp_barrier
#endif
  implicit none
  !
  !  input variables - note irvec is dimensioned with nrr_k 
  !                    (which is assumed to be larger than nrr_q)
  !
  integer :: nbnd, nrk, nmodes, nq, nrr, irvec (3, nrk)
  ! number of electronic bands
  ! number of electronic WS points
  ! number of branches
  ! number of qpoints
  ! number of WS points and coordinates
  complex(kind=DP) :: dynq (nmodes, nmodes, nq), epmatwe (nbnd, nbnd, nrk, nmodes)
  ! dynamical matrix in bloch representation (Cartesian coordinates)
  ! EP matrix in electron-wannier representation and phonon bloch representation
  !   (Cartesian coordinates)
  real(kind=DP) :: xk (3, nq), wslen (nrr) 
  ! kpoint coordinates (cartesian in units of 2piba)
  ! WS vectors length (alat units)
  !
  !  output variables
  !
  ! EP matrix in electron-wannier representation and phonon-wannier representation
  !
  !
  ! work variables 
  !
  real(kind=DP), PARAMETER :: bohr2ang = 0.5291772108, twopi = 6.28318530717959
  complex(kind=DP), PARAMETER :: ci = (0.d0,1.d0), czero = (0.d0, 0.d0), cone = (1.d0,0.d0)
  integer :: ik, ir, ire
  real(kind=DP) :: rdotk, tmp, rvec1(3), rvec2(3), len1, len2
  complex(kind=DP) :: cfac
  real(kind=DP), PARAMETER :: rydcm1 = 13.6058d0 * 8065.5d0
  !

  !
  !----------------------------------------------------------
  !  Fourier transform to go into Wannier basis
  !----------------------------------------------------------
  !
  !  D (R) = (1/nk) sum_k e^{-ikR} D (k)
  !
  ! bring xk in crystal coordinates
  !
  CALL cryst_to_cart (nq, xk, at, -1)
  !
  DO ir = 1, nrr
    !
    epmatwp = czero
    !
    DO ik = 1, nq
       !
       ! direct read of epmatwe for this iq 
       CALL rwepmatw ( epmatwe, nbnd, nrk, nmodes, ik, iunepmatwe, -1)
       !
       rdotk = twopi * dot_product( xk ( :, ik), float(irvec( :, ir) ))
       cfac = exp( -ci*rdotk ) / float(nq)
       epmatwp = epmatwp + cfac * epmatwe
       !
    ENDDO
    !
    ! direct write of epmatwp for this ir 
    CALL rwepmatwp ( nbnd, nrk, nmodes, ir, iunepmatwp, +1)
    !
    !  check spatial decay of e-p matrix elements in wannier basis - electrons
    !  + phonons
    !
    !  we plot: R_e, R_p, max_{m,n,nu} |g(m,n,nu;R_e,R_p)|
    !
#ifdef __PARA
    IF (mpime.eq.ionode_id) THEN
#endif
      IF (ir.eq.1) open(unit=303,file='decay.epmat_wanep',status='unknown')
      DO ire = 1, nrk
        !
        rvec1 = float(irvec(1,ire))*at(:,1) + &
                float(irvec(2,ire))*at(:,2) + &
                float(irvec(3,ire))*at(:,3)
        rvec2 = float(irvec(1,ir))*at(:,1) + &
                float(irvec(2,ir))*at(:,2) + &
                float(irvec(3,ir))*at(:,3)
        len1 = sqrt(rvec1(1)**2.d0+rvec1(2)**2.d0+rvec1(3)**2.d0)
        len2 = sqrt(rvec2(1)**2.d0+rvec2(2)**2.d0+rvec2(3)**2.d0)
        tmp =  maxval ( abs( epmatwp (:, :, ire, :) ) )
        !
        ! rvec1 : electron-electron0 distance
        ! rvec2 : phonon - electron0 distance
        !
        WRITE(303, '(5f15.10)') len1 * celldm (1) * bohr2ang, &
                                len2 * celldm (1) * bohr2ang, tmp
      ENDDO
      IF (ir.eq.nrr) close(303)
#ifdef __PARA
    ENDIF
#endif
    !
  ENDDO
  !
  ! bring xk back into cart coord
  !
  CALL cryst_to_cart (nq, xk, bg, 1)
  !
  END SUBROUTINE ephbloch2wanp
