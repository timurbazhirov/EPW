  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! 
  !---------------------------------------------------------------------------
  subroutine ephwan2blochp ( nmodes, xxq, irvec, ndegen, nrr_q, cuf, epmatf, nbnd, nrr_k )
  !---------------------------------------------------------------------------
  !
  ! even though this is for phonons, I use the same notations
  ! adopted for the electronic case (nmodes->nmodes etc)
  !
  !
#include "f_defs.h"
  USE kinds, only : DP
  USE epwcom, only : iunepmatwp, epf_mem, parallel_k, parallel_q
  use el_phon, only : epmatw17
  USE mp_global, ONLY : my_pool_id, intra_pool_comm, &
       inter_pool_comm,  npool
  use mp,            only : mp_bcast, mp_sum
  USE io_global,     only : ionode_id
  implicit none
  !
  !  input variables
  !
  integer :: nmodes, nrr_q, irvec ( 3, nrr_q), ndegen (nrr_q), nbnd, nrr_k
  ! number of bands (possibly in tyhe optimal subspace)
  ! number of WS points
  ! coordinates of WS points
  ! degeneracy of WS points
  ! n of bands
  ! n of electronic WS points
  complex(kind=DP) :: epmatw ( nbnd, nbnd, nrr_k, nmodes), cuf (nmodes, nmodes)
  ! e-p matrix in Wanner representation
  ! rotation matrix U(k)
  real(kind=DP) :: xxq(3)
  ! kpoint for the interpolation (WARNING: this must be in crystal coord!)
  !
  !  output variables
  !
  complex(kind=DP) :: epmatf (nbnd, nbnd, nrr_k, nmodes)
  ! e-p matrix in Bloch representation, fine grid
  !
  ! work variables 
  !
  integer :: ibnd, jbnd, ir, ire, ir_start, ir_stop, nrst, rest
  real(kind=DP), parameter :: twopi = 6.28318530717959
  real(kind=DP) :: rdotk
  complex(kind=DP) :: cfac, eptmp( nbnd, nbnd, nrr_k, nmodes)
  complex(kind=DP), parameter :: ci = (0.d0,1.d0), &
     czero = (0.d0, 0.d0), cone = (1.d0, 0.d0)
  !
  CALL start_clock('ephW2Bp')
  !----------------------------------------------------------
  !  STEP 3: inverse Fourier transform of g to fine k mesh
  !----------------------------------------------------------
  !
  !  g~ (k') = sum_R 1/ndegen(R) e^{-ik'R} g (R)
  !
  !  g~(k') is epmatf (nmodes, nmodes, ik )
  !  every pool works with its own subset of k points on the fine grid
  !
  IF (parallel_k) THEN
     CALL para_bounds(ir_start, ir_stop, nrr_q)
  ELSEIF (parallel_q) THEN
     ir_start = 1
     ir_stop  = nrr_q
  ELSE 
     CALL errore ('ephwan2blochp', 'Problem with parallel_k/q scheme', nrr_q)
  ENDIF
  !
  eptmp = czero
  !
  DO ir = ir_start, ir_stop
     !   
     !  direct read of epmatwp for this ir 
     IF (.not.epf_mem) call rwepmatw ( epmatw, nbnd, nrr_k, nmodes, ir, iunepmatwp, -1)
     !
     ! note xxq is assumed to be already in cryst coord
     !
     rdotk = twopi * dot_product ( xxq, float(irvec(:, ir)) )
     cfac = exp( ci*rdotk ) / float( ndegen(ir) )
     !
     DO ibnd = 1, nbnd
      DO jbnd = 1, nbnd
       DO ire = 1, nrr_k
        !
          IF (epf_mem) then
             eptmp (ibnd,jbnd,ire,:) = eptmp (ibnd,jbnd,ire,:) + cfac * epmatw17 ( ibnd, jbnd, ire, ir, :)
          ELSE
             eptmp (ibnd,jbnd,ire,:) = eptmp (ibnd,jbnd,ire,:) + cfac * epmatw ( ibnd, jbnd, ire, :)
          ENDIF
        !
       ENDDO
      ENDDO
     ENDDO
     !
  ENDDO
#ifdef __PARA
  IF (parallel_k) CALL mp_sum(eptmp, inter_pool_comm)
#endif
  !
  !----------------------------------------------------------
  !  STEP 4: un-rotate to Bloch space, fine grid
  !----------------------------------------------------------
  !
  ! epmatf(j) = sum_i eptmp(i) * uf(i,j)
  !
  DO ibnd = 1, nbnd
   DO jbnd = 1, nbnd
    DO ire = 1, nrr_k
     !
     CALL zgemv ('t', nmodes, nmodes, cone, cuf, nmodes, eptmp(ibnd,jbnd,ire,:), &
        1, czero, epmatf(ibnd,jbnd,ire,:), 1 )
     !
    ENDDO
   ENDDO
  ENDDO
  !
  CALL stop_clock('ephW2Bp')
  !
  end subroutine ephwan2blochp

