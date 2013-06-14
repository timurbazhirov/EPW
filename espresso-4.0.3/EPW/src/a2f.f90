  !
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino 
  !                                                                           
  ! This file is distributed under the terms of the GNU General Public       
  ! License. See the file `LICENSE' in the root directory of the          
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .        
  !
  !-----------------------------------------------------------------------
  SUBROUTINE eliashberg_a2f( lambda, lambda_v )
  !-----------------------------------------------------------------------
  !
  !  Compute the Eliasberg spectral function
  !  in the Migdal approximation. 
  !  
  !  If the q-points are not on a uniform grid (i.e. a line)
  !  the function will not be correct
  !  
  !  02/2009 works in serial on ionode at the moment.  can be parallelized
  !  03/2009 added transport spectral function -- this involves a v_k dot v_kq term 
  !          in the quantities coming from selfen_phon.f90.  Not fully implemented  
  !  10/2009 the code is transitioning to 'on-the-fly' phonon selfenergies
  !          and this routine is not currently functional
  !  
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE phcom,     ONLY : nmodes
  USE epwcom,    ONLY : degaussq, delta_qsmear, nqsmear, nqstep, eps_acustic
  USE el_phon,   ONLY : nxqf, wf, wqf
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id, npool, mpime
#endif
  USE io_global, ONLY : stdout, ionode_id
  USE io_files,  ONLY : find_free_unit, prefix
  implicit none
  !
  real(kind=DP), PARAMETER :: ryd2mev = 13605.8, ryd2ev = 13.6058, two = 2.d0, zero = 0.d0,  &
     pi = 3.14159265358979, twopi = 6.28318530717958
  integer :: imode, iq, iw, ismear, nqx, iua2ffil, iua2ftrfil, iufildos
  real(kind=DP) :: weight
  real(kind=DP) :: lambda(nmodes, nxqf), lambda_v(nmodes, nxqf), lambda_tot
  real(kind=DP) :: iomega, sigma, a2F_tmp, a2F_tr_tmp, om_max, dw, w0, l, l_tr
  real(kind=DP), allocatable :: a2F(:,:), a2F_tr(:,:), l_a2F(:), dosph(:,:)
  real(kind=DP), external :: w0gauss
  !
  !
  CALL start_clock('a2F')
#ifdef __PARA
  IF (mpime.eq.ionode_id) &
#endif
  !
  iua2ffil   = find_free_unit()
  OPEN (unit = iua2ffil, file = TRIM(prefix)//".a2f", form = 'formatted')
  iua2ftrfil = find_free_unit()
  OPEN (unit = iua2ftrfil, file = TRIM(prefix)//".a2f_tr", form = 'formatted')
  iufildos = find_free_unit()
  OPEN (unit = iufildos, file = TRIM(prefix)//".phdos", form = 'formatted')
  !
  WRITE(stdout,'(/5x,a)') REPEAT('=',67)
  WRITE(stdout,'(5x,"Eliashberg Spectral Function in the Migdal Approximation")') 
  WRITE(stdout,'(5x,a/)') REPEAT('=',67)
  !
  ALLOCATE( a2F(nqstep+1, nqsmear) )
  ALLOCATE( a2F_tr(nqstep+1, nqsmear) )
  ALLOCATE( dosph(nqstep+1, nqsmear) )
  ALLOCATE( l_a2F(nqsmear) )
  !
  !om_max = ( MAXVAL( wf(:,:) ) - MINVAL( wf(:,:) ) ) + 5.d0/ryd2mev
  om_max = MAXVAL( wf(:,:) ) + 1.d0/ryd2mev
  dw = om_max/float(nqstep)
  !
  lambda_tot = zero
  l_a2F(:) = zero
  a2F(:,:) = zero
  a2F_tr(:,:) = zero
  dosph(:,:) = zero
  !
  DO ismear = 1, nqsmear
     !
     sigma = degaussq + (ismear-1) * delta_qsmear
     ! go from meV to Ryd
     sigma = sigma / ryd2mev
     !
     DO iw = 1, nqstep+1  ! loop over points on the a2F(w) graph
        !
        iomega = (float(iw)-1.d0) * dw ! step through the frequncies we wish to plot
        !
        DO iq = 1, nxqf ! loop over q-points 
           !
           DO imode = 1, nmodes ! loop over modes
              !
              w0 = wf(imode,iq)
              l  = lambda(imode,iq)
              l_tr = lambda_v(imode, iq)
              IF (lambda(imode, iq) .lt. 0.d0)  l = 0.d0 ! sanity check
              IF (lambda_v(imode, iq) .lt. 0.d0)  l_tr = 0.d0 ! sanity check
              IF (wf(imode, iq) .lt. 0.d0 ) w0 = 0.d0
              IF (lambda_v(imode,iq) .lt. 0.d0 ) lambda_v(imode,iq) = 0.d0
              !
              a2F_tmp    = wqf(iq) * w0 * l / two
              a2F_tr_tmp = wqf(iq) * w0 * l_tr / two
              IF (a2F_tr_tmp .lt. 0.d0 ) CALL errore ('a2F','something wrong with a2F_tr ' ,-1)
              !
              weight = w0gauss ( (iomega - w0)/sigma, 0 ) / sigma
              a2F(iw,ismear) = a2F(iw,ismear) + a2F_tmp * weight
              a2F_tr(iw,ismear) = a2F_tr(iw,ismear) + a2F_tr_tmp * weight
              dosph(iw,ismear)  = dosph(iw,ismear) + wqf(iq) * weight
              !
           ENDDO
           !
        ENDDO
        !
        !  output a2F
        !
        IF (ismear .eq. nqsmear) WRITE (iua2ffil,   '(f12.6, 15f12.6)') iomega*ryd2mev, a2F(iw,:)
        IF (ismear .eq. nqsmear) WRITE (iua2ftrfil, '(f12.6, 15f12.6)') iomega*ryd2mev, a2F_tr(iw,:)
        IF (ismear .eq. nqsmear) WRITE (iufildos,   '(f12.6, 15f12.6)') iomega*ryd2mev, dosph(iw,:)/ryd2mev
        !
        ! do the integral 2 int (a2F(w)/w dw)
        !
        IF (iomega .gt. eps_acustic) & 
           l_a2F(ismear) = l_a2F(ismear) + two * a2F(iw,ismear) / iomega * dw
        !
     ENDDO
     !
  ENDDO
  !
  WRITE(iua2ffil,*) " #   Integrated el-ph coupling"
  WRITE(iua2ffil,'("  #    ", 15f12.6)') l_a2F(:) 
  !
  CLOSE(iua2ffil)
  CLOSE(iua2ftrfil)
  CLOSE(iufildos)
  !
  DO iq = 1, nxqf ! loop over q-points 
     DO imode = 1, nmodes ! loop over modes
        IF (lambda(imode,iq) .gt. 0.d0 .and. wf(imode,iq) .gt. eps_acustic ) & 
           lambda_tot = lambda_tot + wqf(iq) * lambda(imode,iq)
     ENDDO
  ENDDO
  WRITE (6,'(5x,a,f12.6)') "lambda : ", lambda_tot
  WRITE(6,*)
!
  DEALLOCATE(l_a2F, a2F, a2F_tr, dosph)
!
#ifdef _PARA
  ENDIF
  CALL mp_barrier
#endif
  !
  CALL stop_clock('a2F')
  CALL print_clock('a2F')
  !
  RETURN
  !
  !-----------------------------------------------------------------------
  END SUBROUTINE eliashberg_a2f
  !-----------------------------------------------------------------------
