  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !-----------------------------------------------------------------------
  subroutine spectral_func
  !-----------------------------------------------------------------------
  !
  !  compute the electron spectral function including the  electron-
  !  phonon interaction in the Migdal approximation. 
  !  
  !  We take the trace of the spectral function to simulate the photoemission
  !  intensity. I do not consider the c-axis average for the time being.
  !  The main approximation is constant dipole matrix element and diagonal
  !  selfenergy. The diagonality can be checked numerically. 
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !  12/2009 Might be obsolete
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, only : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  use phcom, only :lgamma, nmodes
  use epwcom, only : nbndsub, lrepmatf, iunepmatf, &
      fsthick, eptemp, ngaussw, degaussw, iuetf, wmin, wmax, nw, nbndskip
  use pwcom, only : nelec, ef, isk
  use el_phon
#ifdef __PARA
  use mp, only : mp_barrier
  use mp_global, only : me_pool
#endif
  implicit none
  !
  real(kind=DP), external :: efermig
  real(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0)
  integer :: iw, ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1), weight
  real(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq, ww, dw, &
     wgauss, dosef, dos_ef, sigmar(nbndsub, nksf, nw),   &
     sigmai(nbndsub, nksf, nw), a(nw, nksf)
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf
  real(kind=DP), allocatable :: xkf_all(:,:) , etf_all(:,:), a_all(:,:)
  !

  WRITE(6,'(/5x,a)') repeat('=',67)
  WRITE(6,'(5x,"Electron Spectral Function in the Migdal Approximation")') 
  WRITE(6,'(5x,a/)') repeat('=',67)
 
  ! 
  ! energy range and spacing for spectral function
  !
  wmin = -fsthick
  wmax =  fsthick
  dw = ( wmax - wmin ) / float (nw-1)
  !
  IF ( fsthick.lt.1.d3 ) &
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
      'Fermi Surface thickness = ', fsthick, ' Ry'
  WRITE(stdout, '(/5x,a,f10.6,a)' ) &
    'Golden Rule strictly enforced with T = ',eptemp, ' Ry'
  !
  ! Fermi level and corresponding DOS
  !
  ! here we take into account that we may skip bands when we wannierize
  ! (spin-unpolarized)
  IF (nbndskip.gt.0) then
     nelec = nelec - two * nbndskip
     WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
     WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec
  ENDIF
  !
!  CALL efermig(etf, nbndsub, nksf, nelec, wkf, degaussw, ngaussw, ef0)
  ef0 = efermig(etf, nbndsub, nksf, nelec, wkf, degaussw, ngaussw,0,isk)
  DOsef = dos_ef (ngaussw, degaussw, ef0, etf, wkf, nksf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin
  DOsef = dosef / two
  !
  WRITE (6, 100) degaussw, ngaussw
  WRITE (6, 101) dosef, ef0 * ryd2ev
  !
  ! make sure q point weights add up to 1
  wqf = wqf/sum(wqf)
  !
  ! loop over all q points of the fine mesh
  ! 
  sigmar = zero
  sigmai = zero
  !
  DO iq = 1, nxqf
    !
    fermicount = 0
    !
    DO ik = 1, nksqf
      !
      IF (lgamma) then
        ikk = ik
        ikq = ik
      ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
      ENDIF
      !
      ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
      !
      nrec = (iq-1) * nksf + ikk
      CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
      nrec = (iq-1) * nksf + ikq
      CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
      !
      ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
      IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
           ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
        !
        fermicount = fermicount + 1
        !
        DO imode = 1, nmodes
          !
          ! the phonon frequency and Bose occupation
          wq = wf (imode, iq)
          wgq = wgauss( -wq/eptemp, -99)
          wgq = wgq / ( one - two * wgq )
          !
          !  we read the e-p matrix
          !
          nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
          CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
          !
          DO ibnd = 1, ibndmax-ibndmin+1
            !
            !  the energy of the electron at k
            ekk = etf (ibndmin-1+ibnd, ikk) - ef0
            !
            DO jbnd = 1, ibndmax-ibndmin+1
              !
              !  the fermi occupation for k+q
              ekq = etf (ibndmin-1+jbnd, ikq) - ef0
              wgkq = wgauss( -ekq/eptemp, -99)  
              !
              ! here we take into account the zero-point sqrt(hbar/2M\omega)
              ! with hbar = 1 and M already contained in the eigenmodes
              ! g2 is Ry^2, wkf must already account for the spin factor
              !
              ! @@@@ not sure if it is i,j or j,i
              g2 = abs(epf (jbnd, ibnd))**two / ( two * wq ) 
              !
              ! Note the +ci and -ci below, this comes from Eq. (28) of my Notes - FG
              ! 
              !  the loop over the energy variable
              !
              DO iw = 1, nw
                !
                ww = wmin + float (iw-1) * dw
                !
                weight =  two * wqf(iq) *  &
                  ( (       wgkq + wgq ) / ( ww - ( ekq - wq ) - ci * degaussw )  +  &
                    ( one - wgkq + wgq ) / ( ww - ( ekq + wq ) + ci * degaussw ) ) 
                !
                sigmar(ibndmin-1+ibnd,ikk,iw) = sigmar(ibndmin-1+ibnd,ikk,iw) + g2 * real  ( weight )
                sigmai(ibndmin-1+ibnd,ikk,iw) = sigmai(ibndmin-1+ibnd,ikk,iw) + g2 * aimag ( weight )
                ! 
              ENDDO
              !
            ENDDO
          ENDDO
          !
        ENDDO
        !
        ! endif tfermi 
      ENDIF 
      !
    ENDDO
    !
    ! no collection from pools here, q points are not distributed
    !
    ! end loop on q
  ENDDO
  !
  ! construct the trace of the spectral function (assume diagonal selfenrgy and constant
  !  matrix elements for dipole transitions)
  !
  DO ik = 1, nksqf
    !
    IF (lgamma) then
      ikk = ik
      ikq = ik
    ELSE
      ikk = 2 * ik - 1
      ikq = ikk + 1
    ENDIF
    !
    DO iw = 1, nw
      !
      ww = wmin + float (iw-1) * dw
      a(iw,ikk) = zero 
      !
      DO ibnd = 1, ibndmax-ibndmin+1
        !
        !  the energy of the electron at k
        ekk = etf (ibndmin-1+ibnd, ikk) - ef0
        !
        a(iw,ikk) = a(iw,ikk) + abs( sigmai(ibndmin-1+ibnd,ikk,iw) ) / pi / &
           ( ( ww - ekk - sigmar(ibndmin-1+ibnd,ikk,iw) )**two + ( sigmai(ibndmin-1+ibnd,ikk,iw) )**two ) 
        !
      ENDDO
      !
    ENDDO
    !
  ENDDO

  !
  WRITE( stdout, '(/5x,a,i5,a,i5/)' ) &
    'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nksqf
  !
  ! The k points are distributed amond pools: here we collect them
  !
  nksqtotf = nkstotf/2 ! even-odd for k,k+q
  !
  allocate ( xkf_all   ( 3,       nkstotf ), &
             etf_all   ( nbndsub, nkstotf ), &
             a_all ( nw, nkstotf)  )
  !
#ifdef __PARA
  !
  ! note that poolgather2 works with the doubled grid (k and k+q)
  ! therefore we need to use the sigma array with both grids, even
  ! though one of them is useless. This should be fixed be modifying
  ! poolgather2 (it's a waste of memory).
  !
  CALL poolgather2 ( 3,       nkstotf, nksf, xkf,   xkf_all  )
  CALL poolgather2 ( nbndsub, nkstotf, nksf, etf,   etf_all  )
  CALL poolgather2 ( nw, nkstotf, nksf, a, a_all)
  !
#else
  !
  xkf_all = xkf
  etf_all = etf
  a_all = a
  !
#endif
  !
  DO ik = 1, nksqtotf
    !
    IF (lgamma) then
      ikk = ik
      ikq = ik
    ELSE
      ikk = 2 * ik - 1
      ikq = ikk + 1
    ENDIF
    !
    WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf_all (:,ikk)
    WRITE(stdout,'(5x,a)') repeat('-',67)
    DO iw = 1, nw
      ww = wmin + float (iw-1) * dw
      WRITE(stdout, 103) ik, ryd2ev * ww, a_all(iw,ikk) / ryd2ev
#ifdef __PARA

      IF (me_pool == 0) & 
#endif
      WRITE(903,'(2x,i7,2x,f10.5,2x,e12.5)') ik, ryd2ev * ww, a_all(iw,ikk) / ryd2ev
    ENDDO
    WRITE(903,*)
    WRITE(stdout,'(5x,a/)') repeat('-',67)
    !
    !
  ENDDO
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
103 format(5x,'ik = ',i7,'  w = ',f10.5,' eV   A(k,w) = ',e12.5,' eV^-1') 
  !
  end subroutine spectral_func
