  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE selfen_elec
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the electron self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the electron linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  USE phcom, ONLY : lgamma, nmodes
  USE epwcom, ONLY : nbndsub, lrepmatf, iunepmatf, &
      fsthick, eptemp, ngaussw, degaussw, iuetf,   &
      nbndskip, ecutse, parallel_k, &
      parallel_q, epf_mem, etf_mem, eig_read, eps_acustic
  USE pwcom, ONLY : nelec, ef, isk
  USE el_phon, ONLY : etf, ibndmin, ibndmax, nksf, etfq, &
      epf17, wkf, nksqf, nxqf, wf, wqf, xkf, nkstotf
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id
#endif
  implicit none
  !
  REAL(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, imode, nrec, iq, fermicount
  complex(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  REAL(kind=DP) :: g2, ekk, ekq, wq, ef0, wgq, wgkq,  &
       weight, wgauss, dosef, dos_ef, sigmar(nbndsub, nksf), &
       sigmai(nbndsub, nksf), zi(nbndsub, nksf), eptemp0
  logical :: already_skipped
  REAL(kind=DP), external :: efermig
  !
  ! variables for collecting data from all pools in parallel case 
  !
  integer :: nksqtotf
  REAL(kind=DP), allocatable :: xkf_all(:,:) , etf_all(:,:), &
                                sigmar_all (:,:), sigmai_all (:,:), zi_all (:,:)
  !
  WRITE(6,'(/5x,a)') repeat('=',67)
  WRITE(6,'(5x,"Electron (Imaginary) Self-Energy in the Migdal Approximation")') 
  WRITE(6,'(5x,a/)') repeat('=',67)
  !
  ! loop over temperatures can be introduced
  !
  eptemp0 = eptemp(1)
  !
  IF ( fsthick .lt. 1.d3 ) &
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
      'Fermi Surface thickness = ', fsthick, ' Ry'
  WRITE(stdout, '(/5x,a,e18.9,a)' ) &
    'Golden Rule strictly enforced with T = ',eptemp0, ' Ry'
  !
  ! here we take into account that we may skip bands when we wannierize
  ! (spin-unpolarized)
  already_skipped = .false.
  IF ( nbndskip .gt. 0 ) THEN
     IF ( .not. already_skipped ) THEN 
        nelec = nelec - two * nbndskip
        already_skipped = .true.
        WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
        WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec
     ENDIF
  ENDIF
  !
  ! Fermi level and corresponding DOS
  !
  ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw,ngaussw,0,isk)
  !
  !   if 'fine' Fermi level differs by more than 250 meV, there is probably
  !   something wrong with the wannier functions
  IF ( abs(ef0 - ef) * ryd2eV .gt. 0.5 .and. (.not.eig_read) ) &
       CALL errore ('selfen_elec', 'Something wrong with Ef, check MLWFs', 1)
  !
  dosef = dos_ef (ngaussw, degaussw, ef0, etf, wkf, nksf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin
  dosef = dosef / two
  !
  WRITE (6, 100) degaussw, ngaussw
  WRITE (6, 101) dosef, ef0 * ryd2ev
  WRITE (6, 101) dosef, ef  * ryd2ev
  WRITE (6,'(a)') ' '
  !
  sigmar = zero
  sigmai = zero
  zi = zero
  !
  ! loop over all k points of the fine mesh
  !
  fermicount = 0 
  DO ik = 1, nksqf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     ! loop over the q-points
     !
     DO iq = 1, nxqf
        !
        ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
        !
        IF (etf_mem) THEN
           etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk, iq)
           etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq, iq)
        ELSE
           nrec = (iq-1) * nksf + ikk
           CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
           nrec = (iq-1) * nksf + ikq
           CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        ENDIF
        !
        ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
        !
        IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
             ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
           !
           fermicount = fermicount + 1
           DO imode = 1, nmodes
              !
              ! the phonon frequency and Bose occupation
              wq = wf (imode, iq)
              wgq = wgauss( -wq/eptemp0, -99)
              wgq = wgq / ( one - two * wgq )
              !
              !  we read the e-p matrix
              !
              IF (etf_mem) THEN
                 epf(:,:) = epf17 ( ik, iq, :, :, imode)
              ELSE
                 nrec = (iq-1) * nmodes * nksqf + (imode-1) * nksqf + ik
                 CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
              ENDIF
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                 !
                 !  the energy of the electron at k (relative to Ef)
                 ekk = etf (ibndmin-1+ibnd, ikk) - ef0
                 !
                 DO jbnd = 1, ibndmax-ibndmin+1
                    !
                    !  the fermi occupation for k+q
                    ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                    wgkq = wgauss( -ekq/eptemp0, -99)  
                    !
                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    IF (wq .gt. eps_acustic) THEN
                       g2 = abs(epf (jbnd, ibnd))**two / ( two * wq )
                    ELSE
                       g2 = 0.d0
                    ENDIF
                    !
                    ! There is a sign error for wq in Eq. 9 of Comp. Phys. Comm. 181, 2140 (2010). - RM
                    ! The sign was corrected according to Eq. (7.282) page 489 from Mahan's book 
                    ! (Many-Particle Physics, 3rd edition)
                    ! 
                    weight = wqf(iq) * real (                                        &
                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) )
!                   ecutse needs to be defined if it's used 
!@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                    !
                    sigmar(ibndmin-1+ibnd,ikk) = sigmar(ibndmin-1+ibnd,ikk) + g2 * weight
                    !
                    weight = wqf(iq) * aimag (                                        &
                         ( (       wgkq + wgq ) / ( ekk - ( ekq - wq ) - ci * degaussw )  +  &
                           ( one - wgkq + wgq ) / ( ekk - ( ekq + wq ) - ci * degaussw ) ) ) 
!@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                    !
                    sigmai(ibndmin-1+ibnd,ikk) = sigmai(ibndmin-1+ibnd,ikk) + g2 * weight
                    !
                    ! Z FACTOR: -\frac{\partial\Re\Sigma}{\partial\omega}
                    !
                    weight = wqf(iq) * &
                         ( (       wgkq + wgq ) * ( (ekk - ( ekq - wq ))**two - degaussw**two ) /       &
                                                  ( (ekk - ( ekq - wq ))**two + degaussw**two )**two +  &
                           ( one - wgkq + wgq ) * ( (ekk - ( ekq + wq ))**two - degaussw**two ) /       &
                                                  ( (ekk - ( ekq + wq ))**two + degaussw**two )**two )  
!@                    if ( abs(ekq-ekk) .gt. ecutse ) weight = 0.d0
                    !
                    zi(ibndmin-1+ibnd,ikk) = zi(ibndmin-1+ibnd,ikk) + g2 * weight
                    !
                 ENDDO !jbnd
                 !
              ENDDO !ibnd
              !
           ENDDO !imode
           !
        ENDIF ! endif  fsthick
        !
     ENDDO  ! iq's
     !
  ENDDO ! end loop on k
  !
  ! The k points are distributed among pools: here we collect them
  !
  nksqtotf = nkstotf/2 ! odd-even for k,k+q
  !
  ALLOCATE ( xkf_all    ( 3,       nkstotf ), &
             etf_all    ( nbndsub, nkstotf ), &
             sigmar_all ( nbndsub, nkstotf ), &  
             sigmai_all ( nbndsub, nkstotf ), &  
             zi_all     ( nbndsub, nkstotf )  )
  !
#ifdef __PARA
  !
  ! note that poolgather2 works with the doubled grid (k and k+q)
  ! therefore we need to use the sigma array with both grids, even
  ! though one of them is useless. This should be fixed by modifying
  ! poolgather2 (it's a waste of memory).
  !
  CALL poolgather2 ( 3,       nkstotf, nksf, xkf,    xkf_all  )
  CALL poolgather2 ( nbndsub, nkstotf, nksf, etf,    etf_all  )
  CALL poolgather2 ( nbndsub, nkstotf, nksf, sigmar, sigmar_all)
  CALL poolgather2 ( nbndsub, nkstotf, nksf, sigmai, sigmai_all)
  CALL poolgather2 ( nbndsub, nkstotf, nksf, zi,     zi_all)
  CALL mp_sum(fermicount, inter_pool_comm)
  !
  ! test output from each pool
  ! DO ik = 1, nksqf
  !    IF (lgamma) THEN
  !       ikk = ik
  !    ELSE
  !       ikk = 2 * ik - 1
  !    ENDIF
  !    WRITE(1000+my_pool_id,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf(:,ikk)
  !    WRITE(1000+my_pool_id,'(6(2x,f12.6))') ( ryd2mev*sigmar(ibnd,ikk), ibnd=ibndmin,ibndmax )
  ! ENDDO
  ! CALL mp_barrier()
  !
#else
  !
  xkf_all = xkf
  etf_all = etf
  sigmar_all = sigmar
  sigmai_all = sigmai
  zi_all     = zi
  !
#endif
  !
  WRITE(6,'(5x,"WARNING: only the eigenstates within the Fermi window are meaningful")') 
  !
  WRITE( stdout, '(/5x,a,i5,a,i5/)' ) &
    'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nksqtotf*nxqf
  !
  DO ik = 1, nksqtotf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     !
     WRITE(stdout,'(/5x,"ik = ",i5," coord.: ", 3f9.5)') ik, xkf_all(:,ikk)
     WRITE(stdout,'(5x,a)') repeat('-',67)
     !
     DO ibnd = ibndmin, ibndmax
        !
        ! note that ekk does not depend on q 
        ekk = etf_all (ibnd, ikk) - ef0
        !
        ! calculate Z = 1 / ( 1 -\frac{\partial\Sigma}{\partial\omega} )
        zi_all (ibnd,ikk) = one / ( one + zi_all (ibnd,ikk) )
        !
!        WRITE(stdout, 102) ibnd, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
!            ryd2mev * sigmai_all (ibnd,ikk),  zi_all (ibnd,ikk)
        WRITE(stdout, 103) ik, ryd2ev * ekk, ryd2mev * sigmar_all (ibnd,ikk), &
              ryd2mev * sigmai_all (ibnd,ikk),  zi_all (ibnd,ikk)
       !
    ENDDO
    WRITE(stdout,'(5x,a/)') repeat('-',67)
    !
  ENDDO
  !
100 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'E( ',i3,' )=',f9.3,' eV   Re[Sigma]=',f9.3,' meV   Im[Sigma]=',f9.3,' meV     Z=',f9.3)
103 format(5x,'k( ',i6,' )=',f10.4,' eV   Re[Sigma]=',f10.4,' meV   Im[Sigma]=',f10.4,' meV     Z=',f9.3)
  !
  RETURN
  !
  END SUBROUTINE selfen_elec
