  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine nesting_fn
  !-----------------------------------------------------------------------
  !
  !  10/2009 Compute the nesting function.  Depends only on electronic
  !          eigenvalues.  Copied from selfen_phon
  !  
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     only : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  use phcom,     only : nmodes
  use epwcom,    only : nbndsub, lrepmatf, iunepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, iuetf,     &
                        wmin, wmax, nw, nbndskip, a2f, epf_mem, etf_mem, & 
                        nsmear, delta_smear
  use pwcom,     only : nelec, ef, isk, nbnd
  use el_phon,   only : epf17, ibndmax, ibndmin, etf, &
                        etfq, wkf, xqf, wqf, nksf, nxqf,   &
                        nksqf, wf, nkstotf, xkf, xqf
#ifdef __PARA
  use mp,        only : mp_barrier,mp_sum
  use mp_global, only : me_pool,inter_pool_comm,my_pool_id
#endif
  !
  implicit none
  !
  real(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0), czero=(0.d0,0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, nrec, iq, fermicount, ismear
  real(kind=DP) :: ekk, ekq, ef0, wgkk, wgkq,  &
     weight, w0g1, w0g2, w0gauss, wgauss, dosef, dos_ef, gamma, &
     degaussw0, eptemp0
  !
  real(kind=DP), external :: efermig
  logical :: already_skipped
  !
  !
  !
  WRITE(6,'(/5x,a)') repeat('=',67)
  WRITE(6,'(5x,"Nesting function using double delta approach")') 
  WRITE(6,'(5x,a/)') repeat('=',67)
  !
  IF ( fsthick.lt.1.d3 ) &
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
      'Fermi Surface thickness = ', fsthick, ' Ry'
  WRITE(stdout, '(/5x,a,f10.6,a)' ) &
    'Golden Rule strictly enforced with T = ',eptemp(1), ' Ry'
  already_skipped = .false.
! here we loop on smearing values
  DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1)*delta_smear+degaussw
     eptemp0 = (ismear-1)*delta_smear+eptemp(1)
  ! 
  !
  ! Fermi level and corresponding DOS
  !
  ! here we take into account that we may skip bands when we wannierize
  ! (spin-unpolarized)
  ! 
  IF (nbndskip.gt.0) then
     IF (.not. already_skipped) then 
        nelec = nelec - two * nbndskip
        ! this is necessary for a loop on smearings, as each time through
        ! we do not want to keep subtracting electrons.
        already_skipped = .true.
        WRITE(stdout, '(/5x,"Skipping the first ",i4," bands:")' ) nbndskip
        WRITE(stdout, '(/5x,"The Fermi level will be determined with ",f9.5," electrons")' ) nelec     
     ENDIF
  ENDIF
  !
  !   Note that the weights of k+q points must be set to zero here
  !   no spin-polarized calculation here
  ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw0,ngaussw,0,isk)
  DOsef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin
  DOsef = dosef / two
  !
  WRITE (6, 100) degaussw0, ngaussw
  WRITE (6, 101) dosef, ef0 * ryd2ev
  !
  ! loop over all q points of the fine mesh (this is in the k-para case 
  ! it should always be k-para for selfen_phon)
  ! 
  DO iq = 1, nxqf
    !
    CALL start_clock('nesting')
    !
    fermicount = 0
    gamma = zero
    !
    DO ik = 1, nksqf
       !
       ikk = 2 * ik - 1
       ikq = ikk + 1
       ! 
       ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
       !
       IF (etf_mem) then
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
       IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
            ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
          !
          fermicount = fermicount + 1
          !
          !
          DO ibnd = 1, ibndmax-ibndmin+1
             !
             !  the fermi occupation for k
             ekk = etf (ibndmin-1+ibnd, ikk) - ef0
             wgkk = wgauss( -ekk/eptemp0, -99)
             !
             DO jbnd = 1, ibndmax-ibndmin+1
                !
                !  the fermi occupation for k+q
                ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                wgkq = wgauss( -ekq/eptemp0, -99)  
                !
                !
                ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                !
                !              weight = wkf (ikk) * (wgkk - wgkq) * &
                !                 aimag ( cone / ( ekq - ekk - ci * degaussw0 ) ) 
                !
                ! the below expression is positive-definite, but also an approximation
                ! which neglects some fine features
                !
                w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
                w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
                weight = pi * wkf (ikk) * w0g1 * w0g2
                !
                gamma  =   gamma  + weight 
                !
             ENDDO ! jbnd
          ENDDO   ! ibnd
          !
          !
       ENDIF ! endif fsthick
       !
       CALL stop_clock('nesting')
       !
    ENDDO ! loop on q
#ifdef __PARA
    !
    ! collect contributions from all pools (sum over k-points)
    ! this finishes the integral over the BZ  (k)
    !
   CALL mp_sum(gamma,inter_pool_comm) 
   CALL mp_sum(fermicount, inter_pool_comm)
    !
#endif
    !
    WRITE(6,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:,iq) , wqf(iq)
    WRITE(6,'(5x,a)') repeat('-',67)
    !
    WRITE(6, 102)  gamma
    WRITE(6,'(5x,a/)') repeat('-',67)
    !
    ! test only
#ifdef __PARA

!    if (me.eq.1) & 
    IF (me_pool == 0) &
#endif
         !
         !
         WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
         'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkstotf/2
 ENDDO
  !
 ENDDO !smears
!  ! generate the Eliashberg spectral function
  !
!  if (a2f) call eliashberg_a2f( gamma_all, lambda_all, gamma_all_v)
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,' Nesting function (q)=',f9.3,' units??')
  !
  end subroutine nesting_fn

  !-----------------------------------------------------------------------
  subroutine nesting_fn_fly (iq )
  !-----------------------------------------------------------------------
  !
  !  compute the imaginary part of the phonon self energy due to electron-
  !  phonon interaction in the Migdal approximation. This corresponds to 
  !  the phonon linewidth (half width). The phonon frequency is taken into
  !  account in the energy selection rule.
  !
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation.  This routine is similar to the one above
  !  but it is only called from within ephwann_shuffle and calculates 
  !  the selfenergy for one phonon at a time.  Much smaller footprint on
  !  the disk
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     only : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout
  use phcom,     only : nmodes
  use epwcom,    only : nbndsub, lrepmatf, iunepmatf, fsthick, &
                        eptemp, ngaussw, degaussw, iuetf,     &
                        wmin, wmax, nw, nbndskip, a2f, epf_mem, etf_mem, &
                        nsmear, delta_smear
  use pwcom,     only : nelec, ef, isk, nbnd
  use el_phon,   only : epf17, ibndmax, ibndmin, etf, &
                        etfq, wkf, xqf, wqf, nksf, nxqf,   &
                        nksqf, wf, nkstotf, xkf, xqf
#ifdef __PARA
  use mp,        only : mp_barrier,mp_sum
  use mp_global, only : me_pool,inter_pool_comm,my_pool_id
#endif
  !
  implicit none
  !
  real(kind=DP), parameter :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
                              two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0), czero=(0.d0,0.d0)
  integer :: ik, ikk, ikq, ibnd, jbnd, nrec, iq, fermicount, ismear
  real(kind=DP) :: ekk, ekq, ef0, wgkk, wgkq, &
     weight, w0g1, w0g2, w0gauss, wgauss, dosef, dos_ef, gamma, &
     degaussw0, eptemp0
  !
  real(kind=DP), external :: efermig
  logical :: already_skipped
  !
  !
  !
  IF (iq.eq.1) then 
     WRITE(6,'(/5x,a)') repeat('=',67)
     WRITE(6,'(5x,"Nesting Function in the double delta approx")')
     WRITE(6,'(5x,a/)') repeat('=',67)
     !
     IF ( fsthick.lt.1.d3 ) &
          WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Fermi Surface thickness = ', fsthick, ' Ry'
     WRITE(stdout, '(/5x,a,f10.6,a)' ) &
          'Golden Rule strictly enforced with T = ',eptemp(1), ' Ry'
     IF (nbndskip.gt.0) then
        nelec = nelec - two * nbndskip
        ! this is necessary for a loop on smearings, as each time through
        ! we do not want to keep subtracting electrons.
     ENDIF
  ENDIF

! here we loop on smearing values
  DO ismear = 1, nsmear
     !
     degaussw0 = (ismear-1)*delta_smear+degaussw
     eptemp0 = (ismear-1)*delta_smear+eptemp(1)
  ! 
  !
  ! Fermi level and corresponding DOS
  !
  ! here we take into account that we may skip bands when we wannierize
  ! (spin-unpolarized)
  ! 
  !
  !   Note that the weights of k+q points must be set to zero here
  !   no spin-polarized calculation here
  ef0 = efermig(etf,nbndsub,nksf,nelec,wkf,degaussw0,ngaussw,0,isk)
  DOsef = dos_ef (ngaussw, degaussw0, ef0, etf, wkf, nksf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin
  DOsef = dosef / two
  !
IF (iq.eq.1) then
  WRITE (6, 100) degaussw0, ngaussw
  WRITE (6, 101) dosef, ef0 * ryd2ev
ENDIF
  !
  !
  CALL start_clock('nesting')
  !
  fermicount = 0
  !
  DO ik = 1, nksqf
     !
     ikk = 2 * ik - 1
     ikq = ikk + 1
     ! 
     ! we read the hamiltonian eigenvalues (those at k+q depend on q!) 
     !
     ! when we see references to iq for file readins, it is always = 1 for on the fly calculations
     IF (etf_mem) then
        etf (ibndmin:ibndmax, ikk) = etfq (ibndmin:ibndmax, ikk,  1)
        etf (ibndmin:ibndmax, ikq) = etfq (ibndmin:ibndmax, ikq,  1)
     ELSE
        nrec = (iq-1) * nksf + ikk
        nrec = ikk
        CALL davcio ( etf (ibndmin:ibndmax, ikk), ibndmax-ibndmin+1, iuetf, nrec, - 1)
        nrec = (iq-1) * nksf + ikq
        nrec = ikq
        CALL davcio ( etf (ibndmin:ibndmax, ikq), ibndmax-ibndmin+1, iuetf, nrec, - 1)
     ENDIF
     !
     ! here we must have ef, not ef0, to be consistent with ephwann_shuffle
     IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .and. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) then
        !
        fermicount = fermicount + 1
        !
           DO ibnd = 1, ibndmax-ibndmin+1
              !
              !  the fermi occupation for k
              ekk = etf (ibndmin-1+ibnd, ikk) - ef0
              wgkk = wgauss( -ekk/eptemp0, -99)
              !
              DO jbnd = 1, ibndmax-ibndmin+1
                 !
                 !  the fermi occupation for k+q
                 ekq = etf (ibndmin-1+jbnd, ikq) - ef0
                 wgkq = wgauss( -ekq/eptemp0, -99)  
                 !
                 ! = k-point weight * [f(E_k) - f(E_k+q)]/ [E_k+q - E_k -w_q +id]
                 ! This is the imaginary part of the phonon self-energy, sans the matrix elements
                 !
!                 weight = wkf (ikk) * (wgkk - wgkq) * &
!                      aimag ( cone / ( ekq - ekk  - ci * degaussw ) ) 
                 !
                 ! the below expression is positive-definite, but also an approximation
                 ! which neglects some fine features
                 !
             w0g1 = w0gauss ( ekk / degaussw0, 0) / degaussw0
             w0g2 = w0gauss ( ekq / degaussw0, 0) / degaussw0
             weight = pi * wkf (ikk) * w0g1 * w0g2
                 !
                 gamma  =   gamma  + weight  
                 !
              ENDDO ! jbnd
           ENDDO   ! ibnd
           !
        !
        !
     ENDIF ! endif fsthick
     !
     CALL stop_clock('nesting')
      !
  ENDDO ! loop on k
#ifdef __PARA
  !
  ! collect contributions from all pools (sum over k-points)
  ! this finishes the integral over the BZ  (k)
  !
  CALL mp_sum(gamma,inter_pool_comm) 
  CALL mp_sum(fermicount, inter_pool_comm)
  !
#endif
  !
  WRITE(6,'(/5x,"iq = ",i5," coord.: ", 3f9.5, " wt: ", f9.5)') iq, xqf(:,iq) , wqf(iq)
  WRITE(6,'(5x,a)') repeat('-',67)
     ! 
  WRITE(6, 102)  gamma
  WRITE(6,'(5x,a/)') repeat('-',67)
  CALL flush(6)
  !
       WRITE( stdout, '(/5x,a,i8,a,i8/)' ) &
       'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nkstotf/2
  !
  !
 ENDDO !smears
  !
100 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 format(5x,' Nesting function (q)=',f9.3,' units??')

end subroutine nesting_fn_fly
!
