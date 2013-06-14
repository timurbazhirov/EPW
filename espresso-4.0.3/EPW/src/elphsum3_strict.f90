  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/elphsum - Quantum-ESPRESSO group                 
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine elphsum3_strict
  !-----------------------------------------------------------------------
  !
  !      Sum over BZ of the electron-phonon matrix elements el_ph_mat
  !      and symmetrize over the star of q
  !
  !      Written by Feliciano Giustino for Wannier interpolation.
  !      Based on the original routine elphsum.f90. See endnote for 
  !      rationale of coefficients and conversion factors.
  !
  !      Here I enforce the exact selection rule on the 
  !      energies of the initial and final electronic states:
  !      \frac{ f( E_{n,k} ) - f( E_{m,k+q} ) }{ E_{m,k+q} - E_{n,k} } &
  !           \delta ( E_{m,k+q} - E_{n,k} -\omega_{q\nu} )
  !      as opposed to
  !      \delta ( E_{n,k} - E_{\rm F} ) \delta ( E_{m,k+q} - E_{\rm F} )
  !
  !      The summation on k is performed only for those states with energy
  !      less than fsthick (Fermi Surface thickness) from the FS.
  !
#include "f_defs.h"
  USE ions_base, ONLY : nat, tau, ityp
  USE io_global, ONLY : stdout
  use pwcom
  USE kinds, only : DP
  use phcom
  use epwcom
  use el_phon
  USE control_flags, ONLY :  modenum, noinv, iverbosity 
#ifdef __PARA
!  use para
   use mp_global, ONLY : me_pool,inter_pool_comm 
   use mp,        ONLY : mp_sum
#endif
  implicit none
  !
  real(kind=DP), external :: efermig
  real(kind=DP), parameter :: pibytwo = 3.1415926 / 2.d0, ryd2ev = 13.6058,  &
                              ryd2ghz = 3.289828d6, two = 2.d0, zero = 0.d0, &
                              ryd2mev = 13605.8
  integer :: ik, ikk, ikq, ibnd, jbnd, ipert, jpert, nu, mu, &
       vu, pu, iuelph, ios, iw
  real(kind=DP) :: weight, w0g1, w0g1_ , w0g2, wgauss,dosef,          &
       DOs_ef, ef1, lambda, gamma, wq, ekk, ekq, dw
  complex(kind=DP) :: cw0g2
  complex(kind=DP) :: el_ph_sum (3*nat,3*nat,nw), gamma_mat(3*nat,3*nat)
  complex(kind=DP), allocatable :: el_ph_mat_f ( :, :, :)
  complex(kind = 8), parameter :: ci = (0.d0, 1.d0), cone = (1.d0, 0.d0), &
                 czero = (0.d0, 0.d0)
  !
  ! local variables
  !
  integer :: nq, isq (48), imq
  ! degeneracy of the star of q
  ! index of q in the star of a given sym.op.
  ! index of -q in the star of q (0 if not present)
  ! counter on atoms
  ! counter on atomic type
  ! counter on modes
  ! counter on representation
  real(kind=DP) :: sxq (3, 48)
  ! list of vectors in the star of q
  !
  ! the following are to avoid the mess by star_q.f90 
  !
  integer :: nsym_ , s_ (3, 3, 48), invs_ (48), irt_ (48, nat)
  real(kind=DP) :: rtau_ (3, 48, nat)
  !
  integer :: nrec, fermicount
  logical :: tfermikk, tfermikq
  external dos_ef
  !
  WRITE(stdout, '(/5x,a)' ) 'Entering elphsum3: e-p matrix for the star of q'
  !
  IF ( elinterp ) then
    WRITE(stdout, '(/5x,a)' ) 'Electron interpolation is ON'
  ELSE
    WRITE(stdout, '(/5x,a)' ) 'Electron interpolation is OFF'
    !
    ! use the non-interpolated values and all the bands
    !
    nksf = nks
    nksqf = nksq
    nbndsub = nbnd
    !
    allocate ( etf ( nbnd, nksf ), wkf (nksf) )
    !
    etf = et
    wkf = wk
    !
    ibndmin = 1
    ibndmax = nbnd
    !
  ENDIF
  !
  allocate ( el_ph_mat_f ( ibndmax-ibndmin+1, ibndmax-ibndmin+1, 3*nat) )
  !
  IF ( fsthick.lt.1.d3 ) &
    WRITE(stdout, '(/5x,a,f10.6,a)' ) &
      'Fermi Surface thickness = ', fsthick, ' Ry'
  WRITE(stdout, '(/5x,a,f10.6,a)' ) &
    'Golden Rule strictly enforced with T = ',eptemp, ' Ry'
  !
  IF ( selfen_type .eq. 1 ) then
     WRITE(stdout, '(/5x,"Calculating  REAL part of the phonon selfenergy")' )
     WRITE(stdout, '( 5x,"-- i.e. the approximate frequency shift --")' )
  ELSEif ( selfen_type .eq. 2 ) then
     WRITE(stdout, '(/5x,"Calculating -IMAG part of the phonon selfenergy")' )
     WRITE(stdout, '( 5x,"-- i.e. the *half* width at half maximum --")' )
  ELSE
     CALL errore ('elphsum3_strict','wrong selfen_type',1)
  ENDIF
  !
  ! the frequency sampling must be finer than the small parameter used in the
  ! phonon selfenergy (degaussw)
  !
  dw = ( wmax - wmin ) / float (nw-1)
  !
  IF (filelph.ne.' ') then
#ifdef __PARA
     ! parallel case: only first node writes
     !if (me.ne.1) then
     IF (me_pool /= 0) then 
        iuelph = 0
     ELSE
#endif
     iuelph = 4
     open (unit = iuelph, file = filelph, status = 'unknown', err = &
          100, iostat = ios)
100  call errore ('elphon', 'opening file'//filelph, abs (ios) )
     rewind (iuelph)
     WRITE (iuelph, 9001) xq, degaussw, ngaussw, 3 * nat, nw, nq
     WRITE (iuelph, '(6e14.6)') (w2 (nu) , nu = 1, nmodes)
#ifdef __PARA
     end if
#endif
  ELSE
     iuelph = 0
  ENDIF

  !
  ! Fermi level and corresponding DOS
  !
  !   Note that the weights of k+q points must be set to zero here
!  CALL efermig (etf, nbndsub, nksf, nelec, wkf, degaussw, ngaussw, ef1)
  ef1 = efermig(etf, nbndsub, nksf, nelec, wkf, degaussw, ngaussw,0,isk)
  DOsef = dos_ef (ngaussw, degaussw, ef1, etf, wkf, nksf, nbndsub)
  !   N(Ef) in the equation for lambda is the DOS per spin 
  DOsef = dosef / two   
  !
  WRITE (6, 9000) degaussw, ngaussw
  WRITE (6, 9005) dosef, ef1 * ryd2ev
  IF (iuelph.ne.0) then
    WRITE (iuelph, *) 
    Write (iuelph, 9000) degaussw, ngaussw
    WRITE (iuelph, 9005) dosef, ef1 * ryd2ev
  ENDIF
  !
  ! Sum over bands with gaussian weights
  !
  el_ph_sum = czero
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
    tfermikk = .false.
    tfermikq = .false.
    !
    DO ibnd = 1, nbndsub
      IF ( abs(etf (ibnd, ikk) - ef) .lt. fsthick ) tfermikk = .true.
      IF ( abs(etf (ibnd, ikq) - ef) .lt. fsthick ) tfermikq = .true.
    ENDDO
    !
    IF (tfermikk.and.tfermikq) then
      !
      fermicount = fermicount + 1
      !
      IF ( elinterp ) then
        !
        ! read from file el_ph_mat_f(ibndmin:ibndmax,ibndmin:ibndmax,:) for this ik
        !
        IF (iverbosity.eq.1) &
        WRITE(stdout,'(5x,a,i6)',advance='no') 'Reading el_ph_mat_f from file for ik =',ik
        !
        DO ipert = 1, 3 * nat
          nrec = (ipert - 1) * nksqf + ik
          CALL dasmio ( el_ph_mat_f(:,:,ipert), ibndmax-ibndmin+1, &
            lrepmatf, iunepmatf, nrec, -1)
        ENDDO
        !
        IF (iverbosity.eq.1) &
        WRITE(stdout,*) ' ... done'
        !
      ELSE
        !
        ! use noninterpolated matrix (all bands)
        !
        DO ipert = 1, 3 * nat
          DO ibnd = 1, nbnd
           DO jbnd = 1, nbnd
              el_ph_mat_f(ibnd, jbnd, ipert) = el_ph_mat(ibnd, jbnd, ik, ipert)
           ENDDO
          ENDDO
        ENDDO
        !
      ENDIF
      !
      IF (iverbosity.eq.1) &
      WRITE(stdout,'(12x,a)',advance='no') 'Fermi surface average '
      !
      DO ibnd = 1, ibndmax-ibndmin+1
        !
        ekk = etf (ibndmin-1+ibnd, ikk) - ef1
        !
        DO jbnd = 1, ibndmax-ibndmin+1
          !
          ekq = etf (ibndmin-1+jbnd, ikq) - ef1
          !
          w0g1_ = wgauss( -ekq/eptemp, -99) - wgauss( -ekk/eptemp, -99) 
          !
          !  Loop over frequency
          !
          DO iw = 1, nw
            !
            wq = wmin + float (iw-1) * dw
            !
            IF ( wq .gt. eps_acustic ) then
              w0g1 = w0g1_ / wq 
            ELSE
              w0g1 = zero
            ENDIF
            !
            cw0g2 =  cone / ( ekq - ekk - wq - ci * degaussw ) / pi
            !
            IF ( selfen_type .eq. 1 ) then
               w0g2 =    real ( cw0g2 ) 
            ELSE 
               w0g2 = - aimag ( cw0g2 ) 
            ENDIF
            ! 
            weight = wkf (ikk) * w0g1 * w0g2
            !
            DO jpert = 1, 3 * nat
              DO ipert = 1, 3 * nat
                 el_ph_sum (ipert, jpert, iw) = el_ph_sum (ipert, jpert, iw)   &
                    + weight * conjg (el_ph_mat_f (jbnd, ibnd, ipert))         &
                             *        el_ph_mat_f (jbnd, ibnd, jpert)
              ENDDO
            ENDDO
            !
          ENDDO
          !
        ENDDO
      ENDDO
      !
      IF (iverbosity.eq.1) &
      WRITE(stdout,*) ' ... done'
      !
      ! endif tfermi
    ENDIF
    !
  ENDDO
  !
  WRITE( stdout, '(/14x,a,i5,a,i5)' ) &
    'Number of (k,k+q) pairs on the Fermi surface: ',fermicount, ' out of ', nksqf
  !
  !
#ifdef __PARA
  !
  ! collect contributions from all pools (sum over k-points)
  !
  CALL mp_sum(el_ph_sum,inter_pool_comm)
!  call poolreduce (2 * 3 * nat * 3 * nat * nw, el_ph_sum)
  !
#endif
  !
  deallocate ( el_ph_mat_f )
  !
  !  Loop over frequency again
  !
  DO iw = 1, nw
    !
    wq = wmin + float (iw-1) * dw
    !
    !   Here el_ph_sum is the matrix in the 'pattern' representation
    !   I need the symmetrized matrix in cartesian basis for all q
    !   in the star
    !
    !   Symmetrizes the gamma matrix w.r.t. the small group of q.
    !   It is a second rank tensor obtained from the tensor product
    !   of two first rank tensors therefore it transforms as the 
    !   dynamical matrix
    ! 
    !   The following call also brings ep_ph_sum on the cartesian basis
    !
    CALL symdyn_munu (el_ph_sum(:,:,iw), u, xq, s, invs, rtau, irt, irgq, at, &
       bg, nsymq, nat, irotmq, minus_q)
    !
    !------------------------------------------------------------------------
    !
    ! ... What follows is only for printing out information needed 
    ! ... when noninterpolated calculations are done (summary and check)
    !
    ! ... go from cartesian to the vibrational eigenmodes
    ! ... dyn(mu,nu) contains the eigenmodes divided by the masses
    ! ... i.e. the true eigendisplacements
    !
    gamma_mat = zero
    !
    DO nu = 1, nmodes
       !
       DO pu = 1, 3 * nat
          !
          DO mu = 1, 3 * nat
          DO vu = 1, 3 * nat
             !
             gamma_mat (nu, pu) = gamma_mat (nu, pu) & 
                + conjg (dyn (mu, nu) ) * el_ph_sum (mu, vu, iw) * dyn (vu, pu) 
             !
          ENDDO
          ENDDO
          !
       ENDDO
       !
    ENDDO
    !
    WRITE(stdout, 9008) iw, wq
    !
    gamma_mat = pibytwo * gamma_mat 
    !        
    DO nu = 1, nmodes
       gamma = real( gamma_mat (nu, nu) )
       IF (sqrt (abs (w2 (nu) ) ) .gt. eps_acustic ) then
          lambda = gamma / pi / w2 (nu) / dosef
       ELSE
          lambda = zero
       ENDIF
       gamma = gamma * ryd2mev 
!      gamma = gamma * ryd2ghz
       WRITE (6, 9011) nu, lambda, gamma 
       !
    ENDDO
    !
    ! end loop on iw
    !
  ENDDO
  !
  ! star_q must be called after symdyn_munu, otherwise we symmetrize with the
  ! wrong irt,s,rtau 
  !
  !   Generates the star of q
  !
  CALL star_q (xq, at, bg, ibrav, symm_type, nat, tau, ityp, nr1, &
    nr2, nr3, nsym_ , s_ , invs_ , irt_ , rtau_ , nq, sxq, isq, imq, noinv, &
    modenum)  
  !
  !   Rotates and writes the gamma matrices of the star of q for this iw
  !
  DO iw = 1, nw
    CALL q2qstar_ph (el_ph_sum(:,:,iw), at, bg, nat, nsym_ , s_ , invs_ , irt_ , rtau_ , &
         nq, sxq, isq, imq, iuelph)
    !
  ENDDO
  !
  9000 format(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
  9001 format(3f15.8,f8.3,4i8)
  9005 format(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
  9006 format(5x,'double delta at Ef =',e10.3)
  9007 format(5x,'Gaussian Broadening: ',f10.6,       &
                 ' --> double delta at Ef =',e10.3,   &
                 ' , double delta / Nf^2 = ', e12.5)
  9008 format(/5x,"Selection rule enforced with w(",i3,") = ",f10.6," Ry")
  9010 format(5x,'lambda(',i2,')=',f10.6,'   gamma=',f10.6,' GHz')
  9011 format(5x,'lambda( ',i3,' )=',e20.10,'   gamma=',e20.10,' meV')
  !
  IF (iuelph.ne.0) close (unit = iuelph)
  !
  ! ---------------------------------------------------------------------
  !
  ! Note for efermig call:
  ! Recalculate the Fermi energy Ef=ef1 and the DOS at Ef, dosef = N(Ef)
  ! for this gaussian broadening
  ! Note that the weights of k+q points must be set to zero for the
  ! following call to yield correct results
  !
  !
  ! el_ph_sum(mu,nu)=\sum_k\sum_{i,j}[ <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}>
  !                                  x <psi_{k+q,j}|dvscf_q(nu)*psi_{k,i}>
  !                                  x \delta(e_{k,i}-Ef) \delta(e_{k+q,j}
  !
  ! gamma = \pi \sum_k\sum_{i,j} \delta(e_{k,i}-Ef) \delta(e_{k+q,j}-Ef)
  !         | \sum_mu z(mu,nu) <psi_{k+q,j}|dvscf_q(mu)*psi_{k,i}> |^2
  ! where z(mu,nu) is the mu component of normal mode nu (z = dyn)
  ! gamma(nu) is the phonon linewidth of mode nu
  !
  ! the factor 2 comes from the factor sqrt(hbar/2/M/omega) that appears
  ! in the definition of the electron-phonon matrix element g
  ! The sqrt(1/M) factor is actually hidden into the normal modes
  !
  ! The factor N(Ef)^2 that appears in most formulations of el-ph interact
  ! is absent because we sum, not average, over the Fermi surface.
  ! The factor 2 is provided by the sum over spins (i.e. sum(w_k) = 2 @FG)
  !
  ! lambda is the adimensional el-ph coupling for mode nu:
  ! lambda(nu)= gamma(nu)/(pi N(Ef) \omega_{q,nu}^2)
  ! 3.289828x10^6 is the conversion factor from Ry to GHz
  ! ---------------------------------------------------------------------
  !
  end subroutine elphsum3_strict
