  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE write_ephmat_fly( iq )
  !-----------------------------------------------------------------------
  !
  !  This subroutine writes the elph matrix elements in a format required 
  !  by Eliashberg equations
  ! 
  !  Use matrix elements, electronic eigenvalues and phonon frequencies
  !  from ep-wannier interpolation
  !
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg
  USE io_global, ONLY : stdout, ionode_id
  USE io_files,  ONLY : find_free_unit, prefix, tmp_dir
  USE phcom,     ONLY : lgamma, nmodes
  USE epwcom,    ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, ngaussw, degaussw, iuetf, &
                        wmin, wmax, nw, nbndskip, ecutse, parallel_k, parallel_q, phonselfen, &
                        epf_mem, etf_mem, eig_read, nkf1, nkf2, nkf3, nqf1, nqf2, nqf3, eps_acustic
  USE pwcom,     ONLY : nelec, ef, isk, bg, at, nkstot
  USE el_phon,   ONLY : etf, ibndmin, ibndmax, nksf, etfq, epf17, &
                        wkf, nksqf, nxqf, wf, wqf, xkf, xqf, nkstotf
  USE pwcom,     ONLY : ibrav
  USE mp_global, ONLY : me_pool, inter_pool_comm
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id, npool, mpime
#endif
  !
  IMPLICIT NONE
  REAL(DP), PARAMETER :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
     two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  INTEGER :: ik, ikk, ikq, ibnd, jbnd, imode, iq, nrec, nksqtotf
  INTEGER :: icount, iufileig, iufilfreq, iufileph
  REAL(DP) :: wq, ef0, dosef, dos_ef, g2
  REAL(kind=DP), ALLOCATABLE :: etf_all(:,:)
  COMPLEX(kind=DP) epf (ibndmax-ibndmin+1, ibndmax-ibndmin+1)
  LOGICAL :: already_skipped
  REAL(DP), EXTERNAL :: efermig
  CHARACTER (len=256) :: tempfile, name1
  CHARACTER (len=3) :: filelab
  !
  IF (iq.eq.1) THEN
     !
     IF (.not.phonselfen) THEN 
        WRITE(6,'(/5x,a)') repeat('=',67)
        WRITE(6,'(5x,"Write electron-phonon matrix elements for Eliashberg equations")') 
        WRITE(6,'(5x,a/)') repeat('=',67)
        !
!        WRITE( stdout, '(/5x,a,i7)' ) 'Nr of pools npool = ', npool
        WRITE( stdout, '(/5x,a,i7)' ) '2 x nr of k points nkstotf = ', nkstotf
        WRITE( stdout, '(/5x,a,i7)' ) 'Nr of k points in the pool nksf = ', nksf
        WRITE( stdout, '(/5x,a,i7)' ) 'Nr of k blocks (k,k+q pairs) in the pool nksqf = ',  nksqf
        WRITE( stdout, '(/5x,a,i7,a,i7)' ) 'Nr of q points nxqf = ', nxqf
        WRITE( stdout, '(/5x,a,i7,a,i7)' ) 'Nr of phonon modes nmodes = ', nmodes
        WRITE( stdout, '(/5x,a,i7,a,i7)' ) 'ibndmin = ', ibndmin, ' ibndmax = ', ibndmax
        IF ( fsthick .lt. 1.d3 ) &
           WRITE(stdout, '(/5x,a,f10.6,a)' ) 'Fermi Surface thickness = ', fsthick, ' Ry'
     ENDIF
     !
     ! here we take into account that we may skip bands when we wannierize
     ! (spin-unpolarized)
     !
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
     ! etf(:,ikk) = etfq(:,ikk,1), etf(:,ikq) = etfq(:,ikk,nxqf) - last q point
     ! since wkf(:,ikq) = 0 these bands do not bring any contribution to ef0 or dosef
     !
     ef0 = efermig(etf, nbndsub, nksf, nelec, wkf, degaussw, ngaussw, 0, isk)
     !
     !   if 'fine' Fermi level differs by more than 250 meV, there is probably
     !   something wrong with the wannier functions
     IF ( abs(ef0 - ef) * ryd2eV .gt. 0.5 .and. (.not.eig_read) ) &
        CALL errore ('selfen_elec', 'Something wrong with Ef, check MLWFs', 1)
     !
     dosef = dos_ef(ngaussw, degaussw, ef0, etf, wkf, nksf, nbndsub)
     !   N(Ef) in the equation for lambda is the DOS per spin
     dosef = dosef / two
     !
     IF ( .not.phonselfen ) THEN 
        WRITE (6, 100) degaussw, ngaussw
        WRITE (6, 101) dosef, ef0 * ryd2ev
        WRITE (6, 101) dosef, ef  * ryd2ev
     ENDIF
     !
     DO ik = 1, nksqf ! loop over k-points
        ikk = 2*ik - 1
        etf(:,ikk) = etfq(:,ikk,1)
     ENDDO
     !
     ALLOCATE ( etf_all( nbndsub, nkstotf) )
#ifdef __PARA
     !
     ! collect contributions from all pools (sum over k-points)
     CALL poolgather2 ( nbndsub, nkstotf, nksf, etf, etf_all )
     !
#endif
     !
     ! write eigenvalues to file
     iufileig = find_free_unit()
     tempfile = trim(tmp_dir) // trim(prefix) // '.egnv'
     OPEN(iufileig, file = tempfile, form = 'formatted')
     WRITE(iufileig,'(i7)') ibndmax-ibndmin+1
     nksqtotf =  nkstotf/2
     DO ik = 1, nksqtotf
        ikk = 2 * ik - 1
        WRITE(iufileig,'(20e18.9)') etf_all(ibndmin:ibndmax,ikk)*ryd2ev
     ENDDO
     CLOSE(iufileig)
     DEALLOCATE( etf_all )
     !
  ENDIF ! iq
  !
  ! write phonon frequencies to file
#ifdef __PARA 
  IF ( my_pool_id .eq. 0 ) THEN
#endif
  IF ( iq .eq. 1 ) THEN
     iufilfreq = find_free_unit()
     name1 = trim(tmp_dir) // trim(prefix) // '.freq'
     OPEN(iufilfreq, file = name1, form = 'formatted')
     WRITE(iufilfreq,'(i7)') nmodes
     WRITE(iufilfreq,'(20e18.9)') wf(:,iq)*ryd2ev
     CLOSE(iufilfreq)
  ELSE
     OPEN(iufilfreq, file = name1, access = 'append', form = 'formatted')
     WRITE(iufilfreq,'(120e18.9)') wf(:,iq)*ryd2ev
     CLOSE(iufilfreq)
  ENDIF
#ifdef __PARA
  ENDIF
#endif
  !
  ! write the e-ph matrix elements in the Bloch representation on the fine mesh
  ! in .ephmat files (one for each pool)
  !
  iufileph = find_free_unit()
  tempfile = trim(tmp_dir) // trim(prefix) // '.ephmat'
#ifdef __PARA
  CALL set_ndnmbr(0,my_pool_id+1,1,npool,filelab)
  tempfile = trim(tmp_dir) // trim(prefix) // '.ephmat' // filelab
#endif
  IF ( iq .eq. 1 ) THEN 
     OPEN(iufileph, file = tempfile, form = 'formatted')
  ELSE
     OPEN(iufileph, file = tempfile, access = 'append', form = 'formatted')
  ENDIF
  !
#ifdef __PARA 
  IF ( my_pool_id .eq. 0 ) THEN
  IF ( iq .eq. 1 ) THEN
     WRITE(iufileph,'(5i7,4e18.9)') npool, nksqtotf, nxqf, nmodes, ibndmax-ibndmin+1, &
                     ef*ryd2ev, ef0*ryd2ev, dosef/ryd2ev, fsthick*ryd2ev
  ENDIF
  ENDIF
#endif
  !
#ifdef __PARA
  IF ( iq .eq. 1 ) WRITE(iufileph,'(2i7)') my_pool_id+1, nksqf
#endif
  DO ik = 1, nksqf ! loop over all k points of the fine mesh
     !
     ikk = 2 * ik - 1
     ikq = ikk + 1
     !
     ! we read the hamiltonian eigenvalues (those at k+q depend on q!)
     !
     ! when we see references to iq for file readinq, it is always = 1 for on the fly calculations
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
     icount = 0
     IF ( ( minval ( abs(etf (:, ikk) - ef) ) .lt. fsthick ) .AND. &
          ( minval ( abs(etf (:, ikq) - ef) ) .lt. fsthick ) ) THEN
        icount = 1
     ENDIF
     IF ( icount .eq. 1 ) THEN
        !
        DO imode = 1, nmodes ! phonon modes
           !
           wq = wf(imode, iq)
           IF ( wq > eps_acustic ) THEN
              !
              !  we read the e-p matrix from disk / memory
              !
              IF (etf_mem) THEN
                 epf(:,:) = epf17 ( ik,  1, :, :, imode)
              ELSE
                 nrec = (imode-1) * nksqf + ik
                 CALL dasmio ( epf, ibndmax-ibndmin+1, lrepmatf, iunepmatf, nrec, -1)
              ENDIF
              !
              DO ibnd = 1, ibndmax-ibndmin+1
                 !
                 DO jbnd = 1, ibndmax-ibndmin+1
                    !
                    ! here we take into account the zero-point sqrt(hbar/2M\omega)
                    ! with hbar = 1 and M already contained in the eigenmodes
                    ! g2 is Ry^2, wkf must already account for the spin factor
                    !
                    g2 = abs( epf(jbnd, ibnd) )**two / ( two * wq )
                    WRITE(iufileph,'(e18.9)') g2*ryd2ev*ryd2ev
                    !
                 ENDDO ! jbnd
                 !
              ENDDO ! ibnd
              !
           ENDIF ! wq
           !
        ENDDO ! imode
        !
     ENDIF ! icount
     !
  ENDDO ! ik's
  CLOSE(iufileph)
  !
!#ifdef __PARA
!  CALL mp_barrier()
!#endif
  !
100 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'E( ',i3,' )=',f9.3,' eV ')
  !
  END SUBROUTINE write_ephmat_fly
  !                                                                            
  !-----------------------------------------------------------------------

