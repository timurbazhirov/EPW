  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE plot_band
  !-----------------------------------------------------------------------
  !
  !  This subroutine writes output files for phonon dispersion and band structure 
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : at, bg, ibrav, symm_type
  USE ions_base, ONLY : nat, ityp, tau
  USE io_global, ONLY : stdout, ionode_id
  USE io_files,  ONLY : find_free_unit, prefix
  USE phcom,     ONLY : lgamma, nmodes
  USE epwcom,    ONLY : nbndsub, lrepmatf, iunepmatf, fsthick, ngaussw, degaussw, iuetf, &
                        wmin, wmax, nw, nqf1, nqf2, nqf3, parallel_k, parallel_q, &
                        epf_mem, etf_mem, eig_read
  USE pwcom,     ONLY : nelec, ef, isk, bg, at
  USE el_phon,   ONLY : nksf, etf, etfq, ibndmin, ibndmax, &
                        wkf, nksqf, nxqf, wf, wqf, xkf, xqf, nkstotf
#ifdef __PARA
  USE mp,        ONLY : mp_barrier, mp_sum
  USE mp_global, ONLY : me_pool, inter_pool_comm, my_pool_id, npool, mpime
#endif
  !
  IMPLICIT NONE
  !
  REAL(DP), PARAMETER :: ryd2mev = 13605.8, one = 1.d0, ryd2ev = 13.6058, &
            two = 2.d0, zero = 0.d0, pi = 3.14159265358979
  INTEGER :: ik, ikk, ikq, ibnd, jbnd, imode, iq, iw
  INTEGER :: iufileig, iufilfreq
  REAL(kind=DP), ALLOCATABLE :: xkf_all(:,:) , etf_all(:,:)
  !
  INTEGER :: nksqtotf
  !
  IF ( .not. etf_mem ) CALL errore ('plot_band', 'etf_mem should be true', 1)
  !
  nksqtotf =  nkstotf/2
  !
#ifdef __PARA
  IF ( my_pool_id .eq. 0 ) THEN
#endif
  !
  iufilfreq = find_free_unit()
  OPEN(iufilfreq, file = "phband.freq", form = 'formatted')
  WRITE(iufilfreq, '(" plot nbnd =",i4," nks =",i4)') nmodes, nxqf
  DO iq = 1, nxqf
     !
     ! crystal to cartesian coordinates
     CALL cryst_to_cart( 1, xqf(:,iq), bg, 1 )
     WRITE(iufilfreq,'(10x,3f10.6)') xqf(:,iq)
     ! back from cartesian to crystal coordinates
     CALL cryst_to_cart( 1, xqf(:,iq), at, -1 )
     WRITE(iufilfreq,'(1000f10.4)') (wf(imode,iq)*ryd2mev, imode=1,nmodes)
     !
  ENDDO
  CLOSE(iufilfreq)
  !
#ifdef __PARA
  ENDIF
#endif
  !
  DO ik = 1, nksqf
     !
     IF (lgamma) THEN
        ikk = ik
        ikq = ik
     ELSE
        ikk = 2 * ik - 1
        ikq = ikk + 1
     ENDIF
     etf(:,ikk) = etfq(:,ikk,1)
     !
  ENDDO
  !
  ALLOCATE ( xkf_all( 3, nkstotf) , etf_all( nbndsub, nkstotf) )
  !
#ifdef __PARA
  !
  CALL poolgather2 ( 3,       nkstotf, nksf, xkf, xkf_all  )
  CALL poolgather2 ( nbndsub, nkstotf, nksf, etf, etf_all  )
  !
#else
  !
  xkf_all = xkf
  etf_all = etf
  !
#endif
  !
  iufileig = find_free_unit()
  OPEN(iufileig, file = "band.eig", form = 'formatted')
  WRITE(iufileig, '(" plot nbnd =",i4," nks =",i4)') nbndsub, nksqtotf
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
     ! crystal to cartesian coordinates
     CALL cryst_to_cart( 1, xkf_all(:,ikk), bg, 1 )
     WRITE(iufileig,'(10x,3f10.6)') xkf_all(:,ikk)
     ! back from cartesian to crystal coordinates
     CALL cryst_to_cart( 1, xkf_all(:,ikk), at, -1 )
     WRITE(iufileig,'(1000f10.4)') (etf_all(ibnd,ikk)*ryd2ev, ibnd=1,nbndsub)
     !
  ENDDO
  CLOSE(iufileig)
  !
  DEALLOCATE( xkf_all )
  DEALLOCATE( etf_all )
  !
100 FORMAT(5x,'Gaussian Broadening: ',f7.3,' Ry, ngauss=',i4)
101 FORMAT(5x,'DOS =',f10.6,' states/spin/Ry/Unit Cell at Ef=',f10.6,' eV')
102 FORMAT(5x,'E( ',i3,' )=',f9.3,' eV ')
  !
  RETURN
  !
  END SUBROUTINE plot_band
