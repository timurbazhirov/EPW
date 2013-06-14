  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------------
  SUBROUTINE wann_run
  !---------------------------------------------------------------------
  !
  !  This is the SUBROUTINE which controls the w90 run.  Primarily,        
  !  we get the phases to remove degeneracies in the wfs, and 
  !  call pw2wan90epw 
  !  
  !---------------------------------------------------------------------
  !
  USE kinds,          ONLY : DP
  USE io_global,      ONLY : stdout, ionode, ionode_id
#ifdef __PARA
  USE mp_global,      ONLY : my_pool_id
#endif
  USE io_files,       ONLY : prefix, iunwfc, find_free_unit, iunigk
  USE wvfct,          ONLY : nbnd, npwx
  USE wavefunctions_module,  ONLY: evc
  USE ions_base,      ONLY : nat, atm, tau
  USE pwcom
  USE phcom
  USE epwcom
  USE el_phon
  USE wannier
  use mp,        only : mp_bcast
  !
  implicit none
  !
  ! work variables
  integer :: i,j,k, counter
  !
  ! intend in variables
  real(kind=DP) :: zero_vect(3)
  real(kind=DP), ALLOCATABLE :: tau_wan(:,:)
  integer :: num_kpts
  logical :: spinors, gamma_only_wan
  !
  !
  real(kind=DP), PARAMETER :: bohr2ang = 0.5291772108
  !
  CALL start_clock( 'WANNIER' )
  ! 
  ! Currently EPW does not incorporate LSDA
  spinors = .false.
  ! For small systems we do not do gamma ONLY.  
  gamma_only_wan = .false.
  !
  zero_vect = 0.0
  !
  mp_grid(1) = nk1
  mp_grid(2) = nk2
  mp_grid(3) = nk3
  num_kpts = mp_grid(1)*mp_grid(2)*mp_grid(3)
  !
  IF (num_kpts .ne. nkstot ) call errore('wannierize','inconsistent nscf and elph k-grids',1) 
  IF (nbnd .lt. n_wannier )  call errore('Must have as many or more bands than Wannier functions',1) 
  ALLOCATE ( kpt_latt(3,npk), tau_wan(3,nat) )
  !
  WRITE(stdout, '(5x,a)') repeat("-",67)
  WRITE(stdout, '(a, i2,a,i2,a,i2,a)') "     Wannierization on ", nk1, " x ", nk2, " x ", nk3 , " electronic grid"
  WRITE(stdout, '(5x,a)') repeat("-",67)
  !
  kpt_latt = xk_cryst
  CALL mp_bcast(kpt_latt,ionode_id)
  !
  ! write the short input file for the wannier90 code
  !
  CALL write_winfil
  !
  ! run the wannier90 code to create MLWFs
  !
  CALL pw2wan90epw
  !
  ! project the Wannier functions onto energy space
  !
!  CALL proj_w90
  !
  WRITE(stdout, '(5x,a)') repeat("-",67)
  CALL print_clock( 'WANNIER' )
  WRITE(stdout, '(5x,a)') repeat("-",67)
  !
  end SUBROUTINE wann_run
  !------------------------------------------------------------
  !
  !------------------------------------------------------------
  SUBROUTINE write_winfil
  !------------------------------------------------------------
  !
  !
  !  This SUBROUTINE write the prefix.win file which wannier90.x
  !  needs to run.  Primarily it contains information about the 
  !  windows USEd for the disentanglement, and the initial projections.
  !  JN - 10/2008  projections now in elph.in file  
  !------------------------------------------------------------
#include "f_defs.h"
  !
  USE io_files,    ONLY : prefix, find_free_unit
#ifdef __PARA
  USE io_global,   ONLY : stdout, ionode
#endif
  USE epwcom
  USE el_phon
  USE wannier,     ONLY : iknum
  USE wvfct,       ONLY : nbnd
  !
  implicit none
  !
  integer :: iuwinfil,i
  !
  logical :: random
  !
  iuwinfil = find_free_unit()
  !
#ifdef __PARA
  IF (ionode) THEN
#endif
  !
  IF (nbndsub .gt. nwanxx) call errore('write_winfil',"Too many wannier bands",nbndsub)
  !
  OPEN (unit = iuwinfil, file = trim(prefix)//".win", form = 'formatted')
  !    
  !  more input and options for interfacing with w90 can/will be added later
  WRITE (iuwinfil,'(a)') "begin projections"
  !
  random = .true.
  DO i = 1, nbndsub+1
     IF (proj(i) .ne. ' ') THEN
        WRITE (iuwinfil,*) trim(proj(i))
        random = .false.
     ENDIF
  ENDDO
  !
  IF (random) WRITE(iuwinfil,*) 'random' 
  !
  WRITE (iuwinfil,'(a)') "end projections"
  !
  WRITE (iuwinfil,'("num_wann ",i3)') nbndsub
  WRITE (iuwinfil,'("iprint ",i3)') iprint
  !
  WRITE (iuwinfil, '("dis_win_min ", f9.3)')  dis_win_min
  WRITE (iuwinfil, '("dis_win_max ", f9.3)')  dis_win_max
  WRITE (iuwinfil, '("dis_froz_min ", f9.3)') dis_froz_min
  WRITE (iuwinfil, '("dis_froz_max ", f9.3)') dis_froz_max
  WRITE (iuwinfil, '("num_iter ", i7)')       num_iter
  !
  ! write any extra PARAMETERs to the prefix.win file
  DO i = 1, nwanxx
     IF (wdata(i) .ne. ' ') write(iuwinfil, *) wdata(i)
  ENDDO
  !
  CLOSE (iuwinfil)
  !
#ifdef __PARA
  ENDIF
#endif
  !
  !
  END SUBROUTINE write_winfil



!------------------------------------------------------------
  SUBROUTINE proj_w90
!------------------------------------------------------------
  !
  ! This SUBROUTINE computes the energy projections of
  ! the computed Wannier functions
  ! 07/2010  Needs work.  Right now this sub is nearly worthless  
  !------------------------------------------------------------
#include "f_defs.h"
  !
  USE io_files,    ONLY : prefix, find_free_unit
#ifdef __PARA
  USE mp_global,   ONLY : inter_pool_comm
  USE io_global,   ONLY : stdout, ionode
  USE mp,              ONLY : mp_sum
#endif
  USE epwcom
  USE wannier,     ONLY : iknum, n_wannier
  USE wvfct,       ONLY : nbnd, et
  USE klist,       ONLY : nks, nkstot
  !
  implicit none
  !
  integer :: iuprojfil, ik, ibnd, ne, ie, iwann
  real(kind=DP), PARAMETER :: ryd2eV = 13.6056923
  complex(kind=DP), PARAMETER :: czero = (0.d0,0.d0)
  real(kind=DP)    :: dE, sigma, argv, en, xxq(3)
  real(kind=DP), ALLOCATABLE    ::  proj_wf(:,:)
  complex(kind=DP), ALLOCATABLE ::  cu(:,:,:), cuq(:,:,:)
  !
  iuprojfil = find_free_unit()
  !
  WRITE(6,'(5x,"Computing energy projections")')
  ! dummy var
  xxq = 0.d0
  !
  ! tmp value
  dE = 0.05
  sigma = 2 * dE
  !
  ! maxvalue = dis_win_max + 1
  ! minvalue = dis_win_min - 1
  ne = int( (dis_win_max - dis_win_min + 1) / dE )
  IF (ne .lt. 1)  CALL errore('proj_wan','Problem with disentanglement window',1)
  !
  ALLOCATE (proj_wf(n_wannier, ne+1))
  proj_wf = 0.d0
  !
  ALLOCATE (cu (nbnd, n_wannier, nks) )
  ALLOCATE (cuq(nbnd, n_wannier, nks) )
  !
  CALL loadumat(nbnd, n_wannier, nks, nkstot, xxq, cu, cuq)
  !
  !
  DO iwann = 1, n_wannier
     !
     DO ie = 1, ne
        en = float(ie)/float(ne) * (dis_win_max - dis_win_min + 1.0) + dis_win_min
        !
        DO ik = 1, nks
           DO ibnd = 1, nbnd
              !
              argv = ( et(ibnd,ik)*ryd2eV - en ) **2/ (2 * sigma **2)
              proj_wf(iwann, ie ) = proj_wf(iwann, ie ) + exp(-argv) * real (  cu(ibnd, iwann,ik ) * conjg( cu(ibnd, iwann,ik) ))
              !
           ENDDO
        ENDDO
     ENDDO
  ENDDO
#ifdef __PARA
  !
  ! sum the contributions from all k-points
  CALL mp_sum(proj_wf, inter_pool_comm)
  !
  IF (ionode) THEN
#endif
     !
     OPEN (unit = iuprojfil, file = trim(prefix)//".projw90", form = 'formatted')
     !
     WRITE(iuprojfil, '(5x,"Wannier energy projections")')
     !
     DO ie = 1, ne
        en =  float(ie)/float(ne) * (dis_win_max - dis_win_min + 1) + dis_win_min
        WRITE(iuprojfil, '(f9.3, 25f8.4)' )  en , proj_wf(:, ie)
     ENDDO
     !
     CLOSE (iuprojfil)
#ifdef __PARA
  ENDIF
#endif
  !
  IF ( ALLOCATED(proj_wf)) DEALLOCATE(proj_wf)
  IF ( ALLOCATED(cu))      DEALLOCATE(cu)
  IF ( ALLOCATED(cuq))     DEALLOCATE(cuq)
  !
!------------------------------------------------------------
  END  SUBROUTINE proj_w90
!------------------------------------------------------------
