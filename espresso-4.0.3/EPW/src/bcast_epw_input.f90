  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Code adapted from PH/bcast_ph_input - Quantum-ESPRESSO group
  ! 09/2009 Very little of this subroutine in necessary.  Many 
  ! excess variables
  !
  !-----------------------------------------------------------------------
  SUBROUTINE bcast_ph_input
  !-----------------------------------------------------------------------
  !
  !     In this routine the first processor sends the input to all
  !     the other processors
  !
  !
#ifdef __PARA
#include "f_defs.h"

  USE pwcom
  USE phcom
  USE epwcom
  USE mp, only: mp_bcast
  USE io_files
  USE control_flags, ONLY: iverbosity
  USE ions_base, ONLY : amass
  USE printout_base, ONLY : title   ! title of the run 


  implicit none
  integer :: root = 0
  !
  ! logicals
  !
  CALL mp_bcast (lgamma, root)
  CALL mp_bcast (epsil, root)
  CALL mp_bcast (trans, root)
  CALL mp_bcast (zue, root)
  CALL mp_bcast (elph, root)
  CALL mp_bcast (lnscf, root)
  CALL mp_bcast (ldisp, root)
  CALL mp_bcast (tshuffle, root)  ! 
  CALL mp_bcast (tshuffle2, root) !
  CALL mp_bcast (elecselfen, root)!
  CALL mp_bcast (phonselfen, root)!
  CALL mp_bcast (ephmatwrite, root)! RM
  CALL mp_bcast (band_plot, root)! RM
  CALL mp_bcast (fly, root)!
  CALL mp_bcast (vme, root)!
  CALL mp_bcast (recover, root)!
  CALL mp_bcast (epbread, root)   !
  CALL mp_bcast (epbwrite, root)  !
  CALL mp_bcast (phinterp, root)  !
  CALL mp_bcast (elinterp, root)  !
  CALL mp_bcast (tphases, root)   !
  CALL mp_bcast (epstrict, root)  !
  CALL mp_bcast (fsthick, root)   !
  CALL mp_bcast (eptemp, root)    !
  CALL mp_bcast (wmin, root)      !
  CALL mp_bcast (wmax, root)      !
  CALL mp_bcast (epwread, root)   !
  CALL mp_bcast (epwwrite, root)  !
  CALL mp_bcast (specfun, root)   !
  CALL mp_bcast (wannierize, root)! JN
  CALL mp_bcast (write_wfn, root) ! 
  CALL mp_bcast (kmaps, root) ! 
  CALL mp_bcast (nest_fn, root) ! 
  CALL mp_bcast (indabs, root) ! 
  CALL mp_bcast (twophoton, root) ! 
  CALL mp_bcast (eig_read, root) ! 
  CALL mp_bcast (parallel_k, root) 
  CALL mp_bcast (parallel_q, root)
  CALL mp_bcast (a2f, root)
  CALL mp_bcast (epf_mem, root)
  CALL mp_bcast (etf_mem, root)
  CALL mp_bcast (rand_q, root)
  CALL mp_bcast (rand_k, root)
  CALL mp_bcast (mp_mesh_q, root)
  CALL mp_bcast (mp_mesh_k, root)
  CALL mp_bcast (wepexst, root)
  CALL mp_bcast (epexst, root)
  !
  ! integers
  !
  CALL mp_bcast (niter_ph, root)
  CALL mp_bcast (nmix_ph, root)
  CALL mp_bcast (maxirr, root)
  CALL mp_bcast (iverbosity, root)
  CALL mp_bcast (ngaussw, root)     ! FG
  CALL mp_bcast (nw, root)          ! 
  CALL mp_bcast (selfen_type, root) ! 
  CALL mp_bcast (nbndsub, root)     ! 
  CALL mp_bcast (nbndskip, root)    ! 
  CALL mp_bcast (nsmear, root)      ! 
  CALL mp_bcast (rand_nq, root)     ! 
  CALL mp_bcast (rand_nk, root)     ! 
  CALL mp_bcast (nkf1, root)
  CALL mp_bcast (nkf2, root)
  CALL mp_bcast (nkf3, root)
  CALL mp_bcast (nqf1, root)
  CALL mp_bcast (nqf2, root)
  CALL mp_bcast (nqf3, root)
  CALL mp_bcast (skf1, root) ! RM
  CALL mp_bcast (skf2, root) !
  CALL mp_bcast (skf3, root) !
  CALL mp_bcast (sqf1, root) !
  CALL mp_bcast (sqf2, root) !
  CALL mp_bcast (sqf3, root) !
  CALL mp_bcast (nqsmear, root ) ! RM
  CALL mp_bcast (nqstep, root) ! RM
  CALL mp_bcast (neptemp, root)
  !
  ! real*8
  !
  CALL mp_bcast (tr2_ph, root)
  CALL mp_bcast (amass, root)
  CALL mp_bcast (alpha_mix, root)
  CALL mp_bcast (xq, root)
  CALL mp_bcast (degaussw, root)  ! FG
  CALL mp_bcast (delta_smear, root)    ! 
  CALL mp_bcast (eminabs, root)    ! 
  CALL mp_bcast (emaxabs, root)    ! 
  CALL mp_bcast (deltaeabs, root)    ! 
  CALL mp_bcast (eps_acustic, root)  ! RM
  CALL mp_bcast (degaussq, root)  ! RM
  CALL mp_bcast (delta_qsmear, root) ! RM 
  !
  ! characters
  !
  CALL mp_bcast (title, root)
  CALL mp_bcast (filelph, root)
  CALL mp_bcast (fildvscf, root)
  CALL mp_bcast (fildrho, root)
  CALL mp_bcast (tmp_dir, root)
  CALL mp_bcast (prefix, root)
  !
  CALL mp_bcast (filkf, root)     ! FG
  CALL mp_bcast (filqf, root)     ! FG
  CALL mp_bcast (filukk, root)    ! FG
  CALL mp_bcast (filukq, root)    ! FG
  CALL mp_bcast (fileig, root)    ! FG
  CALL mp_bcast (fildvscf0, root) !
  CALL mp_bcast (dvscf_dir, root)
#endif
  !
END SUBROUTINE bcast_ph_input
!
!-----------------------------------------------------------------------
SUBROUTINE bcast_ph_input1
  !-----------------------------------------------------------------------
  !
#ifdef __PARA
#include "f_defs.h"

  USE pwcom
  USE phcom
  USE mp, ONLY: mp_bcast
  implicit none
  integer :: root = 0

  !
  ! integers
  !
  CALL mp_bcast (nat_todo, root)
  IF (nat_todo.gt.0) THEN
     CALL mp_bcast (atomo, root)
  ENDIF
#endif
  !  
END SUBROUTINE bcast_ph_input1
