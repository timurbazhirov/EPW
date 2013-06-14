  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Code adapted from PH/deallocate_phq - Quantum-ESPRESSO group
  !
  !----------------------------------------------------------------------
  SUBROUTINE deallocate_epw
  !----------------------------------------------------------------------
  !
  !  deallocates the variables allocated by allocate_epw
  !  this routine is unchanged as of 3/9/09 and should be cleaned and fixed
  !  09/2009 Cleanup still necessary
  !  12/2009 Added variables from elph.f90
  !----------------------------------------------------------------------
  USE epwcom
  USE phcom
  USE becmod, ONLY: becp
  USE wavefunctions_module,  ONLY: evc
  USE el_phon

  !
  IF (lgamma) THEN
     IF(ASSOCIATED(evq)) NULLIFY(evq)
     IF(ASSOCIATED(igkq)) NULLIFY(igkq)
  ELSE
     IF(ASSOCIATED(evq)) DEALLOCATE(evq)
     IF(ASSOCIATED(igkq)) DEALLOCATE(igkq)
  END IF
  !
  IF(ALLOCATED(dvpsi)) DEALLOCATE (dvpsi)    
  IF(ALLOCATED(dpsi)) DEALLOCATE ( dpsi)    
  !
  IF(ALLOCATED(vlocq)) DEALLOCATE (vlocq)
  IF(ALLOCATED(dmuxc)) DEALLOCATE (dmuxc)
  !
  IF(ALLOCATED(eigqts)) DEALLOCATE (eigqts)
  IF(ALLOCATED(rtau)) DEALLOCATE (rtau)
  IF(ASSOCIATED(u)) DEALLOCATE (u)
  IF(ASSOCIATED(ubar)) DEALLOCATE (ubar)
  IF(ALLOCATED(dyn)) DEALLOCATE (dyn)
  IF(ALLOCATED(dyn00)) DEALLOCATE (dyn00)
  IF(ALLOCATED(w2)) DEALLOCATE (w2)
  IF(ASSOCIATED(t)) DEALLOCATE (t)
  IF(ASSOCIATED(tmq)) DEALLOCATE (tmq)
  !
  IF(ALLOCATED(npert)) DEALLOCATE (npert)    
  !
  IF(ALLOCATED(int1)) DEALLOCATE (int1)    
  IF(ALLOCATED(int2)) DEALLOCATE (int2)
  IF(ALLOCATED(int3)) DEALLOCATE (int3)
  IF(ALLOCATED(int4)) DEALLOCATE (int4)
  IF(ALLOCATED(int5)) DEALLOCATE (int5)
  IF(ALLOCATED(dpqq)) DEALLOCATE (dpqq)
  !
  !
  IF(ALLOCATED(alphap))    DEALLOCATE (alphap)    
  IF(ALLOCATED(becp1))     DEALLOCATE(becp1) 
  IF(ALLOCATED(becp))      DEALLOCATE(becp)
  !
  IF(ALLOCATED(drc)) DEALLOCATE(drc)
  !
  IF(ALLOCATED(dvxc_rr)) DEALLOCATE (dvxc_rr)    
  IF(ALLOCATED(dvxc_sr)) DEALLOCATE (dvxc_sr)    
  IF(ALLOCATED(dvxc_ss)) DEALLOCATE (dvxc_ss)    
  IF(ALLOCATED(dvxc_s)) DEALLOCATE (dvxc_s)    
  IF(ALLOCATED(grho)) DEALLOCATE (grho)  
  !
  !  EPW variables
  !
  IF(ALLOCATED(el_ph_mat)) DEALLOCATE (el_ph_mat)    
  IF(ALLOCATED(epmatw17))  DEALLOCATE (epmatw17)    
  IF(ALLOCATED(epf17))     DEALLOCATE (epf17)    
  IF(ALLOCATED(etfq))      DEALLOCATE (etfq)    
  IF(ALLOCATED(etq))       DEALLOCATE (etq)    
  IF(ALLOCATED(etf))       DEALLOCATE (etf)    
  IF(ALLOCATED(wf))        DEALLOCATE (wf)    
  IF(ALLOCATED(xkq))       DEALLOCATE (xkq)    
  IF(ALLOCATED(xkf))       DEALLOCATE (xkf)    
  IF(ALLOCATED(wkf))       DEALLOCATE (wkf)    
  IF(ALLOCATED(xqf))       DEALLOCATE (xqf)    
  IF(ALLOCATED(wqf))       DEALLOCATE (wqf)    
  IF(ALLOCATED(xk_all))    DEALLOCATE (xk_all)    
  IF(ALLOCATED(et_all))    DEALLOCATE (et_all)    
  IF(ALLOCATED(wslen))     DEALLOCATE (wslen)    
  !
  END SUBROUTINE DEALLOCATE_epw
