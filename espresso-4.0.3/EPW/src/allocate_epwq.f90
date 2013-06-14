  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Code adapted from PH/allocate_phq - Quantum-ESPRESSO group
  ! 09/2009 There is a lot of excess in this file.  
  !
  !----------------------------------------------------------------------- 
  subroutine allocate_epwq
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the linear
  ! response problem
  !
#include "f_defs.h"

  USE ions_base, ONLY : nat, ntyp => nsp
  USE pwcom
  USE wavefunctions_module,  ONLY: evc
  USE kinds, only : DP
  USE phcom
  USE el_phon
  USE becmod, ONLY: becp
  USE uspp_param, ONLY: nhm
  implicit none
  !
  !  allocate space for the quantities needed in EPW
  !
  IF (lgamma) THEN
     !
     !  q=0  : evq and igkq are pointers to evc and igk
     !
     evq  => evc
     igkq => igk
  ELSE
     !
     !  q!=0 : evq, igkq are ALLOCATEd and calculated at point k+q
     !
     ALLOCATE (evq ( npwx , nbnd))    
     ALLOCATE (igkq ( npwx))    
  ENDIF
  !
  ALLOCATE (dvpsi ( npwx , nbnd))    
  ALLOCATE ( dpsi ( npwx , nbnd))    
  !
  ALLOCATE (vlocq ( ngm , ntyp))    
  ALLOCATE (dmuxc ( nrxx , nspin , nspin))    
  !
  ALLOCATE (eigqts ( nat))    
  ALLOCATE (rtau ( 3, 48, nat))    
  ALLOCATE (u ( 3 * nat, 3 * nat))    
  ALLOCATE (ubar ( 3 * nat))    
  ALLOCATE (w2 ( 3 * nat))    
  ALLOCATE (dyn00 ( 3 * nat, 3 * nat))
  ALLOCATE (t (max_irr_dim, max_irr_dim, 48,3 * nat))    
  ALLOCATE (tmq (max_irr_dim, max_irr_dim, 3 * nat))    
  ALLOCATE (npert ( 3 * nat))    
  IF (okvan) THEN
     ALLOCATE (int1 ( nhm, nhm, 3, nat, nspin))    
     ALLOCATE (int2 ( nhm , nhm , 3 , nat , nat))    
     ALLOCATE (int3 ( nhm , nhm , max_irr_dim , nat , nspin))    
     ALLOCATE (int4 ( nhm * (nhm + 1)/2,  3 , 3 , nat, nspin))    
     ALLOCATE (int5 ( nhm * (nhm + 1)/2 , 3 , 3 , nat , nat))    
     ALLOCATE (dpqq( nhm, nhm, 3, ntyp))    
     ALLOCATE (alphasum ( nhm * (nhm + 1)/2 , 3 , nat , nspin))    
     ALLOCATE (this_dvkb3_is_on_file(nks))    
     this_dvkb3_is_on_file(:)=.false.
  ENDIF
  ALLOCATE (this_pcxpsi_is_on_file(nks,3))
  this_pcxpsi_is_on_file(:,:)=.false.
  ALLOCATE ( alphap ( nkb , nbnd , 3 , nks))    
  ALLOCATE ( becp1 (nkb, nbnd, nks), becp(nkb, nbnd) )
  IF (elph) ALLOCATE (el_ph_mat( nbnd, nbnd, nks, 3*nat))    
END SUBROUTINE allocate_epwq
