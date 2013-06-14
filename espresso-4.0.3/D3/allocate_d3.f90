!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
subroutine allocate_d3
  !-----------------------------------------------------------------------
  !
  ! dynamical allocation of arrays: quantities needed for the third
  ! derivative of the total energy
  !
  USE ions_base,  ONLY : nat, ntyp => nsp
  use pwcom
  use phcom
  use d3com
  implicit none

  call allocate_phq
  if (lgamma) then
     vlocg0 => vlocq
     npertg0=> npert
     vkb0 => vkb
     ug0 => u
     tg0 => t
  else
     allocate (vlocg0( ngm, ntyp))    
     allocate (ug0( 3*nat, 3*nat))    
     allocate (tg0( max_irr_dim, max_irr_dim, 48, 3*nat))    
     allocate (npertg0( 3*nat))    
     allocate (vkb0( npwx , nkb))    
  endif
  allocate (psidqvpsi( nbnd, nbnd))    
  allocate (d3dyn( 3 * nat, 3 * nat, 3 * nat))    

  if (degauss.ne.0.d0) allocate (ef_sh( 3 * nat))    
  allocate (d3dyn_aux1 ( 3 * nat, 3 * nat, 3 * nat))    
  allocate (d3dyn_aux2 ( 3 * nat, 3 * nat, 3 * nat))    
  allocate (d3dyn_aux3 ( 3 * nat, 3 * nat, 3 * nat))    
  allocate (d3dyn_aux4 ( 3 * nat, 3 * nat, 3 * nat))    
  allocate (d3dyn_aux5 ( 3 * nat, 3 * nat, 3 * nat))    
  allocate (d3dyn_aux6 ( 3 * nat, 3 * nat, 3 * nat))    
  allocate (d3dyn_aux7 ( 3 * nat, 3 * nat, 3 * nat))    
  allocate (d3dyn_aux8 ( 3 * nat, 3 * nat, 3 * nat))    
  allocate (d3dyn_aux9 ( 3 * nat, 3 * nat, 3 * nat))    
  d3dyn_aux1 (:,:,:) = (0.d0, 0.d0)
  d3dyn_aux2 (:,:,:) = (0.d0, 0.d0)
  d3dyn_aux3 (:,:,:) = (0.d0, 0.d0)
  d3dyn_aux4 (:,:,:) = (0.d0, 0.d0)
  d3dyn_aux5 (:,:,:) = (0.d0, 0.d0)
  d3dyn_aux6 (:,:,:) = (0.d0, 0.d0)
  d3dyn_aux7 (:,:,:) = (0.d0, 0.d0)
  d3dyn_aux8 (:,:,:) = (0.d0, 0.d0)
  d3dyn_aux9 (:,:,:) = (0.d0, 0.d0)

  return
end subroutine allocate_d3
