!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------

subroutine tabd (nt, occ_loc)  
  !-----------------------------------------------------------------------
  !
  ! This routine is a table (far from being complete) for the total number
  ! of localized electrons in transition metals or rare earths
  ! (PPs usually are built on non physical configurations)
  !
  USE kinds, ONLY: DP
  USE uspp_param, ONLY: upf
  implicit none
  real(DP) :: occ_loc
  ! output: the total number of d electrons

  integer :: nt
  !
  ! TRANSITION METALS
  !
  if (upf(nt)%psd .eq.'Mn') then
     occ_loc = 5.d0
  elseif (upf(nt)%psd .eq.'Fe') then
     occ_loc = 6.d0
  elseif (upf(nt)%psd .eq.'Co') then
     occ_loc = 7.d0
  elseif (upf(nt)%psd .eq.'Ni') then
     occ_loc = 8.d0
  elseif (upf(nt)%psd .eq.'Cu') then
     occ_loc = 10.d0
     !
     ! RARE EARTHS
     !
  elseif (upf(nt)%psd .eq.'Ce') then
     occ_loc = 2.d0
     !
     ! OTHER ELEMENTS
     !
  elseif (upf(nt)%psd .eq.'C') then
     occ_loc = 2.d0
  elseif (upf(nt)%psd .eq.'O') then
     occ_loc = 4.d0
  elseif (upf(nt)%psd .eq.'H') then
     occ_loc = 1.d0
  else
     occ_loc = 0.d0
     call errore ('tabd', 'pseudopotential not yet inserted', 1)
  endif
  return
end subroutine tabd

