  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/star_q - Quantum-ESPRESSO group   
  !
  !-----------------------------------------------------------------------
  subroutine star_q3 ( xq, at, bg, nsym, s, invs, nq, sxq, isq, imq)
  !-----------------------------------------------------------------------
  !
  ! generate the star of q vectors that are equivalent to the input one
  ! and return their list along with the symmetry ops. needed to obtain
  ! them. symmetry arrays (nsym, s) (in input) are those
  ! appropriate to the crystal symmetry (not to the small-qroup of q).
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE io_global,  ONLY : stdout
  USE kinds, only : DP
  implicit none
  !
  ! input variables
  !
  real(kind=DP) :: xq (3), at (3, 3), bg (3, 3)
  ! q vector
  ! direct lattice vectors
  ! reciprocal lattice vectors
  integer :: nsym, s (3, 3, 48), invs (48), nq, isq (48), imq
  ! number of symmetry operations
  ! the symmetry operations
  ! list of inverse operation indices
  ! degeneracy of the star of q
  ! index of  q in the star for a given sym
  ! index of -q in the star (0 if not present)
  !
  ! output variables
  !
  real(kind=DP) :: sxq (3, 48)
  ! list of vectors in the star of q
  !
  ! Local variables
  !
  integer :: nsq (48), isym, iq, i
  ! number of symmetry ops. of bravais lattice.
  ! counters on symmetry ops.
  ! counter on q-vectors
  ! generic counter
  real(kind=DP) :: rq (3), zero (3)
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! a zero vector: used in eqvect and as dummy q-vector in sgama
  logical, external :: eqvect
  ! function used to compare two vectors
  zero(:) = 0.d0
  !
  ! go to crystal coordinates, rotate, and back to cartesian
  !
  CALL cryst_to_cart ( 1, xq, at, -1)
  DO iq = 1, nq
    CALL irotate ( xq, s(:,:,invs(isq(iq))), sxq(:,iq) )
  ENDDO
  CALL cryst_to_cart ( nq, sxq(:,1:nq), bg, 1)
  CALL cryst_to_cart (  1, xq, bg, 1)
  !
  WRITE ( stdout, * )
  WRITE ( stdout, '(7x,i4,3f14.9)') (iq, (sxq(i,iq), i=1,3), iq=1,nq)
  IF (imq.eq.0) write ( stdout, '(7x,i4,3f12.9)') (iq, (-sxq(i,iq), i=1,3), iq=1,nq)
  !
  end subroutine star_q3

