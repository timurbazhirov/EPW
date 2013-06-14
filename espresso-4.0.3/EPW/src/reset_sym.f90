  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine reset_sym (xq, nsym, s, invs, irt, rtau)
  !-----------------------------------------------------------------------
  !
  ! output values of symmetry arrays (nsym, s, rtau, irt) are those
  ! appropriate to the small-qroup of q.
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE io_global,     ONLY : stdout
  USE kinds,         only : DP
  USE phcom,         only : minus_q, irgq, nsymq, irotmq, t, tmq, max_irr_dim, &
       u, npert, nirr, gi, gimq
  USE pwcom,         only : ibrav, at, sname, symm_type, bg, nr1, nr2, nr3, ftau, &
                    invsym
  USE ions_base,     only : nat, tau, ityp
  USE control_flags, only :  iverbosity, noinv, modenum
  USE epwcom,        only : iswitch 
  !
  implicit none
  real(kind=DP) :: xq (3)
  ! input: q vector

  !-output variables
  integer :: nsym, s (3, 3, 48), invs (48), irt (48, nat)
  ! output: number of symmetry operations
  ! output: the first nq matrices are those that generate the star of q
  !         starting from it
  ! output: list of inverse operation indices
  ! output: for each atom gives the rotated atom

  real(kind=DP) :: rtau (3, 48, nat)
  ! output: for each atom and rotation gives the R vector involved
  ! output: list of vectors in the star of q
  !
  ! Local variables
  !
  integer :: nrot, isym, jsym, table (48, 48), &
       i, j, nks0, npk0, izero
  ! number of symmetry ops. of bravais lattice.
  ! counters on symmetry ops.
  ! index of inverse of isym
  ! group table
  ! counter on q-vectors
  ! generic counter
  ! number of dummy k-points
  ! maximum allowed number of dummy k-points
  ! dummy (zero) value of iswitch passed to sgama
  real(kind=DP) :: xk0 (3), wk0(1), zero (3), mdum(3,nat)
  ! auxiliary list of q (crystal coordinates)
  ! input q in crystal coordinates
  ! rotated q in crystal coordinates
  ! coordinates of fractionary translations
  ! dummy k-points list
  ! a zero vector: used in eqvect and as dummy q-vector in sgama

  logical :: nosym, sym (48)
  ! .t. if the crystal has inversion
  ! dummy output from sgama
  ! input for sgama

  logical, external :: eqvect
  ! function used to compare two vectors
  !
  !  initialize dummy k-point list and zero vector
  !
  izero = 0
  npk0 = 1
  nks0 = 1
  wk0(:) = 1.d0
  xk0(:)= 0.d0
  zero(:) = 0.d0
  !
  !  generate transformation matrices for the bravais lattice
  !
  IF (ibrav == 4 .or. ibrav == 5) then
     CALL hexsym (at, s, sname, nrot)
  ELSEif (ibrav >= 1 .and. ibrav <= 14) then
     CALL cubicsym (at, s, sname, nrot)
  ELSEif (ibrav == 0) then
     IF (symm_type == 'cubic') call cubicsym (at, s, sname, nrot)
     IF (symm_type == 'hexagonal') call hexsym (at, s, sname, nrot)
  ELSE
     CALL errore ('star_q', 'wrong ibrav', 1)
  ENDIF

!  if (noinv) then
  IF (noinv) then
     jsym = 0
     DO isym = 1, nrot
        IF ( s (1, 3, isym) == 0 .and. s (3, 1, isym) == 0 .and. &
             s (2, 3, isym) == 0 .and. s (3, 2, isym) == 0 .and. &
             s (3, 3, isym) == 1) then
           jsym = jsym + 1
           DO i = 1, 3
              DO j = 1, 3
                 s (i, j, jsym) = s (i, j, isym)
              ENDDO
           ENDDO
           sname (jsym) = sname (isym)
        ENDIF
     ENDDO
     nrot = jsym
  ENDIF
  !
  ! extract from it the crystal symmetry group by calling sgama
  !
  nosym = .false.
  CALL sgama (nrot, nat, s, sname, at, bg, tau, ityp, nsym, nr1, &
       nr2, nr3, irt, ftau, npk0, nks0, xk0, wk0, invsym, minus_q, xq, &
       iswitch, modenum, .false., mdum)
  DO isym = 1, nsym
     sym (isym) = .true.
  ENDDO
  CALL sgam_ph (at, bg, nsym, s, irt, tau, rtau, nat, sym)
  !
  ! computes the inverse of each matrix
  !
  CALL multable (nsym, s, table)
  DO isym = 1, nsym
     DO jsym = 1, nsym
        IF (table (isym, jsym) .eq.1) invs (isym) = jsym
     ENDDO
  ENDDO
  !
  IF (nsym.gt.1) then
     CALL set_irr (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
          nirr, gi, gimq, iverbosity)
  ELSE
     CALL set_irr_nosym (nat, at, bg, xq, s, invs, nsym, rtau, irt, &
          irgq, nsymq, minus_q, irotmq, t, tmq, max_irr_dim, u, npert, &
          nirr, gi, gimq, iverbosity)
  ENDIF
  !
  end subroutine reset_sym
