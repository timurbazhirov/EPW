!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine sgam_at_mag (nrot, s, nat, tau, ityp, at, bg, nr1, nr2, &
     nr3, sym, irt, ftau, m_loc, sname, t_rev)
  !-----------------------------------------------------------------------
  !
  !     given a point group, this routine finds the subgroup which is
  !     the point group of the crystal under consideration
  !     non symmorphic groups non allowed, provided that fractional
  !     translations are commensurate with the FFT grid
  !
  !     It sets the array sym, which for each operation of the original
  !     point group is true if this operation is also an operation of the
  !     total point group
  !
#include "f_defs.h"
  USE io_global,  ONLY : stdout
  USE kinds
  implicit none
  !
  !     input variables
  !
  integer :: nrot, s (3, 3, 48), nat, ityp (nat), nr1, nr2, nr3
  REAL(DP) :: m_loc(3,nat), tau (3, nat), at (3, 3), bg (3, 3)
  ! nrot : order of the parent group
  ! s    : symmetry operations of parent group
  ! nat  : number of atoms in the unit cell
  ! ityp : species of each atom in the unit cell
  ! nr*  : dimensions of the FFT mesh
  ! tau  : cartesian coordinates of the atoms
  ! at   : basis of the real-space lattice
  ! bg   :  "   "   "  reciprocal-space lattice
  !
  !     output variables
  !
  integer :: irt (48, nat), ftau (3, 48)
  logical :: sym (48)
  ! irt(isym,na) : sym.op. isym sends atom na into atom irt(isym,na)
  ! ftau(:,isym) : fractional translation associated to sym.op. isym
  !                (in FFT coordinates: crystal axis, multiplied by nr*)
  ! sym(isym)    : flag indicating if sym.op. isym in the parent group
  !                is a true symmetry operation of the crystal
  !
  !    local variables
  !
  integer :: na, kpol, nb, irot, i, j
  ! counters
  real(DP) , allocatable :: xau (:,:), rau (:,:), mxau(:,:), mrau(:,:)
  ! atomic coordinates in crystal axis
  logical :: fractional_translations
  real(DP) :: ft (3), ft1, ft2, ft3
  integer :: t_rev(48)
  character :: sname (48) * 45
  !
  external checksym
  !
  allocate(xau(3,nat))
  allocate(rau(3,nat))
  ALLOCATE(mxau(3,nat))
  ALLOCATE(mrau(3,nat))
  !
  !     Compute the coordinates of each atom in the basis of
  !     the direct lattice vectors
  !
  do na = 1, nat
     do kpol = 1, 3
        xau (kpol, na) = bg (1, kpol) * tau (1, na) + &
                         bg (2, kpol) * tau (2, na) + &
                         bg (3, kpol) * tau (3, na)
        mxau (kpol, na)= bg (1, kpol) * m_loc (1, na) + &
                         bg (2, kpol) * m_loc (2, na) + &
                         bg (3, kpol) * m_loc (3, na)
     enddo
  enddo
  !
  !      check if the identity has fractional translations
  !      (this means that the cell is actually a supercell).
  !      When this happens, fractional translations are disabled,
  !      because there is no guarantee that the generated sym.ops.
  !      form a group
  !
  nb = 1
  irot = 1
  fractional_translations = .true.
  do na = 2, nat
     if (ityp (nb) .eq.ityp (na) ) then
        ft (1) = xau(1,na) - xau(1,nb) - nint( xau(1,na) - xau(1,nb) )
        ft (2) = xau(2,na) - xau(2,nb) - nint( xau(2,na) - xau(2,nb) )
        ft (3) = xau(3,na) - xau(3,nb) - nint( xau(3,na) - xau(3,nb) )


        call checksym_mag (irot, nat, ityp, xau, xau, ft, sym, irt, mxau,&
                           mxau, t_rev(irot))

        if (sym (irot) .and. (abs (ft (1) **2 + ft (2) **2 + ft (3) ** &
             2) ) .lt.1.d-8) call errore ('sgam_at', 'overlapping atoms', na)
        if (sym (irot) ) then
           fractional_translations = .false.
           WRITE( stdout, '(5x,"Found additional translation:",3f10.4)') ft
        endif
     endif

  enddo
  do irot = 1, nrot
     !
     ! check that the grid is compatible with the S rotation
     !
     if ( mod (s (2, 1, irot) * nr1, nr2) .ne.0 .or. &
          mod (s (3, 1, irot) * nr1, nr3) .ne.0 .or. &
          mod (s (1, 2, irot) * nr2, nr1) .ne.0 .or. &
          mod (s (3, 2, irot) * nr2, nr3) .ne.0 .or. &
          mod (s (1, 3, irot) * nr3, nr1) .ne.0 .or. &
          mod (s (2, 3, irot) * nr3, nr2) .ne.0 ) then
        sym (irot) = .false.
        WRITE( stdout, '(5x,"warning: symmetry operation # ",i2, &
             &         " not compatible with FFT grid. ")') irot
        WRITE( stdout, '(3i4)') ( (s (i, j, irot) , j = 1, 3) , i = 1, 3)
        goto 100

     endif
     do na = 1, nat
        do kpol = 1, 3
           ! rau = rotated atom coordinates
           rau (kpol, na) = s (1, kpol, irot) * xau (1, na) + &
                            s (2, kpol, irot) * xau (2, na) + &
                            s (3, kpol, irot) * xau (3, na)
           mrau(kpol,na) = s(1,kpol,irot)*mxau(1,na) + & ! magnetic moment of the
                           s(2,kpol,irot)*mxau(2,na) + & ! atom rotated by the
                           s(3,kpol,irot)*mxau(3,na)   ! present symmerty
                                                       ! operation
        enddo
     enddo
     if (sname(irot)(1:3)=='inv') mrau=-mrau
     !
     !      first attempt: no fractional translation
     !
     do kpol = 1, 3
        ftau (kpol, irot) = 0
        ! input for checksym
        ft (kpol) = 0.d0
     enddo

     call checksym_mag (irot, nat, ityp, xau, rau, ft, sym, irt, mxau, &
                        mrau, t_rev(irot))
     if (.not.sym (irot) .and.fractional_translations) then
        nb = 1
        do na = 1, nat
           if (ityp (nb) .eq.ityp (na) ) then
              !
              !      second attempt: check all possible fractional translations
              !
              ft (1) = rau(1,na) - xau(1,nb) - nint( rau(1,na) - xau(1,nb) )
              ft (2) = rau(2,na) - xau(2,nb) - nint( rau(2,na) - xau(2,nb) )
              ft (3) = rau(3,na) - xau(3,nb) - nint( rau(3,na) - xau(3,nb) )

              call checksym_mag (irot,nat,ityp,xau,rau,ft,sym,irt,mxau, &
                                 mrau,t_rev(irot))
              if (sym (irot) ) then
                 ! convert ft to FFT coordinates
                 ! for later use in symmetrization
                 ft1 = ft (1) * nr1
                 ft2 = ft (2) * nr2
                 ft3 = ft (3) * nr3
                 ! check if the fractional translations are commensurate
                 ! with the FFT grid, discard sym.op. if not
                 if (abs (ft1 - nint (ft1) ) / nr1.gt.1.0d-5 .or. &
                     abs (ft2 - nint (ft2) ) / nr2.gt.1.0d-5 .or. &
                     abs (ft3 - nint (ft3) ) / nr3.gt.1.0d-5) then
                    WRITE( stdout, '(5x,"warning: symmetry operation", &
                         &     " # ",i2," not allowed.   fractional ", &
                         &     "translation:"/5x,3f11.7,"  in crystal", &
                         &     " coordinates")') irot, ft
                    sym (irot) = .false.
                 endif
                 ftau (1, irot) = nint (ft1)
                 ftau (2, irot) = nint (ft2)
                 ftau (3, irot) = nint (ft3)
                 goto 100
              endif
           endif
        enddo

     endif
100  continue
  enddo
  !
  !   deallocate work space
  !
  deallocate (rau)
  deallocate (xau)
  !
  return
END SUBROUTINE sgam_at_mag

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------
subroutine checksym_mag (ir,nat,ityp,xau,rau,ft,sym,irt,mxau,mrau,t_rev)
  !-----------------------------------------------------------------------
  !
  !   This routine receives as input all the atomic positions xau,
  !   and the rotated rau by the symmetry operation ir. It sets to true
  !   sym(ir) if for each atom na, it is possible to find an atom nb
  !   which is of the same type of na, and coincide with it after the
  !   symmetry operation. Fractional translations are allowed.
  !
  !   Revised layout 1 may 1995 by A. Dal Corso
  !
  USE kinds
  implicit none
  !
  !     first the dummy variables
  !
  integer :: nat, ityp (nat), irt (48, nat), ir, t_rev
  ! input: the total number of atoms
  ! input: the type of each atom
  ! output: the rotated of each atom
  ! input: the rotation to be tested
  ! output: time reversal operation is present
  real(DP) :: xau (3, nat), rau (3, nat), ft (3)
  ! input: the initial vectors
  ! input: the rotated vectors
  ! input: the possible fractionary translat
  REAL(DP) :: mxau(3,nat), mrau(3,nat)
  ! input: the rotated vectors
  ! input: the possible fractionary translation
  REAL(DP), PARAMETER :: mt(3) = (/0,0,0/)
  logical :: sym (48), sym_xyz
  ! output: if true this is a symmetry opera
  !
  !  few local variables
  !
  integer :: na, nb, na1, t1, t2
  ! counter on atoms
  ! counter on atoms
  logical :: eqvect
  ! the testing function

  external eqvect

  t1 = 1
  t2 = 1

  do na = 1, nat
     na1 = 0
     do nb = 1, nat
        if(ityp(na).eq.ityp(nb).and. &
           eqvect(rau (1, na),xau(1,nb),ft)) na1 = nb
     enddo
     !
     IF ( na1 /= 0 ) THEN
        !
        if( abs(mrau(1,na) - mxau(1,na1))+       &
            abs(mrau(2,na) - mxau(2,na1))+       &
            abs(mrau(3,na) - mxau(3,na1)).gt.1.0D-5) t1 = 0
        if( abs(mrau(1,na) + mxau(1,na1))+       &
            abs(mrau(2,na) + mxau(2,na1))+       &
            abs(mrau(3,na) + mxau(3,na1)).gt.1.0D-5) t2 = 0
        !
     END IF
     !
     if(na1.eq.0.or.(t1+t2).eq.0) then
       sym(ir) = .false.
       t_rev = 0
       return
     else
       irt(ir,na) = na1
     endif
  enddo

  if(t1+t2.eq.2) then
    sym(ir) = .true.
    t_rev = 0
  elseif(t1.eq.1) then
    sym(ir) = .true.
    t_rev = 0
  else
    sym(ir) = .true.
    t_rev = 1
  endif

  return
END SUBROUTINE checksym_mag

