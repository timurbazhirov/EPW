!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine dielec_test
  !-----------------------------------------------------------------------
  !
  ! Calculates the dielectric tensor using the finite-differences-derivative
  ! of the wavefunctions. This should be used only for testing purposes
  ! while doing a raman calculation
  !
#include "f_defs.h"
  use kinds, only : DP
  use pwcom
  USE io_files, ONLY: iunigk
  USE wavefunctions_module,  ONLY: evc
  use phcom
  USE ramanm
  USE mp_global,            ONLY : inter_pool_comm, intra_pool_comm
  USE mp,                   ONLY : mp_sum

  implicit none

  integer :: ibnd, ipol, jpol, nrec, ik, i1, i2
  real(DP) :: w_, weight, tmp
  complex(DP) :: ZDOTC

  epsilon (:,:) = 0.d0
  if (nksq > 1) rewind (unit=iunigk)
  do ik = 1, nksq
     if (nksq > 1) read (iunigk) npw, igk
     weight = wk (ik)
     w_ = - fpi * weight / omega
     call davcio (evc, lrwfc, iuwfc, ik, -1)
     do ipol = 1, 6
        nrec = (ipol - 1) * nksq + ik
        call davcio (dpsi, lrd2w, iud2w, nrec, -1)
        tmp = 0.d0
        do ibnd = 1, nbnd_occ (ik)
           tmp = tmp + 2.0d0 * w_ *                        &
              real (ZDOTC (npw, evc (1, ibnd), 1, dpsi (1, ibnd), 1))
        enddo
        i1 = a1j (ipol)
        i2 = a2j (ipol)
        epsilon (i1, i2) = epsilon (i1, i2) + tmp
        if (i1.ne.i2 ) epsilon (i2, i1) = epsilon (i2, i1) + tmp
     enddo
  enddo
#ifdef __PARA
  call mp_sum ( epsilon, intra_pool_comm )
  call mp_sum ( epsilon, inter_pool_comm )
#endif
  !
  !  symmetrize
  !
!  write(6,'(/,10x,''Unsymmetrized in crystal axis '',/)')
!  write(6,'(10x,''('',3f15.5,'' )'')') ((epsilon(ipol,jpol), &
!                                       ipol=1,3),jpol=1,3)

  call symtns(epsilon,nsym,s)
  !
  !  pass to cartesian axis
  !
!  write(6,'(/,10x,''Symmetrized in crystal axis '',/)')
!  write(6,'(10x,''('',3f15.5,'' )'')') ((epsilon(ipol,jpol), &
!                                  ipol=1,3),jpol=1,3)

  call trntns(epsilon,at,bg,1)
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     epsilon (ipol, ipol) = epsilon (ipol, ipol) + 1.d0
  end do
  !
  !  and print the result
  !
  write(6,'(/,10x,''Dielectric constant from finite-differences'',/)')
  write(6,'(10x,''('',3f18.9,'' )'')') ((epsilon(ipol,jpol),       &
                                  ipol=1,3),jpol=1,3)

  return
end subroutine dielec_test

