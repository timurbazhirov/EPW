!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine dielec
  !-----------------------------------------------------------------------
  !
  !      calculates the dielectric tensor
  !
#include "f_defs.h"

  USE io_global,  ONLY : stdout
  USE io_files, ONLY: iunigk
  use pwcom
  USE noncollin_module, ONLY : npol
  USE kinds, only : DP
  use phcom
  USE mp_global,        ONLY : inter_pool_comm, intra_pool_comm
  USE mp,               ONLY : mp_sum

  implicit none

  integer :: ibnd, ipol, jpol, nrec, ik
  ! counter on polarizations
  ! counter on records
  ! counter on k points
  real(DP) :: w, weight

  complex(DP) :: ZDOTC

  call start_clock ('dielec')
  epsilon(:,:) = 0.d0
  if (nksq > 1) rewind (unit = iunigk)
  do ik = 1, nksq
     if (nksq > 1) read (iunigk) npw, igk
     weight = wk (ik)
     w = fpi * weight / omega
     do ipol = 1, 3
        nrec = (ipol - 1) * nksq + ik
        call davcio (dvpsi, lrebar, iuebar, nrec, - 1)
        do jpol = 1, 3
           nrec = (jpol - 1) * nksq + ik
           call davcio (dpsi, lrdwf, iudwf, nrec, - 1)
           do ibnd = 1, nbnd_occ (ik)
              !
              !  this is the real part of <DeltaV*psi(E)|DeltaPsi(E)>
              !
              epsilon(ipol,jpol)=epsilon(ipol,jpol)-4.d0*w* DBLE( &
                   ZDOTC(npwx*npol, dvpsi (1, ibnd), 1, dpsi (1, ibnd), 1))
           enddo
        enddo
     enddo
  enddo
#ifdef __PARA
  call mp_sum ( epsilon, intra_pool_comm )
  call mp_sum ( epsilon, inter_pool_comm )
#endif
  !
  !      symmetrize
  !
  !       WRITE( stdout,'(/,10x,"Unsymmetrized in crystal axis ",/)')
  !       WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)

  call symtns (epsilon, nsym, s)
  !
  !    pass to cartesian axis
  !
  !      WRITE( stdout,'(/,10x,"Symmetrized in crystal axis ",/)')
  !      WRITE( stdout,'(10x,"(",3f15.5," )")') ((epsilon(ipol,jpol),
  !     +                                ipol=1,3),jpol=1,3)
  call trntns (epsilon, at, bg, 1)
  !
  ! add the diagonal part
  !
  do ipol = 1, 3
     epsilon (ipol, ipol) = epsilon (ipol, ipol) + 1.d0
  enddo
  !
  !  and print the result
  !
  IF (lnoloc) THEN 
      WRITE( stdout, '(/,10x,"Dielectric constant in cartesian axis (DV_Hxc=0)",/)')
  ELSE IF (lrpa) THEN
     WRITE( stdout, '(/,10x,"RPA dielectric constant in cartesian axis (DV_xc=0)",/)')
  ELSE
      WRITE( stdout, '(/,10x,"Dielectric constant in cartesian axis ",/)')
  ENDIF

 

  WRITE( stdout, '(10x,"(",3f18.9," )")') ((epsilon(ipol,jpol), ipol=1,3), jpol=1,3)
  call stop_clock ('dielec')

  return
end subroutine dielec
