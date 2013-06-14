!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!---------------------------------------------------------------------
subroutine sym_dmag (nper, irr, dmagtosym)
  !---------------------------------------------------------------------
  ! symmetrize the change of the magnetization density
  ! belonging to an irreproducible representation
  !
  USE kinds, only : DP
  USE constants, ONLY: tpi
  USE gvect, ONLY: nr1, nr2, nr3, nrx1, nrx2, nrx3
  USE cell_base, ONLY : at, bg
  USE symme, ONLY : s, ftau, t_rev
  USE char,  ONLY : sname
  USE lsda_mod, ONLY: nspin
  use phcom
  implicit none

  integer :: nper, irr
  ! the number of perturbations
  ! the representation under conside

  complex(DP) :: dmagtosym (nrx1, nrx2, nrx3, nspin, nper)
  ! the magnetization to symmetrize (only 2:4 components)

  integer :: is, ri, rj, rk, i, j, k, ipert, jpert, ipol, isym, &
       irot, kpol
  !  counter on spin polarizations
  !
  !  the rotated points
  !
  !
  !  counter on mesh points
  !
  ! counter on perturbations
  ! counter on perturbations
  ! counter on polarizations
  ! counter on symmetries
  ! the rotation

  real(DP) :: g1 (48), g2 (48), g3 (48), in1, in2, in3
  ! used to construct the phases
  ! auxiliary variables

  complex(DP), allocatable :: dmagsym (:,:,:,:,:), dmags(:,:)
  ! the symmetrized potential
  complex(DP) ::  aux2(3), term (3, 48), phase (48), mag(3), magrot(3)
  ! auxiliary space
  ! the multiplication factor
  ! the phase factor

  if (nsymq == 1.and. (.not.minus_q) ) return
  call start_clock ('sym_dmag')

  allocate (dmagsym(  nrx1 , nrx2 , nrx3 , 3, nper))    
  allocate (dmags( 3, nper))    
  !
  ! if necessary we symmetrize with respect to  S(irotmq)*q = -q + Gi
  !
  in1 = tpi / DBLE (nr1)
  in2 = tpi / DBLE (nr2)
  in3 = tpi / DBLE (nr3)

  if (minus_q) then
     g1 (1) = 0.d0
     g2 (1) = 0.d0
     g3 (1) = 0.d0
     do ipol = 1, 3
        g1 (1) = g1 (1) + gimq (ipol) * in1 * at (ipol, 1)
        g2 (1) = g2 (1) + gimq (ipol) * in2 * at (ipol, 2)
        g3 (1) = g3 (1) + gimq (ipol) * in3 * at (ipol, 3)
     enddo
     term (1, 1) = CMPLX (cos (g1 (1) ), sin (g1 (1) ) )
     term (2, 1) = CMPLX (cos (g2 (1) ), sin (g2 (1) ) )
     term (3, 1) = CMPLX (cos (g3 (1) ), sin (g3 (1) ) )
     phase (1) = (1.d0, 0.d0)
     do k = 1, nr3
        do j = 1, nr2
           do i = 1, nr1
              ri = s (1, 1, irotmq) * (i - 1) + s (2, 1, irotmq) * (j - 1) &
                 + s (3, 1, irotmq) * (k - 1) - ftau (1, irotmq)
              ri = mod (ri, nr1) + 1
              if (ri < 1) ri = ri + nr1
              rj = s (1, 2, irotmq) * (i - 1) + s (2, 2, irotmq) * (j - 1) &
                 + s (3, 2, irotmq) * (k - 1) - ftau (2, irotmq)
              rj = mod (rj, nr2) + 1
              if (rj < 1) rj = rj + nr2
              rk = s (1, 3, irotmq) * (i - 1) + s (2, 3, irotmq) * (j - 1) &
                 + s (3, 3, irotmq) * (k - 1) - ftau (3, irotmq)
              rk = mod (rk, nr3) + 1
              if (rk < 1) rk = rk + nr3

              do ipert = 1, nper
                 aux2 = (0.d0, 0.d0)
                 do jpert = 1, nper
                    do is=2,4
                       aux2(is-1) = aux2(is-1) + tmq (jpert, ipert, irr) * &
                           dmagtosym (ri, rj, rk, is, jpert) * phase (1)
                    enddo
                 enddo
                 do kpol = 1, 3
                    mag(kpol)=bg(1,kpol)*aux2(1) + bg(2,kpol)*aux2(2) + &
                              bg(3,kpol)*aux2(3)
                 enddo
! rotate the magnetic moment
                 do kpol = 1, 3
                    magrot(kpol) = s(1,kpol,invs(irotmq))*mag(1) + &
                                   s(2,kpol,invs(irotmq))*mag(2) + &
                                   s(3,kpol,invs(irotmq))*mag(3)
                 enddo
                 if (sname(irotmq)(1:3)=='inv') magrot=-magrot
                 if(t_rev(irotmq).eq.1) magrot=-magrot
! go back to carthesian coordinates
                 do kpol = 1, 3
                    mag(kpol)=at(kpol,1)*magrot(1) + &
                              at(kpol,2)*magrot(2) + &
                              at(kpol,3)*magrot(3)
                    dmagsym(i,j,k,kpol,ipert)=(dmagtosym(i,j,k,kpol+1,ipert)+&
                         CONJG(mag(kpol)) ) * 0.5d0
                 enddo
              enddo
              phase (1) = phase (1) * term (1, 1)
           enddo
           phase (1) = phase (1) * term (2, 1)
        enddo
        phase (1) = phase (1) * term (3, 1)
     enddo
     do ipert = 1, nper
        do is=2,4
           dmagtosym(:, :, :, is, ipert) = dmagsym (:, :, :, is-1, ipert)
        end do
     enddo
  endif
  !
  ! Here we symmetrize with respect to the small group of q
  !
  do isym = 1, nsymq
     g1 (isym) = 0.d0
     g2 (isym) = 0.d0
     g3 (isym) = 0.d0
     do ipol = 1, 3
        g1 (isym) = g1 (isym) + gi (ipol, isym) * in1 * at (ipol, 1)
        g2 (isym) = g2 (isym) + gi (ipol, isym) * in2 * at (ipol, 2)
        g3 (isym) = g3 (isym) + gi (ipol, isym) * in3 * at (ipol, 3)
     enddo
     term (1, isym) = CMPLX (cos (g1 (isym) ), sin (g1 (isym) ) )
     term (2, isym) = CMPLX (cos (g2 (isym) ), sin (g2 (isym) ) )
     term (3, isym) = CMPLX (cos (g3 (isym) ), sin (g3 (isym) ) )
  enddo

  dmagsym(:,:,:,:,:) = (0.d0, 0.d0)
  do isym = 1, nsymq
     phase (isym) = (1.d0, 0.d0)
  enddo
  do k = 1, nr3
     do j = 1, nr2
        do i = 1, nr1
           do isym = 1, nsymq
              irot = irgq (isym)
              ri = s (1, 1, irot) * (i - 1) + s (2, 1, irot) * (j - 1) &
                 + s (3, 1, irot) * (k - 1) - ftau (1, irot)
              ri = mod (ri, nr1) + 1
              if (ri < 1) ri = ri + nr1
              rj = s (1, 2, irot) * (i - 1) + s (2, 2, irot) * (j - 1) &
                 + s (3, 2, irot) * (k - 1) - ftau (2, irot)
              rj = mod (rj, nr2) + 1
              if (rj < 1) rj = rj + nr2
              rk = s (1, 3, irot) * (i - 1) + s (2, 3, irot) * (j - 1) &
                 + s (3, 3, irot) * (k - 1) - ftau (3, irot)
              rk = mod (rk, nr3) + 1
              if (rk < 1) rk = rk + nr3
              dmags=(0.d0,0.d0)
              do ipert = 1, nper
                 do jpert = 1, nper
                    do is=2,4
                       dmags(is-1,ipert)=dmags(is-1,ipert) + &
                            t (jpert, ipert, irot, irr) * &
                            dmagtosym (ri, rj, rk, is, jpert) * phase (isym)
                    enddo
                 enddo
                 do kpol = 1, 3
                    mag(kpol)=bg(1,kpol)*dmags(1,ipert) + &
                              bg(2,kpol)*dmags(2,ipert) + &
                              bg(3,kpol)*dmags(3,ipert)
                 enddo
! rotate the magnetic moment
                 do kpol = 1, 3
                    magrot(kpol) = s(1,kpol,invs(irot))*mag(1) + &
                                   s(2,kpol,invs(irot))*mag(2) + &
                                   s(3,kpol,invs(irot))*mag(3)
                 enddo
                 if (sname(irot)(1:3)=='inv') magrot=-magrot
                 if(t_rev(irot).eq.1) magrot=-magrot
! go back to carthesian coordinates
                 do kpol = 1, 3
                    mag(kpol)=at(kpol,1)*magrot(1) + &
                              at(kpol,2)*magrot(2) + &
                              at(kpol,3)*magrot(3)
                 enddo
                 dmagsym(i,j,k,1,ipert)=dmagsym(i,j,k,1,ipert)+mag(1)
                 dmagsym(i,j,k,2,ipert)=dmagsym(i,j,k,2,ipert)+mag(2)
                 dmagsym(i,j,k,3,ipert)=dmagsym(i,j,k,3,ipert)+mag(3)
              enddo
           enddo
           do isym = 1, nsymq
              phase (isym) = phase (isym) * term (1, isym)
           enddo
        enddo
        do isym = 1, nsymq
           phase (isym) = phase (isym) * term (2, isym)
        enddo
     enddo
     do isym = 1, nsymq
        phase (isym) = phase (isym) * term (3, isym)
     enddo
  enddo

  do is=2,4
     do ipert = 1, nper
        dmagtosym(:,:,:,is,ipert) = dmagsym(:,:,:,is-1,ipert) / DBLE (nsymq)
     enddo
  enddo

  deallocate (dmags)
  deallocate (dmagsym)

  call stop_clock ('sym_dmag')
  return
end subroutine sym_dmag
