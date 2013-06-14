!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from PH/dvqpsi_us_only
!
!----------------------------------------------------------------------
subroutine dvqpsi_us_only3 (ik, mode, uact, xxk)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector uact.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp, ntyp => nsp
  use pwcom
  USE kinds, only : DP
  USE uspp_param, ONLY: nh
  use phcom
  implicit none
  !
  !   The dummy variables
  !

  integer :: ik, mode
  ! input: the k point
  ! input: the actual perturbation
  real(kind=DP) :: xxk(3)
  ! input: the k point (cartesian coordinates)
  complex(kind=DP) :: uact (3 * nat)
  ! input: the pattern of displacements
  !
  !   And the local variables
  !

  integer :: na, nb, mu, nu, ig, igg, nt, ibnd, ijkb0, &
       ikb, jkb, ih, jh, ipol
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! the point k+q
  ! counter on G vectors
  ! auxiliary counter on G vectors
  ! counter on atomic types
  ! counter on bands
  ! auxiliary variable for counting
  ! counter on becp functions
  ! counter on becp functions
  ! counter on n index
  ! counter on m index
  ! counter on polarizations

  real(kind=DP), parameter :: eps = 1.d-12

  complex(kind=DP), allocatable :: ps1 (:,:), ps2 (:,:,:), aux (:)
  ! work space

  logical :: ok

  CALL start_clock ('dvqpsi_us_on')
  allocate (ps1 ( nkb , nbnd))    
  allocate (ps2 ( nkb , nbnd , 3))    
  allocate (aux ( npwx))    
  IF (lsda) current_spin = isk (ik)
  !
  !   we first compute the coefficients of the vectors
  !
  ps1(:,:)   = (0.d0, 0.d0)
  ps2(:,:,:) = (0.d0, 0.d0)
  ijkb0 = 0
  DO nt = 1, ntyp
     DO na = 1, nat
        IF (ityp (na) .eq.nt) then
           mu = 3 * (na - 1)
           DO ih = 1, nh (nt)
              ikb = ijkb0 + ih
              DO jh = 1, nh (nt)
                 jkb = ijkb0 + jh
                 DO ipol = 1, 3
                    DO ibnd = 1, nbnd
                       IF ( abs (uact (mu + 1) ) + &
                            abs (uact (mu + 2) ) + &
                            abs (uact (mu + 3) ) > eps) then
                          ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                               (deeq (ih, jh, na, current_spin) - &
                                et (ibnd, ik) * qq (ih, jh, nt) ) * &
                                alphap(jkb, ibnd, ipol, ik) * uact (mu + ipol)
                          ps2 (ikb, ibnd, ipol) = ps2 (ikb, ibnd, ipol) +&
                               (deeq (ih,jh, na, current_spin) - &
                                et (ibnd, ik) * qq (ih, jh, nt) ) * &
                                (0.d0, -1.d0) * becp1 (jkb, ibnd, ik) * &
                                uact (mu + ipol) * tpiba
                          IF (okvan) then
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                  (int1 (ih, jh, ipol,na, current_spin) * &
                                  becp1 (jkb, ibnd, ik) ) * uact (mu +ipol)
                          ENDIF
                       ENDIF
                       IF (okvan) then
                          DO nb = 1, nat
                             nu = 3 * (nb - 1)
                             ps1 (ikb, ibnd) = ps1 (ikb, ibnd) + &
                                  (int2 (ih, jh, ipol, nb, na) * &
                                   becp1 (jkb, ibnd, ik) ) * uact (nu + ipol)
                          ENDDO
                       ENDIF
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           ijkb0 = ijkb0 + nh (nt)
        ENDIF
     ENDDO
  ENDDO
  !
  !      This term is proportional to beta(k+q+G)
  !
  IF (nkb.gt.0) call ZGEMM ('N', 'N', npwq, nbnd, nkb, &
       (1.d0, 0.d0) , vkb, npwx, ps1, nkb, (1.d0, 0.d0) , dvpsi, npwx)
  !
  !      This term is proportional to (k+q+G)_\alpha*beta(k+q+G)
  !
  DO ikb = 1, nkb
     DO ipol = 1, 3
        ok = .false.
        DO ibnd = 1, nbnd
           ok = ok.or. (abs (ps2 (ikb, ibnd, ipol) ) .gt.eps)
        ENDDO
        IF (ok) then
           DO ig = 1, npwq
              igg = igkq (ig)
              aux (ig) =  vkb(ig, ikb) * (xxk(ipol) + g(ipol, igg) )
           ENDDO
           DO ibnd = 1, nbnd
              CALL ZAXPY (npwq, ps2(ikb,ibnd,ipol), aux, 1, dvpsi(1,ibnd), 1)
           ENDDO
        ENDIF
     ENDDO

  ENDDO
  deallocate (aux)
  deallocate (ps2)
  deallocate (ps1)

  CALL stop_clock ('dvqpsi_us_on')

end subroutine dvqpsi_us_only3
