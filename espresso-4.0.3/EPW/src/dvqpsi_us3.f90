
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! adapted from PH/dvqpsi_us
!
!----------------------------------------------------------------------
SUBROUTINE dvqpsi_us3 (ik, mode, uact, addnlcc, xxk, xq0)
  !----------------------------------------------------------------------
  !
  ! This routine calculates dV_bare/dtau * psi for one perturbation
  ! with a given q. The displacements are described by a vector u.
  ! The result is stored in dvpsi. The routine is called for each k point
  ! and for each pattern u. It computes simultaneously all the bands.
  !
#include "f_defs.h"
  !
  USE ions_base, ONLY : nat, ityp
  USE pwcom
  ! use atom, only : nlcc
  USE uspp_param, ONLY: upf  
  USE wavefunctions_module,  ONLY: evc
  USE kinds, ONLY : DP
  USE phcom
  implicit none
  !
  !   The dummy variables
  !

  integer :: ik, mode
  ! input: the k point
  ! input: the actual perturbation
  real(kind=DP) :: xxk(3), xq0(3)
  ! input: the k point (cartesian coordinates)   
  ! the first q point of the star (cartesian)
  ! NB: in previous versions I declared xxk as integer (...)
  ! but this was not producing bad results since we only pass
  ! the location to dvqpsi_us_only3 where we declare it as real
  ! (in practice we could even define it as string in this sub...)
  complex(kind=DP) :: uact (3 * nat)
  ! input: the pattern of displacements
  logical :: addnlcc
  !
  !   And the local variables
  !

  integer :: na, mu, ikk, ig, nt, ibnd, ir, is
  ! counter on atoms
  ! counter on modes
  ! the point k
  ! counter on G vectors
  ! the type of atom
  ! counter on bands
  ! counter on real mesh

  complex(kind=DP) :: gtau, gu, fact, u1, u2, u3, gu0
  complex(kind=DP) , ALLOCATABLE, TARGET :: aux (:)
  complex(kind=DP) , ALLOCATABLE :: aux1 (:), aux2 (:)
  complex(kind=DP) , POINTER :: auxs (:)
  ! work space

  CALL start_clock ('dvqpsi_us')
  IF (nlcc_any) THEN
     allocate (aux( nrxx))    
     IF (doublegrid) THEN
        ALLOCATE (auxs( nrxxs))    
     ELSE
        auxs => aux
     ENDIF
  ENDIF
  ALLOCATE (aux1( nrxxs))    
  ALLOCATE (aux2( nrxxs))    
  !
  !    We start by computing the contribution of the local potential.
  !    The computation of the derivative of the local potential is done in
  !    reciprocal space while the product with the wavefunction is done in
  !    real space
  !
  IF (lgamma) THEN
     ikk = ik
  ELSE
     ikk = 2 * ik - 1
  ENDIF
  dvpsi(:,:) = (0.d0, 0.d0)
  aux1(:) = (0.d0, 0.d0)
  DO na = 1, nat
     fact = tpiba * (0.d0, -1.d0) * eigqts (na)
     mu = 3 * (na - 1)
     IF (abs (uact (mu + 1) ) + abs (uact (mu + 2) ) + abs (uact (mu + &
          3) ) .gt.1.0d-12) THEN
        nt = ityp (na)
        u1 = uact (mu + 1)
        u2 = uact (mu + 2)
        u3 = uact (mu + 3)
        gu0 = xq0 (1) * u1 + xq0 (2) * u2 + xq0 (3) * u3
        DO ig = 1, ngms
           gtau = eigts1 (ig1 (ig), na) * eigts2 (ig2 (ig), na) * eigts3 ( &
                ig3 (ig), na)
           gu = gu0 + g (1, ig) * u1 + g (2, ig) * u2 + g (3, ig) * u3
           aux1 (nls (ig) ) = aux1 (nls (ig) ) + vlocq (ig, nt) * gu * &
                fact * gtau
        ENDDO
     ENDIF
  ENDDO
  !
  ! add NLCC when present
  !
   IF (nlcc_any.and.addnlcc) THEN
      aux(:) = (0.d0, 0.d0)
      DO na = 1,nat
         fact = tpiba*(0.d0,-1.d0)*eigqts(na)
         mu = 3*(na-1)
         IF (abs(uact(mu+1))+abs(uact(mu+2))  &
                         +abs(uact(mu+3)).gt.1.0d-12) then
            nt=ityp(na)
            u1 = uact(mu+1)
            u2 = uact(mu+2)
            u3 = uact(mu+3)
            gu0 = xq0(1)*u1 +xq0(2)*u2+xq0(3)*u3
            IF (upf(nt)%nlcc) THEN
               DO ig = 1,ngm
                  gtau = eigts1(ig1(ig),na)*   &
                         eigts2(ig2(ig),na)*   &
                         eigts3(ig3(ig),na)
                  gu = gu0+g(1,ig)*u1+g(2,ig)*u2+g(3,ig)*u3
                  aux(nl(ig))=aux(nl(ig))+drc(ig,nt)*gu*fact*gtau
               ENDDO
            ENDIF
         ENDIF
      ENDDO
      CALL cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,+1)
      IF (.not.lsda) THEN
         DO ir=1,nrxx
            aux(ir) = aux(ir) * dmuxc(ir,1,1)
         END DO
      ELSE
         is=isk(ikk)
         DO ir=1,nrxx
            aux(ir) = aux(ir) * 0.5d0 *  &
                 (dmuxc(ir,is,1)+dmuxc(ir,is,2))
         ENDDO
      ENDIF
      CALL cft3(aux,nr1,nr2,nr3,nrx1,nrx2,nrx3,-1)
      IF (doublegrid) THEN
         auxs(:) = (0.d0, 0.d0)
         DO ig=1,ngms
            auxs(nls(ig)) = aux(nl(ig))
         ENDDO
      ENDIF
      aux1(:) = aux1(:) + auxs(:)
   ENDIF
  !
  ! Now we compute dV_loc/dtau in real space
  !
  CALL cft3s (aux1, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 1)
  DO ibnd = 1, nbnd
     aux2(:) = (0.d0, 0.d0)
     DO ig = 1, npw
        aux2 (nls (igk (ig) ) ) = evc (ig, ibnd)
     ENDDO
     !
     !  This wavefunction is computed in real space
     !
     CALL cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, + 2)
     DO ir = 1, nrxxs
        aux2 (ir) = aux2 (ir) * aux1 (ir)
     ENDDO
     !
     ! and finally dV_loc/dtau * psi is transformed in reciprocal space
     !
     CALL cft3s (aux2, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, - 2)
     DO ig = 1, npwq
        dvpsi (ig, ibnd) = aux2 (nls (igkq (ig) ) )
     ENDDO
  ENDDO
  !
  DEALLOCATE (aux2)
  DEALLOCATE (aux1)
  IF (nlcc_any) THEN
     DEALLOCATE (aux)
     IF (doublegrid) DEALLOCATE (auxs)
  ENDIF
  !
  !   We add the contribution of the nonlocal potential in the US form
  !   First a term similar to the KB case.
  !   Then a term due to the change of the D coefficients in the perturbat
  !
  CALL dvqpsi_us_only3 (ik, mode, uact, xxk)

  CALL stop_clock ('dvqpsi_us')

END SUBROUTINE dvqpsi_us3
