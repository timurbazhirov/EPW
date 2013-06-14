!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------
SUBROUTINE dvanqq2
  !----------------------------------------------------------------------
  !
  ! New
  ! This routine calculates two integrals of the Q functions and
  ! its derivatives with c V_loc and V_eff which are used
  ! to compute term dV_bare/dtau * psi  in addusdvqpsi.
  ! The result is stored in int1,int2. The routine is called
  ! for each q in nqc. 
  ! 
  !
  ! OLD
  ! This routine calculates four integrals of the Q functions and
  ! its derivatives with c V_loc and V_eff which are used
  ! to compute term dV_bare/dtau * psi  in addusdvqpsi and in addusdynmat.
  ! The result is stored in int1,int2,int4,int5. The routine is called
  ! ONLY once. int4 and int5 are deallocated after USE in
  ! addusdynmat, and int1 and int2 saved on disk by that routine.
  !
#include "f_defs.h"
  !
  USE mp_global,        ONLY : my_pool_id, npool
  USE ions_base,        ONLY : nat, ityp, ntyp => nsp
  USE io_files,         ONLY : prefix, tmp_dir
  USE pwcom
  USE scf,              ONLY : v, vltot
  USE noncollin_module, ONLY : noncolin, npol
  USE kinds,            ONLY : DP
  USE phcom
  USE uspp_param,       ONLY: upf, lmaxq, nh
  USE mp_global,        ONLY: intra_pool_comm
  USE mp,               ONLY: mp_sum

  implicit none
  !
  !   And the local variables
  !

  integer :: nt, na, nb, ig, nta, ntb, ir, ih, jh, ijh, ipol, jpol, is, nspin0
  ! counters
  integer :: is1, is2, ijs, lh, kh, find_ijh

  real(DP), ALLOCATABLE :: qmod (:), qmodg (:), qpg (:,:), &
       ylmkq (:,:), ylmk0 (:,:)
  ! the modulus of q+G
  ! the modulus of G
  ! the  q+G vectors
  ! the spherical harmonics

  complex(DP) :: fact, fact1, ZDOTC
  complex(DP), ALLOCATABLE :: aux1 (:), aux2 (:),&
       aux3 (:), aux5 (:), veff (:,:), sk(:)
  ! work space
  complex(DP), ALLOCATABLE, TARGET :: qgm(:)
  ! the augmentation function at G
  complex(DP), POINTER :: qgmq (:)
  ! the augmentation function at q+G
  character (len=256) :: tempfile
  character (len=3) :: filelab
  integer :: iurecover
  logical :: exst

  IF (.not.okvan) RETURN

!  if (recover.and..not.ldisp) return

  nspin0=nspin
  IF (nspin==4.and..not.domag) nspin0=1

  CALL start_clock ('dvanqq2')
  int1(:,:,:,:,:) = (0.d0, 0.d0)
  int2(:,:,:,:,:) = (0.d0, 0.d0)
  ALLOCATE (sk  (  ngm))    
  ALLOCATE (aux1(  ngm))    
  ALLOCATE (aux2(  ngm))    
  ALLOCATE (aux3(  ngm))    
  ALLOCATE (aux5(  ngm))    
  ALLOCATE (qmodg( ngm))    
  ALLOCATE (ylmk0( ngm , lmaxq * lmaxq))    
  ALLOCATE (qgm  ( ngm))    
  ALLOCATE (ylmkq(ngm , lmaxq * lmaxq))    
  ALLOCATE (qmod( ngm))    
  ALLOCATE (qgmq( ngm))    
  !
  !     compute spherical harmonics
  !
  CALL ylmr2 (lmaxq * lmaxq, ngm, g, gg, ylmk0)
  DO ig = 1, ngm
     qmodg (ig) = sqrt (gg (ig) )
  ENDDO
  ALLOCATE (qpg (3, ngm))    
  CALL setqmod (ngm, xq, g, qmod, qpg)
  CALL ylmr2 (lmaxq * lmaxq, ngm, qpg, qmod, ylmkq)
  DEALLOCATE (qpg)
  DO ig = 1, ngm
     qmod (ig) = sqrt (qmod (ig) )
  ENDDO
  !   we start by computing the FT of the effective potential
  !
  ALLOCATE (veff ( nrxx , nspin))    
  DO is = 1, nspin
     IF (nspin.ne.4.or.is==1) THEN
        DO ir = 1, nrxx
           veff (ir, is) = CMPLX (vltot (ir) + v%of_r (ir, is), 0.d0)
        ENDDO
     ELSE
        DO ir = 1, nrxx
           veff (ir, is) = CMPLX (v%of_r (ir, is), 0.d0)
        ENDDO
     ENDIF
     CALL cft3 (veff (1, is), nr1, nr2, nr3, nrx1, nrx2, nrx3, - 1)
  ENDDO
  !
  !     We compute here two of the three integrals needed in the phonon
  !
  fact1 = CMPLX (0.d0, - tpiba * omega)
  !
  tempfile = trim(tmp_dir) // trim(prefix) // '.recover' 
#ifdef __PARA
     CALL set_ndnmbr (0,my_pool_id+1,1,npool,filelab)
  tempfile = trim(tmp_dir) // trim(prefix) // '.recover' // filelab
#endif
  iurecover = 178
  !
  IF (recover) THEN
     WRITE (6,*) "    Using recover mode for USPP"
     inquire(file = tempfile, exist=exst)
     IF (.not. exst ) CALL errore( 'dvanqq2', 'recover files not found ', 1)
     OPEN (iurecover, file = tempfile, form = 'unformatted')
     READ (iurecover) int1, int2
     CLOSE(iurecover)
  ELSE
  DO ntb = 1, ntyp
     IF (upf(ntb)%tvanp ) THEN
        ijh = 0
        DO ih = 1, nh (ntb)
           DO jh = ih, nh (ntb)
              ijh = ijh + 1
              !
              !    compute the augmentation function
              !
              CALL qvan2 (ngm, ih, jh, ntb, qmodg, qgm, ylmk0)
              CALL qvan2 (ngm, ih, jh, ntb, qmod, qgmq, ylmkq)
              !
              !     NB: for this integral the moving atom and the atom of Q
              !     do not necessarily coincide
              !
              DO nb = 1, nat
                 IF (ityp (nb) == ntb) THEN
                    DO ig = 1, ngm
                       aux1 (ig) = qgmq (ig) * eigts1 (ig1 (ig), nb) &
                                             * eigts2 (ig2 (ig), nb) &
                                             * eigts3 (ig3 (ig), nb)
                    ENDDO
                    DO na = 1, nat
                       fact = eigqts (na) * CONJG(eigqts (nb) )
                       !
                       !    nb is the atom of the augmentation function
                       !
                       nta = ityp (na)
                       DO ig=1, ngm
                          sk(ig)=vlocq(ig,nta) * eigts1(ig1 (ig), na) &
                                               * eigts2(ig2 (ig), na) &
                                               * eigts3(ig3 (ig), na) 
                       ENDDO
                       DO ipol = 1, 3
                          DO ig=1, ngm
                            aux5(ig)= sk(ig) * (g (ipol, ig) + xq (ipol) )
                          ENDDO
                          int2 (ih, jh, ipol, na, nb) = fact * fact1 * &
                                ZDOTC (ngm, aux1, 1, aux5, 1)
                       ENDDO
                    ENDDO
                    DO ig = 1, ngm
                       aux1 (ig) = qgm (ig) * eigts1 (ig1 (ig), nb) &
                                               * eigts2 (ig2 (ig), nb) &
                                               * eigts3 (ig3 (ig), nb)
                    ENDDO
                    DO is = 1, nspin0
                       DO ipol = 1, 3
                          DO ig = 1, ngm
                             aux2 (ig) = veff (nl (ig), is) * g (ipol, ig)
                          ENDDO
                          int1 (ih, jh, ipol, nb, is) = - fact1 * &
                               ZDOTC (ngm, aux1, 1, aux2, 1)
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        DO ih = 1, nh (ntb)
           DO jh = ih + 1, nh (ntb)
              !
              !    We use the symmetry properties of the integral factor
              !
              DO nb = 1, nat
                 IF (ityp (nb) == ntb) THEN
                    DO ipol = 1, 3
                       DO is = 1, nspin
                          int1(jh,ih,ipol,nb,is) = int1(ih,jh,ipol,nb,is)
                       ENDDO
                       DO na = 1, nat
                          int2(jh,ih,ipol,na,nb) = int2(ih,jh,ipol,na,nb)
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
  ENDDO
#ifdef __PARA
  CALL mp_sum(  int1, intra_pool_comm )
  CALL mp_sum(  int2, intra_pool_comm )
#endif
     OPEN (iurecover, file = tempfile, form = 'unformatted')
     WRITE (iurecover) int1, int2
     CLOSE(iurecover)
  ENDIF


  DEALLOCATE (veff)
  DEALLOCATE(qgmq)
  DEALLOCATE (qmod)
  DEALLOCATE (ylmkq)
  DEALLOCATE (qgm)
  DEALLOCATE (ylmk0)
  DEALLOCATE (qmodg)
  DEALLOCATE (aux5)
  DEALLOCATE (aux3)
  DEALLOCATE (aux2)
  DEALLOCATE (aux1)
  DEALLOCATE (sk)

  CALL stop_clock ('dvanqq2')
  CALL print_clock('dvanqq2')
  RETURN
END SUBROUTINE dvanqq2
