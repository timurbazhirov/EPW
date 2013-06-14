  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine readmat_shuffle ( q, dyn, w2, iq, nqc)
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, only : DP
  use io_files, only : prefix, tmp_dir
  USE cell_base, ONLY : ibrav, celldm, omega
  USE ions_base, ONLY : amass, tau, nat, ntyp => nsp, ityp
  USE el_phon
  !
  implicit none
  !
  ! Input
  !
  integer :: iq, nqc
  !  the current q point
  !  the total number of qpoints in the list
  !
  ! output
  !
  complex(kind=DP) :: dyn (3 * nat, 3 * nat)
  real(kind=DP) :: w2 (3 * nat)
  !
  ! local (control variables)
  !
  integer :: iudyn
  real(kind=DP) :: q (3), eps
  integer :: ntyp_, nat_, ibrav_, ityp_, ios
  real(kind=DP) :: celldm_ (6), amass_, tau_ (3), q_ (3)
  !
  ! local
  !
  real(kind=DP) :: dynr (2, 3, nat, 3, nat), sumr(2)
  character(len=80) :: line
  character(len=3)  :: atm
  integer :: nt, na, nb, naa, nbb, nu, mu, i, j, ipol, jpol
  !
  character (len=256) :: tempfile
  character (len=3) :: filelab
  ! file label
  !
  eps = 1.d-5
  !
  !  the call to set_ndnmbr is just a trick to get quickly
  !  a file label by exploiting an existing subroutine
  !  (if you look at the sub you will find that the original
  !  purpose was for pools and nodes)
  !
  CALL set_ndnmbr ( 0, iq, 1, nqc, filelab)
  tempfile = trim(tmp_dir) // trim(prefix) // '.dyn_q' // filelab
  !
  iudyn = 81
  open (unit = iudyn, file = tempfile, status = 'unknown', err = 100, iostat = ios)
100  call errore ('openfilepw', 'opening file'//tempfile, abs (ios) )
  rewind (iudyn)
  !
  !
  !
  read (iudyn, '(a)') line
  read (iudyn, '(a)') line
  read (iudyn, * ) ntyp_, nat_, ibrav_, celldm_
  IF (ntyp.ne.ntyp_.or.nat.ne.nat_.or.ibrav_.ne.ibrav.or.abs ( &
     celldm_ (1) - celldm (1) ) .gt.1.0d-5) call errore ('readmat_shuffle', &
     'inconsistent data', 1)
  DO nt = 1, ntyp
     read (iudyn, * ) i, atm, amass_
     IF (nt.ne.i.or.abs (amass_ - amass (nt) ) .gt.1.0d-5) call errore ( &
          'readmat', 'inconsistent data', 1 + nt)
  ENDDO
  DO na = 1, nat
     read (iudyn, * ) i, ityp_, tau_
     IF (na.ne.i.or.ityp_.ne.ityp (na) ) call errore ('readmat', &
        'inconsistent data', 10 + na)
  ENDDO
  !
  !
  read (iudyn, '(///a)') line
  read (line (11:80), * ) (q_ (i), i = 1, 3)
  IF ( abs(q_ (1)-q(1)).gt.eps .or. &
       abs(q_ (2)-q(2)).gt.eps .or. &
       abs(q_ (3)-q(3)).gt.eps ) call errore &
     ('readmat', 'wrong qpoint', 1)
  read (iudyn, '(a)') line
  !
  !
  DO na = 1, nat
     DO nb = 1, nat
        read (iudyn, * ) naa, nbb
        IF (na.ne.naa.or.nb.ne.nbb) call errore &
             ('readmat', 'error reading file', nb)
        read (iudyn, * ) &
            ((dynr (1,i,na,j,nb), dynr (2,i,na,j,nb), j = 1, 3), i = 1, 3)
     ENDDO
  ENDDO
  !
  IF ( abs(q(1)).lt.eps .and. abs(q(2)).lt.eps .and. abs(q(3)).lt.eps ) then
    !
    ! in the case q=0 we impose the acoustic sum rule
    ! [Gonze and Lee, PRB 55, 10361 (1998), Eq. (45)]
    !
    WRITE(6,'(8x,a)') 'Imposing acoustic sum rule on the dynamical matrix'
    !
    DO na = 1,nat
     DO ipol = 1,3
      DO jpol = ipol,3
         !
         sumr(1) = sum ( dynr (1,ipol,na,jpol,:) )
         sumr(2) = sum ( dynr (2,ipol,na,jpol,:) )
         !
         dynr (:,ipol,na,jpol,na) = dynr (:,ipol,na,jpol,na) - sumr
         !
      END DO
     END DO
    END DO
    !
  ENDIF
  !
  !  fill the two-indices dynamical matrix in cartesian coordinates
  !
  DO na = 1,nat
   DO nb = 1,nat
      DO ipol = 1,3
       DO jpol = 1,3
         !
         mu = (na-1)*3+ipol
         nu = (nb-1)*3+jpol
         dynq ( mu, nu, iq) = dcmplx ( dynr (1,ipol,na,jpol,nb), dynr (2,ipol,na,jpol,nb) )
         !
       END DO
      END DO
   END DO
  END DO
  !
  close(iudyn)
  !
  end subroutine readmat_shuffle

  !-----------------------------------------------------------------------
  subroutine readmat_shuffle2 ( iq_irr, nqc_irr, nq, iq_first, sxq, imq)
  !-----------------------------------------------------------------------
  !
  ! read dynamical matrix for the q points  
  ! iq_first, iq_first+1, ... iq_first+nq-1
  !
  !-----------------------------------------------------------------------
#include "f_defs.h"
  USE kinds, only : DP
  use io_files, only : prefix, tmp_dir
  USE cell_base, ONLY : ibrav, celldm, omega
  USE ions_base, ONLY : amass, tau, nat, ntyp => nsp, ityp
  USE el_phon, ONLY : dynq
  use epwcom, only : dvscf_dir
  !
  implicit none
  !
  ! Input
  !
  integer :: iq_irr, nqc_irr, nq, iq_first, imq, current_iq
  !  the index of the irreducible q point 
  !  the number of irreducible qpoints 
  !  the number of q points in the star of q
  !  the index of the first qpoint to be read in the uniform q-grid
  !  flag which tells whether we have to consider the -q vectors
  !  the q index in the dynq matrix
  real(kind=DP) :: sxq(3,48)
  !  the q vectors in the star
  !
  ! output
  !
  ! dynq (nmode,nmodes,nqc) (in el_phon.mod)
  !
  ! local 
  !
  real(kind=DP) :: eps, celldm_ (6), amass_, tau_ (3), q(3), &
       dynr (2, 3, nat, 3, nat), sumr(2)
  integer :: iudyn, ntyp_, nat_, ibrav_, ityp_, ios, iq, &
       nt, na, nb, naa, nbb, nu, mu, i, j, ipol, jpol
  character(len=80) :: line
  character(len=3)  :: atm, filelab
  character (len=256) :: tempfile
  !
  eps = 1.d-5
  !
  !  the call to set_ndnmbr is just a trick to get quickly
  !  a file label by exploiting an existing subroutine
  !  (if you look at the sub you will find that the original
  !  purpose was for pools and nodes)
  !
  CALL set_ndnmbr ( 0, iq_irr, 1, nqc_irr, filelab)
  tempfile = trim(dvscf_dir) // trim(prefix) // '.dyn_q' // filelab
  !
  iudyn = 81
  open (unit = iudyn, file = tempfile, status = 'old', iostat = ios)
  IF (ios /=0)  call errore ('readmat_shuffle2', 'opening file'//tempfile, abs (ios) )
  !
  !  read header and run some checks
  !
  read (iudyn, '(a)') line
  read (iudyn, '(a)') line
  read (iudyn, * ) ntyp_, nat_, ibrav_, celldm_
  IF (ntyp.ne.ntyp_.or.nat.ne.nat_.or.ibrav_.ne.ibrav.or.abs ( &
     celldm_ (1) - celldm (1) ) .gt.1.0d-5) call errore ('readmat2', &
     'inconsistent data', 1)
  ! 
  !  skip reading of cell parameters here
  ! 
  IF (ibrav_ .eq. 0) then
     DO i = 1,4
        read (iudyn, * ) line
     ENDDO
  ENDIF
  DO nt = 1, ntyp
     read (iudyn, * ) i, atm, amass_
     IF (nt.ne.i.or.abs (amass_ - amass (nt) ) .gt.1.0d-2) then
     write (6,*) amass_, amass(nt)
     call errore ('readmat', 'inconsistent data', 0)
  endif
  ENDDO
  DO na = 1, nat
     read (iudyn, * ) i, ityp_, tau_
     IF (na.ne.i.or.ityp_.ne.ityp (na) ) call errore ('readmat2', &
        'inconsistent data (names)', 10 + na)
  ENDDO
  !
  !  read dyn mat for all q in the star
  !
  current_iq = iq_first-1 
  !
 DO iq = 1, nq
    !
    current_iq = current_iq + 1
    !
    read (iudyn, '(///a)') line
    read (line (11:80), * ) (q (i), i = 1, 3)
    !
    IF ( abs(q(1)-sxq(1,iq)).gt.eps .or. &
         abs(q(2)-sxq(2,iq)).gt.eps .or. &
         abs(q(3)-sxq(3,iq)).gt.eps ) call errore &
       ('readmat', 'wrong qpoint', 1)
    read (iudyn, '(a)') line
    !
    !
     DO na = 1, nat
       DO nb = 1, nat
          read (iudyn, * ) naa, nbb
          IF (na.ne.naa.or.nb.ne.nbb) call errore &
               ('readmat', 'error reading file', nb)
          read (iudyn, * ) &
              ((dynr (1,i,na,j,nb), dynr (2,i,na,j,nb), j = 1, 3), i = 1, 3)
       ENDDO
    ENDDO
    !
    IF ( abs(q(1)).lt.eps .and. abs(q(2)).lt.eps .and. abs(q(3)).lt.eps ) then
      !
      ! in the case q=0 we impose the acoustic sum rule
      ! [Gonze and Lee, PRB 55, 10361 (1998), Eq. (45)]
      !
      WRITE(6,'(8x,a)') 'Imposing acoustic sum rule on the dynamical matrix'
      !
      DO na = 1,nat
       DO ipol = 1,3
        DO jpol = ipol,3
          !
          sumr(1) = sum ( dynr (1,ipol,na,jpol,:) )
          sumr(2) = sum ( dynr (2,ipol,na,jpol,:) )
          !
          dynr (:,ipol,na,jpol,na) = dynr (:,ipol,na,jpol,na) - sumr
          !
        END DO
       END DO
      END DO
      !
    ENDIF
    !
    !  fill the two-indices dynamical matrix in cartesian coordinates
    !  the proper index in the complete list is iq_first+iq-1
    !
    DO na = 1,nat
     DO nb = 1,nat
        DO ipol = 1,3
         DO jpol = 1,3
           !
           mu = (na-1)*3+ipol
           nu = (nb-1)*3+jpol
           dynq ( mu, nu, current_iq) = &
             dcmplx ( dynr (1,ipol,na,jpol,nb), dynr (2,ipol,na,jpol,nb) )
           !
         END DO
        END DO
     END DO
   END DO
   !
   IF (imq.eq.0) then
      !
!  To be added back in when imq.eq.0 is tested: jn
!      current_iq = current_iq + 1
      !
      read (iudyn, '(///a)') line
      !
      ! if imq=0 we have to compare with -sxq
      !
!  To be added back in when imq.eq.0 is tested: jn
!      if ( abs(q(1)+sxq(1,iq)).gt.eps .or. &
!           abs(q(2)+sxq(2,iq)).gt.eps .or. &
!           abs(q(3)+sxq(3,iq)).gt.eps ) call errore &
!           ('readmat', 'wrong qpoint', -1)
      read (iudyn, '(a)') line
      !
      DO na = 1, nat
         DO nb = 1, nat
            read (iudyn, * ) naa, nbb
            IF (na.ne.naa.or.nb.ne.nbb) call errore &
                 ('readmat', 'error reading file', nb)
            read (iudyn, * ) &
                ((dynr (1,i,na,j,nb), dynr (2,i,na,j,nb), j = 1, 3), i = 1, 3)
         ENDDO
      ENDDO
      !
      !  fill the two-indices dynamical matrix in cartesian coordinates
      !  the proper index in the complete list is iq_first+iq-1
      !
      DO na = 1,nat
       DO nb = 1,nat
          DO ipol = 1,3
           DO jpol = 1,3
             !
             mu = (na-1)*3+ipol
             nu = (nb-1)*3+jpol
!  To be added back in when imq.eq.0 is tested: jn
!             dynq ( mu, nu, current_iq) = &
!               dcmplx ( dynr (1,ipol,na,jpol,nb), dynr (2,ipol,na,jpol,nb) )
             !
           END DO
          END DO
       END DO
      END DO
      !
    ENDIF
    !
  ENDDO
  !
  close(iudyn)
  !
  end subroutine readmat_shuffle2

