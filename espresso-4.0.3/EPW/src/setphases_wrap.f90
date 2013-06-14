  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  subroutine setphases_wrap
  !---------------------------------------------------------------------
  !
  !  This is the wrapper which is used to set the phases of the wavefunctions  
  !  at k and k+q on the coarse mesh.  It should only be called once.
  !  Note that the phases at k+q are for the input 'q' vector, not
  !  the one in the dynamical coarse list.
  !  
  !  The phases for all wavefunctions are now stored in umat_all
  !  
  !
  !---------------------------------------------------------------------
  USE kinds,           ONLY : DP
  use klist,           only : xk, nkstot
  use io_files,        only : iunigk
  use io_global,       only : ionode
  use mp_global,       only : inter_pool_comm, my_pool_id, nproc_pool
  use mp,              only : mp_sum
  use phcom,           only : npwq, igkq, lgamma, lrwfc,  iuwfc, evq, lrdrho
  use el_phon,         only : umat, umat_all
  use pwcom,           only : npw, igk, nbnd, nks
  use epwcom,           only : tphases, fildvscf0, iudvscf0
  use wavefunctions_module, only : evc
  use io_files,        only : iunwfc, iunigk
  use control_flags,   only : iverbosity
  use pwcom,           only : nrxxs
  use lsda_mod,        only : nspin
  use ions_base,       only : nat
  !
  !
  implicit none
  !
  complex(kind=DP) :: v1(nrxxs,nspin), v2(nrxxs,nspin), v3(nrxxs,nspin)
  ! tmp matrices to build deltav
  real(kind=DP) :: deltav(nrxxs)
  ! the fake (local) perturbation in real space, it is real to guarantee 
  ! hermiticity
  integer :: ik, ios, nkk, nkk_abs, ipool, ibnd, jbnd
  real(kind=DP) :: zero_vect(3)
  !
  IF (nproc_pool>1) call errore &
       ('setphases_wrap', 'only one proc per pool', 1)
  !
  allocate (umat_all (nbnd, nbnd, nkstot))
  allocate (umat(nbnd,nbnd,nks))
  umat_all = (0.d0, 0.d0)
  zero_vect = 0.d0
  !
  IF (tphases) then
     !
     WRITE (6,'(5x,a)') 'Setting the phases on |psi_k>'
     !
     CALL davcio_drho ( v1,  lrdrho, iudvscf0,       1, -1 )
     CALL davcio_drho ( v2,  lrdrho, iudvscf0, 3*nat/2, -1 )
     CALL davcio_drho ( v3,  lrdrho, iudvscf0, 3*nat  , -1 )
     deltav= real ( v1(:,1) + v2(:,1) + v3(:,1))
     deltav=deltav ** 3.d0
     !
     !
     IF (nks.gt.1) rewind (unit = iunigk)
     DO ik=1,nks
        IF (nks.gt.1) then
           read (iunigk, err = 100, iostat = ios) npw, igk
100        call errore ('setphases_wrap', 'reading igk', abs (ios) )
        ENDIF
        !
        !
        IF (nks.gt.1) then
           CALL davcio (evc, lrwfc, iuwfc, ik, - 1)
        ENDIF
        !
        CALL ktokpmq ( xk(:,ik), zero_vect, +1, ipool, nkk, nkk_abs)
        !
        CALL setphases ( 1, ik, npw, umat(:,:,ik),deltav)
        umat_all(:,:,nkk_abs) = umat(:,:,ik)
        !
        !
     END DO
  ELSE ! no phases, rotation matrix is then the identity
     !
     WRITE(6,'(5x,a)') 'No wavefunction gauge setting applied' 
     !
     IF (ionode) then
        DO ik = 1, nkstot
           DO ibnd = 1, nbnd
              DO jbnd = 1, nbnd
                 IF (ibnd .eq. jbnd) then
                    umat_all(ibnd,jbnd,ik) = (1.d0, 0.d0)
                 ELSE
                    umat_all(ibnd,jbnd,ik) = (0.d0,0.d0)
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDIF
     DO ik = 1, nks
        DO ibnd = 1, nbnd
           DO jbnd = 1, nbnd
              IF (ibnd .eq. jbnd) then
                 umat(ibnd,jbnd,ik) = (1.d0, 0.d0)
              ELSE
                 umat(ibnd,jbnd,ik) = (0.d0,0.d0)
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDIF
  !
  ! collect the global phase-setting matrix
  !
  CALL mp_sum(umat_all, inter_pool_comm)
  !
  IF (iverbosity .eq. 1) then
     WRITE (6,* ) "Phase setting matrices:"
     DO ik = 1, nkstot
        DO ibnd = 1, nbnd
           WRITE(6, '(8f8.5)') umat_all(:,ibnd, ik)
        ENDDO
        WRITE(6,*)
     ENDDO
  ENDIF
  !
  end subroutine setphases_wrap
