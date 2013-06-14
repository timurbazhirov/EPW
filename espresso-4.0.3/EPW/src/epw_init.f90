  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PH/phq_init - Quantum-ESPRESSO group                 
  !--------------------------------------------------------------------
  SUBROUTINE epw_init(first_run)
  !----------------------------------------------------------------------------
  !
  !     This initialization is done nqc_irr times from elphon_shuffle_wrap
  !     not all of the following code is necessary.  More adaptation from
  !     phq_init is needed   
  !
#include "f_defs.h"
  USE ions_base,            ONLY : nat, ntyp => nsp, ityp, tau
  USE becmod,               ONLY : calbec
  USE constants,            ONLY : eps8
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : iunigk
  USE pwcom
  USE atom,                 ONLY : msh, rgrid
  USE wavefunctions_module, ONLY : evc
  USE kinds,                ONLY : DP
  USE noncollin_module,     ONLY : noncolin, npol
  USE uspp_param,           ONLY : upf
  USE phcom
  !
  IMPLICIT NONE
  !
  ! ... local variables
  !
  INTEGER :: nt, ik, ipol, ibnd, na, ig
    ! counter on atom types
    ! counter on k points
    ! counter on k+q points
    ! counter on polarizations
    ! counter on bands
    ! index for wavefunctions at k
    ! counter on atoms
    ! counter on G vectors
  REAL(DP) :: arg
    ! the argument of the phase
  COMPLEX(DP), ALLOCATABLE :: aux1(:,:)
    ! used to compute alphap
  logical :: first_run
  !
  !
  !
  CALL start_clock( 'epw_init' )
  !
  ALLOCATE( aux1( npwx*npol, nbnd ) )    
  !
  ! ... initialize structure factor array
  !
  CALL struc_fact( nat, tau, ntyp, ityp, ngm, g, bg, nr1, nr2, nr3, &
                   strf, eigts1, eigts2, eigts3 )
  !                 
  DO na = 1, nat
     !
     ! xq here is the first q of the star
     arg = ( xq(1) * tau(1,na) + &
             xq(2) * tau(2,na) + &
             xq(3) * tau(3,na) ) * tpi
     !        
     eigqts(na) = CMPLX( COS( arg ), - SIN( arg ) )
     !
  END DO
  !
  ! compute rhocore for each atomic-type if needed for nlcc
  !
  IF ( nlcc_any ) CALL set_drhoc( xq )
  !
  ! the fourier components of the local potential for each |G|
  !
  CALL init_vloc()
  !
  ! the fourier components of the local potential at q+G
  !
  vlocq(:,:) = 0.D0
  !
  DO nt = 1, ntyp
     !
     CALL setlocq( xq, rgrid(nt)%mesh, msh(nt), rgrid(nt)%rab, rgrid(nt)%r,&
                   upf(nt)%vloc(1), upf(nt)%zp, tpiba2, ngm, g, omega, &
                   vlocq(1,nt) )
     !
  END DO
  !
  ! the parameters defining the pseudopotential
  !
  ! we compute the denominators of the KB types, or the
  ! parameters which define the non-local pseudopotential and
  ! which are independent of the k point for the US case
  !
  !
! this was in the phq_init routine.  uspp?
  CALL init_us_1()
  !
  REWIND( iunigk )
  !
  DO ik = 1, nks
     !
     !
     IF ( lsda ) current_spin = isk( ik )
     !
     ! g2kin is used here as work space
     !
     CALL gk_sort( xk(1,ik), ngm, g, ( ecutwfc / tpiba2 ), npw, igk, g2kin )
     !
     ! if there is only one k-point evc, evq, npw, igk stay in memory
     !
     WRITE( iunigk ) npw, igk
     !
     npwq = npw
     ! The functions vkb(k+G)
     !
     CALL init_us_2( npw, igk, xk(1,ik), vkb )
     !
     ! ... read the wavefunctions at k
     !
     CALL davcio( evc, lrwfc, iuwfc, ik, -1 )
     !
     ! we compute the becp terms which are used in the rest of
     ! the code
     !
     IF (noncolin) THEN
        CALL calbec (npw, vkb, evc, becp1_nc(:,:,:,ik) )
     ELSE
        CALL calbec (npw, vkb, evc, becp1(:,:,ik) )
     ENDIF
     !
     ! we compute the derivative of the becp term with respect to an
     !   atomic displacement
     !
     DO ipol = 1, 3
        aux1=(0.d0,0.d0)
        DO ibnd = 1, nbnd
           DO ig = 1, npw
              aux1(ig,ibnd) = evc(ig,ibnd) * tpiba * ( 0.D0, 1.D0 ) * & 
                              ( xk(ipol,ik) + g(ipol,igk(ig)) )
           END DO
           IF (noncolin) THEN
              DO ig = 1, npw
                 aux1(ig+npwx,ibnd)=evc(ig+npwx,ibnd)*tpiba*(0.D0,1.D0)*& 
                           ( xk(ipol,ik) + g(ipol,igk(ig)) )
              END DO
           END IF
        END DO
        IF (noncolin) THEN
           CALL calbec (npw, vkb, aux1, alphap_nc(:,:,:,ipol,ik) )
        ELSE
           CALL calbec (npw, vkb, aux1, alphap(:,:,ipol,ik) )
        END IF
     END DO
     !
     !
  END DO
  !
  DEALLOCATE( aux1 )
  !
  IF (.not.first_run) CALL dvanqq2()
  !
  !
  CALL stop_clock( 'epw_init' )
  !
  END SUBROUTINE epw_init
