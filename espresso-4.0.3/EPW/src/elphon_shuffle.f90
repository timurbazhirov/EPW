  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                        
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  SUBROUTINE elphon_shuffle ( iq_irr, nqc_irr, iq, gmapsym, eigv, isym, invs0, xq0, timerev )
  !-----------------------------------------------------------------------
  !
  ! Electron-phonon calculation from data saved in fildvscf
  ! Shuffle2 mode (shuffle on electrons + load all phonon q's)
  !
  ! no ultrasoft yet
  !
  !-----------------------------------------------------------------------
  !
#include "f_defs.h"
  !
#ifdef __PARA
  USE mp
  USE mp_global,   ONLY : nproc, my_pool_id, nproc_pool,npool,kunit,&
                          inter_pool_comm
#endif
  !
  USE ions_base, ONLY : nat, ntyp => nsp, ityp, tau
  !
  USE io_global, ONLY : stdout
  USE pwcom,     ONLY : nbnd, ngm, nkstot, nspin, nrxx, nrxxs, doublegrid, nks, xk
  USE kinds,     ONLY : DP
  USE phcom,     ONLY : nmodes, nirr, npert, u
  USE epwcom,    ONLY : tshuffle, tshuffle2
  USE el_phon,   ONLY : epmatq, el_ph_mat
  implicit none
  complex(kind=DP), PARAMETER :: czero = (0.d0,0.d0), cone = (1.d0,0.d0)
  !
  integer :: irr, imode0, ipert, is, iq, iq_irr, nqc_irr
  ! counter on the representations
  ! counter on the modes
  ! the change of Vscf due to perturbations
  ! the current qpoint in the uniform grid
  ! the current ireducible qpoint
  ! the total number of irreducible qpoints in the list
  complex(kind=DP), POINTER :: dvscfin(:,:,:), dvscfins (:,:,:)
  logical :: timerev
  !  true if we are using time reversal
  !
  integer :: tmp_pool_id, nkl, nkr, ik0, iks, ik, ibnd, jbnd
  !
  integer :: gmapsym ( ngm, 48 ), isym, invs0 (48)
  ! the correspondence G-->S(G)
  ! the symmetry which generates the current q in the star
  ! the index of the inverse operations
  complex(kind=DP) :: eigv (ngm, 48)
  ! e^{ iGv} for 1...nsym ( v the fractional translation)
  real(kind=DP) :: xq0(3)
  ! the first q in the star (cartesian)
  !
  CALL start_clock ('elphon')
  !
  ! tshuffle2 implies tshuffle
  !
  IF (tshuffle2) tshuffle = .true.
  !
  ik0 = 0
  tmp_pool_id = 0
  !
#ifdef __PARA
  !
  npool = nproc / nproc_pool
  IF (npool.gt.1) THEN
    !
    ! number of kpoint blocks, kpoints per pool and reminder
    kunit = 1 
    nkl   = kunit * ( nkstot / npool )
    nkr   = ( nkstot - nkl * npool ) / kunit
    ! the reminder goes to the first nkr pools
    IF ( my_pool_id < nkr ) nkl = nkl + kunit
    !
    iks = nkl * my_pool_id + 1
    IF ( my_pool_id >= nkr ) iks = iks + nkr * kunit
    !
    !  the index of the first k point block in this pool - 1
    !  (I will need the index of ik, not ikk)
    !
    ik0 = ( iks - 1 ) / kunit
    !
  ENDIF
  !
#endif
  !    
  ! read Delta Vscf and calculate electron-phonon coefficients
  !
  imode0 = 0
  DO irr = 1, nirr
     ALLOCATE (dvscfin ( nrxx , nspin , npert(irr)) )
     !
     !   read the <prefix>.dvscf_q[iq] files
     !
     dvscfin = (0.d0,0.d0)
#ifdef __PARA
     IF (my_pool_id.eq.0) THEN
#endif
        DO ipert = 1, npert (irr)
           CALL readdvscf ( dvscfin(1,1,ipert), imode0 + ipert, iq_irr, nqc_irr )
        END DO
#ifdef __PARA
     ENDIF
     CALL mp_sum(dvscfin,inter_pool_comm)
     !
#endif
     !
     IF (doublegrid) THEN
        allocate (dvscfins ( nrxxs , nspin , npert(irr)) )
        DO is = 1, nspin
           DO ipert = 1, npert(irr)
              CALL cinterpolate (dvscfin(1,is,ipert),dvscfins(1,is,ipert),-1)
           ENDDO 
        ENDDO
     ELSE
        dvscfins => dvscfin
     ENDIF
     !
     CALL elphel2_shuffle (npert (irr), imode0, dvscfins, gmapsym, eigv, isym, invs0, xq0, timerev)
     !
     imode0 = imode0 + npert (irr)
     IF (doublegrid) DEALLOCATE (dvscfins)
     DEALLOCATE (dvscfin)
  ENDDO
  !
  !  the output e-p matrix in the pattern representation
  !  must be transformed in the cartesian basis
  !  epmat_{CART} = conjg ( U ) * epmat_{PATTERN}
  !
  !  note it is not U^\dagger ! Have a look to symdyn_munu.f90 
  !  for comparison
  !
  DO ibnd = 1, nbnd
    DO jbnd = 1, nbnd
      DO ik = 1, nks
        ! 
        ! Here is where we calculate epmatq, it appears to be
        ! epmatq = cone * conjug(u) * el_ph_mat + czero  
        CALL zgemv ('n', nmodes, nmodes, cone, CONJG ( u ), nmodes, &
          el_ph_mat (ibnd,jbnd,ik,:), 1, czero, epmatq (ibnd,jbnd,ik,:,iq), 1 )
        !
      ENDDO
    ENDDO
  ENDDO
  !
  CALL stop_clock ('elphon')
  !
  END SUBROUTINE elphon_shuffle
