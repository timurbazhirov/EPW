  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------------
  subroutine loadumat ( nbnd, nbndsub, nks, nkstot, xxq, cu, cuq )
  !-----------------------------------------------------------------------
  !
  !   wannier interpolation of e-p vertex:
  !   load rotation matrix on coarse mesh and distribute
  !
  !   07/2010 JN changed the way this was done.  Previously
  !   a pool scatter call that didn't work on hbar was performed
  !   and some crazy packing scheme was needed to use poolscatter
  !   subsequently it is just a bcast followed by an appropriate assignment
  !-----------------------------------------------------------------------
#include "f_defs.h"
    USE kinds, only : DP
    USE io_global, ONLY : stdout
    use klist_epw, only : kmap   
    use control_ph, only : lgamma
    use epwcom, only : filukk
    use control_flags, ONLY : iverbosity
#ifdef __PARA
    USE io_global, ONLY : ionode_id
    USE mp_global, ONLY : mpime, inter_pool_comm, my_pool_id
    use mp,        only : mp_sum, mp_barrier
#endif
    implicit none
    !
    !  input variables
    ! 
    integer :: nbnd, nbndsub, nks, nkstot    ! number of bands
    ! number of bands in the optimal subspace
    ! number of kpoints 
    ! total number of kpoints across pools
    real(kind=DP) :: xxq(3)
    ! the qpoint for folding of U
    !
    !  output variables
    !   
    complex(kind=DP) :: cu ( nbnd, nbndsub, nks), cuq ( nbnd, nbndsub, nks)
    !
    integer, parameter :: iunukk  = 71               
    ! unit with rotation matrix U(k) from wannier code
    !
    ! work variables 
    !
    integer :: ik, ibnd, jbnd, nsize, ind, ios, ik_start, ik_stop
    real(kind=DP) :: aux( 2*nbnd*nbndsub, nkstot), auxq( 2*nbnd*nbndsub, nkstot) 
    complex(kind=DP) :: cu_big ( nbnd, nbndsub, nkstot), cuq_big ( nbnd, nbndsub, nkstot)
    complex(kind=DP), parameter:: czero = (0.d0, 0.d0)

    !
    cu_big = czero
    cuq_big = czero
#ifdef __PARA
    IF (mpime.eq.ionode_id) then
#endif
      !
      ! first proc read rotation matrix (coarse mesh) from file
      !
      open ( unit = iunukk, file = filukk, status = 'old', form = 'formatted',iostat=ios)
      IF (ios /=0) call errore ('loadumat', 'error opening ukk file',iunukk)

      !
      DO ik = 1, nkstot
        !
        !
        DO ibnd = 1, nbnd
          DO jbnd = 1, nbndsub
            read(iunukk, *) cu_big (ibnd, jbnd, ik)
          ENDDO
        ENDDO
      ENDDO
      !
      close ( iunukk )
      !
      !  generate U(k+q) through the map 
      !
      !
      ! here we create the map k+q --> k
      ! (first proc has a copy of all kpoints)
      !
      !  generates kmap(ik) for this xxq
      !
      CALL createkmap2 ( xxq )
      !
      !  and we generate the matrix for the q-displaced mesh
      !
      DO ik = 1, nkstot
         cuq_big (:, :, ik) = cu_big (:, :, kmap(ik) )  
      ENDDO
      !
#ifdef __PARA
    ENDIF
    CALL mp_sum (cu_big,  inter_pool_comm)
    CALL mp_sum (cuq_big, inter_pool_comm)
    !
#endif
  !
    CALL ckbounds(ik_start, ik_stop)
    IF ( (ik_stop-ik_start+1) .ne. nks) call errore('loadumat',"Improper parallel ukk load",1)
    cu = cu_big (:, :, ik_start:ik_stop)
    cuq = cuq_big (:, :, ik_start:ik_stop)
    !
  end subroutine loadumat

