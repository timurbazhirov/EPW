  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !-----------------------------------------------------------------
  subroutine rwepmatw ( epmatw, nbnd, np, nmodes, nrec, iun, iop)
  !-----------------------------------------------------------------
  !
  ! A simple wrapper to the davcio routine to read/write arrays
  ! instead of vectors 
  !-----------------------------------------------------------------
  USE kinds, only : DP
#ifdef __PARA
  use mp, only : mp_barrier
  use mp_global, only : my_pool_id 
#endif
  implicit none
  integer :: lrec, iun, nrec, iop, i, nbnd, np, nmodes, ibnd, jbnd, imode, ip
  !
  ! np is either nrr_k or nq (epmatwe and epmatwp have the same structure)
  !
  complex(kind=DP):: epmatw(nbnd,nbnd,np,nmodes), &
     aux ( nbnd*nbnd*np*nmodes )
  !
  lrec = 2 * nbnd * nbnd * np * nmodes
  !
  IF ( iop .eq. -1 ) then
    !
    !  read matrix
    !
    CALL davcio ( aux, lrec, iun, nrec, -1 )
    !
    i = 0
    DO ibnd = 1, nbnd
     DO jbnd = 1, nbnd
      DO ip = 1, np
       DO imode = 1, nmodes
         i = i + 1
         epmatw ( ibnd, jbnd, ip, imode ) = aux (i)
         
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
  ELSEif ( iop .eq. 1 ) then
    !
    !  write matrix
    !
    i = 0
    DO ibnd = 1, nbnd
     DO jbnd = 1, nbnd
      DO ip = 1, np
       DO imode = 1, nmodes
         i = i + 1
         aux (i) = epmatw ( ibnd, jbnd, ip, imode ) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
    CALL davcio ( aux, lrec, iun, nrec, +1 )
    !
  ELSE
    !
    CALL errore ('rwepmatw','iop not permitted',1)
    !
  ENDIF
  !
  end subroutine rwepmatw


  !-----------------------------------------------------------------
  subroutine rwepmatwp (nbnd, np, nmodes, nrec, iun, iop)
  !-----------------------------------------------------------------
  !
  ! A simple wrapper to the davcio routine to read/write arrays
  ! instead of vectors 
  !
  !-----------------------------------------------------------------
  USE kinds, only : DP
  use el_phon, only   : epmatwp
#ifdef __PARA
  use mp, only : mp_barrier
  use mp_global, only : my_pool_id 
#endif
  implicit none
  integer :: lrec, iun, nrec, iop, i, nbnd, np, nmodes, ibnd, jbnd, imode, ip
  !
  ! np is either nrr_k or nq (epmatwe and epmatwp have the same structure)
  !
  complex(kind=DP):: aux ( nbnd*nbnd*np*nmodes )
  !
  lrec = 2 * nbnd * nbnd * np * nmodes
  !
  IF ( iop .eq. -1 ) then
    !
    !  read matrix
    !
    CALL davcio ( aux, lrec, iun, nrec, -1 )
    !
    i = 0
    DO ibnd = 1, nbnd
     DO jbnd = 1, nbnd
      DO ip = 1, np
       DO imode = 1, nmodes
         i = i + 1
         epmatwp ( ibnd, jbnd, ip, imode ) = aux (i)
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
  ELSEif ( iop .eq. 1 ) then
    !
    !  write matrix
    !
    i = 0
    DO ibnd = 1, nbnd
     DO jbnd = 1, nbnd
      DO ip = 1, np
       DO imode = 1, nmodes
         i = i + 1
         aux (i) = epmatwp ( ibnd, jbnd, ip, imode ) 
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
    CALL davcio ( aux, lrec, iun, nrec, +1 )
    !
  ELSE
    !
    CALL errore ('rwepmatwp','iop not permitted',1)
    !
  ENDIF
  !
  end subroutine rwepmatwp
