  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------
  subroutine readwfc ( ipool, recn, evc0 )
  !--------------------------------------------------------
  !
  !  open wfc files as direct access, read, and close again
  !
  !-------------------------------------------------------------
#include "f_defs.h"
  use io_files, only : prefix, tmp_dir
  use units_ph, only : lrwfc, iuwfc
  use kinds,    only : DP
  use wvfct,    only : npwx
  use pwcom
#ifdef __PARA
  use mp_global,only : nproc, nproc_pool, me_pool, mpime
  use mp_global,only : npool,my_pool_id 
#endif
  !
  implicit none
  integer :: recn, ipool
  !  kpoint number
  !  poolfile to be read (not used in serial case)
  complex(kind=DP) :: evc0 ( npwx, nbnd )
  character (len=3) :: nd_nmbr0
  ! node number for shuffle
  !
  integer :: unf_recl, ios
  character (len=256) :: tempfile

  !
  !  open the wfc file, read and close
  !
#ifdef __PARA
  CALL set_ndnmbr ( ipool, me_pool, nproc_pool, npool, nd_nmbr0)

#endif
#if defined (__ALPHA)
#  define DIRECT_IO_FACTOR 2    
# else
#  define DIRECT_IO_FACTOR 8
#endif

#ifdef __PARA
  tempfile = trim(tmp_dir) // trim(prefix) // '.wfc' // nd_nmbr0
# else
  tempfile = trim(tmp_dir) // trim(prefix) // '.wfc'
#endif

  unf_recl = DIRECT_IO_FACTOR * lrwfc
  !
  open  (iuwfc, file = tempfile, form = 'unformatted', &
          access = 'direct', iostat=ios, recl = unf_recl)
  IF (ios /= 0) call errore ('readwfc', 'error opening wfc file', iuwfc)
  read  (iuwfc, rec = recn) evc0
  close (iuwfc, status = 'keep')
  !
  !
  end subroutine readwfc

