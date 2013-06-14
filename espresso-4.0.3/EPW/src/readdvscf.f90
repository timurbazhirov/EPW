  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------
  subroutine readdvscf ( dvscf, recn, iq, nqc )
  !--------------------------------------------------------
  !
  !  open dvscf files as direct access, read, ad close again
  !
  !-------------------------------------------------------------
#if defined(__ALPHA)
#  define DIRECT_IO_FACTOR 2
#else
#  define DIRECT_IO_FACTOR 8
#endif

#include "f_defs.h"   
  use io_files, only : prefix, tmp_dir
  use units_ph, only : lrdrho
  use kinds,    only : DP
  USE lsda_mod, ONLY : nspin
  USE gvect,    ONLY : nrxx
  use pwcom
  use epwcom,  only: dvscf_dir
#ifdef __PARA
  use mp_global,only : nproc, nproc_pool, me_pool, mpime
  use mp_global,only : npool,my_pool_id
#endif
  !
  implicit none
  integer :: recn, iq, nqc, iudvscf, statb(13), fstat
  !  perturbation number
  !  the current q point
  !  the total number of qpoints in the list
  !  the temporary unit number
  complex(kind=DP) :: dvscf ( nrxx , nspin) 
  !
  integer :: unf_recl,ios
  character (len=256) :: tempfile
  character (len=3) :: filelab
  ! file label 
  iudvscf = 80
  !
  !  the call to set_ndnmbr is just a trick to get quickly
  !  a file label by exploiting an existing subroutine
  !  (if you look at the sub you will find that the original 
  !  purpose was for pools and nodes)
  !
  CALL set_ndnmbr ( 0, iq, 1, nqc, filelab)
  tempfile = trim(dvscf_dir) // trim(prefix) // '.dvscf_q' // filelab
  !
  unf_recl = DIRECT_IO_FACTOR * lrdrho
  !
  !
  !  open the dvscf file, read and close
  !
  open  (iudvscf, file = tempfile, form = 'unformatted', &
          access = 'direct', iostat=ios,recl = unf_recl,status='old')
  IF (ios /= 0) call errore ('readdvscf','error opening' // tempfile, iudvscf)
  !
  ! check that the binary file is long enough
  ! this is tricky to track through error dumps
  ios = fstat ( iudvscf, statb)
  if (recn * unf_recl .gt. statb(8)) call errore('readdvscf', &
       trim(tempfile)//' too short, check ecut', iudvscf)
  !
  read  (iudvscf, rec = recn) dvscf
  close (iudvscf, status = 'keep')
  !
  !
  end subroutine readdvscf
