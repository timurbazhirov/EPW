!
! Copyright (C) 2003 A. Smogunov
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine summary_tran()
!
! It writes transmission coefficients onto the file tran_file
!
  USE kinds, only : DP
  USE io_files, ONLY : tran_file
  USE cond, ONLY : nenergy, earr, tran_tot
  implicit none

  integer ::  i

!
! Output of T onto the file
!
  open (4,file=trim(tran_file),form='formatted', status='unknown')
  write(4,'("# E-Ef, T")')
  do i=1, nenergy
    write(4,'(2f10.6)') earr(i), tran_tot(i)
  enddo
  close(unit=4)

  return
end subroutine summary_tran


