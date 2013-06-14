!
! Copyright (C) 2003 A. Smogunov 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Main program for conductance calculation. 
! This program generalizes to US-PP the method of Choi and Ihm
! ( PRB 59, 2267 (1999) ).
! The ballistic conductance G is calculated via the
! Landauer formula ( G=e^2/h T ), where the total
! transmission T is obtained by solving the scattering 
! problem. 

program pwcond  

 character :: nodenumber * 3  
 call start_postproc (nodenumber)  
 call do_cond (nodenumber)  

 call stop_pp  
stop  
end program pwcond
