  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !
  !-----------------------------------------------------------------------
  subroutine refold ( ngm_g, mill_g)
  !----------------------------------------------------------------------
  !
  !   Map the indices of G+G_0 into those of G 
  !   this is used to calculate electron-phonon matrix elements by
  !   refolding the k+q points into the first BZ (original k grid)
  !
  !   No parallelization on G-vecs at the moment  
  !   (actually this is done on the global array, but in elphel2.f90
  !   every processor has just a chunk of the array, I may need some
  !   communication)
  !
  !   No ultrasoft now
  !
  !   I use the rule : if not found then gmap = 0 
  !   Note that the map will be used only up to npwx (small sphere), 
  !   while the G-vectors lost in the process are on the surface of 
  !   the large sphere (density set).
  !
  !
  !-----------------------------------------------------------------
#include "f_defs.h"
  use kinds, only: DP
  use io_global, only : stdout
  USE control_flags, ONLY : iverbosity
  USE kfold
#ifdef __PARA
  use mp_global, only : my_pool_id,me_pool 
#endif
  implicit none
  integer :: mill_g( 3, ngm_g ), ngm_g, ig0, ig1, ig2, i, j, k, &
             notfound, iukgmap, indold, indnew
  logical :: tfound

  !
  IF (iverbosity.eq.1) then 
    WRITE(stdout,*) '  There are ',ng0vec,'inequivalent folding G_0 vectors'
    DO ig0 = 1, ng0vec
      WRITE(stdout,'(a,i1,a,3i3)') 'g0vec_all( ',ig0,') = ',g0vec_all(:,ig0)
    ENDDO
  ENDIF

  !
  allocate ( gmap ( ngm_g, ng0vec ) )
  !
  !  loop on the inequivalent G-vectors
  !
  DO ig0 = 1, ng0vec
     !
     IF (ig0.eq.1) then
        WRITE(6,'(/5x,"Progress kgmap: ")',advance='no')
        indold = 0
     ENDIF
     indnew = nint(float(ig0)/float(ng0vec)*40)
     IF (indnew.ne.indold) write(6,'(a)',advance='no') '#'
     indold = indnew
     CALL flush(6) ! only on opteron
     !
     !
    IF (iverbosity.eq.1) &
      WRITE(stdout,'(i3,4x,3i3)') ig0, g0vec_all(:, ig0 )
    notfound = 0
    DO ig1 = 1, ngm_g
      !
      !  the initial G-vector
      !
      i = mill_g(1, ig1)
      j = mill_g(2, ig1)
      k = mill_g(3, ig1)
      IF (iverbosity.eq.1) &
        WRITE(stdout,'(5x,i5,4x,3i5)') ig1, i,j,k 
      !
      !  the final G-vector
      !
      i = i + g0vec_all(1, ig0 )
      j = j + g0vec_all(2, ig0 )
      k = k + g0vec_all(3, ig0 )
      IF (iverbosity.eq.1) &
        WRITE(stdout,'(5x,i5,4x,3i5)') ig1, i,j,k 
      !
      ig2 = 0
      tfound = .false.
      DO while ((.not.tfound).and.(ig2.lt.ngm_g))
        !
        ig2 = ig2 + 1
        tfound = (i.eq.mill_g(1, ig2)).and. & 
                 (j.eq.mill_g(2, ig2)).and. & 
                 (k.eq.mill_g(3, ig2))
        IF (iverbosity.eq.1) &
          WRITE(stdout,'(10x,i5,4x,3i5,2x,i1)') ig2, mill_g(:, ig2), tfound
        !
      ENDDO
      IF (tfound) then
        gmap ( ig1, ig0 ) = ig2
      ELSE
        gmap ( ig1, ig0 ) = 0
        notfound = notfound + 1
      ENDIF
      !
    ENDDO
    !
    IF (iverbosity.eq.1) &
      WRITE(stdout,*) 'ig0 = ',ig0,' not found = ', notfound, ' out of ',ngm_g
    ! 
  ENDDO
  ! 
  !  output on file for electron-phonon matrix elements
  !
  ! note that this unit must be the same used in set_kplusq.f90 
  iukgmap = 96
  !
#ifdef __PARA
  IF (me_pool.ne.0.or.my_pool_id.ne.0) iukgmap = stdout
#endif
  !
  DO ig1 = 1, ngm_g
    WRITE (iukgmap, '(9i10)') (gmap ( ig1, ig0), ig0 = 1, ng0vec)
  ENDDO
  !
  IF (iukgmap.ne.stdout) close (iukgmap)
  WRITE(stdout,*)
  !
  end subroutine refold

