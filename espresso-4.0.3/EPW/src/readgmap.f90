  !                                                                            
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  !--------------------------------------------------------------
  subroutine readgmap ( nkstot, ngxx, ng0vec, g0vec_all_r ) 
  !--------------------------------------------------------------
  !
  !  read map of G vectors G -> G-G_0 for a given q point
  !  (this is used for the folding of k+q into the first BZ) 
  !    
  !
  !--------------------------------------------------------------

#ifdef __PARA
  USE mp_global,   ONLY : my_pool_id, inter_pool_comm, mpime
  USE mp,          ONLY : mp_bcast,mp_max,mp_min
#endif
  USE kinds, only : DP
  use io_global, only : stdout, ionode_id
  use wvfct, only : npwx
  use pwcom, only : ngm, nks
  use el_phon, only : shift, gmap 
  USE io_files, ONLY: iunigk, prefix
  implicit none
  !
  ! variables for folding of k+q grid
  !
  INTEGER :: iukgmap, ng0vec, ngxx, iukmap, nkstot
  !   unit with map of folding G-vector indexes
  !   kpoint blocks (k points of all pools)
  !   number of inequivalent such translations
  !   bound for the allocation of the array gmap
  !   unit with map k --> k+q+G
  REAL(kind=DP) :: g0vec_all_r(3,125)
  !   G-vectors needed to fold the k+q grid into the k grid, cartesian coord.
  !
  !  work variables
  !
  integer :: ik, ik1, ishift, ig0, ig, itmp, npw0, igk0 (npwx)
  real(kind=DP) :: tmp
  integer :: ios
  !
  allocate ( shift(nkstot) )
  !
  !  OBSOLETE: now we read directly the igkq to get the proper ngxx
  !
  !  read only a piece of the map to save time 
  !  the proper allocation bound would be ngxx = max(max(igkq))
  !  where the max is taken over the ig and the ik
  !  Here I use a simpler estimate: take the sphere npwx + two
  !  extra shells. This may not work for strange shapes of the
  !  reciproc latt. In this case just set ngxx = ngm_g
  !
  !  ngxx = nint(4./3.*3.14*(2+(3.0/4.0/3.14*float(npwx))**(1./3.))**3.)
  !

  !  Note that the k+q point below does not correspond to the actual (true) 
  !  k+q, but since we only need to take the max over k and k+q this
  !  does not matter
  !
  ngxx = 0
    rewind (unit = iunigk)
    DO ik = 1, nks
      !
      read (iunigk) npw0, igk0 ! k 
      !
      IF (maxval(igk0(1:npw0)).gt.ngxx) ngxx = maxval(igk0(1:npw0))
      !
    ENDDO
  !
#ifdef __PARA
  tmp = float (ngxx)
  CALL mp_max(tmp,inter_pool_comm)  
  ngxx = nint (tmp)
#endif
  !
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
    !
    iukgmap = 96
    open ( unit = iukgmap, file = trim(prefix)//'.kgmap', form = 'formatted',status='old',iostat=ios)
    IF (ios /=0) call errore ('readgmap', 'error opening kgmap file',iukgmap)
    !
    DO ik = 1, nkstot
      read (iukgmap,*) ik1, shift (ik1)
    ENDDO
    read (iukgmap,*) ng0vec
    !
    !  the following seems crazy but I make it for compatibility
    !  with versions up to 2.1.5:
    !
    !  iukgmap has been created by ../PW/set_kplusq.f90 and has
    !  the correct gmap(), but the wrong shift() (actually the
    !  shift for a specific q-point)
    !
    !  since createkmap.f90 has regenerated the shifts for the
    !  present kpoint I read them again in kmap.dat. The above 
    !  'fake' readin is because the gmap appears *after* the
    !  wrong kmap.
    !
    iukmap = 78
    open ( unit = iukmap, file = trim(prefix)//'.kmap', form = 'formatted',status='old',iostat=ios)
    IF (ios /= 0) call errore ('readgmap', 'error opening kmap file', iukmap)
    DO ik = 1, nkstot
      read (iukmap,*) ik1, itmp, shift (ik1)
    ENDDO
    close (iukmap) 
    !
#ifdef __PARA
  ENDIF
  !
  ! first node broadcasts ng0vec to all nodes for allocation of gmap
  !
  CALL mp_bcast( ng0vec, ionode_id, inter_pool_comm )
#endif
  !
  allocate ( gmap (ngxx * ng0vec) )
  !
#ifdef __PARA
  IF (mpime.eq.ionode_id) then
#endif
     !
    DO ig0 = 1, ng0vec
      read (iukgmap,*) g0vec_all_r (:,ig0)
    ENDDO
    DO ig = 1, ngxx
      ! 
      ! at variance with the nscf calculation, here gmap is read as a vector,
      ! 
      read (iukgmap,*) (gmap ( ng0vec * ( ig - 1 ) + ishift ), ishift = 1, ng0vec)
    ENDDO
    !
    close (iukgmap)
    !
#ifdef __PARA
  ENDIF
  !
  ! first node broadcasts everything to all nodes
  !
  CALL mp_bcast( g0vec_all_r, ionode_id, inter_pool_comm )
  CALL mp_bcast( shift, ionode_id, inter_pool_comm )
  CALL mp_bcast( gmap, ionode_id, inter_pool_comm )
  !
#endif
  !
  end subroutine readgmap

