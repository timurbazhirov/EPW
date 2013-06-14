! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Tone: File adapted from pwi2xsf.f file of XCRYSDEN distribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c     ------------------------------------------------------------------------
      program pwi2xsf
c     Reads pre-procesed (with pwi2xsf.sh) PW-input file
c     and converts to XSF format file
c
c     This program reads the NEWLY formated preprocessed-PW.X input
c     
c     Usage: pwi2xsf.sh < PW-preprocessed file
c     ------------------------------------------------------------------------

      implicit none

c maxtyp : maximum number of types of atoms
c maxatom: maximum number of atoms
c maximage: maximum number of images

      integer
     $     maxtyp,      
     $     maxatom,     
     $     maximage,    
     $     ALAT_UNIT,
     $     BOHR_UNIT,
     $     ANGSTROM_UNIT,
     $     CRYSTAL_UNIT

      real*8
     $     bohr

      parameter (
     $     maxtyp  = 100,
     $     maxatom = 10000,
     $     maximage = 50,
     $     bohr    = 0.5291772108d0,
     $
     $     ALAT_UNIT     = 1,
     $     BOHR_UNIT     = 2,
     $     ANGSTROM_UNIT = 3,
     $     CRYSTAL_UNIT  = 4 )          

      integer
     $     ibrav,               ! label for Bravais lattice
     $     nat,                 ! number of atoms
     $     ntyp,                ! number of pseudopotentials
     $     num_of_images,       ! number of NEB images
     $     inp_num_of_images,   ! number of NEB images in the input
     $     atomic_posunit       ! length-unit of atomic positions

      real*8
     $     celldm(6),           ! cell parameters
     $     omega,               ! cell volume (not used)
     $     alat,                ! lattice parameter
     $     a, b, c, cosab, cosac, cosbc ! lattice parameters

      character
     $     calculation*80,      ! type of calculation
     $     line*120             ! line of input
      character*3
     $     atm(maxatom,maximage) ! atomic symbols

      integer
     $     ityp,                ! type of PP
     $     ounit,               ! output unit     
     $     i, j, ipol,          ! dummies
     $     inat, iim, iim_old,  ! counters
     $     i_trimleft_white_space, ! string whitespace-triming function
     $     len                  ! length of string

      real*8
     $     x,y,z,               ! Cartesian coordinates & weights
     $     w1,w2,               ! linear interpolation weights
     $     dx, dy, dz,          ! auxiliary
     $     tau(3,maxatom,maximage), ! atomic coordinates
     +     pv( 3,3 ),           ! lattice vectors (PRIMITIVE)
     +     cv( 3,3 ),           ! lattice vectors (CONVENTIONAL)
     $     old_total_dist, old_dist(maximage), ! old(=input) inter-image distances
     $     new_total_dist, new_dist ! new(=output) inter-image distances

      logical
     $     ltaucry, matches, last_image

      namelist/system/
     $     ibrav, nat, celldm, a, b, c, cosab, cosac, cosbc,
     $     calculation, num_of_images
      
      ounit=6

c     set default values
      calculation   = 'scf'
      num_of_images = 1
      nat       = 0
      ibrav     = 0
      celldm(1) = 0.0d0
      a         = 0.0D0 
      b         = 0.0D0
      c         = 0.0D0
      cosab     = 0.0D0
      cosac     = 0.0D0
      cosbc     = 0.0D0
c
c     read namelist system
c
      read (5,system)
      if ( nat.eq.0 ) then
         print *,'ERROR: while reading INPUT !!!'
         STOP
      endif

c     was lattice specified in terms of A,B,C,...
      if ( celldm(1) .eq. 0.0D0 .AND. a .ne. 0.0D0 ) THEN
         if ( ibrav .eq. 0 ) ibrav = 14 
         celldm(1) = a / bohr
         celldm(2) = b / a
         celldm(3) = c / a
         celldm(4) = cosab
         celldm(5) = cosac
         celldm(6) = cosbc 
      else if ( celldm(1) .ne. 0.0D0 .AND. a .ne. 0.0D0 ) THEN
         print *, 'ERROR: do not specify both celldm and a,b,c !!!'
      endif
c
c     read the rest of the input
c
 990  continue
      read(5,'(a120)',end=999) line
      len = i_trimleft_white_space(line)

c     
c     CELL_PARAMETERS
c     
      if (     line(1:15) .eq. 'CELL_PARAMETERS' ) then
         read (5,*) ((pv(i,j),i=1,3),j=1,3)
         do j=1,3
            do i=1,3
               cv(i,j) = pv(i,j)
            end do
         end do
c
c     ATOMIC_POSITIONS
c
      elseif ( line(1:16) .eq. 'ATOMIC_POSITIONS' ) then
c     find out the length-unit
         line = line(17:len)
         len  = i_trimleft_white_space(line)
         atomic_posunit = ALAT_UNIT         
         if (len.gt.0 ) then
            if ( matches('ALAT',line) ) then
               atomic_posunit = ALAT_UNIT
            elseif ( matches('BOHR',line) ) then
               atomic_posunit = BOHR_UNIT
            elseif ( matches('CRYSTAL',line) ) then
               atomic_posunit = CRYSTAL_UNIT
            elseif ( matches('ANGSTROM',line) ) then
               atomic_posunit = ANGSTROM_UNIT
            endif
         endif

c     
c     read atoms
c     
         if ( (calculation(1:3) .ne. 'NEB') .and.
     $        (calculation(1:3) .ne. 'SMD') )  then
            iim = 1
            call read_atoms(nat,atm(1,1),tau(1,1,1))
         else
c     
c     path calculation (NEB or SMD): read atoms
c     
            iim = 0
            last_image = .false.
            do while(.not.last_image)               
               iim = iim + 1
               read (5,'(a120)') line ! line: first_image
               if ( matches('LAST_IMAGE',line) ) last_image = .true.
               call read_atoms(nat,atm(1,iim),tau(1,1,iim))
            enddo
         endif
      endif
      inp_num_of_images = iim
      goto 990
      
 999  continue
      
      if ( celldm(1).eq.0.0d0 ) then
         print *,'ERROR while reading INPUT: celldm(1)==0.0d0 !!!'
         STOP
      endif

      
      if ( ibrav.ne.0 ) then
         call latgen( ibrav, celldm,
     $        pv(1,1), pv(1,2), pv(1,3),
     $        cv(1,1), cv(1,2), cv(1,3), omega )
         do j=1,3
            do i=1,3
               pv(i,j) = pv(i,j)/celldm(1)
            end do
         end do
         call latgen_conventional(ibrav, celldm,
     $        pv(1,1), pv(1,2), pv(1,3),
     $        cv(1,1), cv(1,2), cv(1,3))
      endif

      alat = bohr*celldm(1)
      call write_XSF_header (num_of_images,alat, pv, cv, nat, ounit)
      
c     coordinates to ANGSTROMs

      do iim=1,inp_num_of_images
         do inat=1,nat
            if     ( atomic_posunit .eq. BOHR_UNIT ) then            
               tau(1,inat,iim) = bohr * tau(1,inat,iim)
               tau(2,inat,iim) = bohr * tau(2,inat,iim)
               tau(3,inat,iim) = bohr * tau(3,inat,iim)
               
            elseif ( atomic_posunit .eq. ALAT_UNIT ) then
               tau(1,inat,iim) = alat * tau(1,inat,iim)
               tau(2,inat,iim) = alat * tau(2,inat,iim)
               tau(3,inat,iim) = alat * tau(3,inat,iim)
               
            elseif ( atomic_posunit .eq. CRYSTAL_UNIT ) then
               call cryst_to_cart(1, tau(1,inat,iim), pv, 1)
            endif
         enddo
      enddo

      IF ( num_of_images .lt. 2 ) then
c     write atoms for non-PATH calculation
         do inat=1,nat
            write(ounit,'(a3,2x,3f15.10)') atm(inat,1),
     $           tau(1,inat,1), tau(2,inat,1), tau(3,inat,1)
         enddo

      ELSE

c     calculate intermediate images for PATH calculation

         old_total_dist = 0.0d0
         old_dist(1)    = 0.0d0
         do iim = 2, inp_num_of_images
            old_dist(iim) = 0.0
            do inat=1,nat
               dx = tau(1,inat,iim) - tau(1,inat,iim-1)
               dy = tau(2,inat,iim) - tau(2,inat,iim-1)
               dz = tau(3,inat,iim) - tau(3,inat,iim-1)
               old_dist(iim) = old_dist(iim) + dx*dx + dy*dy + dz*dz
            enddo
            old_dist(iim) = sqrt( old_dist(iim) )
            
            old_total_dist = old_total_dist + old_dist(iim)
            old_dist(iim)   = old_total_dist
         enddo
         
         new_dist = old_total_dist / dble(num_of_images-1)
     
c     --------------------------------------------------
c     perform INTERPOLATION
c     --------------------------------------------------
         
         new_total_dist = 0.0
         do iim=1,num_of_images-1
            do iim_old=1,inp_num_of_images-1            
               if ( new_total_dist .ge. old_dist(iim_old)
     $              .and.
     $              new_total_dist .lt. old_dist(iim_old+1) + 1d-10 )
     $              then
               
                  w1 = ( old_dist(iim_old+1) - new_total_dist )
     $                 /
     $                 ( old_dist(iim_old+1) - old_dist(iim_old) )
                  w2 = 1.0d0 - w1
                  
                  write(ounit,'('' PRIMCOORD '',i5)') iim
                  write(ounit,*) nat, 1

                  do inat=1,nat
                     x = w1*tau(1,inat,iim_old)+w2*tau(1,inat,iim_old+1)
                     y = w1*tau(2,inat,iim_old)+w2*tau(2,inat,iim_old+1)
                     z = w1*tau(3,inat,iim_old)+w2*tau(3,inat,iim_old+1)
                     write(ounit,'(a3,2x,3f15.10)')
     $                    atm(inat,iim_old), x, y, z
                  enddo
                  goto 11
               endif
            enddo
 11         continue
            new_total_dist = new_total_dist + new_dist
         enddo
      
c     print last image
         write(ounit,'('' PRIMCOORD '',i5)') iim
         write(ounit,*) nat, 1
         do inat=1,nat
            x = tau(1,inat,inp_num_of_images)
            y = tau(2,inat,inp_num_of_images)
            z = tau(3,inat,inp_num_of_images)
            write(ounit,'(a3,2x,3f15.10)')
     $           atm(inat,inp_num_of_images), x, y, z
         enddo
      endif
      END



c---------------------------------------------------------------------
      subroutine latgen_conventional
     +     ( ibrav, celldm, p1, p2, p3, c1, c2, c3 )
c     Generate convetional lattice
c---------------------------------------------------------------------
c
c   Conventional crystallographic vectors c1, c2, and c3.
c   See "latgen" for the meaning of variables
c
      implicit none
c
c     First the input variables
c
      real*8 
     +     celldm( 6 ),         ! input : the dimensions of the lattice
     +     p1( 3 ),             ! input : first lattice vector (PRIMITIVE)
     +     p2( 3 ),             ! input : second lattice vector
     +     p3( 3 ),             ! input : third lattice vector
     +     c1( 3 ),             ! output: first lattice vector(CONVENTIONAL)
     +     c2( 3 ),             ! output: second lattice vector
     +     c3( 3 )              ! output: third lattice vector
      integer
     +       ibrav          ! input: the index of the Bravais lattice
c
      integer i
c
c
      do i = 1, 3
         c1(i) =0.d0
         c2(i) =0.d0
         c3(i) =0.d0
      end do
c
      if ( ibrav .eq. 2 .or. ibrav .eq.3 ) then
c
c     fcc and bcc lattice
c
         c1( 1 ) = 1.0d0
         c2( 2 ) = 1.0d0
         c3( 3 ) = 1.0d0
c
      else if ( ibrav .eq. 7 ) then
c
c     body centered tetragonal lattice
c
         if ( celldm( 1 ) .le. 0.d0 .or. celldm( 3 ) .le. 0.d0 ) 
     +      call errore( 'latgen', 'wrong celldm', 7 )
         c1(1) = 1.0d0
         c2(2) = 1.0d0
         c3(3) = celldm(3)
c
      else if ( ibrav .eq. 10 ) then
c
c     All face centered orthorombic lattice
c
         if ( celldm( 1 ) .le. 0.d0 .or. celldm( 2 ) .le. 0.d0
     +        .or. celldm( 3 ) .le. 0.d0 ) 
     +      call errore( 'latgen', 'wrong celldm', 10 )
         c1(1) = 1.0d0
         c2(2) = celldm(2)
         c3(3) = celldm(3)
c
      elseif ( ibrav .eq. 11 ) then
c
c     Body centered orthorombic lattice
c
         if ( celldm( 1 ) .le. 0.d0 .or. celldm( 2 ) .le. 0.d0
     +        .or. celldm( 3 ) .le. 0.d0 ) 
     +      call errore( 'latgen', 'wrong celldm', 11 )
         c1(1) = 1.0d0
         c2(2) = celldm(2)
         c3(3) = celldm(3)
      else
c     **********
c     all other cases : just copy p vectors to c vectors !!!
c     **********
         do  i = 1, 3
            c1( i ) = p1( i )
            c2( i ) = p2( i )
            c3( i ) = p3( i )
         enddo
      end if
c
      return
      end


c     ------------------------------------------------------------------------
      subroutine read_atoms(nat,atm,coor)
c     read atomic coordinates
c     ------------------------------------------------------------------------
      implicit none
      integer
     $     nat,                 ! number of atoms
     $     ipol,inat,len,       ! counters
     $     i_trimleft_white_space ! integer-function
      character
     $     line*120             ! line of input
      character*3
     $     atm(*)               ! atomic symbols
      real*8
     $     coor(3,*)
      
      do inat=1,nat
 10      continue
         read (5,'(a120)') line
         len = i_trimleft_white_space(line)
         
         if (len.eq.0) then
c     an empty line, read again
            goto 10
         endif

         read (line,*) atm(inat),(coor(ipol,inat),ipol=1,3)
      enddo
      return
      end


c     ------------------------------------------------------------------------
      subroutine write_XSF_header (num_of_images,alat, p, c, nat, ounit)
c     writes the header for XSF structure file
c     ------------------------------------------------------------------------
      real*8
     $     alat,                ! lattice parameter
     $     p(3,3), c(3,3),      ! lattive vectors (PRIMITIVE & CONVETIONAL)
     $     p1(3,3), c1(3,3)     ! lattive vectors in ANGSTROMS unit
      integer
     $     nat,                 ! number of atoms
     $     num_of_images,       ! number of NEB images
     $     ounit                ! output unit     
      integer
     $     i, j                 ! dummies

      do i=1,3
         do j=1,3
            p1(i,j) = alat*p(i,j)
            c1(i,j) = alat*c(i,j)
         enddo
      enddo

      if (num_of_images .gt. 1)
     $     write(ounit,'('' ANIMSTEPS '',i5)') num_of_images

      write(ounit,'('' CRYSTAL'')')
      write(ounit,'(/,'' PRIMVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((p1(i,j),i=1,3),j=1,3)
      write(ounit,'('' CONVVEC'')')
      write(ounit,'(3(f15.10,2x,f15.10,2x,f15.10,/))')
     $     ((c1(i,j),i=1,3),j=1,3)
      if (num_of_images .eq. 1) then
         write(ounit,'('' PRIMCOORD'')')
         write(ounit,*) nat, 1
      endif
      return
      end

c     -------------------------------------------------
      integer function i_trimleft_white_space(word)
c     trim left white spaces out of word
c     -------------------------------------------------
      character word*(*), auxword*80
 
      ilen=len(word)
      auxword=word
      do i=1,ilen
         if ( word(i:i) .eq. ' ' ) then
            auxword=word(i+1:ilen)
         else
            goto 1
         endif
      enddo
 1    continue
      i_trimleft_white_space=len(word)
      word=auxword(1:i_trimleft_white_space)
      return
      END


c     -----------------------------------------------------------------------
      logical function matches (str1, str2)  
c     .true. if str1 is contained in str2, .false. otherwise
c     This function is taken from PWscf package (www.pwscf.org).
c     -----------------------------------------------------------------------
      implicit none  
      character str1*(*), str2*(*)  
      integer len1, len2, l  
      
      len1 = len(str1)  
      len2 = len(str2)  
      do l = 1, len2 - len1 + 1  
         if ( str1(1:len1) .eq. str2(l:l + len1 - 1) ) then  
            matches = .true.
            return
         endif
      enddo
      matches = .false.  
      return  
      end
