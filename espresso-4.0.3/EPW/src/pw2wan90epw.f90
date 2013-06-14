  !
  ! Copyright (C) 2007-2009 Jesse Noffsinger, Brad Malone, Feliciano Giustino  
  !                                                                            
  ! This file is distributed under the terms of the GNU General Public         
  ! License. See the file `LICENSE' in the root directory of the               
  ! present distribution, or http://www.gnu.org/copyleft.gpl.txt .             
  !                                                                            
  ! Adapted from the code PP/pw2wannier - Quantum-ESPRESSO group               
  !------------------------------------------------------------------------
  subroutine pw2wan90epw 
  !------------------------------------------------------------------------
  ! This is the interface to the Wannier90 code: see http://www.wannier.org
  !
  !
  ! 10/2008  Parellel computation of Amn and Mmn 
  ! 12/2008  Added phase setting of overlap matrix elements
  ! 02/2009  works with standard nk1*nk2*nk3 grids
  ! 12/2009  works with USPP 
  !
  !------------------------------------------------------------------------
#include "f_defs.h"
  use kinds,      only : DP
  use io_global,  only : stdout, ionode, ionode_id
  use mp_global,  only : mpime, inter_pool_comm, my_pool_id
  use mp,         only : mp_bcast
  use cell_base,  only : at, bg
  use lsda_mod,   only : nspin, isk
  use klist,      only : nkstot
  use io_files,   only : nd_nmbr, prefix
  use pwcom
  use epwcom
  use el_phon
  use io_files, only: iunigk
  use wannier
  !
  implicit none
  CHARACTER(LEN=4) :: spin_component
  CHARACTER(len=256) :: outdir
  !
  !
  outdir = './'
  seedname = prefix
  spin_component = 'unpolarized'
  wvfn_formatted = .false.
  reduce_unk= .false.
  !
  !
  !
  WRITE(stdout,*) '     Spin CASE ( default = unpolarized )'
  ispinw = 0
  ikstart = 1
  ikstop  =  nkstot
  iknum   = nkstot
  !
  WRITE(stdout,*) 
  WRITE(stdout,*) '     Initializing Wannier90'
  WRITE(stdout,*) 
  CALL setup_nnkp
  CALL ylm_expansion
  CALL compute_amn_para
  CALL compute_mmn_para
  !
  CALL phases_a_m
  !
  CALL write_band
  !
  IF(write_wfn) call write_plot
  !
  WRITE(stdout,*)
  WRITE(stdout,*) '    Running Wannier90'
  CALL run_wannier
  !
  CALL lib_dealloc
  !
  end subroutine pw2wan90epw
!
!-----------------------------------------------------------------------
subroutine lib_dealloc
  !-----------------------------------------------------------------------
  !
  use wannier
  !
  implicit none
  IF (allocated(m_mat) )     deallocate(m_mat)
  IF (allocated(u_mat) )     deallocate(u_mat)
  IF (allocated(u_mat_opt) ) deallocate(u_mat_opt)
  IF (allocated(a_mat) )     deallocate(a_mat)
  IF (allocated(eigval) )    deallocate(eigval)
  !

end subroutine lib_dealloc
!
!-----------------------------------------------------------------------
subroutine setup_nnkp (  )
  !-----------------------------------------------------------------------
  !
  use io_global, only : stdout, ionode, ionode_id
  use mp_global, only : inter_pool_comm, my_pool_id
  use control_flags, only : iverbosity 
  use kinds,     only : DP
  use constants, only : eps6, tpi
  use cell_base, only : at, bg, alat
  use gvect,     only : g, gg
  use ions_base, only : nat, tau, ityp, atm
  use klist,     only : xk
  use mp,        only : mp_bcast
  use wvfct,     only : nbnd,npwx
  use wannier,   only : num_nnmax, mp_grid, atcart, atsym, kpb, g_kpb, &
                            center_w, alpha_w, l_w, mr_w, r_w, zaxis,  &
                            xaxis, excluded_band, rlatt, glatt, gf,    &
                            csph, ig_, iknum, seedname, kpt_latt, nnb, &
                            num_bands, n_wannier, nexband, nnbx
  use pwcom,    only  :nk1, nk2,nk3
  implicit none
  real(DP) :: g_(3), gg_
  integer  :: ik, ib, ig, iw, ia, indexb, type,i
  real(DP) :: xnorm, znorm, coseno
  integer  :: exclude_bands(nbnd)

  real(DP) :: bohr
  num_nnmax = 32

  ! aam: translations between PW2Wannier90 and Wannier90
  ! pw2wannier90   <==>   Wannier90
  !    nbnd                num_bands_tot
  !    n_wannier           num_wann
  !    num_bands           num_bands
  !    nat                 num_atoms
  !    iknum               num_kpts
  !    rlatt               transpose(real_lattice)
  !    glatt               transpose(recip_lattice)
  !    kpt_latt            kpt_latt
  !    nnb                 nntot
  !    kpb                 nnlist
  !    g_kpb               nncell
  !    mp_grid             mp_grid
  !    center_w            proj_site
  !    l_w,mr_w,r_w        proj_l,proj_m,proj_radial
  !    xaxis,zaxis         proj_x,proj_z
  !    alpha_w             proj_zona
  !    exclude_bands       exclude_bands
  !    atcart              atoms_cart
  !    atsym               atom_symbols

  bohr = 0.5291772108d0
  allocate( atcart(3,nat), atsym(nat) )
  allocate( kpb(iknum,num_nnmax), g_kpb(3,iknum,num_nnmax) )
  allocate( center_w(3,nbnd), alpha_w(nbnd), l_w(nbnd), &
       mr_w(nbnd), r_w(nbnd), zaxis(3,nbnd), xaxis(3,nbnd) )
  allocate( excluded_band(nbnd) )

  ! real lattice (Cartesians, Angstrom)
  rlatt(:,:) = transpose(at(:,:))*alat*bohr
  ! reciprocal lattice (Cartesians, Angstrom)
  glatt(:,:) = transpose(bg(:,:))*tpi/(alat*bohr)
  ! atom co-ordinates in Cartesian co-ords and Angstrom units
  atcart(:,:) = tau(:,:)*bohr*alat
  ! atom symbols
  DO ia=1,nat
     type=ityp(ia)
     atsym(ia)=atm(type)
  ENDDO

#ifdef __PARA
   IF (ionode) then
#endif
      CALL wannier_setup(seedname, mp_grid, iknum,       &  ! input
           rlatt, glatt, kpt_latt, nbnd,                 &  ! input
           nat, atsym, atcart, .false., .false.,         &  ! input
           nnb, kpb, g_kpb, num_bands, n_wannier,        &  ! output
           center_w, l_w, mr_w, r_w, zaxis,              &  ! output
           xaxis, alpha_w, exclude_bands)                   ! output
#ifdef __PARA
   ENDIF
#endif
   
   CALL mp_bcast(nnb,ionode_id)
   CALL mp_bcast(kpb,ionode_id)
   CALL mp_bcast(g_kpb,ionode_id)
   CALL mp_bcast(num_bands,ionode_id)
   CALL mp_bcast(n_wannier,ionode_id)
   CALL mp_bcast(center_w,ionode_id)
   CALL mp_bcast(l_w,ionode_id)
   CALL mp_bcast(mr_w,ionode_id)
   CALL mp_bcast(r_w,ionode_id)
   CALL mp_bcast(zaxis,ionode_id)
   CALL mp_bcast(xaxis,ionode_id)
   CALL mp_bcast(alpha_w,ionode_id)
   CALL mp_bcast(exclude_bands,ionode_id)
   !
   !
   WRITE (stdout,*)
   WRITE (stdout,*) '    Initial Wannier projections'
   WRITE (stdout,*)
   DO i=1,n_wannier
      WRITE (stdout, '(5x,"(",3f10.5,") :  l = ",i3, " mr = ", i3)') center_w(:,i), l_w(i), mr_w(i)
   ENDDO
   !
   WRITE(stdout,'(/,"      - Number of bands is (",i3,")")') num_bands 
   WRITE(stdout,'("      - Number of wannier functions is (",i3,")")') n_wannier 
   !
   allocate( gf(npwx,n_wannier), csph(16,n_wannier) ) 
   !
   excluded_band(1:nbnd)=.false.
   nexband=0
   band_loop: do ib=1,nbnd
      indexb=exclude_bands(ib)
      IF (indexb>nbnd .or. indexb<0) then
         CALL errore('setup_nnkp',' wrong excluded band index ', 1)
      ELSEif (indexb.eq.0) then 
         exit band_loop
      ELSE
         nexband=nexband+1
         excluded_band(indexb)=.true.
      ENDIF
   ENDDO band_loop
 
  IF ( (nbnd-nexband).ne.num_bands ) &
       CALL errore('setup_nnkp',' something wrong with num_bands',1)
 
   DO iw=1,n_wannier
      xnorm = sqrt(xaxis(1,iw)*xaxis(1,iw) + xaxis(2,iw)*xaxis(2,iw) + &
           xaxis(3,iw)*xaxis(3,iw))
      IF (xnorm < eps6) call errore ('setup_nnkp',' |xaxis| < eps ',1)
      znorm = sqrt(zaxis(1,iw)*zaxis(1,iw) + zaxis(2,iw)*zaxis(2,iw) + &
           zaxis(3,iw)*zaxis(3,iw))
      IF (znorm < eps6) call errore ('setup_nnkp',' |zaxis| < eps ',1)
      coseno = (xaxis(1,iw)*zaxis(1,iw) + xaxis(2,iw)*zaxis(2,iw) + &
           xaxis(3,iw)*zaxis(3,iw))/xnorm/znorm
      IF (abs(coseno) > eps6) &
           CALL errore('setup_nnkp',' xaxis and zaxis are not orthogonal !',1)
      IF (alpha_w(iw) < eps6) &
           CALL errore('setup_nnkp',' zona value must be positive', 1)
      ! convert wannier center in cartesian coordinates (in unit of alat)
      CALL cryst_to_cart( 1, center_w(:,iw), at, 1 )
   ENDDO
   WRITE(stdout,*) '     - All guiding functions are given '
   nnbx=0
   nnb=max(nnbx,nnb)
 
   allocate( ig_(iknum,nnb) )
 
   DO ik=1, iknum
      DO ib = 1, nnb
         g_(:) = REAL( g_kpb(:,ik,ib) )
         CALL trnvect (g_, at, bg, 1)
         gg_ = g_(1)*g_(1) + g_(2)*g_(2) + g_(3)*g_(3)
         ig_(ik,ib) = 0
         ig = 1
         DO while  (gg(ig) <= gg_ + eps6) 
            IF ( (abs(g(1,ig)-g_(1)) < eps6) .and.  &
                 (abs(g(2,ig)-g_(2)) < eps6) .and.  &
                 (abs(g(3,ig)-g_(3)) < eps6)  ) ig_(ik,ib) = ig
            ig= ig +1
         ENDDO
      ENDDO
   ENDDO
   WRITE(stdout,*) '     - All neighbours are found '
   WRITE(stdout,*)
   !

end subroutine setup_nnkp
 !
 !-----------------------------------------------------------------------
subroutine run_wannier
  !-----------------------------------------------------------------------
  !
  use io_global, only : ionode, ionode_id
  use ions_base, only : nat
  use mp,        only : mp_bcast
  use cell_base, only : celldm
  use io_files,  only : prefix
  use wannier
  use epwcom,    only : eig_read
  use wvfct,     only : nbnd
  use io_files,  only : find_free_unit

  implicit none

  integer             :: i, ik, ibnd, dummy1, dummy2, ios, iummn
  complex(kind=DP), parameter :: czero = (0.d0, 0.d0), cone = (1.d0, 0.d0), ci = (0.d0, 1.d0)
  character (len=256) :: tempfile
  real(kind=DP), parameter       ::  bohr = 0.5291772108d0, zero = 0.0
  !
  allocate(u_mat(n_wannier,n_wannier,iknum))
  allocate(u_mat_opt(num_bands,n_wannier,iknum))
  allocate(lwindow(num_bands,iknum))
  allocate(wann_centers(3,n_wannier))
  allocate(wann_spreads(n_wannier))
  !
  u_mat_opt = czero
  !
#ifdef __PARA
   IF (ionode) then
#endif
      ! read in external eigenvalues, e.g.  GW
      IF (eig_read) then
         WRITE (6,'(5x,a,i5,a,i5,a)') "Reading electronic eigenvalues (", &
              nbnd, ",", iknum,")"
         tempfile=trim(prefix)//'.eig'
         open(1, file=tempfile, form='formatted', action='read', iostat=ios)
         IF (ios /= 0) call errore ('run_wannier','error opening' // tempfile, 1)
         !
         ! the form should be band, kpt, eigenvalue
         !
         DO ik = 1, iknum
            DO ibnd = 1, nbnd
               read (1,*) dummy1, dummy2, eigval (ibnd,ik)
               IF (dummy1.ne.ibnd) call errore('run_wannier', "Incorrect eigenvalue file", 1)
               IF (dummy2.ne.ik) call errore('run_wannier', "Incorrect eigenvalue file", 1)
            ENDDO
         ENDDO
         close(1)
      ENDIF
      !
      tempfile=trim(prefix)//'.mmn'
      iummn = find_free_unit()
      OPEN(iummn, file=tempfile, iostat=ios, form='unformatted')
      WRITE(iummn) m_mat
      CLOSE(iummn)
     CALL wannier_run(seedname,mp_grid,iknum, &                      ! input
          rlatt, glatt, kpt_latt, num_bands,  &                      ! input
          n_wannier, nnb, nat, atsym,         &                      ! input
          atcart,.false., m_mat,a_mat,eigval, &                      ! input
          u_mat, u_mat_opt, lwindow, wann_centers, &                 ! output
          wann_spreads,spreads)                                      ! output
#ifdef __PARA
   ENDIF
#endif
  !
  CALL mp_bcast(u_mat,ionode_id)
  CALL mp_bcast(u_mat_opt,ionode_id)
  CALL mp_bcast(lwindow,ionode_id)
  CALL mp_bcast(wann_centers,ionode_id)
  CALL mp_bcast(wann_spreads,ionode_id)
  CALL mp_bcast(spreads,ionode_id)
  !
  !
  ! output the results of the wannierization
  !
  WRITE (6,*)
  WRITE (6,*) '    Wannier Function centers (cartesian, alat) and spreads (ang):'
  WRITE (6,*)
  DO i=1,n_wannier
     WRITE (6, '(5x,"(",3f10.5,") :  ",f8.5)') wann_centers(:,i)/celldm(1)/bohr, wann_spreads(i)
  ENDDO
  WRITE (6,*)
  !
  ! store the final minimisation matrix on disk for later use
  !
  CALL write_filukk
  !
end subroutine run_wannier
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine compute_amn_para
!-----------------------------------------------------------------------
!  adapted from compute_amn in pw2wannier90.f90
!  parallelization on k-points has been added
!  10/2008 Jesse Noffsinger UC Berkeley
!
   use io_global,       only : stdout, ionode
#ifdef __PARA
   use mp_global,       only : my_pool_id, npool, intra_pool_comm, inter_pool_comm
   use mp,              only : mp_sum
#endif
   use kinds,           only : DP
   use klist,           only : nkstot, xk, nks
   use wvfct,           only : nbnd, npw, npwx, igk, g2kin
   use units_ph,        only : lrwfc, iuwfc
   use io_files,        only : find_free_unit, nwordwfc, iunwfc
   use gvect,           only : g, ngm, ecutwfc
   use cell_base,       only : tpiba2
   use uspp,            only : nkb, vkb
   use becmod,          only : becp
   use wannier
   use ions_base,       only : nat, ntyp => nsp, ityp, tau
   use uspp_param,      only : upf
   use phcom
   use becmod,          only : calbec 
   implicit none
  
   complex(DP) :: amn, ZDOTC
   complex(DP), allocatable  :: evc_tmp(:)
   complex(DP), allocatable  :: evc(:,:)
   complex(kind=DP), parameter :: czero = (0.d0, 0.d0), cone = (1.d0, 0.d0), ci = (0.d0, 1.d0)
   complex(DP), allocatable :: sgf(:,:)
   integer :: amn_tot, ik, ibnd, ibnd1, iw,nkq, nkq_abs, ipool,ik_g
   logical            :: any_uspp
   real(kind=DP)      :: zero_vect(3)
   !
   any_uspp = ANY( upf(:)%tvanp )
   !
   allocate( a_mat(num_bands,n_wannier,iknum))
   allocate( sgf(npwx,n_wannier), evc_tmp(nbnd))
   allocate( evc(npwx,nbnd))
   !
   ! initialize
   a_mat = czero
   !
   zero_vect = 0.d0
   !
   amn_tot = iknum * nbnd * n_wannier
   WRITE (stdout,'(5x,a)') 'AMN'
   !
   IF (any_uspp) then
      IF (allocated(becp)) deallocate (becp)
      allocate ( becp(nkb,n_wannier))
      CALL init_us_1
   end if
   !
#ifdef __PARA
   WRITE(stdout,'(6x,a,i5,a,i4,a)') 'k points = ',iknum, ' in ', npool, ' pools'
#endif
   DO ik=1,nks
      CALL ktokpmq ( xk(:,ik),zero_vect, +1,ipool,nkq,nkq_abs)
      ik_g = nkq_abs
      !
      WRITE (stdout,'(5x,i8, " of ", i4,a)') ik , nks, ' on ionode'
      CALL flush(6)
      CALL davcio( evc, lrwfc, iuwfc, ik, -1 )
      !
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      CALL generate_guiding_functions(ik)   ! they are called gf(npw,n_wannier)
      !
      !  USPP
      !
      IF(any_uspp) then
         CALL init_us_2 (npw, igk, xk (1, ik), vkb)
         ! below we compute the product of beta functions with trial func.
         CALL calbec(npw, vkb, gf, becp) 
         ! and we use it for the product S|trial_func>
         CALL s_psi (npwx, npw, n_wannier, gf, sgf)  
      ELSE
         sgf(:,:) = gf(:,:)
      ENDIF
      !
      DO iw = 1,n_wannier
         ibnd1 = 0 
         DO ibnd = 1,nbnd
            amn = ZDOTC(npwx,evc(1,ibnd),1,sgf(1,iw),1) 
#ifdef __PARA
            CALL mp_sum(amn, intra_pool_comm)
#endif
            IF (excluded_band(ibnd)) cycle
            ibnd1=ibnd1+1
            a_mat(ibnd1,iw,ik_g) = amn
         END DO !bands
      END DO !wannier fns
   END DO  ! k-points
   deallocate (sgf,csph)
   IF(any_uspp) deallocate (becp)
   !
#ifdef __PARA
   CALL mp_sum(a_mat, inter_pool_comm)
#endif
   !
   !
   !
   WRITE(stdout,*)
   WRITE(stdout,'(5x,a)') 'AMN calculated'

!-----------------------------------------------------------------------
end subroutine compute_amn_para
!-----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
subroutine compute_mmn_para
!-----------------------------------------------------------------------
!
!  adapted from compute_mmn in pw2wannier90.f90
!  parallelization on k-points has been added
!  10/2008 Jesse Noffsinger UC Berkeley
!
   use io_global,       only : stdout, ionode
   use mp_global,       only : my_pool_id, nproc, nproc_pool, &
                               npool, inter_pool_comm, intra_pool_comm
#ifdef __PARA
   use mp,              only: mp_sum
#endif
   use kinds,           only: DP
   use wvfct,           only : nbnd, npw, npwx, igk, g2kin
   use wavefunctions_module, only : psic
   use units_ph,        only : lrwfc, iuwfc
   use gsmooth,         only: nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
   use klist,           only : nkstot, xk, nks
   use io_files,        only : nwordwfc, iunwfc,find_free_unit,prefix
   use gvect,           only : g, ngm, ecutwfc
   use cell_base,       only : tpiba2, omega, alat, tpiba, at, bg
   use ions_base,       only : nat, ntyp => nsp, ityp, tau
   use constants,       only : tpi
   use uspp,            only : nkb, vkb
   use uspp_param,      only : upf, lmaxq, nh
   use becmod,          only : becp, calbec
   use wannier

   implicit none

   integer :: mmn_tot, ik, ikp, ib, npwq, i, m, n, ibnd, jbnd
   integer :: ikb, jkb, ih, jh, na, nt, ijkb0, ind, nbt
   complex(DP), allocatable :: phase(:), aux(:), evc_tmp(:), evc(:,:), evcq(:,:), becp2(:,:), Mkb(:,:)
   complex(DP), allocatable :: qb(:,:,:,:), qgm(:)
   real(DP), allocatable    :: qg(:), ylm(:,:), dxk(:,:)
   integer, allocatable     :: igkq(:)
   real(DP)                 :: xktot(3,nkstot)
   real(DP), dimension(3)   :: zero_vect
   complex(DP)              :: mmn, ZDOTC, phase1
   real(DP)                 :: aa, arg, g_(3)
   logical                  :: any_uspp, exst
   integer                  :: nkq, nkq_abs, ipool
   integer                  :: ik_g, ikp_g, ind0, iummn
   complex(kind=DP), parameter :: czero = (0.d0, 0.d0), cone = (1.d0, 0.d0), ci = (0.d0, 1.d0)

   allocate( phase(nrxxs), aux(npwx), evcq(npwx,nbnd), igkq(npwx) )
   allocate(m_mat(num_bands,num_bands,nnb,iknum))
   allocate(evc(npwx,nbnd), evc_tmp(nbnd))
   
   ! close all the wfc files to allow access for each pool to all wfs
   close (unit = iuwfc,  status = 'keep')
   !
   WRITE (stdout,*)
   WRITE (stdout,'(5x,a)') 'MMN'
   !
   ! Get all the k-vector coords to each pool via xktot
   xktot = 0.d0
   IF (ionode) then
      DO i = 1,nkstot
         xktot(:,i) = xk(:,i)
      ENDDO
   ENDIF
#ifdef __PARA
   CALL mp_sum(xktot, inter_pool_comm)
#endif
   !
   zero_vect = 0.0
   m_mat = czero
   !
   mmn_tot = iknum * nnb * nbnd * nbnd
   !
   !   USPP
   !
   any_uspp = ANY( upf(:)%tvanp )
   !
   IF(any_uspp) then
      CALL init_us_1
      allocate ( becp(nkb,nbnd),becp2(nkb,nbnd))
   end if
   !
   nbt = nnb * iknum
   !
   allocate( qg(nbt) )
   allocate (dxk(3,nbt))
   !
   ind = 0
   DO ik=1,iknum
      DO ib=1,nnb
         ind = ind + 1
         ikp = kpb(ik,ib) 
         !
         g_(:) = REAL( g_kpb(:,ik,ib) )
         CALL trnvect (g_, at, bg, 1)
         dxk(:,ind) = xktot(:,ikp) +g_(:) - xktot(:,ik) 
         qg(ind) = dxk(1,ind)*dxk(1,ind)+dxk(2,ind)*dxk(2,ind)+dxk(3,ind)*dxk(3,ind)
      ENDDO
  ENDDO
  !
  !  USPP
  !
  IF(any_uspp) then
     !
     allocate( ylm(nbt,lmaxq*lmaxq), qgm(nbt) )
     allocate( qb (nkb, nkb, ntyp, nbt) )
     !
     CALL ylmr2 (lmaxq*lmaxq, nbt, dxk, qg, ylm)
     qg(:) = sqrt(qg(:)) * tpiba
     !
     DO nt = 1, ntyp
        IF ( upf(nt)%tvanp ) then
           DO ih = 1, nh (nt)
              DO jh = 1, nh (nt)
                 CALL qvan2 (nbt, ih, jh, nt, qg, qgm, ylm)
                 qb (ih, jh, nt, 1:nbt) = omega * qgm(1:nbt)
              ENDDO
           ENDDO
        ENDIF
     ENDDO
     !
     deallocate (qg, qgm, ylm )
     !
  end if
  !
  !
  allocate( Mkb(nbnd,nbnd) )
  !
#ifdef __PARA
  WRITE(stdout,'(6x,a,i5,a,i4,a)') 'k points = ',iknum, ' in ', npool, ' pools'
#endif
  ! get the first k-point in this pool 
  CALL ktokpmq( xk(:, 1), zero_vect, +1, ipool, nkq, nkq_abs)
  ind0 = (nkq_abs-1) * nnb
  !
  ind = ind0
  DO ik=1,nks 
     CALL ktokpmq( xk(:, ik), zero_vect, +1, ipool, nkq, nkq_abs)
     ik_g = nkq_abs
     !
     WRITE (stdout,'(5x,i8, " of ", i4,a)') ik , nks, ' on ionode'
     CALL flush(6)
     !
     !
     !
     CALL readwfc(my_pool_id+1, ik, evc)
     !
     CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
     !
     !  USPP
     !
     IF(any_uspp) then
        CALL init_us_2 (npw, igk, xk(1,ik), vkb)
        CALL calbec(npw, vkb, evc, becp) 
     end if
     !
     !
     DO ib=1,nnb  ! loop on finite diff vectors
        ind = ind + 1
        !
        ikp = kpb(ik_g,ib)
        ikp_g = ikp 
        !
        CALL ktokpmq( xk(:, ik), xktot(:,ikp_g)-xk(:,ik), +1, ipool, nkq, nkq_abs)
        !
        ! read wfc at k+b
        !
        CALL readwfc( ipool, nkq, evcq)
        !
        CALL gk_sort (xktot(1,ikp_g), ngm, g, ecutwfc / tpiba2, npwq, igkq, g2kin)
        !
        ! compute the phase
        !
        phase(:) = (0.d0,0.d0)
        IF ( ig_(ik_g,ib)>0) phase( nls(ig_(ik_g,ib)) ) = (1.d0,0.d0)
        CALL cft3s (phase, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
        !
        !  USPP
        !
        IF(any_uspp) then
           CALL init_us_2 (npwq, igkq, xktot(1,ikp_g), vkb)
           ! below we compute the product of beta functions with |psi> 
           CALL calbec (npwq, vkb, evcq, becp2)  
        end if
        !
        !
        Mkb(:,:) = (0.0d0,0.0d0) 
        !
        IF (any_uspp) then
           ijkb0 = 0
           DO nt = 1, ntyp
              IF ( upf(nt)%tvanp ) then
                 DO na = 1, nat
                    !
                    arg = DOT_PRODUCT( dxk(:,ind), tau(:,na) ) * tpi 
                    phase1 = CMPLX ( COS(arg), -SIN(arg) )
                    !
                    IF ( ityp(na) == nt ) then
                       DO jh = 1, nh(nt)
                          jkb = ijkb0 + jh
                          DO ih = 1, nh(nt)
                             ikb = ijkb0 + ih
                             !
                             DO m = 1,nbnd
                                DO n = 1,nbnd
                                 Mkb(m,n) = Mkb(m,n) + &
                                      phase1 * qb(ih,jh,nt,ind) * &
                                      conjg( becp(ikb,m) ) * becp2(jkb,n) 
                              ENDDO
                           ENDDO
                        ENDDO !ih
                     ENDDO !jh
                        ijkb0 = ijkb0 + nh(nt)
                     ENDIF  !ityp
                  ENDDO  !nat 
               ELSE  !tvanp
                  DO na = 1, nat
                     IF ( ityp(na) == nt ) ijkb0 = ijkb0 + nh(nt)
                  ENDDO
               ENDIF !tvanp
            ENDDO !ntyp
         end if ! any_uspp
         !
         !
         DO m=1,nbnd  ! loop on band m
            psic(:) = (0.d0, 0.d0)
            psic(nls (igk (1:npw) ) ) = evc (1:npw, m)
            CALL cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
            psic(1:nrxxs) = psic(1:nrxxs) * phase(1:nrxxs)
            CALL cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, -2)
            aux(1:npwq) = psic(nls (igkq(1:npwq) ) )
            aa = 0.d0
            !
            !
            DO n=1,nbnd   ! loop on band n
              !
              mmn = ZDOTC (npwq, aux,1,evcq(1,n),1)
              !
              !  Mkb(m,n) = Mkb(m,n) + \sum_{ijI} qb_{ij}^I * e^-i(b*tau_I)
              !             <psi_m,k1| beta_i,k1 > < beta_j,k2 | psi_n,k2 > 
              !
#ifdef __PARA
              CALL mp_sum(mmn, intra_pool_comm)
#endif
              Mkb(m,n) = mmn + Mkb(m,n)
              !
              aa = aa + abs(mmn)**2
              !
            END DO ! n
         END DO   ! m

         DO n=1,nbnd
            IF (excluded_band(n)) cycle
            DO m=1,nbnd
               IF (excluded_band(m)) cycle
               m_mat(m,n,ib,ik_g)=Mkb(m,n)
               !
            ENDDO
         ENDDO
         !
         !
      END DO !ib
   END DO  !ik
   !
#ifdef __PARA
   CALL mp_sum(m_mat, inter_pool_comm)
#endif
   !
   ! 
   deallocate (Mkb, dxk, phase, aux, evcq, evc, evc_tmp, igkq)
   IF(any_uspp) then 
      deallocate (becp, becp2, qb)
      ! the line from allocate_epwq
      allocate ( becp(nkb, nbnd) )
   ENDIF
   !
   !
#ifdef __PARA
   IF (ionode) then
#endif
   !
#ifdef __PARA
   ENDIF
#endif
   !
   WRITE(stdout,*)
   WRITE(stdout,'(5x,a)') 'MMN calculated'
   !
   ! reopen wfc here, leaving unit=20 in the same state
   iuwfc = 20
   CALL diropn(iuwfc,'wfc',lrwfc,exst)  
   !
end subroutine compute_mmn_para
!
!------------------------------------------------------------------------
!
!-----------------------------------------------------------------------
subroutine compute_pmn_para
!-----------------------------------------------------------------------
!  adapted from compute_amn_para
!  06/2010  Jesse Noffsinger
!  
!
   use io_global,       only : stdout, ionode
#ifdef __PARA
   use mp_global,       only : my_pool_id, npool, intra_pool_comm, inter_pool_comm
   use mp,              only : mp_sum
#endif
   use kinds,           only : DP
   use klist,           only : nkstot, xk, nks, wk
   use wvfct,           only : nbnd, npw, npwx, igk, g2kin, wg
   use units_ph,        only : lrwfc, iuwfc
   use io_files,        only : find_free_unit, nwordwfc, iunwfc
   use gvect,           only : g, ngm, ecutwfc
   use cell_base,       only : tpiba2, tpiba
   use uspp,            only : nkb, vkb
   use becmod,          only : becp
   use wannier
   use ions_base,       only : nat, ntyp => nsp, ityp, tau
   use uspp_param,      only : upf
   use phcom
   use becmod,          only : calbec 
   use el_phon,         only : dmec, et_ks, et_mb
   use epwcom,          only : eig_read
   implicit none
  
   complex(DP), allocatable  :: evc(:,:)
   complex(kind=DP), parameter :: czero = (0.d0, 0.d0), cone = (1.d0, 0.d0), ci = (0.d0, 1.d0)
   integer :: ik, ibnd, ig, jbnd, iband1
   real(kind=DP)      :: zero_vect(3)
   real(kind=DP), parameter   :: ryd2ev = 13.6058
   COMPLEX(DP)  :: dipole_aux(3,nbnd,nbnd), caux 
   !
   !
   allocate( dmec(3,nbnd,nbnd,nks))
   allocate( evc(npwx,nbnd))
   !
   ! initialize
   dmec = czero
   dipole_aux(:,:,:) = (0.0_DP,0.0_DP)
   !
   zero_vect = 0.d0
   !
   DO ik=1,nks
      !
      ! read wfc
      CALL davcio( evc, lrwfc, iuwfc, ik, -1 )
      !
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      !
      dipole_aux = (0.d0, 0.d0)
      DO jbnd = 1,nbnd
         DO ibnd = 1,nbnd
            !
            IF ( ibnd .eq. jbnd ) cycle
            !
            ! taken from PP/epsilon.f90
            DO  ig=1,npwx
               IF (igk(ig) .gt. SIZE(g,2) .or. igk(ig).lt.1) cycle
               !
               caux= conjg(evc(ig,ibnd))*evc(ig,jbnd) 
               !
               dipole_aux(:,ibnd,jbnd) = dipole_aux(:,ibnd,jbnd) + &
                    ( g(:,igk(ig)) ) * caux
               !
            ENDDO
               !
         END DO !bands i
      END DO ! bands j
      ! metal diagonal part
     DO iband1 = 1, nbnd
        DO  ig=1,npwx
           IF (igk(ig) .gt. SIZE(g,2) .or. igk(ig).lt.1) cycle
          !
          caux= conjg(evc(ig,iband1))*evc(ig,iband1) 
          !
          dipole_aux(:,iband1,iband1) = dipole_aux(:,iband1,iband1) + &
                                        ( g(:,igk(ig))+ xk(:,ik) ) * caux
          !
        ENDDO
     ENDDO
     ! need to divide by 2pi/a to fix the units
     dmec(:,:,:,ik) = dipole_aux(:,:,:) * tpiba
     !
  END DO  ! k-points
  !
  !
  WRITE(stdout,'(/5x,a)') 'Dipole matrix elements calculated'
  WRITE(stdout,*)
  !
!-----------------------------------------------------------------------
end subroutine compute_pmn_para
!-----------------------------------------------------------------------
!!
!-----------------------------------------------------------------------
subroutine write_filukk
!-----------------------------------------------------------------------
!
!  Here we compute and write out the final ukk matrix which is used by
!  epw.x to localize the electron wavefuctions (and therefore the ep-matrix 
!  elements)
!  10/2008 Jesse Noffsinger UC Berkeley
!  07/2010 Fixed the rotation for ndimwin when lower bands are not included
!
   use io_files,     only : find_free_unit, prefix
   use kinds,        only: DP
   use io_global,    only : ionode
   use wvfct,        only : nbnd
   use wannier,      only : n_wannier, iknum, u_mat, u_mat_opt, lwindow
   use epwcom,       only : filukk
   use io_files,     only : find_free_unit
   !
   implicit none
   !
   complex(kind=DP), parameter :: czero = (0.d0, 0.d0), cone = (1.d0, 0.d0)
   complex(kind=DP), allocatable :: u_kc(:,:,:)
   integer :: iuukk, jbnd, k_wan, ik, ndimwin(iknum)
   !
   !
   !
#ifdef __PARA
   IF (ionode) then
#endif
      !
      ndimwin(:) = 0
      DO ik = 1, iknum
         DO jbnd = 1, nbnd
            IF (lwindow(jbnd,ik)) ndimwin(ik) = ndimwin(ik) +1
         ENDDO
      ENDDO
      allocate( u_kc(nbnd, n_wannier, iknum) )
      u_kc = czero
      !
      ! get the final rotation matrix, which is the product of the optimal
      ! subspace and the rotation among the n_wannier wavefunctions
      DO ik = 1, iknum
         !
         u_kc(1:ndimwin(ik),1:n_wannier,ik) = &
              matmul (u_mat_opt(1:ndimwin(ik),:,ik), u_mat(:,1:n_wannier,ik))
         !
      ENDDO
      !
      iuukk = find_free_unit()
      !
      open (unit = iuukk, file = filukk, form = 'formatted')
      DO ik = 1,iknum
         DO jbnd = 1, nbnd
            DO k_wan = 1, n_wannier
               WRITE (iuukk,*) u_kc(jbnd,k_wan,ik)
            ENDDO
         ENDDO
      ENDDO
      close (iuukk)
      IF ( allocated(u_kc) ) deallocate(u_kc)
#ifdef __PARA
   ENDIF
#endif
   !
!-----------------------------------------------------------------------
end subroutine write_filukk
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
subroutine phases_a_m
!-----------------------------------------------------------------------
! will set phases here on the matrices.  Should not affect the spreads and
! centers found in w90, but will leave the u_mat_opt and u_mat to reflect the
! known phases
!
!
   use mp_global,       only : inter_pool_comm, my_pool_id
   use mp,              only : mp_sum, mp_barrier
   use kinds,           only : DP
   use io_global,       only : ionode
   use klist,           only : nkstot, xk, nks
   use wvfct,           only : nbnd
   use wannier,         only : a_mat, m_mat, n_wannier, nnb, kpb, iknum
   use el_phon,         only : umat, umat_all

   implicit none

   complex(kind=DP), parameter :: czero = (0.d0, 0.d0), cone = (1.d0, 0.d0), ci = (0.d0, 1.d0)
   integer                      ::  ik, ipool, ib, ikb, i,nkq, &
                                    ik_g
   real(kind=DP)                ::  xktot(3,nkstot)
   complex(kind=DP), allocatable ::  a_mat_tmp(:,:,:), m_mn_tmp1(:,:), &
       m_mn_tmp2(:,:), m_mn_tmp3(:,:,:,:)
   real(DP), dimension(3)   :: zero_vect

   xktot = 0.d0
   IF (ionode) then
      DO i = 1,nkstot
         xktot(:,i) = xk(:,i)
      ENDDO
   ENDIF
   CALL mp_sum(xktot, inter_pool_comm)
   !
   allocate(m_mn_tmp3(nbnd,nbnd,nnb,iknum) )
   allocate(a_mat_tmp(nbnd,n_wannier,iknum),m_mn_tmp1(nbnd,nbnd), &
        m_mn_tmp2(nbnd,nbnd))
   !
   ! zero all temporary/work quantities
   !
   zero_vect = 0.0
   a_mat_tmp = czero
   m_mn_tmp1 = czero
   m_mn_tmp2 = czero
   m_mn_tmp3 = czero
   !
   DO ik=1,nks
      CALL ktokpmq ( xk(:,ik),zero_vect, +1,ipool,nkq,ik_g)
      !
      !  GF_n are the guiding functions which are our initial guesses 
      !  Amn(k) = <psi_k,m|GF_n>.  
      !  We want U(k)^dagger<psi_k,m|GF_m>
      !
       CALL zgemm ('c', 'n', nbnd, n_wannier, nbnd, cone, umat(:,:,ik), & 
             nbnd, a_mat(:,:,ik_g), nbnd, czero, a_mat_tmp(:,:,ik_g), nbnd)
      !
    ENDDO
   CALL mp_sum(a_mat_tmp, inter_pool_comm)
   !
   a_mat(:,:,:) = a_mat_tmp(:,:,:)
   !
   !
   DO ik=1,nks
      CALL ktokpmq ( xk(:,ik),zero_vect, +1,ipool,nkq,ik_g)
      !
      DO ib = 1,nnb
         ikb = kpb(ik_g,ib)
         !
         ! Mmn(k,k+b)  = <psi_k_m| psi_(k+b)_n> so we need
         !  (U(k)^dagger <psi_k_m| ) * (|psi_k+b_n> U(k+b)
         ! = U(k)^dagger (M_mn) = m_mat_tmp, 
         ! Mmn(k,k+b)' = m_mat_tmp*U(k+b) 
         !
         CALL zgemm ('c', 'n', nbnd, nbnd, nbnd, cone, umat(:,:,ik), & 
              nbnd, m_mat(:,:,ib,ik_g), nbnd, czero, m_mn_tmp1(:,:), nbnd)
         CALL zgemm ('n', 'n', nbnd, nbnd, nbnd, cone, m_mn_tmp1(:,:), & 
              nbnd, umat_all(:,:,ikb), nbnd, czero, m_mn_tmp2(:,:), nbnd)
         !         !
         !         m_mn_tmp1 = matmul ( conjg( transpose (umat(:,:,ik) )), m_mat(:,:,ib,ik_g ) )
         !         m_mn_tmp2 = matmul ( m_mn_tmp1, umat_g(:,:,ikb) )
         !
         m_mn_tmp3(:,:,ib,ik_g) = m_mn_tmp2
      ENDDO
   ENDDO
   CALL mp_sum(m_mn_tmp3, inter_pool_comm)

   m_mat(:,:,:,:) = m_mn_tmp3(:,:,:,:)
   !
   deallocate(m_mn_tmp3, m_mn_tmp2, m_mn_tmp1, a_mat_tmp)
   !
   !
!-----------------------------------------------------------------------
end subroutine phases_a_m
!-----------------------------------------------------------------------
!
!---------------------------------------
!
! subroutines below are largely unchanged
! from pw2wannier90.x of esp-4.0.1
!
!---------------------------------------
!
subroutine generate_guiding_functions(ik)
   !
   use mp_global,          only : intra_pool_comm
   use mp,                 only : mp_sum
   use io_global,          only : stdout
   use constants,          only : pi, tpi, fpi, eps8
   use wvfct,              only : npw, g2kin, igk
   use gvect,              only : ig1, ig2, ig3, g
   use cell_base,          only : tpiba2, omega, tpiba
   use wannier
   use klist,              only : xk 
   use cell_base,          only : bg

   implicit none

   integer, parameter :: lmax=3, lmax2=(lmax+1)**2
   integer :: iw, ig, ik, l
   integer :: lm, iig
   real(DP) :: arg, anorm
   complex(DP) :: ZDOTC, lphase
   real(DP), allocatable :: gk(:,:), qg(:), ylm(:,:), radial(:,:)
   complex(DP), allocatable :: sk(:) 
   !
   allocate( gk(3,npw), qg(npw), ylm(npw,lmax2), sk(npw), radial(npw,0:lmax) )
   !
   DO ig = 1, npw
      gk (1,ig) = xk(1, ik) + g(1, igk(ig) )
      gk (2,ig) = xk(2, ik) + g(2, igk(ig) )
      gk (3,ig) = xk(3, ik) + g(3, igk(ig) )
      qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
   ENDDO

   CALL ylmr2 (lmax2, npw, gk, qg, ylm)
   ! define qg as the norm of (k+g) in a.u.
   qg(:) = sqrt(qg(:)) * tpiba

   DO iw = 1, n_wannier
      !
      gf(:,iw) = (0.d0,0.d0)

      CALL radialpart(npw, qg, alpha_w(iw), r_w(iw), lmax, radial) 

      DO lm = 1, lmax2
         IF ( abs(csph(lm,iw)) < eps8 ) cycle
         l = int (sqrt( lm-1.d0))
         lphase = (0.d0,-1.d0)**l
         !
         DO ig=1,npw
            gf(ig,iw) = gf(ig,iw) + csph(lm,iw) * ylm(ig,lm) * radial(ig,l) * lphase
         END DO !ig
      END DO ! lm
      DO ig=1,npw
         iig = igk(ig)
         arg = ( gk(1,ig)*center_w(1,iw) + gk(2,ig)*center_w(2,iw) + &
                                           gk(3,ig)*center_w(3,iw) ) * tpi
         ! center_w are cartesian coordinates in units of alat 
         sk(ig) = CMPLX(cos(arg), -sin(arg) )
         gf(ig,iw) = gf(ig,iw) * sk(ig) 
      END DO
      anorm = REAL(ZDOTC(npw,gf(1,iw),1,gf(1,iw),1))
      CALL mp_sum(anorm, intra_pool_comm)
      gf(:,iw) = gf(:,iw) / dsqrt(anorm)
   END DO
   !
   deallocate ( gk, qg, ylm, sk, radial)

end subroutine generate_guiding_functions

subroutine write_band
   use io_global,  only : stdout, ionode
   use wvfct, only : nbnd, et
   use klist, only : nkstot
   use constants, only: rytoev
   use io_files, only : find_free_unit
   use wannier

   implicit none

   integer ik, ibnd, ibnd1, ikevc

   allocate(eigval(num_bands,iknum))

   DO ik=ikstart,ikstop
      ikevc = ik - ikstart + 1
      ibnd1=0
      DO ibnd=1,nbnd
         IF (excluded_band(ibnd)) cycle
         ibnd1=ibnd1 + 1
         eigval(ibnd1,ikevc) = et(ibnd,ik)*rytoev
      END DO
   END DO

end subroutine write_band

!---------------------------------------
subroutine write_plot
!---------------------------------------
!
! JN 06/2009: 
! added a couple of calls -- now works with multiple
! pools/procs (but one proc per pool)
!
   use io_global,  only : stdout, ionode
   use wvfct, only : nbnd, npw, igk, g2kin
   use wavefunctions_module, only : evc, psic
   use units_ph,        only : iuwfc
   use io_files, only : find_free_unit, nwordwfc
   use wannier
   use gsmooth,         only : nls, nrxxs, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s
   use klist,           only : nkstot, xk, nks
   use gvect,           only : g, ngm, ecutwfc
   use cell_base,       only : tpiba2
   use fft_base,        only : cgather_smooth

   implicit none
   integer ik, ibnd, ibnd1, j, spin, nkq, nkq_abs, ipool
   character*20 wfnname

   ! aam: 1/5/06: for writing smaller unk files 
   integer :: n1by2,n2by2,n3by2,i,k,index,pos
   COMPLEX(DP),allocatable :: psic_small(:)   
   real(kind=DP)      :: zero_vect(3)
   !-------------------------------------------!
#ifdef __PARA
   integer nxxs
   COMPLEX(DP),allocatable :: psic_all(:)
   nxxs = nrx1s * nrx2s * nrx3s
   allocate(psic_all(nxxs) )
#endif
   !
   zero_vect = 0.d0
   WRITE (stdout,*) '    Writing out UNK plot files'
   !
   IF (reduce_unk) then
      WRITE(stdout,'(3(a,i5))') 'nr1s =',nr1s,'nr2s=',nr2s,'nr3s=',nr3s
      n1by2=(nr1s+1)/2;n2by2=(nr2s+1)/2;n3by2=(nr3s+1)/2
      WRITE(stdout,'(3(a,i5))') 'n1by2=',n1by2,'n2by2=',n2by2,'n3by2=',n3by2
      allocate(psic_small(n1by2*n2by2*n3by2))   
   ENDIF
   !
   DO ik=1,nks
      CALL ktokpmq ( xk(:,ik),zero_vect, +1,ipool,nkq,nkq_abs)
      !
      iun_plot = find_free_unit()
      spin=ispinw
      IF(ispinw.eq.0) spin=1
      WRITE(wfnname,200) nkq_abs, spin
200   format ('UNK',i5.5,'.',i1)
      !
      IF(wvfn_formatted) then
         open (unit=iun_plot, file=wfnname,form='formatted')
         IF (reduce_unk) then
            WRITE(iun_plot,*)  n1by2,n2by2,n3by2, nkq_abs, nbnd-nexband
         ELSE
            WRITE(iun_plot,*)  nr1s,nr2s,nr3s, nkq_abs, nbnd-nexband
         ENDIF
      ELSE
         open (unit=iun_plot, file=wfnname,form='unformatted')
         IF (reduce_unk) then
            WRITE(iun_plot)  n1by2,n2by2,n3by2, nkq_abs, nbnd-nexband
         ELSE
            WRITE(iun_plot)  nr1s,nr2s,nr3s, nkq_abs, nbnd-nexband
         ENDIF
      ENDIF
      !
      CALL davcio (evc, nwordwfc, iuwfc, ik, -1 )
      CALL gk_sort (xk(1,ik), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
      !
      ibnd1 = 0
      DO ibnd=1,nbnd
         IF (excluded_band(ibnd)) cycle
         ibnd1=ibnd1 + 1
         psic(:) = (0.d0, 0.d0)
         psic(nls (igk (1:npw) ) ) = evc (1:npw, ibnd)
         CALL cft3s (psic, nr1s, nr2s, nr3s, nrx1s, nrx2s, nrx3s, +2)
         IF (reduce_unk) pos=0
#ifdef __PARA
         CALL cgather_smooth(psic,psic_all)
         !
         IF (reduce_unk) then
            DO k=1,nr3s,2
               DO j=1,nr2s,2
                  DO i=1,nr1s,2
                     index = (k-1)*nr3s*nr2s + (j-1)*nr2s + i
                     pos=pos+1
                     psic_small(pos) = psic_all(index) 
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         !
         IF(wvfn_formatted) then
            IF (reduce_unk) then
               WRITE (iun_plot,'(2ES20.10)') (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot,*) (psic_all(j),j=1,nr1s*nr2s*nr3s)
            ENDIF
         ELSE
            IF (reduce_unk) then
               WRITE (iun_plot) (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot) (psic_all(j),j=1,nr1s*nr2s*nr3s)
            ENDIF
         ENDIF
#else
         IF (reduce_unk) then
            DO k=1,nr3s,2
               DO j=1,nr2s,2
                  DO i=1,nr1s,2
                     index = (k-1)*nr3s*nr2s + (j-1)*nr2s + i
                     pos=pos+1
                     psic_small(pos) = psic(index) 
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         IF(wvfn_formatted) then 
            IF (reduce_unk) then
               WRITE (iun_plot,'(2ES20.10)') (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot,*) (psic(j),j=1,nr1s*nr2s*nr3s)
            ENDIF
         ELSE
            IF (reduce_unk) then
               WRITE (iun_plot) (psic_small(j),j=1,n1by2*n2by2*n3by2)
            ELSE
               WRITE (iun_plot) (psic(j),j=1,nr1s*nr2s*nr3s)
            ENDIF
         ENDIF
#endif
      END DO !ibnd
      !
      close (iun_plot)
      !
   END DO  !ik
   !
   IF (reduce_unk) deallocate(psic_small)   
   !
#ifdef __PARA
   deallocate( psic_all )
#endif

!---------------------------------------
end subroutine write_plot
!---------------------------------------
subroutine ylm_expansion 
   use io_global,  only : stdout
   use kinds, only :  DP
   use random_numbers,       only : rndm
   use wannier
   implicit none
   ! local variables
   integer, parameter :: lmax2=16
   integer ::  i, ir, iw
   real(DP) :: capel
   real(DP), allocatable :: r(:,:), rr(:), rp(:,:), ylm_w(:), ylm(:,:), mly(:,:)
   real(DP) :: u(3,3)

   allocate (r(3,lmax2), rp(3,lmax2), rr(lmax2), ylm_w(lmax2))
   allocate (ylm(lmax2,lmax2), mly(lmax2,lmax2) )

   ! generate a set of nr=lmax2 random vectors
   DO ir=1,lmax2
      DO i=1,3
         r(i,ir) = rndm() -0.5d0
      END DO
   END DO
   rr(:) = r(1,:)*r(1,:) + r(2,:)*r(2,:) + r(3,:)*r(3,:)
   !- compute ylm(ir,lm)
   CALL ylmr2(lmax2, lmax2, r, rr, ylm)
   !- store the inverse of ylm(ir,lm) in mly(lm,ir)
   CALL invmat(lmax2, ylm, mly, capel)
   !- check that r points are independent
   CALL check_inverse(lmax2, ylm, mly)

   DO iw=1, n_wannier

      !- define the u matrix that rotate the reference frame
      CALL set_u_matrix (xaxis(:,iw),zaxis(:,iw),u)
      !- find rotated r-vectors 
      rp(:,:) = matmul ( u(:,:) , r(:,:) )
      !- set ylm funtion according to wannier90 (l,mr) indexing in the rotaterd points
      CALL ylm_wannier(ylm_w,l_w(iw),mr_w(iw),rp,lmax2) 

      csph(:,iw) = matmul (mly(:,:), ylm_w(:))

!      write (stdout,*) 
!      write (stdout,'(2i4,2(2x,3f6.3))') l_w(iw), mr_w(iw), xaxis(:,iw), zaxis(:,iw)
!      write (stdout,'(16i6)')   (lm, lm=1,lmax2)
!      write (stdout,'(16f6.3)') (csph(lm,iw), lm=1,lmax2)

   END DO
   deallocate (r, rp, rr, ylm_w, ylm, mly )

end subroutine ylm_expansion

subroutine check_inverse(lmax2, ylm, mly)
   use kinds, only :  DP
   use constants, only :  eps8
   implicit none
   ! I/O variables
   integer :: lmax2
   real(DP) :: ylm(lmax2,lmax2), mly(lmax2,lmax2)
   ! local variables
   real(DP), allocatable :: uno(:,:)
   real(DP) :: capel
   integer :: lm
   !
   allocate (uno(lmax2,lmax2) )
   uno = matmul(mly, ylm)
   capel = 0.d0
   DO lm = 1, lmax2
      uno(lm,lm) = uno(lm,lm) - 1.d0
   END DO
   capel = capel + SUM ( abs(uno(1:lmax2,1:lmax2) ) )
!   write (stdout,*) "capel = ", capel
   IF (capel > eps8) call errore('ylm_expansion', &
                    ' inversion failed: r(*,1:nr) are not all independent !!',1)
   deallocate (uno)
end subroutine check_inverse
   
subroutine set_u_matrix(x,z,u)
   use kinds, only :  DP
   use constants, only : eps8
   implicit none
   ! I/O variables
   real(DP) :: x(3),z(3),u(3,3)
   ! local variables
   real(DP) :: xx, zz, y(3), coseno

   xx = sqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))
   IF (xx < eps8) call errore ('set_u_matrix',' |xaxis| < eps ',1)
!   x(:) = x(:)/xx
   zz = sqrt(z(1)*z(1) + z(2)*z(2) + z(3)*z(3))
   IF (zz < eps8) call errore ('set_u_matrix',' |zaxis| < eps ',1)
!   z(:) = z(:)/zz

   coseno = (x(1)*z(1) + x(2)*z(2) + x(3)*z(3))/xx/zz
   IF (abs(coseno) > eps8) call errore('set_u_matrix',' xaxis and zaxis are not orthogonal !',1)

   y(1) = (z(2)*x(3) - x(2)*z(3))/xx/zz
   y(2) = (z(3)*x(1) - x(3)*z(1))/xx/zz
   y(3) = (z(1)*x(2) - x(1)*z(2))/xx/zz

   u(1,:) = x(:)/xx
   u(2,:) = y(:)
   u(3,:) = z(:)/zz

!   write (stdout,'(3f10.7)') u(:,:)


end subroutine set_u_matrix

subroutine ylm_wannier(ylm,l,mr,r,nr) 
!
! this routine returns in ylm(r) the values at the nr points r(1:3,1:nr) 
! of the spherical harmonic identified  by indices (l,mr) 
! in table 3.1 of the wannierf90 specification.
! 
! No reference to the particular ylm ordering internal to quantum-espresso
! is assumed. 
!
! If ordering in wannier90 code is changed or extended this should be the 
! only place to be modified accordingly
!
   use kinds, only :  DP
   use constants, only : pi, fpi, eps8
   implicit none
! I/O variables
!
   integer :: l, mr, nr
   real(DP) :: ylm(nr), r(3,nr)
!
! local variables
!
   real(DP), external :: s, p_z,px,py, dz2, dxz, dyz, dx2my2, dxy, &
                        fz3, fxz2, fyz2, fzx2my2, fxyz, fxx2m3y2, fy3x2my2
   real(DP) :: rr, cost, phi
   integer :: ir
   real(DP) :: bs2, bs3, bs6, bs12
   bs2 = 1.d0/sqrt(2.d0)
   bs3=1.d0/sqrt(3.d0)
   bs6 = 1.d0/sqrt(6.d0)
   bs12 = 1.d0/sqrt(12.d0)
!
   IF (l > 3 .OR. l < -5 ) call errore('ylm_wannier',' l out of range ', 1)
   IF (l>=0) then
      IF (mr < 1 .OR. mr > 2*l+1) call errore('ylm_wannier','mr out of range' ,1)
   ELSE
      IF (mr < 1 .OR. mr > abs(l)+1 ) call errore('ylm_wannier','mr out of range',1)
   end if

   DO ir=1, nr
      rr = sqrt( r(1,ir)*r(1,ir) +  r(2,ir)*r(2,ir) + r(3,ir)*r(3,ir) )
      IF (rr < eps8) call errore('ylm_wannier',' rr too small ',1)

      cost =  r(3,ir) / rr
      !
      !  beware the arc tan, it is defined modulo pi
      !
      IF (r(1,ir) > eps8) then
         phi = atan( r(2,ir)/r(1,ir) )
      ELSE if (r(1,ir) < -eps8 ) then
         phi = atan( r(2,ir)/r(1,ir) ) + pi
      ELSE
         phi = sign( pi/2.d0,r(2,ir) )
      end if

    
      IF (l==0) then   ! s orbital
                    ylm(ir) = s()  
      end if
      IF (l==1) then   ! p orbitals
         IF (mr==1) ylm(ir) = p_z(cost) 
         IF (mr==2) ylm(ir) = px(cost,phi)
         IF (mr==3) ylm(ir) = py(cost,phi)
      end if
      IF (l==2) then   ! d orbitals
         IF (mr==1) ylm(ir) = dz2(cost)
         IF (mr==2) ylm(ir) = dxz(cost,phi)
         IF (mr==3) ylm(ir) = dyz(cost,phi)
         IF (mr==4) ylm(ir) = dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = dxy(cost,phi)
      ENDIF
      IF (l==3) then   ! f orbitals
         IF (mr==1) ylm(ir) = fz3(cost)
         IF (mr==2) ylm(ir) = fxz2(cost,phi)
         IF (mr==3) ylm(ir) = fyz2(cost,phi)
         IF (mr==4) ylm(ir) = fzx2my2(cost,phi)
         IF (mr==5) ylm(ir) = fxyz(cost,phi)
         IF (mr==6) ylm(ir) = fxx2m3y2(cost,phi)
         IF (mr==7) ylm(ir) = fy3x2my2(cost,phi)
      ENDIF
      IF (l==-1) then  !  sp hybrids
         IF (mr==1) ylm(ir) = bs2 * ( s() + px(cost,phi) ) 
         IF (mr==2) ylm(ir) = bs2 * ( s() - px(cost,phi) ) 
      end if
      IF (l==-2) then  !  sp2 hybrids 
         IF (mr==1) ylm(ir) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s() +2.d0*bs6*px(cost,phi) 
      end if
      IF (l==-3) then  !  sp3 hybrids
         IF (mr==1) ylm(ir) = 0.5d0*(s()+px(cost,phi)+py(cost,phi)+p_z(cost))
         IF (mr==2) ylm(ir) = 0.5d0*(s()+px(cost,phi)-py(cost,phi)-p_z(cost))
         IF (mr==3) ylm(ir) = 0.5d0*(s()-px(cost,phi)+py(cost,phi)-p_z(cost))
         IF (mr==4) ylm(ir) = 0.5d0*(s()-px(cost,phi)-py(cost,phi)+p_z(cost))
      end if
      IF (l==-4) then  !  sp3d hybrids
         IF (mr==1) ylm(ir) = bs3*s()-bs6*px(cost,phi)+bs2*py(cost,phi)
         IF (mr==2) ylm(ir) = bs3*s()-bs6*px(cost,phi)-bs2*py(cost,phi)
         IF (mr==3) ylm(ir) = bs3*s() +2.d0*bs6*px(cost,phi) 
         IF (mr==4) ylm(ir) = bs2*p_z(cost)+bs2*dz2(cost)
         IF (mr==5) ylm(ir) =-bs2*p_z(cost)+bs2*dz2(cost)
      end if
      IF (l==-5) then  ! sp3d2 hybrids
         IF (mr==1) ylm(ir) = bs6*s()-bs2*px(cost,phi)-bs12*dz2(cost)+.5*dx2my2(cost,phi)
         IF (mr==2) ylm(ir) = bs6*s()+bs2*px(cost,phi)-bs12*dz2(cost)+.5*dx2my2(cost,phi)
         IF (mr==3) ylm(ir) = bs6*s()-bs2*py(cost,phi)-bs12*dz2(cost)-.5*dx2my2(cost,phi)
         IF (mr==4) ylm(ir) = bs6*s()+bs2*py(cost,phi)-bs12*dz2(cost)-.5*dx2my2(cost,phi)
         IF (mr==5) ylm(ir) = bs6*s()-bs2*p_z(cost)+bs3*dz2(cost)
         IF (mr==6) ylm(ir) = bs6*s()+bs2*p_z(cost)+bs3*dz2(cost)
      end if

   END DO


end subroutine ylm_wannier

!======== l = 0 =====================================================================
function s()
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) :: s
   s = 1.d0/ sqrt(fpi)
end function s
!======== l = 1 =====================================================================
function p_z(cost)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::p_z, cost
   p_z =  sqrt(3.d0/fpi) * cost
end function p_z
function px(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::px, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   px =  sqrt(3.d0/fpi) * sint * cos(phi)

end function px
function py(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::py, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   py =  sqrt(3.d0/fpi) * sint * sin(phi)

end function py
!======== l = 2 =====================================================================
function dz2(cost)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dz2, cost
   dz2 =  sqrt(1.25d0/fpi) * (3.d0* cost*cost-1.d0)

end function dz2
function dxz(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dxz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxz =  sqrt(15.d0/fpi) * sint*cost * cos(phi)

end function dxz
function dyz(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dyz =  sqrt(15.d0/fpi) * sint*cost * sin(phi)

end function dyz
function dx2my2(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dx2my2 =  sqrt(3.75d0/fpi) * sint*sint * cos(2.d0*phi)

end function dx2my2
function dxy(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::dxy, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   dxy =  sqrt(3.75d0/fpi) * sint*sint * sin(2.d0*phi)

end function dxy
!======== l = 3 =====================================================================
function fz3(cost)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fz3, cost
   fz3 =  0.25d0*sqrt(7.d0/pi) * ( 5.d0 * cost * cost - 3.d0 ) * cost

end function fz3
function fxz2(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fxz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * cos(phi)

end function fxz2
function fyz2(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fyz2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fyz2 =  0.25d0*sqrt(10.5d0/pi) * ( 5.d0 * cost * cost - 1.d0 ) * sint * sin(phi)

end function fyz2
function fzx2my2(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fzx2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fzx2my2 =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * cos(2.d0*phi)

end function fzx2my2
function fxyz(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fxyz, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxyz =  0.25d0*sqrt(105d0/pi) * sint * sint * cost * sin(2.d0*phi)

end function fxyz
function fxx2m3y2(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fxx2m3y2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fxx2m3y2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * cos(3.d0*phi)

end function fxx2m3y2
function fy3x2my2(cost,phi)
   use kinds, only :  DP
   implicit none
   real(DP), parameter :: pi=3.14159265358979d0, fpi =4.d0*pi
   real(DP) ::fy3x2my2, cost, phi, sint
   sint = sqrt(abs(1.d0 - cost*cost))
   fy3x2my2 =  0.25d0*sqrt(17.5d0/pi) * sint * sint * sint * sin(3.d0*phi)

end function fy3x2my2
!
!
!-----------------------------------------------------------------------
subroutine radialpart(ng, q, alfa, rvalue, lmax, radial)
  !-----------------------------------------------------------------------
  !
  ! This routine computes a table with the radial Fourier transform 
  ! of the radial functions.
  !
  use kinds,      only : dp
  use constants,  only : fpi
  use cell_base,  only : omega
  !
  implicit none
  ! I/O
  integer :: ng, rvalue, lmax
  real(DP) :: q(ng), alfa, radial(ng,0:lmax)
  ! local variables
  real(DP), parameter :: xmin=-6.d0, dx=0.025d0, rmax=10.d0

  real(DP) :: rad_int, pref, x
  integer :: l, ir, ig, mesh_r
  real(DP), allocatable :: bes(:), func_r(:), r(:), rij(:), aux(:)

  mesh_r = nint ( ( log ( rmax ) - xmin ) / dx + 1 )
  allocate ( bes(mesh_r), func_r(mesh_r), r(mesh_r), rij(mesh_r) )
  allocate ( aux(mesh_r))
  !
  !    compute the radial mesh
  !
  DO ir = 1, mesh_r
     x = xmin  + DBLE (ir - 1) * dx 
     r (ir) = exp (x) / alfa
     rij (ir) = dx  * r (ir)
  ENDDO
  !
  IF (rvalue==1) func_r(:) = 2.d0 * alfa**(3.d0/2.d0) * exp(-alfa*r(:))
  IF (rvalue==2) func_r(:) = 1.d0/sqrt(8.d0) * alfa**(3.d0/2.d0) * & 
                     (2.0d0 - alfa*r(:)) * exp(-alfa*r(:)*0.5d0)
  IF (rvalue==3) func_r(:) = sqrt(4.d0/27.d0) * alfa**(2.0d0/3.0d0) * &
                     (1.d0 - 1.5d0*alfa*r(:) + 2.d0*(alfa*r(:))**2/27.d0) * &
                                           exp(-alfa*r(:)/3.0d0)
  pref = fpi/sqrt(omega)
  !
  DO l = 0, lmax
     DO ig=1,ng
       CALL sph_bes (mesh_r, r(1), q(ig), l, bes)
       aux(:) = bes(:) * func_r(:) * r(:)
       CALL simpson (mesh_r, aux, rij, rad_int)
       radial(ig,l) = rad_int * pref
     ENDDO
  ENDDO

  deallocate (bes, func_r, r, rij, aux )

end subroutine radialpart

