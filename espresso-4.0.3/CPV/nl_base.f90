!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!
!-----------------------------------------------------------------------
   subroutine nlsm1 ( n, nspmn, nspmx, eigr, c, becp )
!-----------------------------------------------------------------------

      !     computes: the array becp
      !     becp(ia,n,iv,is)=
      !         = sum_g [(-i)**l beta(g,iv,is) e^(-ig.r_ia)]^* c(g,n)
      !         = delta_l0 beta(g=0,iv,is) c(g=0,n)
      !          +sum_g> beta(g,iv,is) 2 re[(i)**l e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of c*(g)=c(-g)  (g> see routine ggen)
      !     input : beta(ig,l,is), eigr, c
      !     output: becp as parameter
      !
      USE kinds,      ONLY : DP
      USE mp,         ONLY : mp_sum
      USE mp_global,  ONLY : nproc_image, intra_image_comm
      USE ions_base,  only : na, nax, nat
      USE gvecw,      only : ngw
      USE uspp,       only : nkb, nhtol, beta
      USE cvan,       only : ish
      USE uspp_param, only : nh
      !
      USE reciprocal_vectors, ONLY : gstart
!
      implicit none

      integer,   intent(in)  :: n, nspmn, nspmx
      real(DP), intent(in)  :: eigr( 2, ngw, nat ), c( 2, ngw, n )
      real(DP), intent(out) :: becp( nkb, n )
      !
      integer   :: isa, ig, is, iv, ia, l, ixr, ixi, inl, i, nhx
      real(DP)  :: signre, signim, arg
      real(DP), allocatable :: becps( :, : )
      real(DP), allocatable :: wrk2( :, :, : )
      !
      call start_clock( 'nlsm1' )

      allocate( wrk2( 2, ngw, nax ) )

      isa = 0
      do is = 1, nspmn - 1
        isa = isa + na(is)
      end do

      do is = nspmn, nspmx
         !
         IF( nh( is ) < 1 ) THEN
            isa = isa + na(is)
            CYCLE
         END IF
         !
         IF( nproc_image > 1 ) THEN
            nhx = nh( is ) * na( is )
            IF( MOD( nhx, 2 ) /= 0 ) nhx = nhx + 1
            ALLOCATE( becps( nhx, n ) )
            becps = 0.0d0
         END IF
         !
         do iv = 1, nh( is )
            !
            l = nhtol( iv, is )
            !
            if (l == 0) then
               ixr = 1
               ixi = 2
               signre =  1.0d0
               signim =  1.0d0
            else if (l == 1) then
               ixr = 2
               ixi = 1
               signre =  1.0d0
               signim = -1.0d0
            else if (l == 2) then
               ixr = 1
               ixi = 2
               signre = -1.0d0
               signim = -1.0d0
            else if (l == 3) then
               ixr = 2
               ixi = 1
               signre = -1.0d0
               signim =  1.0d0
            endif
!
            do ia=1,na(is)
               !
               !  q = 0   component (with weight 1.0)
               !
               if (gstart == 2) then
                  wrk2( 1, 1, ia ) = signre*beta(1,iv,is)*eigr(ixr,1,ia+isa)
                  wrk2( 2, 1, ia ) = signim*beta(1,iv,is)*eigr(ixi,1,ia+isa)
               end if
               !
               !   q > 0   components (with weight 2.0)
               !
               do ig = gstart, ngw
                  arg = 2.0d0 * beta(ig,iv,is)
                  wrk2( 1, ig, ia ) = signre*arg*eigr(ixr,ig,ia+isa)
                  wrk2( 2, ig, ia ) = signim*arg*eigr(ixi,ig,ia+isa)
               end do
               !
            end do
            
            !
            IF( nproc_image > 1 ) THEN
               inl=(iv-1)*na(is)+1
               CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becps( inl, 1 ), nhx )
            ELSE
               inl=ish(is)+(iv-1)*na(is)+1
               CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becp( inl, 1 ), nkb )
            END IF

         end do


         IF( nproc_image > 1 ) THEN
            !
            inl = ish(is) + 1
            !
            CALL mp_sum( becps, intra_image_comm )

            do i = 1, n
               do iv = inl , ( inl + na(is) * nh(is) - 1 )
                  becp( iv, i ) = becps( iv - inl + 1, i )
               end do
            end do

            DEALLOCATE( becps )

         END IF

         isa = isa + na(is)

      end do

      deallocate( wrk2 )
      call stop_clock( 'nlsm1' )

      return
   end subroutine nlsm1
!-----------------------------------------------------------------------


!-------------------------------------------------------------------------
   subroutine nlsm2( ngw, nkb, n, eigr, c, becdr, tred )
!-----------------------------------------------------------------------

      !     computes: the array becdr
      !     becdr(ia,n,iv,is,k)
      !      =2.0 sum_g> g_k beta(g,iv,is) re[ (i)**(l+1) e^(ig.r_ia) c(g,n)]
      !
      !     routine makes use of  c*(g)=c(-g)  (g> see routine ggen)
      !     input : eigr, c
      !     output: becdr
      !
 
      USE kinds,      ONLY : DP
      use ions_base,  only : nax, nsp, na, nat
      use uspp,       only : nhtol, beta  !, nkb
      use cvan,       only : ish
      use uspp_param, only : nh
      use cell_base,  only : tpiba
      use mp,         only : mp_sum
      use mp_global,  only : nproc_image, intra_image_comm

      use reciprocal_vectors, only: gx, gstart
!
      implicit none
    
      integer,   intent(in)  :: ngw, nkb, n
      real(DP), intent(in)  :: eigr(2,ngw,nat), c(2,ngw,n)
      real(DP), intent(out) :: becdr(nkb,n,3)
      logical,   intent(in)  :: tred
      !
      real(DP), allocatable :: gk(:)
      real(DP), allocatable :: wrk2(:,:,:)
      !
      integer   :: ig, is, iv, ia, k, l, ixr, ixi, inl, isa, i
      real(DP) :: signre, signim, arg
!
      call start_clock( 'nlsm2' )

      allocate( gk( ngw ) )
      allocate( wrk2( 2, ngw, nax ) )
      becdr = 0.d0
!
      do k = 1, 3

         do ig=1,ngw
            gk(ig)=gx(k,ig)*tpiba
         end do
!
         isa = 0

         do is=1,nsp

            do iv=1,nh(is)
               !
               !     order of states:  s_1  p_x1  p_z1  p_y1  s_2  p_x2  p_z2  p_y2
               !
               l=nhtol(iv,is)
               if (l.eq.0) then
                  ixr = 2
                  ixi = 1
                  signre =  1.0d0
                  signim = -1.0d0
               else if (l.eq.1) then
                  ixr = 1
                  ixi = 2
                  signre = -1.0d0
                  signim = -1.0d0
               else if (l.eq.2) then
                  ixr = 2
                  ixi = 1
                  signre = -1.0d0
                  signim =  1.0d0
               else if (l == 3) then
                  ixr = 1
                  ixi = 2
                  signre =  1.0d0
                  signim =  1.0d0
               endif
!    
               do ia=1,na(is)
                  !    q = 0   component (with weight 1.0)
                  if (gstart == 2) then
                     wrk2(1,1,ia) = signre*gk(1)*beta(1,iv,is)*eigr(ixr,1,ia+isa)
                     wrk2(2,1,ia) = signim*gk(1)*beta(1,iv,is)*eigr(ixi,1,ia+isa)
                  end if
                  !    q > 0   components (with weight 2.0)
                  do ig=gstart,ngw
                     arg = 2.0d0*gk(ig)*beta(ig,iv,is)
                     wrk2(1,ig,ia) = signre*arg*eigr(ixr,ig,ia+isa)
                     wrk2(2,ig,ia) = signim*arg*eigr(ixi,ig,ia+isa)
                  end do
               end do
               inl=ish(is)+(iv-1)*na(is)+1
               CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, becdr( inl, 1, k ), nkb )
            end do

            isa = isa + na(is)

         end do

         IF( tred .AND. ( nproc_image > 1 ) ) THEN
            CALL mp_sum( becdr(:,:,k), intra_image_comm )
         END IF
      end do

      deallocate( gk )
      deallocate( wrk2 )

      call stop_clock( 'nlsm2' )
!
      return
   end subroutine nlsm2
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
   real(8) function ennl( rhovan, bec )
!-----------------------------------------------------------------------
      !
      ! calculation of nonlocal potential energy term
      !
      use kinds,          only : DP
      use cvan,           only : ish
      use uspp_param,     only : nhm, nh
      use uspp,           only : nkb, dvan
      use electrons_base, only : n => nbsp, nspin, ispin, f
      use ions_base,      only : nsp, nat, na
      !
      implicit none
      !
      ! input
      !
      real(DP) :: bec( nkb, n )
      real(DP) :: rhovan( nhm*(nhm+1)/2, nat, nspin )
      !
      ! local
      !
      real(DP) :: sum, sums(2)
      integer   :: is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss, i
      !
      ennl = 0.d0
      !
      do is = 1, nsp
         do iv = 1, nh(is)
            do jv = iv, nh(is)
               ijv = (jv-1)*jv/2 + iv
               isa = 0
               do ism = 1, is - 1
                  isa = isa + na(ism)
               end do
               do ia = 1, na(is)
                  inl = ish(is)+(iv-1)*na(is)+ia
                  jnl = ish(is)+(jv-1)*na(is)+ia
                  isa = isa+1
                  sums = 0.d0
                  do i = 1, n
                     iss = ispin(i)
                     sums(iss) = sums(iss) + f(i) * bec(inl,i) * bec(jnl,i)
                  end do
                  sum = 0.d0
                  do iss = 1, nspin
                     rhovan( ijv, isa, iss ) = sums( iss )
                     sum = sum + sums( iss )
                  end do
                  if( iv .ne. jv ) sum = 2.d0 * sum
                  ennl = ennl + sum * dvan( jv, iv, is)
               end do
            end do
         end do
      end do
      !
      return
   end function ennl
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
   subroutine force_nl( fion, bec, becdr )
!-----------------------------------------------------------------------
      !
      !     contribution to fion due to nonlocal part
      !
      use kinds,          only : DP
      use uspp,           only : nkb, dvan, deeq
      use uspp_param,     only : nhm, nh
      use cvan,           only : ish, nvb
      use ions_base,      only : nax, nat, nsp, na
      use electrons_base, only : n => nbsp, ispin, f
      use gvecw,          only : ngw
!
      implicit none

      real(DP), intent(in)  :: bec( nkb, n ), becdr( nkb, n, 3 )
      real(DP), intent(out) :: fion( 3, nat )
!
      integer   :: k, is, ia, isa, iss, inl, iv, jv, i
      real(DP) :: temp
      real(DP) :: tmpbec(nhm,n), tmpdr(nhm,n) ! automatic arrays

      do k = 1, 3
         !
         isa=0
         !
         do is = 1, nsp
            !
            do ia = 1, na(is)
               !
               isa = isa + 1
!
               tmpbec = 0.d0
               tmpdr  = 0.d0
!
               do iv=1,nh(is)
                  do jv=1,nh(is)
                     inl=ish(is)+(jv-1)*na(is)+ia
                     do i=1,n
                        iss=ispin(i)
                        temp=dvan(iv,jv,is)+deeq(jv,iv,isa,iss)
                        tmpbec(iv,i)=tmpbec(iv,i)+temp*bec(inl,i)
                     end do
                  end do
               end do
! 
               do iv=1,nh(is)
                  inl=ish(is)+(iv-1)*na(is)+ia
                  do i=1,n
                     tmpdr(iv,i) = f(i) * becdr(inl,i,k)
                  end do
               end do
!
               do i=1,n
                  do iv=1,nh(is)
                     tmpdr(iv,i)=tmpdr(iv,i)*tmpbec(iv,i)
                  end do
               end do
!
               fion(k,isa) = fion(k,isa) - 2.0d0 * SUM( tmpdr )
!
            end do
         end do
      end do
      !
      return
      !
   end subroutine force_nl
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
   subroutine calbec (nspmn,nspmx,eigr,c,bec)
!-----------------------------------------------------------------------

      !     this routine calculates array bec
      !
      !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
      !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
      !
      !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
      !
      
      USE kinds,          ONLY : DP
      use ions_base,      only : na, nat
      use io_global,      only : stdout
      use cvan,           only : ish
      use electrons_base, only : n => nbsp
      use gvecw,          only : ngw
      use control_flags,  only : iprint, iprsta
      use uspp_param,     only : nh
      use uspp,           only : nkb
!
      implicit none
      !
      integer,      intent(in)  :: nspmn, nspmx
      real(DP),    intent(out) :: bec( nkb, n )
      complex(DP), intent(in)  :: c( ngw, n ), eigr( ngw,nat )

      ! local variables

      integer :: is, ia, i , iv
!
!
      call start_clock( 'calbec' )
      call nlsm1( n, nspmn, nspmx, eigr, c, bec )
!
      if ( iprsta > 2 ) then
         WRITE( stdout,*)
         do is=1,nspmx
            if(nspmx.gt.1) then
               WRITE( stdout,'(33x,a,i4)') ' calbec: bec (is)',is
               WRITE( stdout,'(8f9.4)')                                       &
     &              ((bec(ish(is)+(iv-1)*na(is)+1,i),iv=1,nh(is)),i=1,n)
            else
               do ia=1,na(is)
                  WRITE( stdout,'(33x,a,i4)') ' calbec: bec (ia)',ia
                  WRITE( stdout,'(8f9.4)')                                    &
     &             ((bec(ish(is)+(iv-1)*na(is)+ia,i),iv=1,nh(is)),i=1,n)
               end do
            end if
         end do
      endif
      call stop_clock( 'calbec' )
!
      return
   end subroutine calbec
!-----------------------------------------------------------------------


!-----------------------------------------------------------------------
SUBROUTINE caldbec( ngw, nkb, n, nspmn, nspmx, eigr, c, dbec, tred )
  !-----------------------------------------------------------------------
  !
  !     this routine calculates array dbec, derivative of bec:
  !
  !        < psi_n | beta_i,i > = c_n(0) beta_i,i(0) +
  !                 2 sum_g> re(c_n*(g) (-i)**l beta_i,i(g) e^-ig.r_i)
  !
  !     with respect to cell parameters h
  !
  !     routine makes use of c(-g)=c*(g)  and  beta(-g)=beta*(g)
  !
  USE kinds,      ONLY : DP
  use mp,         only : mp_sum
  use mp_global,  only : nproc_image, intra_image_comm
  use ions_base,  only : na, nax, nat
  use cvan,       only : ish
  use cdvan,      only : dbeta
  use uspp,       only : nhtol
  use uspp_param, only : nh, nhm
  use reciprocal_vectors, only : gstart
  !
  implicit none
  !
  integer,      intent(in)  :: ngw, nkb, n
  integer,      intent(in)  :: nspmn, nspmx
  complex(DP), intent(in)  :: c(ngw,n)
  real(DP),    intent(in)  :: eigr(2,ngw,nat)
  real(DP),    intent(out) :: dbec( nkb, n, 3, 3 )
  logical,      intent(in)  :: tred
  !
  real(DP), allocatable :: wrk2(:,:,:)
  !
  integer   :: ig, is, iv, ia, l, ixr, ixi, inl, i, j, ii, isa
  real(DP) :: signre, signim, arg
  !
  allocate( wrk2( 2, ngw, nax ) )
  !
  !
  do j=1,3
     do i=1,3

        isa = 0
        do is = 1, nspmn - 1
          isa = isa + na(is)
        end do

        do is=nspmn,nspmx
           do iv=1,nh(is)
              l=nhtol(iv,is)
              if (l == 0) then
                 ixr = 1
                 ixi = 2
                 signre =  1.0d0
                 signim =  1.0d0
              else if (l == 1) then
                 ixr = 2
                 ixi = 1
                 signre =  1.0d0
                 signim = -1.0d0
              else if (l == 2) then
                 ixr = 1
                 ixi = 2
                 signre = -1.0d0
                 signim = -1.0d0
              else if (l == 3) then
                 ixr = 2
                 ixi = 1
                 signre = -1.0d0
                 signim =  1.0d0
              else
                 CALL errore(' caldbec  ', ' l not implemented ', ABS( l ) )
              endif
              !
              do ia=1,na(is)
                 if (gstart == 2) then
                    !     q = 0   component (with weight 1.0)
                    wrk2(1,1,ia)= signre*dbeta(1,iv,is,i,j)*eigr(ixr,1,ia+isa)
                    wrk2(2,1,ia)= signim*dbeta(1,iv,is,i,j)*eigr(ixi,1,ia+isa)
                 end if
                 !     q > 0   components (with weight 2.0)
                 do ig = gstart, ngw
                    arg = 2.0d0*dbeta(ig,iv,is,i,j)
                    wrk2(1,ig,ia) = signre*arg*eigr(ixr,ig,ia+isa)
                    wrk2(2,ig,ia) = signim*arg*eigr(ixi,ig,ia+isa)
                 end do
              end do
              inl=ish(is)+(iv-1)*na(is)+1
              CALL DGEMM( 'T', 'N', na(is), n, 2*ngw, 1.0d0, wrk2, 2*ngw, c, 2*ngw, 0.0d0, dbec( inl, 1, i, j ), nkb )
           end do
           if( ( nproc_image > 1 ) .AND. tred ) then
              inl=ish(is)+1
              do ii=1,n
                 call mp_sum( dbec( inl : (inl + na(is)*nh(is) - 1), ii,i,j), intra_image_comm )
              end do
           end if
           isa = isa + na(is)
        end do
     end do
  end do

  deallocate( wrk2 )
  !
  return
end subroutine caldbec
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
subroutine dennl( bec, denl )
  !-----------------------------------------------------------------------
  !
  !  compute the contribution of the non local part of the
  !  pseudopotentials to the derivative of E with respect to h
  !
  USE kinds,      ONLY : DP
  use cvan,       only : ish
  use uspp_param, only : nh
  use uspp,       only : nkb, dvan, deeq
  use cdvan,      ONLY : drhovan, dbec
  use ions_base,  only : nsp, na
  use cell_base,  only : h
  use io_global,  only : stdout
  !
  use electrons_base,     only : n => nbsp, ispin, f, nspin
  use reciprocal_vectors, only : gstart

  implicit none

  real(DP), intent(in)  :: bec( nkb, n )
  real(DP), intent(out) :: denl( 3, 3 )

  real(DP) :: dsum(3,3),dsums(2,3,3), detmp(3,3)
  integer   :: is, iv, jv, ijv, inl, jnl, isa, ism, ia, iss, i,j,k
  !
  denl=0.d0
  do is=1,nsp
     do iv=1,nh(is)
        do jv=iv,nh(is)
           ijv = (jv-1)*jv/2 + iv
           isa=0
           do ism=1,is-1
              isa=isa+na(ism)
           end do
           do ia=1,na(is)
              inl=ish(is)+(iv-1)*na(is)+ia
              jnl=ish(is)+(jv-1)*na(is)+ia
              isa=isa+1
              dsums=0.d0
              do i=1,n
                 iss=ispin(i)
                 do k=1,3
                    do j=1,3
                       dsums(iss,k,j)=dsums(iss,k,j)+f(i)*       &
 &                          (dbec(inl,i,k,j)*bec(jnl,i)          &
 &                          + bec(inl,i)*dbec(jnl,i,k,j))
                    enddo
                 enddo
              end do
              !
              do iss=1,nspin
                 dsum=0.d0
                 do k=1,3
                    do j=1,3
                       drhovan(ijv,isa,iss,j,k)=dsums(iss,j,k)
                       dsum(j,k)=dsum(j,k)+dsums(iss,j,k)
                    enddo
                 enddo
                 if(iv.ne.jv) dsum=2.d0*dsum
                 denl = denl + dsum * dvan(jv,iv,is)
              end do
           end do
        end do
     end do
  end do

!  WRITE(6,*) 'DEBUG enl (CP) = '
!  detmp = denl
!  detmp = MATMUL( detmp(:,:), TRANSPOSE( h ) )
!  WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

  !
  return
end subroutine dennl
!-----------------------------------------------------------------------



!-----------------------------------------------------------------------
subroutine nlfq( c, eigr, bec, becdr, fion )
  !-----------------------------------------------------------------------
  !
  !     contribution to fion due to nonlocal part
  !
  USE kinds,          ONLY : DP
  use uspp,           only : nkb, dvan, deeq
  use uspp_param,     only : nhm, nh
  use cvan,           only : ish, nvb
  use ions_base,      only : nax, nat, nsp, na
  use electrons_base, only : n => nbsp, ispin, f, nspin, iupdwn, nupdwn
  use gvecw,          only : ngw
  use constants,      only : pi, fpi
  use mp_global,      only : me_image, intra_image_comm, nproc_image
  use mp,             only : mp_sum
  USE cp_main_variables, ONLY: nlax, descla, la_proc
  USE descriptors,       ONLY: nlar_ , nlac_ , ilar_ , ilac_ , lambda_node_ , &
                               la_myr_ , la_myc_
  !
  implicit none
  !
  real(DP),    intent(in)  :: bec( nkb, n ), c( 2, ngw, n )
  real(DP),    intent(out) :: becdr( nkb, n, 3 )
  complex(DP), intent(in)  :: eigr( ngw, nat )
  real(DP),    intent(out) :: fion( 3, nat )
  !
  integer   :: k, is, ia, isa, iss, inl, iv, jv, i, ir, nr, nss, istart, ioff
  real(DP) :: temp
  !
  real(DP), allocatable :: tmpbec(:,:), tmpdr(:,:) 
  real(DP), allocatable :: fion_loc(:,:)
  !
  call start_clock( 'nlfq' )
  !
  !
  !     nlsm2 fills becdr
  !
  call nlsm2( ngw, nkb, n, eigr, c, becdr, .true. )
  !
  allocate ( fion_loc( 3, nat ) )
  !
  fion_loc = 0.0d0

  allocate ( tmpbec( nhm, nlax ), tmpdr( nhm, nlax ) )
  !
  DO k=1,3

     isa=0
     DO is=1,nsp
        DO ia=1,na(is)

           isa=isa+1

           DO iss = 1, nspin

              nss = nupdwn( iss )
              istart = iupdwn( iss )

              IF( la_proc .AND. &
                  ( descla( la_myr_ , iss ) == descla( la_myc_ , iss ) ) ) THEN

                 ! only processors on the diagonal of the square proc grid enter here.
                 ! This is to distribute the load among different multi-core nodes,
                 ! and maximize the memory bandwith per core.

                 tmpbec = 0.d0
                 tmpdr  = 0.d0

                 ir = descla( ilar_ , iss )
                 nr = descla( nlar_ , iss )

                 ioff = istart-1+ir-1

                 do iv=1,nh(is)
                    do jv=1,nh(is)
                       inl=ish(is)+(jv-1)*na(is)+ia
                       temp=dvan(iv,jv,is)+deeq(jv,iv,isa,iss)
                       do i=1,nr
                          tmpbec(iv,i)=tmpbec(iv,i)+temp*bec(inl,i+ioff)
                       end do
                    end do
                 end do

                 do iv=1,nh(is)
                    inl=ish(is)+(iv-1)*na(is)+ia
                    do i=1,nr
                       tmpdr(iv,i)=f(i+ioff)*becdr(inl,i+ioff,k)
                    end do
                 end do

                 do i=1,nr
                    do iv=1,nh(is)
                       tmpdr(iv,i)=tmpdr(iv,i)*tmpbec(iv,i)
                    end do
                 end do

                 fion_loc(k,isa) = fion_loc(k,isa)-2.d0*SUM(tmpdr)

              END IF
           END DO
        END DO
     END DO
  END DO
  !
  CALL mp_sum( fion_loc, intra_image_comm )
  !
  fion = fion + fion_loc
  !
  !     end of x/y/z loop
  !
  deallocate ( fion_loc )
  !
  deallocate ( tmpbec, tmpdr )

  call stop_clock( 'nlfq' )
  !
  return
end subroutine nlfq

