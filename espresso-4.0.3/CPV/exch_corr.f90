!
! Copyright (C) 2002-2005 FPMD-CPV groups
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!

#include "f_defs.h"


! ... Gradient Correction & exchange and correlation

!=----------------------------------------------------------------------------=!
        SUBROUTINE v2gc_x( v2xc, grho, rhor, vpot )
!=----------------------------------------------------------------------------=!

          USE kinds,              ONLY: DP
          USE fft_base,           ONLY: dfftp
          USE cell_base,          ONLY: tpiba
          USE reciprocal_vectors, ONLY: gstart, gx
          use grid_dimensions,    only: nnrx
          USE gvecp,              ONLY: ngm
          USE cp_interfaces,      ONLY: fwfft, invfft
!                                                                       
          implicit none
!                                                                       
          REAL(DP) ::  vpot(:,:)
          REAL(DP), intent(in)  ::  v2xc(:,:,:)
          REAL(DP), intent(in)  ::  grho(:,:,:)
          REAL(DP), intent(in)  ::  rhor(:,:)
!                                                                       
          integer :: ig, ipol, is, js, nspin
          COMPLEX(DP), allocatable ::  psi(:)
          COMPLEX(DP), allocatable ::  vtemp(:)
          COMPLEX(DP), allocatable ::  vtemp_pol(:)
          REAL(DP), ALLOCATABLE :: v(:)
          REAL(DP) :: fac
! ...                                                                   
          nspin = SIZE(rhor,2)

          fac = 1.0d0

          ALLOCATE( vtemp( ngm ) )
          ALLOCATE( vtemp_pol( ngm ) )
          ALLOCATE( psi( nnrx ) )

          DO js = 1, nspin
            !
            vtemp = 0.0d0

            DO ipol = 1, 3
              DO is = 1, nspin
                !
                psi( 1:nnrx ) = fac * v2xc( 1:nnrx, js, is ) * grho( 1:nnrx, ipol, is )
                !
                CALL fwfft(   'Dense', psi, dfftp )
                CALL psi2rho( 'Dense', psi, dfftp%nnr, vtemp_pol, ngm )
                !
                DO ig = gstart, ngm
                  vtemp(ig) = vtemp(ig) + vtemp_pol(ig) *  CMPLX( 0.d0, tpiba * gx( ipol, ig ) )
                END DO
                !
              END DO
            END DO
            !
            CALL rho2psi( 'Dense', psi, dfftp%nnr, vtemp, ngm )
            CALL invfft(  'Dense', psi, dfftp )

            vpot( 1:nnrx, js ) = vpot( 1:nnrx, js) - DBLE( psi( 1:nnrx ) )

          END DO

          DEALLOCATE( psi )
          DEALLOCATE( vtemp_pol )
          DEALLOCATE( vtemp )

          RETURN
        END SUBROUTINE v2gc_x


!=----------------------------------------------------------------------------=!
    SUBROUTINE stress_gc_x(grho, v2xc, gcpail, omega)
!=----------------------------------------------------------------------------=!
!
        USE kinds,           ONLY: DP
        use grid_dimensions, only: nr1, nr2, nr3, nnrx

        IMPLICIT NONE
!
        REAL(DP) ::  v2xc(:,:,:)
        REAL(DP) ::  grho(:,:,:)
        REAL(DP) ::  gcpail(6)
        REAL(DP) ::  omega
!
        REAL(DP) :: stre, grhoi, grhoj
        INTEGER :: i, ipol, jpol, ic, is, js, nspin
        INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
        INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)
! ...
        nspin = SIZE(grho,3)

        DO ic = 1, 6
          ipol = alpha(ic)
          jpol = beta(ic)
          stre = 0.0d0
          DO is = 1, nspin
            DO js = 1, nspin
              DO i = 1, nnrx
                stre = stre + v2xc(i,is,js) * grho(i,ipol,js) * grho(i,jpol,is)
              END DO
            END DO
          END DO
          gcpail(ic) = - stre * omega / DBLE(nr1*nr2*nr3)
        END DO

      RETURN
    END SUBROUTINE stress_gc_x


!=----------------------------------------------------------------------------=!
    SUBROUTINE stress_xc_x( dexc, strvxc, sfac, vxc, grho, v2xc, &
        gagb, tnlcc, rhocp, box)
!=----------------------------------------------------------------------------=!

      USE kinds,              ONLY: DP
      USE ions_base,          ONLY: nsp
      USE cell_base,          ONLY: tpiba, boxdimensions
      USE funct,              ONLY: dft_is_gradient
      USE reciprocal_vectors, ONLY: gstart, g
      USE gvecp,              ONLY: ngm
      USE io_global,          ONLY: stdout
      USE cp_interfaces,      ONLY: stress_gc

      IMPLICIT NONE

      ! -- ARGUMENT

      type (boxdimensions), intent(in) :: box
      LOGICAL :: tnlcc(:)
      COMPLEX(DP) :: vxc(:,:)
      COMPLEX(DP), INTENT(IN) :: sfac(:,:)
      REAL(DP) :: dexc(:), strvxc
      REAL(DP) :: grho(:,:,:)
      REAL(DP) :: v2xc(:,:,:)
      REAL(DP) :: gagb(:,:)
      REAL(DP) :: rhocp(:,:)

      INTEGER, DIMENSION(6), PARAMETER :: alpha = (/ 1,2,3,2,3,3 /)
      INTEGER, DIMENSION(6), PARAMETER :: beta  = (/ 1,1,1,2,2,3 /)
      ! ...  dalbe(:) = delta(alpha(:),beta(:))
      REAL(DP),  DIMENSION(6), PARAMETER :: dalbe = &
         (/ 1.0_DP, 0.0_DP, 0.0_DP, 1.0_DP, 0.0_DP, 1.0_DP /)

      COMPLEX(DP) :: tex1, tex2, tex3
      REAL(DP) :: gcpail(6), omega, detmp( 3, 3 )
      REAL(DP) :: dcc( 6 )
      INTEGER :: ig, k, is, ispin, nspin
      INTEGER :: i, j

      omega = box%deth
      nspin = SIZE(vxc, 2)

      DEXC = 0.0d0
      dcc  = 0.0d0

      ! ... computes omega * \sum_{G}[ S(G)*rhopr(G)* G_{alpha} G_{beta}/|G|]
      ! ... (252) Phd thesis Dal Corso. Opposite sign.

      IF ( ANY( tnlcc ) ) THEN

        DO ig = gstart, ngm
          tex1 = ( 0.0d0 , 0.0d0 )
          DO is=1,nsp
            IF ( tnlcc(is) ) THEN
              tex1 = tex1 + sfac( ig, is ) * CMPLX(rhocp(ig,is), 0.d0)
            END IF
          END DO
          tex2 = 0.0d0
          DO ispin = 1, nspin
            tex2 = tex2 + CONJG( vxc(ig, ispin) )
          END DO
          tex3 = DBLE(tex1 * tex2) / SQRT( g( ig ) ) / tpiba
          dcc = dcc + tex3 * gagb(:,ig)
        END DO
        dcc = dcc * 2.0d0 * omega

        ! DEBUG
        !DO k=1,6
        !  detmp(alpha(k),beta(k)) = dcc(k)
        !  detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
        !END DO
        !detmp = MATMUL( detmp(:,:), box%m1(:,:) )
        !WRITE( stdout,*) "derivative of e(xc) - nlcc part"
        !WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
        !WRITE( stdout,*) SUM( rhocp(:,1) )


      END IF

      ! ... (E_{xc} - \int dr v_{xc}(n) n(r))/omega part of the stress
      ! ... this part of the stress is diagonal.

      dexc = strvxc * dalbe

      IF ( dft_is_gradient() ) THEN
        CALL stress_gc(grho, v2xc, gcpail, omega)
        dexc = dexc + gcpail
      END IF

      ! DEBUG
      !DO k=1,6
      !   detmp(alpha(k),beta(k)) = dexc(k)
      !   detmp(beta(k),alpha(k)) = detmp(alpha(k),beta(k))
      !END DO
      !detmp = MATMUL( detmp(:,:), box%m1(:,:) )
      !WRITE( stdout,*) "derivative of e(xc)"
      !WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)

      dexc = dexc + dcc

      RETURN

5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)

      END SUBROUTINE stress_xc_x




!=----------------------------------------------------------------------------=!
     SUBROUTINE exch_corr_energy_x(rhoetr, grho, vpot, exc, vxc, v2xc)
!=----------------------------------------------------------------------------=!

        USE kinds,           ONLY: DP
        use grid_dimensions, ONLY: nnrx
        USE funct,           ONLY: dft_is_gradient
        USE cp_interfaces,   ONLY: v2gc

        implicit none

        REAL (DP) :: rhoetr(:,:)
        REAL (DP) :: grho(:,:,:)
        REAL (DP) :: vpot(:,:)
        REAL (DP) :: exc              ! E_xc   energy ( not multiplied by Omega/nnr )
        REAL (DP) :: vxc              ! SUM ( v(r) * rho(r) )
        REAL (DP) :: v2xc(:,:,:)
        !
        REAL (DP), EXTERNAL :: ddot

        INTEGER :: nspin, ispin
        logical :: is_gradient

        is_gradient =  dft_is_gradient()


        !  vpot = vxc(rhoetr); vpot(r) <-- u(r)

        nspin = SIZE( rhoetr, 2 )
        !
        IF( SIZE( vpot, 1 ) /= nnrx ) &
           CALL errore(" exch_corr_energy ", " inconsistent size for vpot ", 1 )
        !
        CALL exch_corr_wrapper( nnrx, nspin, grho(1,1,1), rhoetr(1,1), exc, vpot(1,1), v2xc(1,1,1) )
        !
        IF( dft_is_gradient() ) THEN
          ! ... vpot additional term for gradient correction
          CALL v2gc( v2xc, grho, rhoetr, vpot )
        END If

        !
        ! vxc = SUM( vpot * rhoetr )
        !
        vxc = 0.0d0
        DO ispin = 1, nspin
           vxc = vxc + DDOT ( nnrx, vpot(1,ispin), 1, rhoetr(1,ispin), 1 )
        END DO


        RETURN
      END SUBROUTINE exch_corr_energy_x




!=----------------------------------------------------------------------------=!
!  CP subroutines
!=----------------------------------------------------------------------------=!

      subroutine exch_corr_h( nspin, rhog, rhor, rhoc, sfac, exc, dxc, self_exc )
!
! calculate exch-corr potential, energy, and derivatives dxc(i,j)
! of e(xc) with respect to to cell parameter h(i,j)
!     
      use funct,           only : dft_is_gradient, dft_is_meta
      use gvecp,           only : ng => ngm
      use gvecs,           only : ngs
      use grid_dimensions, only : nr1, nr2, nr3, nnr => nnrx
      use cell_base,       only : ainv, omega, h
      use ions_base,       only : nsp
      use control_flags,   only : tpre, iprsta
      use derho,           only : drhor
      use core,            only : drhocg, nlcc_any
      use mp,              only : mp_sum
      use metagga,         ONLY : kedtaur
      USE io_global,       ONLY : stdout
      USE mp_global,       ONLY : intra_image_comm
      use kinds,           ONLY : DP
      use constants,       ONLY : au_gpa
      USE sic_module,      ONLY : self_interaction, sic_alpha
      USE cp_interfaces,   ONLY : fillgrad
!
      implicit none

      ! input     
      !
      integer nspin
      !
      ! rhog contains the charge density in G space
      ! rhor contains the charge density in R space
      !
      complex(DP) :: rhog( ng, nspin )
      complex(DP) :: sfac( ngs, nsp )
      !
      ! output
      ! rhor contains the exchange-correlation potential
      !
      real(DP) :: rhor( nnr, nspin ), rhoc( nnr )
      real(DP) :: dxc( 3, 3 ), exc
      real(DP) :: dcc( 3, 3 ), drc( 3, 3 )
      !
      ! local
      !
      integer :: i, j, ir, iss
      real(DP) :: dexc(3,3)
      real(DP), allocatable :: gradr(:,:,:)
      !
      !sic
      REAL(DP) :: self_exc
      REAL(DP), ALLOCATABLE :: self_rho( :,: ), self_gradr(:,:,:)
      complex(DP), ALLOCATABLE :: self_rhog( :,: )
      LOGICAL  :: ttsic
      real(DP) :: detmp(3,3)
      !
      !     filling of gradr with the gradient of rho using fft's
      !
      if ( dft_is_gradient() ) then
         !
         allocate( gradr( nnr, 3, nspin ) )
         call fillgrad( nspin, rhog, gradr )
         ! 
      end if


      ttsic = (self_interaction /= 0 )
      !
      IF ( ttsic ) THEN
         !
         IF ( dft_is_meta() ) CALL errore ('exch_corr_h', &
                               'SIC and metadynamics not together', 1)
         IF ( tpre ) CALL errore( 'exch_corr_h', 'SIC and stress not implemented', 1)

         !  allocate the sic_arrays
         !
         ALLOCATE( self_rho( nnr, nspin ) )
         ALLOCATE( self_rhog( ng, nspin ) )
         IF( dft_is_gradient() ) ALLOCATE( self_gradr( nnr, 3, nspin ) )

         self_rho(:, 1) = rhor( :, 2)
         self_rho(:, 2) = rhor( :, 2)

         IF( dft_is_gradient() ) THEN
            self_gradr(:, :, 1) = gradr(:, :, 2)
            self_gradr(:, :, 2) = gradr(:, :, 2)
         ENDIF

         self_rhog(:, 1) = rhog( :, 2)
         self_rhog(:, 2) = rhog( :, 2)
!
      END IF
!
      self_exc = 0.d0
!
      if( dft_is_meta() ) then
         !
         call tpssmeta( nnr, nspin, gradr, rhor, kedtaur, exc )
         !
      else
         !
         CALL exch_corr_cp(nnr, nspin, gradr, rhor, exc)
         !
         IF ( ttsic ) THEN
            CALL exch_corr_cp(nnr, nspin, self_gradr, self_rho, self_exc)
            self_exc = sic_alpha * (exc - self_exc)
            exc = exc - self_exc
         END IF
!
      end if

      call mp_sum( exc, intra_image_comm )
      IF ( ttsic ) call mp_sum( self_exc, intra_image_comm )

      exc = exc * omega / DBLE( nr1 * nr2 * nr3 )
      IF ( ttsic ) self_exc = self_exc * omega/DBLE(nr1 * nr2 *nr3 )

      !     WRITE(*,*) 'Debug: calcolo exc', exc, 'eself', self_exc
      !
      ! exchange-correlation contribution to pressure
      !
      dxc = 0.0d0
      !
      if ( tpre ) then
         !
         !  Add term: Vxc( r ) * Drhovan( r )_ij - Vxc( r ) * rho( r ) * ((H^-1)^t)_ij
         !
         do iss = 1, nspin
            do j=1,3
               do i=1,3
                  do ir=1,nnr
                     dxc(i,j) = dxc(i,j) + rhor( ir, iss ) * drhor( ir, iss, i, j )
                  end do
               end do
            end do
         end do
         !
         dxc = dxc * omega / ( nr1*nr2*nr3 )
         !
         call mp_sum ( dxc, intra_image_comm )
         !
         do j = 1, 3
            do i = 1, 3
               dxc( i, j ) = dxc( i, j ) + exc * ainv( j, i )
            end do
         end do
         !
         ! DEBUG
         !
         ! write (stdout,*) "derivative of e(xc)"
         ! write (stdout,5555) ((dxc(i,j),j=1,3),i=1,3)
         !
         IF( iprsta >= 2 ) THEN
            DO i=1,3
               DO j=1,3
                  detmp(i,j)=exc*ainv(j,i)
               END DO
            END DO
            WRITE( stdout,*) "derivative of e(xc) - diag - kbar"
            detmp = -1.0d0 * MATMUL( detmp, TRANSPOSE( h ) ) / omega * au_gpa * 10.0d0
            WRITE( stdout,5555) ((detmp(i,j),j=1,3),i=1,3)
         END IF
         !
      end if
      !
      if (dft_is_gradient()) then
         !
         !  Add second part of the xc-potential to rhor
         !  Compute contribution to the stress dexc
         !
         call gradh( nspin, gradr, rhog, rhor, dexc)
         !
         if (tpre) then
            !
            call mp_sum ( dexc, intra_image_comm )
            !
            dxc = dxc + dexc
            !
         end if
         !
      end if
      !

      IF( ttsic ) THEN
!
         IF (dft_is_gradient()) then
      
            call gradh( nspin, self_gradr, self_rhog, self_rho, dexc)
          
            gradr(:,:, 1) = (1.d0 - sic_alpha ) * gradr(:,:, 1)
            gradr(:,:, 2) = (1.d0 - sic_alpha ) * gradr(:,:, 2) + &
                          &  sic_alpha * ( self_gradr(:,:,1) + self_gradr(:,:,2) )
         ENDIF

         rhor(:, 1) = (1.d0 - sic_alpha ) * rhor(:, 1)
         rhor(:, 2) = (1.d0 - sic_alpha ) * rhor(:, 2) + &
                    &  sic_alpha * ( self_rho(:,1) + self_rho(:,2) )

         IF(ALLOCATED(self_gradr)) DEALLOCATE(self_gradr)
         IF(ALLOCATED(self_rhog)) DEALLOCATE(self_rhog)
         IF(ALLOCATED(self_rho)) DEALLOCATE(self_rho)
!
      ENDIF

      IF( tpre ) THEN
         !
         dcc = 0.0d0
         !
         IF( nlcc_any ) CALL  denlcc( nnr, nspin, rhor, sfac, drhocg, dcc )
         !
         ! DEBUG
         !
         !  write (stdout,*) "derivative of e(xc) - nlcc part"
         !  write (stdout,5555) ((dcc(i,j),j=1,3),i=1,3)
         !
         dxc = dxc + dcc
         !
         do iss = 1, nspin
            drc = 0.0d0
            IF( nlcc_any ) THEN
               do j=1,3
                  do i=1,3
                     do ir=1,nnr
                        drc(i,j) = drc(i,j) + rhor( ir, iss ) * rhoc( ir ) * ainv(j,i)
                     end do
                  end do
               end do
               call mp_sum ( drc, intra_image_comm )
            END IF
            dxc = dxc - drc * ( 1.0d0 / nspin ) * omega / ( nr1*nr2*nr3 )
         end do
         !
      END IF
      !

      IF( ALLOCATED( gradr ) ) DEALLOCATE( gradr )

5555  format(1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5/                                &
     &       1x,f12.5,1x,f12.5,1x,f12.5//)
      !
      return
      end subroutine exch_corr_h


!=----------------------------------------------------------------------------=!

      subroutine gradh( nspin, gradr, rhog, rhor, dexc )
!     _________________________________________________________________
!     
!     calculate the second part of gradient corrected xc potential
!     plus the gradient-correction contribution to pressure
!           
      USE kinds,              ONLY: DP
      use control_flags, only: iprint, tpre
      use reciprocal_vectors, only: gx
      use recvecs_indexes, only: np, nm
      use gvecp, only: ng => ngm
      use grid_dimensions, only: nr1, nr2, nr3, nnr => nnrx, nr1x, nr2x, nr3x
      use cell_base, only: ainv, tpiba, omega
      use derho, only: drhog
      USE cp_interfaces, ONLY: fwfft, invfft
      USE fft_base,      ONLY: dfftp
!                 
      implicit none  
! input                   
      integer nspin
      real(DP)    :: gradr( nnr, 3, nspin ), rhor( nnr, nspin ), dexc( 3, 3 )
      complex(DP) :: rhog( ng, nspin )
!
      complex(DP), allocatable:: v(:)
      complex(DP), allocatable:: x(:), vtemp(:)
      complex(DP) ::  ci, fp, fm
      integer :: iss, ig, ir, i,j
!
      allocate(v(nnr))
      allocate(x(ng))
      allocate(vtemp(ng))
      !
      ci=(0.0d0,1.0d0)
      !
      dexc = 0.0d0
      !
      do iss=1, nspin
!     _________________________________________________________________
!     second part xc-potential: 3 forward ffts
!
         do ir=1,nnr
            v(ir)=CMPLX(gradr(ir,1,iss),0.d0)
         end do
         call fwfft('Dense',v, dfftp )
         do ig=1,ng
            x(ig)=ci*tpiba*gx(1,ig)*v(np(ig))
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,ng
                     vtemp(ig) = omega*ci*CONJG(v(np(ig)))*             &
     &                    tpiba*(-rhog(ig,iss)*gx(i,ig)*ainv(j,1)+      &
     &                    gx(1,ig)*drhog(ig,iss,i,j))
                  end do
                  dexc(i,j) = dexc(i,j) + DBLE(SUM(vtemp))*2.0d0
               end do
            end do
         endif
!
         do ir=1,nnr
            v(ir)=CMPLX(gradr(ir,2,iss),gradr(ir,3,iss))
         end do
         call fwfft('Dense',v, dfftp )
!
         do ig=1,ng
            fp=v(np(ig))+v(nm(ig))
            fm=v(np(ig))-v(nm(ig))
            x(ig) = x(ig) +                                             &
     &           ci*tpiba*gx(2,ig)*0.5d0*CMPLX( DBLE(fp),AIMAG(fm))
            x(ig) = x(ig) +                                             &
     &           ci*tpiba*gx(3,ig)*0.5d0*CMPLX(AIMAG(fp),-DBLE(fm))
         end do
!
         if(tpre) then
            do i=1,3
               do j=1,3
                  do ig=1,ng
                     fp=v(np(ig))+v(nm(ig))
                     fm=v(np(ig))-v(nm(ig))
                     vtemp(ig) = omega*ci*                              &
     &                    (0.5d0*CMPLX(DBLE(fp),-AIMAG(fm))*              &
     &                    tpiba*(-rhog(ig,iss)*gx(i,ig)*ainv(j,2)+      &
     &                    gx(2,ig)*drhog(ig,iss,i,j))+                  &
     &                    0.5d0*CMPLX(AIMAG(fp),DBLE(fm))*tpiba*          &
     &                    (-rhog(ig,iss)*gx(i,ig)*ainv(j,3)+            &
     &                    gx(3,ig)*drhog(ig,iss,i,j)))
                  end do
                  dexc(i,j) = dexc(i,j) + 2.0d0*DBLE(SUM(vtemp))
               end do
            end do
         endif
!     _________________________________________________________________
!     second part xc-potential: 1 inverse fft
!
         do ig=1,nnr
            v(ig)=(0.0d0,0.0d0)
         end do
         do ig=1,ng
            v(np(ig))=x(ig)
            v(nm(ig))=CONJG(x(ig))
         end do
         call invfft('Dense',v, dfftp )
         do ir=1,nnr
            rhor(ir,iss)=rhor(ir,iss)-DBLE(v(ir))
         end do
      end do
!
      deallocate(vtemp)
      deallocate(x)
      deallocate(v)
!
      return
      end subroutine gradh

