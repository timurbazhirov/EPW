!
! Copyright (C) 2001-2005 Quantum-ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "f_defs.h"
!-----------------------------------------------------------------------
SUBROUTINE g_tensor_crystal
  !-----------------------------------------------------------------------
  !
  ! ... Compute the g-tensor: PRL 88, 086403 (2002)
  !
  USE kinds,                       ONLY : DP
  USE pwcom,                       ONLY : nrxxs, lmaxx, ap, pi, ecutwfc, &
                                          nl, vkb
  USE io_global,                   ONLY : stdout
  USE io_files,                    ONLY : nwordwfc, iunwfc
  USE cell_base,                   ONLY : at, bg, omega, tpiba, tpiba2
  USE wavefunctions_module,        ONLY : evc
  USE klist,                       ONLY : nks, nkstot, wk, xk, nelec
  USE wvfct,                       ONLY : nbnd, npwx, npw, igk, wg, g2kin, &
                                          current_k
  USE lsda_mod,                    ONLY : current_spin, lsda, isk, nspin
  USE becmod,                      ONLY : becp, calbec
  USE symme,                       ONLY : nsym, s, ftau
  USE scf,                         ONLY : v, vltot, rho
  USE gvect,                       ONLY : ngm, nr1, nr2, nr3, nrx1, nrx2, &
                                          nrx3, nrxx, nlm, g
  USE gipaw_module,                ONLY : j_bare, b_ind, b_ind_r, q_gipaw, &
                                          evq, alpha, nbnd_occ, iverbosity, &
                                          isolve, conv_threshold
  USE buffers,                     ONLY : get_buffer
  USE paw_gipaw,                   ONLY : paw_recon, paw_nkb, paw_vkb, paw_becp
  USE mp_global,                   ONLY : my_pool_id, &
                                          inter_pool_comm, intra_pool_comm
  USE mp,                          ONLY : mp_sum
  USE ions_base,                   ONLY : nat
  
  !-- local variables ----------------------------------------------------
  IMPLICIT NONE
  complex(dp), allocatable, dimension(:,:,:) :: p_evc, vel_evc, g_vel_evc
  complex(dp), allocatable :: aux(:,:), aux1(:)

  ! Q tensor of eq. (65) (pGv => HH in Paratec, vGv => VV in Paratec)
  real(dp) :: q_pGv(3,3,-1:1), q_vGv(3,3,-1:1)

  ! F tensor of eq. (64)
  real(dp) :: f_pGv(3,3,-1:1), f_vGv(3,3,-1:1)

  ! chi_bare tensor of eq. (64)
  real(dp) :: chi_bare_pGv(3,3), chi_bare_vGv(3,3)

  ! f-sum rule
  real(dp) :: f_sum(3,3)

  real(dp) :: q(3), braket, delta_rmc, gipaw_delta_rmc, rmc_gipaw
  integer :: ia, ib, ik, ipol, jpol, i, ibnd, isign, ispin
  complex(dp), external :: ZDOTC

  integer :: s_maj, s_min
  real(dp) :: e_hartree, charge, s_weight, rho_diff, d_omega
  real(dp), allocatable :: grad_vr(:,:), v_local(:,:)
  real(dp), allocatable :: grad_vh(:,:), vh(:,:)
  real(dp), dimension ( 3, 3 ) :: delta_g_rmc, delta_g_bare, delta_g_soo, &
       delta_g_soo_2, delta_g_paramagn, delta_g_diamagn, delta_g_total, &
       delta_g_rmc_gipaw
  real(dp) :: g_e = 2.0023192778_DP, g_prime, units_Ry2Ha = 0.5_DP
  
  logical :: tevaluate_chi, tcalculate_correct_delta_g_soo
  
  real(dp) :: diamagnetic_corr_tensor(3,3)
  real(dp) :: paramagnetic_corr_tensor(3,3)
  real(dp) :: sigma_paramagnetic(3,3)
  
  !-----------------------------------------------------------------------
  
  delta_g_paramagn = 0.0_DP
  delta_g_diamagn = 0.0_DP

  g_prime = 2 * ( g_e - 1 )
  
  !<apsi> TMPTMPTMP Move to input or make it default?
  tcalculate_correct_delta_g_soo = .true.
  !</apsi>
  if ( tcalculate_correct_delta_g_soo ) then
     tevaluate_chi = .true.
  else
     tevaluate_chi = .false.
  end if
  
  ! Select majority and minority spin components
  rho_diff = SUM ( rho%of_r( :, 1 ) - rho%of_r( :, nspin ) )
#ifdef __PARA
  call mp_sum ( rho_diff, intra_pool_comm )
#endif
  if ( rho_diff > +1.0d-3 ) then
     s_maj = 1
     s_min = nspin
  else if ( rho_diff < -1.0d-3 ) then
     s_maj = nspin
     s_min = 1
  else
     write ( stdout, * ) "WARNING: rho_diff zero!"
  end if
  
  ! allocate memory
  allocate ( p_evc(npwx,nbnd,3), vel_evc(npwx,nbnd,3), &
       aux(npwx,nbnd), g_vel_evc(npwx,nbnd,3) )
  
  ! zero the Q tensors
  q_pGv(:,:,:) = 0.0_dp
  q_vGv(:,:,:) = 0.0_dp
  
  ! zero the current and the field
  j_bare(:,:,:,:) = (0.0_dp,0.0_dp)
  b_ind(:,:,:) = (0.0_dp,0.0_dp)
  
  delta_rmc = 0.0_DP
  gipaw_delta_rmc = 0.0_DP
  
  write(stdout, '(5X,''Computing the magnetic susceptibility'',$)')
  write(stdout, '(5X,''isolve='',I1,4X,''ethr='',E10.4)') isolve, conv_threshold
  !====================================================================
  ! loop over k-points
  !====================================================================
  do ik = 1, nks
#ifdef __PARA
    if (my_pool_id == 0) &
    write(*, '(5X,''k-point #'',I5,'' of '',I5,'' (my_pool_id='',I3)') &
      ik, nks, my_pool_id
#else
    write(*, '(5X,''k-point #'',I5,'' of '',I5)') ik, nks
#endif
    current_k = ik
    current_spin = isk(ik)
    
    if ( current_spin == s_maj ) then
       s_weight = +1
    else
       s_weight = -1
    end if
    
    ! initialize at k-point k 
    call gk_sort(xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin)
    g2kin(:) = g2kin(:) * tpiba2
    call init_us_2(npw,igk,xk(1,ik),vkb)
    
    ! read wfcs from file and compute becp
    call get_buffer (evc, nwordwfc, iunwfc, ik)
    
    ! this is the case q = 0 (like the case of the f-sum rule)
    q(:) = 0.0_dp
    !!!write(*,'(''q='',3(F12.4))') q
    call compute_u_kq(ik, q)
    
    !<apsi>
    call init_gipaw_2_no_phase (npw, igk, xk (1, ik), paw_vkb)
    !it was: call ccalbec (paw_nkb, npwx, npw, nbnd, paw_becp, paw_vkb, evc)
    call calbec (npw, paw_vkb, evc, paw_becp)
    diamagnetic_corr_tensor = 0.0d0
    call diamagnetic_correction ( diamagnetic_corr_tensor )
    delta_g_diamagn = delta_g_diamagn + s_weight * diamagnetic_corr_tensor
    !</apsi>
    
    do ibnd = 1, nbnd_occ(ik)
       delta_rmc = delta_rmc - s_weight &
            * wg(ibnd,ik) * SUM ( g2kin(:npw) &
            * conjg(evc(:npw,ibnd)) * evc(:npw,ibnd) )
    end do
    
    CALL relativistic_mass_correction ( rmc_gipaw )
    gipaw_delta_rmc = gipaw_delta_rmc - s_weight * rmc_gipaw
    
    ! compute p_k|evc>, v_k|evc> and G_k v_{k,k}|evc>
    call apply_operators

    !------------------------------------------------------------------
    ! f-sum rule (pGv term only) 
    !------------------------------------------------------------------
    do ia = 1, 3 
      do ib = 1, 3
        do ibnd = 1, nbnd_occ(ik)
          braket = 2.0_dp*real(ZDOTC(npw, p_evc(1,ibnd,ia), 1, &
                                        g_vel_evc(1,ibnd,ib), 1), DP)
          f_sum(ia,ib) = f_sum(ia,ib) + wg(ibnd,ik) * braket
        enddo
      enddo
    enddo

    !------------------------------------------------------------------
    ! pGv and vGv contribution to chi_{bare}
    !------------------------------------------------------------------
    if ( tevaluate_chi ) then
      do i = 1, 3
        call add_to_tensor(q_pGv(:,:,0), p_evc, g_vel_evc)
        call add_to_tensor(q_vGv(:,:,0), vel_evc, g_vel_evc)
      enddo
    end if
    
    !------------------------------------------------------------------
    ! loop over -q and +q
    !------------------------------------------------------------------
    do isign = -1, 1, 2
      
      ! loop over cartesian directions
      do i = 1, 3
        ! set the q vector
        q(:) = 0.0_dp
        q(i) = real(isign,dp) * q_gipaw
        !!!write(*,'(''q='',3(F12.4))') q
        
        ! compute the wfcs at k+q
        call compute_u_kq(ik, q)

        ! compute p_k|evc>, v_k|evc> and G_{k+q} v_{k+q,k}|evc>
        call apply_operators
        
        !<apsi>
        call init_gipaw_2_no_phase (npw, igk, xk(:,ik)+q(:), paw_vkb)
        call paramagnetic_correction ( paramagnetic_corr_tensor )
        paramagnetic_corr_tensor = paramagnetic_corr_tensor * s_weight
        call add_to_sigma_para ( paramagnetic_corr_tensor, delta_g_paramagn )
        !</apsi>
        
        if ( tevaluate_chi ) then
          ! pGv and vGv contribution to chi_bare
          call add_to_tensor(q_pGv(:,:,isign), p_evc, g_vel_evc)
          call add_to_tensor(q_vGv(:,:,isign), vel_evc, g_vel_evc)
        end if
        
        ! now the j_bare term  
        call add_to_current(j_bare(:,:,:,current_spin), evc, g_vel_evc)
      enddo  ! i
      
    enddo  ! isign
  enddo  ! ik

#ifdef __PARA
  call mp_sum( f_sum, intra_pool_comm)
  call mp_sum( q_pGv, intra_pool_comm)
  call mp_sum( q_vGv, intra_pool_comm)
  call mp_sum( delta_rmc, intra_pool_comm)
  !<ceres> other reductions? </ceres>
#endif
  
#ifdef __PARA
  call mp_sum( f_sum, inter_pool_comm )
  call mp_sum( q_pGv, inter_pool_comm )
  call mp_sum( q_vGv, inter_pool_comm )
  call mp_sum( delta_g_diamagn, inter_pool_comm )
  call mp_sum( delta_g_paramagn, inter_pool_comm )
  call mp_sum( delta_rmc, inter_pool_comm )
  call mp_sum( gipaw_delta_rmc, inter_pool_comm )
  call mp_sum( j_bare, inter_pool_comm )
#endif
  
  !====================================================================
  ! print out results
  !====================================================================
  write(stdout,'(5X,''End of magnetic susceptibility calculation'',/)')
  
  ! f-sum rule
  if (iverbosity > 0) then
    write(stdout, '(5X,''f-sum rule:'')')
    write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum
  endif
  call sym_cart_tensor(f_sum)
  write(stdout, '(5X,''f-sum rule (symmetrized):'')')
  write(stdout, '(3(5X,3(F12.6,2X)/))') f_sum
  
  !--------------------------------------------------------------------
  ! now get the current, induced field and delta_g
  !--------------------------------------------------------------------
  chi_bare_pGv(:,:) = chi_bare_pGv(:,:) / omega
  j_bare(:,:,:,:) = j_bare(:,:,:,:) * alpha &
       / (2.0_dp * q_gipaw * tpiba * omega)
  
  ! either you symmetrize the current ...
  do ispin = 1, nspin
#ifdef __PARA
     call psymmetrize_field(j_bare(:,:,:,ispin),1)
#else
     call symmetrize_field(j_bare(:,:,:,ispin),1)
#endif
  end do
  
  !
  ! calculate the susceptibility
  !
  ! F_{ij} = (2 - \delta_{ij}) Q_{ij}
  if ( tevaluate_chi ) then
     do ipol = 1, 3
        do jpol = 1, 3
           f_pGv(ipol,jpol,:) = 2.0_dp*q_pGv(ipol,jpol,:)
           if (ipol == jpol) f_pGv(ipol,jpol,:) = q_pGv(ipol,jpol,:)
           
           f_vGv(ipol,jpol,:) = 2.0_dp*q_vGv(ipol,jpol,:)
           if (ipol == jpol) f_vGv(ipol,jpol,:) = q_vGv(ipol,jpol,:)
        enddo
     enddo
     
     ! compute chi_bare both pGv and vGv terms
     chi_bare_pGv(:,:) = f_pGv(:,:,1) - 2.0_dp*f_pGv(:,:,0) + f_pGv(:,:,-1)
     chi_bare_pGv(:,:) = -0.50_dp * chi_bare_pGv(:,:) * alpha &
          / ( q_gipaw * tpiba )**2
     if (iverbosity > 0) then
        write(stdout, '(5X,''chi_bare pGv (HH) in paratec units:'')')
        write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) / alpha ** 2
     endif
     call sym_cart_tensor(chi_bare_pGv)
     if (iverbosity > 0) then
        write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_pGv(:,:) / alpha ** 2
     endif
     
     chi_bare_vGv(:,:) = f_vGv(:,:,1) - 2.0_dp*f_vGv(:,:,0) + f_vGv(:,:,-1)
     chi_bare_vGv(:,:) = -0.50_dp * chi_bare_vGv(:,:) * alpha &
          / ( q_gipaw * tpiba)**2
     if (iverbosity > 0) then
        write(stdout, '(5X,''chi_bare vGv (VV) in paratec units:'')')
        write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) / alpha ** 2
     endif
     call sym_cart_tensor(chi_bare_vGv)
     if (iverbosity > 0) then
        write(stdout, '(3(5X,3(F12.6,2X)/))') chi_bare_vGv(:,:) / alpha ** 2
     endif
  end if
  
  ! compute induced field
  do ipol = 1, 3
    call biot_savart( ipol )
  enddo
  
  d_omega = omega / REAL ( nr1 * nr2 * nr3, DP )
  
  !
  ! ***************** spin-orbit-bare *******************
  !
  
  ! <apsi> TMPTMPTMP PLEASE CHECK FOR VANDERBILT/HARD GRIDS
  allocate ( grad_vr ( 3, nrxx ), v_local ( nrxx, nspin ) )
  do ispin = 1, nspin
     v_local(:,ispin) = vltot(:) + v%of_r(:,ispin)
  end do
  ! <ceres> which spin channel? </ceres>
  call gradient ( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, v_local, &
       ngm, g, nl, grad_vr )
  grad_vr = grad_vr * units_Ry2Ha
  deallocate ( v_local )
  ! </apsi>
  
  do ipol = 1, 3
     delta_g_bare ( ipol, 1 ) = SUM ( &
          ( j_bare(:,2,ipol,s_maj)-j_bare(:,2,ipol,s_min)) * grad_vr ( 3, : ) &
          -(j_bare(:,3,ipol,s_maj)-j_bare(:,3,ipol,s_min)) * grad_vr ( 2, : ) )
     delta_g_bare ( ipol, 2 ) = SUM ( &
          ( j_bare(:,3,ipol,s_maj)-j_bare(:,3,ipol,s_min)) * grad_vr ( 1, : ) &
          -(j_bare(:,1,ipol,s_maj)-j_bare(:,1,ipol,s_min)) * grad_vr ( 3, : ) )
     delta_g_bare ( ipol, 3 ) = SUM ( &
          ( j_bare(:,1,ipol,s_maj)-j_bare(:,1,ipol,s_min)) * grad_vr ( 2, : ) &
          -(j_bare(:,2,ipol,s_maj)-j_bare(:,2,ipol,s_min)) * grad_vr ( 1, : ) )
  end do
  deallocate ( grad_vr )
  
#ifdef __PARA
  call mp_sum ( delta_g_bare, intra_pool_comm )
#endif
  
  delta_g_bare = delta_g_bare * d_omega
  delta_g_bare = - delta_g_bare * g_prime / 2.0_dp * alpha * 1.0d6
  
  !
  ! ***************** spin-other-orbit *******************
  !
  
  ! This is the form used in the implementation in 'paratec'
  
  ! calculate the spin-other-orbit term a'la paratec:
  !   int_r j_up(r) x v_h[n_unpaired] d^3r
  allocate(grad_vh(3,nrxx), vh(nrxx,nspin), aux1(nrxx))
  ! transform rho to G-space
  aux1(:) = rho%of_r(:,s_maj)-rho%of_r(:,s_min)
  CALL cft3( aux1, nr1, nr2, nr3, nrx1, nrx2, nrx3, -1 )
  rho%of_g(:,1) = aux1(nl(:))
  rho%of_g(:,2) = 0.d0
  vh = 0.d0
  call v_h(rho%of_g, e_hartree, charge, vh)
  ! <ceres> which spin channel? </ceres>
  call gradient ( nrx1, nrx2, nrx3, nr1, nr2, nr3, nrxx, vh, &
       ngm, g, nl, grad_vh )
  grad_vh = grad_vh * units_Ry2Ha
  deallocate(vh, aux1)
  
  ! -2*j_bare_dw(r) because of the self-interaction correction:
  !   j_bare(r) - [j_bare_up(r)-j_bare_dw(r)] = -2*j_bare_dw(r)
  do ipol = 1, 3
     delta_g_soo ( ipol, 1 ) = 2 * SUM &
          ( j_bare(:,2,ipol,s_min) * grad_vh ( 3, : ) &
          - j_bare(:,3,ipol,s_min) * grad_vh ( 2, : ) )
     delta_g_soo ( ipol, 2 ) = 2 * SUM &
          ( j_bare(:,3,ipol,s_min) * grad_vh ( 1, : ) &
          - j_bare(:,1,ipol,s_min) * grad_vh ( 3, : ) )
     delta_g_soo ( ipol, 3 ) = 2 * SUM &
          ( j_bare(:,1,ipol,s_min) * grad_vh ( 2, : ) &
          - j_bare(:,2,ipol,s_min) * grad_vh ( 1, : ) )
  end do
  deallocate ( grad_vh )
  
#ifdef __PARA
  call mp_sum ( delta_g_soo, intra_pool_comm )
#endif
  
  delta_g_soo = delta_g_soo * d_omega
  delta_g_soo = delta_g_soo * 2 * alpha * 1.0d6
  
  ! This is obtained using the equation (7) of Pickard et Mauri, PRL 88/086043
  if ( tcalculate_correct_delta_g_soo ) then
     !<apsi> The G=0 term is still not yet in
     !</apsi>
     do jpol = 1, 3
        do ipol = 1, 3
           delta_g_soo_2 ( ipol, jpol ) = SUM ( b_ind_r(:,ipol,jpol) &
                * ( rho%of_r(:,s_maj)-rho%of_r(:,s_min) ) )
        end do
     end do
#ifdef __PARA
     call mp_sum ( delta_g_soo_2, intra_pool_comm )
#endif
  else
     delta_g_soo_2 = 0.0_dp
  end if
  
  delta_g_soo_2 = delta_g_soo_2 * d_omega
  delta_g_soo_2 = delta_g_soo_2 * 2 * 1.0d6
  
  !
  ! ***************** relativistic-mass-correction *******************
  !
  
  delta_rmc = delta_rmc * alpha ** 2 * g_e * units_Ry2Ha * 1.0d6
  delta_g_rmc = 0.0_DP
  do i = 1, 3
     delta_g_rmc(i,i) = delta_rmc
  end do
  
  ! ***************** relativistic-mass-correction gipaw *******************
  gipaw_delta_rmc = gipaw_delta_rmc * alpha ** 2 * g_e * units_Ry2Ha * 1e6
  delta_g_rmc_gipaw = 0.0_DP
  do i = 1, 3
     delta_g_rmc_gipaw(i,i) = gipaw_delta_rmc
  end do
  
  ! ***************** diamagnetic reconstruction *******************
  delta_g_diamagn(:,:) = delta_g_diamagn(:,:) &
       * g_prime / 4 * units_Ry2Ha * 1e6
  
  ! ***************** paramagnetic reconstruction *******************
  delta_g_paramagn(:,:) = delta_g_paramagn(:,:) * g_prime / 2 * 1e6 &
       * units_Ry2Ha
  
  ! ***************** total delta_g *******************
  delta_g_total = delta_g_rmc + delta_g_rmc_gipaw + delta_g_bare &
       + delta_g_diamagn + delta_g_paramagn + delta_g_soo 
  

  ! print results
  write (stdout,*)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - relativistic-mass-correction'
  write (stdout,*) 
  write (stdout, '(3(5X,3(F12.4,2X)/))' ) delta_g_rmc(:,:)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - relativistic-mass-correction gipaw'
  write (stdout,*) 
  write (stdout, '(3(5X,3(F12.4,2X)/))' ) delta_g_rmc_gipaw(:,:)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-orbit-bare'
  write (stdout,*)
  if (iverbosity > 0) &
    write (stdout, '(3(5X,3(F12.4,2X)/),/)' ) delta_g_bare(:,:)
  call sym_cart_tensor(delta_g_bare) 
  write (stdout, '(3(5X,3(F12.4,2X)/))' ) delta_g_bare(:,:)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-orbit diamagnetic correction (GIPAW)'
  write (stdout,*) 
  if (iverbosity > 0) &
    write ( stdout, '(3(5X,3(F12.4,2X)/),/)' ) delta_g_diamagn(:,:)
  call sym_cart_tensor(delta_g_diamagn) 
  write (stdout, '(3(5X,3(F12.4,2X)/))' ) delta_g_diamagn(:,:)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-orbit paramagnetic correction (GIPAW)'
  write (stdout,*) 
  if (iverbosity > 0) &
    write (stdout, '(3(5X,3(F12.4,2X)/),/)' ) delta_g_paramagn(:,:)
  call sym_cart_tensor(delta_g_paramagn) 
  write (stdout, '(3(5X,3(F12.4,2X)/))' ) delta_g_paramagn(:,:)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-other-orbit'
  write (stdout,*) 
  if (iverbosity > 0) &
    write (stdout, '(3(5X,3(F12.4,2X)/),/)' ) delta_g_soo(:,:)
  call sym_cart_tensor(delta_g_soo) 
  write (stdout, '(3(5X,3(F12.4,2X)/))' ) delta_g_soo(:,:)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - spin-other-orbit, version 2'
  write (stdout,*) 
  if (iverbosity > 0) &
    write (stdout, '(3(5X,3(F12.4,2X)/),/)' ) delta_g_soo_2(:,:)
  call sym_cart_tensor(delta_g_soo_2) 
  write (stdout, '(3(5X,3(F12.4,2X)/))' ) delta_g_soo_2(:,:)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  write (stdout,*) 'Delta g - total'
  write (stdout,*) 
  if (iverbosity > 0) &
    write (stdout, '(3(5X,3(F12.4,2X)/),/)' ) delta_g_total(:,:)
  call sym_cart_tensor(delta_g_total) 
  write (stdout, '(3(5X,3(F12.4,2X)/))' ) delta_g_total(:,:)

  write (stdout,*) '**********************************************'
  write (stdout,*)
  
  deallocate(p_evc, vel_evc, aux, g_vel_evc, j_bare, b_ind)
  
CONTAINS
  
  !====================================================================
  ! compute p_k|evc>, v_k|evc> and G_k v_{k+q,k}|evc>
  !====================================================================
  SUBROUTINE apply_operators
    implicit none
    integer ipol

    do ipol = 1, 3
      call apply_p(evc, p_evc(1,1,ipol), ik, ipol, q)
      call apply_vel(evc, vel_evc(1,1,ipol), ik, ipol, q)
      ! necessary because aux is overwritten by subroutine greenfunction
      aux(:,:) = vel_evc(:,:,ipol)
      call greenfunction(ik, aux, g_vel_evc(1,1,ipol), q)
    enddo
  END SUBROUTINE apply_operators


  !====================================================================
  ! add contribution the Q tensors
  ! Q_{\alpha,\beta} += <(e_i \times ul)_\alpha | (e_i \times ur)_\beta>
  !====================================================================
  SUBROUTINE add_to_tensor(qt, ul, ur)
    implicit none
    real(dp), intent(inout) :: qt(3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd,3), ur(npwx,nbnd,3)
    real(dp) :: braket
    integer :: ibnd, ia, ib, comp_ia, comp_ib, ind(3,3), mult(3,3)
    
    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)

    do ia = 1, 3    ! ia = alpha
      comp_ia = ind(ia,i)
      if (mult(ia,i) == 0) cycle

      do ib = 1, 3    ! ib = beta
        comp_ib = ind(ib,i)
        if (mult(ib,i) == 0) cycle

        do ibnd = 1, nbnd_occ(ik)
          braket = real(ZDOTC(npw, ul(1,ibnd,comp_ia), 1, &
                                   ur(1,ibnd,comp_ib), 1), DP)
          qt(ia,ib) = qt(ia,ib) + wg(ibnd,ik) * &
                      braket * mult(ia,i) * mult(ib,i)
        enddo  ! ibnd

      enddo  ! ib
    enddo  ! ia
  END SUBROUTINE add_to_tensor
  
  
  !====================================================================
  ! add contribution the the current
  ! j(r)_{\alpha,\beta} += <ul|J(r)|(B\times e_i \cdot ur)>
  !====================================================================
  SUBROUTINE add_to_current(j, ul, ur)
    implicit none
    real(dp), intent(inout) :: j(nrxxs,3,3)
    complex(dp), intent(in) :: ul(npwx,nbnd), ur(npwx,nbnd,3)
    real(dp) :: fact
    integer :: ibdir, icomp, ind(3,3), mult(3,3)
    
    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)
    
    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir, i)
      fact = real(mult(ibdir,i)*isign)
      call j_para(fact, ul(1,1), ur(1,1,icomp), ik, q, j(1,1,ibdir))
    enddo
  END SUBROUTINE add_to_current
  
  !====================================================================
  ! ...
  !====================================================================
  SUBROUTINE relativistic_mass_correction ( rmc_gipaw )
    
    USE ions_base,      ONLY : nat, ityp, ntyp => nsp
    USE gipaw_module,   ONLY : radial_integral_rmc
    
    implicit none
    
    ! Arguments
    real(dp), intent(inout):: rmc_gipaw

    integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
    integer :: nt, ibnd, na, lm, nrc, ijkb0
    complex(dp) :: rmc_corr
    complex(dp) :: bec_product
    
    rmc_corr = 0.0_dp
    
    do ibnd = 1, nbnd
       ijkb0 = 0
       do nt = 1, ntyp
          do na = 1, nat
             if (ityp (na) == nt) then
                do ih = 1, paw_recon(nt)%paw_nh
                   ikb = ijkb0 + ih
                   nbs1=paw_recon(nt)%paw_indv(ih)
                   l1=paw_recon(nt)%paw_nhtol(ih)
                   m1=paw_recon(nt)%paw_nhtom(ih)
                   lm1=m1+l1**2
                   do jh = 1, paw_recon(nt)%paw_nh
                      jkb = ijkb0 + jh
                      nbs2=paw_recon(nt)%paw_indv(jh)
                      l2=paw_recon(nt)%paw_nhtol(jh)
                      m2=paw_recon(nt)%paw_nhtom(jh)
                      lm2=m2+l2**2
                      
                      IF ( l1 /= l2 ) CYCLE
                      IF ( m1 /= m2 ) CYCLE
                      
                      bec_product = paw_becp(jkb,ibnd) &
                           * CONJG( paw_becp(ikb,ibnd) )
                      
                      rmc_corr = rmc_corr &
                              + bec_product &
                              * radial_integral_rmc(nbs1,nbs2,nt) &
                              * wg(ibnd,ik)
                      
                   enddo
                enddo
                ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
             endif
          enddo
       enddo
    enddo
    
    rmc_gipaw = REAL ( rmc_corr, dp )
    
  END SUBROUTINE relativistic_mass_correction
  
  
  !====================================================================
  ! ...
  !====================================================================
  SUBROUTINE diamagnetic_correction ( diamagnetic_tensor )
    
    USE ions_base,      ONLY : nat, ityp, ntyp => nsp
    USE gipaw_module,   ONLY : radial_integral_diamagnetic_so
    
    implicit none
    
    ! Arguments
    real(dp), intent(inout):: diamagnetic_tensor(3,3)

    integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
    integer :: nt, ibnd, na, lm, nrc, ijkb0
    complex(dp) , allocatable :: dia_corr(:)
    complex(dp) :: bec_product
    
    allocate ( dia_corr ( lmaxx**2 ) )
    dia_corr = 0.0_dp
    
    do ibnd = 1, nbnd
       ijkb0 = 0
       do nt = 1, ntyp
          do na = 1, nat
             if (ityp (na) == nt) then
                do ih = 1, paw_recon(nt)%paw_nh
                   ikb = ijkb0 + ih
                   nbs1 = paw_recon(nt)%paw_indv(ih)
                   l1 = paw_recon(nt)%paw_nhtol(ih)
                   m1 = paw_recon(nt)%paw_nhtom(ih)
                   lm1 = m1+l1**2
                   do jh = 1, paw_recon(nt)%paw_nh
                      jkb = ijkb0 + jh
                      nbs2 = paw_recon(nt)%paw_indv(jh)
                      l2 = paw_recon(nt)%paw_nhtol(jh)
                      m2 = paw_recon(nt)%paw_nhtom(jh)
                      lm2 = m2 + l2**2
                      
                      bec_product = paw_becp(jkb,ibnd) &
                           * CONJG( paw_becp(ikb,ibnd) )
                      
                      !<apsi> s/non-trace-zero component
                      ! 2/3 to separate the non-trace vanishing component
                      ! 1/(2c^2) from the equation (59) in PM-PRB
                      IF ( l1 == l2 .AND. m1 == m2 ) THEN
                         diamagnetic_tensor(1,1) &
                              = diamagnetic_tensor(1,1) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2
                         diamagnetic_tensor(2,2) &
                              = diamagnetic_tensor(2,2) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2
                         diamagnetic_tensor(3,3) &
                              = diamagnetic_tensor(3,3) &
                              + 2.0_dp / 3.0_dp * bec_product &
                              * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                              * wg(ibnd,ik) * alpha ** 2
                      END IF
                      
                      ! 2/3 to separate the non-trace vanishing component
                      do lm = 5, 9
                         dia_corr(lm) =  dia_corr(lm) &
                              + bec_product / 3.0_dp &
                              * radial_integral_diamagnetic_so(nbs1,nbs2,nt) &
                              * ap(lm,lm1,lm2) * wg(ibnd,ik) * alpha ** 2
                      enddo
                   enddo
                enddo
                ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
             endif
          enddo
       enddo
    enddo
    
    if ( iverbosity > 10 ) then
       write(stdout,'("Debug: gt, dia  ",5F14.8)') dia_corr(5:9)
    end if
    
    !
    !  transform in cartesian coordinates
    !
    
    dia_corr(5:9) = - sqrt(4.0_dp * pi/5.0_dp) * dia_corr(5:9)
    
    diamagnetic_tensor(1,1) = diamagnetic_tensor(1,1) &
         + sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(2,2) = diamagnetic_tensor(2,2) &
         - sqrt(3.0_dp) * dia_corr(8) - dia_corr(5)
    diamagnetic_tensor(3,3) = diamagnetic_tensor(3,3) &
         + dia_corr(5) * 2.0_dp
    diamagnetic_tensor(1,2) = diamagnetic_tensor(1,2) &
         +  dia_corr(9) * sqrt(3.0_dp)
    diamagnetic_tensor(2,1) = diamagnetic_tensor(1,2)
    diamagnetic_tensor(1,3) = diamagnetic_tensor(1,3) &
         - dia_corr(6) * sqrt(3.0_dp)
    diamagnetic_tensor(3,1) = diamagnetic_tensor(1,3)
    diamagnetic_tensor(2,3) = diamagnetic_tensor(2,3) &
         - dia_corr(7) * sqrt(3.0_dp)
    diamagnetic_tensor(3,2) = diamagnetic_tensor(2,3)
    
    ! dia_corr(5) = 3z^2-1
    ! dia_corr(6) = -xz
    ! dia_corr(7) = -yz
    ! dia_corr(8) = x^2-y^2
    ! dia_corr(9) = xy
    
    deallocate ( dia_corr )
    
  END SUBROUTINE diamagnetic_correction
  

  !====================================================================
  ! ...
  !====================================================================
  SUBROUTINE paramagnetic_correction ( paramagnetic_tensor )
    
    USE ions_base,      ONLY : nat, ityp, ntyp => nsp
    USE gipaw_module,   ONLY : lx, ly, lz, paw_becp2, &
                              radial_integral_paramagnetic_so
    USE becmod,         ONLY : calbec
    
    implicit none
    
    ! Arguments
    real(dp), intent(inout):: paramagnetic_tensor(3,3)
    
    integer :: l1, m1, lm1, l2, m2, lm2, ih, ikb, nbs1, jh, jkb, nbs2
    integer :: nt, ibnd, na, lm, j, ijkb0, ipol
    complex(dp) :: bec_product
    complex(dp) , allocatable :: para_corr(:)
    
    integer, parameter :: ng_ = 27, lmax2_ = 16
    integer :: mg, i1, i2, i3
    real(DP) :: g_ (3, ng_), gg_ (ng_)
    real(DP) :: ylm_ (ng_,lmax2_)
    
    !--------------------------------------------------------------------------
    
    allocate (para_corr(3))
    
    !
    !  calculation of the reconstruction part
    !
    
    do ipol = 1, 3 
       
       if ( ipol == i ) cycle !TESTTESTTEST
       
       !call ccalbec ( paw_nkb, npwx, npw, nbnd, paw_becp2, paw_vkb, &
       !     g_vel_evc(1,1,ipol) )
       call calbec (npw, paw_vkb, g_vel_evc(:,:,ipol), paw_becp2)
       
       para_corr = 0.0_dp
       
       do ibnd = 1, nbnd
          ijkb0 = 0
          do nt = 1, ntyp
             do na = 1, nat
                
                if (ityp (na) == nt) then
                   do ih = 1, paw_recon(nt)%paw_nh
                      ikb = ijkb0 + ih
                      nbs1 = paw_recon(nt)%paw_indv(ih)
                      l1 = paw_recon(nt)%paw_nhtol(ih)
                      m1 = paw_recon(nt)%paw_nhtom(ih)
                      lm1 = m1 + l1**2
                      
                      do jh = 1, paw_recon(nt)%paw_nh
                         jkb = ijkb0 + jh
                         nbs2 = paw_recon(nt)%paw_indv(jh)
                         l2 = paw_recon(nt)%paw_nhtol(jh)
                         m2 = paw_recon(nt)%paw_nhtom(jh)
                         lm2 = m2 + l2**2
                         
                         if ( l1 /= l2 ) cycle
                         
                         bec_product = CONJG(paw_becp(ikb,ibnd)) &
                              * paw_becp2(jkb,ibnd)
                         
                         para_corr(1) = para_corr(1) &
                              + bec_product &
                              * radial_integral_paramagnetic_so(nbs1,nbs2,nt) &
                              * lx ( lm1, lm2 ) * wg(ibnd,ik) * alpha ** 2
                         para_corr(2) = para_corr(2) &
                              + bec_product &
                              * radial_integral_paramagnetic_so(nbs1,nbs2,nt) &
                              * ly ( lm1, lm2 ) * wg(ibnd,ik) * alpha ** 2
                         para_corr(3) = para_corr(3) &
                              + bec_product &
                              * radial_integral_paramagnetic_so(nbs1,nbs2,nt) &
                              * lz ( lm1, lm2 ) * wg(ibnd,ik) * alpha ** 2
                      enddo
                   enddo
                   ijkb0 = ijkb0 + paw_recon(nt)%paw_nh
                endif
             enddo
          enddo
       enddo
       
       paramagnetic_tensor ( :, ipol ) = REAL ( para_corr, dp )
       
       if ( iverbosity > 10 ) then
          write(stdout,'("Debug: gt, para ",2I3,3(F16.7,2X))') &
               ipol, i*isign, REAL ( para_corr(1:3) ) * 1e6
       end if
       
    end do
    
    deallocate( para_corr )
    
  END SUBROUTINE paramagnetic_correction
  
  !====================================================================
  ! ...
  !====================================================================
  SUBROUTINE add_to_sigma_para( paramagnetic_correction, sigma_paramagnetic )
    implicit none
    real(dp), intent(in) :: paramagnetic_correction(3,3)
    real(dp), intent(inout) :: sigma_paramagnetic(3,3)
    real(dp) :: fact
    integer :: ibdir, icomp, ipol, ind(3,3), mult(3,3)
    
    ! index for the cross product
    ind(:,1) = (/ 1, 3, 2/);  mult(:,1) = (/ 0,-1, 1 /)
    ind(:,2) = (/ 3, 2, 1/);  mult(:,2) = (/ 1, 0,-1 /)
    ind(:,3) = (/ 2, 1, 3/);  mult(:,3) = (/-1, 1, 0 /)
    
    ! loop over B direction
    do ibdir = 1, 3
      if (i == ibdir) cycle
      icomp = ind(ibdir, i)
      fact = real(mult(ibdir,i)*isign)
      
      do ipol = 1, 3
         sigma_paramagnetic ( ipol, icomp ) &
              = sigma_paramagnetic ( ipol, icomp ) &
              + fact * paramagnetic_correction ( ipol, ibdir ) &
              / ( 2 * q_gipaw * tpiba )
         
!   Do iq=1, 3    ! loop over all q-points
!      ...
!      Do iperm0 =1,-1,-2
!         b0 = Mod(iq+iperm0+2,3) +1 ! gives different b0 for iperm0
!         p0 = Mod(iq-iperm0+2,3) +1 ! gives p0 = abs(q x b0)
!         ...
!         Call take_nonloc_deriv_kq_k(p0,k_gspace, rk(1),u_k(1,i),&
!         ...
!         para_corr_shift_tot(b0,:,:,:) = para_corr_shift_tot(b0,:,:,:)-&
!              Real(iperm0,dp)*Real(iqsign,dp)*para_corr_shift*&
!              kpoints%w(irk)/qmag/dtwo/Real(crys%nspin,dp)
         
      end do
    enddo
  END SUBROUTINE add_to_sigma_para
  
END SUBROUTINE g_tensor_crystal

