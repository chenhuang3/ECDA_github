module comm_data

  implicit  none 
  
  logical :: do_envOF = .false. , &    ! env is treated by OFDFT 
             do_pbc   = .false. ,&     ! PBC is activated
             precond_vks = .false., &  ! precondition vKS
             do_env_correction = .true., & 
             share_atom = .false.  !

  integer :: logOF = 9945       ! unit for OFDFT 
  integer :: vwreg_type = 1     ! 1: rho' = rho + c
                                ! 2: rho' = [1-exp(-rho/c)]*rho
                                ! 3: rho' = [1-exp(-sqrt(rho)/c)]*rho

  real(8) :: vw_lam = 1.d0/9.d0,  &  ! coefficient for VW KEDF 
             vw_reg = 1e-5           ! threshold for density to be regularized

  integer :: dfet_mix_type = 1, & ! 1: pot mix, 2: den mix
             dfet_pen_type = 1    ! 1=> Yakuwa potential 
                                  ! 2=> (rho1+rho2-rho_ref)
end module 
