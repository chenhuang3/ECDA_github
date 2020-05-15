!
! Compute sigma|x> with sigma defined as
!
! sigma = [A_clu'^-1 - I]*\chi_cluster + [A_env'^-1 - I]*\chi_env
!
! Input: vec
! Output: vec_out
!
subroutine sigma_vec(nfft, nspin, j_atom, ucvol, cluster_weight, env_weight, & 
                     dfermi_dvks, vks_clu, vks_env, vec, vec_out)

     implicit none 

     integer :: nfft, nspin, j_atom
     real(8) :: vks_clu(nfft,nspin), & 
                vks_env(nfft,nspin), &
                cluster_weight(nfft), & 
                env_weight(nfft), & 
                dfermi_dvks(nfft,nspin,2), & 
                vec(nfft,nspin), & 
                vec_out(nfft,nspin), & 
                ucvol


     ! local vars 
     integer :: system_type, isp 
     real(8) :: tmpvec2(nfft,nspin), & 
                tmpvec(nfft) 
                


     !>>>>>>>>>>>> function <<<<<<<<<<<<!

     ! get [A_clu'^-1 - I]*\chi_cluster * z 
     system_type = 1 ! cluster
     call apply_chi(-1,system_type,j_atom,nfft,nspin,vks_clu,vec,tmpvec2)
     vec_out = -tmpvec2
     do isp = 1,nspin 
       call apply_Ainv(nfft,ucvol,.true.,cluster_weight,dfermi_dvks(:,isp,1),tmpvec2(:,isp),tmpvec)
       vec_out(:,isp) = vec_out(:,isp) + tmpvec
     enddo 

     ! get [A_env'^-1 - I]*\chi_env * z 
     system_type = 2 ! env 
     call apply_chi(-1,system_type,j_atom,nfft,nspin,vks_env,vec,tmpvec2)
     vec_out = vec_out - tmpvec2
     do isp = 1,nspin 
       call apply_Ainv(nfft,ucvol,.true.,env_weight,dfermi_dvks(:,isp,2),tmpvec2(:,isp),tmpvec)
       vec_out(:,isp) = vec_out(:,isp) + tmpvec
     enddo 

  end subroutine 

