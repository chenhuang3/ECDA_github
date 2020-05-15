

subroutine proj_dump_cluster_env_vks(j_atom,nspin,nfft,nsys,total_fermi, & 
                                     prj_vks,vac,clu_weight,env_weight,global_vks)
 use comm_data 

 implicit none 

 integer :: nsys, nfft, nspin, j_atom
 real(8) :: clu_weight(nfft), & 
            env_weight(nfft), & 
            vac, & 
            total_fermi(nspin), & 
            global_vks(nfft,nspin), & 
            w(nfft), & 
            prj_vks(nfft,nspin,2) ! one for cluster and another one for env
 !            
 integer :: isp 
 character(len=200) :: ss
 




 ! ================= cluster  ==================
 do isp=1,nspin
   prj_vks(:,isp,1) = (global_vks(:,isp)-total_fermi(isp)-vac)*clu_weight
 enddo
 write(ss,'(i4)')j_atom
 ss='./subsys'//trim(adjustl(ss))
 call chdir(ss)

 open(file='cluster_vks.dat',action='write',form='unformatted',unit=111)
 write(111)prj_vks(:,:,1)
 close(111)
 call chdir('../')
 



 ! ============== environment =============
 if (.not. do_envOF) then 
   write(ss,'(i4)')j_atom
   ss='./subsys'//trim(adjustl(ss))//'_env'
   call chdir(ss)
 
   do isp=1,nspin
     prj_vks(:,isp,2) = (global_vks(:,isp)-total_fermi(isp)-vac)*env_weight
   enddo
   open(file='env_vks.dat',action='write',form='unformatted',unit=111)
   write(111)prj_vks(:,:,2)
   close(111)
   call chdir('../')
 endif 
 
end subroutine proj_dump_cluster_env_vks
