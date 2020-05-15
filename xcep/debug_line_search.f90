  subroutine debug_linesearch(nfft,nspin)
     
     integer :: nfft, nspin
     real(8) :: oldg(nfft,nspin),oldu(nfft,nspin),theta

!!!!!!!!!!!!!!!!!! ===============================

       ! run subsystems
       call dump_uemb_samecell(u,nfft,nspin)
       do s=1,nsys
         call run_subsystem(s,1)
       enddo
     
       ! wait the scf of subsystem to finish, get eigenvalues 
       ! recompute fermi level by taking all subsystem into considerations
       ! tell abinit the fermi level of the whole system 
       do s=1,nsys
         call wait_for_system_scf(s)
       enddo
       call fermi_level_across_subsystem(nsys,nspin,nelectr,tsmear)
       
       ! after updating occ, now get all data from subsystems
       do s=1,nsys
         call wait_for_system_end(s)
         call get_subsystem_data(s,nfft,dvol,nspin,sub_etotal(s),&
              sub_nlps(s),sub_TS(s),sub_ke(s),sub_rhor(:,:,s))
       enddo

       call check_subsystem_charges(dvol,nspin,nfft,nsys,ref_rho,sub_rhor)
   
       ! get cgW 
       cgW = - (sum(sub_etotal) + dvol*sum((sum(sub_rhor,3)-ref_rho)*u))
       ! penalty function
       pen = 0.d0
       do s=1,nspin
         call laplacian(nfft,cell_nfft(1),cell_nfft(2),cell_nfft(3),u(:,s),qvec,lap)
         pen = pen - pen_coeff*sum(u(:,s)*lap)*dvol
       enddo
       cgW = cgW + pen
       print *,'(main) cgW: ',cgW,'   penalty: ',pen
       
       ! get gradient 
       grad = ref_rho - sum(sub_rhor,3) 
       ! gradient due to penality 
       do s=1,nspin
         call laplacian(nfft,cell_nfft(1),cell_nfft(2),cell_nfft(3),u(:,s),qvec,lap)
         grad(:,s) = grad(:,s) - 2.d0*pen_coeff*lap
       enddo

!!!!!!!!!!!!!!!!!!

     oldg = grad
     oldu = u
     theta = 0.d0
     print *,'cgW: ',cgW,' penalty: ',pen,' dot(g,g): ',sum(grad*grad)*dvol
    ! open(file='dir_lbfgs.dat',unit=111,action='read',form='unformatted')
    ! read(111) oldg
    ! close(111)
    ! oldg = -oldg

!!!!!!!!!!!!!!!!!!!     
 
     !============ loop to do line search  =============
     do 
     theta = theta + 0.2d0
     u =- oldg * theta + oldu

!! ================================     

       ! run subsystems
       call dump_uemb_samecell(u,nfft,nspin)
       do s=1,nsys
         call run_subsystem(s,1)
       enddo
     
       ! wait the scf of subsystem to finish, get eigenvalues 
       ! recompute fermi level by taking all subsystem into considerations
       ! tell abinit the fermi level of the whole system 
       do s=1,nsys
         call wait_for_system_scf(s)
       enddo
       call fermi_level_across_subsystem(nsys,nspin,nelectr,tsmear)
       
       ! after updating occ, now get all data from subsystems
       do s=1,nsys
         call wait_for_system_end(s)
         call get_subsystem_data(s,nfft,dvol,nspin,sub_etotal(s), & 
              sub_nlps(s),sub_TS(s),sub_ke(s),sub_rhor(:,:,s))
       enddo

       call check_subsystem_charges(dvol,nspin,nfft,nsys,ref_rho,sub_rhor)
   
       ! get cgW 
       cgW = - (sum(sub_etotal) + dvol*sum((sum(sub_rhor,3)-ref_rho)*u))
       ! penalty function
       pen = 0.d0
       do s=1,nspin
         call laplacian(nfft,cell_nfft(1),cell_nfft(2),cell_nfft(3),u(:,s),qvec,lap)
         pen = pen - pen_coeff*sum(u(:,s)*lap)*dvol
       enddo
       cgW = cgW + pen
       print *,'(main) cgW: ',cgW,'   penalty: ',pen
       
       ! get gradient 
       grad = ref_rho - sum(sub_rhor,3) 
       ! gradient due to penality 
       do s=1,nspin
         call laplacian(nfft,cell_nfft(1),cell_nfft(2),cell_nfft(3),u(:,s),qvec,lap)
         grad(:,s) = grad(:,s) - 2.d0*pen_coeff*lap
       enddo


!! ================================

     write(6,'(a,f10.5,a,f16.8,a,f12.4,a,es12.4)') &  
     'theta: ',theta,'  cgW: ',cgW,' penalty: ',pen,' dedtheta: ',-sum(oldg*grad)*dvol
     enddo
 end subroutine
