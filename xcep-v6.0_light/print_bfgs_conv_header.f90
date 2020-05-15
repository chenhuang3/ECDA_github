subroutine print_bfgs_conv_header(file_dfet_out,global_fermi,nspin)
   integer :: file_dfet_out, global_fermi, nspin 
   !
   if (global_fermi<0) then 
     if (nspin==1) then 
       write(file_dfet_out,'(a)') & 
       'nnew          cgW          gmax        gnorm       ek_clu    ek_env   ref_diff '
       write(715,'(a)') & 
       'nnew          cgW          gmax        gnorm       ek_clu    ek_env   ref_diff '
     endif
     if (nspin==2) then 
       write(file_dfet_out,'(a)') & 
         'nnew           cgW          gmax        gnorm    de_ref(up)  de_ref(dn)    ek_clu    ek_env  penalty'
       write(715,'(a)') & 
         'nnew           cgW          gmax        gnorm    de_ref(up)  de_ref(dn)    ek_clu    ek_env  penalty'
     endif
   else
     if (nspin==2) then 
       write(file_dfet_out,'(a)') & 
         'nnew        cgW      gmax      gnorm    |rho-ref_rho|(up,down)    ek_clu   ek_env   penalty'
       write(715,'(a)') & 
         'nnew        cgW      gmax      gnorm    |rho-ref_rho|(up,down)    ek_clu   ek_env   penalty'
     else
       write(file_dfet_out,'(a)') & 
         'nnew          cgW          gmax       gnorm      ref_diff'
       write(715,'(a)') & 
         'nnew          cgW          gmax       gnorm      ref_diff'
     endif
   endif

   call flush(file_dfet_out)
   call flush(715)

end subroutine print_bfgs_conv_header
