subroutine  load_xcart_cluster_env(myrank,natom,natom_clu,natom_env,xcart_clu,xcart_env)

 use comm_data

 implicit none 

 integer :: natom, & 
            natom_clu, & 
            natom_env, & 
            myrank 

 real(8) :: xcart_clu(3,natom), & 
            xcart_env(3,natom)


 ! local vars 
 character (len=500) :: fname, ss
 logical :: file_ext
 integer :: jatom, stat, itmp


 jatom = myrank+1

 write(ss,*)jatom


 ! cluster ==============
 fname = 'subsys'//trim(adjustl(ss))//'/atom_coords.dat'
 do while (.true.)
   open(file=fname,unit=111,action='read',form='unformatted',iostat=stat)
   if (stat/=0) then 
     call sleep(1)
     cycle 
   endif 
   read(111,iostat=stat) xcart_clu(:,1:natom_clu)
   if (stat/=0) then  
     close(111)
     call sleep(1)
     cycle 
   else
     close(111)
     exit
   endif 
 enddo 



 ! env  ======================
 if (do_envOF .eqv. .false.) then 
   fname = 'subsys'//trim(adjustl(ss))//'_env/atom_coords.dat'
   do while(.true.)
     open(file=fname,unit=111,action='read',form='unformatted',IOSTAT=stat)  
     if (stat/=0) then 
       call sleep(1)
       cycle 
     endif 
     read(111,iostat=stat) xcart_env(:,1:natom_env)
     if (stat/=0) then
       close(111)
       call sleep(1)
       cycle 
     else 
       close(111)
       exit
     endif 
   enddo 
 endif 

end subroutine 
