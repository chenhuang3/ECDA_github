!===================================
! dump cluster and env structure 
!
! Created by Chen Huang 1/15/2019
!===================================
subroutine write_cluster_env_structure_XSF(natom,rprimd,xcart,ilist,nlist,mlist,nshare,ishare)

 implicit none 
 integer :: natom,  & 
            nlist(natom), ilist(mlist,natom), mlist, & 
            nshare(natom), & 
            ishare(mlist,natom)

 real(8) :: rprimd(3,3),xcart(3,natom)

 ! local vars 
 character(len=500) :: filename, ss
 logical :: in_list
 integer :: j_atom, s1, k1, shared, kk, ii

 ! loop over all clusters 
 do j_atom = 1,natom 
   write(ss,*)j_atom 
   filename = "cluster_"//trim(adjustl(ss))//".vasp"
   open(file=filename,unit=111,action='write',form='formatted')
   write(111,*)'comment line ---------'
   write(111,*)'1.0 '
   write(111,'(3f16.4)') rprimd(1:3,1)*0.5291772106
   write(111,'(3f16.4)') rprimd(1:3,2)*0.5291772106
   write(111,'(3f16.4)') rprimd(1:3,3)*0.5291772106

   if (nshare(j_atom)>0) then 
      write(111,'(a)') 'N C S O'
      write(111,'(4i4)')1,nlist(j_atom)-nshare(j_atom),nshare(j_atom),natom-nlist(j_atom)-1
   else
      write(111,'(a)') 'N C O'
      write(111,'(3i4)')1,nlist(j_atom),natom-nlist(j_atom)-1
   endif 
   write(111,'(a)') 'Cart'

   ! write atom cooridates for this cluster 
   ! itself
   write(111,'(3f16.8)')xcart(1:3,j_atom)*0.5291772016d0

   ! cluster atoms 
   if (nshare(j_atom)==0) then 
     do kk=1,nlist(j_atom)
       write(111,'(3f16.8)')xcart(1:3,ilist(kk,j_atom))*0.5291772106d0
     enddo 
   endif 
   
   ! we have shared atoms
   if (nshare(j_atom)>0) then 
     ! print non-shared atoms 
     do kk=1,nlist(j_atom) 
        shared = 0
        do s1=1,nshare(j_atom)
           if (ishare(s1,j_atom)==ilist(kk,j_atom)) then 
              shared=1
           endif 
        enddo 
        if (shared==0) write(111,'(3f16.8)')xcart(1:3,ilist(kk,j_atom))*0.5291772106d0
     enddo 

     ! print shared atoms now 
     do kk=1,nlist(j_atom) 
        shared = 0
        do s1=1,nshare(j_atom)
           if (ishare(s1,j_atom)==ilist(kk,j_atom)) then 
              shared=1
           endif 
        enddo 
        if (shared==1) write(111,'(3f16.8)')xcart(1:3,ilist(kk,j_atom))*0.5291772106d0
     enddo 
   endif 
  

   ! env atoms 
   do kk=1,natom
      if (kk==j_atom) cycle
      in_list = .false. 
      do ii=1,nlist(j_atom) 
        if (kk==ilist(ii,j_atom)) in_list = .true. 
      enddo 
      if (.not. in_list) & 
        write(111,'(3f16.8)')xcart(1:3,kk)*0.5291772106d0
   enddo 
   close(111)
 enddo




end subroutine  
