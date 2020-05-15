
subroutine  load_atom_info(natom,xcart,znucl,atom_info,mlist,nlist,ilist,nshare,ishare)

  use mpi
  use comm_data, only: share_atom 

  implicit none 

  character(len=100) :: string
  character(len=20000) :: one_line

  integer :: natom,ii,file_unit=200,iost,int_znucl(natom), & 
             jj,mlist,nlist(natom),ilist(mlist,natom), & 
             nshare(natom), & 
             ishare(mlist,natom)

  real(8) :: xcart(3,natom), znucl(natom), & 
             atom_info(5,natom)

  integer :: itmp,j
  character(len=100) :: stmp
  integer :: myrank, ierr, counter, iatom, idx, i1, i2

  call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )

  nshare = 0
  ishare = 0
  
  open(file='param.in',unit=file_unit,action='read',iostat=iost)

  if (iost/=0) then
    print *,'error on opening param.in file, iostat=',iost
    stop
    close(file_unit)
  endif

  if (myrank==0) print *,'natom: ',natom

  do while (.true.) 

    read(file_unit,*,iostat=iost) string 
    if (iost<0) exit                 ! end of the file 

    ! load information of each atom, which subsystem they are in and their type.
    if (string(1:9) .eq. 'atom_info' ) then 
      if (myrank==0) then 
        print *,''
        print *,'  --------------------------- atom info ------------------------------'
        print *,'  subsystem_ID    atom_type       X           Y           Z         '
      endif 
      !
      do ii=1,natom
        read(file_unit,*)atom_info(1:5,ii),nlist(ii)
        if (nlist(ii)/=0) then 
           backspace  file_unit 
           read(file_unit,*)atom_info(1:5,ii),nlist(ii)
           if (nlist(ii)==999) then 
            !  print *,'**** for atom ',ii,' *** no environmental atoms '
              !
              ! no environment atoms, all atoms belong to cluster
              ! we set the number of buffer atoms to natom - 1
              ! all the atoms in teh system are assigned as buffer atoms 
              !
              nlist(ii) = natom-1
              counter = 0
              ! assign ilist(:)
              do iatom=1,natom 
                 if (iatom/=ii) then 
                    counter = counter + 1
                    ilist(counter,ii) = iatom
                 endif 
              enddo
           else
              backspace file_unit
              read(file_unit,*)atom_info(1:5,ii),nlist(ii),ilist(1:nlist(ii),ii)
              !
              ! check if the list contains the centeral atom 
              !
              do i1=1,nlist(ii) 
                if (ilist(i1,ii)==ii) then 
                   print *,'cluster index : ',ii
                   print *,'the buffer list contains the centeral atom, stop! chekc your param.ini file'
                   stop
                endif 
              enddo 
              !
              ! check if there are two same atom indexes
              !
              do i1=1,nlist(ii)
                do i2=1,nlist(ii) 
                  if (i1==i2) cycle 
                  if (ilist(i1,ii)==ilist(i2,ii)) then 
                    print *,'error, one atom appears twice in the list! stop check your param.in file '
                    print *,'atom index : ',ii
                    print *,'they are ',ilist(i1,ii),' and ',ilist(i2,ii)
                    print *,'ilist: ',ilist(1:,nlist(ii))
                    call flush(6)
                    stop
                  endif 
                enddo
              enddo 
           endif
        endif
        if (myrank==0) then 
          write(6,'(a,i3,a,i3,a,f12.5,a,f12.5,a,f12.5,a)')& 
          '   ',int(atom_info(1,ii)),'               ',int(atom_info(2,ii)), & 
          '   ',atom_info(3,ii),'',atom_info(4,ii),'',atom_info(5,ii),'    Angstrom'
        endif
        !
        znucl(ii) = atom_info(2,ii) 
        xcart(1:3,ii) = atom_info(3:5,ii)/0.52917721067d0
        !
      enddo
      if (myrank==0)  print *,''
    endif 


    !====================
    ! load shared atoms 
    !====================
    if (share_atom) then 
      if (string(1:16) .eq. 'block_share_atom' ) then 
        do j=1,natom
          read(file_unit,*)itmp,stmp,nshare(j),stmp
          backspace file_unit
          read(file_unit,*)itmp,stmp,itmp,stmp,ishare(1:nshare(j),j)
        enddo
      endif 
    endif 

  enddo 



  if (myrank==0) then 
    print *,''
    print *, '---------- buffer atoms ----------'
    write(9944,*) ''
    write(9944,'(a)') '---------- buffer atoms ---------'
    do ii=1,natom
      write(6,'(a,i3,a)',advance='no')'atom: ',ii,' buffer atoms: '
      write(9944,'(a,i3,a)',advance='no')'atom: ',ii,' buffer atoms: '
      if ( nlist(ii)/=0) then 
        do jj=1,nlist(ii)
          write(6,'(i4)',advance='no') ilist(jj,ii) 
          write(9944,'(i4)',advance='no') ilist(jj,ii) 
        enddo
      else
        write(6,'(a)',advance='no') ' none'
        write(9944,'(a)',advance='no') ' none'
      endif 
      write(6,*)''
      write(9944,*)''
    enddo
    print *,''
    write(9944,*)''
  endif 


  !====================== 
  ! shared atoms 
  !====================== 
  if(myrank==0) then 
    print *,''
    print *, '---------- shared atoms ----------'
    write(9944,*) ''
    write(9944,'(a)') '---------- shared atoms ---------'
    do ii=1,natom 
      write(6,'(a,i3,a)',advance='no')'cluster: ',ii,' share atoms: '
      write(9944,'(a,i3,a)',advance='no')'cluster: ',ii,' share atoms: '
      if ( nshare(ii)/=0) then 
        do jj=1,nshare(ii)
          write(6,'(i4)',advance='no')    ishare(jj,ii) 
          write(9944,'(i4)',advance='no') ishare(jj,ii) 
        enddo
      else
        write(6,'(a)',advance='no') ' none'
        write(9944,'(a)',advance='no') ' none'
      endif 
      write(6,*)''
      write(9944,*)''
    enddo
    write(9944,*)''
    write(9944,*)''
  endif 


  close(file_unit)


end subroutine load_atom_info



