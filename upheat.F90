!CS 420: Project

program unsteady_heat

use mpi
implicit none
!variable declaration
integer   nrow_blockto, nrow_blockfrom, srow_blockto, srow_blockfrom,temp
integer:: ierr,tempew,tempns, max_steps,flag,displs
integer:: rank,comm2d,north,south,east,west
integer:: nprocs,remx,remy
integer:: status(MPI_STATUS_SIZE)
integer :: Nx,Ny,i,j,it, maxit,startx,endx,starty,endy,Nxtot,Nytot !no .of grid points in the x and y direction
integer:: global_startx, global_starty, global_endx, global_endy
integer,dimension(2)::dims,coords
integer:: sizes(2),subsizes(2),starts(2)
logical,dimension(2)::period
real(kind(0.d0)), dimension(:,:), allocatable :: Tglobal_old, Tglobal_new, A
real(kind(0.d0)), dimension(30000)::tempor
real(kind(0.d0)):: maxtol, local_error,global_error,c,alpha,dt,dx
double precision start_time, end_time
integer  file,etype,filetype
integer (kind=MPI_OFFSET_KIND):: disp
integer, parameter :: count=3
integer, dimension(3) :: buf


call MPI_INIT(ierr)




!MPI communication logistics
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_size(MPI_COMM_WORLD, nprocs, ierr)
!pre-intialize this array giving full freedom to MPI
dims(1)=0
dims(2)=0
call MPI_Dims_create(nprocs,2,dims,ierr)
period(1)= .false.
period(2)= .false.
call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, .false., comm2d, ierr)
call MPI_Cart_coords(comm2d,rank, 2, coords,ierr)
east = MPI_Proc_null
west = MPI_Proc_null
north = MPI_Proc_null
south = MPI_Proc_null
!identify my neighbours
call MPI_Cart_shift(comm2d, 1, -1, east, west, ierr)
call MPI_Cart_shift(comm2d, 0, -1, south,north, ierr)



!Physical problem logistics
!no. of grid points set by user for the square domain
Nxtot=200
Nytot=Nxtot
dx=1.0/(Nxtot-1)
!coefficient of diffusion
alpha=100.0
!time stepping parameters
max_steps=30000
flag=0
maxtol=1.d-4 !in case running for steady state computation
!choose dt such that the time stepping is stable
dt=(dx*dx)/(4*alpha)
!need to compute c
c= (alpha*dt)/(dx*dx)

global_error=10.0

!no. of interior points exclusing boundaries
Nx=Nxtot-2
Ny=Nytot-2
!global temperature across the domain
allocate(Tglobal_old(Nytot,Nxtot))
allocate(Tglobal_new(Nytot,Nxtot))
!intialize with the intitial conditions at t=0
!Boundary conditions
Tglobal_old(1,:)=100.0
Tglobal_old(Nytot,:)=100.0
Tglobal_old(:,1)=100.0
Tglobal_old(:,Nxtot)=100.0
!next time step matrix
Tglobal_new(1,:)=100.0
Tglobal_new(Nytot,:)=100.0
Tglobal_new(:,1)=100.0
Tglobal_new(:,Nxtot)=100.0
!interior points
do i=2,Nytot-1
  do j=2,Nxtot-1
    Tglobal_old(i,j)=0.0
    Tglobal_new(i,j)=0.0
  end do
end do
!now, we get each processor's portion of the global matrix
!this indexing is in terms of interior numbering and disregards boundaries
IF (MOD(Nx,dims(2)).eq. 0) THEN

    nx=Nx/dims(2)

    startx=(coords(2)*nx)+1
    endx= startx+nx-1

ELSE

  remx=MOD(Nx,dims(2))
  nx=(Nx/dims(2))

  startx=(coords(2)*nx)+1
  endx= startx+nx-1

!only the right edge rank will update endx
  if(coords(2).eq.(dims(2)-1)) then
  endx=endx+remx
  end if

END IF
!same logic for y index
IF (MOD(Ny,dims(1)).eq. 0) THEN

    ny=Ny/dims(1)
    starty=(coords(1)*ny)+1
    endy= starty+ny-1

ELSE

  remy=MOD(Ny,dims(1))
  ny=(Ny/dims(1))
  starty=(coords(1)*ny)+1
  endy= starty+ny-1

!only the right edge rank will update endx
  if(coords(1).eq.(dims(1)-1)) then
  endy=endy+remy
  end if

END IF


!now we regard the boundary and find global indices
global_startx=startx+1
global_starty=starty+1
global_endx=endx+1
global_endy=endy+1

start_time=MPI_Wtime()
tempew=(global_endy-global_starty+1)
tempns= global_endx-global_startx+1




!main time stepping loop
do while(flag<max_steps)

!communication phase
!get the information from other processes and send them my information

!communicate with north neighbour
!define a type vector
   call  MPI_Type_vector(tempns,1, Nxtot, MPI_DOUBLE_PRECISION, nrow_blockto, ierr)
   call  MPI_Type_commit (nrow_blockto, ierr)
   call  MPI_SEND(Tglobal_old(global_starty,global_startx),1,nrow_blockto,  north, 0, MPI_COMM_WORLD, ierr)

   call  MPI_Type_vector(tempns,1, Nxtot, MPI_DOUBLE_PRECISION, nrow_blockto, ierr)
   call  MPI_Type_commit (nrow_blockto, ierr)
   call  MPI_RECV(Tglobal_old(global_starty-1,global_startx),1,nrow_blockto, north, 0, MPI_COMM_WORLD, status, ierr)

!communicate with south neighbour
!define a type vector
   call  MPI_Type_vector(tempns,1, Nxtot, MPI_DOUBLE_PRECISION, srow_blockto, ierr)
   call  MPI_Type_commit (srow_blockto, ierr)
   call  MPI_SEND(Tglobal_old(global_endy,global_startx),1,srow_blockto, south, 0, MPI_COMM_WORLD, ierr)

   call  MPI_Type_vector(tempns,1, Nxtot, MPI_DOUBLE_PRECISION, srow_blockto, ierr)
   call  MPI_Type_commit (srow_blockto, ierr)
   call  MPI_RECV(Tglobal_old(global_endy+1,global_startx),1,srow_blockto, south, 0, MPI_COMM_WORLD, status, ierr)




!communicate with east neighbour
   call MPI_SEND(Tglobal_old(global_starty,global_endx),tempew,MPI_DOUBLE_PRECISION, east, 0, MPI_COMM_WORLD, ierr)
   call MPI_RECV(Tglobal_old(global_starty,global_endx+1),tempew,MPI_DOUBLE_PRECISION, east, 0, MPI_COMM_WORLD, status, ierr)

!communicate with west neighbour
   call MPI_SEND(Tglobal_old(global_starty,global_startx),tempew,MPI_DOUBLE_PRECISION, west, 0, MPI_COMM_WORLD, ierr)
   call MPI_RECV(Tglobal_old(global_starty,global_startx-1),tempew,MPI_DOUBLE_PRECISION, west, 0, MPI_COMM_WORLD, status, ierr)

!computation phase
!now that we have updated the new ghost information from other processes, let's march forward
!advancing one time step
!local_error=0.0
   do i=global_starty,global_endy
    do j=global_startx,global_endx
Tglobal_new(i,j)=Tglobal_old(i,j)+(c*Tglobal_old(i-1,j))+(c*Tglobal_old(i+1,j))+(c*Tglobal_old(i,j-1))+(c*Tglobal_old(i,j+1))&
                 -(4*c*Tglobal_old(i,j))
 !local_error=local_error+((Tglobal_new(i,j)-Tglobal_old(i,j))*(Tglobal_new(i,j)-Tglobal_old(i,j)))
     end do
   end do

   !if (rank .eq. 0) THEN
    !print*,Tglobal_new(global_starty,global_startx)
   !end if

    !call MPI_Barrier(MPI_comm_world, ierr)
    !transfer new to old and repeat loop
     do i=global_starty,global_endy
     do j=global_startx,global_endx
       Tglobal_old(i,j)=Tglobal_new(i,j)
     end do
     end do


!uncomment incase you need to track error
    !call MPI_Allreduce(local_error, global_error,1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_Barrier(MPI_comm_world, ierr)

   !global_error=global_error**0.5


!print midpoint teperature
   !if (rank .eq. 3) THEN
    !print*,Tglobal_new(global_starty,global_startx)
   !end if


flag=flag+1

end do

call MPI_Barrier(MPI_comm_world, ierr)
end_time= MPI_Wtime()

if (rank .eq. 0) THEN
 print*,(end_time-start_time)
 !print*,Tglobal_new(global_starty,global_endx)
end if

!use mpi io to print stuff out to file
!first transfer Tglobal_new to a local_matrix
allocate(A(tempew,tempns))

!portion of the data that each process owns now transferred to a local A
do i=1,tempew
  do j=1,tempns
A(i,j)= Tglobal_new(global_starty-1+i,global_startx-1+j)
!print*,A(i,j)
  end do
end do

! commented out portion
!used for file writing and other minor debugging operations
!if (rank .eq. 1) THEN

!  do i=1,tempew
!    do j=1,tempns
!  print*,A(i,j)
!    end do
!  end do

!  open(1, file="data1.dat",status='new')
  !do i=1,30000
  !  do j=1,Nytot

  !write(1,*) tempor(i)
!write(1,*) global_starty,global_startx
  !  end do
  !end do
  !close(1)


!end if

!call mpi_type_create_subarray(2, [5,5], [5,5], [0,0], mpi_order_fortran, &
!mpi_integer, filetype, ierr)

!call mpi_type_commit(filetype,ierr)

!call MPI_Type_create_subarray (2, [Ny, Nx], [tempew, tempns], &
                              !      [0, 0], &
                              !      MPI_ORDER_FORTRAN, MPI_INTEGER, filetype,  ierr)
!  call MPI_Type_commit(filetype, ierr)



  !call mpi_file_open(mpi_comm_world, 'test.txt', &
  !mpi_mode_wronly + mpi_mode_create, &
  !mpi_info_null, file, ierr)


!  call mpi_file_set_view(file, disp, mpi_integer, filetype, 'native', &
!  mpi_info_null, ierr)

 !disp = 0
 !buf(1)=1
 !buf(2)=2
 !buf(3)=3

!call mpi_file_write_all(file, disp, count, mpi_integer, status, ierr)
!CALL MPI_FILE_SEEK(file,disp,&
!                    MPI_SEEK_SET, ierr)
!call MPI_File_write(file, buf, count, MPI_INTEGER, &
!                          MPI_STATUS_IGNORE, ierr)

!call mpi_file_close(file, ierr)


!displs=[0,]
!counts=[]
!now, we need to get all processes' portion of Tglobal_new into rank 0
!call  MPI_Type_vector(tempns*tempew,tempew, Nxtot, MPI_DOUBLE_PRECISION, temp, ierr)
!call  MPI_Type_commit (temp, ierr)


!call MPI_Gatherv(A(1,1),tempew*tempns,MPI_DOUBLE_PRECISION,  Tglobal_new(global_starty,global_startx), counts,&
                    !displs , MPI_DOUBLE_PRECISION, 0,MPI_COMM_WORLD,ierr)

!open(1, file="data1.dat",status='new')
!do i=1,tempew
  !do j=1,tempns

!write(1,*) A(i,j)
!  end do
!end do
!close(1)

call MPI_FINALIZE(ierr)




end program unsteady_heat
