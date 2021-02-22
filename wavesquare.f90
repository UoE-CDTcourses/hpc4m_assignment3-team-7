program wave_squares
implicit none
include 'mpif.h'
real*8:: dx,dt,r,x,y,final_t,t_start,t_end
integer:: i,j,k,T,M,Jproc,start_row,start_col,n
real*8, allocatable:: u(:,:,:),right(:),left(:),sright(:),sleft(:),solt(:,:),sol(:,:)
integer:: comm,rank,nproc,ierr,num
character(10):: fileid
character(20):: filename
real*8:: lrcor,ulcor,urcor,llcor


comm = MPI_COMM_WORLD
call MPI_INIT(ierr)
call MPI_COMM_RANK(comm, rank, ierr)
call MPI_COMM_SIZE(comm, nproc, ierr)

n=int(sqrt(real(nproc)))           !number of squares per side
M=2306                    !number of spatial intervals (points from 0 to nx)
M=M-mod(M-1,n)      !make nx satisfy: (nx-1)/nproc is integer
final_t=1.d0              !final integration time     

dx=2.d0/M                 !spatial step
dt=0.2d0/M              !temporal step
T=int(final_t/dt)         !number of temporal intervals
r=(dt**2)/(dx**2)               !parameter for Euler  

Jproc=(M-1)/n+2
allocate(u(1:Jproc,1:Jproc,0:2))
allocate(right(1:Jproc),left(1:Jproc),sright(1:Jproc),sleft(1:Jproc))
!allocate(solt(0:M,1:nproc*Jproc),sol(0:M,0:M))
u=0.d0
num=0
!% % % % % % % % % % % % % % % % % % % % % % 

call mpi_barrier(comm,ierr)
t_start=mpi_wtime()

start_row=(Jproc-2)*int(rank/n)
start_col=(Jproc-2)*mod(rank,n)

!set initial condition
do i=1,Jproc
  x=-1+(start_row+i-1)*dx
  do j=1,Jproc
    y=-1+(start_col+j-1)*dx
    u(i,j,0)=exp(-40*((x-0.4d0)**2+y**2))
    u(i,j,1)=u(i,j,0)
  enddo
enddo


do k=1,T

!perform Euler method for points in the interior
  do i=2,Jproc-1
    do j=2,Jproc-1
      u(i,j,2)=r*(u(i+1,j,1)+u(i-1,j,1)+u(i,j+1,1)+u(i,j-1,1))+(2-4*r)*u(i,j,1)-u(i,j,0)
    enddo
  enddo


!message swap between columns of processors
  if (mod(rank,n)==0) then
    
    sright=u(2:Jproc-1,Jproc-1,2)
    call mpi_send(sright,Jproc-2,mpi_real8,rank+1,1,comm,ierr)

    call mpi_recv(right,Jproc-2,mpi_real8,rank+1,1,comm,mpi_status_ignore,ierr)
    u(2:Jproc-1,Jproc,2)=right
    
    if (int(rank/n).ne.n-1) then
      
      lrcor=u(Jproc-1,Jproc-1,2)
      call mpi_send(lrcor,1,mpi_real8,rank+n+1,1,comm,ierr)
      
      call mpi_recv(lrcor,1,mpi_real8,rank+n+1,1,comm,mpi_status_ignore,ierr)
      u(Jproc,Jproc,2)=lrcor
      
    endif
    
  elseif (mod(rank,n)==n-1) then

    call mpi_recv(left,Jproc-2,mpi_real8,rank-1,1,comm,mpi_status_ignore,ierr)
    u(2:Jproc-1,1,2)=left
    
    sleft=u(2:Jproc-1,2,2)
    call mpi_send(sleft,Jproc-2,mpi_real8,rank-1,1,comm,ierr)

    if (int(rank/n).ne.0) then
      
      call mpi_recv(ulcor,1,mpi_real8,rank-n-1,1,comm,mpi_status_ignore,ierr)
      u(1,1,2)=ulcor
      
      ulcor=u(2,2,2)
      call mpi_send(ulcor,1,mpi_real8,rank-n-1,1,comm,ierr)
      
    endif
    
  else 

    call mpi_recv(left,Jproc-2,mpi_real8,rank-1,1,comm,mpi_status_ignore,ierr)    
    u(2:Jproc-1,1,2)=left    

    sleft=u(2:Jproc-1,2,2)
    call mpi_send(sleft,Jproc-2,mpi_real8,rank-1,1,comm,ierr)

    sright=u(2:Jproc-1,Jproc-1,2)
    call mpi_send(sright,Jproc-2,mpi_real8,rank+1,1,comm,ierr)
    
    call mpi_recv(right,Jproc-2,mpi_real8,rank+1,1,comm,mpi_status_ignore,ierr)    
    u(2:Jproc-1,Jproc,2)=right

    if (int(rank/n).ne.0) then
      
      call mpi_recv(ulcor,1,mpi_real8,rank-n-1,1,comm,mpi_status_ignore,ierr)
      u(1,1,2)=ulcor
      
      ulcor=u(2,2,2)
      call mpi_send(ulcor,1,mpi_real8,rank-n-1,1,comm,ierr)
    
    endif
    if (int(rank/n).ne.n-1) then

      lrcor=u(Jproc-1,Jproc-1,2)
      call mpi_send(lrcor,1,mpi_real8,rank+n+1,1,comm,ierr)
      
      call mpi_recv(lrcor,1,mpi_real8,rank+n+1,1,comm,mpi_status_ignore,ierr)
      u(Jproc,Jproc,2)=lrcor
      
    endif
    
  endif

!  call mpi_barrier(comm,ierr)
!message swap between rows of processors

  if (int(rank/n)==0) then
    
    sright=u(Jproc-1,2:Jproc-1,2)
    call mpi_send(sright,Jproc-2,mpi_real8,rank+n,1,comm,ierr)

    call mpi_recv(right,Jproc-2,mpi_real8,rank+n,1,comm,mpi_status_ignore,ierr)
    u(Jproc,2:Jproc-1,2)=right
    
    if (mod(rank,n).ne.0) then
      llcor=u(Jproc-1,2,2)
      call mpi_send(llcor,1,mpi_real8,rank+n-1,1,comm,ierr)
      
      call mpi_recv(llcor,1,mpi_real8,rank+n-1,1,comm,mpi_status_ignore,ierr)
      u(Jproc,1,2)=llcor

    endif
    
  elseif (int(rank/n)==n-1) then

    call mpi_recv(left,Jproc-2,mpi_real8,rank-n,1,comm,mpi_status_ignore,ierr)
    u(1,2:Jproc-1,2)=left
    
    sleft=u(2,2:Jproc-1,2)
    call mpi_send(sleft,Jproc-2,mpi_real8,rank-n,1,comm,ierr)

    if (mod(rank,n).ne.n-1) then

      call mpi_recv(urcor,1,mpi_real8,rank-n+1,1,comm,mpi_status_ignore,ierr)
      u(1,Jproc,2)=urcor

      urcor=u(2,Jproc-1,2)
      call mpi_send(urcor,1,mpi_real8,rank-n+1,1,comm,ierr)

    endif
  
  else 

    call mpi_recv(left,Jproc-2,mpi_real8,rank-n,1,comm,mpi_status_ignore,ierr)    
    u(1,2:Jproc-1,2)=left    

    sleft=u(2,2:Jproc-1,2)
    call mpi_send(sleft,Jproc-2,mpi_real8,rank-n,1,comm,ierr)

    sright=u(Jproc-1,2:Jproc-1,2)
    call mpi_send(sright,Jproc-2,mpi_real8,rank+n,1,comm,ierr)
    
    call mpi_recv(right,Jproc-2,mpi_real8,rank+n,1,comm,mpi_status_ignore,ierr)    
    u(Jproc,2:Jproc-1,2)=right

    if (mod(rank,n).ne.n-1) then

      call mpi_recv(urcor,1,mpi_real8,rank-n+1,1,comm,mpi_status_ignore,ierr)
      u(1,Jproc,2)=urcor

      urcor=u(2,Jproc-1,2)
      call mpi_send(urcor,1,mpi_real8,rank-n+1,1,comm,ierr)

    endif
    
    if (mod(rank,n).ne.0) then
      llcor=u(Jproc-1,2,2)
      call mpi_send(llcor,1,mpi_real8,rank+n-1,1,comm,ierr)
      
      call mpi_recv(llcor,1,mpi_real8,rank+n-1,1,comm,mpi_status_ignore,ierr)
      u(Jproc,1,2)=llcor

    endif
  
  endif


!save solution at 0, 1/3, 2/3, 3/3 of final time
!  if ((k==1).or.(k==int(T/3)).or.(k==int(2*T/3)).or.(k==T)) then
!
!    call mpi_gather(u(:,:,2),Jproc*(M+1),mpi_real8,solt,Jproc*(M+1),mpi_real8,0,comm,ierr)
!    if (rank==0) then
!      do i=1,nproc
!        sol(:,(i-1)*(Jproc-2):i*(Jproc-2)-1)=solt(:,(i-1)*Jproc+1:i*Jproc-2)    
!      enddo
!      sol(:,M-1:M)=solt(:,nproc*Jproc-1:nproc*Jproc)
!      
!      write(fileid,'(i0)') num
!      num=num+1
!      filename='wave'//trim(adjustl(fileid))//'.dat'
!      open(10,file=trim(filename))
!      
!      do i=0,M
!        write(10,*) sol(i,:)
!      enddo
!      close(10)
!
!      solt=0.d0
!      sol=0.d0
!    endif
!  endif


!make current solution the initial one for next iteration
  u(:,:,0)=u(:,:,1)
  u(:,:,1)=u(:,:,2)
enddo

!do i=0,nproc-1
!  if (rank==i) then
!    do j=1,Jproc
!      print *,(u(j,k,2),k=1,Jproc)
!    enddo
!  endif
!  call mpi_barrier(comm,ierr)
!  print *
!  print *
!enddo



call mpi_barrier(comm,ierr)
t_end=mpi_wtime()

if (rank==0) then
  print *,'Squares'
  print *,'Processors :',nproc
  print *,'Squares    :',n
  print *,'Dimension  :',M
  print *,'Time       :',t_end-t_start
endif





deallocate(u,right,left,sright,sleft)
call mpi_finalize(ierr)
end program wave_squares
