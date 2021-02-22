program wave
implicit none
real*8:: dx,dt,dir,r,x,y,t_start,t_end
integer:: i,j,k,T,M
real*8, allocatable:: u(:,:,:)
character(10):: fileid
character(20)::filename

M=2306
dx=2.d0/M         !spatial step per side
dt=0.2d0/M                !temporal step
r=(dt**2)/(dx**2)
T=int(2.d0/dt)          !number of temporal steps
dir=0.d0            !dirichlet boundary condition

allocate(u(0:M,0:M,0:2))

call cpu_time(t_start)

do i=0,M
  x=-1+i*dx
  do j=0,M
    y=-1+j*dx
    u(i,j,0)=exp(-300*((x-0.4d0)**2+y**2))
    u(i,j,1)=u(i,j,0)
  enddo
enddo

 do k=1,T
  u(0,:,2)=dir
  u(M,:,2)=dir
  u(:,0,2)=dir
  u(:,M,2)=dir

  do i=1,M-1
    do j=1,M-1
      u(i,j,2)=r*(u(i+1,j,1)+u(i-1,j,1)+u(i,j+1,1)+u(i,j-1,1))+(2-4*r)*u(i,j,1)-u(i,j,0)
    enddo
  enddo

  u(:,:,0)=u(:,:,1)
  u(:,:,1)=u(:,:,2)

  ! if ((k==1).or.(k==int(T/3)).or.(k==int(2*T/3)).or.(k==T)) then
   !if (mod(k-1,10)==0) then
   !  write(fileid,'(i0)') k-1
   !  filename='wave'//trim(adjustl(fileid))//'.dat'
   !  open(10,file=trim(filename))
   !  do i=0,M
   !    write(10,*) u(i,:,2)
   !  enddo
   !  print *,'time ',k
   !endif


 enddo
call cpu_time(t_end)
print *,'Serial'
print *,'DImension :',M
print *,'Time      :',t_end-t_start
!do i=0,M
!  print *,u(:,i,2)
!enddo

deallocate(u)
end program wave
