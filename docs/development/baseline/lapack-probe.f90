program probe
 implicit none
 integer, parameter :: n=16
 complex(8) :: a(n,n),z(n,n),work(1024)
 real(8) :: w(n),rwork(4096)
 integer :: isuppz(2*n),iwork(4096),m,info,i
 a=(0d0,0d0)
 do i=1,n
  a(i,i)=dble(i)
 end do
 call zheevr('V','A','U',n,a,n,0d0,0d0,1,n,0d0,m,w,z,n,isuppz, &
             work,1024,rwork,4096,iwork,4096,info)
 print *, 'info=',info,'first=',w(1),'last=',w(n)
 if(info/=0 .or. abs(w(1)-1d0)>1d-12 .or. abs(w(n)-16d0)>1d-12) stop 1
end program probe
