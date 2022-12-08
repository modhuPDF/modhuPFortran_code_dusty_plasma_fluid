	program NLsolver
	use module_init
	implicit none
	integer :: i,j
	real*8, dimension(m,n)::xp_old,xo_old
	real*8 :: R_xp_A, R_xo_A, R_xp_B, R_xo_B
	real*8, parameter::cef1=0.0d0,cef2=-1.0d0 
!!!!!   cef1=1, cef2=1 for cylinder,    cef1=0,cef2=-1 for cartesian 
	real*8 :: diffus_t1,diffus_t2,diffus_t3,nonlin_p1,nonlin_p2,nonlin_p3,nonlin_p4
	real*8,dimension(m,n) :: n_fric,i_coll,diffus_t,nonlin_t,sum_t
	call create_arrays

	call create_source
	!call create_source_gaussian	
	!SFv = SFv2 ; SF = SF2;

	call update_psi
	call update_omg
	call update_boundary


	   niter=900000000
 	do iter = 1,niter


 	  do j=2,m-1
  	  do i=2,n-1
         xp_old(j,i)=xp(j,i); 
         xo_old(j,i)=xo(j,i);
   	 enddo
  	 enddo
! 	Solve for psi_xp
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	call solvex_xp
!	ff0 is known	
	call solvey_xp
!	psi is known
	!xp(2:m-1,2:n-1)=psi;
!	Replace old psi with new psi
	call update_psi
!	Replace old omega boundary with new omega boundary using new psi	
	call update_boundary
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
! solve for omg_xo
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@	
	call solvex_xo
!	gg0 is known	
	call solvey_xo
!	omg is known
!	Replace old omega with new omega
	call update_omg
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	!call check_convergence_me
	call check_convergence_sir


	enddo !!! for iteration   		!!!!!!!!!!!!!!

!@@@@@@@@@@@@@@
	
	contains
!000000000000000

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	subroutine check_convergence_sir
!	write(*,*)'Executing iteration no. ',iter
!	resid=xo-xo0
     if (iter.gt.1)then

!	resid=xp-xp0
	resid=(xo0-omg)/minval(omg)
	resid_psi=(xp0-psi)/minval(psi)

	do iter_p =1,900000000,100

	if( iter==iter_p) then
	write(*,*)iter,maxloc(resid),maxval(resid_psi),maxval(resid)
	endif
	end do

!	if (residold .gt. maxval(resid))then
!	  if (abs(maxval(resid)).lt.xtol)then
	  if ((maxval(abs(resid_psi)).lt.xtol).AND.(maxval(abs(resid)).lt.xtol))then

	call save_output
	write(*,*)''
	write(*,*)'used tolerance',xtol
	write(*,*)''
	write(*,*)'SF=', (SF(101:101,101:103))
	!write(*,*) 'xmu=',xmu
	!write(*,*) 'xnu=',xnu
	!write(*,*) 'xxi=',xxi
	!write(*,*) 'h=',h,n,m,dx
	write(*,*) 'zz=',zz
	write(*,*) 'x=',x(1),x(n-10:n-9),x(n)
	write(*,*) 'y=',y(1),y(m-10:m-9),y(m)
	write(*,*) 'xp=',(xp(101:101,101:103))
	  stop 'Convergence achieved_sir'

	  endif
!	else
!	write(*,*)'residold < resid = ',residold,'<',maxval(resid),maxloc(resid)
!	stop 'Diverging'
 !    endif
     


     endif
	xo0=omg
	xp0=psi
	residold=maxval(resid)
	end subroutine check_convergence_sir
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_convergence_me
	real*8 :: R_xp_A, R_xo_A, R_xp_B, R_xo_B
	R_xp_A = 0.0;
	R_xo_A = 0.0;
 	   do j=2,m-1
 	   do i=2,n-1
          R_xp_B = abs(xp(j,i)-xp_old(j,i)); 
          R_xo_B = abs(xo(j,i)-xo_old(j,i));
          R_xp_A = max(R_xp_B,R_xp_A);
          R_xo_A = max(R_xo_B,R_xo_A) ;
  	  enddo
	  enddo

       write(*,*) 'iter== ', iter
       write(*,*)  R_xp_A, R_xo_A
!     if((R_xp_A.lt.Residue).AND.(R_xo_A.lt.Residue)) then
  !    if((R_xp_A.lt.xtol).AND.(R_xo_A.lt.xtol)) then
	if(iter==10) then

	call save_output
	write(*,*)''
	write(*,*)'used tolerance',xtol
	write(*,*)''
	write(*,*) (SF(101:101,101:103))
	write(*,*) 'xmu=',xmu
	write(*,*) 'xnu=',xnu
	write(*,*) 'xxi=',xxi
	write(*,*) 'zz=',zz
	write(*,*) 'h=',h,n,m,dx
	write(*,*) 'x=',x(1),x(n-10:n-9),x(n)
	write(*,*) 'y=',y(1),y(m-10:m-9),y(m)
	write(*,*) (xp(101:101,101:103))
	stop 'Convergence achieved_me'
      endif

end subroutine check_convergence_me
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	subroutine save_output
!00000000000000000000000000
!      c     Output to a file   


!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	  do j=2,m-1
 ! 	  do i=2,n-1
!%%%%%%%%%%%%%%%%%%  11
!	diffus_t1 = (xo(j,i+1)-2.0d0*xo(j,i)+xo(j,i-1))/(dx**2)
!	diffus_t2 = (xo(j,i+1)-xo(j,i-1))/(x(i)*2.0d0*dx) -(xo(j,i)/x(i)**2)
!	diffus_t3 = (xo(j+1,i)-2.0d0*xo(j,i)+xo(j-1,i))/(dy**2)
   	 
!	  diffus_t(j,i) = diffus_t1 + diffus_t2 + diffus_t3
!%%%%%%%%%%%%%%%%%%  22
!	nonlin_p1 = -(xp(j+1,i)-xp(j-1,i))/(2.0d0*dy) 
!	nonlin_p2 =  (xo(j,i+1)-xo(j,i-1))/(2.0d0*dx)
!	nonlin_p3 =  (xp(j,i+1)-xp(j,i-1))/(2.0d0*dx) + xp(j,i)/x(i)
!	nonlin_p4 =  (xo(j+1,i)-xo(j-1,i))/(2.0d0*dy)
	
!	  nonlin_t(j,i) = -Re*(nonlin_p1*nonlin_p2 + nonlin_p3*nonlin_p4) ;
!%%%%%%%%%%%%%%%%%%  33
 !   	n_fric(j,i)= -k1*xo(j,i);
!%%%%%%%%%%%%%%%%%%  44
  !  	i_coll(j,i)= k2*SF(j,i);

!	sum_t(j,i)= diffus_t(j,i) + n_fric(j,i) + i_coll(j,i) + nonlin_t(j,i)
!	enddo
 ! 	 enddo
!	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








	  open(11,file='Xarray1b_r_15_x_505cy_001_cart.dat')
	  open(12,file='Yarray1b_r_15_x_505cy_001_cart.dat')
	open(13,file='Source_v1b_r_15_x_505cy_001_s1_cart.dat')
	open(14,file='Source_w1b_r_15_x_505cy_001_s1_cart.dat')

!	open(15,file='diffusion_in.dat')
!	open(16,file='convection_in.dat')
!	open(17,file='ion_source_in.dat')
!	open(18,file='nutral_coll_in.dat')
!	open(19,file='sum_terms_in.dat') 

	open(20,file='Output_par1_mu_3e6_xi_4_nu_3_r_15_x_505cy_toll5_s1_cart.dat')
	open(21,file='Output_psi1_mu_3e6_xi_4_nu_3_r_15_x_505cy_toll5_s1_cart.dat')
	open(22,file='Output_omg1_mu_3e6_xi_4_nu_3_r_15_x_505cy_toll5_s1_cart.dat')

	write(11,*)(x(i),i=1,n)
	write(12,*)(y(i),i=1,m)
	do i=1,m
	write(13,*)(SFv(i,j),j=1,n)
	enddo
	do i=1,m
	write(14,*)(SF(i,j),j=1,n)
	enddo
!	sum_t(j,i)= diffus_t(j,i) - n_fric(j,i) + i_coll(j,i)   - nonlin_t(j,i)
!	sum_t(j,i)= diffus_t(j,i) - n_fric(j,i) + i_coll(j,i)   - nonlin_t(j,i)
	do j=2,m-1
	write(15,*)(diffus_t(j,i),i=2,n-1)
	enddo
	do j=2,m-1
	write(16,*)(nonlin_t(j,i),i=2,n-1)
	enddo
	do j=2,m-1
	write(17,*)(i_coll(j,i),i=2,n-1)
	enddo
	do j=2,m-1
	write(18,*)(n_fric(j,i),i=2,n-1)
	enddo
	do j=2,m-1
	write(19,*)(sum_t(j,i),i=2,n-1)
	enddo

	!write(17,*) n,h,xmu,xxi,xnu,Residue
	write(20,*) dx,dt,xmu,xxi,xnu,Residue
	do i=1,m
	write(21,*)(xp(i,j),j=1,n)
	enddo
	do i=1,m
	write(22,*)(xo(i,j),j=1,n)
	enddo

	do i=2,m-1
	do j=2,n-1
!	if(i.gt.1.and.i.lt.m)then
!	xo(i,j)=xo0(i-1,j-1)
!	endif
!	if(j.gt.1.and.j.lt.n)then
	xo(i,j)=xo0(i-1,j-1)
!	endif
	enddo
	enddo

!	do i=1,m
!	write(16,*)(xo(i,j),j=1,n)
!	enddo


	close(11)
	close(12)
	close(13)
	close(14)
	close(15)
	close(16)
	close(17)
	close(18)
	close(19)
	close(20)
	close(21)
	close(22)
	close(23)
       write(*,*) 'Core not dumped' 
	end subroutine save_output
!00000000000000000000000000

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine update_psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=2,m-1
	do j=2,n-1
	xp(i,j)=psi(i-1,j-1)
	enddo
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine update_psi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine update_omg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	do i=2,m-1
	do j=2,n-1
	xo(i,j)=omg(i-1,j-1)
	enddo
	enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine update_omg
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine update_boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	do i=1,n
        xp(1,i)=0.0; xp(m,i)=0.0;
        enddo
        do i=1,m
        xp(i,1)=0.0;  xp(i,n)=0.0;
        enddo
        do i=1,n
        !xo(1,i)=0.0;
        xo(1,i)=(-0.5/dy**2)*(8.0*xp(2,i)-xp(3,i))
        !xo(m,i)=0.0;
        xo(m,i)=(-0.5/dy**2)*(8.0*xp(m-1,i)-xp(m-2,i))
        enddo
        do i=2,m-1
        !xo(i,1)=0.0;
        !xo(i,1)=(-0.5/dx**2)*(8.0*xp(i,2)-xp(i,3))-3.d0/dx
        xo(i,1)=(-0.5/dx**2)*(8.0*xp(i,2)-xp(i,3))
        !xo(i,n)=0.0;
        xo(i,n)=(-0.5/dx**2)*(8.0*xp(i,n-1)-xp(i,n-2))
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine update_boundary
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!1111111111111111111
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine solvex_xp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer ::is,i,j,k,ii
	real*8	::am(n-2,n),asm(n-2,n-2),au(n-2,n-2),su(n)
	real*8	::xpdpy(m,n),xpdpydpx(m,n),sm(n),ssm(n),ff(n)
	real*8  :: Ax(n),Bx(n),Cx(n),Ay(n),By(n),Cy(n)
!	function ff0=solvex(xp,xo)
	do i=1,n
	Ax(i)=-dt/dx**2+cef1*(dt/(2.0d0*x(i)*dx))
	Bx(i)=1+2.0*dt/dx**2
	Cx(i)=-dt/dx**2-cef1*(dt/(2.0d0*x(i)*dx))

	Ay(i)=-dt/dy**2
	By(i)=1+2.0*dt/dy**2
	Cy(i)=-dt/dy**2
	end do
!	k=2;
	do k=2,m-1
!	write(*,*)'Solving for ky =', k
	do i=2,n-1
	is=i-1
	am(i,i-1)=Ax(i)
	if (i-1-1.gt.0)then
	asm(is,i-1-1)=Ax(i)
	endif
	am(i,i+0)=Bx(i)
	asm(is,i+0-1)=Bx(i)
	am(i,i+1)=Cx(i)
	if (i+1-1.lt.n-1)then
	asm(is,i+1-1)=Cx(i)
	endif
	enddo

	xpdpy(k,1)=0
	xpdpy(k,n)=0

	do j=2,n-1
	xpdpy(k,j)=Ay(j)*xp(k-1,j)+(By(j)-1.0)*xp(k,j)+Cy(j)*xp(k+1,j)
	enddo
	do j=2,n-1
	xpdpydpx(k,j)=Ax(j)*xpdpy(k,j-1)+(Bx(j)-1.0)*xpdpy(k,j)+Cx(j)*xpdpy(k,j+1)
	enddo

	do j=2,n-1
	   sm(j)=xp(k,j)+dt*xo(k,j)+xpdpydpx(k,j)- cef1*dt*xp(k,j)/(x(j))**2.
	ssm(j-1)=xp(k,j)+dt*xo(k,j)+xpdpydpx(k,j)- cef1*dt*xp(k,j)/(x(j))**2.
!	sm(j)=xp(k,j)+dt*xmu*xo(k,j)+xpdpydpx(k,j)- dt*xp(k,j)/(x(j))**2.
!	ssm(j-1)=xp(k,j)+dt*xmu*xo(k,j)+xpdpydpx(k,j)- dt*xp(k,j)/(x(j))**2.
	enddo

	ff(1)=0.0
	ff(n)=0.0

	sm(2)=sm(2)-Ax(1)*ff(1)
	ssm(2-1)=ssm(2-1)-Ax(1)*ff(1)
	ssm(n-2)=ssm(n-2)-Cx(n)*ff(n)

!Matrix am and asm are tridiaginal matrices, They need to be upper 
!traingularized for solution using back substitution
!	l-u decompose
	
	do i=2,n-2 !i=1s()
	au(i,i)=asm(i,i)-(asm(i-1,i)/asm(i-1,i-1))*asm(i,i-1)
				 !b_i-(c_i-1/b_i-1)*a_i
	su(i)=ssm(i)-(ssm(i-1)/asm(i-1,i-1))*asm(i,i-1)
	asm(i,i)=au(i,i)
	ssm(i)=su(i)
	asm(i,i-1)=0.0*asm(i,i-1) !! Ai* must update to zero 
	enddo

	ff(n-2)=ssm(n-2)/asm(n-2,n-2)
	ff0(k-1,n-2)=ff(n-2)
!	for i=n-3:-1:1;
	do i=n-3,1,-1
	ff(i)=ssm(i)/asm(i,i)
!	ff(i)=ff(i)-sum(asm(i,i+1:n-2).*ff(i+1:n-2))/asm(i,i);
	ff(i)=ff(i)-(asm(i,i+1)*ff(i+1))/asm(i,i)
	ff0(k-1,i)=ff(i)
	enddo

	enddo ! for k
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine solvex_xp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!2222222222
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine solvey_xp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer ::is,i,j,k,ii
	real*8	::am(m-2,m),asm(m-2,m-2),au(m-2,m-2),su(m)
	real*8	::ssm(m),ff(m)
	real*8 :: Ax(m),Bx(m),Cx(m),Ay(m),By(m),Cy(m)
!	function ff0=solvex(xp,xo)
	
	do k=2,n-1
	!write(*,*)'Solving for kx =', k
	!do k=1,n
	do i=2,m-1
	Ax(i)=-dt/dx**2+cef1*(dt/(2.0d0*x(k)*dx))
	Bx(i)=1+2.0*dt/dx**2
	Cx(i)=-dt/dx**2-cef1*(dt/(2.0d0*x(k)*dx))

	Ay(i)=-dt/dy**2
	By(i)=1+2.0*dt/dy**2
	Cy(i)=-dt/dy**2
	enddo
	!enddo
!	
	do i=2,m-1
	is=i-1
	am(i,i-1)=Ay(i)
	if (i-1-1.gt.0)then
	asm(is,i-1-1)=Ay(i)
	endif
	am(i,i+0)=By(i)
	asm(is,i+0-1)=By(i)
	am(i,i+1)=Cy(i)
	if (i+1-1.lt.m-1)then
	asm(is,i+1-1)=Cy(i)
	endif
	enddo
	
	do i=1,m-2
	ssm(i)= ff0(i,k-1)
	enddo

	do i=2,m-2
	au(i,i)=asm(i,i)-(asm(i-1,i)/asm(i-1,i-1))*asm(i,i-1)
	su(i)=ssm(i)-(ssm(i-1)/asm(i-1,i-1))*asm(i,i-1)
	asm(i,i)=au(i,i)
	ssm(i)=su(i)
	asm(i,i-1)=0.0*asm(i,i-1)
	enddo

	ff(m-2)=ssm(m-2)/asm(m-2,m-2)
!	ff0(k-1,m-2)=ff(m-2)
	psi(m-2,k-1)=ff(m-2)
!	for i=m-3:-1:1;
	do i=m-3,1,-1
	ff(i)=ssm(i)/asm(i,i)
!	ff(i)=ff(i)-sum(asm(i,i+1:n-2).*ff(i+1:n-2))/asm(i,i);
	ff(i)=ff(i)-(asm(i,i+1)*ff(i+1))/asm(i,i)
	psi(i,k-1)=ff(i)
	enddo
	enddo ! kx loop
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine solvey_xp
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!%%%%%%%%%%%%%%%%%%%%%%%@@@@@@@@@@@@@@@@@@@
!33333333
	subroutine solvex_xo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer ::is,i,j,k,ii
	real*8	:: am(n-2,n),asm(n-2,n-2),au(n-2,n-2),su(n)
	real*8	:: xodpy(m,n),xodpydpx(m,n),sm(n),ssm(n),ff(n)
	real*8	:: sx(m),sy(m),Ax(n),Bx(n),Cx(n),Ay(n),By(n),Cy(n)

!	
!	k=2; ..............................................
	do k=2,m-1
	!write(*,*)'Solving for ky =', k
	! % u_x at the left boundary
!!!!!!!!!!!!!!
	!sy(k)=0;
	!Ax(1)=-dt/dx**2+(dt/(2.d0*dx))*((1.0d0/x(1))+Re*sy(k))
	!Bx(1)=1+2.0*dt/dx**2
	!Cx(1)=-dt/dx**2-(dt/(2.d0*dx))*((1.0d0/x(1))+Re*sy(k))
!	!!!! % u_x at the right boundary
	!sy(k)=0;
	!Ax(n)=-dt/dx**2 +(dt/(2.d0*dx))*((1.0d0/x(n))+Re*sy(k))
	!Bx(n)=1+2.0*dt/dx**2
	!Cx(n)=-dt/dx**2 -(dt/(2.d0*dx))*((1.0d0/x(n))+Re*sy(k))
!	 % at the interior
	do i=2,n-1
     	sy(k) = (xp(k+1,i)-xp(k-1,i))/(2.0*dy) ; 

	Ax(i)=-dt/dx**2+(dt/(2.d0*dx))*(cef1*(1.0d0/x(i))+Re*sy(k))*cef2
	Bx(i)=1.d0+2.0*dt/dx**2;
	Cx(i)=-dt/dx**2-(dt/(2.d0*dx))*(cef1*(1.0d0/x(i))+Re*sy(k))*cef2
	end do

	!% u_y at the left boundary  !% perfect slip
 	!sx(k)=(xp(k,2)-xp(k,1))/dx; 
	!% u_y at the left boundary  !% no slip
	sx(k)= 0.0;
	Ay(1)=-dt/dy**2-(dt*Re/(2.d0*dy))*sx(k)*cef2!+xp(k,1)/x(1))
	By(1)=1+2.0*dt/dy**2;
	Cy(1)=-dt/dy**2+(dt*Re/(2.d0*dy))*sx(k)*cef2!+xp(k,1)/x(1))
!%%%%%%%%%
	!% u_y at the right boundary   !% perfect slip
 	!sx(k)=(xp(k,N-1)-xp(k,N))/dx; 
	!% u_y at the right boundary   !% no slip
 	sx(k)=0.0; 
	Ay(n)=-dt/dy**2-(dt*Re/(2.d0*dy))*sx(k)*cef2 !+xp(k,n)/x(n))
	By(n)=1+2.0*dt/dy**2;
	Cy(n)=-dt/dy**2+(dt*Re/(2.d0*dy))*sx(k)*cef2 !+xp(k,n)/x(n))
!%%%%%%%%%
!	!% u_y at the right boundary   !% no slip
!	sx(k)=0.0; 
!	Ay(n)=-dt/dy**2!-(dt*Re/(2.d0*dy))*(sx(k)+xp(k,n)/x(n))
!	By(n)=1+2.0*dt/dy**2;
!	Cy(n)=-dt/dy**2!+(dt*Re/(2.d0*dy))*(sx(k)+xp(k,n)/x(n))
!%%%%%%

!	 % inside the interior
	do i=2,n-1! % 
    	sx(k)=(xp(k,i+1)-xp(k,i-1))/(2.0*dx);

	Ay(i)=-dt/dy**2-(dt*Re/(2.d0*dy))*(sx(k)+cef1*xp(k,i)/x(i))*cef2
	By(i)=1+2.0*dt/dy**2;
	Cy(i)=-dt/dy**2+(dt*Re/(2.d0*dy))*(sx(k)+cef1*xp(k,i)/x(i))*cef2
	end do

	!!11111%for k=2:2;

	do i=2,n-1
	is=i-1
	am(i,i-1)=Ax(i)
	if (i-1-1.gt.0)then
	asm(is,i-1-1)=Ax(i)
	endif
	am(i,i+0)=Bx(i)
	asm(is,i+0-1)=Bx(i)
	am(i,i+1)=Cx(i)
	if (i+1-1.lt.n-1)then
	asm(is,i+1-1)=Cx(i)
	endif
	enddo

	!xpdpy(k,1)=0
	!xpdpy(k,n)=0  % 1 to n to be calculated

	do j=1,n
	xodpy(k,j)=Ay(j)*xo(k-1,j)+(By(j)-1.0)*xo(k,j)+Cy(j)*xo(k+1,j);

	enddo
	do j=2,n-1
	xodpydpx(k,j)=Ax(j)*xodpy(k,j-1)+(Bx(j)-1.0)*xodpy(k,j)+Cx(j)*xodpy(k,j+1);
	enddo

	do j=2,n-1
   	   sm(j)=xo(k,j)+xodpydpx(k,j)-k1*dt*xo(k,j) &
			 + k2*dt*SF(k,j)- cef1*dt*xo(k,j)/(x(j))**2;   
  	ssm(j-1)=xo(k,j)+xodpydpx(k,j)-k1*dt*xo(k,j)  &
			 + k2*dt*SF(k,j)- cef1*dt*xo(k,j)/(x(j))**2;
	
	enddo

!	ff(1)=0.0
!	ff(n)=0.0
	ff(1)=(Ay(1)*xo(k-1,1)+By(1)*xo(k,1)+Cy(1)*xo(k+1,1));
	sm(2)=sm(2)-Ax(2)*ff(1);
	ssm(2-1)=ssm(2-1)-Ax(2)*ff(1)

	ff(n)=(Ay(n)*xo(k-1,n)+By(n)*xo(k,n)+Cy(n)*xo(k+1,n));
	sm(n-1)=sm(n-1)-Cx(n-1)*ff(n);
	ssm(n-2)=ssm(n-2)-Cx(n-1)*ff(n);

	do i=2,n-2
	au(i,i)=asm(i,i)-(asm(i-1,i)/asm(i-1,i-1))*asm(i,i-1)
				 !b_i-(c_i-1/b_i-1)*a_i
	su(i)=ssm(i)-(ssm(i-1)/asm(i-1,i-1))*asm(i,i-1)
	asm(i,i)=au(i,i)
	ssm(i)=su(i)
	asm(i,i-1)=0.0*asm(i,i-1) !! Ai* must update to zero 
	enddo

	ff(n-2)=ssm(n-2)/asm(n-2,n-2)
	gg0(k-1,n-2)=ff(n-2)
!	for i=n-3:-1:1;
	do i=n-3,1,-1
	ff(i)=ssm(i)/asm(i,i)
!	ff(i)=ff(i)-sum(asm(i,i+1:n-2).*ff(i+1:n-2))/asm(i,i);
	ff(i)=ff(i)-(asm(i,i+1)*ff(i+1))/asm(i,i)
	gg0(k-1,i)=ff(i)
	enddo

	enddo ! for k

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine solvex_xo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!44444444
	subroutine solvey_xo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	integer ::is,i,j,k,ii
	real*8	::am(m-2,m),asm(m-2,m-2),au(m-2,m-2),su(m)
	real*8	::sm(m),ssm(m),ff(m)
	!real*8	:: sx(m),sy(m),Ax(n),Bx(n),Cx(n),Ay(n),By(n),Cy(n)
	real*8	:: sx(n),sy(n),Ax(m),Bx(m),Cx(m),Ay(m),By(m),Cy(m)
!	function ff0=solvex(xp,xo)

!	k=2; ..............................................
	do k=2,n-1
	!write(*,*)'Solving for ky =', k
	! % u_x at the lower boundary
	!sy(k)=0;
	!Ax(1)=-dt/dx**2 +(dt/(2.d0*dx))*((1.0d0/x(k))+Re*sy(k))
	!Bx(1)=1+2.0*dt/dx**2
	!Cx(1)=-dt/dx**2 -(dt/(2.d0*dx))*((1.0d0/x(k))+Re*sy(k))
!	 % u_x at the upper boundary
	!sy(k)=0;
	!Ax(m)=-dt/dx**2 +(dt/(2.d0*dx))*((1.0d0/x(k))+Re*sy(k))
	!Bx(m)=1+2.0*dt/dx**2
	!Cx(m)=-dt/dx**2 -(dt/(2.d0*dx))*((1.0d0/x(k))+Re*sy(k))
!	 % at the interior
	!do i=2,m-1
     !	sy(k) = (xp(i+1,k)-xp(i-1,k))/(2.0*dy) ; 

	!Ax(i)=-dt/dx**2+(dt/(2.d0*dx))*((1.0d0/x(k))+Re*sy(k))
!	Bx(i)=1+2.0*dt/dx**2;
!	Cx(i)=-dt/dx**2-(dt/(2.d0*dx))*((1.0d0/x(k))+Re*sy(k))
!	end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!% u_y at the lower boundary/confinement
 	  	! sx(k)=(xp(k,2)-xp(k,1))/dx;
!	sx(k)=0.0;
!	Ay(1)=-dt/dy**2!+dt*Re*sx(k)/(2*dy) ;
!	By(1)=1+2.0*dt/dy**2;
!	Cy(1)=-dt/dy**2!-dt*Re*sx(k)/(2*dy);
	!% u_y at the upper boundary
!	sx(k)=0.0; !% no slip
!	Ay(m)=-dt/dy**2!+dt*Re*sx(k)/(2*dy) ;
!	By(m)=1+2.0*dt/dy**2;
!	Cy(m)=-dt/dy**2!-dt*Re*sx(k)/(2*dy);
!	 % at the interior
	do i=2,m-1! % 
    	sx(k)=(xp(i,k+1)-xp(i,k-1))/(2.0*dx);
	Ay(i)=-dt/dy**2-(dt*Re/(2.d0*dy))*(sx(k)+cef1*xp(i,k)/x(k))*cef2
	By(i)=1+2.0*dt/dy**2;
	Cy(i)=-dt/dy**2+(dt*Re/(2.d0*dy))*(sx(k)+cef1*xp(i,k)/x(k))*cef2
	end do
	!!11111%for k=2:2;
	do i=2,m-1
	is=i-1
	am(i,i-1)=Ay(i)
	if (i-1-1.gt.0)then
	asm(is,i-1-1)=Ay(i)
	endif
	am(i,i+0)=By(i)
	asm(is,i+0-1)=By(i)
	am(i,i+1)=Cy(i)
	if (i+1-1.lt.m-1)then
	asm(is,i+1-1)=Cy(i)
	endif
	enddo
!%% source from gg0
		do i=1,m-2
		ssm(i)= gg0(i,k-1)
		enddo
	ff(1)=xo(1,k);
	ff(m)=xo(m,k);
	ssm(1)=ssm(1)-Ay(2)*ff(1);
	ssm(m-2)=ssm(m-2)-Cy(m-1)*ff(m);


	do i=2,m-2
	au(i,i)=asm(i,i)-(asm(i-1,i)/asm(i-1,i-1))*asm(i,i-1)
				 !b_i-(c_i-1/b_i-1)*a_i
	su(i)=ssm(i)-(ssm(i-1)/asm(i-1,i-1))*asm(i,i-1)
	asm(i,i)=au(i,i)
	ssm(i)=su(i)
	asm(i,i-1)=0.0*asm(i,i-1) !! Ai* must update to zero 
	enddo

	ff(m-2)=ssm(m-2)/asm(m-2,m-2)
	!gg0(k-1,m-2)=ff(m-2)
	omg(m-2,k-1)=ff(m-2)
!	for i=n-3:-1:1;
	do i=m-3,1,-1
	ff(i)=ssm(i)/asm(i,i)
!	ff(i)=ff(i)-sum(asm(i,i+1:n-2).*ff(i+1:n-2))/asm(i,i);
	ff(i)=ff(i)-(asm(i,i+1)*ff(i+1))/asm(i,i)
	omg(i,k-1)=ff(i)
	enddo

	enddo ! for kx

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	end subroutine solvey_xo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	end program NLsolver

