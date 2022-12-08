	module module_init
        implicit none
	integer, parameter::n=505,m=(1.5*n)  !!!-1 for odd no eg- m=(0.7*n)-1
	integer ::iter,iter_N,niter,iter_p
       ! real*8, parameter::res=1.0d-3,rest=1.0e-3
	real*8 :: dx,dy,dt,h,x(n),y(m)
	real*8 :: Residue,residold,xtol
	real*8, dimension(m,n):: xo,xp,SF,SFv,SF2,SFv2
	real*8, dimension(m-2,n-2):: psi,omg,xo0,xp0,resid,resid_psi
	real*8 :: cnv,xmu,xxi,xnu,zz,Re,k1,k2
	!double precision :: cnv,xmu,xxi,xnu,Am,zz,shi_z,coeff,pi,Residue
        
        ! Global arrays
	real*8, dimension(m-2,n-2):: ff0,gg0
	!real*8 :: Ax,Bx,Cx,Ay,By,Cy

	contains
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine create_arrays
	integer :: i,j
	real*8 ::x0,xn,alpha,Lr,U0

!	x0=5.001; xn=6.0d0;
	x0=0.001; xn=1.0d0;
 	dx=(xn-x0)/dble(n-1);
	dy=dx;
 	do i=1,n
   	  x(i)= 1.0*x0+dble(i-1)*dx ;   
 	enddo
	!h=1;
 	h=dble(m-1)*(dy/2);
 	do i=1,m
 	    y(i)= -h+dble(i-1)*dy ;   
 	enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  !	dt=alpha*dx*dy;
!	Residue = dt*2.0d0*(1.0d-5/1.0d-1)
!	xtol=Residue
       	alpha=1/4.50d0;
  !     	alpha=1.6d0;
  	dt=alpha*dx*dy;
  !	dt=5.0e-5
	!Residue = dt*(1.0d-6) ! xo= u_d/lr[U_0/L_r]
	xtol=1.0e-5
!	 Lr= 1.0d-1; U0 = 1.0d-5;
!	 cnv=0.0;  xmu = 1.0d-3 ;  xxi=1.0d-2; xnu=1.0d+1;!% xnu=1.0d-6;
	 cnv=1.0;  xmu = 3.0d-6 ;  xxi=1.0d-4; xnu=1.0d-3;!% xnu=1.0d-6;
	 Re = cnv/xmu;  k1 = (xxi+xnu)/xmu;  k2 = xxi/xmu;  
!		Re = Re*Lr*U0;
!		k1= k1*Lr**2 
!		k2= k2*Lr**2; 
 !	k1=0.0;k2=0.0;
!C     Time step for iteration and initial guess

	psi=1.d-13;
	omg=1.d-13;
	!SF=423.0e-2;
	end subroutine create_arrays
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	subroutine create_source_cosine	
	!!subroutine create_source
	!!integer :: i,j
	!!real*8 ::Am,pi,shi_z,coeff
	!!real*8 ::xs(n)
      !!do  j=1,m
     !! do  i=1,n		
	!!	Am = -1.0d-0 ;   pi = 3.14159 ; zz=pi;
	!!	xs(i)= zz*(x(i)-x(1))/(x(n)-x(1))
  !!              SFv(j,i) = Am*cos(xs(i)); 
!                SF(j,i) =  ((zz*x(1))/(x(n)-x(1)))*Am*BESSEL_JN(1,xs(i));  
 !!               SF(j,i) =  -(zz/(x(n)-x(1)))*Am*sin(xs(i));  

!!	end do
!!	end do
!!	end subroutine create_source
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	subroutine create_source
	integer :: i,j
	real*8 ::Am,pi,shi_z,coeff,z_1,z_n,dzz
	real*8 ::xs(n)
      do  j=1,m
      do  i=1,n		
		Am = -1.0d-0 ;   pi = 3.14159 ;
		z_1 = 1.0;   ! x(1) 
 		z_n = 6.0;  
	! J_0 roots.. zz1 = 2.4048;  ! zz2=5.5201 ; zz3= 8.6537 ; 
	! J_1 roots.. zz1 = 3.8317;  zz2= 7.0156 ; zz3=10.1735 ;
! %% case for real cartesian approch...........
		dzz = (z_n-z_1)/dble(n-1);
		xs(i) = z_1+ dble(i-1)*dzz;
		!xs(i) = z_1+ (x(i)-x(1)/dx)*(z_n-z_1)/dble(n-1)
                SFv(j,i) = Am*BESSEL_JN(1,xs(i)); ! zz = 2.4048; 
                SF(j,i) = Am*((z_n-z_1)/(x(n)-x(1)))*(BESSEL_JN(0,xs(i))-BESSEL_JN(2,xs(i)))/2.0;  



! %% case for cartesian approch...........
!!!		xs(i)=zz*(x(i)-x_1)/(x(n)-x_1)
 !!!               SFv(j,i) = Am*BESSEL_JN(0,xs(i)); ! zz = 2.4048; 
!                SF(j,i) =  ((zz*x(1))/(x(n)-x(1)))*Am*BESSEL_JN(1,xs(i));  
 !!!              SF(j,i) =  ((zz)/(x(n)-x_1))*Am*BESSEL_JN(1,xs(i));  

 !! %% case for !vertical streaming

 	!       SFv(j,i) = Am*BESSEL_JN(0,zz*(x(i))/x(n));             ! zz = 2.4048; 
	!       SF(j,i) =  Am*(zz/x(n))*BESSEL_JN(1,zz*(x(i))/x(n));  !vertical streaming

 !! paper2 %% case for  +half, +1, +-2   ! zz = 2.4048; +h ! zz = 3.8317; +1 
   	!  	    shi_z= cos(pi*y(j)/(2.0d0*h)) ;	   !zz  = 7.0156 ; +-2  horizontal
     	!	   coeff = ((x(n)/zz)*(pi/(2*h))**2 +1);
	!        SF(j,i) = shi_z*coeff*Am*BESSEL_JN(1,zz*(x(i))/x(n)); !!

 !! %% case for 2 vortex  ++horizontal    		!! zz =7.0156 ;  
	!  	   shi_z= (cos(pi*y(j)/(2.0d0*h))) ;
     	!	  coeff = ((x(n)/zz)*(pi/(2*h))**2 +1);
	!       SF(j,i) = 1.0*shi_z*coeff*Am*abs(BESSEL_JN(1,zz*(x(i))/x(n))); !! 


 !! %% case for 2 vortex  ++vertical	!! zz = 3.8317; 
	!	   shi_z= abs(sin(pi*y(j)/(1.0d0*h))) ;
     	!	  coeff = ((x(n)/zz)*(pi/(1*h))**2 +1);
	 !       SF(j,i) = 1.0*shi_z*coeff*Am*BESSEL_JN(1,zz*(x(i))/x(n)); !! 		
			

!! %% case for 4 vortex  h+- v+-          !! zz =7.0156; 
	 !         shi_z= sin(pi*y(j)/(1.0d0*h)) ;
     	!	  coeff = ((x(n)/zz)*(pi/(1*h))**2 +1);
	 !      SF(j,i) = 1.0*shi_z*coeff*Am*BESSEL_JN(1,zz*(x(i))/x(n)); !!

!! %% case for 4 vortex  h++ v++          !! zz =7.0156; 
	!          shi_z= abs(sin(pi*y(j)/(1.0d0*h))) ;
     	!	  coeff = ((x(n)/zz)*(pi/(1*h))**2 +1);
	 !      SF(j,i) = 1.0*shi_z*coeff*Am*abs(BESSEL_JN(1,zz*(x(i))/x(n))); !!

!! %% case for 4 vortex  h+- v+-  h++,v++    !! zz =7.0156; 
	!          shi_z= (sin(pi*y(j)/(1.0d0*h))) ;
     	!	  coeff = ((x(n)/zz)*(pi/(1*h))**2 +1);
	 !      SF(j,i) = 1.0*shi_z*coeff*Am*abs(BESSEL_JN(1,zz*(x(i))/x(n))); !!

	end do
	end do
	end subroutine create_source
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!	subroutine create_source_gaussian	
!	integer :: i,j
!	real*8 ::Am,w,x0,shi_z2,shiR,pi
!      do  j=1,m
 !     do  i=1,n		
!		Am = -1.0d-2;  
!		w=x(10); x0= x(1) ;pi = 3.14159 ;
!	        SFv2(j,i) = Am*exp(-((x(i)-x0)/w)**2);
!	        SF2(j,i) = -(-2.0*(x(i)-x0)/w**2)*SFv2(j,i);

!		shi_z2 = cos(pi*y(j)/(2.0d0*h)) ;
 !      shiR =  (Am/x(i))*(0.5*x0*w*sqrt(pi)*erf((x(i)-x0)/w) &
  !                          - (0.5*w**2)*exp(-((x(i)-x0)/w)**2) ) ;
   !    shiR = shiR - (Am/x(i))*(0.5*x0*w*sqrt(pi)*erf((-x0)/w) &
    !                                   - (0.5*w**2)*exp(-((-x0)/w)**2) ) ;      
     !   SF2(j,i) = shi_z2*( 2*((x(i)-x0)/w**2)*SFv2(j,i) + shiR*(pi/(2*h)**2)  );
             


!	end do
!	end do
!	end subroutine create_source_gaussian


	end module module_init
