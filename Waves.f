      program TruncatedStokesTwoModes
c     ----------------------------------------------------------------
c     $Revision: 1.1 $
c     $Date: 2014/7  $
c     ----------------------------------------------------------------
c     Programmer(s): Nick Pizzo  
c     ----------------------------------------------------------------
C      Test for predictor corrector method
c     ----------------------------------------------------------------
c
      implicit none

C      include 'omp_lib.h'
      
      integer ier, globalstrat, maxl, maxlrst,mout,count
      integer nm, neq, ii,niv,niv0
      parameter (nm=1024)
      parameter ( niv= nm)
      parameter (neq = 4*nm)
      parameter (mout=10000)
      integer msbpre
      integer*4 iout(15)
      double precision pp, fnormtol, scsteptol, a0,c2,eps, IL
      double precision rout(2), uu(neq), scale(neq), c
      double precision constr(neq), y(neq), yp(neq), seconds
      double precision t, dt, v2(neq), vv(neq), tout,pi,ainit
      integer  MMax, ipiv(neq),ipiv2(1:2*nm)
      double precision fval(neq), mu,r2n(1:2*nm),r3n(1:2*nm),
     & r4n(1:2*nm) ,fval2(neq),u(neq),ypo(neq),yo(neq),fv2(2*nm),
     * fv(2*nm), Energy
      double precision tres, Ao(1:2*nm,1:2*nm),r2(1:2*nm),err,KE, PE
      integer reserr, No
      integer iter, info       ,l5    
      integer m, j, kk, l, n, l1, l2, l3, i, ll, indsum, indsum2
      integer indsum3, i4, k, k1, k2,k3, k4, k5, k6, l2p, l2pp, Mt
      integer k3a, k3b, kcut1, kcut2, kcut3, nmp, l4
      double precision  po, tcond,  qout(1:nm,-nm:nm), nu, kpr
      double precision  yp2(1:neq), yp3(1:neq),Bmat(1:2*nm,1:2*nm)
      double precision uout(1:neq), Am(1:2*nm,1:2*nm),Cmat(neq,neq),
     & rn(1:2*nm),y4(1:neq), y3(1:neq),uout4(1:neq),yp4(neq)
      complex*16 p(-nm:nm,-nm:nm), V(-nm:nm), resd2a,
     &  Vec(-nm:nm), res2(-nm:nm)
      complex*16 qmat(-nm:nm,-nm:nm,-nm:nm),res(-nm:nm)
      complex*16 sqmat(-nm:nm,-nm:nm,-nm:nm)
      complex*16 sqmatdot(-nm:nm,-nm:nm,-nm:nm)      
      complex*16 sqmatc(-nm:nm,-nm:nm,-nm:nm) ,qo(-nm:nm,-nm:nm)  
      complex*16 y2(neq),restest(1:nm),gamma(-nm:nm,-nm:nm),
     & beta(-nm:nm,-nm:nm), gammad(-nm:nm,-nm:nm) 
      complex*16 Q2(1:2*nm,1:2*nm), Q0, q(-nm:nm,-nm:nm)
      complex*16 Q3(1:2*nm),sbd(-nm:nm),sb(-nm:nm),so,yd,ydd
      complex*16 Q4(1:2*nm), B(-nm:nm), Bd(-nm:nm), Bdd(-nm:nm),
     & Q5(1:nm) 
      complex*16 sqmat1(-nm:nm,-nm:nm,-nm:nm),sbdd(-nm:nm)
      complex*16 sqmat2(-nm:nm,-nm:nm,-nm:nm)       
      complex*16 sdiff(-nm:nm,-nm:nm,-nm:nm)      
      complex*16 qou(-nm:nm,-nm:nm),res4(-nm:nm),betao,CI
      complex*16 work(2*nm), ypout(1:2*nm),res2o(-nm:nm)
      complex*8 swork((2*nm)*(2*nm+1))
      double precision rwork(2*nm), uout2(1:neq),uout3(1:neq),
     & yinitial(2*nm), u5(1:neq),apert
      complex*16 sumdum, sum1, sum2(nm), sum3(nm), qsum,
     &  resd, resd2, scount, pdiff(-nm:nm,-nm:nm)
      complex*16 qmatdiff(-nm:nm,-nm:nm,-nm:nm)
      complex*16 sqmatdiff(-nm:nm,-nm:nm,-nm:nm)
      complex*16 ssum1, ssum2, forcing(nm), qdiff(-nm:nm,-nm:nm)
      complex*16 sum10(nm), sum20(nm)
      double precision RK1(2*nm),Amo(2*nm,2*nm), yn(neq),yn2(neq),
     & RK2(2*nm),RK3(2*nm), RK4(2*nm), r(2*nm),amp, ydiff(neq),
     & Amo2(2*nm,2*nm)
      complex*16 dot, zdotu,
     & smat(-nm:nm,-nm:nm,-nm:nm), alpha   
      CHARACTER(LEN=20) :: fname
C      call omp_set_num_threads ( 4 )      
      do i=1,neq
      y(i) = 0d0
      enddo
      
C      c= 1.0890724934975324d0
C      c= 1.0256832824767739  d0
      c= 1.0890723286321964d0
C     Initial start up
      do i =1,nm
      u(i)= 0d0
      enddo 
      open(unit=11, file='sc_IC_512_0828.txt', status ='unknown')
      read(11,*) u(1:512)
      close(unit = 11)

C      add a perturbation 
      eps = 0.01d0
      u(1)=u(1)+eps

      No = 1
      do i = No, niv, No
      y(i) = u(i/No)
      enddo 
      do i = niv+1,3*niv
      y(i) = 0.0d0
      enddo
      do i = 3*niv+1,neq
      y(i) = y(i-3*niv)*(c/No)*(i-3*niv)
      enddo
      

CC      adding a perturbation 
      y(1) = y(1) 
      y(nm+1) = y(nm+1)
      y(2*nm+1) = y(2*nm+1) 
      y(3*nm+1) = y(3*nm+1)
            
C      niv0=256*2*2
CC      
CC     Restart 
C      do i =1,nm
C      u(i)= 0d0
C      enddo 
C      open(unit=11, file='SW1024b_00060.txt', status ='unknown')
C      read(11,*) u(1:4*niv0)
C      close(unit = 11)      
CC
C      do i = 1,niv0
C      y(i) = u(i)
C      enddo 
C      do i = niv+1,niv+niv0
C      y(i) = u(i-niv+niv0)
C      enddo
C      do i = 2*niv+1,2*niv+niv0
C      y(i) = u(i-2*niv+2*niv0)
C      enddo
C      do i = 3*niv+1,3*niv+niv0
C      y(i) = u(i-3*niv+3*niv0)
C      enddo
C            
C      
C       print*, u        
      u=y
      WRITE(fname,'(A,I5.5,A)') 'SW1024smooth_00001.txt'
       OPEN(14,file=fname,form='formatted')
       WRITE (14,*) (u(i), i = 1,neq)
 25   format (D20.10) 
       close(unit =14)
       
       t=0d0
       dt = 0.001d0
C      
C
      do count=1,mout
C      print*, count

C      tout = (count)*dt
      


C      WRITE(fname,'(A,I5.5,A)') 't.txt'
C       OPEN(42,file=fname,form='formatted')
C       write (42,*) t
C      CLOSE(UNIT=42) 
  
      WRITE(fname,'(A,I5.5,A)') 'SW1024smooth_',count,'.txt'
       OPEN(15,file=fname,form='formatted')
       READ (15,*) uout(1:neq)
       CLOSE(UNIT=15)      
C       print*, u-uout
      
 100   t=t+dt           
       do i = 1,neq
       y(i)=0d0
       yp(i) =0d0
       enddo
       
       do i =1,neq
       y(i) = uout(i)
       enddo
       do i =1,2*nm
       yp(i)=uout(2*nm+i)
       enddo    
       do i =1,2*nm
       ypo(i) = yp(i)
       enddo   
       do i =1,4*nm
       yo(i) = y(i)
       enddo   
C       print*, ypo
C      print*,sqrt(y(1)**2+y(1+nm)**2)  
      CI=dcmplx(0.0d0,1.0d0)
    


      do n=-nm,nm
      so=dcmplx(0d0,0d0)
      do i =-nm,nm
C      do j=-nm,nm
      k1=n-i
      if (k1.eq.0) then
      k1 = 1
      else
      k1 = 0
      endif
      k2=i
      if (k2.eq.0) then
      k2 = 1
      else
      k2=0
      end if
      k3 =abs(n-i)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if

      so = so + 0.25d0*dcmplx(y(2*niv+abs(i)),sign(1,i)*
     & y(3*niv+abs(i)))*(1-k2)*(1-k1)
     & *k3*dcmplx(y(abs(n-i)),sign(1,n-i)*y(niv+abs(n-i)))/
     & ((abs(i)+k2)*(abs(n-i)+k1))
     & *(n-i)*((1-k2)*sign(1,i)-(1-k1)*sign(1,n-i))  
C	  enddo
	  enddo
	  sb(n)= so
	  enddo

C	  
C	  print*, sb-sb2
	  yd=dcmplx(0d0,0d0)
	  do i =1,nm
	  yd=yd-0.5d0*(dcmplx(y(2*niv+i),y(3*niv+i))*
     & dcmplx(y(i),-y(niv+i))
     & +dcmplx(y(2*niv+i),-y(3*niv+i))*dcmplx(y(i),y(niv+i)))/i
      enddo
      
      ydd=dcmplx(0d0,0d0)
      do i =1,nm
	  ydd=ydd+(dcmplx(y(2*niv+i),y(3*niv+i))*
     & dcmplx(y(2*niv+i),-y(3*niv+i)))/i
      enddo
      
C	  print*, yd
C	  print*, ydd
	  
C	  Next, we set up B(n)
	  do n=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else 
	  k1 =0
	  end if  
	  
	   B(n)=dcmplx(0.0d0,1.0d0)/(n+k1)*( sb(n)-
     & (1-k1)*dcmplx(y(2*nm+abs(n)),sign(1,n)* y(3*nm+abs(n)))/
     & (2d0*(abs(n)+k1))-0.5d0*yd*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*(1-k1)-
     &  yd*k1) 
	 
	  enddo
	  
C	   Compute the KE
C	  KE=0d0
C      do n=-nm,nm
C      KE=KE+0.5d0*abs(n)*B(n)*B(-n)
C      enddo 
CC      print*, KE
CC     Compute the PE (only valid for 1 mode!) 
C      PE=0.25d0*(y(1)**2+y(niv+1)**2)
C     & -0.125d0*(y(1)**2+y(niv+1)**2+0.5d0*(y(2)**2+y(niv+2)**2))**2
C     & +0.0625d0*(y(2)**2+y(niv+2)**2)+
C     &1d0*(0.25d0*y(1)**2*y(2)-0.25d0*y(niv+1)**2*y(2)
C     & +0.5d0*y(niv+1)*y(niv+2)*y(1))
CC      print*, PE
C      Energy = 4*(KE+PE)
C      Print*, Energy
      

     
      do n=-nm,nm
      so=dcmplx(0d0,0d0)
      do i =-nm,nm
      k1=n-i
      if (k1.eq.0) then
      k1 = 1
      else
      k1 = 0
      endif
      k2=i
      if (k2.eq.0) then
      k2 = 1
      else
      k2=0
      end if
      k3 =abs(n-i)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if
      k4 = n
      if (k4.eq.0) then
      k4 =0 
      else 
      k4 =1
      end if
      so = so +0.25d0*((1-k1)*
     & dcmplx(y(2*niv+abs(i)),sign(1,i)*y(3*niv+abs(i)))*k3*(1-k2)
     & *dcmplx(y(abs(n-i)+2*niv),sign(1,n-i)*y(3*niv+abs(n-i))))/
     & ((abs(i)+k2)*(abs(n-i)+k1))
     & *(n-i)*(sign(1,i)-sign(1,n-i))*k4
	  enddo
	  sbd(n)=so
	  enddo
C	  print*, sbd
C	  print*, (yp(2*niv+1),yp(3*niv+1))
C	  Next, we set up Bd(n)

	  do n=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else 
	  k1 =0
	  end if
	  Bd(n) = dcmplx(0,1)/(n+k1)*(-0.5d0*yd*
     & (1-k1)*dcmplx(y(abs(n)+2*niv),sign(1,n)*y(3*niv+abs(n))) + 
     & ydd*dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*(1-k1)/2+ 
     & k1*ydd+  sbd(n))
	  enddo
C	  Bd(0)=dcmplx(0d0,0d0)
	  
C	  print*, Bd(1)
C	  -0.5d0*dcmplx(0,1)*(c**2)*dcmplx(y(1),y(niv+1))	  
C	  Then we set up beta(n,m)
	  
	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =abs(n-m)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if
      k4 =m
      if (k4.eq.0) then
      k4 =0
      else
      k4=1
      end if
	  beta(n,m) = dcmplx(0,1d0)/(n+k1)*(0.25d0*k3*(1-k1)*(1-k2)*
     & dcmplx(y(2*niv+abs(n-m)), sign(1,n-m)*y(3*niv+abs(n-m)))
     & /((abs(m)+1-k4)*(abs(n-m)+k2))*m*((1-k2)*sign(1,n-m)-
     & (k4)* sign(1,m))-0.5d0*yd*k2*k4
     & +0.25d0*dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*
     & dcmplx(y(2*niv+abs(m)),-sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k4))*k4+
     & 0.5d0*k1*dcmplx(y(2*niv+abs(m)),sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k4)))
	  enddo
	  enddo 
	  
	  	  
	  
C	  print*, beta
C	  do i =1,nm
C      diff(i)  = beta(-i,i)
C      enddo
C	  print*, beta(0,1)
C	  Then set up gamma(n,m)

	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =n-m
	  if (k3.eq.0) then
	  k3=1
	  else
	  k3=0
	  endif
	  	  k4 =abs(n-m)
      if (k4.gt.nm) then
      k4 =0
      else
      k4=1
      end if
      k5 =m
      if (k5.eq.0) then
      k5 =0
      else
      k5=1
      end if
	  gamma(n,m) = (1-k1)*dcmplx(0,1d0)/(n+k1)*(-0.5d0
     & *dcmplx(k3*1d0,0d0)/(abs(m)+ 1-k5) + k5*(1-k2)*
     &  0.25d0*k4*dcmplx(y(abs(n-m)), sign(1,n-m)*y(niv+abs(n-m)))
     & /((abs(m)+1-k5)
     & *(abs(n-m)+k2))*(n-m)*((k5)*sign(1,m)-(1-k2)*sign(1,n-m))+
     & 0.25d0*dcmplx(y(abs(m)),-sign(1,m)*y(niv+abs(m)))*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))/abs(m+1-k5)
     &*k5*(1-k1)+0.5d0*k1*dcmplx(y(abs(m)),sign(1,m)*y(niv+abs(m)))/
     & abs(m+(1-k5)))
	  enddo
	  enddo 	  

C      print*, gamma
C      print*, dcmplx(0d0,1d0)/1d0*dcmplx(y(1),y(niv+1))
C	   next we compute gammad

	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =n-m
	  if (k3.eq.0) then
	  k3=1
	  else
	  k3=0
	  endif
	  	  k4 =abs(n-m)
      if (k4.gt.nm) then
      k4 =0
      else
      k4=1
      end if
      k5 =m
      if (k5.eq.0) then
      k5 =0
      else
      k5=1
      end if
	  gammad(n,m) = (1-k1)*dcmplx(0,1d0)/(n+k1)*(0.25d0*k5*
     & k4*(1-k3)*dcmplx(y(2*niv+abs(n-m)), 
     & sign(1,n-m)*y(3*niv+abs(n-m)))
     & /((abs(m)+1-k5)
     & *(abs(n-m)+k2))*(n-m)*((k5)*sign(1,m)-(1-k2)*sign(1,n-m))+
     & 0.25d0*dcmplx(y(abs(m)+2*niv),-sign(1,m)*y(3*niv+abs(m)))*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))
     & /abs(m+1-k5)*k5+
     & 0.25d0*dcmplx(y(abs(m)),-sign(1,m)*y(niv+abs(m)))*
     & dcmplx(y(2*niv+abs(n)),sign(1,n)*y(3*niv+abs(n)))
     & /abs(m+1-k5)*k5+
     & 0.5d0*k1*dcmplx(y(2*niv+abs(m)),sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k5)))
	  enddo
	  enddo
C	  print*, gammad   
C	  print*, gammad(1,2) + 0.25d0*dcmplx(y(2*niv+1),-y(3*niv+1))*
C     &dcmplx(0d0,1d0) 


	  do  m = -nm, nm
          do j = -nm, nm
	      kk = m-j
	      l = m-j
	      k1 = m
	      k2 = j
	      k3a = m
	      k3b = j
	      k4 = m-j
	      k5 = j
	      k6 = m-j
	                         
	      if (abs(kk) .gt. nm) then
	      kk = 0
	      else 
	      kk = 1
	      end if
	      if (k1 .eq. 0) then
	      k1 = 0
	      else
	      k1=1
	      end if
	      if (k2.eq.0) then
	      k2 = 0
	      else
	      k2 = 1
	      end if
	      if (k3a .eq. 0) then
	      k3a=1
	      else 
	      k3a=0
	      end if
	      if (k3b .eq. 0) then
	      k3b=1
	      else 
	      k3b=0
	      end if
	      if (k4 .eq. 0) then
	      k4 = 0
	      else
	      k4 = 1
	      end if
	      if (k5 .eq. 0) then
	      k5=2
	      else
	      k5=1
	      end if
	      if (k6 .eq. 0) then
	      k6=2
	      else
	      k6=1
	      end if
	      
	  p(m,j)=(0.5d0*k5*k6*kk*abs(k4*sign(1,m-j)-k2*sign(1,j))*dcmplx
     &(k4*y(abs(m-j))+1-k4,k4*y(nm+abs(m-j))*sign(1,m-j))-
     & 0.5d0*dcmplx(k1
     &*y(abs(m))+(1-k1), k1*y(abs(m)+nm)*sign(1, m))*dcmplx(k2*y(abs(j
     &))+(1-k2),-k2*y(abs(j)+nm)*sign(1,j)))/(DSQRT(1.0d0*(abs(m)+k3a
     & ))*
     & (k3b+abs(j)))
 	 
      enddo  
      enddo	
      
      

CC      Next we set up the leading term coefficient Q

      Mt= 2*nm+1
      alpha = dcmplx(0.25d0,0.0d0)
      betao = dcmplx(0.0d0,0.0d0)
      

C       seconds = omp_get_wtime ( )  

      call ZGEMM('T','N',Mt,Mt,nm,alpha, p(1:nm,-nm:nm),nm,
     &p(-1:-nm:-1,-nm:nm),nm,betao, q ,Mt)
     
C      print*, q


       do m = 1, nm
       ssum1 = dcmplx(0.0d0,0.0d0)
       do n = -nm, nm
        l2=m
        l2p=n
        l2pp=m*n
        l3=(n+m)*n
        l4=(n+m)
        l5=abs(n+m)

         if (l2 .eq. 0) then
          l2 = 1
         else 
          l2 = 0     
         end if 
         if (l2p .eq. 0) then
          l2p = 1
         else 
          l2p = 0     
         end if  
         if (l2pp .eq. 0) then
          l2pp = 1
         else 
          l2pp = 0     
         end if 
         if (l3 .eq. 0) then
          l3 = 0
         else 
          l3 = 1     
         end if 
         if (l4 .eq. 0) then
          l4 = 1
         else 
          l4 = 0     
         end if 
        if (l5 .gt. nm) then
          l5 = 0
         else 
          l5 = 1     
         end if 
    
              
         ssum1 = ssum1 + l5*l3*dcmplx(y(abs(n)),
     &    sign(1,n)*y(abs(n)+nm))*dcmplx(y(abs(n+m)),
     &    -sign(1,n+m)*y(abs(n+m)+nm))/(8d0*(abs(m)+l2)*(abs(n)+l2p)*
     &    (abs(n+m+l4)))*(2d0*abs(n+m)+abs(m))
 
       enddo
       sum2(m) = ssum1
      enddo

      
C
      sum1 = dcmplx(0.0d0,0.0d0)
        
      do i = 1, nm     
      sum1= sum1+dcmplx(y(i),y(i+nm))*dcmplx(y(i),-y(i+nm))/i
      enddo
C
	  do i = 1,nm
	  V(i) = 0.5d0*dcmplx(y(i),-y(i+nm))/(i**2)
     &        -0.5d0*dcmplx(y(i),-y(i+nm))*sum1/i 
     &        + sum2(i) 

	  enddo  
	  
 

C      
Cc$omp parallel  
Cc$omp do
      do m=-nm,nm
      k=m
      if (k.eq.0) then
      k =0
      else
      k=1
      endif
      resd2=dcmplx(0d0,0d0)
      do n=1,nm
      resd2 = resd2+1d0*k*abs(n)*(-beta(n,m)*B(-n)-beta(-n,m)*B(n)+
     &  gamma(n,m)*Bd(-n)+gamma(-n,m)*Bd(n)+gammad(n,m)*B(-n)+
     & gammad(-n,m)*B(n))
      enddo
      res4(m) = resd2
      enddo 
Cc$omp end do
Cc$omp end parallel

C            
      do i =1,nm
      res2(i) = res4(i) + 0.5d0*V(i)
      res2(-nm+i-1) = dconjg(res2(i))
      enddo 
      
      do i =1,nm
      res2o(i) = res2(i)
      res2o(-nm+i-1) =dconjg(res2(i))
      enddo 
   
      do i=1,nm
      rn(i)=dreal(res2(i))
      rn(i+nm) =dimag(res2(i))
      enddo
        
      do m= 1,nm
      resd = dcmplx(0d0,0d0)
      do j =-nm,nm
      if  (j.eq.0) then 
      k1 = 0
      else
      k1 = 1
      end if
      resd = resd + k1*q(m,j)*dcmplx(yp(abs(j) ),
     & 	sign(1,j)*yp(abs(j) + nm)) 
     &  + k1*q(j,m)*dcmplx(yp(abs(j)),
     &  sign(1,j)*yp(abs(j) + nm))
      enddo   
      Q5(m) = resd
      enddo

  
      do i =1,nm
      do j=1,nm


      Am(i,j) =   dreal(q(i,j)+q(j,i)+q(i,-j)+q(-j,i))
      Am(i,nm+j) = -dimag(q(i,j)+q(j,i)-q(i,-j)-q(-j,i))
      Am(nm+i,j) = dimag(q(i,j)+q(j,i)+q(i,-j)+q(-j,i))
      Am(nm+i,nm+j) = dreal(q(i,j)+q(j,i)-q(i,-j)-q(-j,i))

      enddo 
      enddo 
C      
       Ao = Am
       Amo=Am
C       print*, Am
       call dgesv(2*nm, 1, Ao, 2*nm, ipiv2, rn, 2*nm, info)

       do i =1,2*nm
       RK1(i) = -rn(i)
       enddo 

C       print*, -RK1(1)/(c**2)
       
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C ----------------------------------------------------------------
C       Next we find K2

        do i =1, 2*nm
        y(i) = yo(i)+0.5d0*dt*ypo(i)
        y(i+2*nm) = ypo(i) + 0.5d0*dt*RK1(i)
        enddo
        yp = 0d0
        do i =1,2*nm
        yp(i) = y(i+2*nm)
        enddo 
            
            
C            print*, y
            
            
        CI=dcmplx(0.0d0,1.0d0)
      
C      Setting up the sum in B(n) 

      do n=-nm,nm
      so=dcmplx(0d0,0d0)
      do i =-nm,nm
C      do j=-nm,nm
      k1=n-i
      if (k1.eq.0) then
      k1 = 1
      else
      k1 = 0
      endif
      k2=i
      if (k2.eq.0) then
      k2 = 1
      else
      k2=0
      end if
            k3 =abs(n-i)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if

      so = so + 0.25d0*dcmplx(y(2*niv+abs(i)),sign(1,i)*
     & y(3*niv+abs(i)))*(1-k2)*(1-k1)
     & *k3*dcmplx(y(abs(n-i)),sign(1,n-i)*y(niv+abs(n-i)))/
     & ((abs(i)+k2)*(abs(n-i)+k1))
     & *(n-i)*((1-k2)*sign(1,i)-(1-k1)*sign(1,n-i))  
C	  enddo
	  enddo
	  sb(n)= so
	  enddo

C	  
C	  print*, sb-sb2
	  yd=dcmplx(0d0,0d0)
	  do i =1,nm
	  yd=yd-0.5d0*(dcmplx(y(2*niv+i),y(3*niv+i))*
     & dcmplx(y(i),-y(niv+i))
     & +dcmplx(y(2*niv+i),-y(3*niv+i))*dcmplx(y(i),y(niv+i)))/i
      enddo
      
      ydd=dcmplx(0d0,0d0)
      do i =1,nm
	  ydd=ydd+(dcmplx(y(2*niv+i),y(3*niv+i))*
     & dcmplx(y(2*niv+i),-y(3*niv+i)))/i
      enddo
      
C	  print*, yd
C	  print*, ydd
	  
C	  Next, we set up B(n)
	  do n=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else 
	  k1 =0
	  end if  
	  
	   B(n)=dcmplx(0.0d0,1.0d0)/(n+k1)*( sb(n)-
     & (1-k1)*dcmplx(y(2*nm+abs(n)),sign(1,n)* y(3*nm+abs(n)))/
     & (2d0*(abs(n)+k1))-0.5d0*yd*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*(1-k1)-
     &  yd*k1) 
	 
	  enddo

      
      do n=-nm,nm
      so=dcmplx(0d0,0d0)
      do i =-nm,nm
      k1=n-i
      if (k1.eq.0) then
      k1 = 1
      else
      k1 = 0
      endif
      k2=i
      if (k2.eq.0) then
      k2 = 1
      else
      k2=0
      end if
      k3 =abs(n-i)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if
      k4 = n
      if (k4.eq.0) then
      k4 =0 
      else 
      k4 =1
      end if
      so = so +0.25d0*((1-k1)*
     & dcmplx(y(2*niv+abs(i)),sign(1,i)*y(3*niv+abs(i)))*k3*(1-k2)
     & *dcmplx(y(abs(n-i)+2*niv),sign(1,n-i)*y(3*niv+abs(n-i))))/
     & ((abs(i)+k2)*(abs(n-i)+k1))
     & *(n-i)*(sign(1,i)-sign(1,n-i))*k4
	  enddo
	  sbd(n)=so
	  enddo
C	  print*, sbd
C	  print*, (yp(2*niv+1),yp(3*niv+1))
C	  Next, we set up Bd(n)

	  do n=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else 
	  k1 =0
	  end if
	  Bd(n) = dcmplx(0,1)/(n+k1)*(-0.5d0*yd*
     & (1-k1)*dcmplx(y(abs(n)+2*niv),sign(1,n)*y(3*niv+abs(n))) + 
     & ydd*dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*(1-k1)/2+ 
     & k1*ydd+  sbd(n))
	  enddo
C	  Bd(0)=dcmplx(0d0,0d0)
	  
C	  print*, Bd(1)
C	  -0.5d0*dcmplx(0,1)*(c**2)*dcmplx(y(1),y(niv+1))	  
C	  Then we set up beta(n,m)
	  
	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =abs(n-m)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if
      k4 =m
      if (k4.eq.0) then
      k4 =0
      else
      k4=1
      end if
	  beta(n,m) = dcmplx(0,1d0)/(n+k1)*(0.25d0*k3*(1-k1)*(1-k2)*
     & dcmplx(y(2*niv+abs(n-m)), sign(1,n-m)*y(3*niv+abs(n-m)))
     & /((abs(m)+1-k4)*(abs(n-m)+k2))*m*((1-k2)*sign(1,n-m)-
     & (k4)* sign(1,m))-0.5d0*yd*k2*k4
     & +0.25d0*dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*
     & dcmplx(y(2*niv+abs(m)),-sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k4))*k4+
     & 0.5d0*k1*dcmplx(y(2*niv+abs(m)),sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k4)))
	  enddo
	  enddo 
	  
	  	  
	  
C	  print*, beta(15,-32:32)
C	  do i =1,nm
C      diff(i)  = beta(-i,i)
C      enddo
C	  print*, beta(0,1)
C	  Then set up gamma(n,m)

	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =n-m
	  if (k3.eq.0) then
	  k3=1
	  else
	  k3=0
	  endif
	  	  k4 =abs(n-m)
      if (k4.gt.nm) then
      k4 =0
      else
      k4=1
      end if
      k5 =m
      if (k5.eq.0) then
      k5 =0
      else
      k5=1
      end if
	  gamma(n,m) = (1-k1)*dcmplx(0,1d0)/(n+k1)*(-0.5d0
     & *dcmplx(k3*1d0,0d0)/(abs(m)+ 1-k5) + k5*(1-k2)*
     &  0.25d0*k4*dcmplx(y(abs(n-m)), sign(1,n-m)*y(niv+abs(n-m)))
     & /((abs(m)+1-k5)
     & *(abs(n-m)+k2))*(n-m)*((k5)*sign(1,m)-(1-k2)*sign(1,n-m))+
     & 0.25d0*dcmplx(y(abs(m)),-sign(1,m)*y(niv+abs(m)))*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))/abs(m+1-k5)
     &*k5*(1-k1)+0.5d0*k1*dcmplx(y(abs(m)),sign(1,m)*y(niv+abs(m)))/
     & abs(m+(1-k5)))
	  enddo
	  enddo 	  

C      print*, gamma
C      print*, dcmplx(0d0,1d0)/1d0*dcmplx(y(1),y(niv+1))
C	   next we compute gammad

	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =n-m
	  if (k3.eq.0) then
	  k3=1
	  else
	  k3=0
	  endif
	  	  k4 =abs(n-m)
      if (k4.gt.nm) then
      k4 =0
      else
      k4=1
      end if
      k5 =m
      if (k5.eq.0) then
      k5 =0
      else
      k5=1
      end if
	  gammad(n,m) = (1-k1)*dcmplx(0,1d0)/(n+k1)*(0.25d0*k5*
     & k4*(1-k3)*dcmplx(y(2*niv+abs(n-m)), 
     & sign(1,n-m)*y(3*niv+abs(n-m)))
     & /((abs(m)+1-k5)
     & *(abs(n-m)+k2))*(n-m)*((k5)*sign(1,m)-(1-k2)*sign(1,n-m))+
     & 0.25d0*dcmplx(y(abs(m)+2*niv),-sign(1,m)*y(3*niv+abs(m)))*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))
     & /abs(m+1-k5)*k5+
     & 0.25d0*dcmplx(y(abs(m)),-sign(1,m)*y(niv+abs(m)))*
     & dcmplx(y(2*niv+abs(n)),sign(1,n)*y(3*niv+abs(n)))
     & /abs(m+1-k5)*k5+
     & 0.5d0*k1*dcmplx(y(2*niv+abs(m)),sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k5)))
	  enddo
	  enddo
C	  print*, gammad   
C	  print*, gammad(1,2) + 0.25d0*dcmplx(y(2*niv+1),-y(3*niv+1))*
C     &dcmplx(0d0,1d0) 


	  do  m = -nm, nm
          do j = -nm, nm
	      kk = m-j
	      l = m-j
	      k1 = m
	      k2 = j
	      k3a = m
	      k3b = j
	      k4 = m-j
	      k5 = j
	      k6 = m-j
	                         
	      if (abs(kk) .gt. nm) then
	      kk = 0
	      else 
	      kk = 1
	      end if
	      if (k1 .eq. 0) then
	      k1 = 0
	      else
	      k1=1
	      end if
	      if (k2.eq.0) then
	      k2 = 0
	      else
	      k2 = 1
	      end if
	      if (k3a .eq. 0) then
	      k3a=1
	      else 
	      k3a=0
	      end if
	      if (k3b .eq. 0) then
	      k3b=1
	      else 
	      k3b=0
	      end if
	      if (k4 .eq. 0) then
	      k4 = 0
	      else
	      k4 = 1
	      end if
	      if (k5 .eq. 0) then
	      k5=2
	      else
	      k5=1
	      end if
	      if (k6 .eq. 0) then
	      k6=2
	      else
	      k6=1
	      end if
	      
	  p(m,j)=(0.5d0*k5*k6*kk*abs(k4*sign(1,m-j)-k2*sign(1,j))*dcmplx
     &(k4*y(abs(m-j))+1-k4,k4*y(nm+abs(m-j))*sign(1,m-j))-
     & 0.5d0*dcmplx(k1
     &*y(abs(m))+(1-k1), k1*y(abs(m)+nm)*sign(1, m))*dcmplx(k2*y(abs(j
     &))+(1-k2),-k2*y(abs(j)+nm)*sign(1,j)))/(DSQRT(1.0d0*(abs(m)+k3a
     & ))*
     & (k3b+abs(j)))
 	 
      enddo  
      enddo	
      
      

CC      Next we set up the leading term coefficient Q

      Mt= 2*nm+1
      alpha = dcmplx(0.25d0,0.0d0)
      betao = dcmplx(0.0d0,0.0d0)
      


C       seconds = omp_get_wtime ( )  

      call ZGEMM('T','N',Mt,Mt,nm,alpha, p(1:nm,-nm:nm),nm,
     &p(-1:-nm:-1,-nm:nm),nm,betao, q ,Mt)
      

C       seconds = omp_get_wtime ( ) - seconds;
C       print*, seconds 
C      seconds = omp_get_wtime ( )

C       open(unit=10,file='qr.txt',
C     &  status = 'unknown')
C       do i=-nm,nm
C       do j=-nm,nm
C       write(10,*)  dreal(q(i,j))
C       enddo
C       enddo 
C       close(unit =10)
C  111 format (D20.10)  
C  
C       open(unit=10,file='qi.txt',
C     &  status = 'unknown')
C       do i=-nm,nm
C       do j=-nm,nm
C       write(10,*)  dimag(q(i,j))
C       enddo
C       enddo 
C       close(unit =10)
C  111 format (D20.10)    

Cc$omp parallel  
Cc$omp do

       do m = 1, nm
       ssum1 = dcmplx(0.0d0,0.0d0)
       do n = -nm, nm
        l2=m
        l2p=n
        l2pp=m*n
        l3=(n+m)*n
        l4=(n+m)
        l5=abs(n+m)

         if (l2 .eq. 0) then
          l2 = 1
         else 
          l2 = 0     
         end if 
         if (l2p .eq. 0) then
          l2p = 1
         else 
          l2p = 0     
         end if  
         if (l2pp .eq. 0) then
          l2pp = 1
         else 
          l2pp = 0     
         end if 
         if (l3 .eq. 0) then
          l3 = 0
         else 
          l3 = 1     
         end if 
         if (l4 .eq. 0) then
          l4 = 1
         else 
          l4 = 0     
         end if 
        if (l5 .gt. nm) then
          l5 = 0
         else 
          l5 = 1     
         end if 
    
              
         ssum1 = ssum1 + l5*l3*dcmplx(y(abs(n)),
     &    sign(1,n)*y(abs(n)+nm))*dcmplx(y(abs(n+m)),
     &    -sign(1,n+m)*y(abs(n+m)+nm))/(8d0*(abs(m)+l2)*(abs(n)+l2p)*
     &    (abs(n+m+l4)))*(2d0*abs(n+m)+abs(m))
 
       enddo
       sum2(m) = ssum1
      enddo
      
Cc$omp end do
Cc$omp end parallel




      
C
      sum1 = dcmplx(0.0d0,0.0d0)
        
      do i = 1, nm     
      sum1= sum1+dcmplx(y(i),y(i+nm))*dcmplx(y(i),-y(i+nm))/i
      enddo
C
	  do i = 1,nm
	  V(i) = 0.5d0*dcmplx(y(i),-y(i+nm))/(i**2)
     &        -0.5d0*dcmplx(y(i),-y(i+nm))*sum1/i 
     &        + sum2(i) 

	  enddo  
	  
 

C      
Cc$omp parallel  
Cc$omp do
      do m=-nm,nm
      k=m
      if (k.eq.0) then
      k =0
      else
      k=1
      endif
      resd2=dcmplx(0d0,0d0)
      do n=1,nm
      resd2 = resd2+1d0*k*abs(n)*(-beta(n,m)*B(-n)-beta(-n,m)*B(n)+
     &  gamma(n,m)*Bd(-n)+gamma(-n,m)*Bd(n)+gammad(n,m)*B(-n)+
     & gammad(-n,m)*B(n))
      enddo
C      res1(m) = 0.5d0*dreal(resd2)
C      res3(m) = 0.5d0*dimag(resd2)
      res4(m) = resd2
      enddo 
Cc$omp end do
Cc$omp end parallel


C      enddo 
      
      do i =1,nm
      res2(i) = res4(i) + 0.5d0*V(i)
      res2(-nm+i-1) = dconjg(res2(i))
      enddo 
   
      do i=1,nm
      rn(i)=dreal(res2(i))
      rn(i+nm) =dimag(res2(i))
      enddo
      
C     



      do i =1,nm
      do j=1,nm


      Am(i,j) =   dreal(q(i,j)+q(j,i)+q(i,-j)+q(-j,i))
      Am(i,nm+j) = - dimag(q(i,j)+q(j,i)-q(i,-j)-q(-j,i))
      Am(nm+i,j) = dimag(q(i,j)+q(j,i)+q(i,-j)+q(-j,i))
      Am(nm+i,nm+j) = dreal(q(i,j)+q(j,i)-q(i,-j)-q(-j,i))

      enddo 
      enddo 


    

C      
C      do i=1,2*nm
C      do j=1,2*nm
C      Ao(i,j) =Amat(2*nm+i,2*nm+j)
C      enddo
C      enddo 
C      
       Ao = Am
      
       call dgesv(2*nm, 1, Ao, 2*nm, ipiv2, rn, 2*nm, info )
       rn =-rn
       
C       print*, rn
       do i =1,2*nm
       RK2(i) = rn(i)
       enddo
       
C       print*, RK2
       
C        This is K2
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------
C        Next we find K3
C        y=0d0
        do i =1, 2*nm
        y(i) = yo(i)+dt*0.5d0*(ypo(i)+0.5d0*dt*RK1(i))
        y(i+2*nm) = ypo(i) + 0.5d0*dt*RK2(i)
        enddo
        yp = 0d0
        do i =1,2*nm
        yp(i) = y(i+2*nm)
        enddo 
        
C        print*, y
            
        CI=dcmplx(0.0d0,1.0d0)
      
C      Setting up the sum in B(n) 

      do n=-nm,nm
      so=dcmplx(0d0,0d0)
      do i =-nm,nm
C      do j=-nm,nm
      k1=n-i
      if (k1.eq.0) then
      k1 = 1
      else
      k1 = 0
      endif
      k2=i
      if (k2.eq.0) then
      k2 = 1
      else
      k2=0
      end if
            k3 =abs(n-i)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if

      so = so + 0.25d0*dcmplx(y(2*niv+abs(i)),sign(1,i)*
     & y(3*niv+abs(i)))*(1-k2)*(1-k1)
     & *k3*dcmplx(y(abs(n-i)),sign(1,n-i)*y(niv+abs(n-i)))/
     & ((abs(i)+k2)*(abs(n-i)+k1))
     & *(n-i)*((1-k2)*sign(1,i)-(1-k1)*sign(1,n-i))  
C	  enddo
	  enddo
	  sb(n)= so
	  enddo

C	  
C	  print*, sb-sb2
	  yd=dcmplx(0d0,0d0)
	  do i =1,nm
	  yd=yd-0.5d0*(dcmplx(y(2*niv+i),y(3*niv+i))*
     & dcmplx(y(i),-y(niv+i))
     & +dcmplx(y(2*niv+i),-y(3*niv+i))*dcmplx(y(i),y(niv+i)))/i
      enddo
      
      ydd=dcmplx(0d0,0d0)
      do i =1,nm
	  ydd=ydd+(dcmplx(y(2*niv+i),y(3*niv+i))*
     & dcmplx(y(2*niv+i),-y(3*niv+i)))/i
      enddo
      
C	  print*, yd
C	  print*, ydd
	  
C	  Next, we set up B(n)
	  do n=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else 
	  k1 =0
	  end if  
	  
	   B(n)=dcmplx(0.0d0,1.0d0)/(n+k1)*( sb(n)-
     & (1-k1)*dcmplx(y(2*nm+abs(n)),sign(1,n)* y(3*nm+abs(n)))/
     & (2d0*(abs(n)+k1))-0.5d0*yd*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*(1-k1)-
     &  yd*k1) 
	 
	  enddo

C      print*, B
      
      do n=-nm,nm
      so=dcmplx(0d0,0d0)
      do i =-nm,nm
      k1=n-i
      if (k1.eq.0) then
      k1 = 1
      else
      k1 = 0
      endif
      k2=i
      if (k2.eq.0) then
      k2 = 1
      else
      k2=0
      end if
      k3 =abs(n-i)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if
      k4 = n
      if (k4.eq.0) then
      k4 =0 
      else 
      k4 =1
      end if
      so = so +0.25d0*((1-k1)*
     & dcmplx(y(2*niv+abs(i)),sign(1,i)*y(3*niv+abs(i)))*k3*(1-k2)
     & *dcmplx(y(abs(n-i)+2*niv),sign(1,n-i)*y(3*niv+abs(n-i))))/
     & ((abs(i)+k2)*(abs(n-i)+k1))
     & *(n-i)*(sign(1,i)-sign(1,n-i))*k4
	  enddo
	  sbd(n)=so
	  enddo
C	  print*, sbd
C	  print*, (yp(2*niv+1),yp(3*niv+1))
C	  Next, we set up Bd(n)

	  do n=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else 
	  k1 =0
	  end if
	  Bd(n) = dcmplx(0,1)/(n+k1)*(-0.5d0*yd*
     & (1-k1)*dcmplx(y(abs(n)+2*niv),sign(1,n)*y(3*niv+abs(n))) + 
     & ydd*dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*(1-k1)/2+ 
     & k1*ydd+  sbd(n))
	  enddo
C	  Bd(0)=dcmplx(0d0,0d0)
	  
C	  print*, Bd
C	  -0.5d0*dcmplx(0,1)*(c**2)*dcmplx(y(1),y(niv+1))	  
C	  Then we set up beta(n,m)
	  beta=0d0
	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =abs(n-m)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if
      k4 =m
      if (k4.eq.0) then
      k4 =0
      else
      k4=1
      end if
	  beta(n,m) = dcmplx(0,1d0)/(n+k1)*(0.25d0*k3*(1-k1)*(1-k2)*
     & dcmplx(y(2*niv+abs(n-m)), sign(1,n-m)*y(3*niv+abs(n-m)))
     & /((abs(m)+1-k4)*(abs(n-m)+k2))*m*((1-k2)*sign(1,n-m)-
     & (k4)* sign(1,m))-0.5d0*yd*k2*k4
     & +0.25d0*dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*
     & dcmplx(y(2*niv+abs(m)),-sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k4))*k4+
     & 0.5d0*k1*dcmplx(y(2*niv+abs(m)),sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k4)))
	  enddo
	  enddo 
	  
	  	  
	  
C	  print*, beta
C	  do i =1,nm
C      diff(i)  = beta(-i,i)
C      enddo
C	  print*, beta(0,1)
C	  Then set up gamma(n,m)

	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =n-m
	  if (k3.eq.0) then
	  k3=1
	  else
	  k3=0
	  endif
	  	  k4 =abs(n-m)
      if (k4.gt.nm) then
      k4 =0
      else
      k4=1
      end if
      k5 =m
      if (k5.eq.0) then
      k5 =0
      else
      k5=1
      end if
	  gamma(n,m) = (1-k1)*dcmplx(0,1d0)/(n+k1)*(-0.5d0
     & *dcmplx(k3*1d0,0d0)/(abs(m)+ 1-k5) + k5*(1-k2)*
     &  0.25d0*k4*dcmplx(y(abs(n-m)), sign(1,n-m)*y(niv+abs(n-m)))
     & /((abs(m)+1-k5)
     & *(abs(n-m)+k2))*(n-m)*((k5)*sign(1,m)-(1-k2)*sign(1,n-m))+
     & 0.25d0*dcmplx(y(abs(m)),-sign(1,m)*y(niv+abs(m)))*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))/abs(m+1-k5)
     &*k5*(1-k1)+0.5d0*k1*dcmplx(y(abs(m)),sign(1,m)*y(niv+abs(m)))/
     & abs(m+(1-k5)))
	  enddo
	  enddo 	  

C      print*, gamma
C      print*, dcmplx(0d0,1d0)/1d0*dcmplx(y(1),y(niv+1))
C	   next we compute gammad

	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =n-m
	  if (k3.eq.0) then
	  k3=1
	  else
	  k3=0
	  endif
	  	  k4 =abs(n-m)
      if (k4.gt.nm) then
      k4 =0
      else
      k4=1
      end if
      k5 =m
      if (k5.eq.0) then
      k5 =0
      else
      k5=1
      end if
	  gammad(n,m) = (1-k1)*dcmplx(0,1d0)/(n+k1)*(0.25d0*k5*
     & k4*(1-k3)*dcmplx(y(2*niv+abs(n-m)), 
     & sign(1,n-m)*y(3*niv+abs(n-m)))
     & /((abs(m)+1-k5)
     & *(abs(n-m)+k2))*(n-m)*((k5)*sign(1,m)-(1-k2)*sign(1,n-m))+
     & 0.25d0*dcmplx(y(abs(m)+2*niv),-sign(1,m)*y(3*niv+abs(m)))*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))
     & /abs(m+1-k5)*k5+
     & 0.25d0*dcmplx(y(abs(m)),-sign(1,m)*y(niv+abs(m)))*
     & dcmplx(y(2*niv+abs(n)),sign(1,n)*y(3*niv+abs(n)))
     & /abs(m+1-k5)*k5+
     & 0.5d0*k1*dcmplx(y(2*niv+abs(m)),sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k5)))
	  enddo
	  enddo
C	  print*, gammad   
C	  print*, gammad(1,2) + 0.25d0*dcmplx(y(2*niv+1),-y(3*niv+1))*
C     &dcmplx(0d0,1d0) 


	  do  m = -nm, nm
          do j = -nm, nm
	      kk = m-j
	      l = m-j
	      k1 = m
	      k2 = j
	      k3a = m
	      k3b = j
	      k4 = m-j
	      k5 = j
	      k6 = m-j
	                         
	      if (abs(kk) .gt. nm) then
	      kk = 0
	      else 
	      kk = 1
	      end if
	      if (k1 .eq. 0) then
	      k1 = 0
	      else
	      k1=1
	      end if
	      if (k2.eq.0) then
	      k2 = 0
	      else
	      k2 = 1
	      end if
	      if (k3a .eq. 0) then
	      k3a=1
	      else 
	      k3a=0
	      end if
	      if (k3b .eq. 0) then
	      k3b=1
	      else 
	      k3b=0
	      end if
	      if (k4 .eq. 0) then
	      k4 = 0
	      else
	      k4 = 1
	      end if
	      if (k5 .eq. 0) then
	      k5=2
	      else
	      k5=1
	      end if
	      if (k6 .eq. 0) then
	      k6=2
	      else
	      k6=1
	      end if
	      
	  p(m,j)=(0.5d0*k5*k6*kk*abs(k4*sign(1,m-j)-k2*sign(1,j))*dcmplx
     &(k4*y(abs(m-j))+1-k4,k4*y(nm+abs(m-j))*sign(1,m-j))-
     & 0.5d0*dcmplx(k1
     &*y(abs(m))+(1-k1), k1*y(abs(m)+nm)*sign(1, m))*dcmplx(k2*y(abs(j
     &))+(1-k2),-k2*y(abs(j)+nm)*sign(1,j)))/(DSQRT(1.0d0*(abs(m)+k3a
     & ))*
     & (k3b+abs(j)))
 	 
      enddo  
      enddo	
      
      

CC      Next we set up the leading term coefficient Q

      Mt= 2*nm+1
      alpha = dcmplx(0.25d0,0.0d0)
      betao = dcmplx(0.0d0,0.0d0)
      


C       seconds = omp_get_wtime ( )  

      call ZGEMM('T','N',Mt,Mt,nm,alpha, p(1:nm,-nm:nm),nm,
     &p(-1:-nm:-1,-nm:nm),nm,betao, q ,Mt)


C       seconds = omp_get_wtime ( ) - seconds;
C       print*, seconds 
C      seconds = omp_get_wtime ( )

C       open(unit=10,file='qr.txt',
C     &  status = 'unknown')
C       do i=-nm,nm
C       do j=-nm,nm
C       write(10,*)  dreal(q(i,j))
C       enddo
C       enddo 
C       close(unit =10)
C  111 format (D20.10)  
C  
C       open(unit=10,file='qi.txt',
C     &  status = 'unknown')
C       do i=-nm,nm
C       do j=-nm,nm
C       write(10,*)  dimag(q(i,j))
C       enddo
C       enddo 
C       close(unit =10)
C  111 format (D20.10)    

Cc$omp parallel  
Cc$omp do

       do m = 1, nm
       ssum1 = dcmplx(0.0d0,0.0d0)
       do n = -nm, nm
        l2=m
        l2p=n
        l2pp=m*n
        l3=(n+m)*n
        l4=(n+m)
        l5=abs(n+m)

         if (l2 .eq. 0) then
          l2 = 1
         else 
          l2 = 0     
         end if 
         if (l2p .eq. 0) then
          l2p = 1
         else 
          l2p = 0     
         end if  
         if (l2pp .eq. 0) then
          l2pp = 1
         else 
          l2pp = 0     
         end if 
         if (l3 .eq. 0) then
          l3 = 0
         else 
          l3 = 1     
         end if 
         if (l4 .eq. 0) then
          l4 = 1
         else 
          l4 = 0     
         end if 
        if (l5 .gt. nm) then
          l5 = 0
         else 
          l5 = 1     
         end if 
    
              
         ssum1 = ssum1 + l5*l3*dcmplx(y(abs(n)),
     &    sign(1,n)*y(abs(n)+nm))*dcmplx(y(abs(n+m)),
     &    -sign(1,n+m)*y(abs(n+m)+nm))/(8d0*(abs(m)+l2)*(abs(n)+l2p)*
     &    (abs(n+m+l4)))*(2d0*abs(n+m)+abs(m))
 
       enddo
       sum2(m) = ssum1
      enddo
      
Cc$omp end do
Cc$omp end parallel




      
C
      sum1 = dcmplx(0.0d0,0.0d0)
        
      do i = 1, nm     
      sum1= sum1+dcmplx(y(i),y(i+nm))*dcmplx(y(i),-y(i+nm))/i
      enddo
C
	  do i = 1,nm
	  V(i) = 0.5d0*dcmplx(y(i),-y(i+nm))/(i**2)
     &        -0.5d0*dcmplx(y(i),-y(i+nm))*sum1/i 
     &        + sum2(i) 

	  enddo  
	  
 

C      
Cc$omp parallel  
Cc$omp do
      do m=-nm,nm
      k=m
      if (k.eq.0) then
      k =0
      else
      k=1
      endif
      resd2=dcmplx(0d0,0d0)
      do n=1,nm
      resd2 = resd2+1d0*k*abs(n)*(-beta(n,m)*B(-n)-beta(-n,m)*B(n)+
     &  gamma(n,m)*Bd(-n)+gamma(-n,m)*Bd(n)+gammad(n,m)*B(-n)+
     & gammad(-n,m)*B(n))
      enddo
C      res1(m) = 0.5d0*dreal(resd2)
C      res3(m) = 0.5d0*dimag(resd2)
      res4(m) = resd2
      enddo 
Cc$omp end do
Cc$omp end parallel
C      print*, gammad

      do i =1,nm
      res2(i) = res4(i) + 0.5d0*V(i)
      res2(-nm+i-1) = dconjg(res2(i))
      enddo 
   
      do i=1,nm
      rn(i)=dreal(res2(i))
      rn(i+nm) =dimag(res2(i))
      enddo
      
C     
   
      do i =1,nm
      do j=1,nm


      Am(i,j) =   dreal(q(i,j)+q(j,i)+q(i,-j)+q(-j,i))
      Am(i,nm+j) = - dimag(q(i,j)+q(j,i)-q(i,-j)-q(-j,i))
      Am(nm+i,j) = dimag(q(i,j)+q(j,i)+q(i,-j)+q(-j,i))
      Am(nm+i,nm+j) = dreal(q(i,j)+q(j,i)-q(i,-j)-q(-j,i))

      enddo 
      enddo 

       Ao = Am
      
       call dgesv(2*nm, 1, Ao, 2*nm, ipiv2, rn, 2*nm, info )
       rn =-rn
C       print*, rn
       

       do i =1,2*nm
       RK3(i) = rn(i)
       enddo
        
C        print*, RK3
C------------------------------------------------------------------
C------------------------------------------------------------------
C------------------------------------------------------------------
C        Next we find K4
C        do i =1, 2*nm
C        y(i) = yo(i)+dt*0.5d0*(ypo(i)+0.5d0*dt*RK1(i))
C        y(i+2*nm) = ypo(i) + 0.5d0*dt*RK2(i)
C        enddo
        
        do i =1, 2*nm
        y(i) = yo(i)+dt*(ypo(i)+0.5d0*dt*RK2(i))
        y(i+2*nm) = ypo(i) + dt*RK3(i)
        enddo
        yp = 0d0
        do i =1,2*nm
        yp(i) = y(i+2*nm)
        enddo 
        
       
        CI=dcmplx(0.0d0,1.0d0)
      
C      Setting up the sum in B(n) 

      do n=-nm,nm
      so=dcmplx(0d0,0d0)
      do i =-nm,nm
C      do j=-nm,nm
      k1=n-i
      if (k1.eq.0) then
      k1 = 1
      else
      k1 = 0
      endif
      k2=i
      if (k2.eq.0) then
      k2 = 1
      else
      k2=0
      end if
            k3 =abs(n-i)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if

      so = so + 0.25d0*dcmplx(y(2*niv+abs(i)),sign(1,i)*
     & y(3*niv+abs(i)))*(1-k2)*(1-k1)
     & *k3*dcmplx(y(abs(n-i)),sign(1,n-i)*y(niv+abs(n-i)))/
     & ((abs(i)+k2)*(abs(n-i)+k1))
     & *(n-i)*((1-k2)*sign(1,i)-(1-k1)*sign(1,n-i))  
C	  enddo
	  enddo
	  sb(n)= so
	  enddo

C	  
C	  print*, sb-sb2
	  yd=dcmplx(0d0,0d0)
	  do i =1,nm
	  yd=yd-0.5d0*(dcmplx(y(2*niv+i),y(3*niv+i))*
     & dcmplx(y(i),-y(niv+i))
     & +dcmplx(y(2*niv+i),-y(3*niv+i))*dcmplx(y(i),y(niv+i)))/i
      enddo
      
      ydd=dcmplx(0d0,0d0)
      do i =1,nm
	  ydd=ydd+(dcmplx(y(2*niv+i),y(3*niv+i))*
     & dcmplx(y(2*niv+i),-y(3*niv+i)))/i
      enddo
      
C	  print*, yd
C	  print*, ydd
	  
C	  Next, we set up B(n)
	  do n=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else 
	  k1 =0
	  end if  
	  
	   B(n)=dcmplx(0.0d0,1.0d0)/(n+k1)*( sb(n)-
     & (1-k1)*dcmplx(y(2*nm+abs(n)),sign(1,n)* y(3*nm+abs(n)))/
     & (2d0*(abs(n)+k1))-0.5d0*yd*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*(1-k1)-
     &  yd*k1) 
	 
	  enddo

C      print*, B
      
      do n=-nm,nm
      so=dcmplx(0d0,0d0)
      do i =-nm,nm
      k1=n-i
      if (k1.eq.0) then
      k1 = 1
      else
      k1 = 0
      endif
      k2=i
      if (k2.eq.0) then
      k2 = 1
      else
      k2=0
      end if
      k3 =abs(n-i)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if
      k4 = n
      if (k4.eq.0) then
      k4 =0 
      else 
      k4 =1
      end if
      so = so +0.25d0*((1-k1)*
     & dcmplx(y(2*niv+abs(i)),sign(1,i)*y(3*niv+abs(i)))*k3*(1-k2)
     & *dcmplx(y(abs(n-i)+2*niv),sign(1,n-i)*y(3*niv+abs(n-i))))/
     & ((abs(i)+k2)*(abs(n-i)+k1))
     & *(n-i)*(sign(1,i)-sign(1,n-i))*k4
	  enddo
	  sbd(n)=so
	  enddo
C	  print*, sbd
C	  print*, (yp(2*niv+1),yp(3*niv+1))
C	  Next, we set up Bd(n)

	  do n=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else 
	  k1 =0
	  end if
	  Bd(n) = dcmplx(0,1)/(n+k1)*(-0.5d0*yd*
     & (1-k1)*dcmplx(y(abs(n)+2*niv),sign(1,n)*y(3*niv+abs(n))) + 
     & ydd*dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*(1-k1)/2+ 
     & k1*ydd+  sbd(n))
	  enddo
C	  Bd(0)=dcmplx(0d0,0d0)
	  
C	  print*, Bd
C	  -0.5d0*dcmplx(0,1)*(c**2)*dcmplx(y(1),y(niv+1))	  
C	  Then we set up beta(n,m)
	  
	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =abs(n-m)
      if (k3.gt.nm) then
      k3 =0
      else
      k3=1
      end if
      k4 =m
      if (k4.eq.0) then
      k4 =0
      else
      k4=1
      end if
	  beta(n,m) = dcmplx(0,1d0)/(n+k1)*(0.25d0*k3*(1-k1)*(1-k2)*
     & dcmplx(y(2*niv+abs(n-m)), sign(1,n-m)*y(3*niv+abs(n-m)))
     & /((abs(m)+1-k4)*(abs(n-m)+k2))*m*((1-k2)*sign(1,n-m)-
     & (k4)* sign(1,m))-0.5d0*yd*k2*k4
     & +0.25d0*dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))*
     & dcmplx(y(2*niv+abs(m)),-sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k4))*k4+
     & 0.5d0*k1*dcmplx(y(2*niv+abs(m)),sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k4)))
	  enddo
	  enddo 
	  
	  	  
	  
C	  print*, beta
C	  do i =1,nm
C      diff(i)  = beta(-i,i)
C      enddo
C	  print*, beta(0,1)
C	  Then set up gamma(n,m)

	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =n-m
	  if (k3.eq.0) then
	  k3=1
	  else
	  k3=0
	  endif
	  	  k4 =abs(n-m)
      if (k4.gt.nm) then
      k4 =0
      else
      k4=1
      end if
      k5 =m
      if (k5.eq.0) then
      k5 =0
      else
      k5=1
      end if
	  gamma(n,m) = (1-k1)*dcmplx(0,1d0)/(n+k1)*(-0.5d0
     & *dcmplx(k3*1d0,0d0)/(abs(m)+ 1-k5) + k5*(1-k2)*
     &  0.25d0*k4*dcmplx(y(abs(n-m)), sign(1,n-m)*y(niv+abs(n-m)))
     & /((abs(m)+1-k5)
     & *(abs(n-m)+k2))*(n-m)*((k5)*sign(1,m)-(1-k2)*sign(1,n-m))+
     & 0.25d0*dcmplx(y(abs(m)),-sign(1,m)*y(niv+abs(m)))*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))/abs(m+1-k5)
     &*k5*(1-k1)+0.5d0*k1*dcmplx(y(abs(m)),sign(1,m)*y(niv+abs(m)))/
     & abs(m+(1-k5)))
	  enddo
	  enddo 	  

C      print*, gamma
C      print*, dcmplx(0d0,1d0)/1d0*dcmplx(y(1),y(niv+1))
C	   next we compute gammad

	  do n=-nm,nm
	  do m=-nm,nm
	  k1= n
	  if (k1.eq.0) then
	  k1 =1
	  else
	  k1 =0 
	  endif
	  k2 =n-m
	  if (k2.eq.0) then
	  k2 =1
	  else
	  k2=0
	  endif
	  k3 =n-m
	  if (k3.eq.0) then
	  k3=1
	  else
	  k3=0
	  endif
	  	  k4 =abs(n-m)
      if (k4.gt.nm) then
      k4 =0
      else
      k4=1
      end if
      k5 =m
      if (k5.eq.0) then
      k5 =0
      else
      k5=1
      end if
	  gammad(n,m) = (1-k1)*dcmplx(0,1d0)/(n+k1)*(0.25d0*k5*
     & k4*(1-k3)*dcmplx(y(2*niv+abs(n-m)), 
     & sign(1,n-m)*y(3*niv+abs(n-m)))
     & /((abs(m)+1-k5)
     & *(abs(n-m)+k2))*(n-m)*((k5)*sign(1,m)-(1-k2)*sign(1,n-m))+
     & 0.25d0*dcmplx(y(abs(m)+2*niv),-sign(1,m)*y(3*niv+abs(m)))*
     & dcmplx(y(abs(n)),sign(1,n)*y(niv+abs(n)))
     & /abs(m+1-k5)*k5+
     & 0.25d0*dcmplx(y(abs(m)),-sign(1,m)*y(niv+abs(m)))*
     & dcmplx(y(2*niv+abs(n)),sign(1,n)*y(3*niv+abs(n)))
     & /abs(m+1-k5)*k5+
     & 0.5d0*k1*dcmplx(y(2*niv+abs(m)),sign(1,m)*y(3*niv+abs(m)))/
     & abs(m+(1-k5)))
	  enddo
	  enddo
C	  print*, gammad   
C	  print*, gammad(1,2) + 0.25d0*dcmplx(y(2*niv+1),-y(3*niv+1))*
C     &dcmplx(0d0,1d0) 


	  do  m = -nm, nm
          do j = -nm, nm
	      kk = m-j
	      l = m-j
	      k1 = m
	      k2 = j
	      k3a = m
	      k3b = j
	      k4 = m-j
	      k5 = j
	      k6 = m-j
	                         
	      if (abs(kk) .gt. nm) then
	      kk = 0
	      else 
	      kk = 1
	      end if
	      if (k1 .eq. 0) then
	      k1 = 0
	      else
	      k1=1
	      end if
	      if (k2.eq.0) then
	      k2 = 0
	      else
	      k2 = 1
	      end if
	      if (k3a .eq. 0) then
	      k3a=1
	      else 
	      k3a=0
	      end if
	      if (k3b .eq. 0) then
	      k3b=1
	      else 
	      k3b=0
	      end if
	      if (k4 .eq. 0) then
	      k4 = 0
	      else
	      k4 = 1
	      end if
	      if (k5 .eq. 0) then
	      k5=2
	      else
	      k5=1
	      end if
	      if (k6 .eq. 0) then
	      k6=2
	      else
	      k6=1
	      end if
	      
	  p(m,j)=(0.5d0*k5*k6*kk*abs(k4*sign(1,m-j)-k2*sign(1,j))*dcmplx
     &(k4*y(abs(m-j))+1-k4,k4*y(nm+abs(m-j))*sign(1,m-j))-
     & 0.5d0*dcmplx(k1
     &*y(abs(m))+(1-k1), k1*y(abs(m)+nm)*sign(1, m))*dcmplx(k2*y(abs(j
     &))+(1-k2),-k2*y(abs(j)+nm)*sign(1,j)))/(DSQRT(1.0d0*(abs(m)+k3a
     & ))*
     & (k3b+abs(j)))
 	 
      enddo  
      enddo	
      
      

CC      Next we set up the leading term coefficient Q

      Mt= 2*nm+1
      alpha = dcmplx(0.25d0,0.0d0)
      betao = dcmplx(0.0d0,0.0d0)
      


C       seconds = omp_get_wtime ( )  

      call ZGEMM('T','N',Mt,Mt,nm,alpha, p(1:nm,-nm:nm),nm,
     &p(-1:-nm:-1,-nm:nm),nm,betao, q ,Mt)


C       seconds = omp_get_wtime ( ) - seconds;
C       print*, seconds 
C      seconds = omp_get_wtime ( )

C       open(unit=10,file='qr.txt',
C     &  status = 'unknown')
C       do i=-nm,nm
C       do j=-nm,nm
C       write(10,*)  dreal(q(i,j))
C       enddo
C       enddo 
C       close(unit =10)
C  111 format (D20.10)  
C  
C       open(unit=10,file='qi.txt',
C     &  status = 'unknown')
C       do i=-nm,nm
C       do j=-nm,nm
C       write(10,*)  dimag(q(i,j))
C       enddo
C       enddo 
C       close(unit =10)
C  111 format (D20.10)    

Cc$omp parallel  
Cc$omp do

       do m = 1, nm
       ssum1 = dcmplx(0.0d0,0.0d0)
       do n = -nm, nm
        l2=m
        l2p=n
        l2pp=m*n
        l3=(n+m)*n
        l4=(n+m)
        l5=abs(n+m)

         if (l2 .eq. 0) then
          l2 = 1
         else 
          l2 = 0     
         end if 
         if (l2p .eq. 0) then
          l2p = 1
         else 
          l2p = 0     
         end if  
         if (l2pp .eq. 0) then
          l2pp = 1
         else 
          l2pp = 0     
         end if 
         if (l3 .eq. 0) then
          l3 = 0
         else 
          l3 = 1     
         end if 
         if (l4 .eq. 0) then
          l4 = 1
         else 
          l4 = 0     
         end if 
        if (l5 .gt. nm) then
          l5 = 0
         else 
          l5 = 1     
         end if 
    
              
         ssum1 = ssum1 + l5*l3*dcmplx(y(abs(n)),
     &    sign(1,n)*y(abs(n)+nm))*dcmplx(y(abs(n+m)),
     &    -sign(1,n+m)*y(abs(n+m)+nm))/(8d0*(abs(m)+l2)*(abs(n)+l2p)*
     &    (abs(n+m+l4)))*(2d0*abs(n+m)+abs(m))
 
       enddo
       sum2(m) = ssum1
      enddo
      
Cc$omp end do
Cc$omp end parallel




      
C
      sum1 = dcmplx(0.0d0,0.0d0)
        
      do i = 1, nm     
      sum1= sum1+dcmplx(y(i),y(i+nm))*dcmplx(y(i),-y(i+nm))/i
      enddo
C
	  do i = 1,nm
	  V(i) = 0.5d0*dcmplx(y(i),-y(i+nm))/(i**2)
     &        -0.5d0*dcmplx(y(i),-y(i+nm))*sum1/i 
     &        + sum2(i) 

	  enddo  
	  
 

C      
Cc$omp parallel  
Cc$omp do
      do m=-nm,nm
      k=m
      if (k.eq.0) then
      k =0
      else
      k=1
      endif
      resd2=dcmplx(0d0,0d0)
      do n=1,nm
      resd2 = resd2+1d0*k*abs(n)*(-beta(n,m)*B(-n)-beta(-n,m)*B(n)+
     &  gamma(n,m)*Bd(-n)+gamma(-n,m)*Bd(n)+gammad(n,m)*B(-n)+
     & gammad(-n,m)*B(n))
      enddo
C      res1(m) = 0.5d0*dreal(resd2)
C      res3(m) = 0.5d0*dimag(resd2)
      res4(m) = resd2
      enddo 
Cc$omp end do
Cc$omp end parallel

      do i =1,nm
      res2(i) = res4(i) + 0.5d0*V(i)
      res2(-nm+i-1) = dconjg(res2(i))
      enddo 
   
      do i=1,nm
      rn(i)=dreal(res2(i))
      rn(i+nm) =dimag(res2(i))
      enddo
      
C     
   
      do i =1,nm
      do j=1,nm


      Am(i,j) =   dreal(q(i,j)+q(j,i)+q(i,-j)+q(-j,i))
      Am(i,nm+j) = - dimag(q(i,j)+q(j,i)-q(i,-j)-q(-j,i))
      Am(nm+i,j) = dimag(q(i,j)+q(j,i)+q(i,-j)+q(-j,i))
      Am(nm+i,nm+j) = dreal(q(i,j)+q(j,i)-q(i,-j)-q(-j,i))

      enddo 
      enddo 

       Ao = Am
      
       call dgesv(2*nm, 1, Ao, 2*nm, ipiv2, rn, 2*nm, info )
       rn =-rn
C       print*, rn
       
       do i =1,2*nm
       RK4(i) = rn(i)
       enddo
C      print*, RK4
      
      
      do i =1,2*nm
      yn(i) = yo(i) + dt/6d0*(ypo(i)+2d0*(ypo(i)+dt*RK1(i)*0.5d0)+
     & 2d0*(ypo(i)+dt*RK2(i)*0.5d0)+ (ypo(i)+dt*RK3(i)))
      enddo
C      print*, ypo(1)

      do i =1,2*nm
      yn2(i) = yo(i) + dt/2d0*(ypo(i)+(ypo(i)+dt*RK1(i)))
      enddo
      
      
      
  
      do i = 1,nm

     
      fv(i) = dreal(Q5(i)-dt/6d0*res2o(i))
      fv(i+nm) = dimag(Q5(i) - dt/6d0*res2o(i))
      fv2(i) = dreal(Q5(i)-dt/2d0*res2o(i))
      fv2(i+nm) = dimag(Q5(i) - dt/2d0*res2o(i))       
      
      enddo
 
      
     
      Amo2=Amo
        call dgesv(2*nm, 1, Amo, 2*nm, ipiv, fv,2*nm, info ) 
      
C      print*, RK4
            
      do i =1,2*nm
      yn(2*nm+i) =fv(i)+dt/6d0*(2d0*RK2(i)+2d0*RK3(i) + RK4(i)) 
      enddo
      
      call dgesv(2*nm, 1, Amo2, 2*nm, ipiv, fv2,2*nm, info )
      
      do i =1,2*nm
      yn2(2*nm+i) =fv2(i)+dt/2d0*RK2(i) 
      enddo
      do i =1,neq
      ydiff(i)=abs(yn(i)-yn2(i))
      enddo 

C      print*, yn
C      print*, ydiff
      
C      do i =1,neq
      i = 2*nm+1
      if (ydiff(i).ge.1e-7) then 
      t=t-dt
      dt=dt/1.1d0
      print*, "*********Step size Decrease********"
      goto 100
      end if
      if (ydiff(i).le.1e-16) then 
       print*, "*********Step size Increase********"   
       t=t-dt  
      dt=1.1d0*dt
      goto 100
      end if     
      print*, dt
C      enddo 
C     Now we introduce some averaging 
      pi = 3.14159265358d0
C      if (mod(count,50).eq.0) then
CCCC    Simply throw out energy in highest modes
C      do i =1,nm
C      if (i.le.2*nm/3) then 
C      yn(i)  = yn(i)
C      yn(i+nm)=yn(i+nm)
C      yn(i+2*nm)=yn(i+2*nm) 
C      yn(i+3*nm) = yn(i+3*nm)
C      else
C      yn(i)  = 0d0
C      yn(i+nm)=0d0
C      yn(i+2*nm)=0d0
C      yn(i+3*nm) = 0d0
C      end if
C      enddo 
C      end if
C
CCCC    5 point smoothing = LHC (1976), from dommermuth/yue 1986)
C      do i =1,nm
C      yn(i)  = 0.125d0*(5d0+4d0*dcos(pi*i/nm)-dcos(pi*2d0*i/nm))*
C     & yn(i)
C      yn(i+nm)=0.125d0*(5d0+4d0*dcos(pi*i/nm)-
C     & dcos(pi*2d0*i/nm))*yn(i+nm)
C      yn(i+2*nm)=0.125d0*(5d0+4d0*dcos(pi*i/nm)-dcos(pi*2d0*i/nm))*
C     & yn(i+2*nm) 
C      yn(i+3*nm) = 0.125d0*(5d0+4d0*dcos(pi*i/nm)-dcos(pi*2d0*i/nm))*
C     & yn(i+3*nm)
C      enddo
CC      
C      end if

C      print*, yn(1)
C     Compute energy
C      Energy=0d0
CC      c2=-ainit-yn(1)**2-0.5d0*yn(2)**2
C      c2=c**2
C      Energy=0.5d0*c2**2+c2*ainit+0.5d0*ainit**2+yn(1)**2+
C     & 2d0*c2*yn(1)**2+
C     & ainit*yn(1)**2+yn(2)*yn(1)**2+0.25d0*yn(2)**2+c2*yn(2)+
C     & ainit*yn(2) 
     
C      print*, Energy
      WRITE(fname,'(A,I5.5,A)') 'SW1024smooth_',count+1,'.txt'
       OPEN(14,file=fname,form='formatted')
       WRITE (14,*) (yn(i), i = 1,neq)
       close(unit =14) 
 
       WRITE(fname,'(A,I5.5,A)') 't1024smooth_',count,'.txt'
       OPEN(24,file=fname,form='formatted')
       WRITE (24,*) (t)
       close(unit =24)
           
       print*, count
       
       err=0d0
       do i =1,nm
       err=err+y(i)**2+y(i+nm)**2
C       sqrt(abs(yinitial(i)**2+yinitial(i+nm)**2-
C     & (yn(i)**2+ yn(i+nm)**2)))
       enddo
C      print*,err
      if (dt.le.1e-07) then
      print*, "STEP SIZE TOO SMALL"
      stop 
      endif
  

C       
C       
C         
      enddo     



      
C      return
      end

