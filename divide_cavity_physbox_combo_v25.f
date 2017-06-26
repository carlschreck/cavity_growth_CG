      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!
      !!
      !!   Quasi-static growth of cells in rectanular cavity with
      !!   physical outlet via conjugate gradient energy minimization. 
      !! 
      !! 
      !!       Cell type - 1: ellipse
      !!                   2: budding
      !!                   3: disk
      !!
      !!   Division type - 1: -> ->
      !!                   2: -> <-
      !!                   3: <- ->
      !!                   4: random
      !!
      !!        Feedback - growth rate ~ e^(-P/P0)
      !!                   P0=-1: no feedback
      !!
      !!
      !!      Carl  Schreck
      !!      9/17/2015
      !!
      !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      program divide_cavity_physbox_combo

      implicit none
      integer Ntot,Ngen
      parameter(Ntot=8192,Ngen=2000)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,exp,ran2
      double precision ftol,ftol1,fret,alpha0,width,Lx,Ly,ratei,ap
      double precision alpha(Ntot),rate(Ntot),alphar(Ntot),scale(Ntot)
      double precision rate0,desync,phi,flow,tdiv,P,PP(Ntot),D0(Ntot)
      double precision aclone(Ntot*Ngen),dispcm,xa(2),ya(2),PR,PT,P0
      double precision tbirth(Ntot*Ngen),cc,ss,dr(2),dd,corr,att,rat
      double precision xbirth(Ntot*Ngen),ybirth(Ntot*Ngen),distrem
      double precision xdiv(999,2),ydiv(999,2),thdiv(999,2)
      integer N,Nr,seed,iter,steps,i,j,k,kk,c(Ntot,Ngen),m,skip,Nexist
      integer celltype,divtype,nclone(Ntot*Ngen),nclonebox(Ntot*Ngen)
      integer age(Ntot),agehist(Ngen,2),agetot1,agetot2
      integer div,ndiv,idiv(999,2),kdiv
      character file1*80,file2*80,file3*80
      character file4*80,file5*80,file7*80,file8*80
      common /f2com/ width
      common /f3com/ alpha
      common /f4com/ exp,att
      common /f5com/ Lx,Ly,ap
      common /f6com/ P,PP,PT,PR
      common /f8com/ alpha0
      common /f9com/ scale
      common /f10com/ celltype

      ! read geometric parameters
      read(*,*) alpha0
      read(*,*) Lx
      read(*,*) Ly
      read(*,*) ap

      ! read cell parameters
      read(*,*) celltype
      read(*,*) divtype
      read(*,*) P0
      read(*,*) att

      ! read run parameters
      read(*,*) rate0
      read(*,*) desync
      read(*,*) steps
      read(*,*) seed
      read(*,*) distrem
      read(*,*) skip

      ! read output files
      read(*,*) file1
      read(*,*) file2
      read(*,*) file3
      read(*,*) file4
      read(*,*) file5
      read(*,*) file7
      read(*,*) file8

      ! parameters
      D1=1d0       ! Minor axis of particle
      exp=2d0      ! 2 =  LS, 2.5 = Hertzian, >2.9 = RLJ
      width=0.1d0  ! width of neighborlist 
      ftol=1d-16   ! Condition 1 for frprmn: V/N < ftol 
      ftol1=1d-16  ! Condition 2 for frprmn: dV/N < ftol1
      rat=1.5d0    ! ratio of initial cell 2 to cell 1 volume

      open(unit=1,file=TRIM(file1))
      open(unit=2,file=TRIM(file2))
      open(unit=3,file=TRIM(file3))
      open(unit=4,file=TRIM(file4))
      open(unit=5,file=TRIM(file5))
      open(unit=7,file=TRIM(file7))
      open(unit=8,file=TRIM(file8))

      ! initial size & aspect ratios - total vol = (1+rat)*alpha0
      if(celltype.eq.1) then
         d(1)=1d0
         d(2)=1d0
         alpha(1)=alpha0
         alpha(2)=alpha0*rat
      elseif(celltype.eq.2) then
         d(1)=1d0
         d(2)=1d0
         alpha(1)=alpha0
         alpha(2)=dsqrt(alpha0*(1d0+rat)-(alpha0-1d0)**2-2d0)+1d0
      elseif(celltype.eq.3) then
         d(1)=alpha0
         d(2)=dsqrt(alpha0*(1d0+rat)-1d0)
         alpha(1)=1d0
         alpha(2)=1d0
      endif      

      ! random initial config
      N=2
      Nexist=2
      do i=1,N 
         c(i,1)=1
         c(i,2)=i
         x(i)=0d0
         y(i)=(dble(i)-1.5d0)*d(i)*D1
         th(i)=(ran2(seed)-0.5d0)*2d0*pi
         D0(i)=D1
         rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
         age(i)=0
      enddo

      ! calculate # steps until division
      tdiv=dlog10(2d0)/dlog10(1d0+rate0)

      write(*,*) steps, tdiv, 2*int(steps/tdiv)

      ! set histogram of cell ages to 0
      do k=1,2*int(steps/tdiv)
         agehist(k,1)=0
         agehist(k,2)=0
      enddo
      
      ! loop over division
      do k=1,steps 
         ! grow particles
         ndiv=0
         do i=1,N
            if(x(i).gt.-Lx/2d0) then
               ratei=rate(i)
               if(P0.gt.0d0.and.PP(i).gt.0d0) then
                  ratei=ratei*dexp(-PP(i)/P0)
               endif
               if(celltype.eq.1) then
                  alpha(i)=(1d0+ratei)*alpha(i)
               elseif(celltype.eq.2) then
                  alpha(i)=1d0+dsqrt((1d0+ratei)*
     +                 (1d0+(alpha(i)-1d0)**2)-1d0)
               elseif(celltype.eq.3) then
                  d(i)=dsqrt(1d0+ratei)*d(i)
               endif

               if(celltype.le.2.and.alpha(i).gt.2d0*alpha0) then
                  dispcm=alpha0/2d0
                  div=1
                  
                  ! divide into 2 - 1st assigned index N+1
                  N=N+1
                  c(N,1)=c(i,1)+1
                  do j=2,c(N,1)
                     c(N,j)=c(i,j)
                  enddo
                  c(N,c(N,1)+1)=Nexist+1                  
                  D(N)=D(i)
                  x(N)=x(i)+dispcm*dcos(th(i))
                  y(N)=y(i)+dispcm*dsin(th(i))
                  rate(N)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  alpha(N)=alpha0
                  age(N)=0
                  
                  ! divide into 2 - 1st assigned index i
                  c(i,1)=c(i,1)+1
                  c(i,c(i,1)+1)=Nexist+2
                  x(i)=x(i)-dispcm*dcos(th(i))
                  y(i)=y(i)-dispcm*dsin(th(i))
                  rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  alpha(i)=alpha0
                  age(i)=age(i)+1

                  ! initialize # & area of clones
                  nclone(Nexist+1)=0
                  nclone(Nexist+2)=0
                  aclone(Nexist+1)=0d0
                  aclone(Nexist+2)=0d0

                  ! keep track of clones in system
                  nclonebox(Nexist+1)=1
                  nclonebox(Nexist+2)=1
                  do j=2,c(i,1)
                     nclonebox(c(i,j))=nclonebox(c(i,j))+1
                  enddo

                  ! keep track of cell position at birth
                  tbirth(Nexist+1)=dble(k)/tdiv
                  tbirth(Nexist+2)=dble(k)/tdiv
                  xbirth(Nexist+1)=x(N)
                  xbirth(Nexist+2)=x(i)
                  ybirth(Nexist+1)=y(N)
                  ybirth(Nexist+2)=y(i)

                  Nexist=Nexist+2
               elseif(celltype.eq.3.and.D(i).gt.dsqrt(2d0)*D0(i)) then
                  dispcm=1d0/2d0*D0(i)
                  div=1
                  
                  ! divide into 2 - 1st assigned index N+1
                  N=N+1
                  c(N,1)=c(i,1)+1
                  do j=2,c(N,1)
                     c(N,j)=c(i,j)
                  enddo
                  c(N,c(N,1)+1)=Nexist+1                  
                  x(N)=x(i)+dispcm*dcos(th(i))
                  y(N)=y(i)+dispcm*dsin(th(i))
                  rate(N)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  alpha(N)=alpha(i)
                  D(N)=D0(i)
                  D0(N)=D0(i)
                  
                  ! divide into 2 - 1st assigned index i
                  c(i,1)=c(i,1)+1
                  c(i,c(i,1)+1)=Nexist+2
                  x(i)=x(i)-dispcm*dcos(th(i))
                  y(i)=y(i)-dispcm*dsin(th(i))
                  rate(i)=(1d0+(ran2(seed)-0.5d0)*desync)*rate0 
                  D(i)=D0(i)
                  
                  ! initialize # & area of clones
                  nclone(Nexist+1)=0
                  nclone(Nexist+2)=0
                  aclone(Nexist+1)=0d0
                  aclone(Nexist+2)=0d0

                  ! keep track of clones in system
                  nclonebox(Nexist+1)=1
                  nclonebox(Nexist+2)=1
                  do j=2,c(i,1)
                     nclonebox(c(i,j))=nclonebox(c(i,j))+1
                  enddo

                  ! keep track of cell position at birth
                  tbirth(Nexist+1)=dble(k)/tdiv
                  tbirth(Nexist+2)=dble(k)/tdiv
                  xbirth(Nexist+1)=x(N)
                  xbirth(Nexist+2)=x(i)
                  ybirth(Nexist+1)=y(N)
                  ybirth(Nexist+2)=y(i)
                  
                  Nexist=Nexist+2
               else
                  div=0   
               endif
               
               ! types of division: ->->, <-->, -><-, random
               if(div.eq.1) then
                  if(divtype.eq.1) then
                     th(N)=th(i)
                  elseif(divtype.eq.2) then
                     th(N)=th(i)+pi
                  elseif(divtype.eq.3) then
                     th(N)=th(i)
                     th(i)=th(i)+pi
                  elseif(divtype.eq.4) then
                     th(N)=(ran2(seed)-0.5d0)*2d0*pi
                     th(i)=(ran2(seed)-0.5d0)*2d0*pi
                  endif

                  ! temp, remove
                  th(i)=th(i) + 1d-4*(ran2(seed)-0.5d0)
                  th(N)=th(N) + 1d-4*(ran2(seed)-0.5d0)
               endif

               if(div.eq.1) then
                  ndiv=ndiv+1
                  idiv(ndiv,1)=i
                  idiv(ndiv,2)=N
                  xdiv(ndiv,1)=x(i)
                  xdiv(ndiv,2)=x(N)
                  ydiv(ndiv,1)=y(i)
                  ydiv(ndiv,2)=y(N)
                  thdiv(ndiv,1)=th(i)
                  thdiv(ndiv,2)=th(N)
               endif
            endif
         enddo
         
         ! convert from angle to length scale = sqrt(I/m) * angle
         do i=1,N
            if(celltype.eq.1) then
               scale(i)=dsqrt(1d0+alpha(i)**2)/4d0*d(i)         
            elseif(celltype.eq.2) then
               dd=alpha(i)-1d0
               scale(i)=dsqrt(2d0*(1d0+dd**4)/(1+dd**2)+
     +              4d0*(dd*(1d0+dd)/(1+dd**2))**2)/4d0*d(i)
            elseif(celltype.eq.3) then
               scale(i)=dsqrt(2d0)/4d0*d(i)         
            endif
            th(i)=th(i)*scale(i)
         enddo

         ! minimize energy
         call frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
         write(*,*) k,N,fret/dble(N)
         
         ! convert back to angles
         do i=1,N
            th(i)=th(i)/scale(i)
         enddo

         ! output division time & displacement of cells
         do kdiv=1,ndiv
           write(8,'(I,7E)')k,dble(k)/tdiv,x(idiv(ndiv,1))-xdiv(ndiv,1), 
     +      y(idiv(ndiv,1))-ydiv(ndiv,1),th(idiv(ndiv,1))-thdiv(ndiv,1), 
     +      xdiv(ndiv,1),ydiv(ndiv,1),thdiv(ndiv,1)
           write(8,'(I,7E)')k,dble(k)/tdiv,x(idiv(ndiv,2))-xdiv(ndiv,2), 
     +      y(idiv(ndiv,2))-ydiv(ndiv,2),th(idiv(ndiv,2))-thdiv(ndiv,2),
     +      xdiv(ndiv,2),ydiv(ndiv,2),thdiv(ndiv,2)
         enddo
         flush(8)

         ! remove particle in outlet & outside box
         Nr=0d0
         do i=N,1,-1
            if(x(i).lt.-Lx/2d0-distrem) then
               do j=2,c(i,1)+1
                  m=c(i,j)
                  nclone(m)=nclone(m)+1
                  if(celltype.eq.1.or.celltype.eq.3) then
                     aclone(m)=aclone(m)+alpha(i)*d(i)**2
                  elseif(celltype.eq.2) then
                     aclone(m)=aclone(m)+(1d0+(alpha(i)-1d0)**2)*d(i)**2
                  endif
                  nclonebox(m)=nclonebox(m)-1
                  if(nclonebox(m).eq.0) then
                     write(5,'(3I,4F)')m,nclonebox(m),nclone(m),
     +                    pi/4d0*D1**2*aclone(m),tbirth(m),
     +                    xbirth(m),ybirth(m)
                  endif
               enddo
               flush(5)

               if(k.ge.2*int(div)) then
                  agehist(age(i),1)=agehist(age(i),1)+1
               endif

               Nr=Nr+1
               alphar(Nr)=alpha(i)               
               do j=1,c(N,1)+1
                  c(i,j)=c(N,j)
               enddo                  
               D(i)=D(N)
               x(i)=x(N)
               y(i)=y(N)
               th(i)=th(N)
               rate(i)=rate(N)
               alpha(i)=alpha(N)
               age(i)=age(N)
               N=N-1               
            endif
         enddo

         ! write out pressure (bulk, right, top)
         write(2,'(4E)') dble(k)/tdiv, P, PR, PT
         flush(2)

         ! save config         
         if(mod(k,skip).eq.0.and.celltype.eq.1) then
            write(1,*) N
            do i=1,N              
               write(1,'(5E,I)') x(i), y(i), th(i), 
     +              d(i), alpha(i)*d(i), c(i,c(i,1)+1)
            enddo
         elseif(mod(k,skip).eq.0.and.celltype.eq.2) then
            write(1,*) 2*N
            do i=1,N              
               cc=dcos(th(i))
               ss=dsin(th(i))
               dd=alpha(i)-1d0
               dr(1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
               dr(2)=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
               do kk=1,2
                  xa(kk)=x(i)+dr(kk)*cc
                  ya(kk)=y(i)+dr(kk)*ss
               enddo
               write(1,'(3E,I)')xa(1),ya(1),d(i),c(i,c(i,1)+1)
               write(1,'(3E,I)')xa(2),ya(2),d(i)*dd,c(i,c(i,1)+1)
            enddo
         elseif(mod(k,skip).eq.0.and.celltype.eq.3) then
            write(1,*) N
            do i=1,N              
               write(1,'(3E,I)')x(i),y(i),d(i),c(i,c(i,1)+1)
            enddo
         endif
         flush(1)

         ! calculate volume fracion in box (excluding outlet)
         phi=0d0
         do i=1,N
            if(x(i).gt.-Lx/2d0) then
               if(celltype.eq.1) then
                  phi=phi+alpha(i)*D(i)**2
               elseif(celltype.eq.2) then
                  phi=phi+(1d0+(alpha(i)-1d0)**2)*D(i)**2
               elseif(celltype.eq.3) then
                  phi=phi+D(i)**2
               endif
            endif
         enddo
         phi=pi*D1*D1*phi/Lx/Ly/4d0
         write(3,*) dble(k)/tdiv, phi
         flush(3)

         ! calculate flow rate: d(# outflow)/dt & d(mass outflow)/dt
         flow=0d0
         do i=1,Nr
            if(celltype.eq.1) then
               flow=flow+alphar(i)*D(i)**2
            elseif(celltype.eq.2) then
               flow=flow+(1d0+(alphar(i)-1d0)**2)*D(i)**2
            elseif(celltype.eq.3) then
               flow=flow+d(i)**2
            endif
         enddo
         flow=pi*D1*D1*flow/4d0
         write(4,*) dble(k)/tdiv, flow*tdiv, dble(Nr)*tdiv
         flush(4)

         ! calculate age hist of cells in box
         if(mod(k,int(tdiv))) then
            do i=1,N
               if(k.ge.2*int(tdiv)) then
                  agehist(age(i),2)=agehist(age(i),2)+1
               endif
            enddo
         endif
      enddo
      
      do i=1,N
         do j=2,c(i,1)+1
            m=c(i,j)
            nclone(m)=nclone(m)+1
            if(celltype.eq.1.or.celltype.eq.3) then
               aclone(m)=aclone(m)+alpha(i)*d(i)**2
            elseif(celltype.eq.2) then
               aclone(m)=aclone(m)+(1d0+(alpha(i)-1d0)**2)*d(i)**2
            endif 
         enddo
      enddo

      do m=1,Nexist
         if(nclonebox(m).gt.0) then 
            if(nclonebox(m).ne.0) then
               write(5,'(3I,4F)')m,nclonebox(m),nclone(m),
     +              pi/4d0*D1**2*aclone(m),tbirth(m),
     +              xbirth(m),ybirth(m)
            endif
         endif
      enddo      
      flush(5)

      agetot1=0
      agetot2=0
      do k=1,2*int(steps/tdiv)
         agetot1=agetot1+agehist(k,1)
         agetot2=agetot2+agehist(k,2)
      enddo

      do k=1,2*int(steps/tdiv)
         if(agehist(k,1).gt.0.or.agehist(k,2).gt.0) then        
            write(7,'(I,2F)') k, dble(agehist(k,1))/dble(agetot1), 
     +           dble(agehist(k,2))/dble(agetot2)
         endif
      enddo     
      flush(7)
      
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine CG_check(N,x,y,xp,yp,maxdis)
      parameter(Ntot = 8192)
      double precision maxdis,x(Ntot),y(Ntot),xp(Ntot),yp(Ntot)
      integer N

      maxdis=0d0
      do i=1,N
         maxdis=max(dabs(x(i)-xp(i)),maxdis)
         maxdis=max(dabs(y(i)-yp(i)),maxdis)
      enddo
      maxdis=2d0*dsqrt(2d0*maxdis*maxdis)

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makelist(N,x,y,D,D1,xp,yp,countn,nl)
      parameter(Ntot = 8192)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot),D1
      integer countn(Ntot),nl(800,Ntot),N,celltype
      common /f10com/ celltype

      if(celltype.eq.1) then
         call makelist_ellipse(N,x,y,D,D1,xp,yp,countn,nl)
      elseif(celltype.eq.2) then
         call makelist_dimer(N,x,y,D,D1,xp,yp,countn,nl)
      elseif(celltype.eq.3) then
         call makelist_disk(N,x,y,D,D1,xp,yp,countn,nl)
      endif

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makelist_ellipse(N,x,y,D,D1,xp,yp,countn,nl)
      parameter(Ntot = 8192)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot),exp
      double precision D1,xij,yij,rij,dij,rijsq,alpha(Ntot),width,att
      integer countn(Ntot),nl(800,Ntot),N
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att

      do i=1,N
         countn(i)=0
      enddo
      
      do i=1,N-1
         do j=i+1,N
            xij=x(i)-x(j)
            yij=y(i)-y(j)
            rijsq=xij*xij+yij*yij
            dij=dsqrt((alpha(i)**2*D(i)**2+
     +                 alpha(j)**2*D(j)**2)/2d0)*D1 ! major axes
            if(rijsq.lt.(dij+(att+width)*D1)**2) then
               countn(i)=countn(i)+1
               nl(countn(i),i)=j
            end if
         enddo
      enddo
      
      do i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      enddo

      return
      end
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makelist_dimer(N,x,y,D,D1,xp,yp,countn,nl) 
      parameter(Ntot = 8192)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot)
      double precision D1,xij,yij,rij,dij,rijsq,alpha(Ntot),width
      double precision dd,dr1,dr2,dk2,di_up(Ntot),exp,att
      integer countn(Ntot),nl(800,Ntot),N
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att

      do i=1,N
         countn(i)=0
      enddo

      do i=1,N
         dd=alpha(i)-1d0
         dr1=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)/2d0
         dr2=-(1d0+dd)/(1d0+dd**2)*D(i)/2d0
         dk2=dd*D(i)
         di_up(i)=(dk2/2d0-dr2)*2d0
      enddo

      do i=1,N-1
         do j=i+1,N
            xij=x(i)-x(j)
            yij=y(i)-y(j)
            rijsq=xij*xij+yij*yij
            dij=(di_up(i)+di_up(j))/2d0   
            if(rijsq.lt.(dij+(att+width)*D1)**2) then
               countn(i)=countn(i)+1
               nl(countn(i),i)=j
            end if
         enddo
      enddo
      
      do i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      enddo

      return
      end
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine makelist_disk(N,x,y,D,D1,xp,yp,countn,nl) 
      parameter(Ntot = 8192)
      double precision x(Ntot),y(Ntot),xp(Ntot),yp(Ntot),D(Ntot)
      double precision D1,xij,yij,rij,dij,rijsq,width,exp,att
      integer countn(Ntot),nl(800,Ntot),N
      common /f2com/ width
      common /f4com/ exp,att

      do i=1,N
         countn(i)=0
      enddo

      do i=1,N-1
         do j=i+1,N
            xij=x(i)-x(j)
            yij=y(i)-y(j)
            rijsq=xij*xij+yij*yij
            dij=D1*(d(i)+d(j))/2d0
            if(rijsq.lt.(dij+(att+width)*D1)**2) then
               countn(i)=countn(i)+1
               nl(countn(i),i)=j
            end if
         enddo
      enddo
      
      do i=1,N
         xp(i)=x(i)
         yp(i)=y(i)
      enddo

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine func(N,x,y,th,D,D1,V,countn,nl)
      parameter(Ntot = 8192)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V
      integer countn(Ntot),nl(800,Ntot),N,celltype
      common /f10com/ celltype

      if(celltype.eq.1) then
         call func_ellipse(N,x,y,th,D,D1,V,countn,nl)
      elseif(celltype.eq.2) then
         call func_dimer(N,x,y,th,D,D1,V,countn,nl)
      elseif(celltype.eq.3) then
         call func_disk(N,x,y,D,D1,V,countn,nl)
      endif

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dfunc(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      parameter(Ntot = 8192)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1
      double precision fx(Ntot),fy(Ntot),fth(Ntot)
      integer countn(Ntot),nl(800,Ntot),N,celltype
      common /f10com/ celltype

      if(celltype.eq.1) then
         call dfunc_ellipse(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      elseif(celltype.eq.2) then
         call dfunc_dimer(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      elseif(celltype.eq.3) then
         call dfunc_disk(N,x,y,D,D1,fx,fy,fth,countn,nl)
      endif

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine func_ellipse(N,x,y,th,D,D1,V,countn,nl)
      parameter(Ntot = 8192)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V,LJ
      double precision rij,xij,yij,dij,exp,dlnsig,dij_up,sigma
      double precision Lx,Ly,ap,scale(Ntot),rijsq,alpha(Ntot)
      double precision fthi,fthj,fth_c,f_x,f_y,ft,fc,fr,att
      integer countn(Ntot),nl(800,Ntot),N,jj,i,j
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly,ap
      common /f9com/ scale

      ! inter-particle interactions
      V=0d0
      do i=1,N
         if(countn(i).ge.1) then
            do jj=1,countn(i)
               j=nl(jj,i)

               dij_up=D1*dsqrt((alpha(i)**2*D(i)**2+
     +                           alpha(j)**2*D(j)**2)/2d0)
               xij=x(i)-x(j)
               if(dabs(xij).lt.dij_up+att) then
                  yij=y(i)-y(j)
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij_up+att)**2) then 
                     rij=dsqrt(rijsq)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    th(j)/scale(j),D(i),D(j),alpha(i),alpha(j))*D1 
                     if(rij.lt.dij+att) then 
                        if(exp .gt. 2.9) then
                           LJ = (dij/rij)*(dij/rij)
                           LJ = LJ*LJ*LJ
                           V=V+(LJ-1d0)*(LJ-1d0)
                        else
                           V=V+(1d0-rij/dij)**exp/exp-(att/dij)**exp/exp
                        endif
                     end if
                  endif
               endif
            enddo
         end if
      enddo

      ! top wall
      do i=1,N
         dij_up=alpha(i)*D(i)*D1
         yij=y(i)-(Ly-y(i))
         if(dabs(yij).lt.dij_up) then 
            xij=0d0
            rij=dabs(yij)
            if(rij.lt.dij_up) then 
               dij=sigma(rij,xij,yij,th(i)/scale(i),-th(i)/scale(i),
     +              D(i),D(i),alpha(i),alpha(i))*D1 
               if(rij.lt.dij) then
                  if(exp .gt. 2.9) then
                     LJ = (dij/rij)*(dij/rij)
                     LJ = LJ*LJ*LJ
                     V=V+(LJ-1d0)*(LJ-1d0) / 2d0
                  else
                     V=V+(1d0-rij/dij)**exp/exp / 2d0
                  end if
               end if
            endif
         endif
      enddo

      ! bottom wall
      do i=1,N
         dij_up=alpha(i)*D(i)*D1
         yij=y(i)-(-Ly-y(i))
         if(dabs(yij).lt.dij_up) then 
            xij=0d0
            rij=dabs(yij)
            if(rij.lt.dij_up) then 
               dij=sigma(rij,xij,yij,th(i)/scale(i),-th(i)/scale(i),
     +              D(i),D(i),alpha(i),alpha(i))*D1
               if(rij.lt.dij) then
                  if(exp .gt. 2.9) then
                     LJ = (dij/rij)*(dij/rij)
                     LJ = LJ*LJ*LJ
                     V=V+(LJ-1d0)*(LJ-1d0) / 2d0
                  else
                     V=V+(1d0-rij/dij)**exp/exp / 2d0
                  end if
               endif
            end if
         endif
      enddo

      ! right wall
      do i=1,N
         dij_up=alpha(i)*D(i)*D1
         xij=x(i)-(Lx-x(i))
         if(dabs(xij).lt.dij_up) then 
            yij=0d0
            rij=dabs(xij)
            if(rij.lt.dij_up) then 
               dij=sigma(rij,xij,yij,th(i)/scale(i),pi-th(i)/scale(i),
     +              D(i),D(i),alpha(i),alpha(i))*D1
               if(rij.lt.dij) then
                  if(exp .gt. 2.9) then
                     LJ = (dij/rij)*(dij/rij)
                     LJ = LJ*LJ*LJ
                     V=V+(LJ-1d0)*(LJ-1d0) / 2d0
                  else
                     V=V+(1d0-rij/dij)**exp/exp / 2d0
                  end if
               endif
            end if
         endif
      enddo

      ! left wall with outlet
      do i=1,N
         dij_up=alpha(i)*D(i)*D1
         xij=x(i)-(-Lx-x(i))
         if(xij.lt.dij_up) then
            yij=0d0
            rij=dabs(xij)
            dij=sigma(rij,xij,yij,th(i)/scale(i),pi-th(i)/scale(i),
     +           D(i),D(i),alpha(i),alpha(i))*D1            
            if(xij.lt.dij) then ! if touching ghost wall
                  
               ! contact with upper boundary
               xij=x(i)+Lx/2d0
               yij=y(i)-AP/2d0
               xij=2d0*xij
               yij=2d0*yij
               rij=dsqrt(xij*xij+yij*yij)
               call sigma2(rij,xij,yij,th(i)/scale(i),
     +              th(i)/scale(i),D(i),D(i),dij,ft,
     +              fthi,fthj,alpha(i),alpha(i))
               dij=dij*D1
               if(exp .gt. 2.9) then
                  LJ=(dij/rij)*(dij/rij)
                  LJ=LJ*LJ*LJ
                  fc=1d0/rij*LJ*(LJ-1d0)
               else
                  fc=(1d0-rij/dij)**(exp-1d0)/dij     
               endif
               fr=-fc/rij
               fth_c=rij*fc
               f_x=fr*(xij+yij*ft)
               f_y=fr*(yij-xij*ft)               
               if(rij.lt.dij) then
                  if(f_x.lt.0d0.and.f_y.gt.0d0) then ! corner
                     V=V+(1d0-rij/dij)**exp/exp / 2d0
                  elseif(f_x+f_y.gt.0d0) then ! outlet         
                     xij=0d0
                     yij=y(i)-(AP-y(i))
                     rij=dabs(yij)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    -th(i)/scale(i),D(i),D(i),
     +                    alpha(i),alpha(i))*D1
                     V=V+(1d0-rij/dij)**exp/exp / 2d0
                  else  ! wall
                     xij=x(i)-(-Lx-x(i))
                     yij=0d0
                     rij=dabs(xij)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    pi-th(i)/scale(i),D(i),D(i),
     +                    alpha(i),alpha(i))*D1
                     V=V+(1d0-rij/dij)**exp/exp / 2d0       
                  endif
               else
                  if(f_x.lt.0d0.and.f_x+f_y.lt.0d0) then ! outlet       
                     xij=0d0
                     yij=y(i)-(AP-y(i))
                     rij=dabs(yij)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    -th(i)/scale(i),D(i),D(i),
     +                    alpha(i),alpha(i))*D1
                     
                     if(rij.lt.dij) then
                        V=V+(1d0-rij/dij)**exp/exp / 2d0
                     endif                     
                  elseif(f_y.gt.0d0.and.f_x+f_y.gt.0d0) then ! wall     
                     xij=x(i)-(-Lx-x(i))
                     yij=0d0
                     rij=dabs(xij)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    pi-th(i)/scale(i),D(i),D(i),
     +                    alpha(i),alpha(i))*D1                     
                     if(rij.lt.dij) then
                        V=V+(1d0-rij/dij)**exp/exp / 2d0          
                     endif
                  endif
               endif
               
               ! contact with lower boundary
               xij=x(i)+Lx/2d0
               yij=y(i)+AP/2d0
               xij=2d0*xij
               yij=2d0*yij
               rij=dsqrt(xij*xij+yij*yij)
               call sigma2(rij,xij,yij,th(i)/scale(i),
     +              th(i)/scale(i),D(i),D(i),dij,ft,
     +              fthi,fthj,alpha(i),alpha(i))
               dij=dij*D1
               if(exp .gt. 2.9) then
                  LJ=(dij/rij)*(dij/rij)
                  LJ=LJ*LJ*LJ
                  fc=1d0/rij*LJ*(LJ-1d0)
               else
                  fc=(1d0-rij/dij)**(exp-1d0)/dij     
               endif
               fr=-fc/rij
               fth_c=rij*fc
               f_x=fr*(xij+yij*ft)
               f_y=fr*(yij-xij*ft)
               if(rij.lt.dij) then
                  if(f_x.lt.0d0.and.f_y.lt.0d0) then ! corner
                     V=V+(1d0-rij/dij)**exp/exp / 2d0            
                  elseif(f_x-f_y.gt.0d0) then ! outlet
                     xij=0d0
                     yij=y(i)-(-AP-y(i))
                     rij=dabs(yij)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    -th(i)/scale(i),D(i),D(i),
     +                    alpha(i),alpha(i))*D1
                     V=V+(1d0-rij/dij)**exp/exp / 2d0    
                  else  ! wall
                     xij=x(i)-(-Lx-x(i))
                     yij=0d0
                     rij=dabs(xij)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    pi-th(i)/scale(i),D(i),D(i),
     +                    alpha(i),alpha(i))*D1
                     V=V+(1d0-rij/dij)**exp/exp / 2d0  
                  endif
               else      
                  if(f_x-f_y.lt.0d0.and.f_x.lt.0d0) then ! outlet      
                     xij=0d0
                     yij=y(i)-(-AP-y(i))
                     rij=dabs(yij)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    -th(i)/scale(i),D(i),D(i),
     +                    alpha(i),alpha(i))*D1                     
                     if(rij.lt.dij) then                        
                        V=V+(1d0-rij/dij)**exp/exp / 2d0 
                     endif
                  elseif(f_x-f_y.gt.0d0.and.f_y.lt.0d0) then ! wall  
                     xij=x(i)-(-Lx-x(i))
                     yij=0d0
                     rij=dabs(xij)
                     dij=sigma(rij,xij,yij,th(i)/scale(i),
     +                    pi-th(i)/scale(i),D(i),D(i),
     +                    alpha(i),alpha(i))*D1                     
                     if(rij.lt.dij) then
                        V=V+(1d0-rij/dij)**exp/exp / 2d0   
                     endif                     
                  endif
               endif
            endif
         endif
      enddo

      if(exp.gt.2.9) then
         V=V/72d0
      endif

      return				
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dfunc_ellipse(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      parameter(Ntot = 8192)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),sigma,D(Ntot),D1,dij
      double precision fx(Ntot),fy(Ntot),fth(Ntot),rij,xij,yij,fr,exp
      double precision dij_up,alpha(Ntot),LJ,fc,ft,f_x,f_y,scale(Ntot)
      double precision fthi,fthj,fth_c,Lx,Ly,ap,P,rijsq,PP(Ntot),att
      double precision Pij,PT,PR
      integer countn(Ntot),nl(800,Ntot),N,jj,i,j
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly,ap
      common /f6com/ P,PP,PT,PR
      common /f9com/ scale

      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
         PP(i)=0d0
      enddo
      P=0d0

      ! inter-particle interactions
      do i=1,N
         if(countn(i).ge.1) then
            do jj=1,countn(i)
               j=nl(jj,i)
               dij_up=D1*dsqrt((alpha(i)**2*D(i)**2+
     +              alpha(j)**2*D(j)**2)/2d0) ! modify for alpha < 1
               xij=x(i)-x(j)
               if(dabs(xij).lt.dij_up+att) then
                  yij=y(i)-y(j)                  
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij_up+att)**2) then
                     rij=dsqrt(rijsq)
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    th(j)/scale(j),D(i),D(j),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(j))
                     dij=dij*D1
                     if(rij.lt.dij+att) then
                        if(exp .gt. 2.9) then
                           LJ = (dij/rij)*(dij/rij)
                           LJ = LJ*LJ*LJ
                           fc = 1d0/rij*LJ*(LJ-1d0)
                        else
                           fc=(1d0-rij/dij)**(exp-1d0)/dij     
                        endif
                        fr=-fc/rij
                        fth_c=rij*fc
                        f_x=fr*(xij+yij*ft)
                        f_y=fr*(yij-xij*ft)
                        fx(i)=fx(i)+f_x
                        fy(i)=fy(i)+f_y
                        fth(i)=fth(i)+fthi*fth_c
                        fx(j)=fx(j)-f_x
                        fy(j)=fy(j)-f_y
                        fth(j)=fth(j)+fthj*fth_c

                        Pij=-xij*f_x-yij*f_y
                        P=P+2d0*Pij
                        PP(i)=PP(i)+Pij
                        PP(j)=PP(j)+Pij
                     end if
                  end if
               end if
            enddo
         end if
      enddo

      ! top wall
      PT=0d0
      do i=1,N
         dij_up=alpha(i)*D(i)
         yij=y(i)-(Ly-y(i))
         if(dabs(yij).lt.dij_up) then 
            xij=0d0
            rij=dabs(yij)
            if(rij.lt.dij_up) then 
               call sigma2(rij,xij,yij,th(i)/scale(i),-th(i)/scale(i),
     +              D(i),D(i),dij,ft,fthi,fthj,alpha(i),alpha(i))
               dij=dij*D1
               if(rij.lt.dij) then
                  if(exp .gt. 2.9) then
                     LJ = (dij/rij)*(dij/rij)
                     LJ = LJ*LJ*LJ
                     fc = 1d0/rij*LJ*(LJ-1d0)
                  else
                     fc = (1d0-rij/dij)**(exp-1d0)/dij     
                  endif
                  fr=-fc/rij
                  fth_c=rij*fc
                  f_x=fr*(xij+yij*ft)
                  f_y=fr*(yij-xij*ft)
                  fx(i)=fx(i)+f_x
                  fy(i)=fy(i)+f_y
                  fth(i)=fth(i)+fthi*fth_c
                  
                  PT=PT+f_y
                  Pij=-xij*f_x-yij*f_y
                  P=P+Pij
                  PP(i)=PP(i)+Pij
               endif
            endif
         endif
      enddo

      ! bottom wall
      do i=1,N
         dij_up=alpha(i)*D(i)*D1
         yij=y(i)-(-Ly-y(i))
         if(dabs(yij).lt.dij_up) then 
            xij=0d0
            rij=dabs(yij)
            if(rij.lt.dij_up) then                
               call sigma2(rij,xij,yij,th(i)/scale(i),-th(i)/scale(i),
     +              D(i),D(i),dij,ft,fthi,fthj,alpha(i),alpha(i))
               dij=dij*D1
               if(rij.lt.dij) then
                  if(exp .gt. 2.9) then
                     LJ = (dij/rij)*(dij/rij)
                     LJ = LJ*LJ*LJ
                     fc = 1d0/rij*LJ*(LJ-1d0)
                  else
                     fc = (1d0-rij/dij)**(exp-1d0)/dij
                  endif
                  fr=-fc/rij
                  fth_c=rij*fc
                  f_x=fr*(xij+yij*ft)
                  f_y=fr*(yij-xij*ft)
                  fx(i)=fx(i)+f_x
                  fy(i)=fy(i)+f_y
                  fth(i)=fth(i)+fthi*fth_c

                  Pij=-xij*f_x-yij*f_y
                  P=P+Pij
                  PP(i)=PP(i)+Pij
               endif
            end if
         endif
      enddo

      ! right wall
      PR=0d0
      do i=1,N
         dij_up=alpha(i)*D(i)*D1
         xij=x(i)-(Lx-x(i))
         if(dabs(xij).lt.dij_up) then 
            yij=0d0
            rij=dabs(xij)
            if(rij.lt.dij_up) then 
               call sigma2(rij,xij,yij,th(i)/scale(i),pi-th(i)/scale(i),
     +              D(i),D(i),dij,ft,fthi,fthj,alpha(i),alpha(i))
               dij=dij*D1
               if(rij.lt.dij) then
                  if(exp .gt. 2.9) then
                     LJ = (dij/rij)*(dij/rij)
                     LJ = LJ*LJ*LJ
                     fc = 1d0/rij*LJ*(LJ-1d0)
                  else
                     fc = (1d0-rij/dij)**(exp-1d0)/dij     
                  endif
                  fr=-fc/rij
                  fth_c=rij*fc
                  f_x=fr*(xij+yij*ft)
                  f_y=fr*(yij-xij*ft)
                  fx(i)=fx(i)+f_x
                  fy(i)=fy(i)+f_y
                  fth(i)=fth(i)+fthi*fth_c    
              
                  PR=PR+f_x
                  Pij=-xij*f_x-yij*f_y
                  P=P+Pij
                  PP(i)=PP(i)+Pij
               endif
            end if
         endif
      enddo

      !  left wall with outlet
      do i=1,N
         dij_up=alpha(i)*D(i)*D1
         xij=x(i)-(-Lx-x(i))         
         if(xij.lt.dij_up) then
            yij=0d0
            rij=dabs(xij)
            dij=sigma(rij,xij,yij,th(i)/scale(i),pi-th(i)/scale(i),
     +           D(i),D(i),alpha(i),alpha(i))*D1               
            if(xij.lt.dij) then ! if touching ghost wall

               ! contact with upper boundary
               xij=x(i)+Lx/2d0
               yij=y(i)-AP/2d0
               xij=2d0*xij
               yij=2d0*yij
               rij=dsqrt(xij*xij+yij*yij)
               call sigma2(rij,xij,yij,th(i)/scale(i),
     +              th(i)/scale(i),D(i),D(i),dij,ft,
     +              fthi,fthj,alpha(i),alpha(i))
               dij=dij*D1
               if(exp .gt. 2.9) then
                  LJ=(dij/rij)*(dij/rij)
                  LJ=LJ*LJ*LJ
                  fc=1d0/rij*LJ*(LJ-1d0)
               else
                  fc=(1d0-rij/dij)**(exp-1d0)/dij     
               endif
               fr=-fc/rij      
               fth_c=rij*fc
               f_x=fr*(xij+yij*ft)
               f_y=fr*(yij-xij*ft)
               if(rij.lt.dij) then
                  if(f_x.lt.0d0.and.f_y.gt.0d0) then ! corner     
                     fx(i)=fx(i)+f_x
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+fthi*fth_c  

                     Pij=-xij*f_x-yij*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij                   
                  elseif(f_x+f_y.gt.0d0) then ! outlet         
                     xij=0d0
                     yij=y(i)-(AP-y(i))
                     rij=dabs(yij)                     
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    -th(i)/scale(i),D(i),D(i),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(i))
                     dij=dij*D1
                     if(exp .gt. 2.9) then
                        LJ=(dij/rij)*(dij/rij)
                        LJ=LJ*LJ*LJ
                        fc=1d0/rij*LJ*(LJ-1d0)
                     else
                        fc=(1d0-rij/dij)**(exp-1d0)/dij     
                     endif
                     fr=-fc/rij 
                     fth_c=rij*fc
                     f_x=fr*(xij+yij*ft)
                     f_y=fr*(yij-xij*ft)                        
                     fx(i)=fx(i)+f_x
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+fthi*fth_c  

                     Pij=-xij*f_x-yij*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij                
                  else ! wall
                     xij=x(i)-(-Lx-x(i))
                     yij=0d0
                     rij=dabs(xij)                     
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    pi-th(i)/scale(i),D(i),D(i),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(i))
                     dij=dij*D1 
                     if(exp .gt. 2.9) then
                        LJ=(dij/rij)*(dij/rij)
                        LJ=LJ*LJ*LJ
                        fc=1d0/rij*LJ*(LJ-1d0)
                     else
                        fc=(1d0-rij/dij)**(exp-1d0)/dij     
                     endif
                     fr=-fc/rij 
                     fth_c=rij*fc
                     f_x=fr*(xij+yij*ft)
                     f_y=fr*(yij-xij*ft)                        
                     fx(i)=fx(i)+f_x
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+fthi*fth_c    

                     Pij=-xij*f_x-yij*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij               
                  endif
               else
                  if(f_x.lt.0d0.and.f_x+f_y.lt.0d0) then ! outlet   
                     xij=0d0
                     yij=y(i)-(AP-y(i))
                     rij=dabs(yij)                     
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    -th(i)/scale(i),D(i),D(i),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(i))
                     dij=dij*D1                    
                     if(rij.lt.dij) then
                        if(exp .gt. 2.9) then
                           LJ=(dij/rij)*(dij/rij)
                           LJ=LJ*LJ*LJ
                           fc=1d0/rij*LJ*(LJ-1d0)
                        else
                           fc=(1d0-rij/dij)**(exp-1d0)/dij     
                        endif
                        fr=-fc/rij 
                        fth_c=rij*fc
                        f_x=fr*(xij+yij*ft)
                        f_y=fr*(yij-xij*ft)                    
                        fx(i)=fx(i)+f_x
                        fy(i)=fy(i)+f_y
                        fth(i)=fth(i)+fthi*fth_c

                        Pij=-xij*f_x-yij*f_y
                        P=P+Pij
                        PP(i)=PP(i)+Pij
                     endif
                  elseif(f_y.gt.0d0.and.f_x+f_y.gt.0d0) then ! wall     
                     xij=x(i)-(-Lx-x(i))
                     yij=0d0
                     rij=dabs(xij)
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    pi-th(i)/scale(i),D(i),D(i),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(i))
                     dij=dij*D1                     
                     if(rij.lt.dij) then
                        if(exp .gt. 2.9) then
                           LJ=(dij/rij)*(dij/rij)
                           LJ=LJ*LJ*LJ
                           fc=1d0/rij*LJ*(LJ-1d0)
                        else
                           fc=(1d0-rij/dij)**(exp-1d0)/dij     
                        endif
                        fr=-fc/rij 
                        fth_c=rij*fc
                        f_x=fr*(xij+yij*ft)
                        f_y=fr*(yij-xij*ft)                     
                        fx(i)=fx(i)+f_x
                        fy(i)=fy(i)+f_y
                        fth(i)=fth(i)+fthi*fth_c      
                        
                        Pij=-xij*f_x-yij*f_y
                        P=P+Pij
                        PP(i)=PP(i)+Pij
                     endif                     
                  endif
               endif
               
               ! contact with lower boundary
               xij=x(i)+Lx/2d0
               yij=y(i)+AP/2d0
               xij=2d0*xij
               yij=2d0*yij
               rij=dsqrt(xij*xij+yij*yij)
               call sigma2(rij,xij,yij,th(i)/scale(i),
     +              th(i)/scale(i),D(i),D(i),dij,ft,
     +              fthi,fthj,alpha(i),alpha(i))
               dij=dij*D1
               if(exp .gt. 2.9) then
                  LJ=(dij/rij)*(dij/rij)
                  LJ=LJ*LJ*LJ
                  fc=1d0/rij*LJ*(LJ-1d0)
               else
                  fc=(1d0-rij/dij)**(exp-1d0)/dij     
               endif
               fr=-fc/rij
               fth_c=rij*fc
               f_x=fr*(xij+yij*ft)
               f_y=fr*(yij-xij*ft)               
               if(rij.lt.dij) then
                  if(f_x.lt.0d0.and.f_y.lt.0d0) then ! corner       
                     fx(i)=fx(i)+f_x
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+fthi*fth_c          

                     Pij=-xij*f_x-yij*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij             
                  elseif(f_x-f_y.gt.0d0) then ! outlet
                     xij=0d0
                     yij=y(i)-(-AP-y(i))
                     rij=dabs(yij)                      
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    -th(i)/scale(i),D(i),D(i),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(i))
                     dij=dij*D1                     
                     if(exp .gt. 2.9) then
                        LJ=(dij/rij)*(dij/rij)
                        LJ=LJ*LJ*LJ
                        fc=1d0/rij*LJ*(LJ-1d0)
                     else
                        fc=(1d0-rij/dij)**(exp-1d0)/dij     
                     endif
                     fr=-fc/rij 
                     fth_c=rij*fc 
                     f_x=fr*(xij+yij*ft)
                     f_y=fr*(yij-xij*ft)                             
                     fx(i)=fx(i)+f_x
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+fthi*fth_c    

                     Pij=-xij*f_x-yij*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij              
                  else ! wall
                     xij=x(i)-(-Lx-x(i))
                     yij=0d0
                     rij=dabs(xij)
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    pi-th(i)/scale(i),D(i),D(i),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(i))
                     dij=dij*D1
                     if(exp .gt. 2.9) then
                        LJ=(dij/rij)*(dij/rij)
                        LJ=LJ*LJ*LJ
                        fc=1d0/rij*LJ*(LJ-1d0)
                     else
                        fc=(1d0-rij/dij)**(exp-1d0)/dij     
                     endif
                     fr=-fc/rij 
                     fth_c=rij*fc 
                     f_x=fr*(xij+yij*ft)
                     f_y=fr*(yij-xij*ft)                   
                     fx(i)=fx(i)+f_x
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+fthi*fth_c    

                     Pij=-xij*f_x-yij*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij         
                  endif
               else       
                  if(f_x-f_y.lt.0d0.and.f_x.lt.0d0) then ! outlet  
                     xij=0d0
                     yij=y(i)-(-AP-y(i))
                     rij=dabs(yij)
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    -th(i)/scale(i),D(i),D(i),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(i))
                     dij=dij*D1                     
                     if(rij.lt.dij) then
                        if(exp .gt. 2.9) then
                           LJ=(dij/rij)*(dij/rij)
                           LJ=LJ*LJ*LJ
                           fc=1d0/rij*LJ*(LJ-1d0)
                        else
                           fc=(1d0-rij/dij)**(exp-1d0)/dij     
                        endif
                        fr=-fc/rij 
                        fth_c=rij*fc 
                        f_x=fr*(xij+yij*ft)
                        f_y=fr*(yij-xij*ft)                   
                        fx(i)=fx(i)+f_x
                        fy(i)=fy(i)+f_y
                        fth(i)=fth(i)+fthi*fth_c 

                        Pij=-xij*f_x-yij*f_y
                        P=P+Pij
                        PP(i)=PP(i)+Pij             
                     endif                    
                  elseif(f_x-f_y.gt.0d0.and.f_y.lt.0d0) then ! wall   
                     xij=x(i)-(-Lx-x(i))
                     yij=0d0
                     rij=dabs(xij)                     
                     call sigma2(rij,xij,yij,th(i)/scale(i),
     +                    pi-th(i)/scale(i),D(i),D(i),dij,ft,
     +                    fthi,fthj,alpha(i),alpha(i))
                     dij=dij*D1                     
                     if(rij.lt.dij) then
                        if(exp .gt. 2.9) then
                           LJ=(dij/rij)*(dij/rij)
                           LJ=LJ*LJ*LJ
                           fc=1d0/rij*LJ*(LJ-1d0)
                        else
                           fc=(1d0-rij/dij)**(exp-1d0)/dij     
                        endif
                        fr=-fc/rij 
                        fth_c=rij*fc 
                        f_x=fr*(xij+yij*ft)
                        f_y=fr*(yij-xij*ft)                    
                        fx(i)=fx(i)+f_x
                        fy(i)=fy(i)+f_y
                        fth(i)=fth(i)+fthi*fth_c                

                        Pij=-xij*f_x-yij*f_y
                        P=P+Pij
                        PP(i)=PP(i)+Pij
                     endif
                  endif
               endif
            endif
         endif
      enddo

      if(exp.gt.2.9) then
         do i=1, N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
            fth(i)=fth(i)/6d0 
         enddo
         P=P/6d0
      endif

      do i=1,N
         fth(i)=fth(i)/scale(i)
      enddo

      PT=PT/Lx
      PR=PR/Ly
      P=P/4d0/Lx/Ly
      do i=1,N
         PP(i)=PP(i)*dble(N)/4d0/Lx/Ly
      enddo

      return							
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine func_dimer(N,x,y,th,D,D1,V,countn,nl)
      parameter(Ntot = 8192)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),D(Ntot),D1,V,alpha(Ntot)
      double precision rij,xij,yij,dij,exp,dlnsig,dij_up,sigma,LJ
      double precision Lx,Ly,ap,rijsq,dijsq_up,scale(Ntot),c(Ntot),att
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),Vij
      double precision dk(Ntot,2),di_up(Ntot),di1j1,xhit,yhit,yhitout
      integer countn(Ntot),nl(800,Ntot),N
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly,ap
      common /f9com/ scale

      ! convert to from molecules to atoms
      do i=1,N
         c(i)=dcos(th(i)/scale(i))
         s(i)=dsin(th(i)/scale(i))
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
         if(alpha(i).lt.2d0) then 
            di_up(i)=(dk(i,2)/2d0-dr(i,2))*2d0*D1
         else
            di_up(i)=(dk(i,1)/2d0+dr(i,1))*2d0*D1         
         endif
      enddo

      ! inter-particle interactions
      V=0d0
      do i=1,N
         if(countn(i).ge.1) then
            do jj=1,countn(i)
               j=nl(jj,i)               
               dij_up=(di_up(i)+di_up(j))/2d0
               xij=x(i)-x(j)
               if(dabs(xij).lt.dij_up+att) then
                  yij=y(i)-y(j)                  
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij_up+att)**2) then
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     do ki=1,2
                        do kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           yij=ya(i,ki)-ya(j,kj)
                           rijsq=xij**2+yij**2
                           if(rijsq.lt.(dij+att)**2) then
                              rij=dsqrt(rijsq)
                              if(exp .gt. 2.9) then
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 Vij=(LJ-1d0)*(LJ-1d0)
                              else
                                 Vij=(1d0-rij/dij)**exp/exp-
     +                                (att/dij)**exp/exp
                              endif 
                              V=V+Vij*dij**2/di1j1**2
                           endif
                        enddo
                     enddo
                  end if
               end if
            enddo
         end if
      enddo

      ! boundary conditions
      do i=1,N
         do ki=1,2     
            xhit=(Lx-dk(i,ki))/2d0
            yhit=(Ly-dk(i,ki))/2d0
            if(ya(i,ki).gt.yhit) then ! top wall
               V=V+((ya(i,ki)-yhit)/dk(i,1))**exp/exp
            elseif(ya(i,ki).lt.-yhit) then ! bottom wall
               V=V+((ya(i,ki)+yhit)/dk(i,1))**exp/exp
            endif
            if(xa(i,ki).gt.xhit) then ! right wall
               V=V+((xa(i,ki)-xhit)/dk(i,1))**exp/exp
            elseif(xa(i,ki).lt.-xhit) then ! left wall, corner, outlet
               if(dabs(ya(i,ki))+min(xa(i,ki)+Lx/2d0,0d0).gt.AP/2d0)then
                  V=V+((xa(i,ki)+xhit)/dk(i,1))**exp/exp
               elseif(xa(i,ki).gt.-Lx/2d0) then   
                  rijsq=(xa(i,ki)+Lx/2d0)**2+(ya(i,ki)-AP/2d0)**2
                  if(rijsq.lt.dk(i,ki)**2/4d0) then
                     rij=dsqrt(rijsq)
                     V=V+((dk(i,ki)/2d0-rij)/dk(i,1))**exp/exp
                  endif
                  rijsq=(xa(i,ki)+Lx/2d0)**2+(ya(i,ki)+AP/2d0)**2
                  if(rijsq.lt.dk(i,ki)**2/4d0) then
                     rij=dsqrt(rijsq)
                     V=V+((dk(i,ki)/2d0-rij)/dk(i,1))**exp/exp
                  endif
               else                 
                  yhitout=(AP-dk(i,ki))/2d0
                  if(ya(i,ki).gt.yhitout) then
                     V=V+((ya(i,ki)-yhitout)/dk(i,1))**exp/exp
                  endif
                  if(ya(i,ki).lt.-yhitout) then
                     V=V+((ya(i,ki)+yhitout)/dk(i,1))**exp/exp
                  endif
               endif
            endif
         enddo
      enddo

      if(exp.gt.2.9) then
         V=V/72d0
      endif

      return				
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dfunc_dimer(N,x,y,th,D,D1,fx,fy,fth,countn,nl)
      parameter(Ntot = 8192)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),th(Ntot),sigma,D(Ntot),D1,dij
      double precision fx(Ntot),fy(Ntot),fth(Ntot),rij,xij,yij,fr,exp
      double precision dij_up,alpha(Ntot),LJ,fc,ft,f_x,f_y,scale(Ntot)
      double precision fthi,fthj,fth_c,Lx,Ly,ap,P,Pij,rijsq,dijsq_up
      double precision s(Ntot),dd,dr(Ntot,2),xa(Ntot,2),ya(Ntot,2),att
      double precision dk(Ntot,2),di_up(Ntot),di1j1,xhit,yhit,yhitout
      double precision PP(Ntot),c(Ntot),Vij,PT,PR
      integer countn(Ntot),nl(800,Ntot),N
      common /f3com/ alpha ! aspect ratio
      common /f4com/ exp,att
      common /f5com/ Lx,Ly,ap
      common /f6com/ P,PP,PT,PR
      common /f9com/ scale

      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
         PP(i)=0d0
      enddo
      P=0d0

      ! convert to from molecules to atoms
      do i=1,N
         c(i)=dcos(th(i)/scale(i))
         s(i)=dsin(th(i)/scale(i))
         dd=alpha(i)-1d0
         dr(i,1)=(1d0+dd)/(1d0+dd**2)*dd**2*D(i)*D1/2d0
         dr(i,2)=-(1d0+dd)/(1d0+dd**2)*D(i)*D1/2d0
         do k=1,2
            xa(i,k)=x(i)+dr(i,k)*c(i)
            ya(i,k)=y(i)+dr(i,k)*s(i)
         enddo
         dk(i,1)=D(i)*D1
         dk(i,2)=dd*D(i)*D1
         if(alpha(i).lt.2d0) then 
            di_up(i)=(dk(i,2)/2d0-dr(i,2))*2d0*D1
         else
            di_up(i)=(dk(i,1)/2d0+dr(i,1))*2d0*D1         
         endif
      enddo

      ! inter-particle interactions
      do i=1,N
         if(countn(i).ge.1) then
            do jj=1,countn(i)
               j=nl(jj,i)
               dij_up=(di_up(i)+di_up(j))/2d0  
               xij=x(i)-x(j)
               if(dabs(xij).lt.dij_up+att) then
                  yij=y(i)-y(j)                  
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij_up+att)**2) then
                     di1j1=(dk(i,1)+dk(j,1))/2d0
                     do ki=1,2
                        do kj=1,2
                           dij=(dk(i,ki)+dk(j,kj))/2d0
                           xij=xa(i,ki)-xa(j,kj)
                           yij=ya(i,ki)-ya(j,kj)
                           rijsq=xij**2+yij**2
                           if(rijsq.lt.(dij+att)**2) then
                              rij=dsqrt(rijsq)
                              if(exp .gt. 2.9) then
                                 LJ=(dij/rij)*(dij/rij)
                                 LJ=LJ*LJ*LJ
                                 fc=1d0/rij*LJ*(LJ-1d0)
                              else
                                 fc=(1d0-rij/dij)**(exp-1d0)/dij
                              endif    
                              fr=-fc/rij*dij**2/di1j1**2
                              f_x = fr*xij
                              f_y = fr*yij
                              fx(i)=fx(i)+f_x
                              fx(j)=fx(j)-f_x
                              fy(i)=fy(i)+f_y
                              fy(j)=fy(j)-f_y
                              fth(i)=fth(i)+dr(i,ki)*(c(i)*f_y-s(i)*f_x)
                              fth(j)=fth(j)-dr(j,kj)*(c(j)*f_y-s(j)*f_x)

                              Pij=-xij*f_x-yij*f_y
                              P=P+2d0*Pij
                              PP(i)=PP(i)+Pij
                              PP(j)=PP(j)+Pij
                           endif
                        enddo
                     enddo
                  end if
               end if
            enddo
         end if
      enddo

      ! boundary conditions
      PT=0d0
      PR=0d0
      do i=1,N
         do ki=1,2 
            xhit=(Lx-dk(i,ki))/2d0
            yhit=(Ly-dk(i,ki))/2d0
            if(ya(i,ki).gt.yhit) then ! top wall
               f_y=((ya(i,ki)-yhit)/dk(i,1))**(exp-1d0)/dk(i,1)
               fy(i)=fy(i)+f_y                  
               fth(i)=fth(i)+dr(i,ki)*c(i)*f_y

               PT=PT+f_y
               Pij=-2d0*(ya(i,ki)-yhit)*f_y
               P=P+Pij
               PP(i)=PP(i)+Pij
            elseif(ya(i,ki).lt.-yhit) then ! bottom wall
               f_y=((ya(i,ki)+yhit)/dk(i,1))**(exp-1d0)/dk(i,1)
               fy(i)=fy(i)+f_y
               fth(i)=fth(i)+dr(i,ki)*c(i)*f_y

               Pij=-2d0*(ya(i,ki)+yhit)*f_y
               P=P+Pij
               PP(i)=PP(i)+Pij
            endif
            if(xa(i,ki).gt.xhit) then ! right wall
               f_x=((xa(i,ki)-xhit)/dk(i,1))**(exp-1d0)/dk(i,1)
               fx(i)=fx(i)+f_x                  
               fth(i)=fth(i)-dr(i,ki)*s(i)*f_x

               PR=PR+f_x
               Pij=-2d0*(xa(i,ki)-xhit)*f_x
               P=P+Pij
               PP(i)=PP(i)+Pij
            elseif(xa(i,ki).lt.-xhit) then ! left wall, corners, outlet
               if(dabs(ya(i,ki))+min(xa(i,ki)+Lx/2d0,0d0).gt.AP/2d0)then
                  f_x=((xa(i,ki)+xhit)/dk(i,1))**(exp-1d0)/dk(i,1)
                  fx(i)=fx(i)+f_x                  
                  fth(i)=fth(i)-dr(i,ki)*s(i)*f_x

                  Pij=-2d0*(xa(i,ki)+xhit)*f_x
                  P=P+Pij
                  PP(i)=PP(i)+Pij
               elseif(xa(i,ki).gt.-Lx/2d0) then 
                  xij=xa(i,ki)+Lx/2d0
                  yij=ya(i,ki)-AP/2d0
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.dk(i,ki)**2/4d0) then
                     rij=dsqrt(rijsq)
                     fc=((dk(i,ki)/2d0-rij)/dk(i,1))**(exp-1d0)/dk(i,1)
                     f_x=-fc*xij/rij
                     f_y=-fc*yij/rij
                     fx(i)=fx(i)+f_x
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+dr(i,ki)*(c(i)*f_y-s(i)*f_x)   
 
                     Pij=-2d0*xij*f_x-2d0*yij*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij
                  endif
                  yij=ya(i,ki)+AP/2d0
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.dk(i,ki)**2/4d0) then
                     rij=dsqrt(rijsq)
                     fc=((dk(i,ki)/2d0-rij)/dk(i,1))**(exp-1d0)/dk(i,1)
                     f_x=-fc*xij/rij
                     f_y=-fc*yij/rij
                     fx(i)=fx(i)+f_x
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+dr(i,ki)*(c(i)*f_y-s(i)*f_x)    

                     Pij=-2d0*xij*f_x-2d0*yij*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij
                  endif
               else                 
                  yhitout=(AP-dk(i,ki))/2d0
                  if(ya(i,ki).gt.yhitout) then
                     f_y=((ya(i,ki)-yhitout)/dk(i,1))**(exp-1d0)/dk(i,1)
                     fy(i)=fy(i)+f_y                  
                     fth(i)=fth(i)+dr(i,ki)*c(i)*f_y

                     Pij=-2d0*(ya(i,ki)-yhitout)*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij
                  endif
                  if(ya(i,ki).lt.-yhitout) then
                     f_y=((ya(i,ki)+yhitout)/dk(i,1))**(exp-1d0)/dk(i,1)
                     fy(i)=fy(i)+f_y
                     fth(i)=fth(i)+dr(i,ki)*c(i)*f_y

                     Pij=-2d0*(ya(i,ki)+yhitout)*f_y
                     P=P+Pij
                     PP(i)=PP(i)+Pij
                  endif
               endif
            endif
         enddo
      enddo
      
      if(exp .gt. 2.9) then
         do i=1,N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
            fth(i)=fth(i)/6d0 
         enddo
         P=P/6d0
      endif
      
      do i=1,N
         fth(i)=fth(i)/scale(i)
      enddo

      PT=PT/Lx
      PR=PR/Ly
      P=P/4d0/Lx/Ly
      do i=1,N
         PP(i)=PP(i)*dble(N)/4d0/Lx/Ly
      enddo

      return							
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine func_disk(N,x,y,D,D1,V,countn,nl)
      parameter(Ntot = 8192)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),D(Ntot),D1,V,rij,xij,yij,dij,exp
      double precision LJ,Lx,Ly,ap,rijsq,att,Vij,xhit,yhit,yhitout,di
      integer countn(Ntot),nl(800,Ntot),N
      common /f4com/ exp,att
      common /f5com/ Lx,Ly,ap

      ! inter-particle interactions
      V=0d0
      do i=1,N
         if(countn(i).ge.1) then
            do jj=1,countn(i)
               j=nl(jj,i)               
               dij=D1*(d(i)+d(j))/2d0
               xij=x(i)-x(j)
               if(dabs(xij).lt.dij+att) then
                  yij=y(i)-y(j)                  
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij+att)**2) then
                     rij=dsqrt(rijsq)
                     if(exp .gt. 2.9) then
                        LJ=(dij/rij)*(dij/rij)
                        LJ=LJ*LJ*LJ
                        V=V+(LJ-1d0)*(LJ-1d0)
                     else
                        V=V+(1d0-rij/dij)**exp/exp-(att/dij)**exp/exp
                     endif    
                  end if
               end if
            enddo
         end if
      enddo

      ! boundary conditions
      do i=1,N
         di=d(i)*D1
         xhit=(Lx-d(i)*D1)/2d0
         yhit=(Ly-d(i)*D1)/2d0
         if(y(i).gt.yhit) then  ! top wall
            V=V+((y(i)-yhit)/di)**exp/exp
         elseif(y(i).lt.-yhit) then ! bottom wall
            V=V+((y(i)+yhit)/di)**exp/exp
         endif
         if(x(i).gt.xhit) then  ! right wall
            V=V+((x(i)-xhit)/di)**exp/exp
         elseif(x(i).lt.-xhit) then ! left wall, corner, outlet
            if(dabs(y(i))+min(x(i)+Lx/2d0,0d0).gt.AP/2d0)then
               V=V+((x(i)+xhit)/di)**exp/exp
            elseif(x(i).gt.-Lx/2d0) then   
               rijsq=(x(i)+Lx/2d0)**2+(y(i)-AP/2d0)**2
               if(rijsq.lt.di**2/4d0) then
                  rij=dsqrt(rijsq)
                  V=V+((di/2d0-rij)/di)**exp/exp
               endif
               rijsq=(x(i)+Lx/2d0)**2+(y(i)+AP/2d0)**2
               if(rijsq.lt.di**2/4d0) then
                  rij=dsqrt(rijsq)
                  V=V+((di/2d0-rij)/di)**exp/exp
               endif
            else                 
               yhitout=(AP-di)/2d0
               if(y(i).gt.yhitout) then
                  V=V+((y(i)-yhitout)/di)**exp/exp
               endif
               if(y(i).lt.-yhitout) then
                  V=V+((y(i)+yhitout)/di)**exp/exp
               endif
            endif
         endif
      enddo

      if(exp.gt.2.9) then
         V=V/72d0
      endif

      return				
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine dfunc_disk(N,x,y,D,D1,fx,fy,fth,countn,nl)
      parameter(Ntot = 8192)
      double precision pi
      parameter(pi=3.1415926535897932d0)
      double precision x(Ntot),y(Ntot),D(Ntot),D1,dij,fx(Ntot),fy(Ntot)
      double precision fth(Ntot),rij,xij,yij,fr,exp,LJ,fc,ft,f_x,f_y
      double precision Lx,Ly,ap,P,rijsq,PP(Ntot),Pij,Vij,att
      double precision xhit,yhit,yhitout,di,PT,PR
      integer countn(Ntot),nl(800,Ntot),N
      common /f4com/ exp,att
      common /f5com/ Lx,Ly,ap
      common /f6com/ P,PP,PT,PR

      do i=1,N
         fx(i)=0d0
         fy(i)=0d0
         fth(i)=0d0
         PP(i)=0d0
      enddo
      P=0d0

      ! inter-particle interactions
      do i=1,N
         if(countn(i).ge.1) then
            do jj=1,countn(i)
               j=nl(jj,i)
               dij=D1*(d(i)+d(j))/2d0
               xij=x(i)-x(j)
               if(dabs(xij).lt.dij+att) then
                  yij=y(i)-y(j)                  
                  rijsq=xij**2+yij**2
                  if(rijsq.lt.(dij+att)**2) then
                     rij=dsqrt(rijsq)
                     if(exp .gt. 2.9) then
                        LJ=(dij/rij)*(dij/rij)
                        LJ=LJ*LJ*LJ
                        fc=1d0/rij*LJ*(LJ-1d0)
                     else
                        fc=(1d0-rij/dij)**(exp-1d0)/dij
                     endif  
                     fr=-fc/rij
                     f_x = fr*xij
                     f_y = fr*yij
                     fx(i)=fx(i)+f_x
                     fx(j)=fx(j)-f_x
                     fy(i)=fy(i)+f_y
                     fy(j)=fy(j)-f_y

                     Pij=-xij*f_x-yij*f_y
                     P=P+2d0*Pij
                     PP(i)=PP(i)+Pij
                     PP(j)=PP(j)+Pij
                  endif
               end if
            enddo
         end if
      enddo

      ! boundary conditions
      PT=0d0
      PR=0d0
      do i=1,N
         di=D(i)*D1
         xhit=(Lx-di)/2d0
         yhit=(Ly-di)/2d0
         if(y(i).gt.yhit) then  ! top wall
            f_y=((y(i)-yhit)/di)**(exp-1d0)/di
            fy(i)=fy(i)+f_y       
           
            PT=PT+f_y
            Pij=-2d0*(y(i)-yhit)*f_y
            P=P+Pij
            PP(i)=PP(i)+Pij
         elseif(y(i).lt.-yhit) then ! bottom wall
            f_y=((y(i)+yhit)/di)**(exp-1d0)/di
            fy(i)=fy(i)+f_y

            Pij=-2d0*(y(i)+yhit)*f_y
            P=P+Pij
            PP(i)=PP(i)+Pij
         endif
         if(x(i).gt.xhit) then  ! right wall
            f_x=((x(i)-xhit)/di)**(exp-1d0)/di
            fx(i)=fx(i)+f_x        
          
            PR=PR+f_x
            Pij=-2d0*(x(i)-xhit)*f_x
            P=P+Pij
            PP(i)=PP(i)+Pij
         elseif(x(i).lt.-xhit) then ! left wall, corners, outlet
            if(dabs(y(i))+min(x(i)+Lx/2d0,0d0).gt.AP/2d0)then
               f_x=((x(i)+xhit)/di)**(exp-1d0)/di
               fx(i)=fx(i)+f_x   
               
               Pij=-2d0*(x(i)+xhit)*f_x
               P=P+Pij
               PP(i)=PP(i)+Pij
            elseif(x(i).gt.-Lx/2d0) then 
               xij=x(i)+Lx/2d0
               yij=y(i)-AP/2d0
               rijsq=xij**2+yij**2
               if(rijsq.lt.di**2/4d0) then
                  rij=dsqrt(rijsq)
                  fc=((di/2d0-rij)/di)**(exp-1d0)/di
                  f_x=-fc*xij/rij
                  f_y=-fc*yij/rij
                  fx(i)=fx(i)+f_x
                  fy(i)=fy(i)+f_y
                  Pij=-2d0*xij*f_x-2d0*yij*f_y
                  P=P+Pij
                  PP(i)=PP(i)+Pij
               endif
               yij=y(i)+AP/2d0
               rijsq=xij**2+yij**2
               if(rijsq.lt.di**2/4d0) then
                  rij=dsqrt(rijsq)
                  fc=((di/2d0-rij)/di)**(exp-1d0)/di
                  f_x=-fc*xij/rij
                  f_y=-fc*yij/rij
                  fx(i)=fx(i)+f_x
                  fy(i)=fy(i)+f_y

                  Pij=-2d0*xij*f_x-2d0*yij*f_y
                  P=P+Pij
                  PP(i)=PP(i)+Pij
               endif
            else                 
               yhitout=(AP-di)/2d0
               if(y(i).gt.yhitout) then
                  f_y=((y(i)-yhitout)/di)**(exp-1d0)/di
                  fy(i)=fy(i)+f_y    
              
                  Pij=-2d0*(y(i)-yhitout)*f_y
                  P=P+Pij
                  PP(i)=PP(i)+Pij
               endif
               if(y(i).lt.-yhitout) then
                  f_y=((y(i)+yhitout)/di)**(exp-1d0)/di
                  fy(i)=fy(i)+f_y

                  Pij=-2d0*(y(i)+yhitout)*f_y
                  P=P+Pij
                  PP(i)=PP(i)+Pij
               endif
            endif
         endif
      enddo
      
      if(exp .gt. 2.9) then
         do i=1,N
            fx(i)=fx(i)/6d0
            fy(i)=fy(i)/6d0
         enddo
         P=P/6d0
      endif
      
      PT=PT/Lx
      PR=PR/Ly
      P=P/4d0/Lx/Ly
      do i=1,N
         PP(i)=PP(i)*dble(N)/4d0/Lx/Ly
      enddo

      return							
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE frprmn(N,x,y,th,D,D1,ftol,ftol1,iter,fret)
      parameter(Ntot = 8192)
      INTEGER its,iter,ITMAX
      double precision fret,ftol,EPS,ftol1
      PARAMETER (EPS=1d-10,ITMAX=1000000000)
      double precision dgg,fp,gam,gg,gx(Ntot),gy(Ntot),hx(Ntot),hy(Ntot)
      double precision D(Ntot),D1,xix(Ntot),xiy(Ntot),xith(Ntot),width
      double precision x(Ntot),y(Ntot),maxdis,xp(Ntot),yp(Ntot)
      double precision th(Ntot),hth(Ntot),gth(Ntot),exp,att
      integer N,countn(Ntot),nl(800,Ntot)

      ! not needed
      double precision f1,f2,f3,fxe,fye,fthe
      double precision xi,yi,thi,max1,max2,max3,del
      double precision alpha(Ntot),Lx,Ly,ap,alpha0

      ! not needed
      common /f3com/ alpha ! aspect ratio
      common /f5com/ Lx,Ly,ap
      common /f8com/ alpha0

      common /f2com/ width      
      common /f4com/ exp,att

      iter=0

      call makelist(N,x,y,D,D1,xp,yp,countn,nl)
      call func(N,x,y,th,D,D1,fp,countn,nl)

      if (fp.lt.ftol*dble(N).and.att.eq.0d0) then
         fret=fp 
         return
      endif

      call dfunc(N,x,y,th,D,D1,xix,xiy,xith,countn,nl)

      do i=1,N
        gx(i)=-xix(i)
	gy(i)=-xiy(i)
        gth(i)=-xith(i)
        hx(i)=gx(i)
	hy(i)=gy(i)
        hth(i)=gth(i)
        xix(i)=hx(i)
	xiy(i)=hy(i)
        xith(i)=hth(i)
      enddo

      do its=1,ITMAX
         iter=its
         
         call linmin(N,x,y,th,D,D1,xix,xiy,xith,fret,xp,yp,countn,nl)

c         write(*,*) its, fret/dble(N)

         if(att.eq.0d0) then
            if(dabs(fret-fp).lt.ftol1*fp.or.fret.lt.ftol*dble(N)) then
               return
            endif
         else
            if(dabs(fret-fp).lt.ftol1) then
               return
            end if 
         endif
         
         call CG_check(N,x,y,xp,yp,maxdis)	     
         if(maxdis.gt.width*D1) then
            call makelist(N,x,y,D,D1,xp,yp,countn,nl)
         end if
         
         fp=fret
         
         call dfunc(N,x,y,th,D,D1,xix,xiy,xith,countn,nl)
         
         gg=0d0
         dgg=0d0
         
         do i=1,N
            gg=gg+gx(i)*gx(i)+gy(i)*gy(i)+gth(i)*gth(i)
            dgg=dgg+(xix(i)+gx(i))*xix(i)+(xiy(i)+gy(i))*xiy(i)
     +           +(xith(i)+gth(i))*xith(i)
         enddo
         
         if(gg.eq.0d0) then
            return
         endif
         gam=dgg/gg
         do i=1,N
            gx(i)=-xix(i)
            gy(i)=-xiy(i)
            gth(i)=-xith(i)
            hx(i)=gx(i)+gam*hx(i)
            hy(i)=gy(i)+gam*hy(i)
            hth(i)=gth(i)+gam*hth(i)
            xix(i)=hx(i)
            xiy(i)=hy(i)
            xith(i)=hth(i)
         enddo
      enddo
      
c     pause 'frprmn maximum iterations exceeded'
      
      return
      END
C (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE linmin(N,x,y,th,D,D1,xix,xiy,xith,fret,
     +     xpr,ypr,countnr,nlr)
      double precision fret,TOL
      PARAMETER (Ntot=8192,TOL=1d-8)
C     USES dbrent,df1dim,mnbrak
      double precision ax,bx,fa,fb,fx,xmin,xx,xp(Ntot),yp(Ntot),dbrent
      double precision pxcom(Ntot),pycom(Ntot),xixcom(Ntot),xiycom(Ntot)
      double precision x(Ntot),y(Ntot),th(Ntot),xix(Ntot),xiy(Ntot)
      double precision xith(Ntot),Dcom(Ntot),D1com,D(Ntot),D1,width
      double precision pthcom(Ntot),xithcom(Ntot),f1dim,df1dim
      double precision xpr(Ntot),ypr(Ntot)
      integer countn(Ntot),nl(800,Ntot),countnr(Ntot),nlr(800,Ntot)
      integer N,ncom

      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      common /f2com/ width

      EXTERNAL df1dim
      EXTERNAL f1dim

      do i=1,N
        pxcom(i)=x(i)
        pycom(i)=y(i)
        pthcom(i)=th(i)
        xixcom(i)=xix(i)
	xiycom(i)=xiy(i)
        xithcom(i)=xith(i)
        Dcom(i)=D(i)
      enddo
      D1com=D1
      ncom=N

      do i=1,N
         xp(i)=xpr(i)
         yp(i)=ypr(i)
         countn(i)=countnr(i)
         do j=1,countn(i)
            nl(j,i)=nlr(j,i)
         enddo
      enddo

      ax=0d0
      xx=1d0
      call mnbrak(ax,xx,bx,fa,fx,fb,f1dim)
      fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,xmin)
      do i=1,N
        xix(i)=xmin*xix(i)
        xiy(i)=xmin*xiy(i)
        xith(i)=xmin*xith(i)
        x(i)=x(i)+xix(i)
	y(i)=y(i)+xiy(i)
        th(i)=th(i)+xith(i)
      enddo

      return
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE mnbrak(ax,bx,cx,fa,fb,fc,func)
      double precision ax,bx,cx,fa,fb,fc,func,GOLD,GLIMIT,TINY
c      EXTERNAL func
      PARAMETER (GOLD=1.618034d0, GLIMIT=100., TINY=1.d-20)
      double precision dum,fu,q,r,u,ulim
      fa=func(ax)
      fb=func(bx)
      if(fb.gt.fa)then ! was gt
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
      fc=func(cx)
1     if(fb.gt.fc)then ! was ge
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2d0*dsign(max(dabs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0d0)then
          fu=func(u)
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        else if((cx-u)*(u-ulim).gt.0d0)then
          fu=func(u)
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
            fu=func(u)
          endif
        else if((u-ulim)*(ulim-cx).ge.0d0)then ! was ge
          u=ulim
          fu=func(u)
        else
          u=cx+GOLD*(cx-bx)
          fu=func(u)
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION df1dim(x)
      PARAMETER (Ntot=8192)
      double precision df1dim,x,maxdis,fx(Ntot),fy(Ntot),fth(Ntot)
      double precision pxcom(Ntot),pycom(Ntot),pthcom(Ntot)
      double precision xixcom(Ntot),xiycom(Ntot),xithcom(Ntot)
      double precision xt(Ntot),yt(Ntot),tht(Ntot),xp(Ntot),yp(Ntot)
      double precision Dcom(Ntot),D1com,alpha(Ntot),width
      integer ncom,countn(Ntot),nl(800,Ntot)
      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio

      do i=1,ncom
        xt(i)=pxcom(i)+x*xixcom(i)
        yt(i)=pycom(i)+x*xiycom(i)
        tht(i)=pthcom(i)+x*xithcom(i)
      enddo

      call CG_check(ncom,xt,yt,xp,yp,maxdis)
      if(maxdis.gt.width*D1com) then
	call makelist(ncom,xt,yt,Dcom,D1com,xp,yp,countn,nl)
      end if	
      call dfunc(ncom,xt,yt,tht,Dcom,D1com,fx,fy,fth,countn,nl)

      df1dim=0d0
      do i=1,ncom
        df1dim=df1dim+fx(i)*xixcom(i)+fy(i)*xiycom(i)+fth(i)*xithcom(i)
      enddo

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION dbrent(ax,bx,cx,f,df,tol,xmin)
      INTEGER ITMAX
      double precision dbrent,ax,bx,cx,tol,xmin,df,f,ZEPS
      EXTERNAL df,f
      PARAMETER (ITMAX=10000,ZEPS=1.0e-12)
      INTEGER iter
      double precision a,b,d,d1,d2,du,dv,dw,dx,e,fu,fv,fw
      double precision fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm
      LOGICAL ok1,ok2
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
      fx=f(x)
      fv=fx
      fw=fx
      dx=df(x)
      dv=dx
      dw=dx
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*dabs(x)+ZEPS
        tol2=2.*tol1
        if(dabs(x-xm).le.(tol2-0.5d0*(b-a))) goto 3
        if(dabs(e).gt.tol1) then
          d1=2.*(b-a)
          d2=d1
          if(dw.ne.dx) d1=(w-x)*dx/(dx-dw)
          if(dv.ne.dx) d2=(v-x)*dx/(dx-dv)
          u1=x+d1
          u2=x+d2
          ok1=((a-u1)*(u1-b).gt.0d0).and.(dx*d1.le.0d0)
          ok2=((a-u2)*(u2-b).gt.0d0).and.(dx*d2.le.0d0)
          olde=e
          e=d
          if(.not.(ok1.or.ok2))then
            goto 1
          else if (ok1.and.ok2)then
            if(dabs(d1).lt.dabs(d2))then
              d=d1
            else
              d=d2
            endif
          else if (ok1)then
            d=d1
          else
            d=d2
          endif
          if(dabs(d).gt.dabs(0.5d0*olde))goto 1
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=dsign(tol1,xm-x)
          goto 2
        endif
1       if(dx.ge.0d0) then
          e=a-x
        else
          e=b-x
        endif
        d=0.5d0*e
2        if(dabs(d).ge.tol1) then
          u=x+d
           fu=f(u)
        else
          u=x+dsign(tol1,d)
          fu=f(u)
          if(fu.gt.fx)goto 3
        endif
        du=df(u)
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          dv=dw
          w=x
          fw=fx
          dw=dx
          x=u
          fx=fu
          dx=du
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            dv=dw
            w=u
            fw=fu
            dw=du
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
            dv=du
          endif
        endif
11    continue
c      pause 'dbrent exceeded maximum iterations'
3     xmin=x
      dbrent=fx
      return
      END
c  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION f1dim(x)
      PARAMETER (Ntot=8192)
      double precision f1dim,x
      double precision pxcom(Ntot),pycom(Ntot),pthcom(Ntot)
      double precision xixcom(Ntot),xiycom(Ntot),xithcom(Ntot)
      double precision xt(Ntot),yt(Ntot),tht(Ntot),width
      double precision Dcom(Ntot),maxdis,D1com
      double precision xp(Ntot),yp(Ntot),alpha(Ntot)
      integer nl(800,Ntot),countn(Ntot),ncom
      COMMON /f1com/ pxcom,pycom,pthcom,xixcom,xiycom,xithcom,
     +     Dcom,D1com,xp,yp,thp,countn,nl,ncom
      
      common /f2com/ width
      common /f3com/ alpha ! aspect ratio

      do i=1,ncom
        xt(i)=pxcom(i)+x*xixcom(i)
        yt(i)=pycom(i)+x*xiycom(i)
        tht(i)=pthcom(i)+x*xithcom(i)
      enddo

      call CG_check(ncom,xt,yt,xp,yp,maxdis)
      if(maxdis.gt.width*D1com) then
	call makelist(ncom,xt,yt,Dcom,D1com,xp,yp,countn,nl)
      end if
      call func(ncom,xt,yt,tht,Dcom,D1com,f1dim,countn,nl)

      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software Dt+;39.


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      double precision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END
               

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      double precision function sigma(rij,xij,yij,thi,thj,
     +     di,dj,alphai,alphaj)
      
      implicit none

      double precision xij,yij,thi,thj,di,dj
      double precision rmui,rmuj,muij
      double precision rij, li, lj, alphai, alphaj
      double precision chi, alp, sig_0
      double precision n1, n2, d
      double precision dlnsig

      if(alphai.eq.1d0.and.alphaj.eq.1d0) then
         sigma = 0.5d0*(di+dj)
      else
         li = alphai*di
         lj = alphaj*dj
         
         rmui=(xij*dcos(thi)+yij*dsin(thi))/rij
         rmuj=(xij*dcos(thj)+yij*dsin(thj))/rij
         muij=dcos(thi-thj)

         chi = dsqrt((li*li-di*di)*(lj*lj-dj*dj)/
     +      ((lj*lj+di*di)*(li*li+dj*dj)))
         alp = dsqrt(dsqrt((li*li-di*di)*(lj*lj+di*di)/
     +        ((lj*lj-dj*dj)*(li*li+dj*dj))))
         sig_0 = 0.5d0*dsqrt((di*di+dj*dj)*2d0)

         n1=alp*rmui+rmuj/alp
         n2=alp*rmui-rmuj/alp
         d=chi*muij

         sigma=sig_0/dsqrt(1d0-0.5d0*chi*(n1*n1/(1d0+d)+n2*n2/(1d0-d)))
      end if

      return
      end
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma2(rij,xij,yij,thi,thj,di,dj,
     +     sig,dlnsig,fthi,fthj,alphai,alphaj)
      
      double precision xij,yij,thi,thj,di,dj
      double precision rmui,rmuj,muij
      double precision rij, li, lj, alphai, alphaj
      double precision chi, alp, sig_0
      double precision d,n1,n2,n3,n4,n5,n6,n7,n8
      double precision sig,dlnsig,pre
      double precision s2,rmui2,rmuj2,muij2,fthi,fthj

      if(alphai.eq.1d0.and.alphaj.eq.1d0) then
         sig = 0.5d0*(di+dj)
         dlnsig = 0d0
         fthi = 0d0
         fthj = 0d0
      else
         li = alphai*di
         lj = alphaj*dj
         
         rmui=(xij*dcos(thi)+yij*dsin(thi))/rij
         rmuj=(xij*dcos(thj)+yij*dsin(thj))/rij
         muij=dcos(thi-thj)
 
         chi = dsqrt((li*li-di*di)*(lj*lj-dj*dj)/
     +      ((lj*lj+di*di)*(li*li+dj*dj)))
         alp = dsqrt(dsqrt((li*li-di*di)*(lj*lj+di*di)/
     +        ((lj*lj-dj*dj)*(li*li+dj*dj))))
         sig_0 = 0.5d0*dsqrt((di*di+dj*dj)*2d0)
         
         d=chi*muij
         n1=alp*rmui+rmuj/alp
         n2=n1/(1d0+d)
         n3=alp*rmui-rmuj/alp
         n4=n3/(1d0-d)

         sig = sig_0/dsqrt(1d0-0.5d0*chi*(n1*n2+n3*n4))

         s2 = sig/sig_0
         s2 = s2*s2
         rmui2 = (yij*dcos(thi)-xij*dsin(thi))/rij
         rmuj2 = (yij*dcos(thj)-xij*dsin(thj))/rij
         n5 = alp*rmui2+rmuj2/alp
         n6 = alp*rmui2-rmuj2/alp
         pre=0.5d0*chi*s2
         dlnsig = -pre*(n2*n5+n4*n6)

         n7 = n2+n4
         muij2 = dsin(thi-thj)
         n8 = chi*muij2*(n2*n2-n4*n4)
         fthi = 0.5d0*pre*(2d0*alp*rmui2*(n2+n4)+n8)
         fthj = 0.5d0*pre*(2d0/alp*rmuj2*(n2-n4)-n8)                
      end if

      return
      end

