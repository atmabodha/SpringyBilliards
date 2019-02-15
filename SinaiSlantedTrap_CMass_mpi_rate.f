      program main

         include 'mpif.h'
         parameter(send_data_tag = 2000,return_data_tag1 = 2001)
         parameter(return_data_tag2 = 2002,return_data_tag3 = 2003)

         integer my_id, root_process, ierr, status(MPI_STATUS_SIZE)
         integer num_procs, an_id,mpi_cnt,imnt

         integer MAXiter,NI,NIa,NIas
         parameter(MAXiter=500000,NI=10000)
         integer iter,iI,iT,iM,iTmin,iTmins
         integer iTcount(NI),total_iTcount,iTT
         integer bhit,chbrak,barcoll
         integer idum,rtaccept

         double precision start,finish,run_time

c        Parameters to read mass values from rate_m.txt
         parameter(Nm=22)
         integer i,j,k
         character aa
         character*3 bb
         character*4 cc
         integer file_index(Nm),fp
         double precision mass(Nm)

c        Parameters required for Least Squares
         integer arrsize
         double precision arr1(MAXiter),arr2(MAXiter)
         double precision xbar,ybar,xy,xsq
         double precision slope(Nm),intercept


c        Disc is of radius r and centered at (0,L)
         double precision r,L,angle,A0,A,w,phase,phase0,dia

         double precision mp,mw,vw,vwb,mparray(Nm)

         double precision x0,z0,x0s(NI),z0s(NI)
         double precision u0,w0,E0,Eb0,Ep0
         double precision vel,maxvel
         double precision fac

         double precision x1,z1
         double precision x2,z2

         double precision u1,w1
         double precision u2,w2

c       Coordinates of the Trapezium 
c       (-L,0), (L,0), (x4,z4) and (x6,z6)
         double precision m4,m5,m6,c4,c5,c6
         double precision x4,z4,x6,z6
         double precision tanglet,tangleb,h1,h2,hm,trapheight

         double precision ttemp,ttemp1,ttemp2,t1,t2,tHmax
         double precision ts,m,c
         double precision xtemp,ztemp
         double precision tanglei,tangler,tanglec
         double precision tstep,en(MAXiter),entemp
         double precision enw(MAXiter),enwavg(MAXiter),enwlog(MAXiter)
         double precision tA,tB,tC

         double precision total_time,total_time_temp
         double precision broot1,broot2,xacc


         double precision PI         

         double precision bvel,funcheight
         external bvel,funcheight
         double precision rtsafe,ran0
         external rtsafe,funcd,zbrak,ran0,brak_osc


         integer filex
         character*35 fileb
         character*14 filed
         character*70 filee(Nm)
         common /heightparam/ r,L,A,w,PI,phase
         common /funcdparam/ z1,w1,total_time

         common /brakparam/ chbrak
         common /rtparam/ rtaccept

         PI=1.0d0
         PI=4.0d0*datan(PI)

ccc      SINAI BILLIARD PARAMETERS
c        Centered at (0,0)
         
         tangleb=06.d0*PI/180.d0
         tanglet=06.d0*PI/180.d0
         r=20.d0/10.d0
         L=4.0d0         
 
         fileb="data/SinaiSlantedTrap_06-06-20_Mass"
         filed="e8_Eb9_10k.txt"

         open(50,file="data/SinaiSlantedTrap_06-06-20_Eb9_rate.txt"
     *                      ,status="unknown") 

         m4=-1.d0/dtan(tangleb)
         c4=L*m4
         m5=dtan(tanglet)
         c5=2.d0*L+L*(1.d0+2.d0*dtan(tangleb))*dtan(tanglet)
         m6=-m4
         c6=c4

         x4=-L*(1.d0+2.d0*dtan(tangleb))
         z4=2.d0*L
         x6=(c6-c5)/(m5-m6)
c         x6=L*(1.d0+(2.d0+(1.d0+2.d0*dtan(tangleb))
c     *               *dtan(tanglet))*dtan(tangleb))
c     *                   /(1-dtan(tanglet)*dtan(tangleb))
         z6=x6*m5+c5

         h1=z4
         h2=z6
         hm=c5

         dia=dsqrt(z6**2+(x6+L)**2)

c         A0=0.2
         w=2.0d0*PI

         
ccc      INITIAL CONDITIONS

         x0=0.3d0
         z0=0.33d0

         E0=1.0d0
         Eb0=0.9d0
         Ep0=E0-Eb0
         mw=1.0d0
         w0=u0
         A0=dsqrt(2.0d0*Eb0/mw)/w

         tstep=0.1*(2.0*PI/w)

         open(16,file="SinaiSlantedTrap_gam-m-slope.txt"
     *                    ,access="append",status="old")
         open(18,file="data/rate_m.txt",status="unknown")

         do i=1,Nm
            aa=""
            bb=""
            read(18,*) j
            
            if(j.lt.1000) then
               write(aa,10) 0
               write(bb,30) j
               cc=aa//bb
            else
               write(cc,40) j
             endif
            write(*,*) i,j,cc

            mass(i)=j*1.0d-8
            file_index(i)=19+i

            filee(i)=fileb//cc//filed
            write(*,*) 'File Name=',filee(i)

            open(file_index(i),file=filee(i),status="unknown")


         enddo
10       format(I1)
40       format(I4)
30       format(I3)


c MPI BEGIN

         root_process = 0

         call MPI_INIT(ierr)

         call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)
         call MPI_COMM_SIZE(MPI_COMM_WORLD,num_procs,ierr)

c
c ROOT MPI PROCESS
c

         if(my_id.eq.root_process) then

            write(*,*) 'MPI',my_id
            write(*,*) "A0=",A0
            write(*,*) 'Eb/E=',0.5d0*mw*(A0*w)**2/E0,E0
            write(*,*) 'm/M=',mp/mw
            write(*,*) 'u0=',u0
 
           iTmin=MAXiter
            total_iTcount=0

            if(num_procs.lt.(Nm/2+1)) pause 'Not enough processors'
            num_procs=Nm/2+1


            do iI=1,num_procs-1
                  call MPI_RECV(slope(iI*2-1),1,MPI_DOUBLE_PRECISION,iI,
     *             return_data_tag1,MPI_COMM_WORLD,status,ierr)


                  call MPI_RECV(slope(iI*2),1,MPI_DOUBLE_PRECISION,iI,
     *             return_data_tag1,MPI_COMM_WORLD,status,ierr)
            enddo

            do iM=1,Nm
               write(*,*) mass(iM),abs(slope(iM))

               write(50,*) mass(iM),abs(slope(iM))

                  arr1(iM)=log(mass(iM))
                  arr2(iM)=log(abs(slope(iM)))

                  xbar=xbar+arr1(iM)
                  ybar=ybar+arr2(iM)
                  xy=xy+arr1(iM)*arr2(iM)
                  xsq=xsq+(arr1(iM))**2

            enddo

            xbar=xbar/Nm
            ybar=ybar/Nm

            slope(1)=(xy-Nm*xbar*ybar)/(xsq-Nm*xbar**2)
            intercept=(ybar*xsq-xbar*xy)/(xsq-Nm*xbar**2)

            write(*,*) 'tangleb','tanglet','r/L','slope'
            write(*,*) tangleb*180.d0/PI,tanglet*180.d0/PI,r/L,slope(1)

            write(16,*) tangleb*180.d0/PI,tanglet*180.d0/PI,r/L
     *                  ,slope(1),intercept,Eb0
         else

c           open(19,file="SinaiTrap.txt",status="unknown")

            start=MPI_Wtime()
            idum=mod(int(start)+my_id,10000)
c            write(*,*) 'start,idum=',my_id,start,idum

c            write(*,*) my_id,mass(my_id*2-1),mass(my_id*2)


            do iM=1,2
               mp=mass(my_id*2-2+iM)
               u0=dsqrt(Ep0/mp)
               fp=file_index(my_id*2-2+iM)
 
               write(*,*) 'Mass',my_id,iM,mp,u0

            do iter=1,MAXiter
               enwavg(iter)=0.0
            enddo
            NIa=0
            iTmin=MAXiter

            do imnt=1,NI

                x0=r+(L-r)*ran0(idum)
                z0=L+r+(L-r)*ran0(idum)

c                 x0=2.0
c                 z0=4.0

               do iter=1,MAXiter
                  enw(iter)=0.0
               enddo
            
               A=A0
               phase=0.0


            x2=x0
            z2=z0

            u2=u0
            w2=w0

            vel=dsqrt(u2**2+w2**2)

            bhit=0
            total_time=0.0
            maxvel=0.0
            iT=0

            chbrak=0
            iter=0
            rtaccept=1

            do while((iter.le.MAXiter))

              iter=iter+1
               x1=x2
               z1=z2

               u1=u2
               w1=w2

c               if((my_id.eq.1).and.(imnt.eq.1)) 
c     *            write(*,*) my_id,imnt,iter,x1,z1,u1,w1,bhit
 
           
               ts=-x1/u1
               m=w1/u1
               c=z1+w1*ts-L

                barcoll=0
              
               if(u1.gt.0.d0) then
                  x2=(-m*c-dsqrt(r**2*(1+m**2)-c**2))/(1+m**2)
               else
                  x2=(-m*c+dsqrt(r**2*(1+m**2)-c**2))/(1+m**2)
               endif
               t2=(x2-x1)/u1 

               if(t2.gt.0.d0) then
                  bhit=3
               else
                  m=w1/u1
                  c=z1-x1*w1/u1

                  xtemp=(c5-c)/(m-m5)
                  ttemp=(xtemp-x1)/u1
                  ztemp=z1+ttemp*w1

                  if((ttemp.gt.0.d0).and.(bhit.ne.5).and.
     *                  (xtemp.ge.x4).and.(xtemp.le.x6)) then
                     t2=(x1*dtan(tanglet)+hm-z1)/(w1-u1*dtan(tanglet))
                     bhit=5
                  else
                     phase0=w*total_time+phase

                     tHmax=dia/dsqrt(w1**2+u1**2)       

                     call brak_osc(z1,w1,0.d0,A,w,phase0
     *                   ,barcoll,bhit,ttemp1,ttemp2,tHmax)

                      if(barcoll.eq.1) then
                        ttemp=0.d0                   
                        xacc=1.0d-12
                        ttemp=rtsafe(funcd,ttemp1,ttemp2,xacc)

                        xtemp=x1+ttemp*u1
                        ztemp=z1+ttemp*w1

                        if(ztemp.lt.(xtemp*m6+c6)) then
                           bhit=6
                           x2=(c6-c)/(m-m6)
                           t2=(x2-x1)/u1
                        elseif(ztemp.lt.(xtemp*m4+c4)) then
                           bhit=41
                           x2=(c4-c)/(m-m4)
                           t2=(x2-x1)/u1
                         else
                           bhit=10
                           t2=ttemp
                        endif
                      else
                         
                         x2=(c6-c)/(m-m6)
                         t2=(x2-x1)/u1

                         if((t2.gt.0.d0).and.(bhit.ne.6)) then
                            bhit=6
                         else
                            x2=(c4-c)/(m-m4)
                            t2=(x2-x1)/u1
                            bhit=42
                         endif
                      endif
                  endif
                endif

               x2=x1+t2*u1
               z2=z1+t2*w1

               total_time=total_time+t2

               if(bhit.eq.10) then
                  u2=u1
                   vw=bvel(total_time)
                   vwb=vw

                phase0=mod(datan2(w*A*dsin(w*total_time+phase),vw)
     *                             -w*(total_time),2*PI)
                    
                    w2=(2.0d0*mw*vw-(mw-mp)*w1)/(mp+mw)
                    vw=(2.0d0*mp*w1+(mw-mp)*vw)/(mp+mw)

                    A=dsqrt(z2**2+vw**2/w**2)
              phase=dmod(datan2(w*z2,vw)-w*(total_time),2*PI)

                
               elseif((bhit.eq.41).or.(bhit.eq.42)) then
                  tanglei=datan2(w1,u1)
                  tanglec=PI+2.0*tangleb-tanglei
                  
                  vel=dsqrt(u1**2+w1**2)
                  u2=dcos(tanglec)*vel         
                  w2=dsin(tanglec)*vel  
 
               elseif(bhit.eq.6) then
                  tanglei=datan2(w1,u1)
                  tanglec=PI-2.0*tangleb-tanglei
                  
                  vel=dsqrt(u1**2+w1**2)
                  u2=dcos(tanglec)*vel         
                  w2=dsin(tanglec)*vel  
 
               elseif(bhit.eq.5) then
                  tanglei=datan2(w1,u1)
                  tanglec=2.0*tanglet-tanglei
                  
                  vel=dsqrt(u1**2+w1**2)
                  u2=dcos(tanglec)*vel         
                  w2=dsin(tanglec)*vel  
               elseif(bhit.eq.3) then
                  tanglec=datan2(z2-L,x2)
                  tanglei=datan2(w1,u1)

                  tangler=PI+2.0d0*tanglec-tanglei

                  vel=dsqrt(u1**2+w1**2)
                  u2=dcos(tangler)*vel     
                  w2=dsin(tangler)*vel         
               endif
 
               if(total_time.gt.(iT*tstep)) then
                 iTnew=int(total_time/tstep)+1
                 if(iTnew.le.MAXiter) then
                    do iIT=iT+1,iTnew
                        enw(iIT)=A**2
                     enddo
                  endif
                  iT=iTnew

               endif
c               write(19,*) iter,x2,z2,u2,w2,bhit,total_time
            enddo

            if((rtaccept.eq.1).and.(iT.gt.1000)) then
               if(iTmin.gt.iT) iTmin=iT
               
               do iT=1,iTmin
                  enwavg(iT)=enwavg(iT)+enw(iT)
               enddo
               NIa=NIa+1
            endif

            if(iTmin.gt.MAXiter) iTmin=MAXiter

           enddo

           arrsize=0
           xbar=0.d0
           ybar=0.d0
           xy=0.d0
           xsq=0.d0

            do iter=1,iTmin
         enwlog(iter)=log(abs((enwavg(iter)/NIa)*0.5d0*mw*w**2-0.5))
               write(fp,*) iter,iter*tstep*(0.5d0*w/PI)
     *               ,enwlog(iter),mp,(enwavg(iter)/NIa)*0.5d0*mw*w**2

               if(enwlog(iter).ge.-1.9) then
                  arrsize=arrsize+1
                  arr1(arrsize)=iter*tstep*(0.5d0*w/PI)
                  arr2(arrsize)=enwlog(iter)

                  xbar=xbar+arr1(arrsize)
                  ybar=ybar+arr2(arrsize)
                  xy=xy+arr1(arrsize)*arr2(arrsize)
                  xsq=xsq+(arr1(arrsize))**2
               endif

            enddo

            if(enwlog(iTmin).gt.-1.9) pause 'log(E-0.5)>-1.9'

            xbar=xbar/arrsize
            ybar=ybar/arrsize

            slope(iM)=(xy-arrsize*xbar*ybar)/(xsq-arrsize*xbar**2)
            intercept=(ybar*xsq-xbar*xy)/(xsq-arrsize*xbar**2)

           finish=MPI_Wtime()
           write(*,*) 'finish',my_id,finish-start,mp,enwlog(iTmin)
           write(*,60) mp,slope(iM),intercept


           enddo

60         format(f12.9,' ',f12.9,' ',f12.9)

          call MPI_SEND(slope(1),1,MPI_DOUBLE_PRECISION,root_process
     *             ,return_data_tag1,MPI_COMM_WORLD,ierr)

          call MPI_SEND(slope(2),1,MPI_DOUBLE_PRECISION,root_process
     *             ,return_data_tag1,MPI_COMM_WORLD,ierr)


c           call MPI_SEND(NIas,1,MPI_INT,root_process,return_data_tag2
c     *             ,MPI_COMM_WORLD,ierr)
c            

c            call MPI_SEND(enwavgs(1),iTmins,MPI_DOUBLE_PRECISION,
c     *          root_process,return_data_tag3,MPI_COMM_WORLD,ierr)

 
         endif

         call MPI_FINALIZE(ierr)
         stop

      end

      double precision function bvel(t)
         double precision t
         double precision r,L,A,w,PI,phase

         common /heightparam/ r,L,A,w,PI,phase

         bvel=w*A*dcos(w*t+phase)
         return
      endfunction

      subroutine funcd(t,f,fd)
         double precision t,f,fd
         double precision r,L,A,w,PI,phase
         double precision z1,w1,total_time

         common /heightparam/ r,L,A,w,PI,phase
         common /funcdparam/ z1,w1,total_time


         f=z1+t*w1-A*dsin(w*(t+total_time)+phase)
         fd=w1-A*w*dcos(w*(t+total_time)+phase)
         return
      end

      double precision function funcheight(t)
         double precision t
         double precision r,L,A,w,PI,phase
         double precision z1,w1,total_time

         common /heightparam/ r,L,A,w,PI,phase
         common /funcdparam/ z1,w1,total_time

         funcheight=z1+t*w1-A*dsin(w*(t+total_time)+phase)
         return
      end
