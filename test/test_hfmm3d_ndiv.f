      implicit real *8 (a-h,o-z)
      integer ndvals(100) 
      real *8 epsvals(100)
      complex *16 zkvals(100)
      integer ifcdvals(2,3)
      integer ifpghtvals(2,100)
      real *8 ntfac(100)
      integer nsvals(100),ntvals(100)
      integer, allocatable :: ndivs(:,:)
      real *8, allocatable :: src(:,:),targ(:,:)
      complex *16, allocatable :: charges(:),dipvec(:,:)
      complex *16, allocatable :: pot(:),pottarg(:),grad(:,:)
      complex *16, allocatable :: gradtarg(:,:)
      complex *16 hess(100),hesstarg(100)
      complex *16 ima
      character *90 fname
      

      done = 1
      pi = atan(done)*4

      ima = dcmplx(0.0d0,1.0d0)

      nnd = 6
      ndvals(1) = 1
      ndvals(2) = 2
      ndvals(3) = 3
      ndvals(4) = 4
      ndvals(5) = 6
      ndvals(6) = 8

      nep = 5
      epsvals(1) = 0.50001d-2
      epsvals(2) = 0.50001d-3
      epsvals(3) = 0.50001d-6
      epsvals(4) = 0.50001d-9
      epsvals(5) = 0.501d-12

      nndiv = 5
      allocate(ndivs(nndiv,nep))
      do i=1,2
        ndivs(1,i) = 100
        ndivs(2,i) = 200
        ndivs(3,i) = 300
        ndivs(4,i) = 400
        ndivs(5,i) = 500
      enddo

      do i=3,5
        ndivs(1,i) = 300
        ndivs(2,i) = 500
        ndivs(3,i) = 700
        ndivs(4,i) = 900
        ndivs(5,i) = 1100
      enddo

      
      nicd = 3
      ifcdvals(1,1) = 1
      ifcdvals(2,1) = 0

      ifcdvals(1,2) = 0
      ifcdvals(2,2) = 1

      ifcdvals(1,3) = 1
      ifcdvals(2,3) = 1

      ifnear = 1


c
c
c       8 distributions x 2 (pot/pot+grad) to be tested
c       idist = 1 - vol to vol uniform (s to s)
c       idist = 2 - vol to vol unif^2 (s to s)
c       idist = 3 - vol to bdry unif (s to t) 
c       idist = 4 - vol to bdry unif^2 (s to t) 
c       idist = 5 - bdry to bdry (s to s)
c       idist = 6 - bdry to bdry + vol unif (s to st)
c       idist = 7 - bdry to bdry + vol unif^2 (s to st)
c       idist = 8 - vol to vol + bdry unif (s to st)
c       idist = 9 - vol to vol + bdry unif^2 (s to st)
c
            
      nidist = 9
      ifpghtvals(1,1) = 1
      ifpghtvals(2,1) = 0
      ntfac(1) = 0

      ifpghtvals(1,2) = 1
      ifpghtvals(2,2) = 0
      ntfac(2) = 0

      ifpghtvals(1,3) = 0
      ifpghtvals(2,3) = 1
      ntfac(3) = 0.1d0 

      ifpghtvals(1,4) = 0
      ifpghtvals(2,4) = 1
      ntfac(4) = 0.1d0

      ifpghtvals(1,5) = 1
      ifpghtvals(2,5) = 0
      ntfac(5) = 0

      do i=6,9
        ifpghtvals(1,i) = 1
        ifpghtvals(2,i) = 1
        ntfac(i) = 10
      enddo

      ntfac(8) = 0.1d0
      ntfac(9) = 0.1d0

      idivflag = 2
      nns = 5
      nsvals(1) = 4000
      nsvals(2) = 8000 
      nsvals(3) = 16000
      nsvals(4) = 32000
      nsvals(5) = 64000

      nzk = 3
      zkvals(1) = 2*pi
      zkvals(2) = 8*2*pi
      zkvals(3) = 64*2*pi

 1110 format(2(2x,i3),2(2x,e11.5),2(2x,i6),6(2x,i1)) 
      do ind = 1,nnd
        nd = ndvals(ind)
        do iep = 1,nep
          eps = epsvals(iep)
          do izk=2,2
            zk = zkvals(izk)
            do icd = 1,nicd
              ifcharge = ifcdvals(1,icd)
              ifdipole = ifcdvals(2,icd)
              do insvals = 1,nns
                ns = nsvals(insvals)

                do ipg = 1,2
                  write(fname,'(a,i1,a,i1,a,i1,a,i1,a,i1,a,i1,a)') 
     1              "helm-res/ompon_nd",ind,"_iprec",iep,"_zk",izk,
     1              "_cd",icd,"_n",insvals,"_pg",ipg,".txt"
                  open(unit=33,file=trim(fname))
                  do idist = 1,nidist

                    nt = ntfac(idist)*ns
                    ifpgh = ipg*ifpghtvals(1,idist)
                    ifpght = ipg*ifpghtvals(2,idist)

                    allocate(src(3,ns),targ(3,nt),charges(ns))
                    allocate(dipvec(3,ns),pot(ns),grad(3,ns))
                    allocate(pottarg(nt),gradtarg(3,nt))

                    call get_src_targ(idist,ns,nt,src,targ)
                    do i=1,ns
                      charges(i) = hkrand(0) + ima*hkrand(0)
                      dipvec(1,i) = hkrand(0) + ima*hkrand(0)
                      dipvec(2,i) = hkrand(0) + ima*hkrand(0)
                      dipvec(3,i) = hkrand(0) + ima*hkrand(0)
                    enddo

                    do indiv = 1,nndiv
                      ndiv = ndivs(indiv,iep)
                      print *, idist,indiv,ndiv,ns,nt
                      call cpu_time(t1)
C$                      t1 = omp_get_wtime()
                      
                      call hfmm3d_ndiv(nd,eps,zk,ns,src,ifcharge,
     1                 charges,ifdipole,dipvec,ifpgh,pot,grad,hess,nt,
     2                 targ,ifpghtarg,pottarg,gradtarg,hesstarg,ndiv,
     3                 idivflag,ifnear)

                      call cpu_time(t2)
C$                      t2 = omp_get_wtime()                      

                      ttot = t2-t1
                      sp = (ns+nt)/(ttot)
                      write(33,1110) ndiv,idist,ttot,sp,ns,nt,nd,iep,
     1                    ifcharge,ifdipole,ifpgh,ifpght
                    enddo
cc                     end of indiv loop                    

                    deallocate(src,targ,charges,dipvec,pot,grad)
                    deallocate(pottarg,gradtarg)
                  enddo
                  close(33)
cc                  end of idist loop                  
                enddo
cc               end of ipg loop                
              enddo
cc             end of ins loop                
            enddo
cc             icd
          enddo
cc           izk
        enddo
cc         iep
      enddo
cc      ind

      stop
      end



      subroutine get_src_targ(idist,ns,nt,src,targ)
      implicit real *8 (a-h,o-z)
      real *8 src(3,ns),targ(3,nt)

      done = 1
      pi = atan(done)*4
      
      do i=1,ns
        src(1,i) = 2*hkrand(0)-1
        src(2,i) = 2*hkrand(0)-1
        src(3,i) = 2*hkrand(0)-1
      enddo

      do i=1,nt
        targ(1,i) = 2*hkrand(0)-1
        targ(2,i) = 2*hkrand(0)-1
        targ(3,i) = 2*hkrand(0)-1
      enddo


c       idist = 1 - vol to vol uniform (s to s)
c       idist = 2 - vol to vol unif^2 (s to s)
c       idist = 3 - vol to bdry unif (s to t) 
c       idist = 4 - vol to bdry unif^2 (s to t) 
c       idist = 5 - bdry to bdry (s to s)
c       idist = 6 - bdry to bdry + vol unif (s to st)
c       idist = 7 - bdry to bdry + vol unif^2 (s to st)
c       idist = 8 - vol to vol + bdry unif (s to st)
c       idist = 9 - vol to vol + bdry unif^2 (s to st)
c


      if(idist.eq.1.or.idist.eq.3.or.idist.eq.8) then
        do i=1,ns
          src(1,i) = 2*hkrand(0) - 1
          src(2,i) = 2*hkrand(0) - 1
          src(3,i) = 2*hkrand(0) - 1
        enddo

      endif


      if(idist.eq.2.or.idist.eq.4.or.idist.eq.9) then
        do i=1,ns
          src(1,i) = 2*hkrand(0)**2 - 1
          src(2,i) = 2*hkrand(0)**2 - 1
          src(3,i) = 2*hkrand(0)**2 - 1
        enddo

      endif

      if(idist.eq.5.or.idist.eq.6.or.idist.eq.7) then
        do i=1,ns
          thet = pi*hkrand(0)
          phi = thet*hkrand(0)
          src(1,i) = cos(thet)
          src(2,i) = 0.1*sin(thet)*cos(phi)
          src(3,i) = 0.3*sin(thet)*sin(phi)
        enddo
      endif

      if(idist.eq.3.or.idist.eq.4.or.idist.eq.8.or.idist.eq.9) then
        do i=1,nt
          thet = pi*hkrand(0)
          phi = thet*hkrand(0)
          targ(1,i) = cos(thet)
          targ(2,i) = 0.1*sin(thet)*cos(phi)/2
          targ(3,i) = 0.3*sin(thet)*sin(phi)/2
        enddo
      endif

      if(idist.eq.6) then
        do i=1,nt
          targ(1,i) = 2*hkrand(0)-1
          targ(2,i) = 2*hkrand(0)-1
          targ(3,i) = 2*hkrand(0)-1
        enddo
      endif

      if(idist.eq.7) then
        do i=1,nt
          targ(1,i) = 2*hkrand(0)**2-1
          targ(2,i) = 2*hkrand(0)**2-1
          targ(3,i) = 2*hkrand(0)**2-1
        enddo
      endif

      return
      end


       
      


      
