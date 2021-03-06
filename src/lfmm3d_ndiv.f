c
       subroutine lfmm3d(nd,eps,nsource,source,ifcharge,
     $    charge,ifdipole,dipvec,ifpgh,pot,grad,hess,ntarg,
     $    targ,ifpghtarg,pottarg,gradtarg,hesstarg,ndiv,idivflag,
     $    ifnear)
c
c        Laplace FMM in R^{3}: evaluate all pairwise particle
c        interactions (ignoring self-interactions) and interactions
c        with targs.
c
c        We use (1/r) for the Green's function, without the
c        1/(4 \pi) scaling.
c
c
c        Input parameters:
c
c   nd:   number of densities
c
c   eps:  requested precision
c
c   nsource in: integer  
c                number of sources
c
c   source  in: double precision (3,nsource)
c                source(k,j) is the kth component of the jth
c                source locations
c
c   ifcharge  in: integer  
c             charge computation flag
c              ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c 
c   charge    in: double precision (nsource) 
c              charge strengths
c
c   ifdipole   in: integer
c              dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c
c
c   dipvec   in: double precision (3,nsource) 
c              dipole orientation vectors
c
c   ifpgh   in: integer
c              flag for evaluating potential/gradient at the sources
c              ifpgh = 1, only potential is evaluated
c              ifpgh = 2, potential and gradients are evaluated
c
c   ntarg  in: integer  
c                 number of targs 
c
c   targ  in: double precision (3,ntarg)
c               targ(k,j) is the kth component of the jth
c               targ location
c
c   ifpghtarg   in: integer
c              flag for evaluating potential/gradient at the targs
c              ifpghtarg = 1, only potential is evaluated
c              ifpghtarg = 2, potential and gradient are evaluated
c
c
c     OUTPUT parameters:
c
c
c   pot:    out: double precision(nd,nsource) 
c               potential at the source locations
c
c   grad:   out: double precision(nd,3,nsource)
c               gradient at the source locations
c
c   hess    out: double precision(nd,6,nsource)
c               hessian at the source locations
c                 (currently not supported)
c
c   pottarg:    out: double precision(nd,ntarg) 
c               potential at the targ locations
c
c   gradtarg:   out: double precision(nd,3,ntarg)
c               gradient at the targ locations
c
c   hesstarg    out: double precision(nd,6,ntarg)
c                hessian at the target locations - currently not
c                supported
     
       implicit none

       integer nd

       double precision eps

       integer ifcharge,ifdipole
       integer ifpgh,ifpghtarg

       integer ntarg,nsource
       

       double precision source(3,*),targ(3,*)
       double precision charge(nd,*)
       double precision dipvec(nd,3,*)

       double precision pot(nd,*),grad(nd,3,*),hess(nd,6,*)
       double precision pottarg(nd,*),gradtarg(nd,3,*),hesstarg(nd,6,*)

c
cc       tree variables
c
       integer mhung,idivflag,ndiv,isep,nboxes,nbmax,nlevels
       integer *8 ltree
       integer nlmax
       integer mnbors,mnlist1,mnlist2,mnlist3,mnlist4
       integer *8 ipointer(32)
       integer, allocatable :: itree(:)
       double precision, allocatable :: treecenters(:,:),boxsize(:)

c
cc       temporary sorted arrays
c
       double precision, allocatable :: sourcesort(:,:),targsort(:,:)
       double precision, allocatable :: radsrc(:)
       double precision, allocatable :: chargesort(:,:)
       double precision, allocatable :: dipvecsort(:,:,:)

       double precision, allocatable :: potsort(:,:),gradsort(:,:,:)
       double precision, allocatable :: hesssort(:,:,:)
       double precision, allocatable :: pottargsort(:,:)
       double precision, allocatable :: gradtargsort(:,:,:)
       double precision, allocatable :: hesstargsort(:,:,:)
c
cc        temporary fmm arrays
c
       double precision epsfmm
       integer, allocatable :: nterms(:)
       integer *8, allocatable :: iaddr(:,:)
       double precision, allocatable :: scales(:)
       double precision, allocatable :: rmlexp(:)

       integer lmptemp,nmax
       integer *8 lmptot
       double precision, allocatable :: mptemp(:),mptemp2(:)

       integer ifnear

c
cc       temporary variables not main fmm routine but
c        not used in particle code
       double precision expc(3),scjsort(1),radexp
       double complex texpssort(100)
       double precision expcsort(3)
       integer ntj,nexpc,nadd

c
cc         other temporary variables
c
        integer i,iert,ifprint,ilev,idim,ier
        double precision time1,time2,omp_get_wtime,second

c
c
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c       
      ifprint=1

c
cc        figure out tree structure
c
c
cc         set criterion for box subdivision
c


c
cc      set tree flags
c 
       isep = 1
       nlmax = 200
       nlevels = 0
       nboxes = 0
       mhung = 0
       ltree = 0

       nexpc = 0
       nadd = 0
       ntj = 0

       mnlist1 = 0
       mnlist2 = 0
       mnlist3 = 0
       mnlist4 = 0
       nbmax = 0


       allocate(radsrc(nsource))
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)       
       do i=1,nsource
          radsrc(i) = 0
       enddo
C$OMP END PARALLEL DO  

       radexp = 0

c
cc     memory management code for contructing level restricted tree
        iert = 0
        call mklraptreemem(iert,source,nsource,radsrc,targ,ntarg,
     1        expc,nexpc,radexp,idivflag,ndiv,isep,nlmax,nbmax,
     2        nlevels,nboxes,mnbors,mnlist1,mnlist2,mnlist3,
     3        mnlist4,mhung,ltree)

        if(ifprint.ge.1) print *, ltree/1.0d9


        if(iert.ne.0) then
           print *, "Error in allocating tree memory"
           stop
        endif


        allocate(itree(ltree))
        allocate(boxsize(0:nlevels))
        allocate(treecenters(3,nboxes))

c       Call tree code
        call mklraptree(source,nsource,radsrc,targ,ntarg,expc,
     1               nexpc,radexp,idivflag,ndiv,isep,mhung,mnbors,
     2               mnlist1,mnlist2,mnlist3,mnlist4,nlevels,
     2               nboxes,treecenters,boxsize,itree,ltree,ipointer)

c     Allocate sorted source and targ arrays      

      allocate(sourcesort(3,nsource))
      allocate(targsort(3,ntarg))
      if(ifcharge.eq.1) allocate(chargesort(nd,nsource))

      if(ifdipole.eq.1) then
         allocate(dipvecsort(nd,3,nsource))
      endif

      if(ifpgh.eq.1) then 
        allocate(potsort(nd,nsource),gradsort(nd,3,1),hesssort(nd,6,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,nsource),gradsort(nd,3,nsource),
     1       hesssort(nd,6,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,nsource),gradsort(nd,3,nsource),
     1       hesssort(nd,6,nsource))
      else
        allocate(potsort(nd,1),gradsort(nd,3,1),hesssort(nd,6,1))
      endif

      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,1),
     1      hesstargsort(nd,6,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,ntarg),
     1        hesstargsort(nd,6,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,ntarg),gradtargsort(nd,3,ntarg),
     1        hesstargsort(nd,6,ntarg))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,3,1),
     1     hesstargsort(nd,6,1))
      endif


c     scaling factor for multipole and local expansions at all levels
c
      allocate(scales(0:nlevels),nterms(0:nlevels))
      do ilev = 0,nlevels
        scales(ilev) = boxsize(ilev)
        scales(ilev) = 1.0d0
        scales(ilev) = 1.1d0
      enddo
c
cc      initialize potential and gradient at source
c       locations
c
      if(ifpgh.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpgh.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)

        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif


      if(ifpgh.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,nsource
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
            hesssort(idim,1,i) = 0
            hesssort(idim,2,i) = 0
            hesssort(idim,3,i) = 0
            hesssort(idim,4,i) = 0
            hesssort(idim,5,i) = 0
            hesssort(idim,6,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif



c
cc       initialize potential and gradient  at targ
c        locations
c
      if(ifpghtarg.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpghtarg.eq.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifpghtarg.eq.3) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)
        do i=1,ntarg
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
            hesstargsort(idim,1,i) = 0
            hesstargsort(idim,2,i) = 0
            hesstargsort(idim,3,i) = 0
            hesstargsort(idim,4,i) = 0
            hesstargsort(idim,5,i) = 0
            hesstargsort(idim,6,i) = 0
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

c     Compute length of expansions at each level      
      nmax = 0
      do i=0,nlevels
         call l3dterms(eps,nterms(i))
         if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))
      lmptemp = (nmax+1)*(2*nmax+1)*2*nd
      allocate(mptemp(lmptemp),mptemp2(lmptemp))

c
cc     reorder sources 
c
      call dreorderf(3,nsource,source,sourcesort,itree(ipointer(5)))
      if(ifcharge.eq.1) call dreorderf(nd,nsource,charge,chargesort,
     1                     itree(ipointer(5)))

      if(ifdipole.eq.1) then
         call dreorderf(3*nd,nsource,dipvec,dipvecsort,
     1       itree(ipointer(5)))
      endif

c
cc      reorder targs
c
      call dreorderf(3,ntarg,targ,targsort,itree(ipointer(6)))
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
      call mpalloc(nd,itree(ipointer(1)),iaddr,nlevels,lmptot,nterms)
      if(ifprint.ge. 1) print *, "lmptot =",lmptot/1.0d9


      allocate(rmlexp(lmptot),stat=ier)
      if(ier.ne.0) then
         print *, "Cannot allocate mpole expansion workspace"
         print *, "lmptot=", lmptot
         stop
      endif

c     Memory allocation is complete. 
c     Call main fmm routine

      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call lfmm3dmain(nd,eps,
     $   nsource,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipvecsort,
     $   ntarg,targsort,nexpc,expcsort,
     $   epsfmm,iaddr,rmlexp,lmptot,mptemp,mptemp2,lmptemp,
     $   itree,ltree,ipointer,isep,ndiv,nlevels,
     $   nboxes,boxsize,mnbors,mnlist1,mnlist2,mnlist3,mnlist4,
     $   scales,treecenters,itree(ipointer(1)),nterms,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,hesstargsort,ntj,
     $   texpssort,scjsort,ifnear)

      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)



      if(ifpgh.eq.1) then
        call dreorderi(nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
      endif
      if(ifpgh.eq.2) then 
        call dreorderi(nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
        call dreorderi(3*nd,nsource,gradsort,grad,
     1                 itree(ipointer(5)))
      endif

      if(ifpgh.eq.3) then 
        call dreorderi(nd,nsource,potsort,pot,
     1                 itree(ipointer(5)))
        call dreorderi(3*nd,nsource,gradsort,grad,
     1                 itree(ipointer(5)))
        call dreorderi(6*nd,nsource,hesssort,hess,
     1                 itree(ipointer(5)))
      endif


      if(ifpghtarg.eq.1) then
        call dreorderi(nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
        call dreorderi(3*nd,ntarg,gradtargsort,gradtarg,
     1     itree(ipointer(6)))
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(nd,ntarg,pottargsort,pottarg,
     1     itree(ipointer(6)))
        call dreorderi(3*nd,ntarg,gradtargsort,gradtarg,
     1     itree(ipointer(6)))
        call dreorderi(6*nd,ntarg,hesstargsort,hesstarg,
     1     itree(ipointer(6)))
      endif

      return
      end
