c      program gafortran
      subroutine USERSUB ( parentv, funcvalv, iv)
c
c  This is version 1.6.4, last updated on 1/6/97.
c  
c  nce mod start 3 Mar 97.
c  nce mod last   9 Sep 97.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c	USE DFPORT 
c	USE DFLIB
c
      implicit double precision (a-h,o-z)
c      save
c
      include 'params.f'
      dimension parent(indmax,nparmax),child(indmax,nparmax)
      dimension fitness(indmax),nposibl(nparmax),nichflg(nparmax)
      dimension iparent(indmax,nchrmax),ichild(indmax,nchrmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension ibest(nchrmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
c
      REAL*8 parentv(indmax,nparmax),funcvalv(indmax),iv
      INTEGER*4 res
c      INTEGER*2 status
c      CHARACTER*80 buf
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga7   / child,ichild
      common / ga8   / nichflg
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
c  Call the input subroutine.
c      TIME0=SECNDS(0.0)
      call input
c
c  $$$$$ Main generational processing loop. $$$$$
c	CALL GETARG(1, buf, status)
c	read (buf,*) iv
c	WRITE(*,*) iv
c
      if (iv.eq.1) then
c         res=system("del gabig.out")
c         res=system("copy gabig.orig gabig.out")
         res=system("rm gabig.out")
         res=system("touch gabig.out")
         irestrt=0
      else         
         irestrt=1
      endif
c
c  Perform necessary initialization and read the ga.restart file.
      i=iv
      call initial(istart,npossum,ig2sum)
c      kount=0
c      do 20 i=istart,maxgen+istart-1
         write (6,1111) i
         write (24,1111) i
         write(24,1050)
c Skip if iv=-1.  If iv=-1 then we are restarting.
         if (iv.ne.-1) then
c
c  Evaluate the population, assign fitness, establish the best 
c  individual, and write output information. 
      do 31 j=1,npopsiz
         fitness(j)=funcvalv(j)
 31   continue
c
         call evalout(iskip,iend,ibest)
         if(npopsiz.eq.1 .or. iskip.ne.0) then
            close(24)
c            stop
         endif
c
c  Implement "niching".
         if (iniche.ne.0) call niche
c
c  Enter selection, crossover and mutation loop.
         ncross=0
         ipick=npopsiz
         do 45 j=1,npopsiz,nchild
c
c  Perform selection.
            call selectn(ipick,j,mate1,mate2)
c
c  Now perform crossover between the randomly selected pair.
            call crosovr(ncross,j,mate1,mate2)
 45      continue
         write(6,1225) ncross
         write(24,1225) ncross
c
c  Now perform random mutations.  If running micro-GA, skip mutation.
         if (microga.eq.0) call mutate
c
c  Write child array back into parent array for new generation.  Check
c  to see if the best parent was replicated.
         call newgen(ielite,npossum,ig2sum,ibest)
c
c  Implement micro-GA if enabled.
         if (microga.ne.0) call gamicro(i,npossum,ig2sum,ibest)
c
c  Write to restart file.
         call restart(i,istart,kount)
c         CPU = SECNDS(TIME0)
c         write(6,1400) cpu
c         write(24,1400) cpu
c 20   continue
c  $$$$$ End of main generational processing loop. $$$$$
c 999  continue
c
c Decode parameters for CodeV
c
      endif
      do 32 j=1,npopsiz
         call decode(j,parentv,iparent)
 32   continue
c
      CLOSE (24)
c
c      res=system("type ga.out>>gabig.out")
      res=system("cat ga.out>>gabig.out")
c
c
 1050 format(1x,' #      Binary Code',16x,'Param1  Param2  Fitness')
 1111 format(//'#################  Generation',i5,'  #################')
 1225 format(/'  Number of Crossovers      =',i5)
c 1400 format(2x,'CPU time for generation=',e10.4)
c
c      stop
      end
c
c#######################################################################
      subroutine input
c
c  This subroutine inputs information from the ga.inp file. 
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension nposibl(nparmax),nichflg(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga8   / nichflg
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
      namelist / ga   / irestrt,npopsiz,pmutate,maxgen,idum,pcross,
     +                  itourny,ielite,icreep,pcreep,iunifrm,iniche,
     +                  iskip,iend,nchild,nparam,parmin,parmax,nposibl,
     +                  nowrite,nichflg,microga,kountmx
c
      kountmx=5
      irestrt=0
      itourny=0
      ielite=0
      iunifrm=0
      iniche=0
      iskip=0
      iend=0
      nchild=1
      do 2 i=1,nparam
         nichflg(i)=1
 2    continue
      microga=0
c
      OPEN (UNIT=24, FILE='ga.out', STATUS='UNKNOWN')
      rewind 24
      OPEN (UNIT=23, FILE='ga.inp', STATUS='OLD')
      READ (23, NML = ga)
      CLOSE (23)
      itourny=1
c      if (itourny.eq.0) nchild=2
c
c  Check for array sizing errors.
      if (npopsiz.gt.indmax) then
         write(6,1600) npopsiz
         write(24,1600) npopsiz
         close(24)
c         stop
      endif
      if (nparam.gt.nparmax) then
         write(6,1700) nparam
         write(24,1700) nparam
         close(24)
c         stop
      endif
c
c  If using the microga option, reset some input variables
      if (microga.ne.0) then
         pmutate=0.0
         pcreep=0.0
         itourny=1
         ielite=1
         iniche=0
         nchild=1
         if (iunifrm.eq.0) then
            pcross=1.0
         else
            pcross=0.5
         endif
      endif
c
 1600 format(1x,'ERROR: npopsiz > indmax.  Set indmax = ',i6)
 1700 format(1x,'ERROR: nparam > nparmax.  Set nparmax = ',i6)
c
      return
      end
c
c#######################################################################
      subroutine initial(istart,npossum,ig2sum)
c
c  This subroutine sets up the program by generating the g0, g1 and 
c  ig2 arrays, and counting the number of chromosomes required for the
c  specified input.  The subroutine also initializes the random number 
c  generator, parent and iparent arrays (reads the ga.restart file).
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension parent(indmax,nparmax),iparent(indmax,nchrmax)
      dimension nposibl(nparmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
c
      do 3 i=1,nparam
         g0(i)=parmin(i)
         pardel(i)=parmax(i)-parmin(i)
         g1(i)=pardel(i)/dble(nposibl(i)-1)
 3    continue
      do 6 i=1,nparam
         do 7 j=1,30
            n2j=2**j
            if (n2j.ge.nposibl(i)) then
               ig2(i)=j
               goto 8
            endif
            if (j.ge.30) then
               write(6,2000)
               write(24,2000)
               close(24)
c               stop
            endif
 7       continue
 8       continue
 6    continue
c
c  Count the total number of chromosomes (bits) required
      nchrome=0
      npossum=0
      ig2sum=0
      do 9 i=1,nparam
         nchrome=nchrome+ig2(i)
         npossum=npossum+nposibl(i)
         ig2sum=ig2sum+(2**ig2(i))
 9    continue
      if (nchrome.gt.nchrmax) then
         write(6,1800) nchrmax
         write(24,1800) nchrome
         close(24)
c         stop
      endif
c
      if (npossum.lt.ig2sum .and. microga.ne.0) then
         write(6,2100)
         write(24,2100)
      endif
c
c  Initialize random number generator
      call ran3(idum,rand)
c
      if(irestrt.eq.0) then
c  Initialize the random distribution of parameters in the individual
c  parents when irestrt=0.
c	   write (24,*) "NO read ga.restart."
         istart=1
         do 10 i=1,npopsiz
            do 15 j=1,nchrome
               call ran3(1,rand)
               iparent(i,j)=1
               if(rand.lt.0.5) iparent(i,j)=0
 15         continue
 10      continue
         if (npossum.lt.ig2sum) call possibl(parent,iparent)
      else
c  If irestrt.ne.0, read from restart file.
c	   write (24,*) "READ ga.restart."
         OPEN (UNIT=25, FILE='ga.restart', STATUS='OLD')
         rewind 25
         read(25,*) istart,npopsiz
         do 1 j=1,npopsiz
            read(25,*) k,(iparent(j,l),l=1,nchrome)
 1       continue
         CLOSE (25)
      endif
c
      if(irestrt.ne.0) call ran3(idum-istart,rand)
c
 1800 format(1x,'ERROR: nchrome > nchrmax.  Set nchrmax = ',i6)
 2000 format(1x,'ERROR: You have a parameter with a number of '/
     +       1x,'   possibilities > 2**30!  If you really desire this,'/
     +       1x,'   change the DO loop 7 statement and recompile.'//
     +       1x,'   You may also need to alter the code to work with'/
     +       1x,'   REAL numbers rather than INTEGER numbers; Fortran'/
     +       1x,'   does not like to compute 2**j when j>30.')
 2100 format(1x,'WARNING: for some cases, a considerable performance'/
     +       1x,'   reduction has been observed when running a non-'/
     +       1x,'   optimal number of bits with the micro-GA.'/
     +       1x,'   If possible, use values for nposibl of 2**n,'/
     +       1x,'   e.g. 2, 4, 8, 16, 32, 64, etc.  See ReadMe file.')
c
      return
      end
c
c#######################################################################
      subroutine evalout(iskip,iend,ibest)
c
c  This subroutine evaluates the population, assigns fitness, 
c  establishes the best individual, and outputs information. 
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension parent(indmax,nparmax),iparent(indmax,nchrmax)
      dimension fitness(indmax)
      dimension paramsm(nparmax),paramav(nparmax),ibest(nchrmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
c
      fitsum=0.0
      best=-1.0d10
      do 29 n=1,nparam
         paramsm(n)=0.0
 29   continue
      jstart=1
      jend=npopsiz
      if(iskip.ne.0) jstart=iskip
      if(iend.ne.0) jend=iend
      do 30 j=jstart,jend
         call decode(j,parent,iparent)
         if(iskip.ne.0 .and. iend.ne.0 .and. iskip.eq.iend)
     +   write(6,1075) j,(iparent(j,k),k=1,nchrome),
     +                   (parent(j,kk),kk=1,nparam),0.0
c
c  Call function evaluator, write out individual and fitness, and add
c  to the summation for later averaging.
c         call func(j,funcval)
c         fitness(j)=funcval
c         write(24,1075) j,(iparent(j,k),k=1,nchrome),
c     +                  (parent(j,kk),kk=1,nparam),fitness(j)
         write(24,1075) j,
     +                  (parent(j,kk),kk=1,nparam),fitness(j)
         fitsum=fitsum+fitness(j)
         do 22 n=1,nparam
            paramsm(n)=paramsm(n)+parent(j,n)
 22      continue
c
c  Check to see if fitness of individual j is the best fitness.
         if (fitness(j).gt.best) then
            best=fitness(j)
            jbest=j
            do 24 k=1,nchrome
               ibest(k)=iparent(j,k)
 24         continue
         endif
 30   continue
c
c  Compute parameter and fitness averages.
      fbar=fitsum/dble(npopsiz)
      do 23 n=1,nparam
         paramav(n)=paramsm(n)/dble(npopsiz)
 23   continue
c
c  Write output information
      if (npopsiz.eq.1) then
         write(24,1075) 1,(iparent(1,k),k=1,nchrome),
     +                  (parent(1,k),k=1,nparam),fitness(1)
         write(24,*) ' Average Values:'
         write(24,1275) (parent(1,k),k=1,nparam),fbar
      else
         write(24,1275) (paramav(k),k=1,nparam),fbar
      endif
      write(6,1100) fbar
      write(24,1100) fbar
      write(6,1200) best
      write(24,1200) best
c
 1075 format(i3, 1x, 9(1x, f8.4), 2x, f10.0)
 1100 format(1x,'Average Function Value of Generation=',f14.3)
 1200 format(1x,'Maximum Function Value              =',f14.3/)
 1275 format(/' Average Values:'/,4x,9(1x, f8.4),2x,f10.0/
     +     )
      return
      end
c
c#######################################################################
      subroutine niche
c
c  Implement "niching" through Goldberg's multidimensional phenotypic 
c  sharing scheme with a triangular sharing function.  To find the
c  multidimensional distance from the best individual, normalize all
c  parameter differences.
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension parent(indmax,nparmax),iparent(indmax,nchrmax)
      dimension fitness(indmax),nposibl(nparmax),nichflg(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga8   / nichflg
c
c   Variable definitions:
c
c  alpha   = power law exponent for sharing function; typically = 1.0
c  del     = normalized multidimensional distance between ii and all
c            other members of the population
c            (equals the square root of del2)
c  del2    = sum of the squares of the normalized multidimensional 
c            distance between member ii and all other members of 
c            the population
c  nniche  = number of niched parameters
c  sigshar = normalized distance to be compared with del; in some sense,
c            1/sigshar can be viewed as the number of regions over which
c            the sharing function should focus, e.g. with sigshar=0.1, 
c            the sharing function will try to clump in ten distinct  
c            regions of the phase space.  A value of sigshar on the
c            order of 0.1 seems to work best.
c  share   = sharing function between individual ii and j
c  sumshar = sum of the sharing functions for individual ii
c
      alpha=1.0
      sigshar=0.1
      nniche=0
      do 33 jj=1,nparam
         nniche=nniche+nichflg(jj)
 33   continue
      if (nniche.eq.0) then
         write(6,1900)
         write(24,1900)
         close(24)
c         stop
      endif
      do 34 ii=1,npopsiz
         sumshar=0.0
         do 35 j=1,npopsiz
            del2=0.0
            do 36 k=1,nparam
               if (nichflg(k).ne.0) then
                  del2=del2+((parent(j,k)-parent(ii,k))/pardel(k))**2.0
               endif
 36         continue
            del=(dsqrt(del2))/dble(nniche)
            if (del.lt.sigshar) then
               share=1.0-((del/sigshar)**alpha)
            else
               share=0.0
            endif
            sumshar=sumshar+share/dble(npopsiz)
 35      continue
         if (sumshar.ne.0.0) fitness(ii)=fitness(ii)/sumshar
 34   continue
c
 1900 format(1x,'ERROR: iniche=1 and all values in nichflg array = 0'/
     +       1x,'       Do you want to niche or not?')
c
      return
      end
c
c#######################################################################
      subroutine selectn(ipick,j,mate1,mate2)
c
c  Subroutine for selection operator.  Presently, tournament selection
c  is the only option available.
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension parent(indmax,nparmax),child(indmax,nparmax)
      dimension fitness(indmax)
      dimension iparent(indmax,nchrmax),ichild(indmax,nchrmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
c  If tournament selection is chosen (i.e. itourny=1), then
c  implement "tournament" selection for selection of new population.
      if(itourny.eq.1) then
         call select(mate1,ipick)
         call select(mate2,ipick)
c        write(3,*) mate1,mate2,fitness(mate1),fitness(mate2)
         do 46 n=1,nchrome
            ichild(j,n)=iparent(mate1,n)
            if(nchild.eq.2) ichild(j+1,n)=iparent(mate2,n)
 46      continue
      endif
c
      return
      end
c
c#######################################################################
      subroutine crosovr(ncross,j,mate1,mate2)
c
c  Subroutine for crossover between the randomly selected pair.
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension parent(indmax,nparmax),child(indmax,nparmax)
      dimension iparent(indmax,nchrmax),ichild(indmax,nchrmax)
c
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
      if (iunifrm.eq.0) then
c  Single-point crossover at a random chromosome point.
         call ran3(1,rand)
         if(rand.gt.pcross) goto 69
         ncross=ncross+1
         call ran3(1,rand)
         icross=2+dint(dble(nchrome-1)*rand)
         do 50 n=icross,nchrome
            ichild(j,n)=iparent(mate2,n)
            if(nchild.eq.2) ichild(j+1,n)=iparent(mate1,n)
 50      continue
      else
c  Perform uniform crossover between the randomly selected pair.
         do 60 n=1,nchrome
            call ran3(1,rand)
            if(rand.le.pcross) then
               ncross=ncross+1
               ichild(j,n)=iparent(mate2,n)
               if(nchild.eq.2) ichild(j+1,n)=iparent(mate1,n)
            endif
 60      continue
      endif
 69   continue
c
      return
      end
c
c#######################################################################
      subroutine mutate
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension nposibl(nparmax)
      dimension child(indmax,nparmax),ichild(indmax,nchrmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      common / ga7   / child,ichild
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
c  This subroutine performs mutations on the children generation.
c  Perform random jump mutation if a random number is less than pmutate.
c  Perform random creep mutation if a different random number is less
c  than pcreep.  
      nmutate=0
      ncreep=0
      do 70 j=1,npopsiz
         do 75 k=1,nchrome
c  Jump mutation
            call ran3(1,rand)
            if (rand.le.pmutate) then
               nmutate=nmutate+1
               if(ichild(j,k).eq.0) then
                  ichild(j,k)=1
               else
                  ichild(j,k)=0
               endif
               if (nowrite.eq.0) write(6,1300) j,k
               if (nowrite.eq.0) write(24,1300) j,k
            endif
 75      continue
c  Creep mutation (one discrete position away).
         if (icreep.ne.0) then
            do 76 k=1,nparam
               call ran3(1,rand)
               if(rand.le.pcreep) then
                  call decode(j,child,ichild)
                  ncreep=ncreep+1
                  creep=1.0
                  call ran3(1,rand)
                  if (rand.lt.0.5) creep=-1.0
                  child(j,k)=child(j,k)+g1(k)*creep
                  if (child(j,k).gt.parmax(k)) then
                     child(j,k)=parmax(k)-1.0*g1(k)
                  elseif (child(j,k).lt.parmin(k)) then
                     child(j,k)=parmin(k)+1.0*g1(k)
                  endif
                  call code(j,k,child,ichild)
                  if (nowrite.eq.0) write(6,1350) j,k
                  if (nowrite.eq.0) write(24,1350) j,k
               endif
 76         continue
         endif
 70   continue
      write(6,1250) nmutate,ncreep
      write(24,1250) nmutate,ncreep
c
 1250 format(/'  Number of Jump Mutations  =',i5/
     +        '  Number of Creep Mutations =',i5)
 1300 format('*** Jump mutation performed on individual  ',i4,
     +       ', chromosome ',i3,' ***')
 1350 format('*** Creep mutation performed on individual ',i4,
     +       ', parameter  ',i3,' ***')
c
      return
      end
c
c#######################################################################
      subroutine newgen(ielite,npossum,ig2sum,ibest)
c
c  Write child array back into parent array for new generation.  Check
c  to see if the best parent was replicated; if not, and if ielite=1,
c  then reproduce the best parent into a random slot.
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension parent(indmax,nparmax),child(indmax,nparmax)
      dimension iparent(indmax,nchrmax),ichild(indmax,nchrmax)
      dimension ibest(nchrmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga7   / child,ichild
c
      if (npossum.lt.ig2sum) call possibl(child,ichild)
      kelite=0
      do 94 j=1,npopsiz
         jelite=0
         do 95 n=1,nchrome
            iparent(j,n)=ichild(j,n)
            if (iparent(j,n).eq.ibest(n)) jelite=jelite+1
            if (jelite.eq.nchrome) kelite=1
 95      continue
 94   continue
      if (ielite.ne.0 .and. kelite.eq.0) then
         call ran3(1,rand)
         irand=1+dint(dble(npopsiz)*rand)
         do 96 n=1,nchrome
            iparent(irand,n)=ibest(n)
 96      continue
         write(24,1260) irand
      endif
c
 1260 format('  Elitist Reproduction on Individual ',i4)
c
      return
      end
c
c#######################################################################
      subroutine gamicro(i,npossum,ig2sum,ibest)
c
c  Micro-GA implementation subroutine
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension parent(indmax,nparmax),iparent(indmax,nchrmax)
      dimension ibest(nchrmax)
c
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
c
c  First, check for convergence of micro population.
c  If converged, start a new generation with best individual and fill
c  the remainder of the population with new randomly generated parents.
c
c  Count number of different bits from best member in micro-population
      icount=0
      do 81 j=1,npopsiz
         do 82 n=1,nchrome
            if(iparent(j,n).ne.ibest(n)) icount=icount+1
 82      continue
 81   continue
c
c  If icount less than 5% of number of bits, then consider population
c  to be converged.  Restart with best individual and random others.
      diffrac=dble(icount)/dble((npopsiz-1)*nchrome)
      if (diffrac.lt.0.05) then
      do 87 n=1,nchrome
         iparent(1,n)=ibest(n)
 87   continue
      do 88 j=2,npopsiz
         do 89 n=1,nchrome
            call ran3(1,rand)
            iparent(j,n)=1
            if(rand.lt.0.5) iparent(j,n)=0
 89      continue
 88   continue
      if (npossum.lt.ig2sum) call possibl(parent,iparent)
      write(6,1375) i
      write(24,1375) i
      endif
c
 1375 format(//'%%%%%%%  Restart micro-population at generation',
     +       i5,'  %%%%%%%')
c
      return
      end
c
c#######################################################################
      subroutine select(mate,ipick)
c
c  This routine selects the better of two possible parents for mating.
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      dimension parent(indmax,nparmax),iparent(indmax,nchrmax)
      dimension fitness(indmax)
c
      if(ipick+1.gt.npopsiz) call shuffle(ipick)
      ifirst=ipick
      isecond=ipick+1
      ipick=ipick+2
      if(fitness(ifirst).gt.fitness(isecond)) then
         mate=ifirst
      else
         mate=isecond
      endif
c     write(3,*)'select',ifirst,isecond,fitness(ifirst),fitness(isecond)
c
      return
      end
c
c#######################################################################
      subroutine shuffle(ipick)
c
c  This routine shuffles the parent array and its corresponding fitness
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      common / ga4   / fitness
      dimension parent(indmax,nparmax),iparent(indmax,nchrmax)
      dimension fitness(indmax)
c
      ipick=1
      do 10 j=1,npopsiz-1
         call ran3(1,rand)
         iother=j+1+dint(dble(npopsiz-j)*rand)
         do 20 n=1,nchrome
            itemp=iparent(iother,n)
            iparent(iother,n)=iparent(j,n)
            iparent(j,n)=itemp
 20      continue
         temp=fitness(iother)
         fitness(iother)=fitness(j)
         fitness(j)=temp
 10   continue
c
      return
      end
c
c#######################################################################
      subroutine decode(i,array,iarray)
c
c  This routine decodes a binary string to a real number.
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      dimension array(indmax,nparmax),iarray(indmax,nchrmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
c
      l=1
      do 10 k=1,nparam
         iparam=0
         m=l
         do 20 j=m,m+ig2(k)-1
            l=l+1
            iparam=iparam+iarray(i,j)*(2**(m+ig2(k)-1-j))
 20      continue
         array(i,k)=g0(k)+g1(k)*dble(iparam)
 10   continue
c
      return
      end
c
c#######################################################################
      subroutine code(j,k,array,iarray)
c
c  This routine codes a parameter into a binary string.
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      dimension array(indmax,nparmax),iarray(indmax,nchrmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax)
c
c  First, establish the beginning location of the parameter string of
c  interest.
      istart=1
      do 10 i=1,k-1
         istart=istart+ig2(i)
 10   continue
c
c  Find the equivalent coded parameter value, and back out the binary
c  string by factors of two.
      m=ig2(k)-1
      if (g1(k).eq.0.0) return
      iparam=nint((array(j,k)-g0(k))/g1(k))
      do 20 i=istart,istart+ig2(k)-1
         iarray(j,i)=0
         if ((iparam+1).gt.(2**m)) then
            iarray(j,i)=1
            iparam=iparam-2**m
         endif
         m=m-1
 20   continue
c     write(3,*)array(j,k),iparam,(iarray(j,i),i=istart,istart+ig2(k)-1)
c
      return
      end
c
c#######################################################################
c
      subroutine possibl(array,iarray)
c
c  This subroutine determines whether or not all parameters are within
c  the specified range of possibility.  If not, the parameter is
c  randomly reassigned within the range.  This subroutine is only
c  necessary when the number of possibilities per parameter is not
c  optimized to be 2**n, i.e. if npossum < ig2sum.
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga5   / g0,g1,ig2
      common / ga6   / parmax,parmin,pardel,nposibl
      dimension array(indmax,nparmax),iarray(indmax,nchrmax)
      dimension g0(nparmax),g1(nparmax),ig2(nparmax),nposibl(nparmax)
      dimension parmax(nparmax),parmin(nparmax),pardel(nparmax)
c
      do 10 i=1,npopsiz
         call decode(i,array,iarray)
         do 20 j=1,nparam
            n2ig2j=2**ig2(j)
            if(nposibl(j).ne.n2ig2j .and. array(i,j).gt.parmax(j)) then
               call ran3(1,rand)
               irand=dint(dble(nposibl(j))*rand)
               array(i,j)=g0(j)+dble(irand)*g1(j)
               call code(i,j,array,iarray)
               if (nowrite.eq.0) write(6,1000) i,j
               if (nowrite.eq.0) write(24,1000) i,j
            endif
 20      continue
 10   continue
c
 1000 format('*** Parameter adjustment to individual     ',i4,
     +       ', parameter  ',i3,' ***')
c
      return
      end
c
c#######################################################################
      subroutine restart(i,istart,kount)
c
c  This subroutine writes restart information to the ga.restart file. 
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      common / ga1   / npopsiz,nowrite
      common / ga2   / nparam,nchrome
      common / ga3   / parent,iparent
      dimension parent(indmax,nparmax),iparent(indmax,nchrmax)
      common /inputga/ pcross,pmutate,pcreep,maxgen,idum,irestrt,
     +                 itourny,ielite,icreep,iunifrm,iniche,
     +                 iskip,iend,nchild,microga,kountmx
c
      kount=kount+1
c      if(i.eq.maxgen+istart-1 .or. kount.eq.kountmx) then
         OPEN (UNIT=25, FILE='ga.restart', STATUS='OLD')
         rewind 25
         write(25,*) i+1,npopsiz
         do 80 j=1,npopsiz
            write(25,1500) j,(iparent(j,l),l=1,nchrome)
 80      continue
         CLOSE (25)
         kount=0
c      endif
c
 1500 format(i5,3x,30i2)
c
      return
      end
c
c#######################################################################  
      subroutine ran3(idum,rand)
c
c  Returns a uniform random deviate between 0.0 and 1.0.  Set idum to 
c  any negative value to initialize or reinitialize the sequence.
c  This function is taken from W.H. Press', "Numerical Recipes" p. 199.
c
      implicit double precision (a-h,m,o-z)  
      save
c      implicit real*4(m)
      parameter (mbig=4000000.,mseed=1618033.,mz=0.,fac=1./mbig)
c     parameter (mbig=1000000000,mseed=161803398,mz=0,fac=1./mbig)
c
c  According to Knuth, any large mbig, and any smaller (but still large)
c  mseed can be substituted for the above values.
      dimension ma(55)
      data iff /0/
      if (idum.lt.0 .or. iff.eq.0) then
         iff=1
         mj=mseed-dble(iabs(idum))
         mj=dmod(mj,mbig)
         ma(55)=mj
         mk=1
         do 11 i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk.lt.mz) mk=mk+mbig
            mj=ma(ii)
 11      continue
         do 13 k=1,4
            do 12 i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i).lt.mz) ma(i)=ma(i)+mbig
 12         continue
 13      continue
         inext=0
         inextp=31
         idum=1
      endif
      inext=inext+1
      if(inext.eq.56) inext=1
      inextp=inextp+1
      if(inextp.eq.56) inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj.lt.mz) mj=mj+mbig
      ma(inext)=mj
      rand=mj*fac
      return
      end
c
c#######################################################################
c    
      subroutine func(j,funcval)
c
      implicit double precision (a-h,o-z)
      save
c
      include 'params.f'
      dimension parent(indmax,nparmax)
      dimension iparent(indmax,nchrmax)
      real*8 funcvalv(indmax)
c
      common / ga3   / parent,iparent
c      common / ga9   / funcvalv
c
c  This is a two-dimensional version of the multimodal function with
c  decreasing peaks used by Goldberg and Richardson (1987, see ReadMe
c  file for complete reference).  In two dimensions, this function has
c  25 peaks, but only one global maximum.  It is a reasonably tough
c  problem for the GA.  The uniform crossover micro-GA does well, the
c  single-point crossover micro-GA is slower, more conventional GA
c  techniques require reasonably large population sizes.
c
c      pi=4.0*atan(1.d0)
c      f11=(sin(5.1*pi*parent(j,1) + 0.5))**6.0
c      f12=exp(-4.0*log(2.0)*((parent(j,1)-0.0667)**2.0)/0.64)
c      
c      f21=(sin(5.1*pi*parent(j,2) + 0.5))**6.0
c      f22=(exp(-4.0*log(2.0)*((parent(j,2)-0.0667)**2.0)/0.64))
c      
c      funcval=f11*f12*f21*f22
c
c
c  Transfer fitness values from CodeV to here.
      funcval=funcvalv(j)
c
      return
      end
c#######################################################################
