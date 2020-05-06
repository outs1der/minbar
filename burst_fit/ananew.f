      program ananew

c new on April 6, 2016:
c added provision to increase ifitstart by max 2 s to allow for odd deviation
c flip at start decay fit in power law in many short bursts. Find new code
c through variable 'iset2'
c     
c     ananew version 3, jz 2016-02-9
c     modification:
c     - added full output to cumulative file 'fit.results.txt'
c     - columns explained in 'fit.results.readme.txt'

      implicit none
      integer*4 i,j,k,n,ierr,nlc,op,wfc,bn,in,in2,imax,ipk,in3a,in3b,
     >     ipkarray(5),ifl,idu,signal_cygnus_x1,l
      parameter (nlc=100000)
      integer*2 nodata(nlc)
      real*8 t1(nlc), t2(nlc), peaktime,starttime,f1,f2,mjd,mjdc,tss(nlc),t0,
     >     starttime_grad,xmin,burstid_mjd
      real*4 ltf(nlc), x, y, pka(nlc), pkae(nlc), pkas(nlc),pk,pke,nustararea,
     >     pkarray(5),pkearray(5),diff,difs,chi,jemxarea,pcaarea,rms,s1,
     >     stderr, stdsig, decaytime, decaytimee, fitpeak1, fluence,delay,
     >     fitpeake1,fitpeak2,fitpeake2,fluencee,totexp,bg1,bge1,deltat,bgs,
     >     fpowbg,plindex,plindexe,bg2,bge2,bg,bge,chi3,d,fitpeak3,fitpeake3,
     >     d1,d2,pk5,pk5e,fluencet1,fluencete1,tt,tt2,ff1,ff2,pownorm,ft,
     >     tmax,chidif,tau,taue,chi2,rp(3),rpe(3),pi,fac,expon,ff3,wght(nlc),
     >     totpl2,totple2,totgs,totgse,agauss,pl,ple,bgo,bgoe,bg3,bge3,sh,
     >     plindex2,plindexe2,plindex3,plindexe3,chie,fte,grad,grade,gradm,
     >     fluencet2,fluencete2,fluencet3,fluencete3,totpl3,totple3,fracgs,
     >     fracgse,te1,te1e,threshold,te2,te2e,te3,te3e,duration,duratione,
     >     fluencef,fluencefe,tri,itri,hardstart,harddiff
      integer*4 pkas_sort(nlc), iread,ipk5,ist,instr,number(nlc),ifitstart,
     >     ibg,itry,stringindex(5,2),ii,if,ijump,gauss,ibest,ib2,idi2,
     >     ift,ist_grad,nend,irise,iseed,irise2,iset,burstid_id,iset2,irise3
      character*128 ifile, teststring, istring,istring2, source, catstring,burstid_src
      character*4 idi
      character*3 instrs,burstid_instr
      character*1 itype
      character*6 sstr
      common / LCdata / t1, pka, pkae, nodata
      logical jemxr
      common / fp / pownorm

      COMMON /rand/ iseed

      iset2=0
      itype = '-'

c initialize seed value randomizer
      open(29,file='iseed',status='old',err=909)
      read(29,*) iseed
      close(29)
      goto 910
 909  iseed = -10
 910  continue

c 
c      sh = 10.                  ! hardwired start time for fit if soft time is smaller
      ifile = 'anainput.txt'
      jemxr = .false.

      pcaarea = 1400.
      jemxarea = 64.
      nustararea = 1000.
      
c obtain fit settings and instrument
      open(21,file='instr.txt',status='old')
      read(21,*) instr, tmax, chidif, sh, hardstart
      close(21)
      instrs='---'
      if(instr.eq.1) instrs='WFC'
      if(instr.eq.2) instrs='PCA'
      if(instr.eq.3) instrs='JEM'
      if(instr.eq.4) instrs='XRT'
      if(instr.eq.5) instrs='NUS'

c obtain MJD
      open(21,file='mjd.txt',status='old')
      read(21,'(a)') istring
      close(21)
      if(instr.eq.2) then
         i = index(istring,'MJD')+3
         j = i+10
         read(istring(i:j),'(f11.5)') mjd
      else if(instr.eq.3) then
         i = index(istring,'MJD')+3
         j = i+11
         read(istring(i:j),'(f12.6)') mjd
      endif

c obtain source name
      do j = 1, 128
         source(j:j) = ' '
      enddo
      open(21,file='CAT',status='old')
      i = 0
 119  i = i + 1
      read(21,'(a)',end=120) catstring
      k = 0
      stringindex(1,1) = 1
      do j = 1, 64
         if(catstring(j:j).eq.',') then
            k = k + 1
            stringindex(k,2) = j-1
            if(k.lt.5) stringindex(k+1,1) = j+1
         endif
         if(k.eq.5) goto 121
      enddo
 121  stringindex(5,2) = index(catstring,'  ') - 1
      do k = 1, 5
         j = 0
         if(stringindex(k,1).lt.stringindex(k,2)) then
            j = index(istring,catstring(stringindex(k,1):stringindex(k,2)))
         endif
         if(j.gt.0) then
            source = catstring(stringindex(1,1):stringindex(1,2))
            goto 120
         endif
      enddo
      goto 119
 120  close(21)
      signal_cygnus_x1 = index(source,'Cyg')
      write(6,'("Source: ",a,1x,i2)') source, signal_cygnus_x1

c obtain input data
c first, for wfc
      if(instr.eq.1) then
         open(71,file=ifile,status='old')
         iread = 0
         i = 0
 1       read(71,'(a)',end=999) teststring
         i = i + 1
         backspace(71)
         if(teststring(1:5).eq.'No da') then
            read(71,102) teststring, t1(i), t2(i)
            pka (i) =  0.
            pkae(i) = -1.
            pkas(i) =  0.
            nodata(i) = 1
            if(i.gt.1) nodata(i-1) = 1
            if(i.lt.nlc-1) nodata(i+1) = 1
         else if(teststring(1:5).eq.'No at') then
            read(71,103) teststring, t1(i), t2(i)
            pka (i) =  0.
            pkae(i) = -1.
            pkas(i) =  0.
            nodata(i) = 0
         else
            read(71,*) j,k,t1(i),t2(i),ltf(i),x,y,pka(i),
     >           pkae(i),pkas(i)
            nodata(i) = 0
            iread = iread + 1
         endif
 102     format(a25,f12.6,1x,f12.6)
 103     format(a29,f12.6,1x,f12.6)
         goto 1
 999     write(6,*) 'Data points read and with data: ',i,iread
         n = i
         close(71)
         if(iread.eq.0) write(6,*) 'File has no valuable data'

         deltat = sngl((t1(2)-t1(1))*86400.d0)
         deltat=real(nint(deltat))
         mjd = t1(50)           ! JZZ
c     immediately reset time array to seconds
         if(n.gt.15) then
            t0 = t1(50)         ! JZZ
            do i = 1, n
               t1(i) = (t1(i)-t0)*86400.d0
               t2(i) = (t2(i)-t0)*86400.d0
            enddo
         else
            stop 'Time array too short'
         endif
      endif

c now is a good moment to get burstid_id, because only now do we know mjd for all 3 instruments
c obtain burst id
      burstid_id = -1
      open(21,file='burstid.txt',status='old')
      i = 0
 801  i = i + 1
      read(21,'(a)',end=804,err=804) istring2
      j = index(istring2,'" ') - 1
      k = j-1
      burstid_src(1:k) = istring2(2:j)
      do l = k+1, 128
         burstid_src(l:l) = ' '
      enddo
      j = j + 3
      do k = j, 128
         if(istring2(k:k).eq.' ') then
            goto 802
         endif
      enddo
 802  k = k - 1
      read(istring2(j:k),*) burstid_mjd
      j = k + 2
      do k = j, 128
         if(istring2(k:k).eq.' ') then
            goto 803
         endif
      enddo
 803  k = k - 1
      read(istring2(j:k),'(a)') burstid_instr
      j = k + 2
      read(istring2(j:),*) burstid_id
      if(index(burstid_src(1:20),source(1:10)).gt.0.and.86400.d0*dabs(burstid_mjd-mjd).lt.100.d0.and.
     >   ( burstid_instr(1:2).eq.'SW' .and. instr.eq.1 .or.
     >     burstid_instr(1:2).eq.'XP' .and. instr.eq.2 .or.
     >     burstid_instr(1:2).eq.'IJ' .and. instr.eq.3 ) ) then
         goto 805
      endif
      goto 801
 804  burstid_id = -1
 805  close(21)
      if(instr.eq.1) goto 791

c second, get input data for PCA
      if(instr.eq.2) then
         tri = 1.
         iset = 0
         open(24,file='tri.txt',err=932)
         i = 0
 933     i = i + 1
         read(24,'(i4,1x,f5.3)',err=934,end=934) idi2, itri
         if(burstid_id.eq.idi2) then
            tri = itri
            iset = 1
         endif
         goto 933
 934     close(24)
 932     continue
         if(iset.eq.1) write(6,*) 'KKK ',burstid_id,tri
         if(iset.eq.0) write(6,'("No tri found for ",i4)') burstid_id
         open(21,file='anainput.txt')
         i = 0
 21      i = i + 1
         read(21,*,end=299) number(i),t1(i),f1,f2
         pka(i)=sngl(f1)/tri
         pkae(i)=sngl(f2)/tri
         goto 21
 299     close(21)
         n = i-1
         if(n.gt.nlc) stop 'Too many data points..'
         deltat = sngl((t1(2)-t1(1)))
c third, for jemx
      else if(instr.eq.3) then
         open(21,file='anainput.txt')
         i = 0
 31      i = i + 1
         read(21,*,end=399) number(i),t1(i),f1,f2
         pka(i)=sngl(f1)
         pkae(i)=sngl(f2)
         goto 31
 399     close(21)
         n = i-1
         if(n.gt.nlc) stop 'Too many data points..'
         deltat = sngl((t1(2)-t1(1)))
c fourth, for nustar
      else if(instr.eq.4) then
         open(21,file='anainput.txt')
         i = 0
 41      i = i + 1
         read(21,*,end=499) number(i),t1(i),f1,f2
         pka(i)=sngl(f1)
         pkae(i)=sngl(f2)
         goto 41
 499     close(21)
         n = i-1
         if(n.gt.nlc) stop 'Too many data points..'
         deltat = sngl((t1(2)-t1(1)))
c fifth, for swift?
      else if(instr.eq.5) then
         open(21,file='anainput.txt')
         i = 0
 51      i = i + 1
         read(21,*,end=599) number(i),t1(i),f1,f2
         if(number(i).eq.61) t0 = t1(i)
         pka(i)=sngl(f1)
         pkae(i)=sngl(f2)
         goto 51
 599     close(21)
         n = i-1
         do i = 1, n
            t1(i) = t1(i) - t0
         enddo
         if(n.gt.nlc) stop 'Too many data points..'
         deltat = sngl((t1(2)-t1(1)))
      else
         stop 'Invalid instrument spec'
      endif
 791  call ms_insort(pka,pkas_sort,n)

c     delimit array to tmax
      do i = 1, n
         if(t1(i).gt.tmax) then
            j = i-1
            goto 221
         endif
      enddo
      goto 222
 221  n = j
 222  continue


c     renormalize some jem-x light curves
      stderr = 0.
      j = 0
      do i = 1, n
         if(pkae(i).gt.0.) then
            j = j + 1
            stderr = stderr + pkae(i)
         endif
      enddo
      stderr = stderr / j
      if(instr.eq.3) then
         jemxr = .false.
         j = 0
         if(stderr.lt.1.) then
            do i = 1, n
               pka(i) = pka(i) * 100.
               pkae(i) = pkae(i) * 100.
            enddo
            jemxr = .true.
         endif
      endif

c     normalize to per cm2, like for wfc
      if(instr.eq.2) then
         do i = 1, n
            pka(i) = pka(i) / pcaarea
            pkae(i) = pkae(i) / pcaarea
         enddo
      else if(instr.eq.3) then
         do i = 1, n
            pka(i) = pka(i) / jemxarea
            pkae(i) = pkae(i) / jemxarea
         enddo
      else if(instr.eq.4) then
         do i = 1, n
            pka(i) = pka(i) / nustararea
            pkae(i) = pkae(i) / nustararea
         enddo
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DETERMINE PEAK
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     determine peak
      pk = -1.E16
      peaktime = 0.d0
      do i = 1, n
         if(pka(i).gt.pk) then 
            ipk = i
            peaktime = t1(ipk)
            pk = pka(i)
            pke = pkae(i)
         endif
      enddo
      call getpeak(pka,pkae,n,ipk5,pk5,pk5e)
      call getpeakscales(pka,pkae,n,ipkarray,pkarray,pkearray)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DETERMINE BACKGROUND AND SUBTRACT FROM DATA AND PEAK VALUES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     subtract background
      bg = 0.
      ibg = 0
      rms = 0.
      s1 = stderr
      if(instr.eq.2) stderr = 0.
      do i = 1, n
         if(t1(i).lt.-15.d0.and.t1(i).gt.-50.) then
            ibg = ibg + 1
            bg = bg + pka(i)
            rms = rms + pka(i)**2
            if(instr.eq.2) stderr = stderr + pkae(i)
         endif
      enddo
      if(ibg.gt.1) then
         bg = bg / ibg
         rms = sqrt((rms/ibg - bg**2))
         stderr = stderr / ibg
         bge = stderr
      else
         bg = 0.
         rms = 0.
         stderr = s1
         bge = -1.
      endif
      if(ibg.le.0) then
         open(31,file='fit.error.txt',status='unknown',access='append')
         write(31,310) mjd, instr, 'no pre-burst data'
         close(31)
         itype = 'p'
      endif

c subtract background from data arrays and peak values
      do i = 1, n
         pka(i) = pka(i) - bg
      enddo
      pk = pk - bg
      do i = 1, 5, 2
         pkarray(i) = pkarray(i) - bg
      enddo
      
c     mend peak time if 5-sec average is very different (could result if data is very noisy)
      if(dabs(t1(ipk)-t1(ipk5)).gt.6.d0.and.abs(pk5/pk).gt.0.5) then ! was 0.4 JJJZZZ
         pk = -1.E16
         do i = ipk5-2, ipk5+2
            if(pka(i).gt.pk) then
               ipk = i
               pk = pka(i)
            endif
         enddo
      endif

      goto 881 ! skip this for now JJJZZZ
ccc the following is SKIPPED!
c     redefine peak if a better value comes out of binned result
      if(pk/pke.gt.5.) goto 881 ! only do this if first value is bad
      diff = abs(pkarray(3)-pk)
      difs = sqrt(pkearray(3)**2+pke**2)
      if(diff/difs.lt.2.) then
         ipk = ipkarray(3)
         pk = pkarray(3)
         pke = pkearray(3)
      endif
      diff = abs(pkarray(5)-pk)
      difs = sqrt(pkearray(5)**2+pke**2)
      if(diff/difs.lt.2.) then
         ipk = ipkarray(5)
         pk = pkarray(5)
         pke = pkearray(5)
      endif
 881  continue

      peaktime = t1(ipk)

c     define threshold for duration calculation later on
      threshold = 0.05*pk
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DETERMINE BURST START TIME
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     determine burst start time as first data point where flux is negative
      ist = 1
      starttime = t1(1)
      itry = 0
      do i = ipk, 1, -1
c         if(itry.eq.0.and.(pka(i)-bg).lt.0.0) then ! JJJZZZ  lt pkae(i)) then
         if(instr.ne.2.and.(itry.eq.0.and.(pka(i)-bg).lt.0.00*pk ) .or.
     >      instr.eq.2.and.(itry.eq.0.and.(pka(i)-bg).lt.0.05*pk)) then ! JJJZZZ  lt pkae(i)) then
            ist = i
            itry = 1
         else if(itry.eq.1) then
            if((ist-i)*deltat.lt.15.) then ! do not allow to have precursor 15 s from main burst
               if(pkae(i).gt.0.) then
c                  if((pka(i)-bg)/pkae(i).gt.3.5) then ! JJJZZZ
                   if((pka(i)-bg)/pkae(i).gt.3.5 .and.
     >                (pka(i)-bg).gt.0.1*pk) then ! JJJZZZ
                     itry = 0
                  endif
               endif
            endif
         endif
      enddo
c     5     continue
      if(ist.gt.1) then         ! JNEW
         if(abs(t1(ist)).lt.20.) then
            if(instr.ne.2) ist = ist + 1
            starttime = t1(ist)+0.5d0 * deltat
         else
            do i = 1, n
               if(dabs(t1(i)).lt.xmin) then
                  xmin = dabs(t1(i))
                  ist = i
               endif
            enddo
            starttime = 0.d0
         endif
      endif
c     alternative method to determine starttime - search for point with
c     largest gradient
      starttime_grad = starttime
      ist_grad = ist
      gradm = -1.e16
      do i = ipk-1, max(1,ipk-20), -1
         grad  = (pka(i+1)-pka(i))
         grade = sqrt( pkae(i+1)**2 + pkae(i)**2 )
         if(grad.gt.gradm) then
            gradm = grad
            starttime_grad = t1(i)
            ist_grad = i
         endif
      enddo
      if(instr.eq.21) then ! JZZ SEP 2018 - this was for RXTE, now use old algortihm also for XTE
            starttime = starttime_grad
            d = 1.e16
            ist = 1
            do i = 1, n
               if(abs(t1(i)-starttime).lt.d) then
                  d = abs(t1(i)-starttime)
                  ist = i
               endif
            enddo
c        endif
      endif
      
      if(hardstart.gt.-1.e3) then
         d = 1.e16
         ist = 1
         do i = 1, n
            if(abs(t1(i)-(hardstart+starttime)).lt.d) then
               d = abs(t1(i)-(hardstart+starttime))
               ist = i
            endif
         enddo
         starttime = hardstart+starttime
      endif

      do i = 1, n
         tss(i) = t1(i)-starttime
      enddo

      write(6,*) 'Starttime,peaktime=',Starttime,peaktime,ist

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INVESTIGATE SLOW RISE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
c     check whether there is some rise before start time
      irise = 0
      j = ist
 955  j = j - 1
      if(j.eq.1) then
         goto 956
      endif
      if((pka(j)-bg)/pkae(j).gt.5.) then
         irise = irise + 1
         goto 955
      else
         goto 956
      endif
 956  continue
      
c     check whether there is some rise before start time
      irise2 = 0
      j = ist
 855  j = j - 1
      if(j.eq.1) then
         goto 856
      endif
      if(pka(j).gt.threshold) then
         irise2 = irise2 + 1
         goto 855
      else
         goto 856
      endif
 856  continue

c     check whether there is some rise before peak time ! JZ Aug 16
      irise3 = 0
      j = ipk
 755  j = j - 1
      if(j.eq.1) then
         goto 756
      endif
      if(pka(j).gt.threshold) then
         irise3 = irise3 + 1
         goto 755
      else
         goto 756
      endif
 756  continue

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DETERMINE START POINT FOR FITS OF DECAY PORTIONS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      ifitstart = ipk
      if(pkarray(1)/pkearray(1).gt.6.) then
         do i = n, ipk+1, -1 ! JZJZJZ Mar16 changed from ipk to ipk+1
            if(pka(i).gt.0.50*pk) then
               if(pka(i-1)-pka(i).lt.-5.*pkae(i)) then
                  nend = i
                  goto 920
               endif
            endif
         enddo
         nend = ipk+2
 920     continue               
         do i = n, nend, -1
            if(pka(i).gt.0.75*pk.and.tss(i).lt.tmax) then
               ifitstart = i    !+ 1
               goto 921
            endif
         enddo
         ifitstart = nend
 921     continue
      else
         ifitstart = ipk+1
      endif
      if(pke.gt.0..and.(pk/pke).lt.6.) ifitstart = ipk +1
      if(signal_cygnus_x1.gt.0) ifitstart=ipk
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     SET NORMALIZATION POINT FOR POWER LAW FITS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set pownorm
      pownorm = tss(ifitstart)
      if(pownorm.le.0.) pownorm = 1.


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DO FITS
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c set some starting values for fit parameters
      i = n
      d1 = 10.*deltat
      if(tmax.gt.200.) d1 = tmax / 10.
c      if(source(1:5).eq.'GX 17') then
c         d1 = 300.
c         d2 = 1.
c      endif
      d2 = 1000.                ! skip data beyond 100 s to bring all data under same denomer
      if(instr.gt.3) d2=100.    ! to accommodate GX17+2
      if(i.gt.ifitstart+1) then ! JZJZ
 778     call getLCparams (instr,tmax,chidif,t1,pka,pkae,starttime,
     >        ifitstart,ipk,i,fitpeak1,fitpeake1,decaytime,
     >        decaytimee,bg1,bge1,deltat,d1,d2,chi,wght,sh)
 22      call getLCparams4(instr,tmax,chidif,t1,pka,pkae,starttime,
     >        ifitstart,ipk,i,fitpeak2,fitpeake2,plindex2,
     >        plindexe2,bg2,bge2,deltat,d1,d2,chi2,wght,sh,iset2)
         if(iset2.gt.0) then
            iset2 = -iset2
            goto 778
         endif
         fitpeak3 = fitpeak2    ! start values for fit with added gauss
         plindex3 = plindex2
         bg3 = bg2
         call getLCparams7(instr,tmax,chidif,decaytime,t1,pka,pkae,starttime,
     >        ifitstart,ipk,i,fitpeak3,fitpeake3,plindex3,
     >        plindexe3,bg3,bge3,rp,rpe,deltat,d1,d2,chi3,wght,sh)
         if(chi3.lt.chi2-1.) then
            tau = pownorm*(exp(1.)**(-1/plindex3)-1.)
            taue= tau-pownorm*(exp(1.)**(-1/(plindex3-plindexe3))-1.)
         else
            tau = pownorm*(exp(1.)**(-1/plindex2)-1.)
            taue= tau-pownorm*(exp(1.)**(-1/(plindex2-plindexe2))-1.)
         endif
      else
         fitpeak1 = 0.
         fitpeake1 = -1.
         fitpeak2 = 0.
         fitpeake2 = -1.
         fitpeak3 = 0.
         fitpeake3 = -1.
         decaytime = 0.
         decaytimee = -1.
         plindex2 = 0.
         plindexe2 = -1.
         plindex3 = 0.
         plindexe3 = -1.
         bg1 = 0.
         bge1 = -1.
         bg2 = 0.
         bge2 = -1.
         bg3 = 0.
         bge3 = -1.
      endif

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DECIDE WHICH FUNCTION IS BEST
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ibest = -1
      if(chi.lt.chi2.and.chi.lt.chi3) ibest = 1
      if(chi2.lt.chi.and.chi2.lt.chi3) ibest = 2
      if(chi3.gt.0..and.chi3.lt.chi.and.chi3.lt.chi2) ibest = 3
      if(chi.lt.1.1) then
         ibest = 1
      else if(chi2.lt.1.1) then
         ibest = 2
      else if(ibest.eq.3.and.chi3.gt.0.95*chi2) then
         ibest = 2
      else if(ibest.eq.3.and.chi2.lt.1.5.and.chi3.gt.0.85*chi2) then
         ibest = 2
      endif
      if(ibest.eq.2.and.plindex2+plindexe2.gt.-1.0.and.
     >   plindex3+plindexe3.lt.-1.0.and.abs(chi3-chi2)/chi2.lt.0.1) then
         ibest = 3
      endif
      if(signal_cygnus_x1.gt.0) ibest=1
      if(ibest.le.0) then
         if(chi2.le.chi ) ibest = 2
         if(chi .lt.chi2) ibest = 1
      endif
      if(ibest.le.0) then
         if(chi3.le.chi ) ibest = 3
         if(chi .lt.chi3) ibest = 1
      endif
      if(ibest.le.0) then
         if(chi3.le.chi2) ibest = 2
         if(chi2.lt.chi3) ibest = 3
      endif
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DETERMINE FLUENCE NUMBERS FOR 3 FITS AND DIRECTLY FROM DATA
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     determine straightforward fluence
      ft = 0.
      fte= 0.
      do i = ist-irise, n
         ft = ft + pka(i)
         fte= fte+ pkae(i)**2
      enddo
      fte = sqrt(fte)

c     determine model fluences for 3 model components
      totexp = fitpeak1 * exp(-tss(ifitstart)/decaytime)*decaytime
      totpl2 = -fitpeak2*pownorm/(plindex2+1)
      totple2= sqrt( (fitpeake2/fitpeak2)**2 + (plindexe2/(plindex2+1))**2 )
      totple2= totple2*totpl2
      if(chi3.lt.50..and.chi3.gt.0.01) then ! maart2016
         totpl3 = fitpeak3*pownorm/(-plindex3-1)
         totple3= sqrt( (fitpeake3/fitpeak3)**2 + (plindexe3/(-plindex3-1))**2 )
         totple3= totple3*totpl3
         if(pownorm.lt.rp(2)) then
            totgs = rp(1)*(1.+agauss(pownorm,rp(2),rp(3)))/2.
            totgse= rpe(1)*(1.+agauss(pownorm,rp(2),rp(3)))/2.
         else
            totgs = rp(1)*(1.-agauss(pownorm,rp(2),rp(3)))/2.
            totgse= rpe(1)*(1.-agauss(pownorm,rp(2),rp(3)))/2.
         endif
      else
         totpl3 = 0.
         totple3 = -1.
         totgs = 0.
         totgse = -1.
      endif

c     determine fluence of non-fitted part
      fluence = 0.
      fluencee = 0.
      do i = ist-irise, ifitstart
         fluence  = fluence  + pka(i) - bg ! JJJZZZ
         fluencee = fluencee + pkae(i)**2
      enddo

c     determine total model fluences
      ift = 0
      fluencet1 = fluence + totexp
      if(decaytimee.lt.1.E-16) then
         fluencet1  = 0.
         fluencete1 = 0.
      else
         fluencete1 = sqrt(fluencee + totexp**2 * ((fitpeake1/fitpeak1)**2 + 
     >        (decaytimee/decaytime)**2))
      endif

      fluencet2  = fluence + totpl2
      fluencete2 = sqrt(fluencee + totple2**2)
      
      fluencet3  = fluence + totpl3 + totgs
      fluencete3 = sqrt(fluencee**2 + totple3**2 + totgse**2)

      if(totgse.gt.0.) then
         fracgs  = totgs / ft
         fracgse = fracgs * sqrt( (totgse/totgs)**2 + (fte/ft)**2 )
      else
         fracgs  = 0.
         fracgse = fracgs * sqrt( (fte/ft)**2 )
      endif

      if(instr.gt.1) then
         i = index(istring,'MJD')-2
      else
         i = index(istring,'.txt') - 1
      endif
      if(chi .gt.9999.) chi3 = -1.
      if(chi2.gt.9999.) chi3 = -1.
      if(chi3.gt.9999.) chi3 = -1.

c     determine duration with a very simple method
c     follow fitted functions untill they decay to 5%
c     calculate errors with a monte carlo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     DETERMINE BURST END TIMES, DEFINED AS TIME WHEN MODEL PREDICTS A
c     DECAY TO 5% OF PEAK VALUE
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      call tend1(fitpeak1,fitpeake1,decaytime,decaytimee,bg1,threshold,te1,te1e)
      call tend2(fitpeak2,fitpeake2,pownorm,plindex2,plindexe2,bg2,threshold,te2,te2e)
      call tend3(fitpeak3,fitpeake3,pownorm,plindex3,plindexe3,rp,rpe,bg3,threshold,te3,te3e)

c     set final duration and fluence to those for best best fit model unless
c     power law index is larger than -1 when the fluence is set to straight one
      duration  = -1.
      duratione = -1.
      fluencef  = ft
      fluencefe = fte
      ifl=0
      idu=0
      if(ibest.eq.1) then
         if(te1.gt.0.) then
            duration = irise2 + te1
            duratione = sqrt(1.+te1e**2)
            idu=1
         endif
         if(decaytime.gt.0.) then
            fluencef  = fluencet1
            fluencefe = fluencete1
            ifl=1
         endif
      else if(ibest.eq.2) then
         if(te2.gt.0.) then
            duration = irise2 + te2
            duratione = sqrt(1.+te2e**2)
            idu=2
         endif
         if(plindex2.lt.-1.) then
            fluencef  = fluencet2
            fluencefe = fluencete2
            ifl=2
         endif
      else if(ibest.eq.3) then
         if(te3.gt.0.) then
            duration = irise2 + te3
            duratione = sqrt(1.+te3e**2)
            idu=3
         endif
         if(plindex3.lt.-1.2) then
            fluencef  = fluencet3
            fluencefe = fluencete3
            ifl=3
         endif
      endif
      if(ibest.eq.2.and.plindex2.gt.-1.2.and.(abs(chi-chi2)/chi2.lt.0.1.or.abs(chi-chi2).lt.1.)) then
        fluencef  = fluencet1
        fluencefe = fluencete1
        ifl=1
      endif
      if(ibest.eq.3.and.plindex3.gt.-1.2.and.(abs(chi-chi3)/chi3.lt.0.1.or.abs(chi-chi3).lt.1.)) then
        fluencef  = fluencet1
        fluencefe = fluencete1
        ifl=1
      endif
      if(ibest.eq.3.and.plindex3.gt.-1.2.and.(abs(chi2-chi3)/chi3.lt.0.1.or.abs(chi2-chi3).lt.1.)) then
        fluencef  = fluencet2
        fluencefe = fluencete2
        ifl=2
      endif
      if(ibest.eq.3.and.te3.lt.0..and.te2.gt.0..and.(abs(chi2-chi3)/chi3.lt.0.1.or.abs(chi2-chi3).lt.1.)) then
        fluencef  = fluencet2
        fluencefe = fluencete2
        ifl=2
        duration  = irise2 + te2
        duratione = sqrt(1.+te2e**2)
        idu=2
      endif
      if(duration.lt.duratione.and.idu.gt.1) then
         duration = irise2 + te1
         duratione = sqrt(1.+te1e**2)
         idu=1
      endif
      if(fluencef.lt.fluencefe) then
         fluencef  = ft
         fluencefe = fte
         ifl=0
         if(fluencef.lt.fluencefe) then
            fluencef = fluencet1
            fluencefe = fluencete1
            ifl = 1
         endif
      endif
      if(duration.lt.0..or.(duration.gt.300.d0.and.duration.lt.2.*duratione
     >   .and.duration.gt.10.*te1)) then ! JZJZ
         duration  = irise2 + te1
         duratione = sqrt(1.+te1e**2)
         idu=1
      endif
      if(duration.gt.tss(n).or.duration.le.1..or.idu.eq.0) then ! JZJZ Mar16 was2*tss
         open(31,file='fit.error.txt',status='unknown',access='append')
         write(31,310) mjd, instr,'data stretch too short'
         close(31)
         itype = 'd'
      endif
      if(fitpeak1.lt.0.) itype = 'n'

 310  format(f11.5,2x,i1,': ',a)
      
      open(31,file='ana.tri.out',access='append',status='unknown')
      if(burstid_id.le.0) then
         write(31,*) 'No burstid_id for ',istring(1:64)
      else
         write(31,'(i4,2x,a15,2x,f12.5,2x,f5.3)') burstid_id,burstid_src(1:15),burstid_mjd,tri
      endif
      close(31)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      mjd = starttime/86400.d0 + mjd
      delay = t1(ifitstart) - starttime

      if(chi2.gt.0..and.chi3.gt.0..and.chi3.lt.chi2-1.) then
         pl = plindex3
         ple = plindexe3
         bgo = bg3
         bgoe = bge3
         chie = chi3
         ib2 = 3
      else
         pl = plindex2
         ple = plindexe2
         bgo = bg2
         bgoe = bge2
         chie = chi2
         ib2 = 2
      endif
      open(24,file='fit.results.readme.txt',status='unknown')
      write(24,*) 'Columns:'
      write(24,*) ' 1 - burst-id'
      write(24,*) ' 2 - Source name'
      write(24,*) ' 3 - Instrument (1=WFC,2=PCA,3=JEMX)'
      write(24,*) ' 4 - MJD start of burst'
      write(24,*) ' 5 - Exponential fit Peak flux'
      write(24,*) ' 6 - Exponential fit Peak flux 1-sigma uncertainty'
      write(24,*) ' 7 - Exponential fit e-folding decay time'
      write(24,*) ' 8 - Exponential fit e-folding decay time uncert.'
      write(24,*) ' 9 - Exponential fit background'
      write(24,*) '10 - Exponential fit background uncert.'
      write(24,*) '11 - Exponential fit fluence'
      write(24,*) '12 - Exponential fit fluence uncert.'
      write(24,*) '13 - Exponential fit Chi2-red'
      write(24,*) '14 - Power law fit peak flux'
      write(24,*) '15 - Power law fit peak flux uncert.'
      write(24,*) '16 - Power law fit decay index'
      write(24,*) '17 - Power law fit decay index uncert.'
      write(24,*) '18 - Power law fit background'
      write(24,*) '19 - Power law fit background uncert.'
      write(24,*) '20 - Power law fit fluence'
      write(24,*) '21 - Power law fit fluence uncert.'
      write(24,*) '22 - Power law fit decay index Chi2-red'
      write(24,*) '23 - Power law + gaussian fit peak flux'
      write(24,*) '24 - Power law + gaussian fit peak flux uncert.'
      write(24,*) '25 - Power law + gaussian fit decay index'
      write(24,*) '26 - Power law + gaussian fit decay index uncert.'
      write(24,*) '27 - Power law + gaussian fit Gauss normalization'
      write(24,*) '28 - Power law + gaussian fit Gauss norm. uncrt.'
      write(24,*) '29 - Power law + gaussian fit Gauss centroid'
      write(24,*) '30 - Power law + gaussian fit Gauss centr. uncrt.'
      write(24,*) '31 - Power law + gaussian fit Gauss st. dev.'
      write(24,*) '32 - Power law + gaussian fit Gauss st. dev. uncrt.'
      write(24,*) '33 - Power law + gaussian fit background'
      write(24,*) '34 - Power law + gaussian fit background uncert.'
      write(24,*) '35 - Power law + gaussian fit fluence'
      write(24,*) '36 - Power law + gaussian fit fluence uncert.'
      write(24,*) '37 - Fraction gauss to modeled fluence'
      write(24,*) '38 - Fraction gauss to modeled fluence, uncrt.'
      write(24,*) '39 - Power law + gaussian fit chi2-red'
      write(24,*) '40 - Absolute peak flux'
      write(24,*) '41 - Absolute peak flux uncrt.'
      write(24,*) '42 - File containing data fitted'
      close(24)
      ii = index(istring,'.txt') + 3
      open(24,file='fit.results.txt',status='unknown',access='append')
      if(decaytimee.gt.0..and.chi2.gt.0.) decaytimee = decaytimee*sqrt(chi2) ! JZ Aug 16
      call correctmjd(instr,mjd,mjdc)
      write(24,1101)  burstid_id, source,instr,mjdc,itype,
     >     fitpeak1,fitpeake1,decaytime,decaytimee,bg1,bge1,fluencet1,fluencete1,chi,
     >     fitpeak2,fitpeake2,plindex2,plindexe2,bg2,bge2,fluencet2,fluencete2,chi2,
     >     fitpeak3,fitpeake3,plindex3,plindexe3,(rp(i),rpe(i),i=1,3),bg3,bge3,
     >     fluencet3,fluencete3,fracgs,fracgse,chi3,pk,pke,duration,duratione,pownorm,
     >     istring(1:ii)
 1101 format(i4,1x,a20,1x,i1,1x,f12.6,1x,a1, ! format mjd improved
     >     2(1x,f8.2),2(1x,f7.3),2(1x,f9.2),2(1x,f9.3),1x,f7.2,
     >     2(1x,f8.2),2(1x,f7.3),2(1x,f9.2),2(1x,f9.3),1x,f7.2,
     >     2(1x,f8.2),2(1x,f7.3),1x,3(1x,f7.2,1x,f7.2),2(1x,f9.2),2(1x,f9.3),
     >     2(1x,f6.3),1x,f7.2,2(1x,f8.2),3(1x,f7.2),1x,a)
      close(24)
c specifics about rise times and end times 
      open(24,file='fit.times.txt',status='unknown',access='append')
      write(24,1102)  burstid_id,source,instr,mjdc, irise3,irise2,te1,te1e,te2,te2e,
     >     te3,te3e,ibest,decaytime,decaytimee,duration,duratione,idu,istring(1:ii) ! JZ Aug 16
 1102 format(i4,1x,a20,1x,i1,1x,f12.6,2(1x,i3),6(1x,f8.2),1x,i2,4(1x,f8.2),1x,i1,1x,a)
      close(24)
c specifics about rise times and end times 
      open(24,file='fit.fluences.txt',status='unknown',access='append')
      write(24,1105) burstid_id,source,instr,mjdc,ft,fte,fluencet1,fluencete1,fluencet2,fluencete2,fluencet3,fluencete3,fluencef,fluencefe,ifl,istring(1:ii)
 1105 format(i4,1x,a20,1x,i1,1x,f12.6,10(1x,f8.2),1x,i1,1x,a)
      close(24)
c concise output, intended for catalog
      open(24,file='fit.results.concise.txt',status='unknown',
     >     access='append')
      write(24,1103) burstid_id,source,instrs,mjdc,itype,ibest,duration,duratione,decaytime,decaytimee,pk,pke,
     >     fluencef,fluencefe,irise3,istring(1:ii)
 1103 format(i4,1x,a20,1x,a3,1x,f12.6,1x,a1,1x,i2,8(1x,f10.3),1x,i3,1x,a) ! JZ Aug 16
      close(24)
c message when ifitstart was changed in power law fit
      if(iset2.ne.0) then
         open(24,file='fit.changestart.txt',status='unknown',
     >        access='append')
         write(24,1113) burstid_id,source,instrs,mjdc,-iset2,istring(1:ii)
 1113    format(i4,1x,a20,1x,a3,1x,f12.6,1x,i2,1x,a)
         close(24)
      endif
c html output, intended for catalog
      open(24,file='index0.html',status='unknown',
     >     access='append')
      write(24,1104) burstid_id,source,instrs,mjdc,itype,ibest,duration,duratione,decaytime,decaytimee,pk,pke,
     >     fluencef,fluencefe,irise3,istring(1:ii)
 1104 format('<th> ',i4,1x,'<th> ',a20,1x,'<th> ',a3,1x,'<th> ',f12.6,1x,'<th> ',a1,1x,'<th> ',i2,8(1x,'<th> ',f10.3),1x,'<th> ',i3,1x,'<th> <a href="',a,'_fit.pdf" target="plot tab">plot</a><tr>')
      close(24)
c screen output      
      write(6,110) source,instr,mjdc,pownorm,pl,ple,tau,taue,
     >     decaytime,decaytimee,pk,pke,fluencet1, fluencete1, ift,
     >     bg1,bge1,(starttime), ist, chie, ib2, ibest, chi, istring
c gauss output
      open(28,file='ana.gauss.out',status='unknown',access='append')
      write(28,112) source,instr,mjdc,pownorm,(rp(i),rpe(i),i=1,3), istring(1:ii)
      close(28)
c fit parameters to put on the plot
      open(28,file='ana.plot.out',status='unknown')
      write(28,*) (rp(i),rpe(i),i=1,3),mjdc,fracgs,fracgse ! JJJZ was geformatteerd met 1122 maar dat ging soms bij GX17+2 fout
      write(28,*) plindex2, plindexe2, plindex3, plindexe3, chi2, chi3
      write(28,*) duration, duratione, pk, pke, fluencef, fluencefe, ibest, ifl,idu
      write(28,*) burstid_id
      close(28)
c     endif
      
 110  format('Results: ',a20,2x,i1,2x,f12.6,2x,f6.1,6(2x,f12.3),2(2x,f9.3),
     >     2(2x,f10.3),2x,i1,2(2x,f9.4),2x,f5.1,2x,i5,2x,f7.2,
     >     2(2x,i1),2x,f7.2,2x,a)
 112  format(a20,2x,i1,2x,f11.5,2x,f6.1,3(2x,f7.2,2x,f7.2),2x,a40)
 1122 format(3(2x,f10.2,2x,f10.2),2x,f11.5,2(2x,f6.3))

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     OUTPUT MODEL LIGHT CURVES, FOR LATER PLOTTING
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      open(23,file='tt.out',status='unknown')
      write(23,*) n
      do i = 1, n
         if = 0
         if(wght(i).gt.0.) if = 1
         if(sh.gt.-99.) then
            tt = t1(i)-sh
         else
            tt = t1(i)-t1(ifitstart)
         endif
         tt2 = t1(i)-starttime
         if(tt.gt.-10..and.if.eq.1) then ! was tt.gt.0
            ff1 = fitpeak1*exp(-tt2/decaytime)+bg1
         else
            ff1 = bg1
         endif
         ff2=0.
         if(plindex2.lt.0.) then
            if(tt.gt.-10..and.if.eq.1) then
               ff2 = fitpeak2 * (tt2/pownorm)**plindex2 + bg2 
            else
               ff2 = bg2
            endif
         endif
         pi = 4. * atan(1.)
         fac = 1./(sqrt(2.*pi)*rp(3))
         expon = exp(-0.5*((tt2-rp(2))/rp(3))**2)
         ff3=0.
         if(plindex3.lt.0.) then
            if(tt.gt.-10..and.if.eq.1) then
               ff3 = fitpeak3 * (tt2/pownorm)**plindex3 + bg3 
     >              + rp(1)*fac*expon
            else
               ff3 = bg3
            endif
         endif
         write(23,*) tss(i),pka(i),pkae(i),ff1,ff2,ff3,if
      enddo
      close(23)

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine getLCparams(instr,tmax,chidif,t1,flux,fluxe,starttime,ifitstart,ipk,n,fitpeak,fitpeake,
     >                       time,timee,bg,bge,deltat,d1,d2,chisqr,wght,shard)

      implicit none
      integer*4 i,j,k,n,ipk,ieval,itry,nmax,instr,ifitstart
      parameter (nmax=100000)
      real*8 t1(*),starttime
      real*4 flux(*), fluxe(*), time, fluence, tt(nmax),chisqr,
     >       yfit(nmax),a(3),sigmaa(3),flamda,chisqrold,timee,tmax,chidif,
     >       wght(*),fitpeak,fitpeake,bg,bge,deltat,d1,d2,d3,fexp,shard
      real*8   scratch(3,3)
      external fexp, fexpderiv

      if(n.gt.nmax) stop 'getLCparams arrays too short..'
      d3 = -5. ! skip data below d3 sigma from zero
      itry = 0
      ieval = 0
      flamda = 0.001
      chisqrold = 0.
      if(n-ifitstart.le.0) then
         do i = 1, 3
            a(i) = 0.
            sigmaa(i) = -1.
         enddo
         return
      endif
      a(1) = d1
 3    a(2) = -1.e16
      do i = 1, n
         if(flux(i).gt.a(2)) a(2) = flux(i)
      enddo
      if(a(2).lt.0.) a(2)=-a(2)
      a(3) = 0. ; ! JNEW
      a(2) = a(2) - a(3)
      do i = 1, 3
         sigmaa(i) = 1.
      enddo
      do i = 1, n
         tt(i) = (t1(i)-starttime)
      enddo
      do i = 1, n
         if(fluxe(i).lt.1.E-6) then
            wght(i) = 0.
         else
            wght(i) = 1./fluxe(i)**2
         endif
         if(tt(i).gt.tmax.or.tt(i).lt.-50.) wght(i)=0.
         if(shard.gt.-99.) then
            if(tt(i).ge.-20..and.tt(i).le.shard) wght(i)=0. ! LLL
         else
            if(i.lt.ifitstart.and.tt(i).ge.-20.d0) wght(i)=0. ! KLM
         endif
      enddo
 2    call ma_lsqfits(tt,flux,wght,1,n,3,a,sigmaa,flamda,yfit,chisqr,
     >     ieval,scratch,fexp,fexpderiv)
      if(ieval.gt.10000) goto 1
      if((a(1).lt.2.).and.itry.eq.0) then
         a(1) = d1*2 ! JJJ changed 22-1-2016 
         a(2) = d2/2
         a(3) = 0. ! JNEW
         itry = 1
         goto 3
      endif
      if(chisqr.gt.1.e10) goto 1
      if(ieval.ge.2) then
         if(abs((chisqr-chisqrold)/chisqr).le.chidif) goto 1
      endif
      chisqrold = chisqr
      if(flamda.lt.1.e-5) flamda = 0.01
      goto 2
 1    time=a(1)
      timee=sigmaa(1)
      fitpeak=a(2)
      fitpeake=sigmaa(2)
      bg=a(3)
      bge=sigmaa(3)
      if(time.gt.2.*(tt(n)-tt(1)) .and. timee.lt.time/10.) then ! JZZ nov 2017
         time = 0.
         timee = 0.
         fitpeak = 0.
         fitpeake = 0.
      endif

      write(6,*) 'fit exponential', a(1),a(2),a(3),chisqr
      return
      end

c ======================================================================

      subroutine getLCparams4(instr,tmax,chidif,t1,flux,fluxe,starttime,ifitstart,ipk,n,fitpeak,fitpeake,
     >                       time,timee,bg,bge,deltat,d1,d2,chisqr,wght,shard,iset2)

      implicit none
      integer*4 i,j,k,n,ipk,ieval,itry,nmax,instr,ii,ifitstart,ie,iset2
      parameter (nmax=100000)
      real*8 t1(*),starttime
      real*4 flux(*), fluxe(*), time, fluence, tt(nmax),chisqr,chidif,
     >       yfit(nmax),a(3),sigmaa(3),flamda,chisqrold,timee,fpowbg,tmax,
     >       wght(*),fitpeak,fitpeake,deltat,d1,d2,d3,bg,bge,test,shard,
     >       y1,y2,y3
      real*8   scratch(3,3)
      real*4   pownorm
      common / fp / pownorm
      external fpowbg, fpowbgderiv

      if(n.gt.nmax) stop 'getLCparams arrays too short..'
      d3 = -5. ! skip data below d3 sigma from zero
 777  itry = 0
      ieval = 0
      flamda = 0.001
      chisqrold = 0.
      if(n-ifitstart.le.0) then
         do i = 1, 3
            a(i) = 0.
            sigmaa(i) = -1.
         enddo
         return
      endif
      a(1) = -1.
      a(2) = 1.
      a(3) = 0.
 3    do i = 1, 3
         sigmaa(i) = 1.
      enddo
      do i = 1, n
         tt(i) = (t1(i)-starttime)
      enddo
      do i = 1, n
         if(fluxe(i).lt.1.E-6) then
            wght(i) = 0.
         else
            wght(i) = 1./fluxe(i)**2
         endif
         if(tt(i).gt.tmax.or.tt(i).lt.-50.) wght(i)=0.
         if(shard.gt.-99.) then
            if(tt(i).ge.-20..and.tt(i).le.shard) wght(i)=0. ! LLL
         else
            if(i.lt.ifitstart.and.tt(i).ge.-20.d0) wght(i)=0. ! KLM
         endif
      enddo
 2    call ma_lsqfits(tt,flux,wght,1,n,3,a,sigmaa,flamda,yfit,chisqr,
     >     ieval,scratch,fpowbg,fpowbgderiv)
      if(ieval.gt.10000) goto 1
      ie = 0
      if((a(1).gt.0.or.a(1).lt.-3.).and.itry.lt.100) then
         a(1) = -1.
         ie = 1
      endif
      if((a(2).lt.0.).and.itry.lt.100) then
         a(2) = 2.
         ie = 1
      endif
c      if(a(3).lt.-3.*sigmaa(3).and.itry.lt.100) then
c         a(3) = 0.
c         ie = 1
c      endif
      if(ie.eq.1) then
         itry = itry+ 1
         goto 3
      endif
      if(chisqr.gt.1.e10) goto 1
      if(ieval.ge.2) then
         if(abs((chisqr-chisqrold)/chisqr).le.chidif) goto 1
      endif
      chisqrold = chisqr
      if(flamda.lt.1.e-5) flamda = 0.01
      goto 2
 1    time=a(1)
      timee=sigmaa(1)
      fitpeak=a(2)
      fitpeake=sigmaa(2)
      bg=a(3)
      bge=sigmaa(3)
      if(time.gt.2.*(tt(n)-tt(1))) then
         time = 0.
         timee = 0.
         fitpeak = 0.
         fitpeake = 0.
         bg = 0.
         bge = 0.
      endif

      y1 = (flux(ifitstart  )-yfit(ifitstart  ))/fluxe(ifitstart  )
      y2 = (flux(ifitstart+1)-yfit(ifitstart+1))/fluxe(ifitstart+1)
      y3 = (flux(ifitstart+2)-yfit(ifitstart+2))/fluxe(ifitstart+2)
      if(y1.lt.-5..and.y2.gt.0..and.y3.gt.0..and.(iset2.eq.0.or.iset2.eq.1.or.iset2.eq.2)) then
         if(ifitstart.lt.n-3) ifitstart = ifitstart + 1
         iset2=iset2+1
         goto 777
      endif

      return
      end

c ======================================================================

      subroutine getLCparams7(instr,tmax,chidif,dt,t1,flux,fluxe,starttime,ifitstart,ipk,n,fitpeak,fitpeake,
     >                       time,timee,bg,bge,rp,rpe,deltat,d1,d2,chisqr,wght,shard)

      implicit none
      integer*4 i,j,k,n,ipk,ieval,itry,nmax,instr,ii,ifitstart,ie
      parameter (nmax=100000)
      real*8 t1(*),starttime
      real*4 flux(*), fluxe(*), time, fluence, tt(nmax),chisqr,chidif,
     >       yfit(nmax),a(5),sigmaa(5),flamda,chisqrold,timee,fpowbg,tmax,
     >       wght(*),fitpeak,fitpeake,deltat,d1,d2,d3,bg,bge,test,
     >       rp(*),rpe(*),tau,dt,shard
      real*8   scratch(5,5)
      real*4   pownorm
      common / fp / pownorm
      external fpowgs3, fpowgs3deriv

      if(n.gt.nmax) stop 'getLCparams arrays too short..'
      itry = 0
      ieval = 0
      flamda = 0.001
      chisqrold = 0.
      if(n-ifitstart.le.0) then
         do i = 1, 5
            a(i) = 0.
            sigmaa(i) = -1.
         enddo
         return
      endif
      a(1) = -1.6
      a(2) = fitpeak
      a(3) = bg
      a(4) = 10.
      a(5) = 50. ! dt/2.
 3    do i = 1, 5
         sigmaa(i) = 1.
      enddo
      do i = 1, n
         tt(i) = (t1(i)-starttime)
      enddo
      do i = 1, n
         if(fluxe(i).lt.1.E-6) then
            wght(i) = 0.
         else
            wght(i) = 1./fluxe(i)**2
         endif
         if(tt(i).gt.tmax.or.tt(i).lt.-50.) wght(i)=0.
         if(shard.gt.-99.) then
            if(tt(i).ge.-20..and.tt(i).le.shard) wght(i)=0. ! LLL
         else
            if(i.lt.ifitstart.and.tt(i).ge.-20.d0) wght(i)=0. ! KLM
         endif
      enddo
 2    call ma_lsqfits(tt,flux,wght,1,n,5,a,sigmaa,flamda,yfit,chisqr,
     >     ieval,scratch,fpowgs3,fpowgs3deriv)
      if(ieval.gt.10000) goto 1
      ie = 0
      if((a(1).gt.-0.5.or.a(1).lt.-6.5).and.itry.lt.1000) then !JZJ
         a(1) = -1.5
         a(2) = 2.
         ie = 1
      endif
      if((a(2).lt.0.).and.itry.lt.1000) then
         a(2) = 2.
         ie = 1
      endif
      if((a(4).lt.0.).and.itry.lt.1000) then
         a(4) = 10.
         ie = 1
      endif
      if((a(4).gt.200..and.a(5).gt.100.).and.itry.lt.1000) then
         a(4) = 2.
         a(5) = 100.
         ie = 1
      endif
      if((a(5).lt.0.).and.itry.lt.1000) then
         a(5) =-a(5)
         ie = 1
      endif
      if(a(5).gt.150..and.itry.lt.1000) then
         a(5) = 10.
         ie = 1
      endif
      if((a(5).lt.5..or.a(5).gt.100.).and.itry.lt.1000) then
         a(5) = 40.
         ie = 1
      endif
         itry = itry+ 1
      if(ie.eq.1) then
         if(flamda.lt.1.e-5) flamda = 0.01
         goto 3
      endif
      if(chisqr.gt.1.e10) goto 1
      if(ieval.ge.2) then
         if(abs((chisqr-chisqrold)/chisqr).le.chidif) goto 1
      endif
      chisqrold = chisqr
      if(flamda.lt.1.e-5) flamda = 0.01
      goto 2
 1    time=a(1)
      timee=sigmaa(1)
      fitpeak=a(2)
      fitpeake=sigmaa(2)
      bg=a(3)
      bge=sigmaa(3)
      rp(1) = sqrt(a(4)**2+1.)-1.
      rpe(1) = rp(1)*abs(sigmaa(4)/a(4))
      rp(2) = 0.
      rpe(2) = -1.
      rp(3) = sqrt(a(5)**2+1.)-1.
      rpe(3) = rp(3)*abs(sigmaa(5)/a(5))

      return
      end

cc==========

      real*4 function fexp(x,i,a,n)
      implicit none
      integer*4 i, n
      real*4 a(*),x(*),expon

      if(-x(i)/a(1).lt.172.) then
         expon = exp(-x(i)/a(1))
      else
         expon = 0.
      endif
      if(x(i).gt.0.) then
         fexp = a(2) * expon + a(3)
      else
         fexp = a(3)
      endif
      return
      end

      subroutine fexpderiv(x,i,a,n,deriv)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),deriv(*),expon

      if(-x(i)/a(1).lt.172.) then
         expon = exp(-x(i)/a(1))
      else
         expon = 0.
      endif
      
      if(x(i).gt.0.) then
         deriv(1) = a(2) * (x(i)/a(1)**2) * expon
         deriv(2) = expon
      else
         deriv(1) = 0.
         deriv(2) = 0.
      endif
      deriv(3) = 1.

      return
      end

cc==========

      real*4 function fpow(x,i,a,n)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),pownorm
      common / fp / pownorm

      if(x(i).gt.0.) then
         fpow = a(2) * (x(i)/pownorm)**a(1)
      else
         fpow = 0.
      endif

      return
      end

      subroutine fpowderiv(x,i,a,n,deriv)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),deriv(*),pownorm
      common / fp / pownorm
      
      if(x(i).gt.0.) then
         deriv(1) = alog(x(i)/pownorm) * ((x(i)/pownorm)**a(1)) * a(2)
         deriv(2) = (x(i)/pownorm)**a(1)
      else
         deriv(1) = 0.
         deriv(2) = 0.
      endif
      
      return
      end

cc==========

      real*4 function fpowbg(x,i,a,n)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),pownorm
      common / fp / pownorm

      if(x(i).gt.0.) then
         fpowbg = a(2) * (x(i)/pownorm)**a(1) + a(3)
      else
         fpowbg = a(3)
      endif

      return
      end

      subroutine fpowbgderiv(x,i,a,n,deriv)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),deriv(*),pownorm
      common / fp / pownorm
      
      if(x(i).gt.0.) then
         deriv(1) = alog(x(i)/pownorm) * ((x(i)/pownorm)**a(1)) * a(2)
         deriv(2) = (x(i)/pownorm)**a(1)
         deriv(3) = 1.
      else
         deriv(1) = 0.
         deriv(2) = 0.
         deriv(3) = 1.
      endif

      return
      end

cc==========

      real*4 function fpowgs(x,i,a,n)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),pownorm,expon,fac,pi
      common / fp / pownorm

      pi = 4. * atan(1.)
      fac = 1./(sqrt(2.*pi)*a(6))
      expon = exp(-0.5*((x(i)-a(5))/a(6))**2)
      if(x(i).gt.0.) then
         fpowgs = a(2) * (x(i)/pownorm)**a(1) + a(3) +
     >            a(4)*fac*expon
      else
         fpowgs = a(3)
      endif

      return
      end

      subroutine fpowgsderiv(x,i,a,n,deriv)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),deriv(*),pownorm,expon,fac,pi
      common / fp / pownorm
      
      pi = 4. * atan(1.)
      fac = 1./(sqrt(2.*pi)*a(6))
      expon = exp(-0.5*((x(i)-a(5))/a(6))**2)
      if(x(i).gt.0.) then
         deriv(1) = alog(x(i)/pownorm) * ((x(i)/pownorm)**a(1)) * a(2)
         deriv(2) = (x(i)/pownorm)**a(1)
         deriv(3) = 1.
         deriv(4) = fac*expon
         deriv(5) = a(4) * fac * ( (x(i)-a(5))/a(6)**2    ) * expon
         deriv(6) = a(4) * fac * ( (x(i)-a(5))**2/a(6)**3 ) * expon
     >             -a(4) * fac * expon / a(6)
      else
         deriv(1) = 0.
         deriv(2) = 0.
         deriv(3) = 1.
         deriv(4) = 0.
         deriv(5) = 0.
         deriv(6) = 0.
      endif
      return
      end

cc==========

      real*4 function fpowgs2(x,i,a,n)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),pownorm,expon,fac,pi
      common / fp / pownorm

      pi = 4. * atan(1.)
      fac = 1./(sqrt(2.*pi)*a(5))
      expon = exp(-0.5*(x(i)/a(5))**2)
      if(x(i).gt.0.) then
         fpowgs2 = a(2) * (x(i)/pownorm)**a(1) + a(3) +
     >            a(4)*fac*expon
      else
         fpowgs2 = a(3)
      endif

      return
      end

      subroutine fpowgs2deriv(x,i,a,n,deriv)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),deriv(*),pownorm,expon,fac,pi
      common / fp / pownorm
      
      pi = 4. * atan(1.)
      fac = 1./(sqrt(2.*pi)*a(5))
      expon = exp(-0.5*(x(i)/a(5))**2)
      if(x(i).gt.0.) then
         deriv(1) = alog(x(i)/pownorm) * ((x(i)/pownorm)**a(1)) * a(2)
         deriv(2) = (x(i)/pownorm)**a(1)
         deriv(3) = 1.
         deriv(4) = fac*expon
         deriv(5) = a(4) * fac * ( x(i)**2/a(5)**3 ) * expon
     >             -a(4) * fac * expon / a(5)
      else
         deriv(1) = 0.
         deriv(2) = 0.
         deriv(3) = 1.
         deriv(4) = 0.
         deriv(5) = 0.
         deriv(6) = 0.
      endif
      return
      end

cc==========

      real*4 function fpowgs3(x,i,a,n)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),pownorm,expon,fac,pi
      common / fp / pownorm

      pi = 4. * atan(1.)
      fac = 1./(sqrt(2.*pi)* (sqrt(a(5)**2+1.)-1.))
      expon = exp(-0.5*(x(i)/(sqrt(a(5)**2+1.)-1.))**2)

      if(x(i).gt.0.) then
         fpowgs3 = a(2) * (x(i)/pownorm)**a(1) + a(3) +
     >        (sqrt(a(4)**2+1)-1.)*fac*expon
      else
         fpowgs3 = a(3)
      endif

      return
      end

      subroutine fpowgs3deriv(x,i,a,n,deriv)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),deriv(*),pownorm,expon,fac,pi
      common / fp / pownorm
      
      pi = 4. * atan(1.)
      fac = 1./(sqrt(2.*pi)*(sqrt(a(5)**2+1.)-1.))
      expon = exp(-0.5*(x(i)/(sqrt(a(5)**2+1.)-1.))**2)
      if(x(i).gt.0.) then
         deriv(1) = alog(x(i)/pownorm) * ((x(i)/pownorm)**a(1)) * a(2)
         deriv(2) = (x(i)/pownorm)**a(1)
         deriv(3) = 1.
         deriv(4) = fac*expon*a(4)/sqrt(a(4)**2+1.)
         deriv(5) = (sqrt(a(4)**2+1.)-1.) * fac * expon * (-0.5*x(i)**2  ) * (-2.*a(5)/(sqrt(a(5)**2+1.)-1.)**3./sqrt(a(5)**2+1.))
     >             +(sqrt(a(4)**2+1.)-1.) *       expon * (1./sqrt(2.*pi)) * (   -a(5)/(sqrt(a(5)**2+1.)-1.)**2./sqrt(a(5)**2+1.))
      else
         deriv(1) = 0.
         deriv(2) = 0.
         deriv(3) = 1.
         deriv(4) = 0.
         deriv(5) = 0.
      endif
      return
      end

cc==========

      real*4 function fexp2(x,i,a,n)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),expon
	
      if(-x(i)/a(1).lt.172.) then
         expon = exp(-x(i)/a(1))
      else
         expon = 0.
      endif
      fexp2 = a(2) * expon
      return
      end

      subroutine fexpderiv2(x,i,a,n,deriv)
      implicit none
      integer*4 i,n
      real*4 a(*),x(*),deriv(*),expon

      if(-x(i)/a(1).lt.172.) then
         expon = exp(-x(i)/a(1))
      else
         expon = 0.
      endif
      
      deriv(1) = a(2) * (+x(i)/a(1)**2) * expon
      deriv(2) = expon
      
      return
      end

c
c revisioned subroutine ma_lsqfit to be able to choose functions
c functn and fderiv
c
      SUBROUTINE MA_LSQFITS(X,Y,WEIGHT,IXS,IXE,NTERMS,A,SIGMAA,
     &     FLAMDA,YFIT,CHISQR,IEVAL,ARRAY,FUNCTN,FDERIV)
C
      IMPLICIT NONE
      INTEGER*4 I,J,NTERMS,IXS,IXE,NFREE,NTERM,IEVAL,K,NN
      REAL*4 X(*),Y(*),WEIGHT(*),A(*),SIGMAA(*),YFIT(*),
     >       ALPHA(30,30),BETA(30),DERIV(30),B(30),FLAMDA,CHISQ1,FUNCTI,
     >       FUNCTN,FU_CHISQ,CHISQR

      REAL*8		ARRAY(NTERMS,*), DET
      INTEGER*4 IKS1(30),IKS2(30)
      EXTERNAL FUNCTN,FDERIV
      DATA NFREE/0/
C
C** INITIALISE
C
      IF(IEVAL .GT. 0)GOTO 30
      IF (NTERMS.GT.30) STOP 'TOO MANY VARIABLES : MAXIMUM IS 30'
      NFREE=IXE-IXS+1-NTERMS
      NN=0
      DO I = IXS, IXE
         IF(WEIGHT(I).gt.0.) NN=NN+1
      ENDDO
      NFREE=NN-NTERMS
      IF (NFREE) 10,10,30
   10 CHISQR=0.
      RETURN
C
C** EVALUATE ALPHA AND BETA MATRICES
C
   30 DO 31 J=1,NTERMS
        BETA(J)=0.
        DO 31 K=1,NTERMS
         ALPHA(J,K)=0.
   31 CONTINUE
      DO 41 I=IXS,IXE
        CALL FDERIV(X,I,A,NTERMS,DERIV)
        YFIT(I)=FUNCTN(X,I,A,NTERMS)
        FUNCTI=WEIGHT(I)*(Y(I)-YFIT(I))
        DO 40 J=1,NTERMS
          BETA(J)=BETA(J)+FUNCTI*DERIV(J)
          DO 40 K=1,J
            ALPHA(J,K)=ALPHA(J,K)+WEIGHT(I)*DERIV(J)*DERIV(K)
   40   CONTINUE
        DO J=1,3
           DO K=1,3
           ENDDO
        ENDDO
   41 CONTINUE
      DO 42 J=1,NTERMS
      DO 42 K=1,J
        ALPHA(K,J)=ALPHA(J,K)
   42 CONTINUE
C
C** CHECK DIAGONAL OF ALPHA MATRIX
C
      DO 50 J = 1, NTERMS
	IF(ALPHA(J,J).LT.1E-32)THEN
		CHISQR = 1E+32		! IF DIAGONAL = 0.
		RETURN
	ENDIF
   50 CONTINUE
C
C** EVALUATE CHI SQUARE AT STARTING POINT
C
      CHISQ1=FU_CHISQ(Y,YFIT,WEIGHT,IXS,IXE,NFREE)
C
C** INVERT MODIFIED CURVATURE MATRIX TO FIND NEW PARAMETERS
C
   60 DO 62 J=1,NTERMS
        DO 61 K=1,NTERMS
          ARRAY(J,K)=ALPHA(J,K)/SQRT(ALPHA(J,J))/SQRT(ALPHA(K,K))
   61   CONTINUE
        ARRAY(J,J)=1.+FLAMDA
   62 CONTINUE
      CALL MM_MATINV(ARRAY,NTERMS,DET,IKS1,IKS2)
      DO 63 J=1,NTERMS
        B(J)=A(J)
        DO 63 K=1, NTERMS
          B(J)=B(J)+BETA(K)*ARRAY(J,K)/SQRT(ALPHA(J,J))/SQRT(ALPHA(K,K))
   63 CONTINUE
C
C** IF CHI SQUARE INCREASED, INCREASE FLAMDA AND TRY AGAIN
C
      DO 70 I=IXS,IXE
        YFIT(I)=FUNCTN(X,I,B,NTERMS)
   70 CONTINUE
      CHISQR=FU_CHISQ(Y,YFIT,WEIGHT,IXS,IXE,NFREE)
      IEVAL=IEVAL+1
c      WRITE(*,*) 'IEVAL ',IEVAL,'  CHISQR ',CHISQR,chisq1
      IF(ieval.GT.1500) GOTO 100 !! NEWNEW JZ, oct 95
      IF (CHISQR.GT.1E+32) GOTO 100
      IF (CHISQ1-CHISQR) 71,100,100
   71 FLAMDA=10.*FLAMDA
      GOTO 60
C
C** EVALUATE PARAMETERS AND UNCERTAINTIES
C
  100 DO 101 J=1,NTERMS
        A(J)=B(J)
        SIGMAA(J)=SQRT(ABS(ARRAY(J,J)/ALPHA(J,J)))
  101 CONTINUE
      FLAMDA=FLAMDA/10.
      RETURN
      END
C======================================================================
C  MM_MATINV.F77		LROLIB			  1982-OCT-26
C======================================================================

      SUBROUTINE MM_MATINV(ARRAY,NORDER,DET,IK,JK)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION ARRAY(NORDER,NORDER),IK(NORDER),JK(NORDER)
C
      DET=1.
      DO 100 K=1,NORDER
C
C** FIND LARGEST ELEMENT ARRAY(I,J) IN REST OF MATRIX
        RMAX=0.
   21   DO 30 I=K,NORDER
        DO 30 J=K,NORDER
          IF (ABS(RMAX)-ABS(ARRAY(I,J))) 24,24,30
   24     RMAX=ARRAY(I,J)
          IK(K)=I
          JK(K)=J
   30   CONTINUE
C
C** INTERCHANGE ROWS AND COLUMNS TO PUT RMAX IN ARRAY(K,K)
   31   IF (RMAX) 41,32,41
   32   DET=0.
        GOTO 140
   41   I=IK(K)
        IF (I-K) 21,51,43
   43   DO 50 J=1, NORDER
          SAVE=ARRAY(K,J)
          ARRAY(K,J)=ARRAY(I,J)
          ARRAY(I,J)=-SAVE
   50   CONTINUE
   51   J=JK(K)
        IF (J-K) 21,61,53
   53   DO 60 I=1,NORDER
          SAVE=ARRAY(I,K)
          ARRAY(I,K)=ARRAY(I,J)
          ARRAY(I,J)=-SAVE
   60   CONTINUE
C
C** ACCUMULATE ELEMENTS OF INVERSE MATRIX
   61   DO 70 I=1, NORDER
          IF (I-K) 63,70,63
   63     ARRAY(I,K)=-ARRAY(I,K)/RMAX
   70   CONTINUE
   71   DO 80 I=1,NORDER
        DO 80 J=1,NORDER
          IF (I-K) 74,80,74
   74     IF (J-K) 75,80,75
   75     ARRAY(I,J)=ARRAY(I,J)+ARRAY(I,K)*ARRAY(K,J)
   80   CONTINUE
   81   DO 90 J=1,NORDER
          IF (J-K) 83,90,83
   83     ARRAY(K,J)=ARRAY(K,J)/RMAX
   90   CONTINUE
        ARRAY(K,K)=1./RMAX
        DET=DET*RMAX
  100 CONTINUE
C
C** RESTORE ORDERING OF MATRIX
  101 DO 130 L=1,NORDER
        K=NORDER-L+1
        J=IK(K)
        IF (J-K) 111,111,105
  105   DO 110 I=1,NORDER
          SAVE=ARRAY(I,K)
          ARRAY(I,K)=-ARRAY(I,J)
          ARRAY(I,J)=SAVE
  110   CONTINUE
  111   I=JK(K)
        IF (I-K) 130,130,113
  113   DO 120 J=1,NORDER
          SAVE=ARRAY(K,J)
          ARRAY(K,J)=-ARRAY(I,J)
          ARRAY(I,J)=SAVE
  120   CONTINUE
  130 CONTINUE
  140 RETURN
      END

C======================================================================
C  FU_CHISQ.F77			LROLIB			  1982-OCT-29
C======================================================================

      FUNCTION FU_CHISQ ( Y, YFIT, WEIGHT, IXS, IXE, NFREE)
      DIMENSION Y(IXE),YFIT(IXE),WEIGHT(IXE)
C
      CHISQ=0.
      DO 10 I=IXS,IXE
        IF (YFIT(I).GT.1E+36) GOTO 30
        CHISQ=CHISQ+WEIGHT(I)*(Y(I)-YFIT(I))**2.
        IF (CHISQ.GT.1E+32) GOTO 30
   10 CONTINUE
      FU_CHISQ=CHISQ / NFREE
      RETURN
   30 FU_CHISQ=1E+32
      RETURN
      END
c ======================================================================

      subroutine getpeak(pka,pkae,n,ipk,pk,pke)

      implicit none
      integer*4 i,j,k,n,ipk,nlc
      parameter (nlc=3000)
      real*4 pka(*), pkae(*), pkas(nlc),pk,pke,ps,pss
      
       pk = -1.E16
       if(n.gt.nlc) stop 'problem nlc'
       do i = 3, n - 2
          ps = 0.
          pss = 0.
          do j = -2, 2
             if(pkae(i+j).gt.1.e-8) then
                ps = ps  + pka(i+j)/pkae(i+j)**2
                pss= pss + 1.      /pkae(i+j)**2
             else
                goto 9
             endif
          enddo
          if(pss.gt.1.e-8) then
             ps = ps/pss
             pss = 1./sqrt(pss)
          else
             ps = 0.
             pss = 0.
          endif
          if(ps.gt.pk.and.pss.gt.1.e-8) then
             ipk = i
             pk  = ps
             pke = pss
          endif
 9     enddo

       return
       end

c ======================================================================

      subroutine getpeakscales(pka,pkae,n,ipk,pk,pke)

      implicit none
      integer*4 i,j,k,n,ipk(5),nlc,side,iscale
      parameter (nlc=3000)
      real*4 pka(*), pkae(*), pkas(nlc),pk(5),pke(5),ps,pss

      do iscale = 1, 5
         pk(iscale) = -1.E16
         pke(iscale) = -1.
         ipk(iscale) = -1
      enddo
      do iscale = 1, 5, 2
         side = (iscale - 1 ) / 2
         do i = 1+side, n-side
            ps = 0.
            pss = 0.
            do j = -side, side
               if(pkae(i+j).gt.1.e-8) then
                  ps = ps  + pka(i+j)/pkae(i+j)**2
                  pss= pss + 1.      /pkae(i+j)**2
               endif
            enddo
            if(pss.gt.1.e-8) then
               ps = ps/pss
               pss = 1./sqrt(pss)
            else
               ps = 0.
               pss = 0.
            endif
            if(ps.gt.pk(iscale).and.pss.gt.1.e-8) then
               ipk(iscale) = i
               pk(iscale)  = ps
               pke(iscale) = pss
            endif
         enddo
      enddo
      
      return
      end
      
c ======================================================================

C  MS_INSORT.F77		LROLIB			  1982-MAR-16
C======================================================================

      SUBROUTINE MS_INSORT(A, IX, JJ)
C*********************************************************************
C THIS SUBPROGRAM IS PORTABLE TO OTHER TYPES OF COMPUTER
C ALTHOUGH ORIGINALLY WRITTEN IN FORTRAN IV, THE PROGRAM IS
C ALSO USABLE WITH FORTRAN 77
C*********************************************************************
      DIMENSION A(JJ), IX(JJ)
      DIMENSION IU(16), IL(16)
	DO 5 I = 1,JJ
	IX(I) = I
 5	CONTINUE
      II = 1
      M = 1
      I = II
      J = JJ
10    IF (I.GE.J) GO TO 80
20    K = I
      IJ = (J+I)/2
      T = A(IX(IJ))
	T1 = IX(IJ)
      IF (A(IX(I)).LE.T) GO TO 30
      IX(IJ) = IX(I)
      IX(I) = T1
      T = A(IX(IJ))
	T1 = IX(IJ)
30    L = J
      IF (A(IX(J)).GE.T) GO TO 50
      IX(IJ) = IX(J)
      IX(J) = T1
      T = A(IX(IJ))
	T1 = IX(IJ)
      IF (A(IX(I)).LE.T) GO TO 50
      IX(IJ) = IX(I)
      IX(I) = T1
      T = A(IX(IJ))
	T1 = IX(IJ)
      GO TO 50
40    IX(L) = IX(K)
      IX(K) = TT1
50    L = L - 1
      IF (A(IX(L)).GT.T) GO TO 50
      TT = A(IX(L))
	TT1 = IX(L)
60    K = K + 1                                                         
      IF (A(IX(K)).LT.T) GO TO 60
      IF (K.LE.L) GO TO 40
      IF (L-I.LE.J-K) GO TO 70                                          
      IL(M) = I                                                         
      IU(M) = L                                                         
      I = K                                                             
      M = M + 1                                                         
      GO TO 90                                                          
70    IL(M) = K                                                         
      IU(M) = J                                                         
      J = L                                                             
      M = M + 1                                                         
      GO TO 90                                                          
80    M = M - 1                                                         
      IF (M.EQ.0) RETURN                                                
      I = IL(M)                                                         
      J = IU(M)                                                         
90    IF (J-I.GE.11) GO TO 20                                           
      IF (I.EQ.II) GO TO 10                                             
      I = I - 1                                                         
100   I = I + 1                                                         
      IF (I.EQ.J) GO TO 80
      T = A(IX(I+1))
	T1 = IX(I+1)
      IF (A(IX(I)).LE.T) GO TO 100
      K = I                                                             
110   IX(K+1) = IX(K)
      K = K - 1
      IF (T.LT.A(IX(K))) GO TO 110
      IX(K+1) = T1
      GO TO 100
      END                                                               

c function agauss.f
c
c source
c   Bevington, page 48.
c
c purpose
c   evaluate integral of gaussian probability function
c
c usage
c   result = agauss (x, averag, sigma)
c
c description of parameters
c   x      - limit for integral
c   averag - mean of distribution
c   sigma  - standard deviation of distribution
c   integration range is averag +/- z*sigma
c      where z = abs(x-averag)/sigma
c
c subroutines and function subprograms required
c   none
c
	real*4 function agauss (x,averag,sigma)
	double precision z,y2,term,sum,denom
11	z=abs(x-averag)/sigma
        if(z.gt.10.) then
           agauss=1.
           return
        endif
	agauss=0.
	if (z) 42,42,21
21	term=0.7071067812*z
22	sum=term
	y2=(z**2)/2.
	denom=1.
c
c accumulate sums of terms
c
31	denom=denom+2.
32	term=term*(y2*2./denom)
33	sum=sum+term
	if (term/sum-1.e-10) 41,41,31
41	agauss=1.128379167*sum*dexp(-y2)
42	return
	end

c     ============================================================
      
      subroutine tend1(norm,norme,tau,taue,bg,threshold,te1,te1e)

      implicit none
      integer*4 i,j,k,n,iseed,try
      real*4 norm,norme,tau,taue,threshold,te1,te1e,normt,taut,tt(1000),gauss,bg

      COMMON /rand/ iseed

      if(norm.lt.0..or.tau.lt.0.1) then
         te1  = -1.
         te1e = -1.
         return
      endif
      
      try = 0
      te1  = 0.
      te1e = 0.
      k    = 0
 1    k = k + 1
      if(try.eq.10000) then
         te1  = 0.
         te1e = -1.
         return
      endif
      normt = gauss(norm,norme)
      if(normt.le.0.) then
         try = try + 1
         k = k - 1
         goto 1
      endif
      taut  = gauss(tau ,taue )
      if(taut.le.0.) then
         try = try + 1
         k = k - 1
         goto 1
      endif
      tt(k) = -taut*alog((threshold)/normt)
      te1   = te1 + tt(k)
      if(k.lt.1000) goto 1
      te1 = te1 / 1000.
      te1e = 0.
      do i = 1, 1000
         te1e = te1e + (tt(i)-te1)**2.
      enddo
      te1e = sqrt(te1e/(k-1))
      return
      end

c     ============================================================
      
      subroutine tend2(norm,norme,tnorm,pl,ple,bg,threshold,te1,te1e)

      implicit none
      integer*4 i,j,k,l,n,iseed,try
      real*4 norm,norme,normt,pl,ple,plt,threshold,te1,te1e,tt(1000),
     >       gauss,tnorm,ff,bg,rp(3),pltime

      COMMON /rand/ iseed

      if(threshold.lt.bg.or.norm.lt.0..or.pl.gt.0.) then
         te1  = -1.
         te1e = -1.
         return
      endif

      rp(1) = 0.
      rp(2) = 0.
      rp(3) = 0.
      
      te1  = 0.
      te1e = 0.
      k    = 0
      try = 0
      l    = 0
 1    k = k + 1
      if(try.eq.10000) then
         te1  = 0.
         te1e = -1.
         return
      endif
      normt = gauss(norm,norme)
      if(normt.le.0.) then
         try = try + 1
         k = k - 1
         goto 1
      endif
      plt   = gauss(pl ,ple )
      if(plt.gt.0.) then
         try = try + 1
         k = k - 1
         goto 1
      endif

c     determine tend empirically
      tt(k) = pltime(plt,normt,tnorm,bg,threshold,rp)

      if(tt(k).gt.0.) then
         l = l + 1
         te1   = te1 + tt(k)
      endif
      if(k.lt.1000) goto 1


      if(l.gt.1) then
         te1 = te1 / l
         te1e = 0.
         do i = 1, 1000
            if(tt(i).gt.0.) te1e = te1e + (tt(i)-te1)**2.
         enddo
         te1e = sqrt(te1e/(k-1))
      endif
      te1 = pltime(pl,norm,tnorm,bg,threshold,rp)

      return
      end

c     ============================================================
      
      subroutine tend3(norm,norme,tnorm,pl,ple,rp,rpe,bg,threshold,te1,te1e)

      implicit none
      integer*4 i,j,k,l,n,iseed,try
      real*4 norm,norme,normt,pl,ple,plt,threshold,te1,te1e,pltime,
     >       tt(1000),gauss,tnorm,rp(3),rpe(3),rpt(3),expon,fac,pi,ff,bg

      COMMON /rand/ iseed

      if(threshold.lt.bg.or.norm.lt.0..or.pl.gt.0.) then ! .or.rp(1).lt.0.) then JZJZ
         te1  = -1.
         te1e = -1.
         return
      endif

      
      pi = 4. * atan(1.)
      te1  = 0.
      te1e = 0.
      l = 0
      try =0
      k    = 0
 1    k = k + 1
      if(try.eq.10000) then
         te1  = 0.
         te1e = -1.
         return
      endif
      normt = gauss(norm,norme)
      plt   = gauss(pl ,ple )
      do i = 1, 3
         rpt(i) = gauss(rp(i),rpe(i))
      enddo
      if(normt.le.0..or.plt.gt.-0.1) then
         try = try + 1
         k = k - 1
         goto 1
      endif

c     determine tend empirically
      tt(k) = pltime(plt,normt,tnorm,bg,threshold,rpt)
      
      if(tt(k).gt.0.) then
         l = l + 1
         te1   = te1 + tt(k)
      endif
      if(k.lt.1000) goto 1

      if(l.gt.1) then
         te1 = te1 / l
         te1e = 0.
         do i = 1, 1000
            if(tt(i).gt.0.) te1e = te1e + (tt(i)-te1)**2.
         enddo
         te1e = sqrt(te1e/(l-1))
         te1 = pltime(pl,norm,tnorm,bg,threshold,rp)
      else
         te1  = -1.
         te1e = -1.
      endif
      
      return
      end
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      real*4 function pltime(pl,norm,tnorm,bg,threshold,rp)
c probes up to 10000 times tnorm      
      implicit none
      integer*4 i,j
      real*4   pl,norm,rp(3),fac,expon,tnorm,ff,threshold,bg,pi

      pi = 4. * atan(1.)

      pltime = -1.
      j = 0
 2    j = j + 1
      if(abs(rp(1)).gt.1E-16) then
         expon = exp(-0.5*(((1.+j/10.)*tnorm-rp(2))/rp(3))**2)
         fac = 1./(sqrt(2.*pi)*rp(3))
      else
         expon = 0.
         fac = 0.
      endif
      ff = norm * ((1.+j/10.))**pl
     >     + rp(1)*fac*expon + bg
      if(ff.lt.threshold) then
         pltime = (1.+j/10.)*tnorm
      else
         if(j.lt.100000) goto 2
      endif

      return
      end
      
      REAL*4 FUNCTION  gauss ( mn, sgm )
*
************************ C O M I S ******************SRU***Utrecht******
*                                      *
* Function   : GAUSS               Rev.: 1.00    may 1985  author: jjs *
*                                      *
* -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  *
*
* Real*4 function, returning one sample from a Normal distributed
* variable, with mean mn and standard deviation sgm.
*      Based on lrolib:FU_NORMAL by willem van dijk
*
*   literature: knuth, seminumerical algorithms, page 104
*
*                 - Parameters -
* I/O   name         type              description
*
*  O   normal         R*4   Function, returns sample.
*  I   mn, sgm        R*4   Normal distribution parameters.
*
       IMPLICIT none
       INTEGER*4 iseed 
       REAL*4    s, u1, u2, v1, v2, mn, sgm, ran2 
       COMMON /rand/ iseed 

       u2 = ran2 ( iseed ) 
1      CONTINUE                ! REPEAT
       u1 = u2 
       u2 = ran2 ( iseed ) 
       v1 = u1 + u1 - 1.0 
       v2 = u2 + u2 - 1.0 
       s  = v1 * v1 + v2 * v2 
       IF ( s .GT. 1. ) GOTO 1         ! UNTIL s <= 1
       gauss = mn + v1 * sgm * sqrt ( -2.0 * alog ( s ) / s ) 

       RETURN 
       END

      function ran2(idum) 
      integer idum, im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv 
      real    ran2,am,eps,rnmx 
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,ia1
     $     =40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,ntab
     $     =32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps) 
      integer idum2,j,k,iv(ntab),iy 
      save iv,iy,idum2 
      data idum2/123456789/, iv/ntab*0/, iy/0/ 
      if(idum.le.0) then 
         idum=max(-idum,1) 
         idum2=idum 
         do j=ntab+8,1,-1 
            k=idum/iq1 
            idum=ia1*(idum-k*iq1)-k*ir1 
            if(idum.lt.0) idum=idum+im1 
            if(j.le.ntab) iv(j)=idum 
         enddo 
         iy=iv(1) 
      endif 
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if(idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if(idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1) iy=iy+imm1
      ran2=min(am*iy,rnmx)
      return
      end

      subroutine correctmjd(i,mjdi,mjdo)

c incorporate varying leap seconds
      
      implicit none
      integer*4 i
      real*8 mjdi,mjdo,leapsecond

      mjdo=mjdi

      if(i.eq.2) mjdo=mjdi+(leapsecond(mjdi)+32.184d0-60.184d0-3.37D0)/86400.d0
      if(i.eq.3) mjdo=mjdi+(leapsecond(mjdi)+32.184d0-64.184d0)/86400.d0

      return
      end

      real*8 function leapsecond(mjd)

c leap second, to go from TAI to UTC

      implicit none
      real*8 mjd,jd

      jd=mjd+2400000.5

      if(jd.gt.2437300.5) leapsecond=1.4228180 
      if(jd.gt.2437512.5) leapsecond=  1.3728180 
      if(jd.gt.2437665.5) leapsecond=  1.8458580 
      if(jd.gt.2438334.5) leapsecond=  1.9458580 
      if(jd.gt.2438395.5) leapsecond=  3.2401300 
      if(jd.gt.2438486.5) leapsecond=  3.3401300 
      if(jd.gt.2438639.5) leapsecond=  3.4401300 
      if(jd.gt.2438761.5) leapsecond=  3.5401300 
      if(jd.gt.2438820.5) leapsecond=  3.6401300 
      if(jd.gt.2438942.5) leapsecond=  3.7401300 
      if(jd.gt.2439004.5) leapsecond=  3.8401300 
      if(jd.gt.2439126.5) leapsecond=  4.3131700 
      if(jd.gt.2439887.5) leapsecond=  4.2131700 
      if(jd.gt.2441317.5) leapsecond= 10.0       
      if(jd.gt.2441499.5) leapsecond= 11.0       
      if(jd.gt.2441683.5) leapsecond= 12.0       
      if(jd.gt.2442048.5) leapsecond= 13.0       
      if(jd.gt.2442413.5) leapsecond= 14.0       
      if(jd.gt.2442778.5) leapsecond= 15.0       
      if(jd.gt.2443144.5) leapsecond= 16.0       
      if(jd.gt.2443509.5) leapsecond= 17.0       
      if(jd.gt.2443874.5) leapsecond= 18.0       
      if(jd.gt.2444239.5) leapsecond= 19.0       
      if(jd.gt.2444786.5) leapsecond= 20.0       
      if(jd.gt.2445151.5) leapsecond= 21.0       
      if(jd.gt.2445516.5) leapsecond= 22.0       
      if(jd.gt.2446247.5) leapsecond= 23.0       
      if(jd.gt.2447161.5) leapsecond= 24.0       
      if(jd.gt.2447892.5) leapsecond= 25.0       
      if(jd.gt.2448257.5) leapsecond= 26.0       
      if(jd.gt.2448804.5) leapsecond= 27.0       
      if(jd.gt.2449169.5) leapsecond= 28.0       
      if(jd.gt.2449534.5) leapsecond= 29.0       
      if(jd.gt.2450083.5) leapsecond= 30.0       
      if(jd.gt.2450630.5) leapsecond= 31.0       
      if(jd.gt.2451179.5) leapsecond= 32.0       
      if(jd.gt.2453736.5) leapsecond= 33.0       
      if(jd.gt.2454832.5) leapsecond= 34.0       
      if(jd.gt.2456109.5) leapsecond= 35.0       
      if(jd.gt.2457204.5) leapsecond= 36.0       
      if(jd.gt.2457754.5) leapsecond= 37.0       

      return
      end
