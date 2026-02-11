      program ROP2009MA
c
c     Code to supply regional optical model potential (OMP) parameters for alpha
c     -particles at incident energies E<50 MeV, on target nuclei with 50<A<209,
c     (M. Avrigeanu et al., At. Data Nucl. Data Tables 95(2009)501; EFFDOC-1081,
c     OECD/NEA, Paris, Nov. 2009).
c     V. Avrigeanu, IFIN-HH Bucharest, P.O. Box MG6, E-mail: vavrig@nipne.ro
c     
c     Code compiled with: MS Fortran (WINDOWS)
c	g95 -c -fsloppy-char -ftrace=frame ROP2009MA.for
c	g95 -o ROP2009MA -fsloppy-char ROP2009MA.o
c     and run with:
c	echo ROP2009MA.in ROP2009MAza | ROP2009MA
c	
c	INPUT FILES
c	ROP2009MA.in = input data (see below)       -       FORMAT: character*12
c
c	OUTPUT FILE
c	ROP2009MAza  = plain OMP parameters for (z,a) target nucleus, in tabular  
c		       form as well as in RIPL-2 format  -  FORMAT: character*12
c
c	INPUT DATA (ROP2009MA.in file) besides IN/OUT file names on command line
c	[ALL INPUT PARAMETERS ARE READ IN  RIPL-2  FREE FORMAT READ STATEMENTS,
c	e.g. http://www-nds.iaea.org/RIPL-2/optical/om-data/om-parameter.readme]
c	------------------------------------------------
c	iref
c	author      [1 line of author names]
c	reference   [1 line of reference information]
c	summary     [4 lines of descriptive information]
c	emin,emax
c	izmin,izmax
c	iamin,iamax
c	imodel,izproj,iaproj,irel,idr
c	(1) *****	
c	iztarget,iatarget	 [the only different by above RIPL-2 parameters]
c	    ***** STOP ***** IF END OF INPUT DATA FILE
c	    ***** SKIP TO (1)
c	------------------------------------------------
c
c     om-utility.cmb common block
c     Common blocks and declarations for om-summary.f and
c     om-modify.f codes (1/16/2002)
c
c     Parameter statement
c
      parameter (ndim1=10, ndim2=13, ndim3=25, ndim4=30, ndim5=10,
     1 ndim6=10,ndim7=120)
c
c     Common block and declarations for RIPL optical parameters
c
      character*1 author,refer,summary
      common/lib/author(80),refer(80),summary(320),iref,emin,emax,
     + izmin,izmax,iamin,iamax,imodel,jrange(6),epot(6,ndim1),
     + rco(6,ndim1,ndim2),aco(6,ndim1,ndim2),pot(6,ndim1,ndim3),
     + ncoll(ndim4),nvib(ndim4),nisotop,iz(ndim4),ia(ndim4),
     + lmax(ndim4),bandk(ndim4),def(ndim4,ndim5),idef(ndim4),
     + izproj,iaproj,exv(ndim7,ndim4),iparv(ndim7,ndim4),irel,idr,
     + nph(ndim7,ndim4),defv(ndim7,ndim4),thetm(ndim7,ndim4),
     + beta0(ndim4),gamma0(ndim4),xmubeta(ndim4),ex(ndim6,ndim4),
     + spin(ndim6,ndim4),ipar(ndim6,ndim4),spinv(ndim7,ndim4),
     + jcoul,ecoul(ndim1),rcoul(ndim1),rcoul0(ndim1),beta(ndim1),
     + rcoul1(ndim1),rcoul2(ndim1)
      common/io/ki,ko,ieof
c
      character*12 filein,fileou1,fileou2
C     read(*,'(A)')filein
C     read(*,'(A)')fileout
      filein ='ROP2009MA.in'
      fileou1='ROP2009MA.ou'
      fileou2='ROP2009MA.ri'
	ki=5
	ko=8			
      OPEN(ki,FILE=filein,STATUS='OLD')			
      OPEN(6,FILE=fileou1,STATUS='UNKNOWN')
      OPEN(ko,FILE=fileou2,STATUS='UNKNOWN')

      read(ki,*,end=999) iref0
	iref = iref0-1
      read(ki,9) (author(m),m=1,80)
      read(ki,9) (refer(m),m=1,80)
      read(ki,9) (summary(m),m=1,320)
      read(ki,*) emin,emax
      read(ki,*) izmin,izmax
      read(ki,*) iamin,iamax
      read(ki,*) imodel,izproj,iaproj,irel,idr
	z1=izproj	
	a1=iaproj	
  100 read(5,*,end=999)iztarget,iatarget
	z2=iztarget
	a2=iatarget
c
	iref = iref +1
c
      do 110 i=1,6
	jrange(i)=0
	ecoul(i)=0.
	rcoul0(i)=0.
	rcoul(i)=0.
	rcoul1(i)=0.
	rcoul2(i)=0.
	beta(i)=0.
	do 110 j=1,6
	epot(i,j)=0.
	do 110 k=1,11
	rco(i,j,k)=0.
	aco(i,j,k)=0.
	pot(i,j,k)=0.
	do 110 kV=12,25
	pot(i,j,kV)=0.
  110 continue
c
	a213=a2**(1./3.)
	a13=a1**(1./3.)+a2**(1./3.)
	RB=1.36*a13+0.5				! W.N.'80, p.6,8
	VClb=1.44*z1*z2/RB
	E2=0.9*VClb*(a1+a2)/a2			!E2=(2.59+10.4/a2)*z2/RB
	if(a2.le.130.)then
	  E1=-3.12-0.762*a2**(1./3.)+1.24*E2
	  else
	  E1=-3.45-0.762*a2**(1./3.)+1.24*E2
	endif
	E3 = 23.6+0.181*z2/a2**(1./3.)
	E3rR=25.
	E3aR=29.1-0.22*z2/a2**(1./3.)
	E4 = 50.
c
	jrange(1)=5
	jrange(2)=1
	jrange(4)=3
	jcoul=1
c
	epot(1,1)=E2
	if(E3.lt.E3rR)then
	  epot(1,2)=E3
	  epot(1,3)=E3rR
	  epot(1,4)=E3aR
	  elseif(E3.lt.E3aR)then
	    epot(1,2)=E3rR
	    epot(1,3)=E3
	    epot(1,4)=E3aR
	    else
	    epot(1,2)=E3rR
	    epot(1,3)=E3aR
	    epot(1,4)=E3
	endif
	epot(1,5)=E4
	epot(2,1)=E4
	epot(4,1)=E1
	epot(4,2)=E2
	epot(4,3)=E4
	ecoul(1)=E4
C
	rco(1,1,1)=1.18
	rco(1,1,2)=0.012
	rco(1,2,1)=rco(1,1,1)
	rco(1,2,2)=rco(1,1,2)
	if(E3.lt.E3rR)then
	  rco(1,3,1)=rco(1,1,1)
	  rco(1,3,2)=rco(1,1,2)
	  rco(1,4,1)=1.48
	  rco(1,5,1)=rco(1,4,1)
	  else
	  rco(1,3,1)=1.48
	  rco(1,4,1)=rco(1,3,1)
	  rco(1,5,1)=rco(1,3,1)
	endif
	aco(1,1,1)=0.631+0.016*z2/a213+(-0.001*z2/a213)*E2
	aco(1,2,1)=0.631+0.016*z2/a213
	aco(1,2,2)=-0.001*z2/a213
	aco(1,3,1)=aco(1,2,1)
	aco(1,3,2)=aco(1,2,2)
	if(E3.lt.E3aR)then
	  aco(1,4,1)=aco(1,2,1)
	  aco(1,4,2)=aco(1,2,2)
	  aco(1,5,1)=0.684-0.016*z2/a213
	  aco(1,5,2)=-0.0026+0.00026*z2/a213
	  else
	  aco(1,4,1)=0.684-0.016*z2/a213
	  aco(1,4,2)=-0.0026+0.00026*z2/a213
	  aco(1,5,1)=aco(1,4,1)
	  aco(1,5,2)=aco(1,4,2)
	endif
	pot(1,1,1)=168.0+0.733*z2/a213
	pot(1,1,2)=-2.64
	pot(1,2,1)=pot(1,1,1)
	pot(1,2,2)=pot(1,1,2)
	if(E3.lt.E3rR)then
	  pot(1,3,1)=116.5+0.337*z2/a213
	  pot(1,3,2)=-0.453
	  pot(1,4,1)=pot(1,3,1)
	  pot(1,4,2)=pot(1,3,2)
	  pot(1,5,1)=pot(1,3,1)
	  pot(1,5,2)=pot(1,3,2)
	  elseif(E3.lt.E3aR)then
	    pot(1,3,1)=pot(1,1,1)
	    pot(1,3,2)=pot(1,1,2)
	    pot(1,4,1)=116.5+0.337*z2/a213
	    pot(1,4,2)=-0.453
	    pot(1,5,1)=pot(1,4,1)
	    pot(1,5,2)=pot(1,4,2)
	    else
	    pot(1,3,1)=pot(1,1,1)
	    pot(1,3,2)=pot(1,1,2)
	    pot(1,4,1)=pot(1,1,1)
	    pot(1,4,2)=pot(1,1,2)
	    pot(1,5,1)=116.5+0.337*z2/a213
	    pot(1,5,2)=-0.453
	endif
c
	rco(2,1,1)=1.34
	aco(2,1,1)=0.50
	pot(2,1,1)=2.73-2.88*a213
	pot(2,1,2)=1.11
	rco(4,1,1)=1.52
	rco(4,2,1)=rco(4,1,1)
	rco(4,3,1)=rco(4,1,1)
	aco(4,1,1)=0.729-0.074*a213
	aco(4,2,1)=aco(4,1,1)
	aco(4,3,1)=aco(4,1,1)
	if(a2.le.130.)then
	  pot(4,1,1)=3.5
	  else
	  pot(4,1,1)=1.5
	endif
	pot(4,2,1)=22.2+4.57*a2**(1./3.)-7.446*epot(4,2)
	pot(4,2,2)=6.
	pot(4,3,1)=22.2+4.57*a2**(1./3.)
	pot(4,3,2)=-1.446
	rcoul(1)=1.3
c
      iz(1)=ifix(z2)
      ia(1)=ifix(a2)
      write(6,10)iref,author,refer
      write(6,11)(summary(m),m=1,240)
	
      izmin = iztarget
	izmax = iztarget
	iamin = iatarget-1
	iamax = iatarget+1

      write(6,12)emin,emax,izmin,izmax,iamin,iamax
      write(6,14)ifix(z2),ifix(a2)
	write(6,16)E1,E2,E3,E3aR
c
      if(E3.lt.E3rR)then
	write(6,21)rco(1,1,1),rco(1,1,2),epot(1,3),rco(1,4,1),epot(1,3)   	
	else
	write(6,21)rco(1,1,1),rco(1,1,2),epot(1,2),rco(1,4,1),epot(1,2)   	
      endif
      if(E3.lt.E3aR)then
	  if(aco(1,5,2).lt.0.0)then
	    write(6,22)aco(1,1,1),           epot(1,1),
     *	      	       aco(1,2,1),aco(1,2,2),epot(1,1),epot(1,4),
     *	               aco(1,5,1),aco(1,5,2),epot(1,4)
	    else
	    write(6,23)aco(1,1,1),           epot(1,1),
     *	      	       aco(1,2,1),aco(1,2,2),epot(1,1),epot(1,4),
     *	               aco(1,5,1),aco(1,5,2),epot(1,4)
	  endif
	  elseif(aco(1,4,2).lt.0.0)then
            write(6,22)aco(1,1,1),epot(1,1),
     *		       aco(1,2,1),aco(1,2,2),epot(1,1),epot(1,3),
     *	               aco(1,4,1),aco(1,4,2),epot(1,3)
	    else
            write(6,23)aco(1,1,1),epot(1,1),
     *		       aco(1,2,1),aco(1,2,2),epot(1,1),epot(1,3),
     *	               aco(1,4,1),aco(1,4,2),epot(1,3)
      endif
      if(E3.lt.E3rR)then
        write(6,24)pot(1,1,1),pot(1,1,2),epot(1,2),
     *		   pot(1,3,1),pot(1,3,2),epot(1,2) 
	elseif(E3.lt.E3aR)then
	  write(6,24)pot(1,1,1),pot(1,1,2),epot(1,3),
     *		     pot(1,4,1),pot(1,4,2),epot(1,3) 
	  else
	  write(6,24)pot(1,1,1),pot(1,1,2),epot(1,4),
     *		     pot(1,5,1),pot(1,5,2),epot(1,4) 
      endif
      write(6,26)rco(2,1,1)   	
      write(6,27)aco(2,1,1)   	
      write(6,28)pot(2,1,1),pot(2,1,2) 
      write(6,29)rco(4,1,1)   	
      write(6,30)aco(4,1,1)   	
      write(6,31)pot(4,1,1),epot(4,1),pot(4,2,1),pot(4,2,2),epot(4,1),
     *			    epot(4,2),pot(4,3,1),pot(4,3,2),epot(4,2)
c
      write(6,32)
      write(6,34)
      write(6,*)'E(L.S.)   V_R    r_R    a_R',
     *		     '     W_V   r_WV   a_WV',
     *		     '     W_D   r_WD   a_WD    r_C'
      write(6,*)' (MeV)   (MeV)   (fm)   (fm)',
     *		      '   (MeV)  (fm)   (fm)',
     *		     '    (MeV)  (fm)   (fm)    (fm)'
      write(6,34)
c
	do 200 l=1,50
	E=float(l)				!*0.1
	if(E.lt.E3rR)then
	  rR=rco(1,1,1)+rco(1,1,2)*E
	  else
	  rR=rco(1,4,1)
	endif	 
	if(E.lt.epot(1,1)) then
	  aR=aco(1,1,1)
	  elseif(E3.lt.E3aR)then
	    if(E.lt.epot(1,4)) then
	      aR=aco(1,2,1)+aco(1,2,2)*E
	      else
	      aR=aco(1,5,1)+aco(1,5,2)*E
	    endif
	    else
	    if(E.lt.epot(1,3)) then
	      aR=aco(1,2,1)+aco(1,2,2)*E
	      else
	      aR=aco(1,4,1)+aco(1,4,2)*E
	    endif
	endif	  
c
	if(E.lt.E3) then
	  VR=pot(1,1,1)+pot(1,1,2)*E
	  else
	  VR=pot(1,5,1)+pot(1,5,2)*E
	endif
c
	rWV=rco(2,1,1)
	aWV=aco(2,1,1)
	WV=amax1(pot(2,1,1)+pot(2,1,2)*E,0.)
	rWD=rco(4,1,1)
	aWD=aco(4,1,1)
	if(E.lt.epot(4,1)) then
	  WD=pot(4,1,1)
	  elseif(E.lt.epot(4,2)) then
	    WD=amax1(pot(4,2,1)+pot(4,2,2)*E,0.)
	    else
	    WD=amax1(pot(4,3,1)+pot(4,3,2)*E,0.)
	endif
c
	rC=rcoul(1)
  200 write(6,36)E,VR,rR,aR,WV,rWV,aWV,WD,rWD,aWD,rC
c
      write(6,34)
      write(6,38)
c
      call omout11
      goto 100
c
  999 STOP
    1 format(10i5)
    4 format(f10.3,5f10.4)
    9 format(80a)
c
   10 FORMAT(30X,'Optical Model Potential'//
     *'iref      = ',I5/
     *'authors   = ',80a/
     *'reference = ',80a/
     *'summary:')
   11 FORMAT(5x,80a)
   12 FORMAT('E range   =',F6.3,' - ',F8.4/
     *	     'Z range   =',I6,' - ',I3/
     *	     'A range   =',I6,' - ',I3/)
   14 FORMAT('Target: Z =',I3,'  A =',I4)
   16 FORMAT('        E1=',F6.3,'  E2=',F6.3,'  E3=',F6.3,'  E3a=',F6.3/
     *)
   21 FORMAT(' r_R(fm) =',F7.3,'+',F5.3'*E,     E <',F5.1,' MeV',/
     *       '         =',F7.3,  ',             E >',F5.1,' MeV'/)   
   22 FORMAT(' a_R(fm) =',F7.3,',                        E <',F7.3,' MeV
     *'     /'         =',F7.3,F8.5,'*E,'  ,F7.3,' MeV < E <',F7.3,' MeV
     *'     /'         =',F7.3,F8.5,'*E,              E >',F7.3,' MeV'/)   
   23 FORMAT(' a_R(fm) =',F7.3,',                        E <',F7.3,' MeV
     *'     /'         =',F7.3,F8.5,'*E,'  ,F7.3,' MeV < E <',F7.3,' MeV
     *' /'         =',F7.3,'+',F7.5,'*E,              E >',F7.3,' MeV'/)   
   24 FORMAT(' V_R(MeV)=',F6.2,F5.2,'*E,','    E <',F7.3,' MeV'/
     *       '         =',F6.2,F5.2,'*E,','    E >',F7.3,' MeV'/)   
   25 FORMAT(' V_R(MeV)=',F6.2,'+321/E**0.5')
   26 FORMAT(' r_WV(fm)=',F6.2)   
   27 FORMAT(' a_WV(fm)=',F6.2)   
   28 FORMAT(' W_V(MeV)=',F6.2,'+',F4.2,'*E'/)   
   29 FORMAT(' r_WD(fm)=',F6.2)   
   30 FORMAT(' a_WD(fm)=',F7.3)   
   31 FORMAT(' W_D(MeV)=',F6.2,',                       E <',F7.3,' MeV
     *'  /'         =',F7.2,'+',F5.3,'*E,',F7.3,' MeV < E <',F7.3,' MeV
     *'  /'         =',F7.2,F6.3,   '*E,              E >',F7.3,' MeV'/)
   32 FORMAT('----- Calculated Potential Parameters -----'/)
   34 FORMAT(80('-'))
   36 FORMAT(F6.1,F9.2,2F7.3,F8.3,2F7.3,F8.3,2F7.3,F7.2)
   38 FORMAT(/'----- RIPL Library Format -----'/)
      END
c
      subroutine omout11
c     ******************************************************************
c
c     routine to write optical model parameters into RIPL-2 library
c
c     ******************************************************************
      parameter (ndim1=10, ndim2=13, ndim3=25, ndim4=30, ndim5=10,
     1 ndim6=10,ndim7=120)
c
c     Common block and declarations for RIPL optical parameters
c
      character*1 author,refer,summary
      common/lib/author(80),refer(80),summary(320),iref,emin,emax,
     + izmin,izmax,iamin,iamax,imodel,jrange(6),epot(6,ndim1),
     + rco(6,ndim1,ndim2),aco(6,ndim1,ndim2),pot(6,ndim1,ndim3),
     + ncoll(ndim4),nvib(ndim4),nisotop,iz(ndim4),ia(ndim4),
     + lmax(ndim4),bandk(ndim4),def(ndim4,ndim5),idef(ndim4),
     + izproj,iaproj,exv(ndim7,ndim4),iparv(ndim7,ndim4),irel,idr,
     + nph(ndim7,ndim4),defv(ndim7,ndim4),thetm(ndim7,ndim4),
     + beta0(ndim4),gamma0(ndim4),xmubeta(ndim4),ex(ndim6,ndim4),
     + spin(ndim6,ndim4),ipar(ndim6,ndim4),spinv(ndim7,ndim4),
     + jcoul,ecoul(ndim1),rcoul(ndim1),rcoul0(ndim1),beta(ndim1),
     + rcoul1(ndim1),rcoul2(ndim1)
      common/io/ki,ko,ieof
c
      character*8 ldum
      dimension rdum(ndim2),adum(ndim2),pdum(ndim3)
c
  1   format(10i5)
  2   format(10a8)
  4   format(f10.3,7f10.4)
  5   format(f12.5,(5(1x,1pe11.4)),/12x,(5(1x,1pe11.4)))
  6   format(5i5,f5.1,(4(1x,1pe10.3)),/30x,(4(1x,1pe10.3)))
  7   format(f12.8,f7.1,2i4,1p,2(1x,e11.4))
  8   format(2i5,1p,3(1x,e11.4))
  9   format(80a1)
 11   format('Target:  Z = ',i2,'   A = ',i3)
c
      data ldum/'++++++++'/
c
      write(ko,1) iref
      write(ko,9) author
      write(ko,9) refer
      write(ko,9) (summary(m),m=1,240)
      write(ko,11)iz(1),ia(1)
      write(ko,4) emin,emax
      write(ko,1) izmin,izmax
      write(ko,1) iamin,iamax
      write(ko,1) imodel,izproj,iaproj,irel,idr
      do 100 i=1,6
        write(ko,1) jrange(i)
        if(jrange(i).eq.0) go to 100
        krange=iabs(jrange(i))
        do 98 j=1,krange
          write(ko,4) epot(i,j)
c         write(ko,5) (rco(i,j,n),n=1,ndim2)
c         write(ko,5) (aco(i,j,n),n=1,ndim2)
c         write(ko,5) (pot(i,j,n),n=1,ndim3)
          do 82 n=1,ndim2
          rdum(n)=rco(i,j,n)
  82      adum(n)=aco(i,j,n)
          do 84 n=1,ndim3
  84      pdum(n)=pot(i,j,n)
          call linew3(rdum,ndim2,ko)
          call linew3(adum,ndim2,ko)
          call linew3(pdum,ndim3,ko)
  98    continue
 100  continue
      write(ko,1)jcoul
      if(jcoul.le.0) go to 110
      do 108 j=1,jcoul
 108  write(ko,4)ecoul(j),rcoul0(j),rcoul(j),rcoul1(j),rcoul2(j),
     >           beta(j), 0., 0.
 110  if(imodel.ne.1) go to 130
      write(ko,1) nisotop
      do 120 n=1,nisotop
        write(ko,6) iz(n),ia(n),ncoll(n),lmax(n),idef(n),bandk(n),
     1   (def(n,k),k=2,idef(n),2)
        do 124 k=1,ncoll(n)
          write(ko,7) ex(k,n),spin(k,n),ipar(k,n)
 124    continue
 120  continue
      go to 200
 130  if(imodel.ne.2) go to 150
      write(ko,1) nisotop
      do 140 n=1,nisotop
        write(ko,1) iz(n),ia(n),nvib(n)
        do 138 k=1,nvib(n)
          write(ko,7) exv(k,n),spinv(k,n),iparv(k,n),nph(k,n),defv(k,n),
     1      thetm(k,n)
 138    continue
 140  continue
      go to 200
 150  if(imodel.ne.3)go to 200
      write(ko,1) nisotop
      do 160 n=1,nisotop
        write(ko,8) iz(n),ia(n),beta0(n),gamma0(n),xmubeta(n)
 160  continue
 200  continue
      write(ko,2)(ldum,l=1,10)
      return
      end
c
      subroutine linew3(ah,nah,ko)
c     ******************************************************************
c
c     utility for omout RIPL single line writes
c     requires subroutine cxfp
c
c     ******************************************************************
c      dimension ah(25),f(25),is(25),j(25)
      dimension ah(nah), f(nah), is(nah), j(nah)
      character*8 start1(2),start2(2),fmt(15),f8(2),f7(2),finish
      character*1 is
      data start1/'(f12.5, ','1x      '/
      data start2/'(13x,   ','        '/
      data f8/'f8.5,  ','a1,i1,1x'/
      data f7/'f7.4,  ','a1,i2,1x'/
      data finish/')'/
c
c     Write first line
c
      l=0
      do 110 k=1,2
        l=l+1
        fmt(l)=start1(k)
 110  continue
      do 120 k=2,7
        call cxfp(ah(k),f(k),is(k),j(k))
        l=l+1
        if(iabs(j(k)).ge.10)go to 115
        fmt(l)=f8(1)
        l=l+1
        fmt(l)=f8(2)
        go to 120
 115    fmt(l)=f7(1)
        l=l+1
        fmt(l)=f7(2)
 120  continue
      l=l+1
      fmt(l)=finish
      write(ko,fmt)ah(1),(f(k),is(k),j(k),k=2,7)
c
c     Write remaining full lines
c
      nleft=nah-7
      nloop=nleft/6
      nremain=nleft-nloop*6
      nold=1
      if(nloop.le.0)go to 165
      do 160 nlp=1,nloop
        nold=nold+6
        l=0
        do 140 k=1,2
          l=l+1
          fmt(l)=start2(k)
 140    continue
        do 150 k=1,6
          call cxfp(ah(k+nold),f(k),is(k),j(k))
          l=l+1
          if(iabs(j(k)).ge.10)go to 145
          fmt(l)=f8(1)
          l=l+1
          fmt(l)=f8(2)
          go to 150
 145      fmt(l)=f7(1)
          l=l+1
          fmt(l)=f7(2)
 150    continue
        l=l+1
        fmt(l)=finish
        write(ko,fmt)(f(k),is(k),j(k),k=1,6)
 160  continue
 165  if(nremain.le.0) go to 200
c
c     Write remaining partial line
c
      nold=nold+6
      l=0
      do 170 k=1,2
        l=l+1
        fmt(l)=start2(k)
 170  continue
      do 180 k=1,nremain
        call cxfp(ah(k+nold),f(k),is(k),j(k))
        l=l+1
        if(iabs(j(k)).ge.10)go to 175
        fmt(l)=f8(1)
        l=l+1
        fmt(l)=f8(2)
        go to 180
 175    fmt(l)=f7(1)
        l=l+1
        fmt(l)=f7(2)
 180  continue
      l=l+1
      fmt(l)=finish
      write(ko,fmt)(f(k),is(k),j(k),k=1,nremain)
 200  return
      end
c
      subroutine cxfp (x,f,s,n)
c     ******************************************************************
c     convert x for punching.
c        x - floating point number = f*10.0**n
c        f - 0.99995 le f lt 9.999995
c        s  sign (hollerith + or -) of exponent
c        n - exponent
c     ******************************************************************
cibm
c     real*8 xx,ff
cibm
      character*1 s,sp,sm
      data sp/'+'/,sm/'-'/
      if (x.ne.0.0) go to 100
      f=0.0
      s=sp
      n=0
      return
ccdc
  100 n=alog10(abs(x))
      if (abs(x).lt.1.0) go to 140
      f=x/10.0**n
      s=sp
      if (iabs(n).lt.10.and.abs(f).lt.9.9999995) go to 170
      if (iabs(n).ge.10.and.abs(f).lt.9.999995) go to 170
      f=f/10.0
      n=n+1
      go to 170
  140 n=1-n
      f=x*10.0**n
      s=sm
      if (iabs(n).lt.10.and.abs(f).lt.9.9999995) go to 170
      if (iabs(n).ge.10.and.abs(f).lt.9.999995) go to 170
      f=f/10.0
      n=n-1
      if (n.gt.0) go to 170
      s=sp
  170 continue
ccdc
cibm
c 100 xx=x
c     n=dlog10(dabs(xx))
c     if (dabs(xx).lt.1.0d0) go to 140
c     ff=xx/10.d0**n
c     s=sp
c     if (iabs(n).lt.10.and.dabs(ff).lt.9.9999995d0) go to 170
c     if (iabs(n).ge.10.and.dabs(ff).lt.9.999995d0) go to 170
c     ff=ff/10.0d0
c     n=n+1
c     go to 170
c 140 n=1-n
c     ff=xx*10.0d0**n
c     s=sm
c     if (iabs(n).lt.10.and.dabs(ff).lt.9.9999995d0) go to 170
c     if (iabs(n).ge.10.and.dabs(ff).lt.9.999995d0) go to 170
c     ff=ff/10.0d0
c     n=n-1
c     if (n.gt.0) go to 170
c     s=sp
c 170 f=ff
cibm
      return
      end
 
