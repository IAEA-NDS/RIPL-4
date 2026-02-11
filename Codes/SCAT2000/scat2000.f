      program SCAT2000
c***********************************************************************
c                         *                 *                          *
c                         *  PROGRAM SCAT2  *                          *
c                         *                 *                          *
c                         *******************                          *
c                                                                      *
c                          JULY      1977                              *
c                  REV.1   JULY      1979                              *
c                  REV.2   JULY      1982                              *
c                  REV.3   OCTOBER   1991                              *
c                  REV.4   NOVEMBER  1993                              *
c                  REV.5   JUNE      1994                              *
c                  REV.6   FEBRUARY  1999 Fortran 90 compatibility     *
c                  REV.7   OCTOBER   2000 Heavy ion option from (1)    *
c                  REV.7a  OCTOBER   2000 Relativistic kinematics      *
c                  REV.7b  OCTOBER   2000 Surface real potential       *
c                  REV.7c  OCTOBER   2000 Imaginary spin-orbit pot.    *
c                                                                      *
c  (1) Anisoara CONSTANTINESCU, Bucarest University, Romania           *
c      private communication                                           *
c***********************************************************************
c  AUTHOR:                                                             *
c     O.BERSILLON                                                      *
c     SERVICE DE PHYSIQUE NUCLEAIRE                                    *
c     CENTRE D'ETUDES DE BRUYERES-LE-CHATEL                            *
c     B.P. 12                                                          *
c     91680 BRUYERES-LE-CHATEL       FRANCE                            *
c     TEL.: +33 1 69.26.54.14                                          *
c     E-MAIL: olivier.bersillon@cea.fr                                 *
c***********************************************************************
c  REFERENCE:                                                          *
c     O.BERSILLON                                                      *
c     SCAT2 : UN PROGRAMME DE MODELE OPTIQUE SPHERIQUE                 *
c     CEA-N-2227, NEANDC(E) 220'L', INDC(FR) 49/L                      *
c     OCTOBER 1981                                                     *
c***********************************************************************
c  Input files:                                                        *
c  ------------                                                        *
c  ie1      14  Names of the input and output files      seq.   80   f *
c  ie       15  Input data file                          seq.   80   f *
c                                                                      *
c  Output files:                                                       *
c  -------------                                                       *
c  is       16  Output listing file                      seq.  133   f *
c  is1      17  Trans.coeff.   in GNASH or STAPRE format seq.   80   f *
c  is2      18  Temporary file in GNASH or STAPRE format seq.   80   f *
c  is3      19  Summary output                           seq.   80   f *
c----------------------------------------------------------------------*
c=====  The names of the five files (15 to 19) must be read from  =====*
c=====      the SCAT2000.DAT file (ie1=14)                        =====*
c***********************************************************************
c  INPUT DATA:                                                         *
c  ===========                                                         *
c   read(ie,*) ipr,ida,iba,ipu,idr,irelat                              *
c   ****                                                               *
c       IPR   Print option                                             *
c           =  1  print transmission coefficients T(l,j)               *
c           =  0  no print of T(l,j)                                   *
c       IDA   Angular distribution calculation option                  *
c           =  1  calculate shape elastic angular distribution for     *
c                 uniformly spaced angles from 0 to 180 degrees by     *
c                 step of 5 degrees                                    *
c           =  2  calculate shape elastic angular distribution for     *
c                 uniformly spaced cosines from 1 to -1 by             *
c                 step of 0.02                                         *
c           =  0  no angular distribution calculation                  *
c           = -1  calculate also Legendre polynomial coefficients      *
c                 (neutrons only)                                      *
c           = -2  calculate also Legendre polynomial coefficients      *
c                 (neutrons only)                                      *
c       IBA   Output option on the file IS1                            *
c           =  1  print optical parameters, cross sections and         *
c                 transmission coefficients in the GNASH format        *
c           =  2  print transmission coefficients in the STAPRE format *
c           =  0  don't print                                          *
c       IPU   Output option on the file IS3                            *
c           =  1  print cross sections and angular distributions       *
c           =  0  don't print                                          *
c       IDR   Usage of dispersion relations                            *
c           =  2  use exact dispersion correction                      *
c                    (volume + surface  real potential)                *
c           =  1  use equivalent volume real potential                 *
c           =  0  don't use                                            *
c       IRELAT = 0 classical    kinematics                             *
c                1 relativistic kinematics                             *
c                                                                      *
c 10 continue                                                          *
c   read(ie,*) ne                                                      *
c   ****                                                               *
c       NE    Number of incident energies                              *
c                                                                      *
c   read(ie,*) (en(j),j=1,ne)                                          *
c   ****                                                               *
c       EN(J) Incident energies (MeV)                                  *
c             if EN(1) > 0.  center of mass energies                   *
c                      < 0.  laboratory energies                       *
c                                                                      *
c 20 continue                                                          *
c   read(ie,*) izt,imt                                                 *
c   ****                                                               *
c       IZT   Atomic number of the target                              *
c       IMT   Mass   number of the target                              *
c                                                                      *
c 30 continue                                                          *
c   read(ie,*) ip,ipot                                                 *
c   ****                                                               *
c       IP    Type of incident particle                                *
c       IPOT  Type of built-in parameters                              *
c                                                                      *
c       IP   =  1 neutron      IPOT = 1 --> Wilmore - Hodgson          *
c                                     2 --> Becchetti - Greenlees      *
c                                     3 --> Ferrer - Carlson - Rapaport*
c                                     4 --> Bersillon - Cindro         *
c                                     5 --> Madland (actinides)        *
c            =  2 proton       IPOT = 1 --> Perey                      *
c                                     2 --> Beccheti - Greenless       *
c            =  3 deuteron     IPOT = 1 --> Lohr - Haeberli            *
c                                     2 --> Perey                      *
c            =  4 triton       IPOT = 1 --> Becchetti - Greenlees      *
c            =  5 helium-3     IPOT = 1 --> Becchetti - Greenlees      *
c            =  6 alpha        IPOT = 1 --> averaged potential         *
c            =  7 heavy ion option                                     *
c    if IP   =  7 then                                                 *
c       read(ie,*) mi,zi,si                                            *
c       ****                                                           *
c           MI = heavy ion mass number                                 *
c           ZI = heavy ion charge                                      *
c           SI = heavy ion spin (restricted to 0, 1/2 or 1)            *
c       read(ie,*) tetaf,tetastep                                      *
c       ****                                                           *
c           TETAF    = maximum angle up to which the                   *
c                      angular distribution is calculated              *
c           TETASTEP = angle step for the angular distribution         *
c                                                                      *
c    if IPOT =  0  optical parameters have to be read                  *
c                                                                      *
c VOLUME  REAL POTENTIAL: Woods-Saxon                    potential # 1 *
c ----------------------                                               *
c   read(ie,*) r(1),re(1),a(1),ae(1),npzen(1)                          *
c   ****                                                               *
c        if NPZEN(1) > 0 strength parameters are read                  *
c                    < 0 volume integrals    are read                  *
c                                                                      *
c        for each energy range j=1,npzen(1)                            *
c     read(ie,*) epot(1,j),(pot(1,j,k),k=1,npcofx)                     *
c     ****                                                             *
c                                                                      *
c       R (1)      = radius (fm)                                       *
c       RE(1)      = linear energy dependence of the radius            *
c       A (1)      = diffuseness (fm)                                  *
c       AE(1)      = linear energy dependence of the diffuseness       *
c       NPZEN(1)   = number of energy range for the potential          *
c       EPOT(1,J)  = upper limit of the j-th energy range              *
c       POT(1,J,K) = strength parameters (k=1,npcofx; npcofx=7)        *
c                                                                      *
c       if POT(1,J,7) .eq. 0                                           *
c           v =   pot(1,j,1)                                           *
c               + pot(1,j,2)*e                                         *
c               + pot(1,j,3)*e*e                                       *
c               + pot(1,j,4)*e*e*e                                     *
c               + pot(1,j,5)*LN(e)                                     *
c               + pot(1,j,6)*SQRT(e)                                   *
c       if POT(1,J,7) .ne. 0                                           *
c           v =   pot(1,j,7)*EXP( pot(1,j,2)*u )                       *
c                 where u = e - ef                                     *
c                                                                      *
c SURFACE REAL POTENTIAL                                 potential # 2 *
c ----------------------                                               *
c   read(ie,*) r(2),re(2),a(2),ae(2),npzen(2)                          *
c   ****                                                               *
c        for each energy range j=1,npzen(2)                            *
c     read(ie,*) epot(2,j),(pot(2,j,k),k=1,npcofx)                     *
c     ****                                                             *
c       if POT(2,J,7) .eq. 0                                           *
c           v =   pot(2,j,1)                                           *
c               + pot(2,j,2)*e                                         *
c               + pot(2,j,3)*e*e                                       *
c               + pot(2,j,4)*e*e*e                                     *
c               + pot(2,j,5)*LN(e)                                     *
c               + pot(2,j,6)*SQRT(e)                                   *
c                                                                      *
c VOLUME  IMAGINARY POTENTIAL: Woods-Saxon               potential # 3 *
c ---------------------------                                          *
c   read(ie,*) r(3),re(3),a(3),ae(3),npzen(3)                          *
c   ****                                                               *
c        for each energy range j=1,npzen(3)                            *
c     read(ie,*) epot(3,j),(pot(3,j,k),k=1,npcofx)                     *
c     ****                                                             *
c       if POT(3,J,7) .eq. 0                                           *
c           v =   pot(3,j,1)                                           *
c               + pot(3,j,2)*e                                         *
c               + pot(3,j,3)*e*e                                       *
c               + pot(3,j,4)*e*e*e                                     *
c               + pot(3,j,5)*LN(e)                                     *
c               + pot(3,j,6)*SQRT(e)                                   *
c       if POT(3,J,7) .ne. 0                                           *
c           v =   pot(3,j,1)*(u**n1) / (u**n1 + pot(3,j,2)**n1)        *
c                 where u  = e - ef                                    *
c                       n1 = pot(3,j,3)                                *
c                                                                      *
c SURFACE IMAGINARY POTENTIAL: if r(4) > 0., Woods-Saxon derivative    *
c ---------------------------  if r(4) < 0., Gaussian    potential # 4 *
c   read(ie,*) r(4),re(4),a(4),ae(4),npzen(4)                          *
c   ****                                                               *
c        for each energy range j=1,npzen(4)                            *
c     read(ie,*) epot(4,j),(pot(4,j,k),k=1,npcofx)                     *
c     ****                                                             *
c       if POT(4,J,7) .eq. 0                                           *
c           v =   pot(4,j,1)                                           *
c               + pot(4,j,2)*e                                         *
c               + pot(4,j,3)*e*e                                       *
c               + pot(4,j,4)*e*e*e                                     *
c               + pot(4,j,5)*LN(e)                                     *
c               + pot(4,j,6)*SQRT(e)                                   *
c       if POT(4,J,7) .ne. 0                                           *
c           v =   pot(4,j,1)*(u**n1) / (u**n1 + pot(4,j,2)**n1) +      *
c                 pot(4,j,4)*(u**n2) / (u**n2 + pot(4,j,5)**n2)        *
c                 where u  = e - ef                                    *
c                       n1 = pot(4,j,3)                                *
c                       n2 = pot(4,j,6)                                *
c                                                                      *
c REAL      SPIN-ORBIT POTENTIAL                         potential # 5 *
c ------------------------------                                       *
c   read(ie,*) r(5),re(5),a(5),ae(5),npzen(5)                          *
c   ****                                                               *
c        for each energy range j=1,npzen(4)                            *
c     read(ie,*) epot(5,j),(pot(5,j,k),k=1,npcofx)                     *
c     ****                                                             *
c       if POT(5,J,7) .eq. 0                                           *
c           v =   pot(5,j,1)                                           *
c               + pot(5,j,2)*e                                         *
c               + pot(5,j,3)*e*e                                       *
c               + pot(5,j,4)*e*e*e                                     *
c               + pot(5,j,5)*LN(e)                                     *
c               + pot(5,j,6)*SQRT(e)                                   *
c       if POT(5,J,7) .ne. 0                                           *
c           v =   pot(5,j,7)*EXP( pot(5,j,2)*u )                       *
c                 where u = e - ef                                     *
c                                                                      *
c IMAGINARY SPIN-ORBIT POTENTIAL                         potential # 6 *
c ------------------------------                                       *
c   read(ie,*) r(6),re(6),a(6),ae(6),npzen(6)                          *
c   ****                                                               *
c        for each energy range j=1,npzen(6)                            *
c     read(ie,*) epot(6,j),(pot(6,j,k),k=1,npcofx)                     *
c     ****                                                             *
c       if POT(6,J,7) .eq. 0                                           *
c           v =   pot(6,j,1)                                           *
c               + pot(6,j,2)*e                                         *
c               + pot(6,j,3)*e*e                                       *
c               + pot(6,j,4)*e*e*e                                     *
c               + pot(6,j,5)*LN(e)                                     *
c               + pot(6,j,6)*SQRT(e)                                   *
c       if POT(6,J,7) .ne. 0                                           *
c           v =   pot(6,j,1)*(u**n1) / (u**n1 + pot(6,j,2)**n1)        *
c                 where u  = e - ef                                    *
c                       n1 = pot(6,j,3)                                *
c                                                                      *
c COULOMB RADIUS                                                       *
c --------------                                                       *
c   read(ie,*) rcoul,beta                                              *
c   ****                                                               *
c       RCOUL  = Coulomb radius                                        *
c       BETA   = nonlocality range                                     *
c                if BETA .ne. 0. the imaginary potential is pure DWS   *
c                                                                      *
c RECYCLE OPTION                                                       *
c --------------                                                       *
c   read(ie,*) isuit                                                   *
c   ****                                                               *
c     ISUIT = 0 end                                                    *
c           = 1 new complete case                             (goto 10)*
c           = 2 the energy grid is conserved for a new system (goto 20)*
c           = 3 the projectile or the parameters are changed  (goto 30)*
c                                                                      *
c  If IBA = 2 only                                                     *
c   read(ie,*) titre                                                   *
c   ****                                                               *
c     TITRE = comment card read in the TLSAP subroutine                *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'angles2.i'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'ener.i'
      INCLUDE 'fact.i'
      INCLUDE 'inout.i'
      INCLUDE 'physcon.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
      INCLUDE 'tce.i'
      INCLUDE 'xs.i'
c
      character        ref*10
      double precision ampi
      integer          iba,ida,idr,iperm,ipr,ipu,irelat,isuit
      integer          i,imt,izt,ip,ipl,ipot,j,k
      integer          lmax,lmaxim,lmaxtl,n,na,ne
c
      dimension        iperm(nemax)
c
c  Incident particles data (mass, charge, spin)
      double precision ami   ,azi   ,asi
      dimension        ami(6),azi(6),asi(6)

c     to call energy2() subroutine
      real             tarmas,tarspi,tarex,zatar
      character*1      tarpar
c                n        p        d        t       He3      He4
c     data ami /1.0     ,1.0     ,2.0     ,3.0     ,3.0     ,4.0     /
      data ami /1.008665,1.007825,2.014102,3.016050,3.016030,4.002603/
      data asi /0.5     ,0.5     ,1.0     ,0.5     ,0.5     ,0.0     /
      data azi /0.      ,1.      ,1.      ,1.      ,2.      ,2.      /
c
      double precision zero,one,two,three
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
      data             three /3.0d+00/
c=======================================================================
c
c-----------------------------------------------------------------------
c  Open files
c-----------------------------------------------------------------------
c
      call FILES(1)
c
c-----------------------------------------------------------------------
c  Physical  constants
c-----------------------------------------------------------------------
c
      ampipm = 1.395688d+02
      ampi0  = 1.349645d+02
      ampi   = (two*ampipm + ampi0)/three
c     amu0c2 = 9.3149386d+02
      e2     = 1.4399652d+00
c     hbarc  = 1.97327053d+02
c-------------------------------------------------------------------------------
c  Physical  constants (following ECIS03 definitions and ENDF manual 2001)
c-------------------------------------------------------------------------------
C CONSTANTS COMPUTED FROM THE FUNDAMENTAL CONSTANTS, ATOMIC MASS, HBAR*CCALC-118
C AND ALPHA, AS GIVEN IN THE EUROPEAN PHYSICAL JOURNAL, PAGE 73, VOLUME CALC-119
C 15 (2000) REFERRING FOR THESE VALUES TO THE 1998 CODATA SET WHICH MAY CALC-120
C BE FOUND AT http://physics.nist.gov/constants                         CALC-121
C     amu0c2=931.494013 +/- 0.000037 MeV                                CALC-122
      amu0c2=931.494013D0                                               CALC-123
C     hbarc=197.3269601 +/- 0.0000078 (*1.0E-9 eV*cm)                   CALC-124
      hbarc=197.3269601D0                                               CALC-125
C     CZ=137.03599976 +/- 0.00000050 without dimension                  CALC-126
C     CZ=137.03599976D0                                                 CALC-127
C
      ck2    = (two*amu0c2)/(hbarc**2)
      ceta   = e2*DSQRT(amu0c2/two)/hbarc
      cso    = hbarc/ampi
      cso    = cso*cso
c
      epotmx = 999.0d+00
      lmaxim = 0
c
c-----------------------------------------------------------------------
c  Numerical constants
c-----------------------------------------------------------------------
c
      pi     = 3.141592654d+00
c
c-----------------------------------------------------------------------
c  Input data
c-----------------------------------------------------------------------
c
      read(ie,*) ipr,ida,iba,ipu,idr,irelat
c     ****
c
   10 continue
      read(ie,*) ne
c     ****
      read(ie,*) (en(j),j=1,ne)
c     ****
c
   20 continue
      read(ie,*) izt,imt

      zatar=float(1000*izt+imt)
      call energy2(zatar,tarmas,tarspi,tarpar,tarex,amu0c2)

c     read(ie,*) izt, mt
c     ****
      zt = DBLE(izt)
c     mt = DBLE(imt)
      mt = dble(tarmas) 
c
   30 continue
      read(ie,*) ip,ipot
c     ****
c
c  Special heavy ion option
      if(ip .eq. 7) then
        read(ie,*) mi,zi,si
c       ****
        read(ie,*) tetaf,tetastep
c       ****
        ipot = 0
      else
        mi = ami(ip)
        zi = azi(ip)
        si = asi(ip)
      endif
c
c  Angular distribution option
      if(ida .ne. 0) then
        call FCT(lfmax,g)
        call PREANG(ida,ip,na)
      endif
c
      ipl = IDNINT(two*si + one)
c
c-----------------------------------------------------------------------
c  Transformation of energies from laboratory to center-of-mass
c-----------------------------------------------------------------------
c
      if(en(1) .lt. zero) then
        do n=1,ne
          call KINEMA(DABS(en(n)),ein(n),mi,mt,amu,ak2,1,irelat)
          en (n) = DABS(en(n))
        end do
      else
        do n=1,ne
          ein(n) = en(n)
        end do
      endif
c
c  Sort center-of-mass energies by increasing values
c
      call SORTR(ein,ne,iperm)
c
c  Potential parameters initialisations
c
      do i=1,nptypx
        a (i) = zero
        ae(i) = zero
        r (i) = zero
        re(i) = zero
        intpot(i)  = 0
        npzen(i)   = 1
        epot (i,0) = zero
        do j=1,npzenx
          do k=1,npcofx
            pot(i,j,k) = zero
          end do
        end do
      end do
      beta  = zero
      eferm = zero
      rcoul = zero
c
      do n=1,nemax
        do j=1,ltcmax
          tc(n,j) = zero
        end do
      end do
c
c-----------------------------------------------------------------------
c  Read optical potentials parameters
c-----------------------------------------------------------------------
c
      if(ipot .ne. 0) then
        call SYSPOT(ip,ipot,ref)
      else
c
        ref = 'USER'
        do i=1,nptypx
          read(ie,*) r(i),re(i),a(i),ae(i),npzen(i)
c         ****
          if(npzen(i) .lt. 0) then
            npzen(i)  = IABS(npzen(i))
            intpot(i) = 1
          endif
          do j=1,npzen(i)
            read(ie,*) epot(i,j),(pot(i,j,k),k=1,npcofx)
c           ****
          end do
          epot(i,npzen(i)) = epotmx
        end do
c
c  Coulomb radius, Fermi energy, non-locatily range
        read(ie,*) rcoul,eferm,beta
c       ****
        if(eferm .gt. zero) then
          print *,'E-Fermi must be negative'
          stop    'E-Fermi'
        endif
c
      endif
c
      read(ie,*) isuit
c     ****
c
      call PRIPOT(idr,irelat,ref)
      call WCHK(3)
      call WCHK(4)
c
c  End of input data
c
      lmaxtl = 0
c
c=======================================================================
c  Loop on incident energies
c=======================================================================
c
      do n=1,ne
        e1 = ein(n)
c  E1 is the current CM energy
        call KINEMA(el,e1,mi,mt,amu,ak2,2,irelat)
        if(ipr .eq. 1) write(is,fmt='(1h1)')
c
        call SCAT(ip,lmax,idr)
c
        if(ipl .eq. 1) then
          call SPIN0 (n,lmax,ipr,ipu)
        else if(ipl .eq. 2) then
          call SPIN05(n,lmax,ipr,ipu)
        else if(ipl .eq. 3) then
          call SPIN1 (n,lmax,ipr,ipu)
        endif
c
        if(lmax .gt. lmaxtl) lmaxtl = lmax
c
        if(ida .lt. 0) then
          call SHAPEL2(lmax,ida,na,ipl)
        else if(ida .gt. 0) then
          call SHAPEC (lmax,ida,na,ipl,ipu)
        endif
c
        if(iba .eq. 1) call TPUN(ip,ne,ipl,ipot,0)
      end do
c
c-----------------------------------------------------------------------
c  Print final results
c-----------------------------------------------------------------------
c
      call PRITC(ne)
c
      if(iba .eq. 1) then
        call TPUN(ip,ne,ipl,ipot,1)
      else if(iba .eq. 2) then
        call TLSTAP(ne,lmaxtl,0)
      endif
      if(lmaxtl .gt. lmaxim) lmaxim = lmaxtl
c
      if(isuit .eq. 1) then
        go to 10
      else if(isuit .eq. 2) then
        go to 20
      else if(isuit .eq. 3) then
        go to 30
      endif
c
      if(iba .eq. 2) call TLSTAP(ne,lmaxtl,lmaxim)
c
c-----------------------------------------------------------------------
c  Close files
c-----------------------------------------------------------------------
c
      call FILES(0)
c
      STOP 'End of SCAT2'
c
      end


      complex*16 function CGAMMA(z)
c***********************************************************************
c  Calculates complex gamma function                                   *
c----------------------------------------------------------------------*
c  Reference:   Y.L.Luke                                               *
c               The special functions and their approximations         *
c               vol.2, Academic Press, New York and London             *
c               (1969) p.304-305                                       *
c***********************************************************************
c
      implicit none
c
      complex*16       z
c
      complex*16       h,s,u,v
      double precision fk,fk1,x
      integer          k,l
c
      double precision g
      dimension        g(16)
      data             g/
     * 41.62443 69164 39068,-51.22424 10223 74774,+11.33875 58134 88977,
     * -0.74773 26877 72388, +0.00878 28774 93061, -0.00000 18990 30264,
     * +0.00000 00019 46335, -0.00000 00001 99345, +0.00000 00000 08433,
     * +0.00000 00000 01486, -0.00000 00000 00806, +0.00000 00000 00293,
     * -0.00000 00000 00102, +0.00000 00000 00037, -0.00000 00000 00014,
     * +0.00000 00000 00006/
c
      double precision zero,half,one,const
      data             zero  /0.0d+00/
      data             half  /0.5d+00/
      data             one   /1.0d+00/
      data             const /2.506628274631001d+00/
      double precision pi
      data             pi    /3.14159 26535 89793d+00/
c=======================================================================
c
      u = z
      x = DBLE(u)
      if(x .ge. one ) go to 3
      if(x .ge. zero) go to 2
      v = one - u
      l = 1
      go to 11
    2 v = u + one
      l = 2
      go to 11
    3 v = u
      l = 3
c
   11 h = one
      s = g(1)
      do k=2,16
        fk  = k - 2
        fk1 = fk + one
        h   = ((v - fk1)/(v + fk))*h
        s   = s + g(k)*h
      end do
      h = v + 4.5d+00
c90   cgamma = const*CDEXP((v - half)*CDLOG(h) - h)*s
      cgamma = const*  EXP((v - half)*  LOG(h) - h)*s
c
      if(l .eq. 1) then
c90     cgamma = pi/(CDSIN(pi*u)*cgamma)
        cgamma = pi/(  SIN(pi*u)*cgamma)
        return
      else if(l .eq. 2) then
        cgamma = cgamma/u
        return
      else if(l .eq. 3) then
        return
      endif
c
      end
      double precision function CLEBG(aj1,aj2,aj3,am1,am2,am3)
c***********************************************************************
c  Calculate Clebsch-Gordan coefficients                               *
c                                                                      *
c  Attention:  cg(j1,j2,j3,;m1,m2,m3) = (-1)**(j1+j2-m3)*              *
c                                       3-J(j1,j2,j3;m1,m2,-m3)        *
c                                       ---                            *
c                                                                      *
c  from    John.G. Wills     ORNL-TM-1949 (August 1967)                *
c                         et Comp.Phys.Comm. 2(1971)381                *
c                                                                      *
c  O.Bersillon     August 1977                                         *
c                                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'fact.i'
c
      double precision aj1,aj2,aj3,am1,am2,am3
c
      double precision a,b,c,d,e,f,h,q,s,t,x
      integer          i
      integer          i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,il
      integer          j,j1,j2,j3,k,l,la,lb,m,m1,m2,m3,n
c
      dimension        i(11)
      equivalence
     *(i(1),i1),(i(2),i2),(i(3),i3),(i( 4), i4),(i( 5), i5),(i(6),i6),
     *(i(7),i7),(i(8),i8),(i(9),i9),(i(10),i10),(i(11),i11)
c
      double precision zero,eps,one,two
      data             zero /0.0d+00/
      data             eps  /1.0d-03/
      data             one  /1.0d+00/
      data             two  /2.0d+00/
c=======================================================================
c
      clebg = zero
c
c  Convert the arguments to integer
c
      j1 = IDINT(two*aj1 + eps)
      j2 = IDINT(two*aj2 + eps)
      j3 = IDINT(two*aj3 + eps)
      m1 = IDINT(two*am1 + SIGN(eps,am1))
      m2 = IDINT(two*am2 + SIGN(eps,am2))
      m3 = IDINT(two*am3 + SIGN(eps,am3))
c
c  Test m1 + m2 = m3
c
      if(m1 + m2 - m3 .ne. 0) return
c
c  Test table size
c
      i(10) = (j1 + j2 + j3)/2 + 2
      n     = i(10)
      i(11) = j3 + 2
      if(i(10) .gt. lfmax) then
        print 9010, i(10),lfmax,aj1,aj2,aj3,am1,am2,am3
        return
      endif
c
      i(1) = j1 + j2 - j3
      i(2) = j2 + j3 - j1
      i(3) = j3 + j1 - j2
      i(4) = j1 - m1
      i(5) = j1 + m1
      i(6) = j2 - m2
      i(7) = j2 + m2
      i(8) = j3 - m3
      i(9) = j3 + m3
c
c  Check i(j) = even, triangular inequality, m less than j,
c    find number of terms
c
      do j=1,9
        k = i(j)/2
        if(i(j) .ne. 2*k) return
        if(k .lt. 0) return
        if(k .lt. n) then
          n = k
        endif
        i(j) = k + 1
      end do
c
      if(m3 .ne. 0 .or. m1 .ne. 0 .or. m1 .ne. 1) then
        il = 0
        la = i1 - i5
        lb = i1 - i6
        if(il .lt. la) then
          il = la
        endif
        if(il .lt. lb) then
          il = lb
        endif
c
c  Form coefficients of sum
c
        c  = (g(i11) - g(i11-1) + g(i1) + g(i2) + g(i3) - g(i10) +
     *        g(i4) + g(i5) + g(i6) + g(i7) + g(i8) + g(i9))/two
        j1 = i1 - il
        j2 = i4 - il
        j3 = i7 - il
        m1 = il + 1
        m2 = il - la + 1
        m3 = il - lb + 1
        c  = c - g(j1) - g(j2) - g(j3) - g(m1) - g(m2) - g(m3)
        c  = DEXP(c)
        if((il - 2*(il/2)) .ne. 0) c = -c
        if(n .lt. 0) return
        if(n .eq. 0) then
          clebg = c
          return
        else
c
c  Form sum
c
          a = j1 - 1
          b = j2 - 1
          h = j3 - 1
          d = m1
          e = m2
          f = m3
          s = one
          q = n - 1
          do j=1,n
            t = (a-q)/(d+q)*(b-q)/(e+q)*(h-q)/(f+q)
            s = one - s*t
            q = q - one
          end do
          clebg = c*s
          return
        endif
      else
c
c  Special formula for m3 = 0 and m1 = 0 or 1/2
c
        k = i10/2
        if(i10 .eq. 2*k) then
          k = 0
        else
          k = 1
        endif
c
        if(m1 .eq. 0) then
          l = 0
          if(k .ne. 0) return
        else if(m1 .eq. 1) then
          l = 1
        endif
c
        x  = l
        m  = i3 + (i1 + k + 1)/2 - l
        m1 = i10/2 + k
        m2 = i4 + i5
        m3 = i6 + i7
        j1 = (i1 + 1 - k    )/2
        j2 = (i2 + 1 + k - l)/2
        j3 = (i3 + 1 + k - l)/2
        clebg = DEXP((g(i11)-g(i11-1)+g(i1)+g(i2)+g(i3)-g(i10))/two +
     *          g(m1)-g(j1)-g(j2)-g(j3)+
     *          x*(g(3)-(g(m2)-g(m2-1)+g(m3)-g(m3-1))/two))
        if((m - 2*(m/2)) .ne. 0) clebg = -clebg
        return
      endif
c
c  Formats
c
 9010 format(' ','FACTORIAL TABLE SIZE IN CLEBG ',2i5,6f5.1)
c
      end
      subroutine DEPCAL(i,el,depth)
c***********************************************************************
c  Calculate the depth of the I-th potential                           *
c----------------------------------------------------------------------*
c  I     = type of the potential                                       *
c  EL    = lab energy                                                  *
c  DEPTH = depth of the potential                                      *
c----------------------------------------------------------------------*
c  Type of potential:  1  real      volume                             *
c                      2  real      surface                            *
c                      3  imaginary volume                             *
c                      4  imaginary surface                            *
c                      5  real      spin-orbit                         *
c                      6  imaginary spin-orbit                         *
c----------------------------------------------------------------------*
c  Called by SCAT                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
c
      integer          i
      double precision el,depth,dept1,dept2,u
c
      integer          nz,n1,n2
c
      double precision zero
      data             zero /0.0d+00/
c=======================================================================
c
      do nz=1,npzen(i)
        depth = zero
        if(el .lt. epot(i,nz)) then
          if(pot(i,nz,7) .eq. zero) then
            depth = pot(i,nz,1)
     *            + pot(i,nz,2)*el
     *            + pot(i,nz,3)*el*el
     *            + pot(i,nz,4)*el*el*el
     *            + pot(i,nz,5)*DLOG(el)
     *            + pot(i,nz,6)*DSQRT(el)
          else
            if(i .eq. 1) then
              u     = el - eferm
              depth = DEXP( pot(1,nz,2)*u )
              depth = pot(1,nz,7)*depth
            elseif(i .eq. 3) then
              u     = el - eferm
              n1    = NINT(pot(3,nz,3))
              depth = (u**n1) / ( u**n1 + pot(3,nz,2)**n1 )
              depth = pot(3,nz,1)*depth
            elseif(i .eq. 4) then
              u     = el - eferm
              n1    = NINT(pot(4,nz,3))
              dept1 = (u**n1) / ( u**n1 + pot(4,nz,2)**n1 )
              dept1 = pot(4,nz,1)*dept1
              n2    = NINT(pot(4,nz,6))
              dept2 = (u**n2) / ( u**n2 + pot(4,nz,5)**n2 )
              dept2 = pot(4,nz,4)*dept2
              depth = dept1 + dept2
            elseif(i .eq. 5) then
              u     = el - eferm
              depth = DEXP( pot(5,nz,2)*u )
              depth = pot(5,nz,7)*depth
            elseif(i .eq. 6) then
              u     = el - eferm
              n1    = NINT(pot(6,nz,3))
              depth = (u**n1) / ( u**n1 + pot(6,nz,2)**n1 )
              depth = pot(6,nz,1)*depth
            endif
          endif
          return
        endif
      end do
c
      return
      end
      subroutine DISPER(idr,av,rv,pote,el)
c***********************************************************************
c  Calculate dispersion relations                                      *
c----------------------------------------------------------------------*
c  IDR  = 2  use real volume + surface  potential                      *
c         1  use equivalent real volume potential                      *
c  AV   =                                                              *
c  RV   =                                                              *
c  POTE =                                                              *
c  EL   = lab energy                                                   *
c----------------------------------------------------------------------*
c  Called by SCAT                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
      INCLUDE 'poten4.i'
c
      double precision av,rv,pote,el
      integer          idr
c
      double precision sqradm,sqrhf,sqrsur,sqrvol
      double precision volhf,volsur,voltot,volvol
c
      dimension        av(nptypx),rv(nptypx),pote(nptypx)
c
      double precision zero
      data             zero  /0.0d+00/
c=======================================================================
c
c-----------------------------------------------------------------------
c  Imaginary volume  corrective term
c-----------------------------------------------------------------------
c
      if(pot(3,1,7) .eq. zero) then
        call DISPIN1(3,el,deltvv)
      else
        call DISPIN2(3,el,deltvv)
      endif
c
c-----------------------------------------------------------------------
c  Imaginary surface corrective term
c-----------------------------------------------------------------------
c
      if(pot(4,1,7) .eq. zero) then
        call DISPIN1(4,el,deltvs)
      else
        call DISPIN2(4,el,deltvs)
      endif
c
      if(idr .eq. 2) return
c
c=======================================================================
c  Real potential = equivalent volume term
c=======================================================================
c
c  Volume integrals
c
      call VOLIN(1,pote(1),rv(1),av(1),volhf )
      volhf  = volhf  / mt
      call VOLIN(1,deltvv ,rv(3),av(3),volvol)
      volvol = volvol / mt
      call VOLIN(2,deltvs ,rv(4),av(4),volsur)
      volsur = volsur / mt
      voltot = volhf + volvol + volsur
c
c  Mean-square radius
c
      call RMS(1,rv(1),av(1),sqrhf)
      call RMS(1,rv(3),av(3),sqrvol)
      call RMS(2,rv(4),av(4),sqrsur)
      sqradm = sqrhf*volhf + sqrvol*volvol + sqrsur*volsur
      sqradm = sqradm / voltot
c
c  Equivalent real radius
c
      call RMS(-1,rv(1),av(1),sqradm)
c
c  Equivalent real depth
c
      call VOLIN(-1,pote(1),rv(1),av(1),voltot*mt)
c
      return
      end
      subroutine DISPIN2(i,el,deltv)
c***********************************************************************
c  Calculate dispersion integrals for potentials whose energy          *
c  dependence is given by BRJM form                                    *
c----------------------------------------------------------------------*
c  I     = type of potential                                           *
c               3  imaginary volume                                    *
c               4  imaginary surface                                   *
c  EL    = lab energy                                                  *
c  DELTV = corrective term                                             *
c----------------------------------------------------------------------*
c  Called by DISPER                                                    *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const2.i'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
c
      double precision el,deltv
      integer          i
c
      double precision denom
      double precision f
      double precision som2,som4,som6
      double precision term2,term4,term6
c
      double precision zero,one,two,three
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
      data             three /3.0d+00/
c=======================================================================
c
      deltv   = zero
      f       = el - eferm
c
      som2    = zero
      som4    = zero
      som6    = zero
c
c  First  component
      if(IDNINT(pot(i,1,3)) .eq. 2) then
        denom =  f**2 + pot(i,1,2)**2
        term2 = one / denom
        term2 = pot(i,1,1)*f*pot(i,1,2)*term2
        som2  = som2 + term2
      else if(IDNINT(pot(i,1,3)) .eq. 4) then
        denom =  f**4 + pot(i,1,2)**4
        term4 = (f**2 + pot(i,1,2)**2) / denom
        term4 = pot(i,1,1)*f*pot(i,1,2)*term4
        term4 = term4 / DSQRT(two)
        som4  = som4 + term4
      else if(IDNINT(pot(i,1,3)) .eq. 6) then
        denom =  f**6 + pot(i,1,2)**6
        term6 = (two*(f**4 + pot(i,1,2)**4) + (f**2)*(pot(i,1,2)**2)) /
     *           denom
        term6 = pot(i,1,1)*f*pot(i,1,2)*term6
        term6 = term6/three
      endif
c
c  Second component
      if(IDNINT(pot(i,1,6)) .eq. 2) then
        denom =  f**2 + pot(i,1,2)**2
        term2 = one / denom
        term2 = pot(i,1,4)*f*pot(i,1,5)*term2
        som2  = som2 + term2
      else if(IDNINT(pot(i,1,6)) .eq. 4) then
        denom =  f**4 + pot(i,1,5)**4
        term4 = (f**2 + pot(i,1,5)**2) / denom
        term4 = pot(i,1,4)*f*pot(i,1,5)*term4
        term4 = term4 / DSQRT(two)
        som4  = som4 + term4
      else if(IDNINT(pot(i,1,6)) .eq. 6) then
        denom =  f**6 + pot(i,1,5)**6
        term6 = (two*(f**4 + pot(i,1,5)**4) + (f**2)*(pot(i,1,5)**2)) /
     *           denom
        term6 = pot(i,1,4)*f*pot(i,1,5)*term6
        term6 = term6/three
      endif
c
      deltv   = som2 + som4 + som6
c
      return
      end
      subroutine FCT(idim,fact)
c***********************************************************************
c  Calculate factorial logarithms from 0! ( =1.) up to (idim-1)!       *
c                                                                      *
c                  k! = fact(k+1)                                      *
c***********************************************************************
c
      implicit none
c
      double precision fact
      integer          idim
      dimension        fact(idim)
c
      integer          k
c
      double precision zero
      data             zero /0.0d+00/
c=======================================================================
c
      fact(1) = zero
      fact(2) = zero
c
      do k=3,idim
        fact(k) = fact(k-1) + DLOG(DBLE(k-1))
      end do
c
      return
      end
      subroutine FILES(iopt)
c***********************************************************************
c  Open and close files                                                *
c----------------------------------------------------------------------*
c  IOPT = 1   open  files                                              *
c         0   close files                                              *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'inout.i'
c
      integer          iopt
c
      character        file*80
c=======================================================================
c
      ie  = 15
      ie1 = 14
      is  = 16
      is1 = 17
      is2 = 18
      is3 = 19
c
c  Close files
c
      if(iopt .eq. 0) then
        ENDFILE is
        ENDFILE is1
        ENDFILE is2
        ENDFILE is3
c
        CLOSE(unit=is )
        CLOSE(unit=is1)
        CLOSE(unit=is2)
        CLOSE(unit=is3)
        return
      endif
c
c  Open  files
c
      OPEN(unit=ie1,file='scat2000.dat',status='old')
      OPEN(unit=is4,file='scat2000.log')
c
      read(ie1,fmt='(a)') file
      OPEN(unit=ie ,file=file,status='old')
c
      read(ie1,fmt='(a)') file
      OPEN(unit=is ,file=file)
c
      read(ie1,fmt='(a)') file
      OPEN(unit=is1,file=file)
c
      read(ie1,fmt='(a)') file
      OPEN(unit=is2,file=file)
c
      read(ie1,fmt='(a)') file
      OPEN(unit=is3,file=file)
c
      return
      end
      subroutine INTEG(npt,h,sl,spo)
c***********************************************************************
c  Integrate the Schroedinger equation                                 *
c----------------------------------------------------------------------*
c  See: A.C.Allison                                                    *
c       Journal of Computational Physics  6(1970)378-391               *
c----------------------------------------------------------------------*
c  Called by SCAT                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'inout.i'
      INCLUDE 'potn.i'
      INCLUDE 'psi.i'
c
      double precision h,sl,spo
      integer          npt
c
      double precision det,f1,f2,f3,h212
      double precision r,r2,s1,s2
      double precision u,w,y
      integer          i,j,n,npt3,npt5,m
c
      dimension        y(2,305),u(305),w(305)
      dimension        f1(2,2),f2(2,2),f3(2,2),s1(2),s2(2)
c
      double precision zero,eps,one,two,ten
      data             zero /0.0d+00/
      data             eps  /1.0d-20/
      data             one  /1.0d+00/
      data             two  /2.0d+00/
      data             ten  /1.0d+01/
c=======================================================================
c
c=======================================================================
c  Integration formula
c
c               2         -1           2               2
c     y   = (i+h *f   /12)  *((2*i-10*h *f /12)*y -(i+h *f   /12)*y   ))
c      n+1         n+1                    n      n        n-1      n-1
c
c     f1 = f(n-1) * h**2/12
c     f2 = f( n ) * h**2/12
c     f3 = f(n+1) * h**2/12
c     s1 = i + f1
c     s2 = 2*i - 10*f2
c     det = determinant ( I + f3 )
c=======================================================================
c
      h212 = h*h/1.2d+01
      npt3 = npt + 3
      npt5 = npt + 5
c
      do n=1,npt3
        u(n) = potr(n) + spo*vso(n)
        w(n) = poti(n) + spo*wso(n)
      end do
c
c  y(1,n) =   real    part of the wave function
c  y(2,n) = imaginary part of the wave function
c
c  Initial conditions
c
      y(1,1)  =  zero
      y(2,1)  =  zero
      y(1,2)  =  eps
      y(2,2)  =  eps
c
      f1(1,1) =  zero
      f1(1,2) =  zero
      f1(2,1) =  zero
      f1(2,2) =  zero
c
      f2(1,1) = -one
      f2(1,2) =  zero
      f2(2,1) =  zero
      f2(2,2) = -one
c
c  Start integration
c
      do n=3,npt5
        m  = n - 2
c       r  = DFLOAT(m)*h
        r  = DBLE  (m)*h
        r2 = r*r
c
        f3(1,1) = (ak2 - sl/r2 - u(m))*h212
        f3(1,2) = (            + w(m))*h212
        f3(2,1) = (            - w(m))*h212
        f3(2,2) = (ak2 - sl/r2 - u(m))*h212
c
        s1(1)   = (one + f1(1,1))*y(1,n-2) + f1(1,2)*y(2,n-2)
        s1(2)   = (one + f1(2,2))*y(2,n-2) + f1(2,1)*y(1,n-2)
c
        s2(1)   = (two - ten*f2(1,1))*y(1,n-1) - ten*f2(1,2)*y(2,n-1)
        s2(2)   = (two - ten*f2(2,2))*y(2,n-1) - ten*f2(2,1)*y(1,n-1)
c
        det     = (f3(1,1)+one)*(f3(2,2)+one) - f3(1,2)*f3(2,1)
c
        y(1,n)  = (f3(2,2)+one)*(s2(1)-s1(1)) - f3(1,2)*(s2(2)-s1(2))
        y(2,n)  = (f3(1,1)+one)*(s2(2)-s1(2)) - f3(2,1)*(s2(1)-s1(1))
        y(1,n)  =  y(1,n) / det
        y(2,n)  =  y(2,n) / det
c
        do i=1,2
          do j=1,2
            f1(i,j) = f2(i,j)
            f2(i,j) = f3(i,j)
          end do
        end do
c
      end do
c
c  Calculation of the derivatives
c
      n     = npt + 2
      psir  = y(1,n)
      psii  = y(2,n)
      psirp = (y(1,n+3) - y(1,n-3) + 9.*(y(1,n-2) - y(1,n+2))
     *        + 45.*(y(1,n+1) - y(1,n-1)))/ (60.*h)
      psiip = (y(2,n+3) - y(2,n-3) + 9.*(y(2,n-2) - y(2,n+2))
     *        + 45.*(y(2,n+1) - y(2,n-1)))/ (60.*h)
c
      return
      end
      subroutine KINEMA(el,e1,mi,mt,amu,ak2,iopt,irelat)
c***********************************************************************
c  Kinematics:   lab  <===>  CM                                        *
c    With relativistic kinematics, the reduced mass is replaced by     *
c    the reduced total energy                                          *
c----------------------------------------------------------------------*
c  EL     = current lab kinetic energy                                 *
c  E1     = current  CM kinetic energy                                 *
c  MI     = incident particle rest mass (in a.m.u.)                    *
c  MT     = target   nucleus  rest mass (in a.m.u.)                    *
c  AMU    = reduced mass                                               *
c  AK2    = CM wave number                                             *
c  IOPT   = 1   from lab to CM                                         *
c           2   from CM  to lab                                        *
c  IRELAT = 0   classical    kinematics                                *
c           1   relativistic kinematics                                *
c----------------------------------------------------------------------*
c  AMU0C2 = a.m.u. in MeV                                              *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'inout.i'
      INCLUDE 'physcon.i'
c
      double precision el,e1,mi,mt,amu,ak2
      integer          iopt,irelat
c
      double precision ck2,mtot,p2,w2
      double precision etoti,etott
c
      double precision one,two
      data             one  /1.0d+00/
      data             two  /2.0d+00/
c=======================================================================
c
      ck2  = (two*amu0c2)/(hbarc**2)
c
      mtot = mi + mt
c
      if(iopt .eq. 1) then
c
c***********************************************************************
c  From lab to CM (the input quantity is Elab)
c***********************************************************************
c
        if(irelat .eq. 0) then
c
c-----------------------------------------------------------------------
c  Classical    kinematics
c-----------------------------------------------------------------------
c
          amu = mi*mt / mtot
          e1  = el*mt / mtot
          w2  = ck2*amu
          ak2 = w2*e1
        else if(irelat .eq. 1) then
c
c-----------------------------------------------------------------------
c  Relativistic kinematics
c-----------------------------------------------------------------------
c
          e1  = amu0c2*mtot*
     *         (DSQRT(one + two*el/(amu0c2*mt*((one+mi/mt)**2))) - one)
          p2  = (el*(el + two*amu0c2*mi)) /
     *          ((one+mi/mt)**2 + two*el/(amu0c2*mt))
          ak2 = p2 / (hbarc*hbarc)
          etoti = DSQRT((amu0c2*mi)**2 + p2)
          etott = DSQRT((amu0c2*mt)**2 + p2)
          amu   = etoti*etott / (etoti + etott)
          amu   = amu / amu0c2
        else
          stop 9873
        endif
      else if(iopt .eq. 2) then
c
c***********************************************************************
c  From  CM to lab (the input quantity is Ecm)
c***********************************************************************
c
        if(irelat .eq. 0) then
c
c-----------------------------------------------------------------------
c  Classical    kinematics
c-----------------------------------------------------------------------
c
          amu = mi*mt / mtot
          el  = e1*mtot / mt
          ak2 = ck2*amu*e1
        else if(irelat .eq. 1) then
c
c-----------------------------------------------------------------------
c  Relativistic kinematics
c-----------------------------------------------------------------------
c
          el  = e1*(e1+two*amu0c2*mtot) / (two*amu0c2*mt)
          p2  = e1*(e1+two*amu0c2*mtot)*(e1+two*amu0c2*mi)*
     *             (e1+two*amu0c2*mt  ) / ( (two*(e1+amu0c2*mtot))**2 )
          ak2 = p2 / (hbarc*hbarc)
          etoti = DSQRT((amu0c2*mi)**2 + p2)
          etott = DSQRT((amu0c2*mt)**2 + p2)
          amu   = etoti*etott / (etoti + etott)
          amu   = amu / amu0c2
        else
          stop 9873
        endif
      else
        stop 9872
      endif
c
      return
      end
      subroutine PREANG(ida,ip,na)
c***********************************************************************
c  Calculate Legendre and associated Legendre polynomials              *
c----------------------------------------------------------------------*
c  If IDA = 1 equidistant angles between 0 and 180 degrees             *
c             by step of 2.5                                           *
c         = 2 angles whose cosines are equidistant between             *
c             -1 and 1 by step of 0.02                                 *
c  IP = type of incident particle (7 = heavy ion)                      *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'angles1.i'
      INCLUDE 'angles2.i'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'inout.i'
      INCLUDE 'poleg.i'
c
      integer          ida,ip,na
c
      double precision d,x,xl
      double precision rad
      integer          i,idaa,l
c
      double precision zero,eps,one,two,three,deg
      data             zero  /0.0d+00/
      data             eps   /1.0d-06/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
      data             three /3.0d+00/
      data             deg   /1.8d+02/
c=======================================================================
c
      rad   = deg / pi
c
      idaa  = IABS(ida)
c
c  Calculate cosines of equally spaced angles
c
      if(idaa .eq. 1) then
        if(ip .ne. 7) then
          na = 73
          d  = deg / DBLE(na-1)
        else
          na = tetaf / tetastep
          if(na .le. nangmx) then
            d  = tetastep
          else
            na = nangmx
            d  = tetaf / DBLE(na-1)
          endif
        endif
c
        do i=1,na
           ang(i) = DBLE(i-1)*d
          cang(i) = DCOS(ang(i)/rad)
          if(DABS(cang(i)) .lt. eps) cang(i) = zero
        end do
c
c  Calculate equally spaced cosines
c
      else if(idaa .eq. 2) then
        na = 101
        d  = two / DBLE(na-1)
        do i=1,na
          cang(i) = one - DBLE(i-1)*d
          if(DABS(cang(i)) .lt. eps) cang(i) = zero
        end do
      endif
c
c  Calculate Legendre polynomials
c
      do i=1,na
        x        = cang(i)
        pl(1,i)  = one
        pl(2,i)  = x
        do l=3,ltcmax
          xl      = DBLE(l-1)
          pl(l,i) = ((two*xl-one)*x*pl(l-1,i)-(xl-one)*pl(l-2,i))
     *             /xl
        end do
      end do
c
c  Calculate associated Legendre polynomials
c
      do i=1,na
        x        = cang(i)
        pl1(1,i) = zero
        pl1(2,i) = DSQRT(one - x*x)
        do l=3,ltcmax
          xl       = DBLE(l-1)
          pl1(l,i) = ((two*xl-one)*x*pl1(l-1,i)-xl*pl1(l-2,i))
     *              /(xl-one)
        end do
      end do
c
      do i=1,na
        x        = cang(i)
        pl2(1,i) = zero
        pl2(2,i) = zero
        pl2(3,i) = three*(one - x*x)
        do l=4,ltcmax
          xl       = DBLE(l-1)
          pl2(l,i) = ((two*xl-one)*x*pl2(l-1,i) - (xl+one)*pl2(l-2,i))
     *              /(xl-two)
        end do
      end do
c
      return
      end
      subroutine PRIPOT(idr,irelat,ref)
c***********************************************************************
c  Print optical potential parameters                                  *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
c
      character        ref*(*)
      integer          idr,irelat
c
      character*10     type(nptypx)
      integer          i,j,k
c=======================================================================
c
      type(1) = 'REAL VOL. '
      type(2) = 'REAL SURF.'
      type(3) = 'IMAG.VOL. '
      type(4) = 'IMAG.SURF.'
      type(5) = 'REAL  S.O.'
      type(6) = 'IMAG. S.O.'
c
      write(is,9000)
      write(is,9010)
      write(is,9020) zi,mi,zt,mt
      write(is,9030) ref
      write(is,9034)
      write(is,9035)
      write(is,9040)
      do i=1,nptypx
        do j=1,npzen(i)
          write(is,9100) type(i),epot(i,j-1),epot(i,j),
     *                  (pot(i,j,k),k=1,npcofx),
     *                   r(i),re(i),a(i),ae(i)
        end do
      end do
      write(is,9050)
      write(is,9200) rcoul,eferm,beta
c
      if(irelat .eq. 0) then
        write(is,9210) 'Kinematics: non-relativistic'
      elseif(irelat .eq. 1) then
        write(is,9210) 'Kinematics:     relativistic'
      endif
c
      if(idr .eq. 0) then
        write(is,9050)
        write(is,9300)
      else if(idr .eq. 1) then
        write(is,9050)
        write(is,9310)
      else if(idr .eq. 2) then
        write(is,9050)
        write(is,9320)
        write(is,9050)
      endif
c
c  Formats
c
 9000 format('1',14x,'TRANSMISSION COEFFICIENTS CALCULATED FROM THE'
     *          ,    ' FOLLOWING OPTICAL MODEL PARAMETERS')
 9010 format(' ',14x,80('*'),///)
 9020 format(' ',34x,'CHARGE',29x,'MASS',//
     *          ,15x,'PROJECTILE',10x,   f6.1,20x,1p,e13.6,/
     *          ,15x,'TARGET    ',10x,0p,f6.1,20x,1p,e13.6,//)
 9030 format(' ','PARAMETER SET : ',a)
 9034 format(' ',57x,'STRENGTH',39x,'RADIUS',14x,'DIFFUSENESS')
 9035 format(' ',15x,'Emin',6x,'Emax',6x,'Cst',9x
     *          ,'E',8x,'E*E',6x,'E*E*E',5x,'LN(E)',5x,'SR(E)'
     *          ,4x,'EXP(E)',17x,'E',19x,'E')
 9040 format(' ',14x,16('-'),3x,67('-'),4x,16('-'),4x,16('-'))
 9050 format(' ',56x)
 9100 format(' ',a10,13(2x,f8.4))
 9200 format(' ','COUL.RAD.',4x,f7.4,14x
     *          ,'E-Fermi =',   f9.4, 2x
     *          ,'BETA =   ',   f7.4)
 9210 format(' ',a)
 9300 format(' ','DISP.REL.: NO')
 9310 format(' ','DISP.REL.: REAL POT. = EQUIVALENT  VOL.TERM')
 9320 format(' ','DISP.REL.: REAL POT. = VOL.TERM + SURF.TERM')
 9400 format(' ','E-FERMI   ',2x,f9.4)
c
      return
      end
      subroutine PRITC(ne)
c***********************************************************************
c  Print a summary table of: - the J-averaged transmission coefficients*
c                            - the compound nucleus cross sections     *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'ener.i'
      INCLUDE 'inout.i'
      INCLUDE 'tce.i'
      INCLUDE 'xs.i'
c
      integer          ne
c
      integer          k,kmax,l,lmax,lmin,lpmax,m
c=======================================================================
c
      write(is,9010) mi,zi,mt,zt
c
c  Find LMAX for printing
c
      do l=1,ltcmax
        if(tc(ne,l) .lt. 5.0d-07) go to 2
      end do
c
    2 lmax = l - 1
      kmax = lmax/10 + 1
c
      do k=1,kmax
        lmin  = 10*(k-1)
        lpmax = 10* k    - 1
        lpmax = MIN0(lpmax,lmax)
        write(is,9020) (l,l=lmin,lpmax)
        write(is,9030)
        lmin  = lmin  + 1
        lpmax = lpmax + 1
        do m=1,ne
          write(is,9040) ein(m),sr(m),(tc(m,l),l=lmin,lpmax)
        end do
        write(is,9060)
      end do
c
c  Formats
c
 9010 format('1',29x,'PROJECTILE  A=',f4.0,3x,'Z=',f4.0,6x,
     *               '    TARGET  A=',f5.0,3x,'Z=',f5.0,//)
 9020 format(' ','  ENERGY',3x,'SIGMA-CN',1x,11(5x,i2,3x))
 9030 format(56x)
 9040 format(' ',f8.4,2x,f9.2,1x,11(2x,f8.6))
 9060 format(/)
c
      return
      end
      double precision function RACAH(a,b,c,d,e,f)
c***********************************************************************
c  Calculate Racah coefficients      w(a,b,c,d;e,f)                    *
c                                                                      *
c  Attention:   w(a,b,c,d;e,f) = (-1)**(a+b+c+d)*6-J(a,b,e;d,c,f)      *
c                                                ---                   *
c                                                                      *
c  from John.G.Wills     ORNL-TM-1949 (August 1967)                    *
c                    and Comp.Phys.Comm. 2(1971)381                    *
c                                                                      *
c  O.Bersillon     August 1977                                         *
c                                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'fact.i'
c
      double precision a,b,c,d,e,f
c
      double precision h,o,p,q,r,s,t,v,w,x,y,z
      integer          i
      integer          i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11,i12,i13
     *                ,i14,i15,i16,il
      integer          ja,jb,jc,jd,je,jf
      integer          j1,j2,j3,j4,j5,j6,j7
      integer          j,k,n
c
      dimension        i(16)
      equivalence
     *(i( 1), i1),(i( 2), i2),(i( 3), i3),(i( 4), i4),(i( 5), i5),
     *(i( 6), i6),(i( 7), i7),(i( 8), i8),(i( 9), i9),(i(10),i10),
     *(i(11),i11),(i(12),i12),(i(13),i13),(i(14),i14),(i(15),i15),
     *(i(16),i16)
c
      double precision zero,eps,one,two
      data             zero /0.0d+00/
      data             eps  /1.0d-03/
      data             one  /1.0d+00/
      data             two  /2.0d+00/
c=======================================================================
c
      racah = zero
c
c  Convert arguments to integer and make useful combinations
c
      ja  = IDINT(two*a + eps)
      jb  = IDINT(two*b + eps)
      jc  = IDINT(two*c + eps)
      jd  = IDINT(two*d + eps)
      je  = IDINT(two*e + eps)
      jf  = IDINT(two*f + eps)
      i1  = ja + jb - je
      i2  = jb + je - ja
      i3  = je + ja - jb
      i4  = jc + jd - je
      i5  = jd + je - jc
      i6  = je + jc - jd
      i7  = ja + jc - jf
      i8  = jc + jf - ja
      i9  = jf + ja - jc
      i10 = jb + jd - jf
      i11 = jd + jf - jb
      i12 = jf + jb - jd
      i13 = ja + jb + je
      i14 = jc + jd + je
      i15 = ja + jc + jf
      i16 = jb + jd + jf
c
c  Check triangular inequalities, find no. of terms in sum,
c    divide I's by 2
c
      n = i16
      do j=1,12
        k = i(j)/2
        if(i(j) .ne. 2*k) return
        if(k .lt. 0) return
        if(k .lt. n) then
          n = k
        endif
        i(j) = k + 1
      end do
c
c  Find minimum value of summation index
c
      il = 0
      do j=13,16
        i(j) = i(j)/2
        if(il .lt. i(j)) then
          il = i(j)
        endif
      end do
      j1 = il  - i13 + 1
      j2 = il  - i14 + 1
      j3 = il  - i15 + 1
      j4 = il  - i16 + 1
      j5 = i13 + i4  - il
      j6 = i15 + i5  - il
      j7 = i16 + i6  - il
      h  = -      DEXP((g(i1)+g(i2)+g(i3)-g(i13+2)+g(i4)+g(i5)+g(i6)-
     *g(i14+2)+g(i7)+g(i8)+g(i9)-g(i15+2)+g(i10)+g(i11)+g(i12)-g(i16+2))
     */two+g(il+2)-g(j1)-g(j2)-g(j3)-g(j4)-g(j5)-g(j6)-g(j7))
      if((j5 - 2*(j5/2)) .ne. 0) h = -h
      if(n .lt. 0) return
      if(n .eq. 0) then
        racah = h
        return
      else
c
        s = one
        q = n  - 1
        p = il + 2
        r = j1
        o = j2
        v = j3
        w = j4
        x = j5 - 1
        y = j6 - 1
        z = j7 - 1
        do j=1,n
          t = (p+q)/(r+q)*(x-q)/(o+q)*(y-q)/(v+q)*(z-q)/(w+q)
          s = one - s*t
          q = q   - one
        end do
        racah = h*s
      endif
c
      return
      end
      subroutine RCWFN(rho,eta,minl,maxl,fc,fcp,gc,gcp,accur,step)
c***********************************************************************
c  Coulomb wave functions calculated at r = rho                        *
c  by the continued-fraction method of J.W.Steed                       *
c  MINL, MAXL are actual l-values                                      *
c----------------------------------------------------------------------*
c  See: A.R.Barnett, D.H.Feng, J.W.Steed and L.J.B.Goldfarb            *
c       Computer Physics Communications  8(1974)377-395                *
c----------------------------------------------------------------------*
c  Called by SCAT                                                      *
c***********************************************************************
      implicit double precision (a-h,o-z)
c
      INCLUDE 'scat2000.parameter'
c
      double precision k,k1,k2,k3,k4,m1,m2,m3,m4
c
      dimension        fc(ltcmax),fcp(ltcmax),gc(ltcmax),gcp(ltcmax)
c=======================================================================
c
      zero  = 0.0d+00
      one   = 1.0d+00
      two   = 2.0d+00
      gpmax = 1.0d+75
c
      pace  = step
      acc   = accur
      if(pace .lt. 1.0d+02) pace = 1.0d+02
      if(acc  .lt. 1.0d-15 .or. acc .gt. 1.0d-06) acc = 1.0d-06
      r     = rho
      ktr   = 1
      lmax  = maxl
      lmin1 = minl+1
c90   xll1  = DFLOAT(minl*lmin1)
      xll1  = DBLE  (minl*lmin1)
      eta2  = eta*eta
      turn  = eta + DSQRT(eta2 + xll1)
      if(r .lt. turn .and. DABS(eta) .ge. 1.0d-06) ktr=-1
      ktrp  = ktr
      go to 2
    1 r     = turn
      tf    = f
      tfp   = fp
      lmax  = minl
      ktrp  = 1
    2 etar  = eta*r
      rho2  = r*r
c90   pl    = DFLOAT(lmax+1)
      pl    = DBLE  (lmax+1)
      pmx   = pl + 0.5d+00
c
c  Continued fraction for fp(maxl)/f(maxl), xl is f, xlprime is fp
c
      fp    = eta/pl + pl/r
      dk    = etar*two
      del   = zero
      d     = zero
      f     = one
      k     = (pl*pl - pl + etar)*(two*pl - one)
      if(pl*pl + pl + etar .ne. zero) go to 3
      r     = r + 1.0d-06
      go to 2
    3 h     = (pl*pl + eta2)*(one - pl*pl)*rho2
      k     = k + dk + pl*pl*6.0d+00
      d     = one/(d*h+k)
      del   = del*(d*k - one)
      if(pl .lt. pmx) del = -r*(pl*pl + eta2)*(pl + one)*d/pl
      pl    = pl + one
      fp    = fp + del
      if(d .lt. zero) f = -f
      if(pl .gt. 2.0d+04) go to 11
      if(DABS(del/fp) .ge. acc) go to 3
      fp    = f*fp
      if(lmax .eq. minl) go to 5
      fc(lmax+1)  = f
      fcp(lmax+1) = fp
c
c  Downward recursion to minl for f and fp, arrays gc, gcp are storage
c
      l     = lmax
      do 4 lp = lmin1,lmax
c90   pl      = DFLOAT(l)
      pl      = DBLE  (l)
      gc(l+1) = eta/pl + pl/r
      gcp(l+1)= DSQRT(eta2 + pl*pl)/pl
      fc(l)   = (gc(l+1)*fc(l+1) + fcp(l+1))/gcp(l+1)
      fcp(l)  =  gc(l+1)*fc(l)   - gcp(l+1)*fc(l+1)
    4 l     = l - 1
      f     = fc(lmin1)
      fp    = fcp(lmin1)
    5 if(ktrp .eq. -1) go to 1
c
c  Repeat for r = turn if rho lt turn
c  Now obtain p + i.q for minl from continued fraction (32)
c  Real arithmetic to facilitate conversion to ibm using double precision
c
      p     = zero
      q     = r - eta
      pl    = zero
      ar    = -(eta2 + xll1)
      ai    = eta
      br    = two*q
      bi    = two
      wi    = two*eta
      dr    =  br/(br*br + bi*bi)
      di    = -bi/(br*br + bi*bi)
      dp    = -(ar*di + ai*dr)
      dq    =   ar*dr - ai*di
    6 p     = p  + dp
      q     = q  + dq
      pl    = pl + two
      ar    = ar + pl
      ai    = ai + wi
      bi    = bi + two
      d     = ar*dr - ai*di + br
      di    = ai*dr + ar*di + bi
      t     = one/(d*d + di*di)
      dr    = d*t
      di    = -t*di
      h     = br*dr - bi*di - one
      k     = bi*dr + br*di
      t     = dp*h  - dq*k
      dq    = dp*k  + dq*h
      dp    = t
      if(pl .gt. 4.6d+04) go to 11
      if(DABS(dp)+DABS(dq) .ge. (DABS(p)+DABS(q))*acc) go to 6
      p     = p/r
      q     = q/r
c
c  Solve for fp, g, gp and normalise f at l = minl
c
      g     = (fp - p*f)/q
      gp    = p*g - q*f
      w     = one/DSQRT(fp*g - f*gp)
      g     = w*g
      gp    = w*gp
      if(ktr .eq. 1) go to 8
      f     = tf
      fp    = tfp
      lmax  = maxl
c
c  Runge-Kutta integration of g(minl) and gp(minl) inwards from turn
c
      if(rho .lt. 0.2d+00*turn) pace=999.0
      r3    = one/3.0d+00
      h     = (rho - turn)/(pace + one)
      h2    = 0.5d+00*h
      i2    = IDINT(pace + 1.0d-03)
      etah  = eta*h
      h2ll  = h2*xll1
      s     = (etah + h2ll/r)/r - h2
    7 rh2   = r + h2
      t     = (etah + h2ll/rh2)/rh2 - h2
      k1    = h2*gp
      m1    = s*g
      k2    = h2*(gp + m1)
      m2    = t *(g  + k1)
      k3    = h *(gp + m2)
      m3    = t *(g  + k2)
      m3    = m3 + m3
      k4    = h2*(gp + m3)
      rh    = r + h
      s     = (etah + h2ll/rh)/rh - h2
      m4    = s*(g + k3)
      g     = g  + (k1 + k2 + k2 + k3 + k4)*r3
      gp    = gp + (m1 + m2 + m2 + m3 + m4)*r3
      r     = rh
      i2    = i2 - 1
      if(DABS(gp) .gt. gpmax) go to 11
      if(i2 .ge. 0) go to 7
      w     = one/(fp*g-f*gp)
c
c  Upward recursion from gc(minl) and gcp(minl), stored values are r,
c  Renormalise fc, fcp for each l-value
c
    8 gc(lmin1)  = g
      gcp(lmin1) = gp
      if(lmax .eq. minl) go to 10
      do 9 l=lmin1,lmax
      t       = gc(l+1)
      gc(l+1) = (gc(l)*gc(l+1) - gcp(l))/gcp(l+1)
      gcp(l+1)= gc(l)*gcp(l+1) - gc(l+1)*t
      fc(l+1) = w*fc(l+1)
    9 fcp(l+1)= w*fcp(l+1)
      fc(lmin1)  = fc(lmin1)*w
      fcp(lmin1) = fcp(lmin1)*w
      go to 12
   10 fc(lmin1)  = w*f
      fcp(lmin1) = w*fp
      go to 12
   11 w     = zero
      g     = zero
      gp    = zero
      go to 8
   12 return
      end
      subroutine RMS(itype,rad,dif,sqradm)
c***********************************************************************
c  Calculate the mean-square radius for the ITYPE potential            *
c----------------------------------------------------------------------*
c  ITYPE  = 1 or -1            Woods-Saxon potential                   *
c           2 or -2 derivative Woods-Saxon potential                   *
c  ITYPE  > 0       calculate SQRADM from RAD and DIF                  *
c         < 0       calculate RAD from SQRADM and DIF (iteratively)    *
c  RAD    = radius      of the potential                               *
c  DIF    = diffuseness of the potential                               *
c  SQRADM = mean-square radius                                         *
c----------------------------------------------------------------------*
c  Called by DISPER                                                    *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'const2.i'
c
      double precision rad,dif,sqradm
      integer          itype
c
      double precision df,dr,fn,fprev,r,sign
c
      double precision eps,one,two
      data             eps   /1.0d-06/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
c
      external         RMSSUR,RMSVOL
      double precision RMSSUR,RMSVOL
c=======================================================================
c
c
c-----------------------------------------------------------------------
c  Volume  potential
c-----------------------------------------------------------------------
c
      if(itype .eq. 1) then
        sqradm = RMSVOL(rad,dif)
      else if(itype .eq. -1) then
        fprev = RMSVOL(rad,dif)
        sign  = one
        r     = rad
        dr    = one
        df    = one
        do while (df .gt. eps)
          r  = r + sign*dr
          fn = RMSVOL(r,dif)
          df = DABS((fn - fprev) / fprev)
          if(fn .gt. sqradm) then
            sign = -one
          else
            sign = one
            dr   = dr / two
          endif
          fprev  = fn
        end do
        rad = r
c
c-----------------------------------------------------------------------
c  Surface potential
c-----------------------------------------------------------------------
c
      else if(itype .eq.  2) then
        sqradm = RMSSUR(rad,dif)
      else if(itype .eq. -2) then
        fprev = RMSSUR(rad,dif)
        sign  = one
        r     = rad
        dr    = one
        df    = one
        do while (df .gt. eps)
          r  = r + sign*dr
          fn = RMSSUR(r,dif)
          df = DABS((fn - fprev) / fprev)
          if(fn .gt. sqradm) then
            sign = -one
          else
            sign = one
            dr   = dr / two
          endif
          fprev  = fn
        end do
        rad = r
      endif
c
      return
      end
      double precision function RMSVOL(rad,dif)
c***********************************************************************
c  Calculate the mean-square radius for a volume  potential            *
c----------------------------------------------------------------------*
c  Called by RMS                                                       *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'const2.i'
c
      double precision rad,dif
c
      double precision arg,arg2,rad2
c
      double precision one,three,five,seven,ten
      data             one    /1.0d+00/
      data             three  /3.0d+00/
      data             five   /5.0d+00/
      data             seven  /7.0d+00/
      data             ten    /1.0d+01/
c=======================================================================
c
      arg    = (pi*dif/rad)**2
      arg2   = arg*arg
      rad2   = rad*rad
c
      rmsvol = three*rad2*(one + seven*arg/three + ten*arg2/three) /
     *         (five*(one + arg))
c
      return
      end
      subroutine SCAT(ip,lmax,idr)
c***********************************************************************
c  Calculate: - transmission coefficients T(l,j)                       *
c             - matrix elements         ETA(l,j)                       *
c----------------------------------------------------------------------*
c  E1  = current  CM energy                                            *
c  EL  = current lab energy                                            *
c----------------------------------------------------------------------*
c  Type of potential:  1  real      volume                             *
c                      2  real      surface                            *
c                      3  imaginary volume                             *
c                      4  imaginary surface                            *
c                      5  real      spin-orbit                         *
c                      6  imaginary spin-orbit                         *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'ener.i'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
      INCLUDE 'poten4.i'
      INCLUDE 'potn.i'
      INCLUDE 'psi.i'
      INCLUDE 'tlj.i'
c
      integer          ip,lmax,idr
c
      character*10     type(nptypx)
      character*10     typd(2)
      double precision a1,a2,a3,a4,a5
      double precision accur,ak,av
      double precision c1,c2
      double precision d1,depth,dt1,dt2,dt3,dt4,dt5,dt6
c     double precision el,eti,etr
      double precision    eti,etr
      double precision fc,fcp,fl,fj,fjmin,gc,gcp,h
      double precision mi13,mt13
      double precision one3
      double precision p1,p2,p3,p4,p5,pote,psim
      double precision r1,r3,r4,rap3,rap4,rcoulb,rho,rint,rm,rv
      double precision redrad,sl,spo,step
      double precision t1,t2,t3,t4,t5,t6
      double precision u1,u2,vcl,w2,y1,y2,yy,zz
      double precision volhf,volint,volsur,volvol
      integer          i,ipl,j,k,l,lm1,lmaxc,npt,npt3
      logical          idbg
c
      dimension        av(nptypx),rv(nptypx),pote(nptypx)
      dimension        fc(ltcmax),fcp(ltcmax),gc(ltcmax),gcp(ltcmax)
      dimension        u1(7),y1(7)
      dimension        volint(nptypx)
c
      double precision argmin
      data             argmin /-1.745d+02/
      double precision zero,epstl,one,two,three,four,seven
      data             zero   /0.0d+00/
      data             epstl  /1.0d-10/
      data             one    /1.0d+00/
      data             two    /2.0d+00/
      data             three  /3.0d+00/
      data             four   /4.0d+00/
      data             seven  /7.0d+00/
c=======================================================================
c
      type(1) = 'REAL VOL. '
      type(2) = 'REAL SURF.'
      type(3) = 'IMAG.VOL. '
      type(4) = 'IMAG.SURF.'
      type(5) = 'REAL  S.O.'
      type(6) = 'IMAG. S.O.'
c
      typd(1) = 'DELT V-VOL'
      typd(2) = 'DELT V-SUR'
c
      idbg    = .false.
c
      one3    =  one / three
c
      do i=1,3
        do j=1,ltcmax
          br(i,j) = zero
          bi(i,j) = zero
          t (i,j) = zero
        end do
      end do
c
c  Physical constants
c    W2  = (2*(amu*c**2)) / ((hbar*c)**2)
c    AK2 = k**2
c    VCL = constant for the Coulomb potential
c
      w2     = ck2*amu
      ak     = DSQRT(ak2)
      zz     = zi*zt
      eta    = ceta*zz*DSQRT(amu/e1)
      vcl    = ak*eta/w2
      mt13   = mt**one3
c
c  Radius
c
      if(ip .ne. 7) then
        mi13 = zero
      else
        mi13 = mi**one3
      endif
c
      do i=1,nptypx
        rv(i) = (DABS(r(i)) + re(i)*el)*(mt13+mi13)
      end do
      rcoulb  = rcoul*(mt13+mi13)
c
c  Diffuseness
c
      do i=1,nptypx
        av(i) = a(i) + ae(i)*el
        if(av(i) .le. zero) av(i) = one
      end do
c
c  Strength
c
c  Real      volume
      call DEPCAL(1,el,pote(1))
      if(intpot(1) .eq. 1) then
        call VOLIN(-1,depth,rv(1),av(1),pote(1))
        pote(1) = depth
      endif
c
c  Real      surface
      call DEPCAL(2,el,pote(2))
      if(intpot(2) .eq. 1) then
        call VOLIN(-2,depth,rv(2),av(2),pote(2))
        pote(2) = depth
      endif
c
c  Imaginary volume
      call DEPCAL(3,el,pote(3))
      if(intpot(3) .eq. 1) then
        call VOLIN(-1,depth,rv(3),av(3),pote(3))
        pote(3) = depth
      endif
c
c  Imaginary surface
      call DEPCAL(4,el,pote(4))
      if(intpot(4) .eq. 1) then
        call VOLIN(-2,depth,rv(4),av(4),pote(4))
        pote(4) = depth
      endif
c
c  Real spin-orbit
      call DEPCAL(5,el,pote(5))
c
c  Imaginary spin-orbit
      call DEPCAL(6,el,pote(6))
c
c  Dispersion relation
c
      if(idr .ne. 0) call DISPER(idr,av,rv,pote,el)
c
c  Matching radius
c
      r1  = rv(1) + seven*a(1)
      r3  = rv(3) + seven*a(3)
      r4  = rv(4) + seven*a(4)
      rm  = 1.5*DMAX1(r1,r3,r4)
      rho = ak*rm
c
c  Calculate significant l-values only
c
      lmaxc = 2*IDINT(2.5*(rho**0.8) + 3.5)
      lmaxc = MIN0(lmaxc,ltcmax)
c
c  Coulomb functions at the matching radius
c
      accur = 1.0d-14
      step  = 999.0
      call RCWFN(rho,eta,0,lmaxc-1,fc,fcp,gc,gcp,accur,step)
c
c=======================================================================
c  Calculate potentials for each integration step (n*h)
c=======================================================================
c
      npt  = 200
      npt3 = npt + 3
      h    = rm / DBLE(npt)
      if(beta .ne. zero) then
        pote(3) = zero
        c2      = beta*beta/16.0d+00
        c1      = four*w2*c2
        d1      = DEXP(four*c2*ak2)
      endif
c
      t1  = one / DEXP(rv(1)/av(1))
      t2  = one / DEXP(rv(2)/av(2))
      t3  = one / DEXP(rv(3)/av(3))
      t4  = one / DEXP(rv(4)/av(4))
      t5  = one / DEXP(rv(5)/av(5))
      t6  = one / DEXP(rv(6)/av(6))
c
      dt1 = DEXP(h/av(1))
      dt2 = DEXP(h/av(2))
      dt3 = DEXP(h/av(3))
      dt4 = DEXP(h/av(4))
      dt5 = DEXP(h/av(5))
      dt6 = DEXP(h/av(6))
c
      do i=1,npt3
c
c  Geometry
        rint = DBLE(i)*h
        t1   = t1*dt1
        t2   = t2*dt2
        t3   = t3*dt3
        t4   = t4*dt4
        t5   = t5*dt5
        t6   = t6*dt6
c
c  Potentials
        potr(i) = zero
        poti(i) = zero
c
c  Real volume  potential
        potr(i) = potr(i) + pote(1) / (one + t1)
c
c  Real surface potential      (Woods-Saxon derivative)
        potr(i) = potr(i) + 4.0*pote(2) * t2/((one + t2)**2)
c
c  Dispersion relations
        if(idr .eq. 2) then
          potr(i) = potr(i) + deltvv / (one + t3) +
     *              four*deltvs * t4/((one + t4)**2)
        endif
c
c  Imaginary volume  potential
        poti(i) = poti(i) + pote(3) / (one + t3)
c
c  Imaginary surface potential (Woods-Saxon derivative)
        if(r(4) .gt. zero) then
          poti(i) = poti(i) + 4.0*pote(4) * t4/((one + t4)**2)
c
c  Imaginary surface potential (Gaussian)
        else if(r(4) .lt. zero) then
          yy = -(((rint - rv(4))/av(4))**2)
          if(yy .gt. argmin) then
            poti(i) = poti(i) + pote(4) * DEXP(yy)
          endif
        endif
c
        if(beta .gt. zero) then
c
c  Local equivalent of a non local potential
c
          p1 = c2/(potr(i)*potr(i) + poti(i)*poti(i))
          p4 = t1/(one + t1)
          p5 = t4/(one + t4)
          p2 =-potr(i)*p4/av(1)
          p3 = poti(i)*(one - 2.*p5)/av(4)
          p4 = p2     *(one - 2.*p4)/av(1)
          p5 = poti(i)*(1.-6.*p5*(1.-p5))/(av(4)*av(4))
          u2 = p1*((potr(i)*p2 + poti(i)*p3)*2./rint +
     *         potr(i)*p4      + poti(i)*p5)
          y2 = p1*((potr(i)*p3 - poti(i)*p2)*2./rint +
     *         potr(i)*p5      - poti(i)*p4)
          y1(1) =  poti(i)/(d1*d1 + 2.*d1*c1*potr(i))
          u1(1) = (potr(i) + poti(i)*c1*y1(1))/(d1 + c1*potr(i))
          do k=1,6
            p1 = c1*y1(k) - y2
            p2 = DSIN(p1)
            p1 = DCOS(p1)
            p3 = one/(d1*DEXP(c1*u1(k) - u2))
            u1(k+1) = (potr(i)*p1 + poti(i)*p2)*p3
            y1(k+1) = (poti(i)*p1 - potr(i)*p2)*p3
          end do
c
          p1 = u1(7) - 2.*u1(6) + u1(5)
          if(p1 .eq. zero) then
            potr(i) = u1(7)
          else
            potr(i) = u1(7) - ((u1(7)-u1(6))**2)/p1
          endif
c
          p2 = y1(7) - 2.*y1(6) + y1(5)
          if(p2 .eq. zero) then
            poti(i) = y1(7)
          else
            poti(i) = y1(7) - ((y1(7) - y1(6))**2)/p2
          endif
        endif
c
c  Real + Coulomb potentials
c
        if(zz .ne. zero) then
          if(rint .gt. rcoulb) then
            potr(i) = potr(i) - 2.0*vcl/rint
          else
            potr(i) = potr(i) - vcl*(3.0 - (rint/rcoulb)**2)/rcoulb
          endif
        endif
c
c  Spin-orbit potentials
c
        vso(i)  = +cso*pote(5) * t5/(av(5)*rint*((one + t5)**2))
        wso(i)  = +cso*pote(6) * t6/(av(6)*rint*((one + t6)**2))
c
        potr(i) = -potr(i)*w2
        poti(i) = -poti(i)*w2
        vso(i)  = -vso(i) *w2
        wso(i)  = -wso(i) *w2
c
      end do
c
c  Print potentials for each grid point
c
      if(idbg) then
        write(is,9110)
        do i=1,npt3
          rint = DBLE(i)*h
          write(is,9100) i,rint,-potr(i)/w2,-poti(i)/w2
     *                   ,-vso(i)/w2,-wso(i)/w2
        end do
        write(is,9080)
      endif
c
      ipl = IDNINT(two*si + one)
c
c=======================================================================
c  Solve equation for each l-value
c=======================================================================
c
      do l=1,lmaxc
        lmax = l
        fl    = DBLE(l-1)
        sl    = fl*(fl + one)
        fj    = fl - si - one
        fjmin = DABS(fl - si)
        do j=1,ipl
          fj  = fj + one
          if(fj .ge. fjmin) then
            spo = fj*(fj + one) - fl*(fl + one) - si*(si + one)
            if(si .eq. one) spo = spo/two
c
            call INTEG(npt,h,sl,spo)
c
c  PSIR  =                     real    part of the internal function
c  PSIRP = derivative of the   real    part of the internal function
c  PSII  =                   imaginary part of the internal function
c  PSIIP = derivative of the imaginary part of the internal function
c
            t1    = fc (l)
            t2    = fcp(l)*ak
            t3    = gc (l)
            t4    = gcp(l)*ak
            psim  = DMAX1(DABS(psir),DABS(psirp),DABS(psii),DABS(psiip))
            psir  = psir  / psim
            psirp = psirp / psim
            psii  = psii  / psim
            psiip = psiip / psim
            a1    = -(t1*psirp - t2*psir + t3*psiip - t4*psii)
            a2    = -(t1*psiip - t2*psii - t3*psirp + t4*psir)
            a3    = +(t1*psirp - t2*psir - t3*psiip + t4*psii)
            a4    = +(t1*psiip - t2*psii + t3*psirp - t4*psir)
            a5    =   a1*a1 + a2*a2
c
c  A1  =    real    part of the denominator
c  A2  =  imaginary part of the denominator
c  A3  =    real    part of the numerator
c  A4  =  imaginary part of the numerator
c  A5  =   modulus       of the denominator
c  ETR = 1.0 - real part of the matrix element
c  ETI =  imaginary part of the matrix element
c
            rap3    =   a3 / a1 - one
            rap4    =   a4 / a2 - one
            etr     = -(rap3*a1*a1 + rap4*a2*a2) / a5
            eti     =   a1*a2*(rap4 - rap3)      / a5
            br(j,l) = etr
            bi(j,l) = eti
            t (j,l) = one - (a3*a3 + a4*a4)/a5
          endif
        end do
c
        if(br(ipl,l) .le. zero) go to 205
        if(DABS(t(ipl,l)/t(ipl,1)) .lt. epstl) go to 205
      end do
      go to 210
c
  205 do j=1,ipl
        br(j,lmax) = zero
        bi(j,lmax) = zero
        t (j,lmax) = zero
      end do
      lmax = lmax - 1
c
  210 lm1 = lmax - 1
c
c  Volume integrals
c
      if(idr .le. 1) then
        call VOLIN(1,pote(1),rv(1),av(1),volint(1))
      else
        call VOLIN(1,pote(1),rv(1),av(1),volhf )
        call VOLIN(1,deltvv ,rv(1),av(1),volvol)
        call VOLIN(2,deltvs ,rv(4),av(4),volsur)
        volint(1) = volhf + volvol + volsur
      endif
      volint(1) = volint(1) / mt
c
      if(r(4) .gt. zero) then
        call VOLIN(2,pote(4),rv(4),av(4),volint(4))
        volint(4) = volint(4) / mt
        call VOLIN(2,pote(2),rv(4),av(4),volint(2))
        volint(2) = volint(2) / mt
      endif
c
      call VOLIN(1,pote(3),rv(3),av(3),volint(3))
      volint(3) = volint(3) / mt
      volint(5) = zero
      volint(6) = zero
c
c-----------------------------------------------------------------------
c  Print potential parameters
c-----------------------------------------------------------------------
c
      write(is,9010) el,amu,ak,eta,rm,h,npt,lm1
      write(is,9020)
      do j=1,nptypx
        redrad = rv(j) / mt13
        write(is,9030) type(j),pote(j),rv(j),redrad,av(j),volint(j)
      end do
c
      if(idr .ne. 0) then
        write(is,9030) typd(1),deltvv
        write(is,9030) typd(2),deltvs
      endif
c
      write(is,9070) rcoulb
      write(is,9080)
c
c  Formats
c
 9010 format(' ','E-l = ',1p,e12.5,3x
     *          ,'AMU =',e12.5,3x
     *          ,'K =',1p,e12.5,3x
     *          ,'ETA =',e12.5,6x
     *          ,'RM =',e12.5,2x
     *          ,'DR =',e12.5,6x
     *          ,i4,' POINTS',8x,'LMAX =',i3,//)
 9020 format(' ','POTENTIAL ',2x,'STRENGTH',2x,'RADIUS',2x
     *          ,'RED. RAD.',2x,'DIFFUSENESS',2x,'VOLUME INTEG.')
 9030 format(' ',a10,1x,f9.4,1x,f7.4,3x,f6.3,5x,f7.4,6x,f8.2)
 9070 format(' ','COULOMB   ',11x,f7.4)
 9080 format(' ',56x)
 9100 format(' ',i3,2x,f8.5,2x,1p,4e14.5)
 9110 format(' ',' I ',3x,' radius ',2x
     *          ,' real central ',' imag central '
     *          ,' real spin-or ',' imag spin-or ')
c
      return
      end
      subroutine SHAPEC(lmax,ida,na,ipl,ipu)
c***********************************************************************
c  Shape elastic differential cross-section (charged-particles)        *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'angles1.i'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'inout.i'
      INCLUDE 'poleg.i'
      INCLUDE 'tlj.i'
c
      integer          lmax,ida,na,ipl,ipu
c
      complex*16       a(nangmx),b(nangmx),c(nangmx),d(nangmx),e(nangmx)
      complex*16       CGAMMA
c     complex*16       chk
      complex*16       csigl(ltcmax),etac(3,ltcmax),fc(nangmx)
      complex*16       c1(ltcmax),c2(ltcmax),c3(ltcmax),c4(ltcmax)
     *                ,c5(ltcmax)
      complex*16       z,z1,z2,z3,z4,z5
      double precision ak,al,arg
      double precision cost,sint
      double precision sigl,dsig,dsigr
      double precision rap,tmpang
      double precision polar,pol20,pol21,pol22
      double precision sig0,som
      double precision twosr
      integer          i,j,jmax,jmin,l,naa
c
      dimension        sigl(ltcmax),dsig(nangmx),dsigr(nangmx)
      dimension        rap(nangmx),tmpang(nangmx)
      dimension        polar(nangmx),pol20(nangmx),pol21(nangmx)
     *                ,pol22(nangmx)
c
      double precision zero,one,two,three,ten
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
      data             three /3.0d+00/
      data             ten   /1.0d+01/
c
      complex*16       ai,onec,zeroc
      data             ai     /(0.0,1.0)/
      data             onec   /(1.0,0.0)/
      data             zeroc  /(0.0,0.0)/
c
      character        angtyp(2)*9
      data             angtyp /'     TETA','COS(TETA)'/
c=======================================================================
c
      twosr = DSQRT(two)
c
      ak    = DSQRT(ak2)
c
c  Initialisations
c
      do i=1,na
        rap(i)   = one
        dsig(i)  = zero
        dsigr(i) = one
        polar(i) = zero
        fc(i)    = zeroc
      end do
c
c  Store angles or cosines into a temporary array
c
      if(ida .eq. 1) then
        do i=1,na
          tmpang(i) =  ang(i)
        end do
      else if(ida .eq. 2) then
        do i=1,na
          tmpang(i) = cang(i)
        end do
      endif
c
      do l=1,lmax
        csigl(l) = onec
        do i=1,ipl
c90       etac(i,l) = DCMPLX(one - br(i,l) , bi(i,l))
          etac(i,l) =  CMPLX(one - br(i,l) , bi(i,l),8)
        end do
      end do
c
      if(eta .gt. zero) then
c
c  Coulomb phase shifts
c
c90     z = DCMPLX(one,eta)
        z =  CMPLX(one,eta,8)
        z = CGAMMA(z)
c90     sig0 = DIMAG(CDLOG(z))
        sig0 = AIMAG(  LOG(z))
c
        som = zero
        do l=1,lmax
          al       = DBLE  (l)
          sigl(l)  = sig0 + som
c90       csigl(l) = CDEXP(two*ai*sigl(l))
          csigl(l) =   EXP(two*ai*sigl(l))
          som      = som + DATAN2(eta,al)
        end do
c
c  Coulomb diffusion amplitude
c
        do i=1,na
          arg = (one - cang(i))/two
          if(arg .ne. zero) then
            z        = -ai*eta*DLOG(arg) + two*ai*sig0
c90         fc(i)    = -eta*CDEXP(z)/(two*ak*arg)
            fc(i)    = -eta*  EXP(z)/(two*ak*arg)
c90         dsigr(i) = CDABS(fc(i))**2
            dsigr(i) =   ABS(fc(i))**2
            dsigr(i) = ten*dsigr(i)
          endif
        end do
      endif
c
      if(ipl .eq. 1) then
c
c  spin 0 particles
c  ****************
        do l=1,lmax
          al    = DBLE  (l-1)
          c1(l) = csigl(l)*(two*al + one)*(onec - etac(1,l))
        end do
c
        do i=1,na
          z1 = zeroc
          do l=1,lmax
            z1 = z1 + c1(l)*pl(l,i)
          end do
          a(i)    = fc(i) + ai*z1/(two*ak)
c90       dsig(i) = CDABS(a(i))**2
          dsig(i) =   ABS(a(i))**2
          dsig(i) = ten*dsig(i)
        end do
c
        do i=2,na
          rap(i) = dsig(i)/dsigr(i)
        end do
c
      else if(ipl .eq. 2) then
c
c  spin 1/2 particles
c  ******************
        do l=1,lmax
          al    = DBLE  (l-1)
          c1(l) = csigl(l)*(
     *            (al+one)*(onec - etac(2,l)) + al*(onec - etac(1,l)) )
          c2(l) = csigl(l)*(etac(2,l) - etac(1,l))
        end do
c
        do i=1,na
          z1 = zeroc
          z2 = zeroc
          do l=1,lmax
            z1 = z1 + c1(l)*pl (l,i)
            z2 = z2 + c2(l)*pl1(l,i)
          end do
          a(i)     = fc(i) + ai*z1/(two*ak)
          b(i)     =       - ai*z2/(two*ak)
c90       dsig(i)  = CDABS(a(i))**2 + CDABS(b(i))**2
          dsig(i)  =   ABS(a(i))**2 +   ABS(b(i))**2
          dsig(i)  = ten*dsig(i)
c90       polar(i) = 2000.0*dimag( a(i)*dconjg(b(i)) )/dsig(i)
          polar(i) = 2000.0*AIMAG( a(i)* CONJG(b(i)) )/dsig(i)
        end do
c
        if(eta .gt. zero) then
          do i=2,na
            rap(i) = dsig(i)/dsigr(i)
          end do
        endif
c
      else if(ipl .eq. 3) then
c
c  spin 1 particles
c  ****************
        do l=1,lmax
          al    = DBLE  (l-1)
          c1(l) = -csigl(l)*(
     *            (    al+one)     *(onec - etac(3,l)) +
     *            (    al    )     *(onec - etac(1,l)) )
          c2(l) = -csigl(l)*(
     *            (    al+two)     *(onec - etac(3,l)) +
     *            (two*al+one)     *(onec - etac(2,l)) +
     *            (    al-one)     *(onec - etac(1,l)) )
          c3(l) = -csigl(l)*(
     *                              (onec - etac(3,l)) -
     *                              (onec - etac(1,l)) )
          c4(l) = -csigl(l)*(
     *            (al*(al+two))    *(onec - etac(3,l)) -
     *            (two*al+one)     *(onec - etac(2,l)) -
     *            (al-one)*(al+one)*(onec - etac(1,l)) )
          c5(l) = -csigl(l)*(
     *            (    al    )     *(onec - etac(3,l)) -
     *            (two*al+one)     *(onec - etac(2,l)) +
     *            (    al+one)     *(onec - etac(1,l)) )
        end do
c
        do i=2,na
          cost = pl(2,i)
          sint = DSQRT(one - cost*cost)
          z1 = zeroc
          z2 = zeroc
          z3 = zeroc
          z4 = zeroc
          z5 = zeroc
          do l=1,lmax
            al = DBLE  (l - 1)
            z1 = z1 + c1(l)*pl (l,i)
            z2 = z2 + c2(l)*pl (l,i)/two
            z3 = z3 + c3(l)*pl1(l,i)/twosr
          end do
c
          do l=2,lmax
            al = DBLE  (l - 1)
            z4 = z4 + c4(l)*pl1(l,i)/(al*(al+one)*twosr)
            z5 = z5 + c5(l)*pl2(l,i)/(al*(al+one)*two)
          end do
c
          a(i) = fc(i) + z1/(two*ai*ak)
          b(i) = fc(i) + z2/(two*ai*ak)
          c(i) =         z3/(two*ai*ak)
          d(i) =         z4/(two*ai*ak)
          e(i) =         z5/(two*ai*ak)
c         chk  = b(i) - a(i) - e(i) + twosr*(d(i) - c(i))*cost/sint
c90       dsig(i)  = CDABS(a(i))**2 + two*( CDABS(b(i))**2 +
c90  *               CDABS(c(i))**2 + CDABS(d(i))**2 + CDABS(e(i))**2 )
          dsig(i)  =   ABS(a(i))**2 + two*(   ABS(b(i))**2 +
     *                 ABS(c(i))**2 +   ABS(d(i))**2 +   ABS(e(i))**2 )
          dsig(i)  = ten*dsig(i)/three
c90       polar(i) = 1.0d+03*two*twosr*DIMAG( a(i)*dconjg(c(i)) +
c90  *               b(i)*DCONJG(d(i)) + d(i)*DCONJG(e(i)) ) /
c90  *               (three*dsig(i))
          polar(i) = 1.0d+03*two*twosr*AIMAG( a(i)* CONJG(c(i)) +
     *               b(i)* CONJG(d(i)) + d(i)* CONJG(e(i)) ) /
     *               (three*dsig(i))
c90       pol20(i) = 1.0d+02*(one - ten*( CDABS(a(i))**2 +
c90  *               two*CDABS(d(i))**2 ) / dsig(i) ) / twosr
          pol20(i) = 1.0d+02*(one - ten*(   ABS(a(i))**2 +
     *               two*  ABS(d(i))**2 ) / dsig(i) ) / twosr
c90       pol21(i) =-1.0d+03*twosr*DREAL( a(i)*DCONJG(c(i)) -
c90  *                b(i)*DCONJG(d(i)) + d(i)*DCONJG(e(i)) ) /
c90  *                                 (dsig(i)*DSQRT(three))
          pol21(i) =-1.0d+03*twosr*DBLE ( a(i)* CONJG(c(i)) -
     *                b(i)* CONJG(d(i)) + d(i)* CONJG(e(i)) ) /
     *                                 (dsig(i)*DSQRT(three))
c90       pol22(i) = 1.0d+03*( two*dreal(b(i)*DCONJG(e(i))) -
c90  *               CDABS(c(i))**2) / (dsig(i)*DSQRT(three))
          pol22(i) = 1.0d+03*( two*DBLE (b(i)* CONJG(e(i))) -
     *                 ABS(c(i))**2) / (dsig(i)*DSQRT(three))
        end do
c
        do i=2,na
          rap(i) = dsig(i)/dsigr(i)
        end do
c
      endif
c
c  Print angular distributions
c
      write(is,9000)
c
      if(ipl .ne. 3) then
c
c  Spin 0 and 1/2 particles
        naa = na/2 + 1
        write(is,9010) angtyp(ida),angtyp(ida)
        do i=1,naa
          jmin = 2*(i-1) + 1
          jmax = MIN0(2*i,na)
          write(is,9030)
     *    (tmpang(j),dsig(j),dsigr(j),rap(j),polar(j),j=jmin,jmax)
        end do
      else
c
c  Spin 1 particles
        naa = na
        write(is,9015) angtyp(ida),angtyp(ida)
        do i=1,naa
          write(is,9035)
     *    tmpang(i),dsig(i),dsigr(i),rap(i),polar(i),pol20(i),
     *    pol21(i),pol22(i),ang(i)
        end do
      endif
c
      if(ipu .ne. 0) then
        if(ipl .ne. 3) then
          do i=1,na
            write(is3,9040) tmpang(i),dsig(i),dsigr(i),rap(i),polar(i)
          end do
        else
          do i=1,na
            write(is3,9040) tmpang(i),dsig(i),dsigr(i),rap(i),polar(i),
     *      pol20(i),pol21(i),pol22(i)
          end do
        endif
      endif
c
c  Formats
c
 9000 format('1',46x,'SHAPE ELASTIC DIFFERENTIAL CROSS-SECTION',/
     *      ,' ',46x,40('*'),/
     *      ,' ',56x,'CENTER-OF-MASS SYSTEM',///)
 9010 format(' ',4x,a9         ,3x,'DSIG(TETA)',2x,'DSIGR(TETA)'
     *          ,3x,'DSIG/DSIGR',1x,'POLARIZATION',2x,
     *           4x,a9         ,3x,'DSIG(TETA)',2x,'DSIGR(TETA)'
     *          ,3x,'DSIG/DSIGR',1x,'POLARIZATION',/)
 9015 format(' ',4x,a9         ,3x,'DSIG(TETA)',2x,'DSIGR(TETA)'
     *          ,3x,'DSIG/DSIGR',2x,'POLAR. VECT',4x,'POLAR. 20'
     *          ,4x,'POLAR. 21' ,4x,'POLAR. 22'  ,4x,a9         ,/)
 9030 format(' ',1p,5e13.5,2x,5e13.5)
 9035 format(' ',1p,9e13.5)
 9040 format(1p,8e10.3)
c
      return
      end
      subroutine SHAPEL2(lmax,ida,na,ipl)
c***********************************************************************
c  Shape elastic differential cross-section for neutrons               *
c----------------------------------------------------------------------*
c  BB are Blatt and Biedenharn coefficients, see tabulated values      *
c       in L.C.Biedenharn  ORNL-1501 (1953)                            *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'angles1.i'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'inout.i'
      INCLUDE 'poleg.i'
      INCLUDE 'tlj.i'
c
      integer          lmax,ida,na,ipl
c
      double precision aj1,aj2,aj2max,aj2min,al,al1,al2
      double precision b1,bb,bl,br1,bi1,da,z1,z2
      double precision cg,rac
      double precision rad,s,som,sym
      integer          idaa,j,kpar,l,lmax2,naa,ncoup,nl
      integer          i,imax,imin
      integer          ij1,ij2,il,il1,il2,il2max,il2min,js1,js2
      integer          lm,lma,lmi
c
      dimension        bl(ldamax),da(nangmx)
c
      double precision zero,eps,one,two,four,eight,ten
      data             zero  /0.0d+00/
      data             eps   /1.0d-03/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
      data             four  /4.0d+00/
      data             eight /8.0d+00/
      data             ten   /1.0d+01/
c
      double precision CLEBG,RACAH
      external         CLEBG,RACAH
c=======================================================================
c
      rad   = 1.8d+02/pi
c
      idaa  = IABS(ida)
      lmax2 = 2*lmax - 1
      ncoup = 0
c
      do il=1,lmax2
        al  = DBLE(il-1)
        som = zero
c
        do il1=1,lmax
          al1 = DBLE(il1-1)
          aj1 = al1 - si - one
          il2min = MAX0(il1,IABS(il1-il)+1)
          kpar   = il + il1 + il2min
          if(2*(kpar/2) .eq. kpar) then
            il2min = il2min + 1
          endif
          il2max = MIN0(il1 + il - 1,lmax)
          if(il2min .le. il2max) then
c
            do ij1=1,ipl
              aj1 = aj1 + one
              js1 = 100*il1 + ij1
              if(aj1 .ge. DABS(al1-si) ) then
                br1 = br(ij1,il1)
                bi1 = bi(ij1,il1)
                z1  = (two*al1 + one)*(two*aj1 + one)
                aj2min = DMAX1(aj1,DABS(aj1 - al)) - eps
                aj2max = aj1 + al + eps
c
                do il2=il2min,il2max,2
                  al2 = DBLE(il2-1)
                  aj2 = al2 - si - one
c
                  do ij2=1,ipl
                    aj2 = aj2 + one
                    js2 = 100*il2 + ij2
                    if((aj2-aj2min)*(aj2-aj2max) .le. zero) then
                      if(js1 - js2 .ne. 0) then
                        sym = 2.0
                      else
                        sym = 1.0
                      endif
                      cg  = CLEBG(al1,al2,al,zero,zero,zero)
                      rac = RACAH(al1,aj1,al2,aj2,si,al)
                      z2  = (two*al2 + one)*(two*aj2 + one)
                      bb  = z1*z2*cg*cg*rac*rac
                      s   = br1*br(ij2,il2) + bi1*bi(ij2,il2)
                      som = som + bb*sym*s
                      ncoup = ncoup + 1
                    endif
                  end do
                end do
              endif
            end do
          endif
        end do
        bl(il) = ten*som/(eight*ak2)
      end do
c
c  Angular distribution
c
      do i=1,na
        som = zero
        do l=1,lmax2
          som = som + bl(l)*pl(l,i)
        end do
        da(i) = som
      end do
c
c  Integral of the angular distribution
c
      som = zero
      if(idaa .eq. 1) then
        do i=2,na
          som = som + (da(i)  *DSIN(ang(i)/rad) +
     *                 da(i-1)*DSIN(ang(i-1)/rad)) / two
        end do
        som = som*2.5/rad
      else if(idaa .eq. 2) then
        do i=1,na
          som = som + da(i)
        end do
        som = som - (da(1) + da(na)) / two
        som = som*0.02
      endif
c
      som = two*pi*som
c
c  Print angular distribution
c
      write(is,9000)
      naa = na/4 + 1
      if(idaa .eq. 1) then
        write(is,9010)
        do i=1,naa
          imin = 4*(i-1) + 1
          imax = MIN0(4*i,na)
          write(is,9030) (ang(j),da(j),j=imin,imax)
        end do
      else if(idaa .eq. 2) then
        write(is,9020)
        do i=1,naa
          imin = 4*(i-1) + 1
          imax = MIN0(4*i,na)
          write(is,9030) (cang(j),da(j),j=imin,imax)
        end do
      endif
c
      b1 = four*pi*bl(1)
      write(is,9040) som,b1
      write(is,9050) ncoup
c
c  Print Legendre coefficients (ENDF prescription)
c
      do l=2,lmax2
        al = DBLE  (l-1)
        bl(l) = bl(l)/((two*al + one)*bl(1))
      end do
      bl(1) = one
      write(is,9060)
      nl = (lmax2 + 4)/5
      do l=1,nl
        lmi = 5*(l-1)
        lma = MIN0(5*l,lmax2) - 1
        write(is,9070) (lm,bl(lm+1),lm=lmi,lma)
      end do
c
c  Formats
c
 9000 format('1',46x,'SHAPE ELASTIC DIFFERENTIAL CROSS-SECTION',/
     *      ,' ',46x,40('*'),/
     *      ,' ',56x,'CENTER-OF-MASS SYSTEM',///)
 9010 format(' ',5x,4('    TETA ',2x,'D.SIGMA/D.OMEGA',6x),/)
 9020 format(' ',5x,4('COS(TETA)',2x,'D.SIGMA/D.OMEGA',6x),/)
 9030 format(' ',5x,4(1p,e12.5,2x,e12.5,6x))
 9040 format(//,' ','INTEGRAL =',1p,e12.5,' MB',10x
     *             ,'BL(0) ='      ,e12.5,' MB')
 9050 format(' ',90x,'NCOUP =',i8)
 9060 format(///,' ',7x,5(2x,'L',7x,'BL(L)',10x),/)
 9070 format(' ',7x,5(i3,3x,1p,e14.7,5x))
c
      return
      end
      subroutine SORTR(a,n,iperm)
c***********************************************************************
c  Sort by increasing order the array A (dimension N)                  *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      double precision a
      integer          iperm,n
      dimension        a(n),iperm(n)
c
      double precision atemp
      integer          iptemp
      integer          i,i1,j,n1
c=======================================================================
c
      if(n .le. 1) return
c
c  Initialize the permutation array
c
      do i=1,n
        iperm(i) = i
      end do
c
c  Sort the array
c
      n1 = n - 1
      do i=1,n1
        i1 = i + 1
        do j=i1,n
          if(a(j) .lt. a(i)) then
            atemp    = a(j)
            a(j)     = a(i)
            a(i)     = atemp
            iptemp   = iperm(j)
            iperm(j) = iperm(i)
            iperm(i) = iptemp
          endif
        end do
      end do
c
      return
      end
      subroutine SPIN0(n,lmax,ipr,ipu)
c***********************************************************************
c  Print transmission coefficients T(l,j)                              *
c  Calculate cross sections: - compound nucleus                        *
c                              shape elastic                           *
c                              total                                   *
c  for spin 0 incident particles                                       *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'ener.i'
      INCLUDE 'inout.i'
      INCLUDE 'tce.i'
      INCLUDE 'tlj.i'
      INCLUDE 'xs.i'
c
      integer          n,lmax,ipr,ipu
c
      double precision al2,tm
      integer          k,l
c
      double precision zero,one,two,ten
      data             zero /0.0d+00/
      data             one  /1.0d+00/
      data             two  /2.0d+00/
      data             ten  /1.0d+01/
c=======================================================================
c
      se(n) = zero
      sr(n) = zero
      st(n) = zero
c
      do l=1,lmax
        k   = l - 1
        al2 = two*DBLE(l-1) + one
        tm  = t(1,l)
        tc(n,l) = tm
        sr(n) = sr(n) +
     *    al2*(two*br(1,l) - br(1,l)*br(1,l) - bi(1,l)*bi(1,l))
        se(n) = se(n) + al2*(br(1,l)*br(1,l) + bi(1,l)*bi(1,l))
        st(n) = st(n) + al2*two*br(1,l)
      end do
c
      if(ipr .ne. 0) then
        write(is,9010)
        do l=1,lmax
          k = l - 1
          write(is,9020) k,tc(n,l),br(1,l),bi(1,l)
        end do
      endif
c
c-----------------------------------------------------------------------
c  Cross sections (in mb)
c-----------------------------------------------------------------------
c
      se(n) = ten*se(n)*pi/ak2
      sr(n) = ten*sr(n)*pi/ak2
      st(n) = ten*st(n)*pi/ak2
c
c     el    = e1*(mi + mt)/mt
      write(is,9030) e1,sr(n),el,se(n),st(n)
c
      if(ipu .ne. 0) then
        write(is3,'(1x,f7.3,1x,3(f10.2,1x))') el,st(n),sr(n),se(n)
      endif
c
c  Formats
c
 9010 format(' ',3x,'L',7x,'TC(L)',3x,'1 - ETA R',7x,'ETA I',/)
 9020 format(' ',i4,1p,3e12.4)
 9030 format(//,' ',   f8.4,' MEV CM' ,5x
     *                 ,'COMPOUND NUCLEUS CROSS SECTION',1p,e14.7
     *       ,/,' ',0p,f8.4,' MEV LAB',7x
     *                 ,'SHAPE ELASTIC CROSS SECTION',   1p,e14.7
     *       ,/,' ',31x,'TOTAL CROSS SECTION',              e14.7
     *                 ,' MB',//)
 9050 format(1p,4e10.3)
c
      return
      end
      subroutine SPIN05(n,lmax,ipr,ipu)
c***********************************************************************
c  Print transmission coefficients T(l,j)                              *
c  Calculate cross sections: - compound nucleus                        *
c                              shape elastic                           *
c                              total                                   *
c  for spin 1/2 incident particles                                     *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'ener.i'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'tce.i'
      INCLUDE 'tlj.i'
      INCLUDE 'xs.i'
c
      integer          n,lmax,ipr,ipu
c
      double precision al,al1,al2,tm
      double precision b1,b2
      double precision p1,r2,rp,s0,s1
      integer          i,k,l
c
      double precision zero,one,two,four,ten
      data             zero /0.0d+00/
      data             one  /1.0d+00/
      data             two  /2.0d+00/
      data             four /4.0d+00/
      data             ten  /1.0d+01/
c=======================================================================
c
      se(n) = zero
      sr(n) = zero
      st(n) = zero
c
      do l=1,lmax
        k   = l - 1
        al  = DBLE(k)
        al1 = al + one
        al2 = two*al + one
        tm  = (al1*t(2,l) + al*t(1,l)) / al2
        tc(n,l) = tm
        sr(n)   = sr(n) +
     *            al1*(two*br(2,l) - br(2,l)*br(2,l) - bi(2,l)*bi(2,l))
     *          + al *(two*br(1,l) - br(1,l)*br(1,l) - bi(1,l)*bi(1,l))
        se(n)   = se(n) +
     *            al1*(br(2,l)*br(2,l) + bi(2,l)*bi(2,l))
     *          + al *(br(1,l)*br(1,l) + bi(1,l)*bi(1,l))
        st(n)   = st(n) +
     *            al1*(two*br(2,l))
     *          + al *(two*br(1,l))
      end do
c
      if(ipr .ne. 0) then
        write(is,9010)
        do l=1,lmax
          k   = l - 1
          al  = DBLE  (k)
          al1 = al + one
          al2 = two*al + one
          b1  = (al1*br(2,l) + al*br(1,l)) / al2
          b2  = (al1*bi(2,l) + al*bi(1,l)) / al2
          write(is,9020) k,tc(n,l),b1,b2,(t(i,l),br(i,l),bi(i,l),i=1,2)
        end do
      endif
c
c-----------------------------------------------------------------------
c  Cross sections (in mb)
c-----------------------------------------------------------------------
c
      se(n) = ten*se(n)*pi/ak2
      sr(n) = ten*sr(n)*pi/ak2
      st(n) = ten*st(n)*pi/ak2
c
c     el    = e1*(mi + mt)/mt
      write(is,9030) e1,sr(n),el,se(n),st(n)
c
      if(ipu .ne. 0) then
        write(is3,'(1x,f7.3,1x,3(f10.2,1x))') el,st(n),sr(n),se(n)
      endif
c
c-----------------------------------------------------------------------
c  Strength functions, scattering radius
c-----------------------------------------------------------------------
c
      if(e1 .lt. 1.0d-01) then
        s0 = tc(n,1)/(two*pi*DSQRT(1.0d+06*e1))
c
        rp = (r(1) + re(1)*el)*(mt**0.333333333)
        r2 = rp*rp
        p1 = (ak2*r2)/(one + ak2*r2)
        s1 = tc(n,2)/(two*pi*p1*DSQRT(1.0d+06*e1))
c
        rp = DSQRT(se(n)/(four*ten*pi))
        write(is,9040) s0,s1,rp
      endif
c
c  Formats
c
 9010 format(' ',3x,'L',7x,'TC(L)',3x
     *             ,'1 - ETA R',7x,'ETA I',9x,'T(L,L-1/2)',3x
     *             ,'1 - ETA R',7x,'ETA I',9x,'T(L,L+1/2)',3x
     *             ,'1 - ETA R',7x,'ETA I',/)
 9020 format(' ',i4,1p,3e12.4,2(7x,3e12.4))
 9030 format(//,' ',   f8.4,' MEV CM' ,5x
     *                 ,'COMPOUND NUCLEUS CROSS SECTION',1p,e14.7
     *       ,/,' ',0p,f8.4,' MEV LAB',7x
     *                 ,'SHAPE ELASTIC CROSS SECTION',   1p,e14.7
     *       ,/,' ',31x,'TOTAL CROSS SECTION',              e14.7
     *                 ,' MB',//)
 9040 format(' ',27x,'STRENGTH FUNCTIONS   S0',1p,e14.7,/
     *      ,' ',48x,'S1',                        e14.7,/
     *      ,' ',33x,'SCATTERING RADIUS',         e14.7,//)
 9050 format(1p,4e10.3)
c
      return
      end
      subroutine SPIN1(n,lmax,ipr,ipu)
c***********************************************************************
c  Print transmission coefficients T(l,j)                              *
c  Calculate cross sections: - compound nucleus                        *
c                              shape elastic                           *
c                              total                                   *
c  for spin 1 incident particles                                       *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'ener.i'
      INCLUDE 'inout.i'
      INCLUDE 'tce.i'
      INCLUDE 'tlj.i'
      INCLUDE 'xs.i'
c
      integer          n,lmax,ipr,ipu
c
      double precision al,al1,al2,al3,tm
      double precision b1,b2
      integer          i,k,l
c
      double precision zero,one,two,three,ten
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
      data             three /3.0d+00/
      data             ten   /1.0d+01/
c=======================================================================
c
      se(n) = zero
      sr(n) = zero
      st(n) = zero
c
      do l=1,lmax
        k   = l - 1
        al  = DBLE(l-1)
        al1 = two*al - one
        al2 = two*al + one
        al3 = two*al + three
        tm  = (al3*t(3,l) + al2*t(2,l) + al1*t(1,l)) / (three*al2)
        tc(n,l) = tm
        sr(n)   = sr(n) +
     *            al3*(two*br(3,l) - br(3,l)*br(3,l) - bi(3,l)*bi(3,l))
     *          + al2*(two*br(2,l) - br(2,l)*br(2,l) - bi(2,l)*bi(2,l))
     *          + al1*(two*br(1,l) - br(1,l)*br(1,l) - bi(1,l)*bi(1,l))
        se(n)   = se(n) +
     *            al3*(br(3,l)*br(3,l) + bi(3,l)*bi(3,l))
     *          + al2*(br(2,l)*br(2,l) + bi(2,l)*bi(2,l))
     *          + al1*(br(1,l)*br(1,l) + bi(1,l)*bi(1,l))
        st(n)   = st(n) +
     *            al3*two*br(3,l)
     *          + al2*two*br(2,l)
     *          + al1*two*br(1,l)
      end do
c
      if(ipr .ne. 0) then
        write(is,9010)
        do l=1,lmax
          k   = l - 1
          al  = DBLE(k)
          al1 = two*al - one
          al2 = two*al + one
          al3 = two*al + three
          b1 = (al3*br(3,l) + al2*br(2,l) + al1*br(1,l))/(three*al2)
          b2 = (al3*bi(3,l) + al2*bi(2,l) + al1*bi(1,l))/(three*al2)
          write(is,9020) k,tc(n,l),b1,b2,(t(i,l),br(i,l),bi(i,l),i=1,3)
        end do
      endif
c
c-----------------------------------------------------------------------
c  Cross sections (in mb)
c-----------------------------------------------------------------------
c
      se(n) = se(n)/three
      sr(n) = sr(n)/three
      st(n) = st(n)/three
c
      se(n) = ten*se(n)*pi/ak2
      sr(n) = ten*sr(n)*pi/ak2
      st(n) = ten*st(n)*pi/ak2
c
c     el    = e1*(mi + mt)/mt
      write(is,9030) e1,sr(n),el,se(n),st(n)
c
      if(ipu .ne. 0) then
        write(is3,'(1x,f7.3,1x,3(f10.2,1x))') el,st(n),sr(n),se(n)
      endif
c
c  Formats
c
 9010 format(' ',1x,'L',7x,'TC(L)',3x
     *             ,'1 - ETA R',7x,'ETA I',3x,'T(L,L-1)',1x
     *             ,'1 - ETA R',5x,'ETA I',5x,'T(L,L)'  ,1x
     *             ,'1 - ETA R',5x,'ETA I',3x,'T(L,L+1)',1x
     *             ,'1 - ETA R',5x,'ETA I',/)
 9020 format(' ',i2,1p,3e12.4,3(1x,3e10.3))
 9030 format(//,' ',   f8.4,' MEV CM' ,5x
     *                 ,'COMPOUND NUCLEUS CROSS SECTION',1p,e14.7
     *       ,/,' ',0p,f8.4,' MEV LAB',7x
     *                 ,'SHAPE ELASTIC CROSS SECTION',   1p,e14.7
     *       ,/,' ',31x,'TOTAL CROSS SECTION',              e14.7
     *                 ,' MB',//)
 9050 format(1p,4e10.3)
c
      return
      end
      subroutine SYSPOT(ip,ipot,ref)
c***********************************************************************
c  Potentials from the compilation by C.M.Perey et F.G.Perey           *
c  Atomic Data and Nuclear Data Tables   17 (1976) 1-101               *
c----------------------------------------------------------------------*
c                      new # (2000 version)         old #              *
c  Type of potential:  1  real      volume          1                  *
c                      2  real      surface                            *
c                      3  imaginary volume          3                  *
c                      4  imaginary surface         2                  *
c                      5  real      spin-orbit      4                  *
c                      6  imaginary spin-orbit      5                  *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
c
      character        ref*(*)
      integer          ip,ipot
c
      double precision mt13,mt23,mt43,mt2,mt3,nmzsa
      double precision one3,two3
      integer          i,j,k
c
      double precision zero,one,two,three
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
      data             three /3.0d+00/
c=======================================================================
c
      one3  = one / three
      two3  = two / three
c
c  Powers of A
c
      mt13  = mt**one3
      mt23  = mt**two3
      mt43  = mt**(one + one3)
      mt2   = mt*mt
      mt3   = mt*mt*mt
      nmzsa = (mt - two*zt)/mt
c
c  Potential parameters initialisations
c
      do i=1,nptypx
        a (i) = zero
        ae(i) = zero
        r (i) = zero
        re(i) = zero
        npzen(i)  = 1
        epot(i,0) = zero
        epot(i,1) = epotmx
        do j=1,npzenx
          do k=1,npcofx
            pot(i,j,k) = zero
          end do
        end do
      end do
      beta   = zero
      rcoul  = zero
c
      if(ip .eq. 1) then
c*************
c  neutrons  *
c*************
c
        if(ipot .eq. 1) then
          rcoul = zero
c
c     WILMORE - HODGSON parameters                    40 < A
c                                                          E < 10
          ref  = 'WILM.HODGS'
          r(1) = 1.322 - 7.6e-04*mt + 4.0e-06*mt2 - 8.0e-09*mt3
          r(4) = 1.266 - 3.7e-04*mt + 2.0e-06*mt2 - 4.0e-09*mt3
          r(5) = r(1)
          a(1) = 0.66
          a(4) = 0.48
          a(5) = a(1)
          pot(1,1,1) =  47.01
          pot(1,1,2) = -0.267
          pot(1,1,3) = -0.00118
          pot(4,1,1) =  9.52
          pot(4,1,2) = -0.053
          pot(5,1,1) =  7.0
        else if(ipot .eq. 2) then
c
c     BECCHETTI - GREENLEES parameters                40 < A
c                                                     10 < E < 50
          ref  = 'BECC.GREEN'
          r(1) = 1.17
          r(3) = 1.26
          r(4) = 1.26
          r(5) = 1.01
          a(1) = 0.75
          a(3) = 0.58
          a(4) = 0.58
          a(5) = 0.75
          pot(1,1,1) =  56.3 - 24.0*nmzsa
          pot(1,1,2) = -0.32
          pot(3,1,1) = -1.56
          pot(3,1,2) =  0.22
          pot(4,1,1) =  13.0 - 12.0*nmzsa
          pot(4,1,2) = -0.25
          pot(5,1,1) =  6.2
        else if(ipot .eq. 3) then
c
c     FERRER et al. parameters                        24 < A < 209
c     NUCL.PHYS. A275(1977)325-341                         E = 11
c
          ref  = 'FERRER'
          r(1) = 1.27
          r(4) = 1.27
          r(5) = 1.08
          a(1) = 0.71
          a(4) = 0.434
          a(5) = 0.71
          pot(1,1,1) = 47.14 - 22.50*nmzsa
          pot(4,1,1) = 12.16 -  2.03*nmzsa
          pot(5,1,1) = 4.55
        else if(ipot .eq. 4) then
c
c     BERSILLON - CINDRO parameters
c     CONTRIBUTION TO THE 5TH INTERNATIONAL SYMPOSIUM ON
c     INTERACTIONS OF FAST NEUTRONS WITH NUCLEI,
c     GAUSSIG, DDR, 17-21 NOV 1975
c
          ref  = 'CIND.BERS.'
          beta = 0.85
          r(1) = 1.182 + 1.93e-04*mt
          r(4) = 1.21
          r(5) = r(1)
          a(1) = 0.65
          a(4) = 0.47
          a(5) = a(1)
          pot(1,1,1) = 71.
          pot(4,1,1) = 7.
          pot(4,1,2) = 0.4
          pot(5,1,1) = 7.
        else if(ipot .eq. 5) then
c
c     MADLAND parameters                         ACTINIDES
c     HARWELL CONFERENCE                         E < 10 MEV
c     SEPTEMBER 25-29, 1978
c     -----  TEMPORARY VALUES  -----
c
          ref   = 'MADLAND'
          r(1)  = 1.264
          r(4)  = 1.256
          r(5)  = 1.01
          a(1)  = 0.612
          a(4)  = 0.553
          a(5)  = 0.75
          ae(4) = 0.0144
          pot(1,1,1) = 50.378 - 27.073*nmzsa
          pot(1,1,2) = -0.354
          pot(4,1,1) =  9.265 - 12.666*nmzsa
          pot(4,1,2) = -0.232
          pot(4,1,3) =  0.03318
          pot(5,1,1) =  6.2
        endif
c
      else if(ip .eq. 2) then
c************
c  protons  *
c************
        if(ipot .eq. 1) then
c
c     PEREY parameters                           30 < A < 100
c                                                     E < 20
          ref  = 'PEREY'
          rcoul= 1.25
          r(1) = 1.25
          r(4) = 1.25
          r(5) = 1.25
          a(1) = 0.65
          a(4) = 0.47
          a(5) = 0.47
          pot(1,1,1) =  53.3 + 27.0*nmzsa + 0.4*zt/mt13
          pot(1,1,2) = -0.55
          pot(4,1,1) =  13.5
          pot(5,1,1) =  7.5
        else if(ipot .eq. 2) then
c
c     BECCHETTI - GREENLEES parameters            A < 40
c                                            20 < E < 50
c     RCOUL from L.R.B.ELTON, "Nuclear Sizes", Oxford Univ.Press,
c                              London, 1961, p.36
          ref  = 'BECC.GREEN'
          rcoul= 1.123 + 2.35141/mt23 - 2.06789/mt43 + 13.78687/mt2
          r(1) = 1.17
          r(3) = 1.32
          r(4) = 1.32
          r(5) = 1.01
          a(1) = 0.75
          a(3) = 0.51 + 0.7*nmzsa
          a(4) = 0.51 + 0.7*nmzsa
          a(5) = 0.75
          pot(1,1,1) =  54.0 + 24.0*nmzsa + 0.4*zt/mt13
          pot(1,1,2) = -0.32
          pot(3,1,1) = -2.7
          pot(3,1,2) =  0.22
          pot(4,1,1) =  11.8 + 12.0*nmzsa
          pot(4,1,2) = -0.25
          pot(5,1,1) =  6.2
        endif
      else if(ip .eq. 3) then
c**************
c  deuterons  *
c**************
c
        if(ipot .eq. 1) then
c
c     LOHR - HAEBERLI parameters                  40 < A
c                                                  8 < E < 13
          ref  = 'LOHR-HAEB.'
          rcoul= 1.3
          r(1) = 1.05
          r(4) = 1.43
          r(5) = 0.75
          a(1) = 0.86
          a(4) = 0.50 + 0.013*mt23
          a(5) = 0.50
          pot(1,1,1) = 91.13 + 2.2*zt/mt13
          pot(4,1,1) = 218.0/mt23
          pot(5,1,1) = 7.0
        else if(ipot .eq. 2) then
c
c     PEREY parameters                                12 < E < 25
c
          ref  = 'PEREY'
          rcoul= 1.15
          r(1) = 1.15
          r(4) = 1.34
          a(1) = 0.81
          a(4) = 0.68
          pot(1,1,1) =  81.0 + 2.0*zt/mt13
          pot(1,1,2) = -0.22
          pot(4,1,1) =  14.4
          pot(4,1,2) =  0.24
        endif
      else if(ip .eq. 4) then
c************
c  tritons  *
c************
        if(ipot .eq. 1) then
c
c     BECCHETTI - GREENLEES parameters                40 < A
c                                                          E < 40
          ref  = 'BECC.GREEN'
          rcoul= 1.30
          r(1) = 1.20
          r(3) = 1.40
          r(5) = 1.20
          a(1) = 0.72
          a(3) = 0.84
          a(5) = 0.72
          pot(1,1,1) =  165.0 - 6.4*nmzsa
          pot(1,1,2) = -0.17
          pot(3,1,1) =  46.0  - 110.0*nmzsa
          pot(3,1,2) = -0.33
          pot(5,1,1) =  2.5
        endif
      else if(ip .eq. 5) then
c*************
c  helium-3  *
c*************
        if(ipot .eq. 1) then
c
c     BECCHETTI - GREENLEES parameters                40 < A
c                                                          E < 40
          ref  = 'BECC.GREEN'
          rcoul= 1.30
          r(1) = 1.20
          r(3) = 1.40
          r(5) = 1.20
          a(1) = 0.72
          a(3) = 0.88
          a(5) = 0.72
          pot(1,1,1) =  151.9 + 50.0*nmzsa
          pot(1,1,2) = -0.17
          pot(3,1,1) =  41.7 - 44.0*nmzsa
          pot(3,1,2) = -0.33
          pot(5,1,1) =  2.5
        endif
      else if(ip .eq. 6) then
c***********
c  alphas  *
c***********
c
        if(ipot .eq. 1) then
c
c     Averaged parameters
c     MAC FADDEN ET SATCHLER  NUCL.PHYS. 84(1966)177
c
          ref  = 'MAC FADDEN'
          rcoul= 1.40
          r(1) = 1.40
          r(3) = 1.40
          a(1) = 0.52
          a(3) = 0.52
          pot(1,1,1) = 185.
          pot(3,1,1) =  25.
        endif
      endif
c
      return
      end
      subroutine TLSTAP(ne,lmaxtl,lmaxim)
c***********************************************************************
c  Print the j-averaged transmission coefficients T(l)                 *
c  in the STAPRE format                                                *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'ener.i'
      INCLUDE 'inout.i'
      INCLUDE 'tce.i'
c
      integer          ne,lmaxtl,lmaxim
c
      character*80     card,titre
      integer          k,l,lmaxt,n,n1,n2,ne1,nen,nzero
c
      double precision zero
      data             zero /0.0d+00/
c=======================================================================
c
      write(is,fmt='(6hTLSTAP,3i3)') ne,lmaxtl,lmaxim
c
      if(lmaxim .eq. 0) then
c
c  First pass
c
        do l=1,lmaxtl
          lmaxt = l
          nzero = 0
          do n=1,ne
            if(tc(n,l) .le. zero) then
              nzero = nzero + 1
            endif
          end do
          nen = ne - nzero
          if(nen .lt. 4) go to 40
        end do
   40   lmaxtl = lmaxt - 1
        write(is2,fmt='(i3)') lmaxtl
c
        do l=1,lmaxtl
          nzero = 0
          do n=1,ne
            if(tc(n,l) .le. zero) then
              nzero = nzero + 1
            endif
          end do
          nen = ne - nzero
          write(is2,fmt='(i3)') nen
          ne1 = (nen + 3) / 4
          do n=1,ne1
            n1 = 4*(n-1)  + nzero + 1
            n2 = min0(4*n + nzero,ne)
            write(is2,9010) (ein(k),tc(k,l),k=n1,n2)
          end do
        end do
c
      else
c
c  Second pass
c
        ENDFILE is2
        REWIND  is2
        read (ie ,fmt='(a)')  titre
        write(is1,fmt='(a)')  titre
        write(is1,fmt='(i3)') lmaxim
  110   read(is2,fmt='(a)',end=130) card
        write(is1,fmt='(a)') card
        go to 110
  130   ENDFILE is1
      endif
c
c  Formats
c
 9010 format(4(0p,f6.3,1p,e14.8))
c
      return
      end
      subroutine TPUN(ip,ne,ipl,ipot,ipass)
c***********************************************************************
c  Print: - potential parameters                                       *
c         - cross-sections                                             *
c         - transmission coefficients in GNASH format                  *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const1.i'
      INCLUDE 'const2.i'
      INCLUDE 'ener.i'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
      INCLUDE 'tce.i'
      INCLUDE 'tlj.i'
      INCLUDE 'xs.i'
c
      integer          ip,ne,ipl,ipot,ipass
c
      character*80     card
      integer          iat,izt,izat
      integer          i,ic1,ic2,ip1,iu,j,ju,k,l,n,nl
c
      double precision zero
      data             zero /0.0d+00/
c=======================================================================
c
      nl  = 20
      k   = 0
      ip1 = 6/ipl
c
c  First pass
c
      if(ipass .eq. 0) then
        do j=1,nl,ip1
          ju = j + ip1 - 1
          ju = MIN0(ju,nl)
          write(is2,9010) ((t(n,l),n=1,ipl),l=j,ju)
        end do
        return
      endif
c
c  Second pass
c
      REWIND is2
c
      iat  = IDNINT(mt)
      izt  = IDNINT(zt)
      izat = 1000*izt + iat
      k = k + 1
      write(is1,9110) izat,ip,ipot,ne,nl,ipl,k
c
      do i=1,nptypx
        k = k + 1
        write(is1,9120) (pot(i,1,j),j=1,6),k
c
        k = k + 1
        write(is1,9130) a(i),ae(i),r(i),re(i),k
      end do
      k = k + 1
      write(is1,9140) rcoul,zero,zero,beta,k
c
      do i=1,ne
        k = k + 1
        write(is1,9130) ein(i),se(i),sr(i),st(i),k
      end do
c
      k = k + 1
      write(is1,9150) k
c
      do i=1,ne,6
        k  = k + 1
        iu = i + 5
        iu = MIN0(iu,ne)
        card = ' '
        ic2 = 0
        do j=i,iu
          ic1 = ic2 +  1
          ic2 = ic1 + 11
          write(card(ic1:ic2),fmt='(1p,e12.5)') ein(j)
        end do
        write(card(73:80),fmt='(i8)') k
        write(is1,fmt='(a)') card
      end do
c
      do i=1,ne
        do j=1,nl,ip1
          k  = k + 1
          ju = j + ip1 - 1
          ju = MIN0(ju,nl)
          read(is2,9010) ((t(n,l),n=1,ipl),l=j,ju)
          card = ' '
          ic2 = 0
          do l=j,ju
            do n=1,ipl
              ic1 = ic2 +  1
              ic2 = ic1 + 11
              write(card(ic1:ic2),fmt='(1p,e12.5)') t(n,l)
            end do
          end do
          write(card(73:80),fmt='(i8)') k
          write(is1,fmt='(a)') card
        end do
      end do
c
      REWIND is2
c
c  Formats
c
 9010 format(1p,6e12.5)
 9110 format(6i12,i8)
 9120 format(1p,6e12.5,    0p,i8)
 9130 format(1p,4e12.5,24x,0p,i8)
 9140 format(1p,4e12.5,24x,0p,i8)
 9150 format(72x,i8)
c
      return
      end
      subroutine VOLIN(itype,depth,rad,dif,volint)
c***********************************************************************
c  Calculate volume integrals                                          *
c----------------------------------------------------------------------*
c  ITYPE  = 1 or -1 Woods-Saxon            potential                   *
c           2 or -2 Woods-Saxon derivative potential                   *
c  ITYPE  > 0       calculate VOLINT from DEPTH, RAD and DIF           *
c         < 0       calculate DEPTH from VOLINT, RAD and DIF           *
c  DEPTH  = depth       of the potential                               *
c  RAD    = radius      of the potential                               *
c  DIF    = diffuseness of the potential                               *
c  VOLINT = volume integral                                            *
c----------------------------------------------------------------------*
c  Called by DISPER, SCAT                                              *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'const2.i'
c
      double precision depth,rad,dif,volint
      integer          itype
c
      double precision arg
c
      double precision zero,one,three,four
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             three /3.0d+00/
      data             four  /4.0d+00/
c=======================================================================
c
      if(rad .eq. zero) return
c
      arg = (pi*dif/rad)**2
c
      if(itype .eq. 1) then
        volint = four*pi*(rad**3)*depth*(one + arg) / three
      else if(itype .eq. -1) then
        depth  = volint / (four*pi*(rad**3)*(one + arg) / three)
      else if(itype .eq.  2) then
        volint = four*four*pi*dif*(rad**2)*depth*(one + arg/three)
      else if(itype .eq. -2) then
        depth  = volint / (four*four*pi*dif*(rad**2)*(one + arg/three))
      endif
c
      return
      end
      subroutine WCHK(i)
c***********************************************************************
c  Check imaginary potential continuity                                *
c----------------------------------------------------------------------*
c  Called by MAIN                                                      *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
c
      integer          i
c
      double precision el
      double precision depth1,depth2,diff
      integer          nz,nz1
c
      double precision zero,eps
      data             zero /0.0d+00/
      data             eps  /1.0d-03/
c=======================================================================
c
      do nz=1,npzen(i)-1
        el = epot(i,nz)
        if(pot(i,nz,7) .eq. zero) then
          depth1 = pot(i,nz,1) + pot(i,nz,2)*el + pot(i,nz,3)*el*el  +
     *             pot(i,nz,4)*el*el*el       + pot(i,nz,5)*DLOG(el) +
     *             pot(i,nz,6)*DSQRT(el)
        else
          depth1 =
     *         pot(i,nz,7)*DEXP(pot(i,nz,2)*(el - eferm)/pot(i,nz,7))
        endif
c
        nz1 = nz + 1
        if(pot(i,nz1,7) .eq. zero) then
          depth2 = pot(i,nz1,1) + pot(i,nz1,2)*el + pot(i,nz1,3)*el*el +
     *             pot(i,nz1,4)*el*el*el       + pot(i,nz1,5)*DLOG(el) +
     *             pot(i,nz1,6)*DSQRT(el)
        else
          depth2 =
     *         pot(i,nz1,7)*DEXP(pot(i,nz1,2)*(el - eferm)/pot(i,nz1,7))
        endif
c
        if(depth1 .gt. eps .and. depth2 .gt. eps) then
          diff = (depth1 - depth2) / depth1
          if(DABS(diff) .gt. eps) then
            write(is,9100) i,el
          endif
        endif
      end do
c
c  Formats
c
 9100 format(' ','POTENTIAL # ',i1
     *          ,' DISCONTINUITY AT E = ',f8.4,' MeV')
c
      return
      end
      double precision function RMSSUR(rad,dif)
c***********************************************************************
c  Calculate the mean-square radius for a surface potential            *
c----------------------------------------------------------------------*
c  Called by RMS                                                       *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'const2.i'
c
      double precision rad,dif
c
      double precision arg,arg2,rad2
c
      double precision one,two,three,seven,fifte
      data             one    /1.0d+00/
      data             two    /2.0d+00/
      data             three  /3.0d+00/
      data             seven  /7.0d+00/
      data             fifte  /1.5d+01/
c=======================================================================
c
      arg    = (pi*dif/rad)**2
      arg2   = arg*arg
      rad2   = rad*rad
c
      rmssur = rad2*(one + two*arg + seven*arg2/fifte) /
     *         (one + arg/three)
c
      return
      end
      subroutine DISPIN1(i,el,deltv)
c***********************************************************************
c  Calculate dispersion integrals for potentials whose energy          *
c  dependence is given by piecewise linear segments                    *
c----------------------------------------------------------------------*
c  I     = type of potential                                           *
c               3  imaginary volume                                    *
c               4  imaginary surface                                   *
c  EL    = lab energy                                                  *
c  DELTV = corrective term                                             *
c----------------------------------------------------------------------*
c  Called by DISPER                                                    *
c***********************************************************************
c
      implicit none
c
      INCLUDE 'scat2000.parameter'
      INCLUDE 'const2.i'
      INCLUDE 'inout.i'
      INCLUDE 'poten1.i'
      INCLUDE 'poten2.i'
      INCLUDE 'poten3.i'
c
      double precision el,deltv
      integer          i
c
      double precision term,z
      integer          icas,n
c
      double precision zero,one,two
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             two   /2.0d+00/
c=======================================================================
c
      deltv     = zero
      epot(i,0) = eferm
c
      z         = el - two*eferm
c
      do n=1,npzen(i)
c
c  The singularity is at the lower bound of the interval
        if(el .eq. epot(i,n-1)) then
          if(epot(i,n) .eq. epotmx .and. pot(i,n,2) .eq. zero) then
            term =  pot(i,n,1)*
     *              DLOG(DABS( (epot(i,n-1) + z)                    ))
            icas =  11
          else
            term = (pot(i,n,1) + pot(i,n,2)*el)*
     *              DLOG(DABS( (epot(i,n) - el)                    ))
     *           - (pot(i,n,1) - pot(i,n,2)*z )*
     *              DLOG(DABS( (epot(i,n) + z )/(epot(i,n-1) + z ) ))
            icas =  12
          endif
c
c  The singularity is at the upper bound of the interval
        else if(el .eq. epot(i,n)) then
          if(epot(i,n) .eq. epotmx .and. pot(i,n,2) .eq. zero) then
            term =  pot(i,n,1)*
     *              DLOG(DABS( (epot(i,n-1) + z)/(epot(i,n-1) - el) ))
            icas =  21
          else
            term = (pot(i,n,1) + pot(i,n,2)*el)*
     *              DLOG(DABS( one             /(epot(i,n-1) - el) ))
     *           - (pot(i,n,1) - pot(i,n,2)*z )*
     *              DLOG(DABS( (epot(i,n) + z )/(epot(i,n-1) + z ) ))
            icas =  22
          endif
c
c  The singularity is inside or outside the interval
        else
          if(epot(i,n) .eq. epotmx .and. pot(i,n,2) .eq. zero) then
            term =  pot(i,n,1)*
     *              DLOG(DABS( (epot(i,n-1) + z)/(epot(i,n-1) - el) ))
            icas =  31
          else
            term = (pot(i,n,1) + pot(i,n,2)*el)*
     *              DLOG(DABS( (epot(i,n) - el)/(epot(i,n-1) - el) ))
     *           - (pot(i,n,1) - pot(i,n,2)*z )*
     *              DLOG(DABS( (epot(i,n) + z )/(epot(i,n-1) + z ) ))
            icas =  32
          endif
        endif
        deltv = deltv + term
      end do
c
      deltv   = deltv / pi
c
      return
      end


      subroutine energy2(za,exactmas,spin,parity,excess,amu0c2)
c
c     Routine energy2 looks up values of ground-state mass excess
c     (in MeV), spin, and parity. Missing data produces a flag.
c     Note that negative za omits Duflo for missing data.
c
      common/xener/nmass,inpgrd,lza(9200),ener(9200),spnpar(9200)
      character*4 bcd
      character*1 parity
      dimension bcd(120)
      real*8 amu0c2

      data k13/13/
c
  1   format(28h0***** ground-state data for,i6,19h not in table *****)
  5   format(5x,'+++++ ground state of ',f6.0,' is incompletely describe
     +d,     spin,parity = ',f6.2,a1,2x,'+++++')
  6   format(5x,'+++++',28x,'assignments changed to,   spin,parity =
     +',f6.2,a1,2x,'+++++')
  8   format(///' M A S S   D A T A   I N P U T   T O   S U B R O U T I
     +N E   E N E R G Y',/(20a4))
c
c     FIRST CALL CAUSES DATA TO BE READ IN FROM UNIT 13
      if(inpgrd.eq.12345) go to 10
      open(unit=13,status='unknown',file='gs-mass-sp.dat')
      open(unit=44,status='unknown',file='massinfo.out')
      read(k13,'(20a4)')bcd
      write(44,8)bcd
      read(k13,'(i7)')nmass
      read(k13,'(5(i7,f11.6,f7.1))')(lza(n),ener(n),spnpar(n),n=1,nmass)
      rewind k13
      inpgrd = 12345
      close (unit=k13)
c
c     SET FLAG FOR DUFLO OPTION
 10   imiss=1
      if(za) 12,15,20
 12   za=abs(za)
      imiss=0
      go to 20
c
c     Z=0,A=0 IS CONSIDERED A PHOTON
   15 excess=0.
      spin=0.
      parity='-'
      exactmas=0.
      return
c
c     FIND REQUESTED NUCLEUS IN TABLE
   20 iza=ifix(za)
      ia=mod(iza,1000)
      amass=float(ia)
      if(iza.lt.lza(1).or.iza.gt.lza(nmass))go to 40
      in=1
      if(iza.eq.lza(1))go to 50
      il=1
      iu=nmass
 122  ii=(il+iu)/2
      if(lza(ii)-iza)124,130,125
 124  il=ii
      go to 126
 125  iu=ii
 126  if(iu-il-1)135,140,122
 130  in=ii
      go to 50
 135  in=il
      go to 150
 140  in=iu
 150  if(lza(in).ne.iza)go to 40
      go to 50
c
c     REQUESTED ISOTOPE IS NOT IN TABLES
 40   write(44,1) iza
      excess= -999999.
      if(imiss.eq.0) stop 7776
      npro=iza/1000
      nneu=ia-npro
      call duflo(nneu,npro,excess,amu0c2)
      exactmas=amass+excess/amu0c2
      go to 200
c
c     PREPARE RETRIEVED DATA FOR EXIT
 50   continue
      excess=ener(in)
      exactmas=amass+excess/amu0c2
      if(spnpar(in).ge.9900.)spin=spnpar(in)-9900.
      if(spnpar(in).ge.9900.)parity='u'
      if(spnpar(in).ge.9900.)go to 200
      if(spnpar(in).ge.100.)parity='+'
      if(spnpar(in).ge.100.)spin=spnpar(in)-100.
      if(spnpar(in).lt.0.)parity='-'
      if(spnpar(in).lt.0.)spin=spnpar(in)+100.
 200  kspin=spin+0.01
      if((parity.ne.'u').and.(kspin.ne.99))return
      write(44,5) za,spin,parity
      if(parity.eq.'u')parity='+'
      if(kspin.eq.99)spin=.25*(1.-(-1.)**ia)
      write(44,6) spin,parity
      return
      end

      subroutine duflo(nn,nz,excess,amu0c2)
c
c                --------------------------------------------
c                |Nuclear mass formula of Duflo-Zuker (1992)|
c                --------------------------------------------
c
      real*8 amu0c2
      dimension xmag(6),zmag(5),a(21)
      data zmag/14.,28.,50.,82.,114./
      data xmag/14.,28.,50.,82.,126.,184./
      data a/16.178,18.422,120.146,202.305,12.454,0.73598,5.204,1.0645,
     +1.4206,0.0548,0.1786,.6181,.0988,.0265,-.1537,.3113,-.6650,-.0553,
     +-.0401,.1774,.4523/
c
      x=nn
      z=nz
      t=abs(x-z)*.5
      v=x+z
      s=v**(2./3.)
      u=v**(1./3.)
c     a5=a(5)
c     if(z.gt.x) a5=0.
      E0=a(1)*v-a(2)*s-a(3)*t*t/v+a(4)*t*t/u/v-a(5)*t/u-a(6)*z*z/u
      esh=0.
      esh1=0.
      do 10 k=2,5
        f1=zmag(k-1)
        f2=zmag(k)
        dfz=f2-f1
        if(z.ge.f1.and.z.lt.f2) then
          roz=(z-f1)/dfz
          pz=roz*(roz-1)*dfz
          do 20 l=2,6
            f3=xmag(l-1)
            f4=xmag(l)
            dfn=f4-f3
            if(x.ge.f3.and.x.lt.f4) then
              ron=(x-f3)/dfn
              pn=ron*(ron-1)*dfn
              esh=(pn+pz)*a(8)+a(10)*pn*pz
              xx=2.*ron-1.
              zz=2.*roz-1.
              txxx=pn*xx
              tzzz=pz*zz
              txxz=pn*zz
              tzzx=pz*xx
              kl=l-k
              if(kl.eq.0) esh1=a(k+10)*(txxx+tzzz)+a(k+15)*(txxz+tzzx)
              if(kl.eq.1)
     +          esh1=a(k+11)*txxx-a(k+16)*txxz+a(k+10)*tzzz-a(k+15)*tzzx
              if(kl.eq.2)
     +          esh1=a(k+12)*txxx+a(k+17)*txxz+a(k+10)*tzzz+a(k+15)*tzzx
                edef=a(9)*(pn+pz)+a(11)*pn*pz
                if(esh.lt.edef) esh=edef
            end if
   20     continue
        end if
   10 continue
      ebin=e0+esh+esh1
      nn2=nn/2
      nz2=nz/2
      nn2=2*nn2
      nz2=2*nz2
      if(nn2.ne.nn) ebin=ebin-a(7)/u
      if(nz2.ne.nz) ebin=ebin-a(7)/u
      rneumas=1.008665
      rpromas=1.007825
      amass=nz*rpromas+nn*rneumas-ebin/amu0c2
      excess=(amass-nz-nn)*amu0c2
      return
      end
