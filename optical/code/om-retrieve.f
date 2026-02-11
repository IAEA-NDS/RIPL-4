      program omretrieve
c
c     Code to retrieve optical model potentials from the RIPL optical
c     model potential (omp) library and to format them for input into
c     the SCAT2000 and ECIS03(ECIS96) computer codes.
c------------------------------------------------------------------------
c     New
c
c     Version date: September 23 2004, RCN
c
c     1. Format expanded to allow for:
c        -Energy dependent radius up to E**2    (r(i,j,12))
c        -Dispersive-like energy dependent radius ( r(i,j,13)>0 )
c        -Soukhovitski et al energy dependence of the real potential
c         See J.Phys.G. Nucl.Part.Phys. 30(2004) 905-920)
c        -Buck and Perey energy dependence of the local real potential
c         (As coded and employed by B. Morillon and P. Romain,
C         See Phys.Rev.C70(2004) 014601)
c
c     2. Analytical dispersive integrals are included
c         See Quesada JM et al, Comp. Phys. Comm. 153(2003) 97
c                               Phys. Rev. C67(2003) 067601
c
c     3. General numerical solution of the dispersive integral is
c        kept to deal with special Ws(E) and Wd(E) form factors.
c         See Capote et al, J.Phys.G. Nucl.Part.Phys. 27(2001) B15-B19
c
c        Numerical integration is used for Nagadi et al dispersive OMP
c         See Nagadi et al, Phys. Rev. C68 (2003) 044610
c
c     4. Relativistic kinematics is used to calculate the kinematical
c        conversion factor xkine and reduced mass amu (for irel>0)
c        Kinematical conversion factor for discrete levels is assumed
c        non-relativistic.
c
c     5. The isimple=1 option to send energy block to SCAT code was
c        eliminated (it could be reproduced by sending energies one
c        by one. The speed gain is negligible today)
c
c     Previous Version date:  13 March 2003, RCN
c
c------------------------------------------------------------------------
c
c     P.G.Young, Group T-2, Los Alamos National Laboratory
c     Mail Stop B243, email address pgy@lanl.gov
c
c     R. Capote Noy, IAEA Nuclear Data Section
c     email address R.CapoteNoy@iaea.org or rcapotenoy@yahoo.com
c
c     Code compiled using:
c         f77 om-retrieve.f -C -O3 -o omretrieve
c     That is, -r8 should not be turned on in order for the dispersive
c     corrections to be computed accurately.
c
c     INPUT FILES
c       omp-parameter-u.dat = RIPL3 Optical Model Potential (OMP) file
c       gs-mass-sp.dat = ground-state mass, spin-parity file
c       ominput.inp = input instructions, defined below
c
c     OUTPUT FILES
c       sc2.inp   = input file for the scat2000 code (modtyp=1)
c       ecis.inp  = standard input file for the ecis03 code (modtyp=2)
c       ecistc.inp= alternate input file for the ecis03 code.
c                   (modtyp=2).  Useful for generating transmission
c                   transmission coefficients.
c       ecisvib   = input file for the ecis03 code with vibrational
c                   model activated (modtyp=2).
c       ecisdw.inp= input file for the ecis03 code with DWBA
c                   model activated (modtyp=3).
c       massinfo.out = descriptive remarks about the ground-state mass
c                      data used in the gs-mass-sp.dat file.
c       omp-table.dat = a table of numerical values of each potential made
c                   at each incident energy when modtyp=2 (ECIS input).
c
c     INPUT DETAILS (ominput.inp file)
c
c     read(5,*) ne
c
c       ne = number of incident energies
c          = 0 to use a built-in array of incident energies
c
c     if(ne.gt.0) read(5,*) (en(n),n=1,ne)
c
c       en(n) = incident energies in MeV in the laboratory system
c
c     read(5,*,end=990) iztar,iatar,irefget,modtyp
c
c       abs(iztar) = Z of target nucleus
c                    Set iztar negative to provide integer projectile
c                    mass in input decks.
c       abs(iatar) = A of target nucleus
c                    Set iatar negative to provide integer target
c                    mass in input decks.
c
c       irefget = reference number of optical model potential to be
c                 retrieved from the RIPL library
c         Set irefget = negative izaproj (projectile) to retrieve all
c                       spherical potentials in the RIPL library for
c                       this izaproj and the inputted izatar (target).
c
c       iabs(modtyp)= 1 to generate SCAT2 input file (sc2.inp)
c       iabs(modtyp)= 2 to generate ECIS03 input files (ecis.inp,
c                       ecistc.inp)
c       iabs(modtyp)= 3 to generate ECIS03 DWBA input files
c                       (ecisdw.inp), using structure information from
c                       an external file (deform.dat), which may be the
c                       om-deformations.dat*file or a user-provided
c                       deform.dat file.
c
c     if(iabs(modtyp).eq.3)read(5,*)kdef
c
c       Set modtyp negative to force imaginary volume and surface
c       potentials to be zero or positive.  If iabs(modtyp)=1, then
c       all SCAT2000 inputs are single energy, i.e., the compact
c       energy representation of SCAT2000 is bypassed.
c
c
c         kdef=1 use JENDL-3.2 betas from om-deformations.dat file
c             =2 use ENSDF(Q)
c             =3 use ENSDF(BE2) and ENSDF(BE3)
c             =4 use Raman and Spear for BE2,BE3
c             =5 use deform.dat file provided by user
c
c       See om-deformations.readme file for more information on the
c       kdef options and for the format of deform.dat for using kdef=5.
c       For kdef=1-4, the om-deformations.dat file must be copied to
c       deform.dat.  Note that where multiple entries occur for a single
c       state in om-deformations.dat, the code makes an ECIS96 input
c       for each entry.
c
      include "om-retrieve.cmb"
c
      character*3 ifin
      data isuit,ifin/0,'FIN'/
      data nee/21/
c
      open(unit=1,file='om-parameter-u.dat')
      open(unit=5,file='ominput.inp')
      open(unit=35,file='omp-table.dat')     
      open(unit=37,file='vint-table.dat')
c
c     Read input information
c
      nsc2=0
      read(5,*) ne
      if(ne.gt.0) go to 88
      nee=0
      go to 90
  88  read(5,*) (en(n),n=1,ne)
c
c     Statement 90 ==> Main Read-Cycle Return
c
      nread = 0
  90  read(5,*,end=990) iztar,iatar,irefget,modtyp
      nonegw=0
      if(modtyp.lt.0) nonegw=1
      modtyp=iabs(modtyp)
      if(modtyp.eq.3)read(5,*)kdef
c
      irndp=0
      if(iztar.lt.0)irndp=1
      iztar=iabs(iztar)
      irndt=0
      if(iatar.lt.0)irndt=1
      iatar=iabs(iatar)
      izatar=1000*iztar+iatar
c
c     output unit
c
      ko=15
      if(modtyp.eq.1) open(unit=15,file='sc2.inp')
      if(modtyp.eq.2) open(unit=15,file='ecis.inp')
      if(modtyp.eq.3) open(unit=15,file='ecisdw.inp')
      open(unit=16,file='ecistc.inp')
c
c     Get optical model parameters
c
  98  rewind ki
 100  call omin11
      if(ieof.eq.1) go to 90

      if(idr.GT.0) open(unit=36,file='DOM-table.dat')
      if(iztar.GT.0) open(unit=38,file='coul-vint-table.dat')
    
c
c     Set energy grid up to emax if nee=0
c
      if(nee.eq.0)call egrid(ne,emin,emax,en)
      call optname
c
c     Check for desired reference and Z,A range
 108  if(irefget.lt.0.and.imodel.lt.1) go to 110
      if(irefget.eq.iref) go to 112
      go to 100
 110  izapr=iabs(irefget)
      izpr=izapr/1000
      iapr=mod(izapr,1000)
      if(izpr.ne.izproj.or.iapr.ne.iaproj) go to 100
 112  if(izmin.eq.0.and.izmax.eq.0) go to 114
      if(iztar.lt.izmin.or.iztar.gt.izmax)
     + write(6,'("Inputted target Z outside range of RIPL potential")')
 114  if(iamin.eq.0.and.iamax.eq.0) go to 120
      if(iatar.lt.iamin.or.iatar.gt.iamax)
     + write(6,'("Inputted target A outside range of RIPL potential")')
c
c     Check if input file can be made for this potential
c
 120  if(imodel.le.2)go to 122
      write(6,'("imodel=",i2," in RIPL library not implemented yet.")')
      go to 1020
c
 122  if(modtyp.ne.1)go to 126
      if(imodel.eq.0)go to 126
      write(6,'("Coupled-channel model not implemented for SCAT2000")')
      go to 1020

 126  if(imodel.eq.2.and.modtyp.eq.3)write(6,'("You are making an ECIS03
     + DWBA input from a vibrational potential.")')
c
      if(imodel.eq.1.and.modtyp.eq.3)write(6,'("You are making an ECIS03
     + DWBA input from a rotational potential.")')
c
c     Initialize program
 130  call setup
      kecis=0

c     Call output routines
      if(modtyp.ne.1)go to 150

c     Send one energy at a time to scatip2
      isimple=2
      nesc=1
      do 146 n=1,ne
      ensc(1)=en(n)
      call scatip2
 146  continue
      go to 200
c
 150  continue
c
c     Check for Gaussian form factor for WD
      krange=iabs(jrange(4))
      if(krange.eq.0)go to 153
      do 151 j=1,krange
      if(rco(4,j,1).lt.0.)go to 152
 151  continue
      go to 153
 152  write(6,'("No ECIS03 inputs made when Gaussian form factor ",
     + "used for WD",/,"iref=",i5)')iref
c
 153  if(modtyp.ne.2)go to 160
c
c     Set up for table print
      call tableset
c
      if(imodel.eq.2)go to 154
c
c     rotational model being used
      kecis=1
      ko=15
      call ecisip
      kecis=2
      ko=16
      call ecisip
      go to 200
c
c     vibrational model being used
 154  kecis=3
      call ecisdwip
      go to 200
 160  if(modtyp.ne.3) stop 5
c
c     DWBA model being used
      kecis=4
      call ecisdwip
c
 200  if(irefget.le.0) go to 100
      go to 90
c
c     End computation
 990  if(modtyp.eq.1) write(15,'(i5)')isuit
      if(modtyp.eq.2) write(15,'(a3)') ifin
      if(modtyp.eq.3) write(15,'(a3)') ifin
      close(15)
      backspace (16,ERR=1010)
      read(16,*,END=1010)
      write(16,'(a3)') ifin
      close(16)
      if(idr.GT.0) close(36)
	close(35)
	close(37)
      if(iztar.GT.0) close(38)
 1020 stop
 1010 close(16,status='delete')
      end program omretrieve

      subroutine omin11
c
c     routine to retrieve optical model parameters from RIPL library
      include "om-retrieve.cmb"
c
  1   format(80a1)
c
c     Zero arrays
      call setr(0.,epot,5*ndim1)
      call setr(0.,rco,5*ndim1*ndim2)
      call setr(0.,aco,5*ndim1*ndim2)
      call setr(0.,pot,5*ndim1*ndim3)
c
      ieof=0
      read(ki,*,end=999) iref
      read(ki,1) (author(m),m=1,80)
      read(ki,1) (refer(m),m=1,80)
      read(ki,1) (summary(m),m=1,320)
      read(ki,*) emin,emax
      read(ki,*) izmin,izmax
      read(ki,*) iamin,iamax
      read(ki,*) imodel,izproj,iaproj,irel,idr
c
      do 100 m=1,6
        read(ki,*) jrange(m)
        if(jrange(m).eq.0) go to 100
        krange=iabs(jrange(m))
        do 98 j=1,krange
          read(ki,*) epot(m,j)
          read(ki,*) (rco(m,j,n),n=1,ndim2)
          read(ki,*) (aco(m,j,n),n=1,ndim2)
          read(ki,*) (pot(m,j,n),n=1,ndim3)
  98    continue
 100  continue
      read(ki,*) jcoul
      if(jcoul.le.0) go to 110
      do 108 j=1,jcoul
 108  read(ki,*) ecoul(j),rcoul0(j),rcoul(j),rcoul1(j),rcoul2(j),
     >            beta(j),acoul(j)
 110  if(imodel.ne.1) go to 130
      read(ki,*) nisotop
      do 120 n=1,nisotop
        read(ki,*) iz(n),ia(n),ncoll(n),lmax(n),idef(n),bandk(n),
     +    (def(n,k),k=2,idef(n),2)
        do 124 k=1,ncoll(n)
          read(ki,*) ex(k,n),spin(k,n),ipar(k,n)
 124    continue
 120  continue
      go to 200
 130  if(imodel.ne.2) go to 150
      read(ki,*) nisotop
      do 140 n=1,nisotop
        read(ki,*) iz(n),ia(n),nvib(n)
        do 138 k=1,nvib(n)
          read(ki,*) exv(k,n),spinv(k,n),iparv(k,n),nph(k,n),defv(k,n),
     +    thetm(k,n)
 138    continue
 140  continue
      go to 200
 150  if(imodel.ne.3)go to 200
      read(ki,*) nisotop
      do 160 n=1,nisotop
        read(ki,*) iz(n),ia(n),beta0(n),gamma0(n),xmubeta(n)
 160  continue
 200  continue
      read(ki,1,end=999)idum
      return
 999  ieof=1
      write(6,*)'End of RIPL-2 optical potential library reached.'
      return
      end subroutine omin11

      subroutine setr (a,b,n)
c     ******************************************************************
c     set all elements of array b(n) to real number a.
c     ******************************************************************
      dimension b(n)
      do 100 k=1,n
  100 b(k)=a
      return
      end subroutine setr

      subroutine setup
c
c     Initialize program
      include "om-retrieve.cmb"
c
      becis1='FFFFFFFFFFFFFFFFFFFFFFFFFFFTFFFFFFFFFFFFFFFFFFFFFF'
      becis2='FFFFFFFFFFFFFTFFTTTFFTTFTFFFFFFFFFFFFFFFFFFFFFFFFF'
c
      ndisx=ndim8
      ifirst1=0
      ifirst2=0
      ifirst3=0
      ifirst4=0
      ifirst5=0
c
      izaproj=1000*izproj + iaproj
      if(izaproj.eq.1   ) ipsc=1
      if(izaproj.eq.1001) ipsc=2
      if(izaproj.eq.1002) ipsc=3
      if(izaproj.eq.1003) ipsc=4
      if(izaproj.eq.2003) ipsc=5
      if(izaproj.eq.2004) ipsc=6
      zproj=real(izproj)
      atar=float(iatar)
      ztar=float(iztar)
      eta=1.-(2.*ztar/atar)
      encoul=0.4*ztar/atar**(1./3.)
c
c     Get Fermi energy, masses and kinematic factor (amu)
      call masses
      if(irndt.eq.1) tarmas=atar
      if(irndp.eq.1) projmas=iaproj
c
c     Internal Flags
      koptom=2
c     Set koptom=1 to compute excited state potentials at incident
c         energy only. (Shortcut that is not recommended.)
c     Set koptom=2 to compute excited state potentials at incident
c         energy minus the excitation energy (corrected to lab).
c         This is the recommended option.
c
      return
      end subroutine setup

c*******************
c     SCAT2000 input
c*******************
      subroutine scatip2
c     routine to write optical model parameters in parameter lib
      include "om-retrieve.cmb"
      include "om-retrieve-cts.cmb"
c
      parameter(npcofx=7)
      dimension rsc(6),resc(6),asc(6),aesc(6),potsc(6,ndim1,npcofx),
     + npzen(6),epotsc(6,ndim1),b(6,ndim1,15)
c
c     To use dispersive optical model package
c
c             functions
      real*8 DOM_int_Ws,DOM_int_Wv,DOM_int_T1,DOM_int_T2,Vhf
      real*8 DOM_int,WVf,WDf,Delta_WD,Delta_WV
c             variables
      real*8 As,Bs,Cs,AAv,Bv,AAvso,Bvso,EEE,Ep,Ea,Ef
      real*8 alpha_PB,beta_PB,gamma_PB,Vnonl
      real*8 DWS,DWV,DWVso,DerDWV,DerDWS,dtmp
      real*8 WDE,WVE
      integer n,iq,nns
c     Common blocks for numerical integration
      common /energy/EEE,Ef,Ep
      common /Wenerg/WDE,WVE
      common /pdatas/As,Bs,Cs,nns,iq
      common /pdatav/AAv,Bv,n

c     Nonlocality constants alpha fixed according to Mahaux 1991
      real*8 AlphaV
      data AlphaV/1.65d0/
c     real*8 AlphaS
c     data AlphaS/0.235d0/
c
      data isuit/1/
      data ipr,ida,iba,ipu/1,0,0,1/

      external Delta_WV,WVf,Delta_WD,WDf
c
c     if(ifirst1.ne.25)open(unit=4,file='sc2.inp')
c     ifirst1=25
c     ko=4
c
      call setr(0.,rsc,6)
      call setr(0.,resc,6)
      call setr(0.,asc,6)
      call setr(0.,aesc,6)
      call setr(0.,potsc,6*npcofx*ndim1)
c

      VDcoul = 0.d0
      el=ensc(1)

c     Get coulomb radius, if needed
      rcoulsc=0.
      if(jcoul.lt.1)go to 100
      jc=1
c     if(isimple.ne.2)go to 90
      do 80 j=1,jcoul
      if(el.gt.ecoul(j))jc=j+1
  80  continue
      jc=min0(jc,jcoul)
  90  rcoulsc=rcoul(jc) + rcoul0(jc)*atar**(-1./3.) +
     +       rcoul1(jc)*atar**(-2./3.) + rcoul2(jc)*atar**(-5./3.)
      betasc=beta(jc)
 100  encoul2=0.
      if(rcoulsc.gt.0.) encoul2=1.73*ztar/(rcoulsc*atar**(1./3.))
      ipot=0


      do 150 i=1,6
      jab=iabs(jrange(i))
      ii=i
      if(i.eq.2)ii=3
      if(i.eq.3)ii=2
      npzen(ii)=jrange(i)
      if(jrange(i).eq.0)go to 150
      do 104 j=1,jab
 104  epotsc(ii,j)=epot(i,j)

      vc = 0.d0
      VVcoul = 0.d0
      VScoul = 0.d0
      DerDWV = 0.d0
      DerDWS = 0.d0
c
c     Find energy range that energy el falls in.
 110  jp=1
      do 112 j=1,jab
      if(el.gt.epot(i,j))jp=j+1
 112  continue
      j=min0(jp,jab)
      js=1
      epotsc(ii,js)=epot(i,j)
      if(npzen(ii).gt.0)npzen(ii)=+1
      if(npzen(ii).lt.0)npzen(ii)=-1
c
      Ef=dble(int(100000*efermi))/100000
      if(pot(i,j,18).ne.0.) Ef=pot(i,j,18) + pot(i,j,19)*atar
      elf = el - Ef
c
c     Calculate radius and diffuseness parameters
      if(rco(i,j,13).eq.0.) then
        rsc(ii)=abs(rco(i,j,1)) + rco(i,j,3)*eta
     *       + rco(i,j,4)/atar + rco(i,j,5)/sqrt(atar)
     *       + rco(i,j,6)*atar**(2./3.) + rco(i,j,7)*atar
     *       + rco(i,j,8)*atar**2  + rco(i,j,9)*atar**3
     *       + rco(i,j,10)*atar**(1./3.)
     *       + rco(i,j,11)*atar**(-1./3.)
C--------------------------------------------------------------------
C     RCN, 08/2004, to handle new extension to the OMP RIPL-2 format
     *       + rco(i,j,2)*el
     *       + rco(i,j,12)*el*el
      else
C     RCN, 09/2004, to handle new extension to the OMP RIPL-2 format
        nn = int(rco(i,j,7))
        rsc(ii)= ( abs(rco(i,j,1)) + rco(i,j,2)*atar ) *
     *           ( 1.d0 - ( rco(i,j,3) + rco(i,j,4)*atar ) * elf**nn/
     *           ( elf**nn + ( rco(i,j,5) + rco(i,j,6)*atar )**nn ) )
      endif
      resc(ii)=0.
C--------------------------------------------------------------------
      asc(ii)=abs(aco(i,j,1)) + aco(i,j,3)*eta
     *       + aco(i,j,4)/atar + aco(i,j,5)/sqrt(atar)
     *       + aco(i,j,6)*atar**(2./3.) + aco(i,j,7)*atar
     *       + aco(i,j,8)*atar**2 + aco(i,j,9)*atar**3
     *       + aco(i,j,10)*atar**(1./3.) + aco(i,j,11)*atar**(-1./3.)
C--------------------------------------------------------------------
C     RCN, 08/2004, to handle new extension to the OMP RIPL-2 format
     *       + aco(i,j,2)*el
      aesc(ii)=0.
C--------------------------------------------------------------------
c
      if(pot(i,j,24).eq.0.) go to 120
c
c     Special Koning-type potential formulas
c
      if (pot(i,j,24).eq.1.) then
c
c       Koning-type formulas
c
        if(i.eq.1) call bcoget(b,j)
c
        elf=el-Ef
        vc=0.
        if(i.eq.1 .and. b(1,j,5).ne.0.) then
          vc = b(1,j,1)*encoul2*( b(1,j,2) - 2.*b(1,j,3)*elf +
     +    3.*b(1,j,4)*elf**2 + b(i,j,14)*b(i,j,13)*exp(-b(i,j,14)*elf) )
	      VDcoul = b(i,j,5)*vc
  	  endif

        nn = int(pot(i,j,13))
c       Retrieving average energy of the particle states Ep
        Ep=Ef
        if( (i.eq.2) .or. (i.eq.4) )
     >     Ep=dble(int(100000*pot(i,j,20)))/100000
        if(Ep.eq.0.) Ep=Ef
        elf = el - Ep

        iq=1
        if(i.eq.4 .and. b(4,j,12).gt.0.) iq=nint(b(4,j,12))

        potsc(ii,js,1)=
     +     b(i,j,1)*( b(i,j,15) - b(i,j,2)*elf + b(i,j,3)*elf**2 -
     +     b(i,j,4)*elf**3 + b(i,j,13)*exp(-b(i,j,14)*elf) ) +
     +     b(i,j,5)*vc + b(i,j,6)*(elf**nn/(elf**nn + b(i,j,7)**nn)) +
     +     b(i,j,8)*exp(-b(i,j,9)*elf**iq)*(elf**nn/
     +     (elf**nn + b(i,j,10)**nn)) + b(i,j,11)*exp(-b(i,j,12)*elf)

      endif

      if (pot(i,j,24).eq.2.) then
c
c       Morillon-Romain formulas
c
        if(i.eq.1) call bcoget(b,j)
c
        elf=el-Ef

        nn = int(pot(i,j,13))
c
c       Vhf(E) calculated from nonlocal approximation
c          as suggested by Perey and Buck
c
        alpha_PB = dble(int(100000*b(i,j,1)))/100000
        beta_PB  = dble(int(100000*b(i,j,2)))/100000
        gamma_PB = dble(int(100000*b(i,j,3)))/100000
        EEE = dble(int(1000000*el))/1000000

        iq=1
        if(i.eq.4 .and. b(4,j,12).gt.0.) iq=nint(b(4,j,12))

        vc=0.
        Vnonl = 0.
        if(i.eq.1) then
	    Vnonl = -Vhf(EEE,alpha_PB,beta_PB,gamma_PB)
C         Numerical derivative of the Vhf
          if(b(1,j,5).ne.0.) then 
	      Vnonlm = -Vhf(EEE-0.05,alpha_PB,beta_PB,gamma_PB)
	      Vnonlp = -Vhf(EEE+0.05,alpha_PB,beta_PB,gamma_PB)
C           Coulomb correction for Hartree-Fock potential
            vc = encoul2*(Vnonlm-Vnonlp)*10.d0
            VDcoul = b(i,j,5)*vc
          endif
        endif
        
        potsc(ii,js,1)=
     +    Vnonl + b(i,j,5)*vc +
     +    b(i,j,6)*(elf**nn/(elf**nn + b(i,j,7)**nn)) +
     +    b(i,j,8)*exp(-b(i,j,9)*elf**iq)*(elf**nn/
     +    (elf**nn + b(i,j,10)**nn)) +
     +    b(i,j,11)*exp(-b(i,j,12)*elf)

      endif

c
c     Nonlocality consideration
c
c     Retrieving energy above which nonlocality in the volume absorptive
c               potential is considered (Ea)
c
      Ea=dble(int(100000*pot(i,j,21)))/100000
      if(Ea.eq.0.) Ea=1000.1d0

      if(i.eq.2 .and. Ea.lt.1000. .and. el.gt.(Ef+Ea) )
     +         potsc(ii,js,1) = potsc(ii,js,1) +
     +       AlphaV*(sqrt(el)+(Ef+Ea)**1.5/(2.*el)-1.5*sqrt(Ef+Ea))
c
      go to 150
c
 120  if(pot(i,j,23).eq.0)go to 130
c
c     Special Varner-type potential formulas
      potsc(ii,js,1)= (pot(i,j,1) + pot(i,j,2)*eta)/
     *    (1.+ exp((pot(i,j,3) - el + pot(i,j,4)*encoul2)/pot(i,j,5)))
      if(pot(i,j,6).eq.0.)go to 150
      potsc(ii,js,1)= potsc(ii,js,1)
     *    + pot(i,j,6)*exp((pot(i,j,7)*el - pot(i,j,8))/pot(i,j,6))
      go to 150
c
 130  if(pot(i,j,22).eq.0)go to 140
c
c     Special Smith-type potential formulas
      pi=acos(-1.d0)
      potsc(ii,js,1)=pot(i,j,1) + pot(i,j,2)*eta
     *    + pot(i,j,6)*exp(pot(i,j,7)*el + pot(i,j,8)*el*el)
     *    + pot(i,j,9)*el*exp(pot(i,j,10)*el**pot(i,j,11))
      if(pot(i,j,5).ne.0.)potsc(ii,js,1)=potsc(ii,js,1)
     *    + pot(i,j,3)*cos(2.*pi*(atar - pot(i,j,4))/pot(i,j,5))
      go to 150
c
c     Standard potential formulas
 140  potsc(ii,js,1)=pot(i,j,1) + pot(i,j,7)*eta + pot(i,j,8)*encoul
     *    + pot(i,j,9)*atar + pot(i,j,10)*atar**(1./3.)
     *    + pot(i,j,11)*atar**(-2./3.) + pot(i,j,12)*encoul2
     *    + (pot(i,j,2) + pot(i,j,13)*eta + pot(i,j,14)*atar)*el
     *    + pot(i,j,3)*el*el + pot(i,j,4)*el*el*el + pot(i,j,6)*sqrt(el)
     *    + (pot(i,j,5) + pot(i,j,15)*eta + pot(i,j,16)*el)*alog(el)
     *    + pot(i,j,17)*encoul/el**2
c
 150  continue
c
      if(rco(4,1,1).lt.0.) rsc(4)=-rsc(4)
      if(nonegw.eq.0) go to 152
      if(potsc(3,js,1).lt.0.0) potsc(3,js,1)=0.
      if(potsc(4,js,1).lt.0.0) potsc(4,js,1)=0.
c
c     To calculate dispersion relation contribution
c
 152  if(abs(idr).ge.2) then
c
c       Exact calculation of the dispersive contribution
c
        EEE = dble(int(1000000*el))/1000000

        i=2
c       Only one energy range
        j=1
c       Real volume contribution from Dispersive relation
        DWV=0.d0
c
        if(pot(2,1,24).ne.0) then

          AAv=dble(int(100000*b(i,j,6)))/100000
          Bv =dble(int(100000*b(i,j,7)))/100000
          n = nint( pot(i,j,13) )
          if(n.eq.0 .or. mod(n,2).eq.1)
     +      stop 'Zero or odd exponent in Wv(E) for dispersive OMP'

c         Retrieving average energy of the particle states Ep
          Ep=dble(int(100000*pot(i,j,20)))/100000
          if(Ep.eq.0.) Ep=Ef

c         analytical DOM integral
          DWV=DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE,n,DerDWV)

c         Coulomb correction for real volume potential 
          DerDWV = -b(1,1,5)*encoul2*DerDWV
C         numerical DOM derivative (not needed for a time being) 
C         DWVp = DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE+0.1d0,n,dtmp)
C         DWVm = DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE-0.1d0,n,dtmp)
C         DerDWV = -b(1,1,5)*encoul2*(DWVp-DWVm)*5.d0
c         if(idr.le.-2) then
c           numerical DOM integral (not needed for a time being)
c           WVE=WVf(AAv,Bv,Ep,Ef,EEE,n)
c           DWV=2*DOM_int(Delta_WV,WVf,Ef,Ef+5.*Bv,150000.d0,EEE,0.d0)
c         endif
c
c         Nonlocality correction to the DOM integral
c           (only used if Ea is non-zero)
          Ea=dble(int(100000*pot(i,j,21)))/100000
          if(Ea.eq.0.) Ea=1000.1d0
          Dwplus = 0.d0
          Dwmin = 0.d0
          T12der = 0.d0  
          if(Ea.lt.1000.) THEN
             Dwplus = AlphaV*DOM_INT_T2(Ef,Ea,EEE)
             dtmp1 = Wvf(AAv,Bv,Ep,Ef,Ef+Ea,n)
             Dwmin = dtmp1*DOM_INT_T1(Ef,Ea,EEE)
             DWV = DWV + Dwplus + Dwmin
c            Coulomb correction for nonlocal dispersive contribution  
c                to real volume potential 
 	       if(b(1,1,5).ne.0.d0) then
               if(eee.ne.0.05d0) then
                 T2p = DOM_INT_T2(Ef,Ea,EEE+0.05d0)  
                 T2m = DOM_INT_T2(Ef,Ea,EEE-0.05d0)  
	 	    T2der = AlphaV*(T2p-T2m)*10.d0 
	           T1p = DOM_INT_T1(Ef,Ea,EEE+0.05d0)
	           T1m = DOM_INT_T1(Ef,Ea,EEE-0.05d0)
	           T1der = dtmp1*(T1p-T1m)*10.d0
	           T12der =  -b(1,1,5)*encoul2* ( T1der + T2der )
		  else
                 T2p = DOM_INT_T2(Ef,Ea,EEE+0.1d0)  
                 T2m = DOM_INT_T2(Ef,Ea,EEE-0.1d0)
		    T2der = AlphaV*(T2p-T2m)*5.d0 
	           T1p = DOM_INT_T1(Ef,Ea,EEE+0.1d0)
	           T1m = DOM_INT_T1(Ef,Ea,EEE-0.1d0)
	           T1der = dtmp1*(T1p-T1m)*5.d0
	           T12der =  -b(1,1,5)*encoul2* ( T1der + T2der )
               endif 
	       endif

          endif
	    VVcoul = DerDWV + T12der 
        endif

        i=4
c       Only one energy range
        j=1
c       Real surface contribution from Dispersive relation
        DWS=0.d0

        if(pot(4,1,24).ne.0) then

          As=dble(int(100000*b(i,j,8)))/100000
          Bs=dble(int(100000*b(i,j,10)))/100000
          Cs=dble(int(100000*b(i,j,9)))/100000
          n = nint( pot(i,j,13) )
          if(n.eq.0 .or. mod(n,2).eq.1)
     +      stop 'Zero or odd exponent in Wd(E) for dispersive OMP'
            iq=1
          if(b(4,j,12).gt.0.) iq=nint(b(4,j,12))

c         Retrieving average energy of the particle states Ep
          Ep=dble(int(100000*pot(i,j,20)))/100000
          if(Ep.eq.0.) Ep=Ef

          if(idr.ge.2) then
c           analytical DOM integral
            DWS = DOM_INT_Ws(Ef,Ep,As,Bs,Cs,EEE,n,DerDWS)
		  VScoul = -b(1,1,5)*encoul2*DerDWS
          endif

          if(idr.le.-2) then
c           numerical DOM integral
            nns=n
            WDE=WDf(As,Bs,Cs,Ep,EEE,n,iq)
            DWS = 2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE,WDE)
c           Coulomb correction for real surface potential 
	      if(b(1,1,5).ne.0.d0) then
              WDE=WDf(As,Bs,Cs,Ep,EEE+0.1d0,n,iq)
              DWSp = 
     >	 	 2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE+0.1d0,WDE)
              WDE=WDf(As,Bs,Cs,Ep,EEE-0.1d0,n,iq)
              DWSm =
     >		 2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE-0.1d0,WDE)
c             Numerical derivative
	        VScoul = -b(1,1,5)*encoul2*(DWSp-DWSm)*5.d0
	      endif
          endif

        endif

        i=6
c       Only one energy range
        j=1
c       Real spin orbit contribution from Dispersive relation
        DWVso=0.d0
c
        if(pot(6,1,24).ne.0 .and. abs(idr).eq.3) then

          AAvso=dble(int(100000*b(i,j,6)))/100000
          Bvso=dble(int(100000*b(i,j,7)))/100000
          n = nint( pot(i,j,13) )
          if(n.eq.0 .or. mod(n,2).eq.1)
     +      stop 'Zero or odd exponent in Wso(E) for dispersive OMP'

c         analytical DOM integral
          DWVso=DOM_INT_Wv(Ef,Ef,AAvso,Bvso,EEE,n,dtmp)

        endif

c       Adding real volume dispersive and Coulomb contribution to the real potential
c       Geometry parameters are the same as for the volume potential(imag and real)
        potsc(1,1,1) = potsc(1,1,1) + DWV + VVcoul
c       Including real surface and Coulomb dispersive contribution
c       Geometry parameters are the same as for the imaginary surface potential
        potsc(2,1,1) = DWS  + VScoul
        rsc(2)=rsc(4)
        asc(2)=asc(4)
        npzen(2)=1
c       Adding real spin orbit dispersive contribution to the real spin orbit potential
c       Geometry parameters are the same as for the imaginary spin orbit potential(imag and real)
        potsc(5,1,1) = potsc(5,1,1) + DWVso

      endif
c
c     Write SCAT2 input
c
      if(nsc2.ge.1)go to 156
 154  itmp=irel
      if(irel.eq.2) itmp=1
c
c     Dispersion relations are considered explicitly in this code
c     so SCAT2 does not need to make any DR related calculation
c         => we are setting idr = 0 always
      write(ko,'(6i5)') ipr,ida,iba,ipu, 0 ,itmp
      nsc2=nsc2+1
 156  write(ko,'(i5)') nesc
      ensc(1)=-ensc(1)
      write(ko,'(7f10.5)') (ensc(j),j=1,nesc)
      ensc(1)=-ensc(1)
      write(ko,'(2i5)') iztar,iatar
      write(ko,'(2i5)') ipsc,ipot
c
c     Factor coming from Dirac equation reduction (relativistic approximation).
c     Using this factor force us to employ relativistic kinematics as well.
      gamma=1.d0
      if(irel.eq.2) then
c       Target system mass in MeV
        EMtar=tarmas*amu0c2
c       Total system mass in MeV
        EMtot=(tarmas+projmas)*amu0c2
c       Total kinetic energy in cm
        Tcm=sqrt(2*EMtar*el + EMtot**2) - EMtot
c       Relativistic correction to the potential (non relativistic target!!).
        gamma=1.d0+Tcm/(Tcm+2*projmas*amu0c2)
      endif
      do 200 i=1,6
      write(ko,'(4f9.5,i5)') rsc(i),resc(i),asc(i),aesc(i),npzen(i)
      do 160 j=1,iabs(npzen(i))
      if(i.le.4) then
c
c     gamma should not be used for SCAT2000 (it must be equal to 1 ).
c     It is left to be used by another OM code
c
      write(ko,5)epotsc(i,j),(gamma*potsc(i,j,k),k=1,npcofx)
      else
      write(ko,5)epotsc(i,j),(      potsc(i,j,k),k=1,npcofx)
      endif
  5   format(f9.3,1p,4(1x,e12.5),/,9x,3(1x,e12.5))
 160  continue
 200  continue
      write(ko,'(3f10.5)') rcoulsc,efermi,betasc
      write(ko,'(i5)') isuit
      return
      end subroutine scatip2

c*************************
c     ECIS96-ECIS03 input
c*************************
      subroutine ecisip
c
c     ROUTINE TO PREPARE AND WRITE ECIS INPUT
      include "om-retrieve.cmb"
      include "om-retrieve-cts.cmb"
c
  42  astep=1.
C 42  astep=10.
      ncol=1
      iterm=20
      npp=1
      rmatch=0.
      convg=1.0e-8
      ecis1=becis1
      ecis2=becis2

C     RCN 08/2004  Output of the cross sections requested
      ecis2(9:9)='T'

      ecis2(15:15)='T'
      ecis2(20:20)='T'
      if(imodel.eq.0)ecis2(40:40)='T'
      if (izaproj.ne.1) ecis2(10:10)='T'
      if (irel.ge.1) ecis1(8:8)='T'
      if(imodel.ne.1)go to 20
      call couch11
      ncol=1
      if(ntar.gt.0)ncol=ncoll(ntar)
      npp=ncol
      if(koptom.eq.1)npp=1
  20  if(kecis.eq.1) go to 50
      ecis2(13:13)='T'
      ecis2(14:14)='F'
      ecis2(15:15)='F'
      astep=180.
  50  do 90 n=1,ne
      el=en(n)
      if(kecis.eq.2)el=en(n)/xkine(el)
c
c     Check for gaussian ff for WD
      call gaussff(igauflg)
      if(igauflg.eq.1)go to 90
c
      if(imodel.ne.1.or.ntar.eq.0)go to 54
      if(spin(1,ntar).ge.2.5.and.el.lt.0.1)convg=1.0e-10
  54  ldwmax=2.4*1.25*iatar**(1./3.)*0.22*sqrt(projmas*en(n))
      njmax=max(2*ldwmax,20)
      jchk=1
      if(kecis.ne.1)jchk=0
      call optmod
      call optname
      write(ko,'(f7.3," MeV ",a8," on ",i3,a2," - Optical model: ",
     * 14a1," REF#=",i5)') el,parname(ipsc),iatar,nuc(iztar),opname,iref
      write(ko,'(a50)') ecis1
      write(ko,'(a50)') ecis2
      write(ko,'(4i5)') ncol,njmax,iterm,npp
      write(ko,'(10x,f10.5,10x,1p,3(2x,e8.1))') rmatch,convg,convg,convg
      if (ecis2(15:15).eq.'T') write(ko,'( )')
      if(imodel.ne.1)go to 60
      call couch22
      go to 62
  60  tgts=0.
      tgtp='+'
      write(ko,'(f5.2,2i2,a1,5f10.5)') tgts,0,1,tgtp,el,
     +  projspi,projmas,tarmas,ztar*zproj
      write(ko,'( )')
  62  do 70 j=1,npp
      if(npp.eq.1)go to 68
      el=en(n)
      if(kecis.eq.2)el=el/xkine(el)
      xratio = tarmas / (tarmas+projmas)
      if(ntar.gt.0) el=el-ex(j,ntar)/xratio
      if(el.lt.0.) el=0.0001
      jchk=0
      if(j.gt.1)call optmod
  68  continue
      write(ko,'(3f10.5)') v,rv,av
      write(ko,'(3f10.5)') w,rw,aw
      write(ko,'(3f10.5)') vd,rvd,avd
      write(ko,'(3f10.5)') wd,rwd,awd
      write(ko,'(3f10.5)') vso,rvso,avso
      write(ko,'(3f10.5)') wso,rwso,awso
C     COULOMB POTENTIAL       
C      1-10   VAL(20)  REDUCED COULOMB RADIUS IN FERMIS.                ECIS-724
C     11-20   VAL(21)  DIFFUSENESS OF A WOODS-SAXON CHARGE DISTRIBUTION.ECIS-725
      write(ko,'(3f10.5)') rc,ac,0.
      write(ko,'(3f10.5)') 0.,0.,0.
      if(j.eq.npp)write(ko,'(3f10.5)') 0.,astep,180.
  70  continue
  90  continue
c     close (unit=ko)
      return
      end subroutine ecisip

      subroutine ecisdwip
c
c     Routine to prepare and write ECIS inputs for vibrational and DWBA cases
c
      include "om-retrieve.cmb"
      include "om-retrieve-cts.cmb"
      character*4 lmodel
      character*1 tgtpx

  14  ncol=2
      iterm=20
      npp=2
C     astep=10.
      astep=1.
      if(koptom.eq.1)npp=1
      rmatch=0.
      convg=1.0e-8
      ecis1=becis1
      ecis2=becis2
C     RCN 08/2004  Output of the cross sections requested
      ecis2(9:9)='T'
      ecis1(12:12)='T'
      if (izaproj.ne.1) ecis2(10:10)='T'
      if (irel.ge.1) ecis1(8:8)='T'
      if(modtyp.eq.2) call vib1
      kexit=0
      if(modtyp.eq.3) iterm=1
      if(modtyp.eq.3) call dwba1(kexit)
      if(kexit.eq.1)return
      lmodel='VIBR'
      if(modtyp.eq.3)lmodel='DWBA'
c
      ecis2(14:14)='F'
      ecis2(15:15)='F'
      ecis2(30:30)='T'
      tgts=sdis(1)
      tgtp='+'
      if(tgts.eq.0.0.and.ipdis(1).ge.0)go to 60
      write(6,'("Nucleus izatar=",i6,": VIB and DWBA options only",
     +  " set up for even-even targets.")') izatar
      stop 16
  60  do 90 i=1,ne
        eninc=en(i)
        ecm=xkine(eninc)*eninc
c
c       Check for gaussian ff for WD
        el=eninc
        call gaussff(igauflg)
        if(igauflg.eq.1)go to 90
c
c       if(tgts.ge.2.5.and.eninc.lt.0.1) convg=1.0e-10
        ldwmax=2.4*1.25*iatar**(1./3.)*0.22*sqrt(projmas*eninc)
        njmax=max(2*ldwmax,20)
        do 88 j=2,ndis
          if (ecm.le.edis(j)+0.1) go to 88
          tgtpx='+'
          if(ipdis(j).lt.0)tgtpx='-'
        write(ko,'(f6.2," MeV ",a8," on ",$)') eninc,parname(ipsc)
        write(ko,'(i3,a2,": ",a4," for ",$)')iatar,nuc(iztar),lmodel
        write(ko,'(f4.1,a1," level at ",f5.3," MeV")')
     +   sdis(j),tgtpx,edis(j)
        write(ko,'(a50)') ecis1
        write(ko,'(a50)') ecis2
        write(ko,'(4i5)') ncol,njmax,iterm,npp
        write(ko,'(10x,f10.5,10x,1p,3(2x,e8.1))') rmatch,convg,convg,
     + convg
        if (ecis2(15:15).eq.'T') write(ko,'( )')
        write(ko,'(f5.2,2i2,a1,5f10.5)') tgts,0,1,tgtp,eninc,
     +    projspi,projmas,tarmas,ztar*zproj
        write(ko,'(2i5)') 0,0
          jj=2
          if(koptom.eq.1)jj=1
c         if(typdis(j).eq.'D')bdis(j)=bdis(j)/(rv*atar**(1./3.))
          write(ko,'(f5.2,2i2,a1,5f10.5)') sdis(j),0,jj,tgtpx,edis(j),
     +      projspi,projmas,tarmas,ztar*zproj
          write(ko,'(2i5)') 1,1
          xjdis=abs(sdis(1)-sdis(j)) + 0.001
          jdis=int(xjdis)
          if(modtyp.eq.3)jdis=langmo(j)
          write(ko,'(i5,5x,f10.5)') jdis,bdis(j)
        do 80 k=1,npp
          if (k.gt.1)go to 74
          el=eninc
          jchk=0
          go to 76
   74     xratio = tarmas / (tarmas+projmas)
          el=eninc-edis(j)/xratio
          if(el.lt.0.) el=0.0001
   76     call optmod
          write(ko,'(3f10.5)') v,rv,av
          write(ko,'(3f10.5)') w,rw,aw
          write(ko,'(3f10.5)') vd,rvd,avd
          write(ko,'(3f10.5)') wd,rwd,awd
          write(ko,'(3f10.5)') vso,rvso,avso
          write(ko,'(3f10.5)') wso,rwso,awso
C         COULOMB POTENTIAL       
C      1-10   VAL(20)  REDUCED COULOMB RADIUS IN FERMIS.                ECIS-724
C     11-20   VAL(21)  DIFFUSENESS OF A WOODS-SAXON CHARGE DISTRIBUTION.ECIS-725
          write(ko,'(3f10.5)') rc,ac,0.
          write(ko,'(3f10.5)') 0.,0.,0.
   80   continue
        write(ko,'(3f10.5)') 0.,astep,180.
   88 continue
   90 continue
      return
      end subroutine ecisdwip

      subroutine vib1
c
c     Routine to setup for vibrational model calculation
c
      include "om-retrieve.cmb"
      include "om-retrieve-cts.cmb"
c
      ecis2(20:20)='T'
c
      ntar=0
      do 44 nis=1,nisotop
      if(iz(nis).ne.iztar.or.ia(nis).ne.iatar)go to 44
      ntar=nis
      go to 48
  44  continue
      write(6,'("Nucleus izatar=",i6," does not have vibrational",
     +  " data.")') izatar
      stop 6
c
  48  ndis=nvib(ntar)
      do 54 k=1,ndis
      edis(k)=exv(k,ntar)
      sdis(k)=spinv(k,ntar)
      ipdis(k)=iparv(k,ntar)
  54  bdis(k)=defv(k,ntar)
      return
      end subroutine vib1

      subroutine dwba1(kexit)
c
c     Routine to setup for vibrational model calculation
c
      include "om-retrieve.cmb"
      include "om-retrieve-cts.cmb"

      character*2 sym,idummy
      character*10 authz,audef1,audef2
c
      if(kdef.eq.1) audef1='JENDL-3.2 '
      if(kdef.eq.1) audef2='JENDL-3.2 '
      if(kdef.eq.2) audef1='ENSDF(Q)  '
      if(kdef.eq.2) audef2='ENSDF(Q)  '
      if(kdef.eq.3) audef1='ENSDF(BE2)'
      if(kdef.eq.3) audef2='ENSDF(BE3)'
      if(kdef.eq.4) audef1='Raman2    '
c     if(kdef.eq.4) audef1='Raman     '
      if(kdef.eq.4) audef2='Kibedi    '
c     if(kdef.eq.4) audef2='Spear     '
      if(kdef.eq.5) audef1='User      '
      if(kdef.eq.5) audef2='User      '
c
      ecis2(20:20)='T'
      if(modtyp.eq.3)ecis2(42:42)='T'
c
      open(unit=16,file='deform.dat')
c
      klev=1
      edis(klev)=0.
      sdis(klev)=0.
      ipdis(klev)=0
      langmo(klev)=0
      bdis(klev)=0.
c
c     Position deform.dat file
      do 25 i=1,4
  25  read(16,'(a2)') idummy
c
  30  read(16,'(2i4,1x,a2,1x,f10.6,1x,f4.1,i3,i2,1x,f10.6,2x,a10)',
     + end=900) kz,ka,sym,edisz,sdisz,ipdisz,langmoz,bdisz,authz
      kza=1000*kz+ka
      if(kza.gt.izatar)go to 900
      if(kza.ne.izatar)go to 30
      if(kdef.eq.5)go to 40
      if(authz.ne.audef1.and.authz.ne.audef2) go to 30
c
  40  klev=klev+1
      edis(klev)=edisz
      sdis(klev)=sdisz
      ipdis(klev)=ipdisz
      langmo(klev)=langmoz
      bdis(klev)=bdisz
      go to 30
c
 900  ntar=0
      ndis=klev
      if(ndis.eq.1) go to 990
      ndis=min0(ndisx,ndis)
      go to 999
c
 990  write(6,'("DWBA deformations not found for izatar=",i7)') izatar
      kexit=1
 999  close (unit=17)
      return
      end subroutine dwba1

      subroutine gaussff(igauflg)
c
c     Routine to see if Gaussian form factor required for WD(el)
c
      include "om-retrieve.cmb"
c
      igauflg=0
      jab=iabs(jrange(4))
      if(jrange(4).lt.1)go to 40
c
      jp=1
      do 30 j=1,jab
      if(el.gt.epot(4,j)) jp=j+1
  30  continue
      j=min0(jp,jab)
      if(rco(4,j,1).lt.0.)igauflg=1
c
  40  return
      end subroutine gaussff

      subroutine optmod
c
c     Routine to generate input for ECIS96 from RIPL-2 library
c
      include "om-retrieve.cmb"
      include "om-retrieve-cts.cmb"

      dimension rlib(6),alib(6),vlib(6),b(6,ndim1,15)
c
c     To use dispersive optical model package
c
c             functions
      real*8 DOM_int_Ws,DOM_int_Wv,DOM_int_T1,DOM_int_T2,Vhf
      real*8 DOM_int,WVf,WDf,Delta_WD,Delta_WV
c             variables
      real*8 As,Bs,Cs,AAv,Bv,AAvso,Bvso,EEE,Ep,Ea,Ef
      real*8 alpha_PB,beta_PB,gamma_PB,Vnonl,VVcoul,VScoul,VDcoul
      real*8 DWS,DWV,DWVnonl,DWVso,DerDWV,DerDWS,dtmp
      real*8 WDE,WVE
	real*8 rvolint,wvolint,rsurint,wsurint,DWVintNONL,DWVint
      real*8 vcoldisp, vcolint,scolint 

      integer n,iq,nns

      common /energy/EEE,Ef,Ep
      common /Wenerg/WDE,WVE
      common /pdatas/As,Bs,Cs,nns,iq
      common /pdatav/AAv,Bv,n

c     Nonlocality constants alpha fixed according to Mahaux 1991
      real*8 AlphaV
      data AlphaV/1.65d0/
c     real*8 AlphaS
c     data AlphaS/0.235d0/
      external Delta_WV,WVf,Delta_WD,WDf

cc
c     Generate optical model parameters for ECIS
c
      pi=acos(-1.d0)
      VDcoul = 0.d0
      rc=0.d0
      if(jcoul.lt.1)go to 194
      jc=1
      do 190 j=1,jcoul
      if(el.gt.ecoul(j)) jc=j+1
 190  continue
      jc=min0(jc,jcoul)
      rc=rcoul0(jc)*atar**(-1./3.) + rcoul(jc) +
     +   rcoul1(jc)*atar**(-2./3.) + rcoul2(jc)*atar**(-5./3.)
	ac=acoul(jc)
 194  encoul2=0.
      if(rc.gt.0.) encoul2=1.73*ztar/(rc*atar**(1./3.))
c
      do 300 i=1,6

      vc = 0.0
      VVcoul = 0.d0
      VScoul = 0.d0
      DerDWV = 0.d0
      DerDWS = 0.d0
      DWVnonl = 0.d0

      rlib(i)=0.
      alib(i)=0.
      vlib(i)=0.
      jab=iabs(jrange(i))
      if(jrange(i).lt.1)go to 300
      jp=1
      do 204 j=1,jab
      if(el.gt.epot(i,j)) jp=j+1
 204  continue
      j=min0(jp,jab)
      Ef=dble(int(100000*efermi))/100000
      if(pot(i,j,18).ne.0.) Ef=pot(i,j,18) + pot(i,j,19)*atar
      elf = el - Ef
c
c     Calculate radius and diffuseness parameters
      if(rco(i,j,13).eq.0.) then
        rlib(i)=abs(rco(i,j,1)) + rco(i,j,3)*eta
     *       + rco(i,j,4)/atar + rco(i,j,5)/sqrt(atar)
     *       + rco(i,j,6)*atar**(2./3.) + rco(i,j,7)*atar
     *       + rco(i,j,8)*atar**2  + rco(i,j,9)*atar**3
     *       + rco(i,j,10)*atar**(1./3.)
     *       + rco(i,j,11)*atar**(-1./3.)
C--------------------------------------------------------------------
C     RCN, 08/2004, to handle new extension to the OMP RIPL-2 format
     *       + rco(i,j,2)*el
     *       + rco(i,j,12)*el*el
      else
C     RCN, 09/2004, to handle new extension to the OMP RIPL-2 format
        nn = int(rco(i,j,7))
        rlib(i)= ( abs(rco(i,j,1)) + rco(i,j,2)*atar ) *
     *           ( 1.d0 - ( rco(i,j,3) + rco(i,j,4)*atar ) * elf**nn/
     *           ( elf**nn + ( rco(i,j,5) + rco(i,j,6)*atar )**nn ) )
      endif

      alib(i)=abs(aco(i,j,1)) + aco(i,j,2)*el + aco(i,j,3)*eta
     *        + aco(i,j,4)/atar + aco(i,j,5)/sqrt(atar)
     *        + aco(i,j,6)*atar**(2./3.) + aco(i,j,7)*atar
     *        + aco(i,j,8)*atar**2 + aco(i,j,9)*atar**3
     *        + aco(i,j,10)*atar**(1./3.) + aco(i,j,11)*atar**(-1./3.)
c
c
      if (pot(i,j,24).eq.0.) go to 210
c
c     Special Koning-type potential formulas
c
      if (pot(i,j,24).eq.1.) then
c
c       Koning-type formulas
c
        elf = el - Ef

        if(i.eq.1) call bcoget(b,j)
        if(i.eq.1 .and. b(1,j,5).ne.0.) then
	    vc = b(1,j,1)*encoul2*( b(1,j,2) - 2.*b(1,j,3)*elf +
     +    3.*b(1,j,4)*elf**2 + b(i,j,14)*b(i,j,13)*exp(-b(i,j,14)*elf) )
	    VDcoul = b(i,j,5)*vc
	  endif 

        nn = int(pot(i,j,13))
c       Retrieving average energy of the particle states Ep
        Ep=Ef
        if( (i.eq.2) .or. (i.eq.4) )
     >     Ep=dble(int(100000*pot(i,j,20)))/100000
        if(Ep.eq.0.) Ep=Ef
        elf = el - Ep

        iq=1
        if(i.eq.4 .and. b(4,j,12).gt.0.) iq=nint(b(4,j,12))

        vlib(i)=
     +    b(i,j,1)*( b(i,j,15) - b(i,j,2)*elf + b(i,j,3)*elf**2 -
     +    b(i,j,4)*elf**3 + b(i,j,13)*exp(-b(i,j,14)*elf) ) +
     +    b(i,j,5)*vc + b(i,j,6)*(elf**nn/(elf**nn + b(i,j,7)**nn)) +
     +    b(i,j,8)*exp(-b(i,j,9)*elf**iq)*(elf**nn/
     +    (elf**nn + b(i,j,10)**nn)) +
     +    b(i,j,11)*exp(-b(i,j,12)*elf)

      endif

      if (pot(i,j,24).eq.2.) then
c
c       Morillon-Romain formulas
c
        elf = el - Ef
        if(i.eq.1) call bcoget(b,j)

c       Vhf(E) calculated from nonlocal approximation
c          as suggested by Perey and Buck
c
        alpha_PB = dble(int(100000*b(i,j,1)))/100000
        beta_PB  = dble(int(100000*b(i,j,2)))/100000
        gamma_PB = dble(int(100000*b(i,j,3)))/100000

        EEE = dble(int(1000000*el))/1000000

        iq=1
        if(i.eq.4 .and. b(4,j,12).gt.0.) iq=nint(b(4,j,12))

        Vnonl = 0.
        if(i.eq.1) then
	    Vnonl = -Vhf(EEE,alpha_PB,beta_PB,gamma_PB)
C         Numerical derivative of the Vhf
          if(b(1,j,5).ne.0.) then
	      Vnonlm = -Vhf(EEE-0.05,alpha_PB,beta_PB,gamma_PB)
	      Vnonlp = -Vhf(EEE+0.05,alpha_PB,beta_PB,gamma_PB)
C           Coulomb correction for Hartree-Fock potential
            vc = encoul2*(Vnonlm-Vnonlp)*10.d0
	      VDcoul = b(i,j,5)*vc
          endif
        endif
        vlib(i)=
     +    Vnonl + b(i,j,5)*vc +
     +    b(i,j,6)*(elf**nn/(elf**nn + b(i,j,7)**nn)) +
     +    b(i,j,8)*exp(-b(i,j,9)*elf**iq)*(elf**nn/
     +    (elf**nn + b(i,j,10)**nn)) +
     +    b(i,j,11)*exp(-b(i,j,12)*elf)
      endif
c
c     Nonlocality consideration
c
c     Retrieving energy above which nonlocality in the volume absorptive
c               potential is considered (Ea)
c
      Ea=dble(int(100000*pot(i,j,21)))/100000
      if(Ea.eq.0.) Ea=1000.1d0
      if(i.eq.2 .and. Ea.lt.1000. .and. el.gt.(Ef+Ea) )
     +         vlib(i)= vlib(i) +
     +       AlphaV*(sqrt(el)+(Ef+Ea)**1.5/(2.*el)-1.5*sqrt(Ef+Ea))
c
      go to 300
c
 210  if (pot(i,j,23).eq.0.) go to 220
c
c     Special Varner-type potential formulas
      vlib(i)= (pot(i,j,1) + pot(i,j,2)*eta)/
     *    (1.+ exp((pot(i,j,3) - el + pot(i,j,4)*encoul2)/pot(i,j,5)))
      if(pot(i,j,6).eq.0.)go to 300
      vlib(i) = vlib(i)
     *    + pot(i,j,6)*exp((pot(i,j,7)*el - pot(i,j,8))/pot(i,j,6))
      go to 300
c
 220  if (pot(i,j,22).eq.0.) go to 230
c
c     Special Smith-type potential formulas
      vlib(i)=pot(i,j,1) + pot(i,j,2)*eta
     *    + pot(i,j,6)*exp(pot(i,j,7)*el + pot(i,j,8)*el*el)
     *    + pot(i,j,9)*el*exp(pot(i,j,10)*el**pot(i,j,11))
      if(pot(i,j,5).ne.0.)vlib(i)=vlib(i)
     *    + pot(i,j,3)*cos(2.*pi*(atar - pot(i,j,4))/pot(i,j,5))
      go to 300
c
c     Standard potential formulas
 230  vlib(i)=pot(i,j,1) + pot(i,j,7)*eta + pot(i,j,8)*encoul
     *    + pot(i,j,9)*atar + pot(i,j,10)*atar**(1./3.)
     *    + pot(i,j,11)*atar**(-2./3.) + pot(i,j,12)*encoul2
     *    + (pot(i,j,2) + pot(i,j,13)*eta + pot(i,j,14)*atar)*el
     *    + pot(i,j,3)*el*el + pot(i,j,4)*el*el*el + pot(i,j,6)*sqrt(el)
      if(el.gt.0.) vlib(i) = vlib(i) + pot(i,j,17)*encoul/el**2
     *    + (pot(i,j,5) + pot(i,j,15)*eta + pot(i,j,16)*el)*alog(el)
 300  continue
c
c     To calculate dispersion relation contribution
c
 152  if(abs(idr).ge.2) then
c
c       Exact calculation of the dispersive contribution
c
        EEE = dble(int(1000000*el))/1000000

        i=2
c       Only one energy range
        j=1
c       Real volume contribution from Dispersive relation
        DWV=0.d0
c
        if(pot(2,1,24).ne.0) then
          AAv=dble(int(100000*b(i,j,6)))/100000
          Bv=dble(int(100000*b(i,j,7)))/100000
          n = nint( pot(i,j,13) )
          if(n.eq.0 .or. mod(n,2).eq.1)
     +      stop 'Zero or odd exponent in Wv(E) for dispersive OMP'

c         Retrieving average energy of the particle states Ep
          Ep=dble(int(100000*pot(i,j,20)))/100000
          if(Ep.eq.0.) Ep=Ef

c         Analytical DOM integral
          DWV=DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE,n,DerDWV)

c         Coulomb correction for real volume potential 
 	    DerDWV = -b(1,1,5)*encoul2*DerDWV
C         numerical DOM derivative (not needed for a time being) 
C         DWVp = DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE+0.1d0,n,dtmp)
C         DWVm = DOM_INT_Wv(Ef,Ep,AAv,Bv,EEE-0.1d0,n,dtmp)
C         DerDWV = -b(1,1,5)*encoul2*(DWVp-DWVm)*5.d0

c         if(idr.le.-2) then
c           numerical DOM integral (not needed for a time being)
c           WVE=WVf(AAv,Bv,Ep,Ef,EEE,n)
c           DWV=2*DOM_int(Delta_WV,WVf,Ef,Ef+5.*Bv,150000.d0,EEE,0.d0)
c         endif

c
c         Nonlocality correction to the DOM integral
c           (only used if Ea is non-zero)
          Ea=dble(int(100000*pot(i,j,21)))/100000
          if(Ea.eq.0.) Ea=1000.1d0
          T12der = 0.d0  
          if(Ea.lt.1000.) THEN
             Dwplus = AlphaV*DOM_INT_T2(Ef,Ea,EEE)
             dtmp1 = Wvf(AAv,Bv,Ep,Ef,Ef+Ea,n)
             Dwmin = dtmp1*DOM_INT_T1(Ef,Ea,EEE)
             DWV = DWV + Dwplus + Dwmin
	       DWVnonl = Dwplus + Dwmin
c            Coulomb correction for nonlocal dispersive contribution  
c                to real volume potential 
 	       if(b(1,1,5).ne.0.d0) then
               if(eee.ne.0.05d0) then
                 T2p = DOM_INT_T2(Ef,Ea,EEE+0.05d0)  
                 T2m = DOM_INT_T2(Ef,Ea,EEE-0.05d0)  
		 	   T2der = AlphaV*(T2p-T2m)*10.d0 
	           T1p = DOM_INT_T1(Ef,Ea,EEE+0.05d0)
	           T1m = DOM_INT_T1(Ef,Ea,EEE-0.05d0)
	           T1der = dtmp1*(T1p-T1m)*10.d0
	           T12der =  -b(1,1,5)*encoul2* ( T1der + T2der )
			 else
                 T2p = DOM_INT_T2(Ef,Ea,EEE+0.1d0)  
                 T2m = DOM_INT_T2(Ef,Ea,EEE-0.1d0)
			   T2der = AlphaV*(T2p-T2m)*5.d0 
	           T1p = DOM_INT_T1(Ef,Ea,EEE+0.1d0)
	           T1m = DOM_INT_T1(Ef,Ea,EEE-0.1d0)
	           T1der = dtmp1*(T1p-T1m)*5.d0
	           T12der =  -b(1,1,5)*encoul2* ( T1der + T2der )
               endif 
	       endif

          endif
	    VVcoul = DerDWV + T12der 
        endif

        i=4
c       Only one energy range
        j=1
c       Real surface contribution from Dispersive relation
        DWS=0.d0

        if(pot(4,1,24).ne.0) then
          As=dble(int(100000*b(i,j,8)))/100000
          Bs=dble(int(100000*b(i,j,10)))/100000
          Cs=dble(int(100000*b(i,j,9)))/100000
          n = nint( pot(i,j,13) )
          if(n.eq.0 .or. mod(n,2).eq.1)
     +      stop 'Zero or odd exponent in Wd(E) for dispersive OMP'

          iq=1
          if(b(4,j,12).gt.0.) iq=nint(b(4,j,12))

c         Retrieving average energy of the particle states Ep
          Ep=dble(int(100000*pot(i,j,20)))/100000
          if(Ep.eq.0.) Ep=Ef
C         Ea=dble(int(100000*pot(i,j,21)))/100000
C         if(Ea.eq.0.) Ea=1000.1d0

          if(idr.ge.2) then
c           analytical DOM integral
            DWS = DOM_INT_Ws(Ef,Ep,As,Bs,Cs,EEE,n,DerDWS)
c           Coulomb correction for real surface potential 
		  VScoul = -b(1,1,5)*encoul2*DerDWS
          endif

          if(idr.le.-2) then
c           numerical DOM integral
            nns=n
            WDE=WDf(As,Bs,Cs,Ep,EEE,n,iq)
            DWS = 2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE,WDE)
c           Coulomb correction for real surface potential 
	      if(b(1,1,5).ne.0.d0) then
              WDE=WDf(As,Bs,Cs,Ep,EEE+0.1d0,n,iq)
              DWSp = 
     >	 	 2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE+0.1d0,WDE)
              WDE=WDf(As,Bs,Cs,Ep,EEE-0.1d0,n,iq)
              DWSm =
     >		 2*DOM_int(Delta_WD,WDf,Ef,Ef+30.d0,2000.d0,EEE-0.1d0,WDE)
c             Numerical derivative
	        VScoul = -b(1,1,5)*encoul2*(DWSp-DWSm)*5.d0
	      endif
          endif

        endif

        i=6
c       Only one energy range
        j=1
c       Real spin orbit contribution from Dispersive relation
        DWVso=0.d0
c
        if(pot(6,1,24).ne.0 .and. abs(idr).eq.3) then

          AAvso=dble(int(100000*b(i,j,6)))/100000
          Bvso=dble(int(100000*b(i,j,7)))/100000
          n = nint( pot(i,j,13) )
          if(n.eq.0 .or. mod(n,2).eq.1)
     +      stop 'Zero or odd exponent in Wso(E) for dispersive OMP'

c         analytical DOM integral
          DWVso=DOM_INT_Wv(Ef,Ef,AAvso,Bvso,EEE,n,dtmp)

        endif

c       Adding real volume dispersive contribution to the real potential
c       Geometry parameters are the same as for the volume potential(imag and real).
        Vhfnum = vlib(1)
        vlib(1)= vlib(1) + DWV + VVcoul
c       Including real surface and Coulomb dispersive contribution
c       Geometry parameters are the same as for the imaginary surface potential.
        vlib(3)= DWS + VScoul
        alib(3)= alib(4)
        rlib(3)= rlib(4)
c       Adding real spin orbit dispersive contribution to the real spin orbit potential
c       Geometry parameters are the same as for the imaginary spin orbit potential(imag and real)
        vlib(5) = vlib(5) + DWVso

      endif

c======= Capote 2001 =======
c
c     Factor coming from Dirac equation reduction (relativistic approximation).
c     Using this factor force us to employ relativistic kinematics as well.
      gamma=1.d0
      if(irel.eq.2) then
c       Target system mass in MeV
        EMtar=tarmas*amu0c2
c       Total system mass in MeV
        EMtot=(tarmas+projmas)*amu0c2
c       Total kinetic energy in cm
        Tcm=sqrt(2*EMtar*el + EMtot**2) - EMtot
c       Relativistic correction to the potential (non relativistic target!!).
        gamma=1.d0+Tcm/(Tcm+2*projmas*amu0c2)
      endif
c
c     gamma should not be used for ECIS (it must be equal to 1 ).
c     It is left to be used by another OM code
c
      v=vlib(1)*gamma
      rv=rlib(1)
      av=alib(1)
      w =vlib(2)*gamma
      if(w.lt.0.0.and.nonegw.eq.1) w=0.
      rw=rlib(2)
      aw=alib(2)
      vd=vlib(3)*gamma
      rvd=rlib(3)
      avd=alib(3)
      wd=vlib(4)*gamma
      if(wd.lt.0.0.and.nonegw.eq.1) wd=0.
      rwd=rlib(4)
      awd=alib(4)
      vso=vlib(5)
      rvso=rlib(5)
      avso=alib(5)
      wso=vlib(6)
      rwso=rlib(6)
      awso=alib(6)

      if(jchk.eq.1) then
C
C       SPHERICAL VOLUME INTEGRALS 
C
        call VOLIN(1,dble(v) ,dble(rv *atar**(1./3.)),dble(av),rvolint)
        rvolint = rvolint  / atar
        
	  call VOLIN(1,dble(DWV),dble(rv *atar**(1./3.)),dble(av),DWVint)
        DWVint = DWVint  / atar

	  call VOLIN(1,dble(DWVnonl),dble(rv *atar**(1./3.)),
     >               dble(av),DWVintNONL)
        DWVintNONL = DWVintNONL  / atar
        
	  call VOLIN(1,dble(w) ,dble(rw *atar**(1./3.)),dble(aw) ,wvolint)
        wvolint = wvolint  / atar

        call VOLIN(2,dble(vd),dble(rvd*atar**(1./3.)),dble(avd),rsurint)
	  rsurint = rsurint / atar

        call VOLIN(2,dble(wd),dble(rwd*atar**(1./3.)),dble(awd),wsurint)
        wsurint = wsurint / atar

        vcolint =	0.d0
        scolint = 0.d0
	  vcoldisp = 0.d0
	  
	  if(ztar.gt.0) then
	  vc = (VDcoul + VVcoul)*gamma
        call VOLIN(1,dble(vc),dble(rv *atar**(1./3.)),dble(av),vcolint)
	  vcolint = vcolint / atar

	  vc = VVcoul*gamma
        call VOLIN(1,dble(vc),dble(rv *atar**(1./3.)),dble(av),vcoldisp)
        vcoldisp = vcoldisp / atar
	  
	  vc = VScoul*gamma 
        call VOLIN(2,dble(vc),dble(rvd*atar**(1./3.)),dble(avd),scolint)
        scolint = scolint / atar
	  endif
C
C       Printing potential values
C
        write(35,'(f7.3,6(f7.3,f6.3,f6.3),13(1x,f6.3,1x))')
     +    el,v,rv,av,w,rw,aw,vd,rvd,avd,wd,rwd,awd,
     +    vso,rvso,avso,wso,rwso,awso
        write(36,'(f7.3,15(1x,f8.3,1x))')
     +    el,DWV,DWV-Dwplus-Dwmin,Dwplus,Dwmin,Dwvso,
     +    v,Vhfnum,v-Vhfnum,VDCoul,VVcoul,T12der,VScoul,DWS,vlib(3)
        write(37,'(f7.3,2x,9(f8.3,2x))')
     +    el,rvolint+rsurint,rvolint,rsurint,
     +       DWVint,DWVintNONL,DWVint-DWVintNONL,
     +       wvolint+wsurint,wvolint,wsurint
        write(38,'(f7.3,2x,3(f8.3,2x),2x,4(f8.3,3x))')
     +    el,rvolint+rsurint,rvolint,rsurint,
     +       vcolint+scolint,vcolint,scolint,vcoldisp
	endif

      return
      end subroutine optmod

      subroutine bcoget(b,j)
c
c     Routine to compute b coefficients for Koning global potential
c
c     Modified by R.Capote to extend RIPL format
c
      include "om-retrieve.cmb"
      dimension b(6,ndim1,15)
c
       call setr(0.,b,90*ndim1)

c      Original Koning dependence
       b(1,j,1)  =  pot(1,j,1) + pot(1,j,2)*atar + pot(1,j,8)*eta

       if((pot(1,j,20).ne.0.) .and.
     +    (pot(1,j,14) + pot(1,j,15)*atar + pot(1,j,16)).ne.0. ) then
c        Soukhovitski dependence
         b(1,j,1)  =  pot(1,j,1) + pot(1,j,2)*atar + pot(1,j,8)*eta +
     +                pot(1,j,20)*eta/
     +               (pot(1,j,14) + pot(1,j,15)*atar + pot(1,j,16))
       endif
       b(1,j,2)  =  pot(1,j,3) + pot(1,j,4)*atar
       b(1,j,3)  =  pot(1,j,5) + pot(1,j,6)*atar
       b(1,j,4)  =  pot(1,j,7)
       b(1,j,5)  =  pot(1,j,9)
       b(1,j,11) =  pot(1,j,10) + pot(1,j,11)*atar
       b(1,j,12) =  pot(1,j,12)

c      b coefficients from 13 to 15 added for Soukhovitski potential
c      V^DISP_R
       b(1,j,13)  =  pot(1,j,16)
c      Lambda_R
       b(1,j,14)  =  pot(1,j,17)
c      V^0_R + V^A_R*(A-232)
       b(1,j,15)  =  pot(1,j,14) + pot(1,j,15)*atar
c      To preserve compatibility with RIPL-2 Koning database
c      b(i,j,15) must be equal to 1. !!! for Koning OMP
       if((pot(1,j,14) + pot(1,j,15)*atar + pot(1,j,16)).eq.0.)
     >          b(1,j,15) = 1.
C      if(abs(b(1,j,15)).lt.1.e-8) b(1,j,15) = 1.

c      Wv( Av )
       b(2,j,6)  =  pot(2,j,1) + pot(2,j,2)*atar
c      Wv( Bv )
       b(2,j,7)  =  pot(2,j,3) + pot(2,j,4)*atar

c      Wd
c
c      added A dependence for As parameter (RCN, 09/2004), i.e.  pot(4,j,7)<>0
c      Wd( As )
c      b(4,j,8)  =  pot(4,j,1) + pot(4,j,8)*eta
       b(4,j,8)  =  pot(4,j,1) + pot(4,j,8)*eta + pot(4,j,7)*atar

c      Wd( Cs )
       if(pot(4,j,3).ne.0.) then
         b(4,j,9)  =  pot(4,j,2) +
     +                pot(4,j,3)/(1.+ exp((atar-pot(4,j,4))/pot(4,j,5)))
       else
         b(4,j,9)  =  pot(4,j,2)
       endif
c      Wd( Bs )
       b(4,j,10) =  pot(4,j,6)
C      Wd( q )
       b(4,j,12) =  pot(4,j,12)

c      Vso
       b(5,j,11) =  pot(5,j,10) + pot(5,j,11)*atar
       b(5,j,12) =  pot(5,j,12)
c      Wso
       b(6,j,6)  =  pot(6,j,1)
       b(6,j,7)  =  pot(6,j,3)

      return
      end subroutine bcoget

      subroutine optname
c
c     Routine to load opname with first author name (or first
c     14 characters).
c
      include "om-retrieve.cmb"
      character*1 opnamex
      dimension opnamex(14)
c
      do 20 i=1,14
  20  opnamex(i)=author(i)
      do 22 i=2,14
      j=i
      if(opnamex(i).eq.' '.and.opnamex(i-1).ne.'.')go to 24
      if(opnamex(i).eq.',')go to 24
  22  continue
  24  do 26 ii=1,j
      i=j-ii+1
      if(opnamex(i).ne.' '.and.opnamex(i).ne.',')go to 28
  26  continue
  28  jj=i
      do 30 i=1,14
  30  opname(i)=' '
      do 32 i=1,jj
  32  opname(i)=opnamex(i)
c
c     Set names for table print
      if(imodel.eq.0)modelz='spherical  '
      if(imodel.eq.1)modelz='rotational '
      if(imodel.eq.2)modelz='vibrational'
      if(imodel.eq.3)modelz='non-axial  '
      if(idr.eq.0) dispz='non-dispersive'
      if(idr.lt.0) dispz='dispersive   (Num.Integr.) '
      if(idr.gt.0) dispz='dispersive   (Anal.Integr.)'
      if(idr.eq.-3)dispz='dispersive+SO(Num.Integr.) '
      if(idr.eq. 3)dispz='dispersive+SO(Anal.Integr.)'
      if(irel.eq.0)relz='non-relativistic'
      if(irel.eq.1)relz='relativ.Schoedr.'
      if(irel.eq.2)relz='relativ + gamma '
      return
      end subroutine optname

      subroutine couch11
c
c     First setup for coupled channels parameters
c
      include "om-retrieve.cmb"
c
      ecis1(1:1)='T'
      ecis1(12:12)='T'
      ecis1(21:21)='T'
      ntar=0
      do 44 nis=1,nisotop
      if(iz(nis).ne.iztar.or.ia(nis).ne.iatar)go to 44
      ntar=nis
      go to 48
  44  continue
      write(6,'("Nucleus izatar=",i6," does not have coupled-channel",
     +  " data.")') izatar
      write(6,'("RIPL-2 reference number: ",i6)')iref
c     stop 6
  48  return
      end subroutine couch11

      subroutine couch22
c
c     Second setup for coupled channels parameters
c
      include "om-retrieve.cmb"
c
      if(ntar.ne.0)go to 54
      write(ko,'("PUT ROTATIONAL BAND STATES IN HERE")')
      write(6,'("Nucleus izatar=",i6," does not have rotational",
     +  " data.")') izatar
      go to 62
  54  do 56 k=1,ncol
      elx=ex(k,ntar)
      if(k.eq.1) elx=el
      tarspi=spin(k,ntar)
      tarpar='+'
      if(ipar(k,ntar).eq.-1)tarpar='-'
      kom=k
      if(koptom.eq.1)kom=1
      write(ko,'(f5.2,2i2,a1,5f10.5)') tarspi,0,kom,tarpar,elx,
     +  projspi,projmas,tarmas,ztar*zproj
  56  continue
      iqm=idef(ntar)
      iqmax=lmax(ntar)
      aspin=bandk(ntar)
      ik=1
      write(ko,'(2i5,f10.5,i5)') iqm,iqmax,aspin,ik
c     iq=iqm/2
      write(ko,'(7f10.5)') (def(ntar,j),j=2,iqm,2)
  62  return
      end subroutine couch22

      subroutine masses
c
c     Routine to retrieve masses and compute separation energies
c     and Fermi energy
c
      include "om-retrieve.cmb"
      include "om-retrieve-cts.cmb"
      character*1 dp
      real*8 ck2
c-------------------------------------------------------------------------------
c  Physical  constants (following ECIS03 definitions and ENDF manual 2001)
c-------------------------------------------------------------------------------
C CONSTANTS COMPUTED FROM THE FUNDAMENTAL CONSTANTS, ATOMIC MASS, HBAR*C
C AND ALPHA, AS GIVEN IN THE EUROPEAN PHYSICAL JOURNAL, PAGE 73, VOLUME 
C 15 (2000) REFERRING FOR THESE VALUES TO THE 1998 CODATA SET WHICH MAY 
C BE FOUND AT http://physics.nist.gov/constants                         
C     amu0c2=931.494013 +/- 0.000037 MeV                                
      amu0c2=931.494013D0                                               
C     hbarc=197.3269601 +/- 0.0000078 (*1.0E-9 eV*cm)                   
      hbarc=197.3269601D0                                               
C
      ck2    = (2.d0*amu0c2)/(hbarc**2)
c
c     Compute lab to cm factor
      zatar=float(izatar)
      call energy2(zatar,tarmas,tarspi,tarpar,tarex)
      zaproj=float(izaproj)
      call energy2(zaproj,projmas,projspi,projpar,projex)
c     Compute Fermi energy
      call energy2(zatar-zaproj,dm,ds,dp,tarexm1)
      call energy2(zatar+zaproj,dm,ds,dp,tarexp1)
      efermi=-0.5*(tarexm1-tarexp1+2.*projex)
      return
      end subroutine masses

      subroutine energy2(za,exactmas,spin,parity,excess)
c
      include "om-retrieve-cts.cmb"
c
c     Routine energy2 looks up values of ground-state mass excess
c     (in MeV), spin, and parity. Missing data produces a flag.
c     Note that negative za omits Duflo for missing data.
c
      common/xener/nmass,inpgrd,lza(9200),ener(9200),spnpar(9200)
      character*4 bcd
      character*1 parity
      dimension bcd(120)
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
      call duflo(nneu,npro,excess)
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
      end subroutine energy2

      subroutine duflo(nn,nz,excess)
      include "om-retrieve-cts.cmb"
c
c                --------------------------------------------
c                |Nuclear mass formula of Duflo-Zuker (1992)|
c                --------------------------------------------
c
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
      end subroutine duflo

      subroutine egrid(ne,eminx,emaxx,en)
c     Routine to generate auto energy grid.
c     Currently setup to produce grid for A. Koning
c
      dimension en(500),enx(500)
c
      enx(1)=0.
      enx(2)=0.001
      enx(3)=0.002
c
c     energy range 0-0.010 MeV, delta E = 0.002 MeV
      do 10 i=4,7
  10  enx(i)=enx(i-1)+0.002
c
c    energy range 0.01-0.10 MeV, delta E = 0.005 MeV
      do 20 i=8,25
  20  enx(i)=enx(i-1)+0.005
c
c    energy range 0.1-1.0 MeV, delta E = 0.05 MeV
      do 30 i=26,43
  30  enx(i)=enx(i-1)+0.05
c
c    energy range 1.0-10.0 MeV, delta E = 0.2 MeV
      do 40 i=44,88
  40  enx(i)=enx(i-1)+0.2
c
c    energy range 10-200 MeV, delta E = 0.5 MeV
      do 50 i=89,468
      nex=i
  50  enx(i)=enx(i-1)+0.5
c
c     Now put energy grid between potential limits.
      ne=1
      do 80 n=2,nex
      if(enx(n).lt.eminx.or.enx(n).gt.emaxx)go to 80
      ne=ne+1
      en(ne)=enx(n)
  80  continue
      return
      end subroutine egrid

      subroutine tableset
c
c     Routine to setup for table print
c
      include "om-retrieve.cmb"
c
      write(35,'(/"Reference Number =",i5," Incident Particle: ",a8)')
     +  iref,parname(ipsc)
      write(35,'("Target: Z=",i3," A=",i4," Number of energies=",i4)')
     +     iztar,iatar,ne
      write(35,'("Type potential: ",a11,2x,a27,2x,a16)')
     +              modelz,dispz,relz
      write(35,'("First Author: ",14a1)') (opname(i),i=1,14)
      write(35,'("Z-Range=",i2,"-",i3,"  A-Range=",i3,"-",i3,"  E-Range="
     +",f8.3,"-",f8.3," MeV         Ef =",f7.3)')
     +  izmin,izmax,iamin,iamax,emin,emax,Efermi
      write(35,'(/,"En(MeV)     V    RV    AV      W    RW    AW ",
     +                    "    VD   RVD   AVD     WD   RWD   AWD ",
     +                    "   VSO   RVSO  AVSO   WSO   RWSO  AWSO")')

	if(idr.GT.0) then
      write(36,'("Reference Number =",i5," Incident Particle: ",a8,
	+ " Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +  iref,parname(ipsc),iztar,iatar,Efermi
      write(36,'("En(MeV)     DWV    DWV(symm)     DW>      DW<     ",
     +             " DWvso      V         Vhf       Vdisp    VDcoul   ",
     +             " VVcoul    T12der    VScoul +   DWS   =   Vs")')
	endif

      write(37,'("Reference Number =",i5," Incident Particle: ",a8,
	+ " Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +  iref,parname(ipsc),iztar,iatar,Efermi
      write(37,'("En(MeV)   RealPot   RVolPot   RSurPot   DVolPot",
     +"   DVolNonL  DVolSym    ImPot    ImVolPot  ImSurPot")')

	if(iztar.GT.0) then
      write(38,'("Reference Number =",i5," Incident Particle: ",a8,
	+ " Target: Z=",i3," A=",i4," Ef =",f7.3)')
     +  iref,parname(ipsc),iztar,iatar,Efermi
      write(38,'("En(MeV)   RealPot   RVolPot    RSurPot ",
     +"  Int(Coul)  VInt(Coul) SInt(Coul) VIntDisp(Coul)")')
	endif

      return
      end subroutine tableset

      subroutine block data
      include "om-retrieve.cmb"
      data ki,idr/1,0/
      data (nuc(i),i=1,103) /'H ','He','Li','Be','B ','C ','N ','O ',
     +  'F ','Ne','Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     +  'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
     +  'As','Se','Br','Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru',
     +  'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ','Xe','Cs','Ba',
     +  'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er',
     +  'Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     +  'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th','Pa','U ',
     +  'Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr'/
      data (parname(i),i=1,7) /'neutron ','proton  ','deuteron',
     +  'triton  ','he-3    ','alpha   ','gamma   '/
      end	subroutine block data

      REAL*8 FUNCTION Vhf(einp,alpha_PB,beta_PB,gamma_PB)
c
c     According to Morillon B, Romain P, PRC70(2004)014601
c
c     Originally coded in c++ by Morillon B. and Romain P.
c
c     Coded in FORTRAN and tested by RCN, August 2004.
c
      real*8 einp,alpha_PB,beta_PB,gamma_PB
      real*8 Vtmp,Etmp,miu_sur_hbar2, coef1, coef2
      integer niter
      include "om-retrieve-cts.cmb"

c     getting amu
      xtmp=xkine(sngl(einp))
      miu_sur_hbar2 = amu / hbarc**2
      coef1 = -0.5d0 * beta_PB**2 * miu_sur_hbar2
      coef2 =  4.0d0 * (gamma_PB * miu_sur_hbar2)**2

      niter = 0
      Vhf = -45.d0
10    niter = niter + 1
      Vtmp = Vhf
      Etmp = einp - Vtmp
      Vhf = alpha_PB * dexp(coef1 * Etmp + coef2 * Etmp**2)
      if( abs(Vhf - Vtmp) .GT. 0.0001 .AND.  niter.LT.10000) goto 10
      return
      end FUNCTION Vhf

      REAL FUNCTION xkine(ei)
c***********************************************************************
c     From lab to CM (the input quantity is el = Elab)
c***********************************************************************
c     RCN 08/2004, xkine calculated by relativistic kinematics when needed

      include "om-retrieve.cmb"
      include "om-retrieve-cts.cmb"

      mtot = (tarmas+projmas)
      xratio = tarmas / mtot

      if(irel .eq. 0) then
c
c-----------------------------------------------------------------------
c  Classical    kinematics (energy independent amu and xkine)
c-----------------------------------------------------------------------
c
          amu = projmas*tarmas / mtot * amu0c2
          xkine = tarmas / mtot
c         e1  = el*xkine
c         w2  = ck2*amu
c         ak2 = w2*e1
      else
c
c-----------------------------------------------------------------------
c  Relativistic kinematics
c-----------------------------------------------------------------------
c
c         e1  = amu0c2*mtot*
c    * (        DSQRT(1.d0 +
c    *       2.d0*el/(amu0c2*tarmas*((1.d0+projmas/tarmas)**2))) - 1.d0)
          p2  = (ei*(ei + 2.d0*amu0c2*projmas)) /
     *          ((1.d0+projmas/tarmas)**2 + 2.d0*ei/(amu0c2*tarmas))
c         ak2 = p2 / (hbarc*hbarc)
          etoti = DSQRT((amu0c2*projmas)**2 + p2)
          etott = DSQRT((amu0c2*tarmas)**2  + p2)
          amu   = etoti*etott / (etoti + etott)
c         amu   = amu / amu0c2
          xkine = etott / (etoti + etott)
      endif
      return
      end FUNCTION xkine
C
C==========================================================================
C     AUTHOR: Dr. Roberto Capote Noy
c
C     e-mail: rcapotenoy@yahoo.com; r.capotenoy@iaea.org
C
C     DISPERSIVE OPTICAL MODEL POTENTIAL PACKAGE
C
c     Analytical dispersive integrals are included
c     see Quesada JM et al, Computer Physics Communications 153(2003) 97
c                           Phys. Rev. C67(2003) 067601
C
      REAL*8 FUNCTION DOM_INT_Wv(Ef,Ep,Av,Bv,E,n,DerivIntWv)
C
C      Analytical dispersive integral and its derivative for
C      Wv(E)=Av*(E-Ep)**n/( (E-Ep)**n + Bv**n )  for E>Ep
C      Wv(E)=Wv(2*Ef-E)                          for E<2Ef-Ep
C      Wv(E)=0                                     OTHERWISE
C
      IMPLICIT NONE
      REAL*8 Ef,Ep,Av,Bv,E,pi
      REAL*8 E0,Ex,Eplus,Emin,Rs,ResEmin,ResEplus
      REAL*8 DerEmin, DerEplus, Rds, DerivIntWv
      COMPLEX*16 Pj,I,Zj,Ztmp
      COMPLEX*8 Fs,Ds


      INTEGER N,j

      DATA I/(0.d0,1.d0)/

      pi=4.d0*atan(1.d0)

      E0 = Ep - Ef
      Ex = E  - Ef
      Eplus = Ex + E0
      Emin  = Ex - E0
      DOM_INT_Wv = 0.d0
	DerivIntWv = 0.d0

      ResEmin  =  Emin**n / (Emin**n + Bv**n)

      DerEmin  =  Emin**(n-1) *
	>           ( Emin**n + Bv**n*(1.d0 + n*log(dabs(Emin)) ) ) 
	>           / (Emin**n + Bv**n)**2 

      ResEplus = -Eplus**n / (Eplus**n + Bv**n)

      DerEplus = -Eplus**(n-1) *
	>           ( Eplus**n + Bv**n*(1.d0+n*log(Eplus)) ) 
	>           / (Eplus**n + Bv**n)**2 

C----------------------------------
C     Complex arithmetic follows
      Fs = (0.d0,0.d0)
	Ds = (0.d0,0.d0)
      do j=1,n
       Ztmp = I*(2*j-1)/dble(n)*pi
       Pj = Bv*exp(Ztmp)
       Zj = Pj * (2*Pj +Eplus -Emin) * Ex
       Zj = Zj / ( (Pj+E0) * (Pj+Eplus) * (Pj-Emin) )
       Fs = Fs + Zj*log(-Pj)
	 Ds = Ds + 2*Pj*(Ex*Ex + (Pj+E0)**2)*log(-Pj)
     >           /( (Pj+Eplus)**2 * (Pj-Emin)**2 )
      enddo

      IF(ABS(IMAG(Fs)).gt.1.e-4) STOP '(F) Too big imag part in Wv'
      Rs  = REAL(Fs)
      IF(ABS(IMAG(Ds)).gt.1.e-4) STOP '(D) Too big imag part in Wv'
	Rds = REAL(Ds)
C----------------------------------

      DOM_INT_Wv = -Av/pi*
     &  (Rs/n  + ResEplus*log(Eplus) + ResEmin*log(dabs(Emin)))

      DerivIntWv = -Av/pi*( Rds/n + DerEplus + DerEmin)

      RETURN
      END FUNCTION DOM_INT_Wv

      REAL*8 FUNCTION DOM_INT_Ws(Ef,Ep,As,Bs,Cs,E,m,DerivIntWs)
C
C      Analytical dispersive integral and its derivative for
C      Ws(E)=As*(E-Ep)**m/( (E-Ep)**m + Bs**m ) * exp(-Cs*(E-Ep)) for E>Ep
C      Ws(E)=Ws(2*Ef-E)                                           for E<2Ef-Ep
C      Ws(E)=0                                                    OTHERWISE
C
      IMPLICIT NONE
      REAL*8 Ef,Ep,As,Bs,Cs,E,EIn
      COMPLEX*16 I,Pj,Zj,Ztmp,zfi
      COMPLEX*8 Fs,Ds
      REAL*8 E0,Ex,Eplus,Emin,pi
      REAL*8 Rs,ResEmin,ResEplus
	REAL*8 DerivIntWs,DerEmin,DerEplus,Rds
      INTEGER m,j

      DATA I/(0.d0,1.d0)/

      pi=4.d0*atan(1.d0)

      E0 = Ep - Ef
      Ex = E  - Ef
      Eplus = Ex + E0
      Emin  = Ex - E0
      DOM_INT_Ws = 0.d0
	DerivIntWs = 0.d0

      ResEmin  =  Emin**m / (Emin**m + Bs**m)

      DerEmin  = -Emin**(m-1) *
	>          ( Emin**m + Bs**m + ( -Cs*Emin**(m+1) +
     >            Bs**m *(-Cs*Emin+m) ) * exp(-Cs*Emin)*EIn(Cs*Emin) ) 
	>           / (Emin**m + Bs**m)**2 

      ResEplus = -Eplus**m / (Eplus**m + Bs**m)

      DerEplus =  Eplus**(m-1) *
	>          ( Eplus**m + Bs**m + ( Cs*Eplus**(m+1) +
     >            Bs**m *(Cs*Eplus+m) ) * exp(Cs*Eplus)*EIn(-Cs*Eplus) ) 
	>           / (Eplus**m + Bs**m)**2 

C----------------------------------
C     Complex arithmetic follows
      Fs = (0.d0,0.d0)
      Ds = (0.d0,0.d0)
      do j=1,m
       Ztmp = I*(2*j-1)/dble(m)*pi
       Pj = Bs*exp(Ztmp)
       Zj = Pj * (2*Pj +Eplus -Emin) * Ex
       Zj = Zj / (Pj+E0) / (Pj+Eplus) / (Pj-Emin)
       Fs = Fs + Zj* zfi(-Pj*Cs)
	 Ds = Ds + 2*Pj*(Ex*Ex + (Pj+E0)**2)*zfi(-Pj*Cs)
     >           /( (Pj+Eplus)**2 * (Pj-Emin)**2 )
      enddo

      IF(ABS(IMAG(Fs)).gt.1.e-4) STOP '(F) Too big imag part in Wv'
      Rs  = REAL(Fs)
      IF(ABS(IMAG(Ds)).gt.1.e-4) STOP '(D) Too big imag part in Wv'
      Rds = REAL(Ds)
C----------------------------------

      DOM_INT_Ws = As/pi*(Rs/m
     &                  - ResEplus*exp(Cs*Eplus)*EIn(-Cs*Eplus)
     &                  - ResEmin*exp(-Cs*Emin)*EIn(Cs*Emin) )

      DerivIntWs = As/pi*( Rds/m + DerEplus + DerEmin)

      RETURN
      END FUNCTION DOM_INT_Ws

      real*8 function WV(A,B,Ep,Ef,E,n)
      IMPLICIT NONE
      real*8 A,B,Ep,Ef,E,ee
      integer n

      WV=0.d0
      if(E.LE.Ef) E=2.d0*Ef-E
      if(E.LT.Ep) return

      ee=(E-Ep)**n
      WV=A*ee/(ee+B**n)

      return
      end function WV

      real*8 function WDD(A,B,C,Ep,Ef,E,m)
      IMPLICIT NONE
      real*8 A,B,C,Ep,Ef,E,ee,arg
      integer m

      WDD=0.d0
      if(E.LE.Ef) E=2.d0*Ef-E
      if(E.LT.Ep) return

      arg=C*(E-Ep)
      IF(arg.GT.15) return
      ee=(E-Ep)**m
      WDD=A*ee/(ee+B**m)*EXP(-arg)
      return
      end function WDD


      REAL*8 FUNCTION DOM_int_T1(Ef,Ea,E)
C
C     Integral over E' corresponding to nonlocal additions T1(E'<<0)
C
      IMPLICIT NONE

      real*8 E,Ea,Ef,Ex,Ea2,Eax,Pi,T11,T12,T13
      Pi=4.d0*ATAN(1.d0)

      Ex=E-Ef
      Ea2=Ea**2
      Eax=Ex+Ea

      T11 = 0.5d0*log(Ea)/Ex
      T12 =  ( (2*Ea+Ex)*log(Ea)+0.5d0*pi*Ex )
     >      /(2.*(Eax**2 + Ea2))
      T13 = -Eax**2*log(Eax)/(Ex*(Eax**2+Ea2))

      DOM_int_T1 = Ex/Pi*(T11+T12+T13)
C
      RETURN
      END FUNCTION DOM_int_T1
C
      REAL*8 FUNCTION DOM_int_T2_JM(Ef,Ea,E)
C
C     Integral over E' corresponding to nonlocal additions T2(E'>>0)
C
      IMPLICIT NONE
      real*8 E,Ea,Ef,EL,Pi,DOM_int

      Pi=4.d0*ATAN(1.d0)
      EL=Ef+Ea
      DOM_int= 1.d0 / Pi * (
     >      sqrt(abs(Ef)) * atan( (2*sqrt(EL*abs(Ef)))/(EL-abs(Ef)) )
     > +    EL**1.5d0/(2*Ef)*log(Ea/EL) )

      IF(E.GT.EL) THEN

      DOM_int = DOM_int + 1.d0/Pi* (
     >  sqrt(E) * log( (sqrt(E)+sqrt(EL)) / (sqrt(E)-sqrt(EL)) ) +
     >  1.5d0*sqrt(EL)*log((E-EL)/Ea) + EL**1.5d0/(2*E)*log(EL/(E-EL)) )

      ELSEIF(E.EQ.EL) THEN

      DOM_int = DOM_int + 1.d0/Pi*1.5d0*sqrt(EL)
     > *log((2**(4.d0/3.d0)*EL)/Ea)

      ELSEIF(E.GT.0.d0 .AND. E.LE.EL) THEN

      DOM_int = DOM_int + 1.d0/Pi * (
     > sqrt(e) * log( (sqrt(E)+sqrt(EL)) / (sqrt(EL)-sqrt(E)) ) +
     > 1.5d0*sqrt(EL)*log((EL-E)/Ea)+EL**1.5d0/(2.d0*E)*log(EL/(EL-E)) )

      ELSEIF(E.EQ.0.d0) THEN

C     CPC formula 
C     DOM_int = DOM_int + 1.d0/Pi*( 0.5*EL**(1./3.)
C    > + log(EL/Ea) + 0.5d0*sqrt(EL) )

      DOM_int = DOM_int + 1.d0/Pi*1.5d0*sqrt(EL)
     > *log(EL/Ea)+0.5d0*sqrt(EL)

      ELSE

      DOM_int = DOM_int + 1.d0/Pi * (
     > -sqrt(abs(E))*atan( 2*(sqrt(EL-abs(E))) / (EL-abs(E)) ) +
     > 1.5d0*sqrt(EL)*log((EL-E)/Ea)+EL**1.5d0/(2.d0*E)*log(EL/(EL-E)) )

      ENDIF

      DOM_int_T2_JM = DOM_int

      RETURN
      END FUNCTION DOM_int_T2_JM

      REAL*8 FUNCTION DOM_int_T2(Ef,Ea,E)
C
C     Integral over E' corresponding to nonlocal additions T2(E'>>0)
C
      IMPLICIT REAL*8(A-H,O-Z)

      Pi=4.d0*ATAN(1.d0)
      El=Ef+Ea

      R1=1.5*DSQRT(El)*dLOG(abs((El-E)/Ea))

      IF(E.eq.0.d0) THEN
        R2=0.5*El**1.5d0*(1.d0/El-dlog(abs(El/Ea))/Ef)
      ELSE
        R2=0.5*El**1.5d0/(E*Ef)*
     >            ( Ef*dlog(dabs(El/(El-E))) -E*dlog(dabs(El/Ea)) )
      ENDIF

      R3=2*DSQRT(dABS(Ef))*( 0.5d0*Pi - atan( DSQRT(El/dabs(Ef)) ) )

      IF(E.GE.0.d0) THEN
        R4=DSQRT(E)*
     >           dlog(dabs( (dsqrt(El)+dsqrt(E))/(dsqrt(El)-dsqrt(E)) ))
      ELSE
        R4=-2.d0*DSQRT(dabs(E))*( 0.5d0*Pi - atan( DSQRT(dabs(El/E)) ) )
      ENDIF

      DOM_int_T2 =  1.d0/Pi*(R1+R2+R3+R4)

      RETURN
      END FUNCTION DOM_int_T2

C
C-----FUNCTION TO EVALUATE exp(Z)*E1(Z)
C
      complex*16 function zfi(za)
C
C Complex exponential integral function multiplied by exponential
C
C AUTHOR: J. Raynal
C
      IMPLICIT NONE
      real*8 aj
      complex*16 za,y
      integer m,i
      zfi=0.d0
      if (za.eq.0.) return
      if (dabs(dreal(za)+18.5d0).ge.25.d0) go to 3
      if (dsqrt(625.d0-(dreal(za)+18.5d0)**2)/1.665d0.lt.dabs(dimag(za))
     1) go to 3
      zfi=-.57721566490153d0-cdlog(za)
      y=1.d0
      do 1 m=1,2000
      aj=m
      y=-y*za/aj
      if (cdabs(y).lt.1.d-15*cdabs(zfi)) go to 2
    1 zfi=zfi-y/aj
    2 zfi=cdexp(za)*zfi
      return
    3 do 4 i=1,20
      aj=21-i
      zfi=aj/(za+zfi)
    4 zfi=aj/(1.d0+zfi)
      zfi=1.d0/(zfi+za)
      return
      end function zfi

C
C-----FUNCTION TO EVALUATE Ei(X)
C
      REAL*8 FUNCTION EIn(X)
      IMPLICIT NONE
      REAL*8 FAC, H, X
      INTEGER N
      EIn = 0.57721566490153d0+LOG(ABS(X))
      FAC = 1.0
      DO N = 1,100
      H = FLOAT(N)
      FAC = FAC*H
      EIn = EIn + X**N/(H*FAC)
      ENDDO
      RETURN
      END FUNCTION EIn

c*******************************************
      real*8 function DELTA_WV(WVf,y)
      real*8 E,Ef,A,B,Ep,y,WVf,WDE,WVE
      integer n
      common /energy/E,Ef,Ep
      common /Wenerg/WDE,WVE
      common /pdatav/A,B,n
C
      DELTA_WV=(WVf(A,B,Ep,Ef,y,n) - WVE)
     >         /((y-Ef)**2-(E-Ef)**2)
      return
      end function DELTA_WV
C
      real*8 function WVf(A,B,Ep,Ef,E,n)
      real*8 A,B,Ep,E,Ef
      integer n
C
      WVf=0.d0
      if(E.LE.Ef) E=2.d0*Ef-E
      if(E.LE.Ep) return
C
      ee=(E-Ep)**n
      WVf=A*ee/(ee+B**n)
C
      return
      end function WVf
C
      real*8 function DELTA_WD(WDf,y)
      real*8 E,Ef,A,B,C,Ep,y,WDf,WDE,WVE
      integer m,iq
      common /energy/E,Ef,Ep
      common /Wenerg/WDE,WVE
      common /pdatas/A,B,C,m,iq
C
      DELTA_WD=(WDf(A,B,C,Ep,y,m,iq) - WDE)
     >            /((y-Ef)**2-(E-Ef)**2)
C
      return
      end function DELTA_WD
C
      real*8 function WDf(A,B,C,Ep,E,m,iq)
      real*8 A,B,C,Ep,E,ee,arg
      integer m,iq
C
      WDf=0.d0
      if(E.Lt.Ep) return
      arg=C*(E-Ep)**iq
      IF(arg.GT.15) return
      ee=(E-Ep)**m
      WDf=A*ee/(ee+B**m)*EXP(-arg)
      return
      end function WDf
C
      REAL*8 FUNCTION DOM_int(DELTAF,F,Ef,Eint,Ecut,E,WE_cte)
C
C     DOM integral (20 points Gauss-Legendre)
C
C     Divided in two intervals for higher accuracy
C     The first interval corresponds to peak of the integrand
C
      DOUBLE PRECISION Eint,Ef,Ecut,WE_cte,E
      DOUBLE PRECISION F,WG,XG,WWW,XXX,DELTAF
      DOUBLE PRECISION ABSC1,CENTR1,HLGTH1,RESG1
      DOUBLE PRECISION ABSC2,CENTR2,HLGTH2,RESG2
      INTEGER J
      EXTERNAL F
      DIMENSION XG(10),WG(10)
C
C     THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C     BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C     CORRESPONDING WEIGHTS ARE GIVEN.
C
C     XG - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
C     WG - WEIGHTS OF THE 20-POINT GAUSS RULE
C
C GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
C AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
C BELL LABS, NOV. 1981.
C
      DATA WG  (  1) / 0.0176140071 3915211831 1861962351 853 D0 /
      DATA WG  (  2) / 0.0406014298 0038694133 1039952274 932 D0 /
      DATA WG  (  3) / 0.0626720483 3410906356 9506535187 042 D0 /
      DATA WG  (  4) / 0.0832767415 7670474872 4758143222 046 D0 /
      DATA WG  (  5) / 0.1019301198 1724043503 6750135480 350 D0 /
      DATA WG  (  6) / 0.1181945319 6151841731 2377377711 382 D0 /
      DATA WG  (  7) / 0.1316886384 4917662689 8494499748 163 D0 /
      DATA WG  (  8) / 0.1420961093 1838205132 9298325067 165 D0 /
      DATA WG  (  9) / 0.1491729864 7260374678 7828737001 969 D0 /
      DATA WG  ( 10) / 0.1527533871 3072585069 8084331955 098 D0 /
C
      DATA XG( 1) / 0.9931285991 8509492478 6122388471 320 D0 /
      DATA XG( 2) / 0.9639719272 7791379126 7666131197 277 D0 /
      DATA XG( 3) / 0.9122344282 5132590586 7752441203 298 D0 /
      DATA XG( 4) / 0.8391169718 2221882339 4529061701 521 D0 /
      DATA XG( 5) / 0.7463319064 6015079261 4305070355 642 D0 /
      DATA XG( 6) / 0.6360536807 2651502545 2836696226 286 D0 /
      DATA XG( 7) / 0.5108670019 5082709800 4364050955 251 D0 /
      DATA XG( 8) / 0.3737060887 1541956067 2548177024 927 D0 /
      DATA XG( 9) / 0.2277858511 4164507808 0496195368 575 D0 /
      DATA XG(10) / 0.0765265211 3349733375 4640409398 838 D0 /
C
      CENTR1 = 0.5D+00*(Ef+Eint)
      HLGTH1 = 0.5D+00*(Eint-Ef)
      CENTR2 = 0.5D+00*(Ecut+Eint)
      HLGTH2 = 0.5D+00*(Ecut-Eint)
C
C     COMPUTE THE 20-POINT GAUSS-KRONROD APPROXIMATION
C     TO THE INTEGRAL in TWO INTERVALS (Ef - Eint, Eint - Ecut)
C
      RESG1 = 0.0D+00
      RESG2 = 0.0D+00
      DO J=1,10
      XXX=XG(J)
      WWW=WG(J)
        ABSC1 = HLGTH1*XXX
      RESG1=RESG1 + WWW*(DELTAF(F,CENTR1-ABSC1)+DELTAF(F,CENTR1+ABSC1))
        ABSC2 = HLGTH2*XXX
      RESG2=RESG2 + WWW*(DELTAF(F,CENTR2-ABSC2)+DELTAF(F,CENTR2+ABSC2))
      ENDDO
c
      CORR=0.5d0*WE_cte/(E-Ef)*dlog((Ecut-(E-Ef))/(Ecut+(E-Ef)))
      DOM_int = ( RESG1*HLGTH1+RESG2*HLGTH2 + CORR) *
     >           (E-Ef)/(acos(-1.d0))
c
      RETURN
      END FUNCTION DOM_int

      subroutine VolIn(itype,depth,rad,dif,volint)	
c**********************************************************
c  Calculate volume integrals (from SCAT2000, O.Bersillon)*
c---------------------------------------------------------*
c  ITYPE  = 1 Woods-Saxon potential                       *
c           2 Woods-Saxon derivative potential            *
c  DEPTH  = depth       of the potential                  *
c  RAD    = radius      of the potential                  *
c  DIF    = diffuseness of the potential                  *
c  VOLINT = volume integral                               *
c**********************************************************
c
      implicit none
c
c
      double precision depth,rad,dif,volint
      integer          itype
c
      double precision arg
c
      double precision zero,one,three,four, pi
      data             zero  /0.0d+00/
      data             one   /1.0d+00/
      data             three /3.0d+00/
      data             four  /4.0d+00/
	data             pi    /3.1415926d0/ 
c=======================================================================
c
      volint = 0.d0 
      if(rad.eq.zero .or. depth.eq.zero) return
c
      arg = (pi*dif/rad)**2
c
      if(itype .eq. 1) then
        volint = four*pi*(rad**3)*depth*(one + arg) / three
      else if(itype .eq.  2) then
        volint = four*four*pi*dif*(rad**2)*depth*(one + arg/three)
      endif
c
      return
      end subroutine VolIn
