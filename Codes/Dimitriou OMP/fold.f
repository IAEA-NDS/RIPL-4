      program alphapot
c
c alpha potential = real double folding potential + imaginary WS potential
c dispersive relations also included if ipot=3
c nztar: Z of target 
c natar: A of target 
c e: incident energy
c
      implicit double precision(a-h,o-z)
      parameter(ngrid=1000,ndis=400,nzmax=130,nnmax=300)
      common/calpha/excst(nzmax,nnmax)
      common/calpha1/e2exp
      dimension va(ngrid),rva(ngrid)
      dimension vopr(ndis),wopr(ndis),vsopr(ndis)
      dimension v(10)
      dimension cdef(nzmax,nnmax),hdef(nzmax,nnmax),b2th(nzmax,nnmax)
      dimension ronth(nzmax,nnmax),anth(nzmax,nnmax)
      dimension rozth(nzmax,nnmax),azth(nzmax,nnmax)
      dimension cnth(nzmax,nnmax),czth(nzmax,nnmax)

      nztar=62
      natar=144
      nntar=natar-nztar
      e=3.82988
      e=10.
      ee=e
      ipot=1
      pi=3.141592654
      a=natar
      z=nztar
      a13=a**(1./3.)
      a2=a*a
      a3=a2*a
c
c---------------------------------------------------------------------------
c  determination of the density parameters
c---------------------------------------------------------------------------
      open(19,file='/nuc/Masses/hfb14-conf')
      read(19,*)
  650 continue
      read(19,660,end=670,err=650) iz,ia,c,h,b2,b4,cn1,tn1,
     &  an1,cz1,tz1,az1
  660 format(2i4,4f6.2,38x,6f9.4)
      tn1=tn1/float(ia)**(1./3.)
      tz1=tz1/float(ia)**(1./3.)
      in=ia-iz
      if (iz.gt.nzmax.or.in.gt.nnmax) goto 650
      b2th(iz,in)=b2
      cnth(iz,in)=cn1
      czth(iz,in)=cz1
      cdef(iz,in)=c
      hdef(iz,in)=h
      ronth(iz,in)=tn1
      anth(iz,in)=an1
      rozth(iz,in)=tz1
      azth(iz,in)=az1
      goto 650
  670 continue
      close(19)
      rn=ronth(nztar,nntar)*a13
      bn=anth(nztar,nntar)
      rz=rozth(nztar,nntar)*a13
      bz=azth(nztar,nntar)
      cn=cnth(nztar,nntar)
      cz=czth(nztar,nntar)
c
c---------------------------------------------------------------------------
c  first excited state determined from nld calculation for the Energy dependence of 
c     the imaginary component of the alpha-potential
c---------------------------------------------------------------------------
      do 800 kz=1,nzmax
      do 800 kn=1,nnmax
  800 excst(kz,kn)=0.
c
      open(8,file='/home/sgoriely/Most/Excstate/excstate.hfb14')
      read(8,*)
      read(8,*)
  810 continue
      read(8,820,end=830) kz,ka,e1
  820 format(2i4,f8.4)
      kn=ka-kz
      if (kn.gt.nnmax) goto 810
      if (kz.gt.nzmax) goto 810
      excst(kz,kn)=e1
      goto 810
  830 continue
      close(8)
      e2exp=excst(nztar,natar-nztar)
c
      goto (350,360,370) ipot
c-----------------------------------------------------------------------
c potential fitted on experimental data  --> only Imaginary part is read
c Double-folding real potential
c Woods-Saxon e-dependent volume imaginary potential
c grama & demetriou & goriely, double folding potential 
c-----------------------------------------------------------------------
  350 continue
      v(1)=1.25
      e0=e2exp
c version 23/01/02
      rw=0.85281+0.02202*a-2.14551e-4*a2+7.92942e-7*a3-
     1     9.94658e-10*a2*a2
      aw=-0.13526+0.02029*a-1.98441e-4*a2+7.35104e-7*a3-
     2     9.15272e-10*a2*a2
      if (aw.le.1.e-2) aw=1.e-2
      as=7.65867-7.5669*e0+2.50486*e0**2
      es=0.1024*a+1.1307*as
      te=1./(1.+exp(-(e-es)/as))
      expj0=77.
      if (a.le.90.) expj0=135.-0.644*a
      ww=3./3.141592/rw**3*expj0*te
      v(3)=rw
      v(4)=aw
      v(8)=ww
      v(5)=0.
      dwv=0.
      dws=0.
      goto 2000
c----------------------------------------------------------------------
c potential fitted on experimental data  --> only Imaginary part is read
c Double-folding real potential
c Woods-Saxon e-dependent volume+ surface imaginary potential
c demetriou & grama & goriely, double folding potential 
c-----------------------------------------------------------------------
  360 continue
      v(1)=1.25
      e0=e2exp
cv version 23/01/02 for volume+surface term
      rw=1.47385-0.00134615*a
      aw=0.29
      v(5)=0.9
      as=7.65867-7.5669*e0+2.50486*e0**2
      es=0.1024*a+1.1307*as
      te=1./(1.+exp(-(e-es)/as))
      expj0=77.
      if(a.le.90.) expj0=135.-0.644*a
cv 7.603=4*1.6*1.09**2
      ww=expj0*te/(rw**3/3.+7.603*v(5)*aw*rw**2/a13)/3.141592
      v(3)=rw
      v(4)=aw
      v(8)=ww
      dwv=0.
      dws=0.
      goto 2000
c-----------------------------------------------------------------------
c potential fitted on experimental data  --> only Imaginary part is read
c Double-folding real potential+dispersive contributions
c using prescription and functions of Capote et al. JPG:NPP 27 (2001) B15
c Woods-Saxon e-dependent volume+ surface imaginary potential
c demetriou & grama & goriely, double folding potential
c-----------------------------------------------------------------------
  370 continue
      v(1)=1.25
      e0=e2exp
cv version 13/12/01 for volume+surface term
      rw=1.47385-0.00134615*a
      aw=0.32
      v(5)=0.5
      as=7.65867-7.5669*e0+2.50486*e0**2
      es=0.1024*a+1.1307*as
      te=1./(1.+exp(-(e-es)/as))
      expj0=77.
      if(a.le.90.) expj0=135.-0.644*a
      ww=expj0*te/(rw**3/3.+7.603*v(5)*aw*rw**2/a13)/3.141592
      v(3)=rw
      v(4)=aw
      v(8)=ww
cv dispersive contributions 
cv volume part
      dw=0.
      Ei=efermi
      Ecut=1000.
      Eint=2.*es
      Eint=efermi+Eint
      expj=ww
      WVFE=WVF(expj,as,es,efermi,ee)
      dw=2*DOM_int(Delta_WVF,WVF,efermi,Eint,Ecut,ee,Ei,WVFE)
      dwv=dw
cv surface part incl. damping function
      dw=0.
      C=0.01951-0.00049571*ee
      if(ee.lt.14.) C=0.49236-0.03427*ee
      expj=v(5)*ww
      WDFE=WDF(expj,as,C,es,efermi,ee)
      dw=2*DOM_int(Delta_WDF,WDF,efermi,Eint,Ecut,ee,Ei,WDFE)
      dws=dw
cv add damping term to surface potential
      v(5)=v(5)*exp(-C*dabs(e-efermi))
 2000 continue

c

      rww=v(3)*a13
      aww=v(4)
      www=v(8)

      rb=20.
      nu=500

      call fold(z,a,ee,rb,nu,rva,va,rz,rn,bz,bn,cz,cn)

c
c interpolation of the alpha optical potential on a given radial grid
c
      m4=200
      d1=rb/dble (m4)
      do 460 k=1,m4
      rr=k*d1
      klo=1
      khi=nu
      if (rr.le.rva(klo)) then
         khi=klo+1
         goto 465
      endif
      if (rr.ge.rva(khi)) then
         klo=nu-1
         goto 465
      endif
  470 if (khi-klo.gt.1) then
        kk=(khi+klo)/2.
        if (rva(kk).gt.rr) then
          khi=kk
        else
          klo=kk
        endif
        goto 470
      endif
  465 continue
      hh=rva(khi)-rva(klo)
      if (abs(hh).lt.1.e-8) then
        write(*,*) ' Problem in interpolation'
        stop
      endif
      aa1=(rva(khi)-rr)/hh
      bb1=(rr-rva(klo))/hh
      vopr(k)=aa1*va(klo)+bb1*va(khi)
      wopr(k)=-www/(1.+dexp((rr-rww)/aww))
cv if v(5).eq.0 only volume; if v(5).ne. 0 volume+surface 
      ws=0. 
      rws=1.09*rww
      aws=1.6*aww
      ws=-4.*www*v(5)*dexp((rr-rws)/aws)/
     1(1.+dexp((rr-rws)/aws))**2
      wopr(k)=wopr(k)+ws
cv dispersive contributions to real potential if dwv.ne.0 and dws.ne.0
      vdiv=0.
      vdis=0.
      vdiv=-dwv/(1.+dexp((rr-rww)/aww))
      vdis=-4.*dws*dexp((rr-rws)/aws)/
     1(1.+dexp((rr-rws)/aws))**2
      vopr(k)=vopr(k)+vdiv+vdis
      vsopr(k)=0.
      write(*,461) ee,k,rr,vopr(k),wopr(k)
  461 format(f9.5,i4,3f10.3)
  460 continue
      end
c
c---------------------------------------------------------------------------
c subroutine for double folding potential for alpha projectiles
c---------------------------------------------------------------------------
c
      SUBROUTINE fold(Z,A,E,rb,mp,rv,U,rz,rn,bz,bn,cz,cn)
c
c Z: Z of target
c A: A of target
c E: incident energy
c rb: maximum radius value, typically rb=20.
c mp= nbr of points for double folding (maximum=ngrid=1000; typicaly mp=500)
c rv: radius grid for the alpha potential
c U: alpha potential in a vector form at radius rv
c rz,rn,bz,bn,cz,cn: density parameters (Fermi function) for protons and neutrons in target A
c
      implicit double precision(a-h,o-z)
      parameter(ngrid=1000)
      common/calpha1/e2exp
      COMMON /SP/ DQ,D1,D2,NR(2),NQ,NU,DU,FRO(ngrid,2),RO(ngrid,2),
     *AM(4),ROE(ngrid,2),BETA
      COMMON/CONSTF/ PI,SQRPI,PI4,PI32                                  
      COMMON/ISO/ ISO 
      common/voli/volj
      DIMENSION BCOF(ngrid),BINT(9),B(1),UE(ngrid)
      DIMENSION U(ngrid),rv(ngrid)
      DIMENSION A1(6),A2(6),B1(6),B2(6)                      
      logical pass
      save
      EQUIVALENCE (AM(1),AMP),(AM(2),AMT),(NR(1),N1),(NR(2),N2)                                           
      DATA A1/1577.1071D0,34751.04D0,19839.156D0,25129.599D0,25129.599D0
     *,-15348.25D0/,B1/0.16,16.0,16.,16.,16.,16./,A2/-1239.9272D0,      
     * -12754.865D0,-9857.0604D0,-10726.653D0,-10726.653D0,5908.71D0/,  
     * B2/0.0625,6.2500,6.25,6.25,6.25,6.25/                            
      data pass/.false./
      PI=3.1415926536D0                                                 
      SQRPI=DSQRT(PI)                                                   
      PI4=4.*3.1415926536D0                                             
      PI32=SQRPI*PI
      nu=mp
      KEY=5                                                       
c     ALPHA=3.55D-03*E+4.01
c     BETA=-0.0305*E+11.5
c     CE=-1.203D-03*E+0.477
c e-independent parameters
      ce=0.44073
      alpha=4.3529
      beta=10.639
      AM(1)=4
      AM(2)=A
      AM(3)=2
      AM(4)=Z 
      R1=rb
      D1=R1/nu
      R2=rb
      D2=R2/nu
c id=1 if d1=d2
      id=1
      RQ=10.
      DQ=RQ/nu
      RU=rb                                        
      NR(1)=nu
      NR(2)=nu
      NQ=nu
      DU=RU/nu
      DO 34 I=1,NU                                                      
      U(I)=0.D0                                                         
      UE(I)=0.D0                                                        
  34  CONTINUE                                                          
      AKOFF=0.5/(PI*PI)                                                 
      IS1=1                                                             
      IS2=1
      ANMZP=DABS(AMP-2.*AM(3))
      ANMZT=DABS(AMT-2.*AM(4))
c     IF(ANMZP.LT.0.00001D0.OR.ANMZT.LT.0.00001D0) IS2=1
      DO 33 ISO=IS1,IS2
      IF(ISO.EQ.2) KEY=6                                                
      DE=0.                                                             
      GOTO (20,20,20,20,15,16),KEY                                      
   15 DE=-276.*(1.-0.005*E/AMP)                                         
c     de=-272.5
      GOTO 20                                                           
   16 DE=227*(1.-0.004*E/AMP)   
   20 CONTINUE  
      IF (pass) GO TO 25                                                        
      NO=1                                                              
      IF(AMP.GT.1.) CALL DENSE(NO)
      pass=.true.
   25 continue
      NO=2                                                              
      R=0.                      
      DO 42 I=1,nu
      RO(I,2)=(cn/(1.+DEXP((R-rn)/bn))+cz/(1.+DEXP((R-rz)/bz)))*r*r*pi4
      ROE(I,2)=RO(I,NO)*DEXP(-RO(I,NO)*BETA)   
   42 R=R+D2
      NINT=MAX0(NR(1),NR(2))                                            
      IF(AMP.EQ.1.) NINT=NR(2)                                          
      IFLG=0                                                            
      Q=0.                                                              
      KF=1                                                       
      DO 200 K=1,NQ                                                     
      QQ=Q*Q                                            
      SUM1=0.                               
      SUM2=0.   
      S1=0.
      S2=0. 
      i1=1
      z0=0.
      ndd=nu
      IF(id.eq.1.AND.IFLG.EQ.0)            
     *         CALL FTRANS(BCOF,i1,NINT,z0,D1,Q,BINT,B,ndd)               
      IF(id.ne.1.AND.IFLG.EQ.0)            
     *         CALL FTRANS(BCOF,i1,N1,z0,D1,Q,BINT,B,ndd)                 
      IF(KEY.NE.5)GOTO 101
      DO 100 J=1,N1                         
  100 S1=S1+BCOF(J)*ROE(J,1)   
  101 DO 102 J=1,N1                         
  102 SUM1=SUM1+BCOF(J)*RO(J,1)
  126    IF(id.ne.1.AND.IFLG.EQ.0)
     *      CALL FTRANS(BCOF,i1,N2,z0,D2,Q,BINT,B,ndd)
      IF(KEY.NE.5) GOTO 131
      DO 130 J=1,N2
  130 S2=S2+BCOF(J)*ROE(J,2)
  131 DO 132 J=1,N2
  132 SUM2=SUM2+BCOF(J)*RO(J,2)
  172 C=A1(KEY)/(QQ+B1(KEY))+A2(KEY)/(QQ+B2(KEY))+DE
      FRO(K,1)=(SUM1*SUM2+S1*S2*ALPHA)*C*QQ
  200 Q=Q+DQ
      R=0.
      DO 300 I=1,NU
      SUM=0.
      CALL FTRANS(BCOF,i1,NQ,z0,DQ,R,BINT,B,ndd)
      DO 400 K=1,NQ
  400 SUM=SUM+BCOF(K)*FRO(K,1)
      U(I)=AKOFF*SUM*CE+U(I)
      rv(i)=r
  300 R=R+DU
  33  CONTINUE
      R=0
      DO 520 I=1,NU
      UD=U(I)
      U(I)=UD+UE(I)  
  520 R=R+DU
      CALL INTU(U,NU,DU,VOLJ,RMS,vmi33,3)
      VOLJ=VOLJ/(AMP*AMT)
      RMS=DSQRT(RMS)
c possible renormalization on the Volume Integral
      rap=1.
      DO 521 I=1,NU
      U(I)=rap*U(I)
  521 continue
      RETURN
      END
C
      SUBROUTINE INTU(U,NR,DU,VOLJ,RMS,RIN,KEY)
      implicit double precision(a-h,o-z)
      parameter(ngrid=1000)
      DIMENSION U(ngrid)
      FACT=0.375*DU*12.5663706144D0
      VOLJ=0.
      RMS=0.
      RIN=0.
      SUM0=0.
      SUM1=0.
      SUM2=0.
      A0=0.
      A1=0.
      A2=0.
      R=0.
      NR2=NR-1
      NR1=(NR2/3)*3
      DO 400 KR=3,NR1,3
      R1=R+DU
      R2=R1+DU
      R3=R2+DU
   2  RR1=R1*R1
      RR2=R2*R2
      RR3=R3*R3
      B1=U(KR-1)*RR1+U(KR)*RR2
      C1=U(KR+1)*RR3
      SUM1=SUM1+A1+3.*B1+C1
      A1=C1
      RRR1=RR1*RR1
      RRR2=RR2*RR2
      RRR3=RR3*RR3
      B2=U(KR-1)*RRR1+U(KR)*RRR2
      C2=U(KR+1)*RRR3
      SUM2=SUM2+A2+3.*B2+C2
      A2=C2
   1  R=R3
  400 CONTINUE
      IF(NR2.LE.NR1) GOTO 3
      NR1=NR1+1
      C5=DU*4.D0/3.D0
      DO 4 KR=NR1,NR2
      SUM0=SUM0+C5*(U(KR)*(KR-1)+U(KR+1)*KR)  
      SUM1=SUM1+C5*DU*(U(KR)*(KR-1)**2+U(KR+1)*KR**2)     
      SUM2=SUM2+C5*DU**3*(U(KR)*(KR-1)**4+U(KR+1)*KR**4)   
    4 CONTINUE
    3 RIN=FACT*SUM0
      VOLJ=FACT*SUM1
      IF(SUM1.NE.0.) RMS=SUM2/SUM1
      RETURN
      END  
C
      SUBROUTINE DENSE(N0)
      implicit double precision(a-h,o-z)
      parameter(ngrid=1000)                                         
      COMMON /SP/ DQQ,D1,D2,NR(2),NQQ,NU,DU,FRO(ngrid,2),RO(ngrid,2)
     *,AM(4),ROE(ngrid,2),BETA
      COMMON /CONSTF/ PI,SQRPI,PI4,PI32                        
      COMMON /ISO/ ISO  
      dimension a(12),ri(12),sq(12)
      DIMENSION FR(2*ngrid),F(2*ngrid)
      DIMENSION ROH(ngrid),FH(ngrid)
      EQUIVALENCE (FRO(1,1),F(1)),(RO(1,1),FR(1))
      NO=N0                    
      NQ=NQQ                                                            
      DQ=DQQ                                                            
      NH=0                                                              
      NFH=0                                                             
      AH=0.                                                             
      RMSH=0.                                                           
      TISO=1.                                    
      IF(ISO.GT.1) TISO=-1.   
      DO 1 I=1,ngrid
      RO(I,NO)=0.                       
      ROE(I,NO)=0.                  
      ROH(I)=0. 
      FRO(I,NO)=0.                                                      
   1  FH(I)=0.                                                          
      MF=1+(NO-1)*nu
      MR=1+(NO-1)*nu
      N=NR(NO)                
  700 AMH=AM(NO)-AH                                                     
      IF(NO.EQ.2)GO TO 50       
      uu=2./3.
      g=dsqrt(uu)
      ri(1)=0.2
      ri(2)=0.6
      ri(3)=0.9
      ri(4)=1.4
      ri(5)=1.9
      ri(6)=2.3
      ri(7)=2.6
      ri(8)=3.1
      ri(9)=3.5
      ri(10)=4.2
      ri(11)=4.9
      ri(12)=5.2
      sq(1)=0.034724
      sq(2)=0.430761
      sq(3)=0.203166
      sq(4)=0.192986
      sq(5)=0.083866
      sq(6)=0.033007
      sq(7)=0.014201
      sq(8)=0.
      sq(9)=0.006860
      sq(10)=0.
      sq(11)=0.000438
      sq(12)=0.
      do 11 i=1,12
   11 a(i)=4.*sq(i)/(2.*pi32*g**3*(1.+2.*(ri(i)/g)**2))
      r=0.
      sum=0.
      do 3 j=1,nu
      s=0.
      do 22 i=1,12
   22 s=s+a(i)*(dexp(-((r-ri(i))/g)**2)+dexp(-((r+ri(i))/g)**2))
      ro(j,1)=s 
      r=r+D1 
    3 continue        
      R=0.                                                              
      ROE(1,1)=0.                                                      
      RO(1,1)=0.                                                       
      DO 80 I=2,N                                                       
      RR=R*R                                                            
      ROE(I,1)=PI4*RR*RO(I,NO)*DEXP(-RO(I,NO)*BETA)   
      RO(I,1)=RO(I,NO)*RR*PI4                                          
   80 R=R+D1 
      RETURN                        
   50 CONTINUE
      RETURN                                                            
      END                                                               
C
      SUBROUTINE FTRANS(BCOF,LM,NPOINT,R,H,P,BINT,B,NDIM)
      implicit double precision(a-h,o-z)
      DIMENSION BCOF(NDIM,LM),BINT(LM,3,3),F(2),BFAC(2),LL(2),
     1 FI(3),B(LM),BDIF(3)
      DATA EPS,EPS1,EPS4 /1.D-3,1.D-13,1.D-1/
      L1=LM-1
      L2=LM-2
      D=H*P
      DFLM=DFAC(L1)
c      dflm1=1.
      IF(LM.GT.1) DFLM1=DFAC(L2)
      XEPS=0.
      IF(L1.GT.0.) XEPS=(EPS*DFLM)**(1./DFLOAT(L1))
      IF(DABS(D).LT.EPS4) GOTO 300
      D2=D+D
      DD=D*D
      XFAC=1./(P*DD)
      X=R*P
      N=(NPOINT+1)/2
      DO 5 IN=1,N
      XX=X*X
      M=MOD(IN,3)+1
      M1=MOD(IN-1,3)+1
      M2=MOD(IN-2,3)+1
      IF(X.EQ.0) GOTO 8
      IF(LM.EQ.1.OR.DABS(X).GE.XEPS) GOTO 100
      IF(DABS(X).GT.EPS1) GOTO 11
    8 DO 9 J=1,3
      DO 9 I=1,LM
    9 BINT(I,J,M)=0.
      GOTO 200
   11 DO 10 I=1,2
      LL(I)=LM-I
      F(I)=1.
   10 BFAC(I)=1.
      DO 15 I2=1,3
   15 FI(I2)=1./(LL(1)+I2)
      IND=0
   20 IND=IND+1
      DO 30 I=1,2
      BFAC(I)=-BFAC(I)/(4.*IND*(LL(I)+IND+0.5))*XX
   30 F(I)=F(I)+BFAC(I)
      DO 35 I2=1,3
   35 FI(I2)=FI(I2)+BFAC(1)/(LL(1)+I2+IND+IND)
      IF(DABS(BFAC(2)).GT.EPS1) GOTO 20
      XEX=X**L2*XFAC
      B(LM)=XEX*X*F(1)/DFLM
      B(L1)=XEX  *F(2)/DFLM1
      DO 60 I2=1,3
   60 BINT(LM,I2,M)=XEX*X**(I2+1)*FI(I2)/DFLM
      BINT(L1, 1,M)=(BINT(LM,2,M)+X*B(L1))/L1
      BINT(L1, 2,M)=L1*BINT(LM,1,M)+X*B(LM)
      BINT(L1, 3,M)=L2*BINT(LM,2,M)+XX*B(LM)
      IF(LM.LE.2) GOTO 200
      LI=L1
      DO 40 L=1,L2
      LI=LI-1
      B(LI)=(LI+LI+1)/X*B(LI+1)-B(LI+2)
      BINT(LI,1,M)=((LI+1)*BINT(LI+2,1,M)+(LI+LI+1)
     * *B(LI+1))/LI
      BINT(LI, 2,M)=LI*BINT(LI+1,1,M)+X*B(LI+1)
   40 BINT(LI,3,M)= (LI-1)*BINT(LI+1,2,M)+XX*B(LI+1)
      GOTO 200
  100 S=DSIN(X)
      C=DCOS(X)
      B(1)=S/X*XFAC
      B(2)=(B(1)-C*XFAC)/X
      BINT(1,1,M)=SININT(X)*XFAC
      BINT(1,2,M)=(-C+1.)*XFAC
      BINT(1,3,M)=(S-X*C)*XFAC
      IF(LM.EQ.1) GOTO 200
      BINT(2,1,M)=-B(1)+XFAC
      BINT(2,2,M)=BINT(1,1,M)-S*XFAC
      BINT(2,3,M)=(-2.*C-X*S+2.)*XFAC
      IF(LM.EQ.2) GOTO 200
      DO 110 L=3,LM
      B(L)=(L+L-3)/X*B(L-1)-B(L-2)
      BINT(L,1,M)=((L-2)*BINT(L-2,1,M)-(L+L-3)*B(L-1))/(L-1)
      BINT(L,2,M)= (L-1)*BINT(L-1,1,M)-X*B(L-1)
  110 BINT(L,3,M)= L*BINT(L-1,2,M)-XX*B(L-1)
  200 CONTINUE
      IF(IN.EQ.1) GOTO 5
      X1=X-D
      X2=X-D2
      DO 230 L=1,LM
      DO 210 I2=1, 3
  210 BDIF(I2)=BINT(L,I2,M)-BINT(L,I2,M1)
      BCOF(2*IN-2,L)=-BDIF(3)+2.*X1*BDIF(2)
     *  -(X1*X1-DD)*BDIF(1)
      IF(IN.GT.2) GOTO 240
      BCOF(1,L)=0.5*BDIF(3)-(X2+1.5*D)*BDIF(2)
     *  +((.5*X2+1.5*D)*X2+DD)*BDIF(1)
      GOTO 230
  240 BCOF(2*IN-3,L)=.5*(BINT(L,3,M)-BINT(L,3,M2))
     1  -X2*(BINT(L,2,M)-BINT(L,2,M2))
     2  -1.5*D*(BINT(L,2,M)-2.*BINT(L,2,M1)+BINT(L,2,M2))
     3  +(.5*X2*X2+DD)*(BINT(L,1,M)-BINT(L,1,M2))
     4  +1.5*X2*D*(BINT(L,1,M)-2.*BINT(L,1,M1)+BINT(L,1,M2))
      IF(IN.EQ.N) BCOF(NPOINT,L)=.5*BDIF(3)-(X-1.5*D)
     1  *BDIF(2)+((.5*X-1.5*D)*X+DD)*BDIF(1)
  230 CONTINUE
    5 X=X+D2
      RETURN
  300 CONTINUE
      XFAC=H/3.
      X=P*R
      IVF=4
      DO 340 IN=1,NPOINT
      IF(X.EQ.0.0) GOTO 301
      IF(LM.EQ.1.OR.DABS(X).GE.XEPS) GOTO 400
  301 IVF=6-IVF
      IF(DABS(X).GE.EPS1) GOTO 360
      BCOF(IN,1)=IVF*XFAC
      IF(LM.EQ.1) GOTO 340
      DO 370 I1=2,LM
  370 BCOF(IN,I1)=0.
      GOTO 340
  360 XX=X*X
      DO 310 I=1,2
      LL(I)=LM-I
      F(I)=1.
  310 BFAC(I)=1.
      IND=0
  320 IND=IND+1
      DO 330 I=1,2
      BFAC(I)=-BFAC(I)/(4.*IND*(LL(I)+IND+.5))*XX
  330 F(I)=F(I)+BFAC(I)
      IF(DABS(BFAC(2)).GT.EPS1) GOTO 320
      XEX=X**L2*XFAC
      BCOF(IN,LM)=XEX*X*F(1)*IVF/DFLM
      BCOF(IN,L1)=XEX  *F(2)*IVF/DFLM1
      IF(LM.LE.2) GOTO  340
      LI=L1
      DO 345 L=1,L2
      LI=LI-1
  345 BCOF(IN,LI)=(LI+LI+1)/X*BCOF(IN,LI+1)-BCOF(IN,LI+2)
  340 X=X+D
      GOTO 410
  400 N1=IN
      DO 351 IN=N1,NPOINT
      IVF=6-IVF
      S=DSIN(X)*XFAC
      C=DCOS(X)*XFAC
      BCOF(IN,1)=S/X*IVF
      IF(LM.EQ.1) GOTO 351
      BCOF(IN,2)=(BCOF(IN,1)-C*IVF)/X
      IF(LM.EQ.2) GOTO 351
      DO 350 L=3,LM
  350 BCOF(IN,L)=(L+L-3)/X*BCOF(IN,L-1)-BCOF(IN,L-2)
  351 X=X+D
  410 DO 420 L=1,LM
      BCOF(1,L)=.5*BCOF(1,L)
  420 BCOF(NPOINT,L)=.5*BCOF(NPOINT,L)
      RETURN
      END
C
      FUNCTION DFAC(L)
      implicit double precision(a-h,o-z)
      DFAC=1.
      IF(L.EQ.0) RETURN
      X=1.
      DO 10 I=1,L
      X=X+2.
   10 DFAC=DFAC*X
      RETURN
      END
      FUNCTION SININT(X)
c$large
      implicit double precision(a-h,o-z)
      IF(DABS(X).GE.16) GOTO 1
      Z=X/16.
      Y=4.*Z**2-2.
      B=       0.00000 00000 00007D0
      A= Y*B  -0.00000 00000 00185D0
      B= Y*A-B+0.00000 00000 04185D0
      A= Y*B-A-0.00000 00000 84710D0
      B= Y*A-B+0.00000 00015 22370D0
      A= Y*B-A-0.00000 00241 00076D0
      B= Y*A-B+0.00000 03329 88589D0
      A= Y*B-A-0.00000 39729 08746D0
      B= Y*A-B+0.00004 04202 71419D0
      A= Y*B-A-0.00034 54691 55569D0
      B= Y*A-B+0.00243 62214 04749D0
      A= Y*B-A-0.01386 74455 89417D0
      B= Y*A-B+0.06203 36794 32003D0
      A= Y*B-A-0.21126 37809 76555D0
      B= Y*A-B+0.53014 88479 16522D0
      A= Y*B-A-0.96832 22369 87086D0
      B= Y*A-B+1.38930 87711 71888D0
      A= Y*B-A-1.92656 50911 50656D0
      B= Y*A-B+2.77875 63817 42663D0
      A= Y*B-A-4.06398 08449 11986D0
      A= Y*A-B+8.10585 29553 61245D0
      SININT=Z*.5*(A-B)
      RETURN
    1 Z=16./X
      G=Z**2
      Y=4.*G-2.
      B=       0.00000 00000 00002D0
      A= Y*B  -0.00000 00000 00014D0
      B= Y*A-B+0.00000 00000 00107D0
      A= Y*B-A-0.00000 00000 00964D0
      B= Y*A-B+0.00000 00000 10308D0
      A= Y*B-A-0.00000 00001 36096D0
      B= Y*A-B+0.00000 00023 56196D0
      A= Y*B-A-0.00000 00586 70317D0
      B= Y*A-B+0.00000 24537 55677D0
      A= Y*B-A-0.00023 37560 41393D0
      A= Y*A-B+0.12452 74580 57854D0
      F= Z*.5*(A-B)
      B=       0.00000 00000 00002D0
      A= Y*B  -0.00000 00000 00012D0
      B= Y*A-B+0.00000 00000 00087D0
      A= Y*B-A-0.00000 00000 00717D0
      B= Y*A-B+0.00000 00000 06875D0
      A= Y*B-A-0.00000 00000 79604D0
      B= Y*A-B+0.00000 00011 69202D0
      A= Y*B-A-0.00000 00234 68225D0
      B= Y*A-B+0.00000 07249 95950D0
      A= Y*B-A-0.00004 26441 82622D0
      A= Y*A-B+0.00772 57121 93407D0
      G= G*.5*(A-B)
      B= 1.57079 63267 94897D0
      IF(X.LT.0.) B=-B
      SININT=B-F*DCOS(X)-G*DSIN(X)
      RETURN
      END
C
      FUNCTION FINT(Y,R,R0,H,NMAX)
      implicit double precision(a-h,o-z)
      DIMENSION Y(1)
      N1=IDINT((R-R0+1.D-8)/H)
      IF(N1.EQ.0) N1=1
      IF((N1+2).EQ.NMAX ) N1=N1-1
      N2=N1+1
      N3=N2+1
      N4=N3+1
      X1=R0+DFLOAT((N1-1))*H
      X2=X1+H
      X3=X2+H
      X4=X3+H
      T1=R-X1
      T2=R-X2
      T3=R-X3
      T4=R-X4
      S=-T2*T3*T4*Y(N1)/3.+T1*T3*T4*Y(N2)-T1*T2*T4*Y(N3)+               
     * T1*T2*T3*Y(N4)/3.
      HH=H*H
      H3=HH*H
      FINT=0.5*S/H3
      RETURN
      END
cv function calculating integrand with Fermi-type volume W(E)
      real*8 function DELTA_WVF(WVF,y)
      real*8 E,Ef,Ei,expj,es,as,C,y,WVF,WDFE,WVFE,dwv,dws
      common /Wenergy/E,es,Ei
      common /eferm/Ef
      common /Wenerg/WDFE,WVFE,dwv,dws
      common /pdatav/expj,as,C
      external WVF

      DELTA_WVF=(WVF(expj,as,es,Ef,y) - WVFE)
     >            /((y-Ei)**2-(E-Ei)**2)

C     DELTA_WD=(WD(A,B,C,Ep,y) -
C    >          WD(A,B,C,Ep,E))
C    >          /((y-Ef)**2-(E-Ef)**2)

      return
      end
cv function for Fermi-type volume W(E)
      real*8 function WVF(A,B,Ep,Ef,E)
      real*8 A,B,E,Ep,Ef,ee

      WVF=0.d0
      if(E.Lt.Ef) return
      ee=E-Ep
      WVF=A/(1.+dexp(ee/B))
      return
      end
cv function calculating integrand for Fermi-type surface W(E)
      real*8 function DELTA_WDF(WDF,y)
      real*8 E,Ef,Ei,A,B,C,Ep,y,WDF,WDFE,WVFE,dwv,dws
      common /Wenergy/E,Ep,Ei
      common /eferm/Ef
      common /Wenerg/WDFE,WVFE,dwv,dws
      common /pdatav/A,B,C
      external WDF

      DELTA_WDF=(WDF(A,B,C,Ep,Ef,y) - WDFE)
     >            /((y-Ei)**2-(E-Ei)**2)

C     DELTA_WD=(WD(A,B,C,Ep,y) -
C    >          WD(A,B,C,Ep,E))
C    >          /((y-Ef)**2-(E-Ef)**2)

      return
      end
cv function for Fermi-type surface W(E)
      real*8 function WDF(A,B,C,Ep,Ef,E)
      real*8 A,B,C,Ep,Ef,E,ee,arg

      WDF=0.d0
      if(E.le.Ef) return
      arg=C*dabs(E-Ef)
      IF(arg.GT.15) return
      ee=E-Ep
      WDF=A*EXP(-arg)/(1.+dexp(ee/B))
      return
      end
c function performing Gaussian integration
      REAL*8 FUNCTION DOM_int(DELTAF,F,Ef,Eint,Ecut,E,Ei,WE_cte)
C
C     DOM integral (20 points Gauss-Legendre)
C
C     Divided in two intervals for higher accuracy
C     The first interval corresponds to peak of the integrand
C
      DOUBLE PRECISION Eint,Ef,Ecut,WE_cte,E,Ei
      DOUBLE PRECISION F,WG,XG,WWW,XXX,DELTAF
      DOUBLE PRECISION ABSC1,CENTR1,HLGTH1,RESG1
      DOUBLE PRECISION ABSC2,CENTR2,HLGTH2,RESG2
      INTEGER J
      EXTERNAL F,DELTAF
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

	CORR=0.5d0*WE_cte/(E-Ei)*dlog((Ecut-(E-Ei))/(Ecut+(E-Ei))) 
      DOM_int = ( RESG1*HLGTH1+RESG2*HLGTH2 + CORR) * 
     >           (E-Ei)/(4.d0*ATAN(1.d0))

      RETURN
      END


