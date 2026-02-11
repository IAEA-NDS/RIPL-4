      MODULE ext_functions
      
      PUBLIC
      
      CONTAINS
      
      REAL*8 FUNCTION DELTA_WV(WVf,y)

      REAL*8  E,Ef,A,B,Ep,y,WVf,WDE,WVE,Ea
      INTEGER n
      COMMON /energy/E,Ef,Ep,Ea
      COMMON /Wenerg/WDE,WVE
      COMMON /pdatav/A,B,n

      DELTA_WV = (WVf(A,B,Ep,Ef,y,n) - WVE)
     >           /((y-Ef)**2-(E-Ef)**2)
      RETURN
      END FUNCTION DELTA_WV
      
      REAL*8 FUNCTION WVf(A,B,Ep,Ef,E,n)

      REAL*8  A,B,Ep,E,Ef
      INTEGER n

      WVf = 0.d0
      IF (E.LE.Ef) E = 2.d0*Ef-E
      IF (E.LE.Ep) RETURN

      ee = (E-Ep)**n
      WVf = A*ee/(ee+B**n)

      RETURN
      END FUNCTION WVf      
      
      REAL*8 FUNCTION DELTA_WD(WDf,y)

      REAL*8 E,Ef,A,B,C,Ep,y,WDf,WDE,WVE,Ea
      INTEGER m,iq
      COMMON /energy/E,Ef,Ep,Ea
      COMMON /Wenerg/WDE,WVE
      COMMON /pdatas/A,B,C,m,iq

      DELTA_WD = (WDf(A,B,C,Ep,y,m,iq) - WDE)
     >            /((y-Ef)**2-(E-Ef)**2)

      RETURN
      END FUNCTION DELTA_WD
      
      REAL*8 FUNCTION WDf(A,B,C,Ep,E,m,iq)
      REAL*8 A,B,C,Ep,E,ee,arg
      INTEGER m,iq

      WDf = 0.d0
      IF (E.Lt.Ep) RETURN
      arg = C*(E-Ep)**iq
      IF (arg.GT.15) RETURN
      ee = (E-Ep)**m
      WDf = A*ee/(ee+B**m)*EXP(-arg)
      
      END FUNCTION WDf      

c==   START integrands in r for the various terms
      real*8 function VolumeR(r) 
      IMPLICIT real*8 (A-H,O-Z)
      common/parameters/vdepth,radius,diffuss
      common/armon/LANDA,LANDAM
      VolumeR=vdepth/(1+exp((r-radius)/diffuss))*r**(landa+2)
      return
      end   function VolumeR

      real*8 function SurfaceR(r) 
      IMPLICIT real*8 (A-H,O-Z)
      common/parameters/vdepth,radius,diffuss
      common/armon/LANDA,LANDAM
      SurfaceR = 4./(1+exp((r-radius)/diffuss))**2*
     > vdepth*exp((r-radius)/diffuss)*r**(landa+2)
      return
      end   function SurfaceR

      real*8 function SpinOrbitR(r) 
      IMPLICIT real*8 (A-H,O-Z)
      common/parameters/vdepth,radius,diffuss
      common/armon/LANDA,LANDAM
C     SpinOrbitR =  -vdepth/(diffuss*radius)*exp((r-radius)/diffuss)
C    >             /  (1+exp((r-radius)/diffuss))**2*r**(landa+2)
      SpinOrbitR =  -vdepth/diffuss*exp((r-radius)/diffuss)
     >             /  (1+exp((r-radius)/diffuss))**2*r**(landa+1)
      return
      end   function SpinOrbitR
      
      END MODULE ext_functions
