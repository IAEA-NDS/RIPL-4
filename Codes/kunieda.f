Ccc   * $Rev: 2537 $
Ccc   * $Author: rcapote $
Ccc   * $Date: 2012-02-10 14:07:34 +0100 (Fri, 10 Feb 2012) $

      SUBROUTINE KUNIEDA(izat)
      IMPLICIT REAL*8(A-H,O-Z)
      common /omp/VR0,VR1,VR2,VR3,VRLA,ALAVR,
     +            WD0,WD1,WDA1,WDBW,WDWID,ALAWD,
     +            WC0,WC1,WCA1,WCBW,WCWID,BNDC,
     +            VS,ALASO,WSO,WS1,WSBW,WSWID,
     +            RR,RRBWC,RRWID,PDIS,ARO,AR1,
     +            RD,AD0,AD1,
     +            RC,AC0,AC1,
     +            RW,AW0,AW1,
     +            RS,ASO,AS1,
     +            RZ,RZBWC,RZWID,AZ,CCOUL,ALF,
     +            CISO,WCISO

c     target mass and charge
c     -----------------------
      AT    = mod(izat,1000)
      AT    = dble(AT)
      AT3   = 1.d0/(AT**0.3333333333)
      IZT   = izat/1000

c     real volume depth : V_R^0
c     --------------------------
      VR0 = -34.40d0 - 6.18d-2*AT + 3.37d-4*AT**2 - 6.52d-7*AT**3
      if(IZT.eq.12) VR0   = -39.0d0
      if(IZT.eq.13) VR0   = -39.0d0
      if(IZT.eq.14) VR0   = -39.0d0
      if(IZT.eq.20) VR0   = -37.5d0
      if(IZT.eq.26) VR0   = -37.0d0
      if(IZT.eq.28) VR0   = -37.0d0
      if(IZT.eq.29) VR0   = -37.0d0
      if(IZT.eq.39) VR0   = -37.5d0
      if(IZT.eq.40) VR0   = -37.5d0
      if(IZT.eq.41) VR0   = -38.0d0
      if(IZT.eq.42) VR0   = -38.0d0
      if(IZT.eq.50) VR0   = -38.0d0
      if(IZT.eq.58) VR0   = -38.5d0
      if(IZT.eq.73) VR0   = -38.5d0
      if(IZT.eq.74) VR0   = -38.5d0
      if(IZT.eq.79) VR0   = -38.5d0
      if(IZT.eq.82) VR0   = -38.5d0
      if(IZT.eq.83) VR0   = -38.5d0
      if(IZT.eq.90) VR0   = -38.5d0
      if(IZT.eq.92) VR0   = -39.0d0

c     real volume depth : V_R^1
c     --------------------------
      VR1   = 0.027d0

c     real volume depth : V_R^2
c     --------------------------
      VR2   = 0.000120d0

c     real volume depth : V_R^3
c     --------------------------
      VR3   = 0.00000035d0

c     real volume depth : V_R^{DISP}
c     -------------------------------
      VRLA  = 94.88d0

c     lambda for real volume depth : lambda_R
c     ----------------------------------------
      ALAVR = 1.05d-2 + 1.63d-4*VR0
      if(IZT.eq.13) ALAVR = 0.00430d0
      if(IZT.eq.20) ALAVR = 0.00450d0
      if(IZT.eq.26) ALAVR = 0.00450d0
      if(IZT.eq.28) ALAVR = 0.00450d0
      if(IZT.eq.29) ALAVR = 0.00450d0
      if(IZT.eq.73) ALAVR = 0.00423d0
      if(IZT.eq.74) ALAVR = 0.00423d0
      if(IZT.eq.79) ALAVR = 0.00423d0
      if(IZT.eq.82) ALAVR = 0.00425d0
      if(IZT.eq.83) ALAVR = 0.00425d0
      if(IZT.eq.90) ALAVR = 0.00425d0
      if(IZT.eq.92) ALAVR = 0.00417d0
c     ------------------------------------------

c     imaginary surface depth : W_D^{DISP}
c     --------------------------------------
      WDBW  = 11.08d0 + 6.83d-5*AT**2
      if(IZT.eq.12) WDBW  = 11.5d0
      if(IZT.eq.13) WDBW  = 11.0d0
      if(IZT.eq.14) WDBW  = 11.0d0
      if(IZT.eq.20) WDBW  = 11.0d0
      if(IZT.eq.26) WDBW  = 11.0d0
      if(IZT.eq.28) WDBW  = 11.0d0
      if(IZT.eq.29) WDBW  = 11.5d0
      if(IZT.eq.39) WDBW  = 11.5d0
      if(IZT.eq.40) WDBW  = 11.5d0
      if(IZT.eq.41) WDBW  = 12.0d0
      if(IZT.eq.42) WDBW  = 12.0d0
      if(IZT.eq.50) WDBW  = 12.8d0
      if(IZT.eq.58) WDBW  = 12.5d0
      if(IZT.eq.73) WDBW  = 13.0d0
      if(IZT.eq.74) WDBW  = 13.0d0
      if(IZT.eq.79) WDBW  = 13.5d0
      if(IZT.eq.82) WDBW  = 14.0d0
      if(IZT.eq.83) WDBW  = 14.0d0
      if(IZT.eq.90) WDBW  = 15.0d0
      if(IZT.eq.92) WDBW  = 15.0d0

c     imaginary surface depth : WID_D
c     --------------------------------------
      WDWID = 12.72d0 - 1.54d-2*AT + 7.14d-5*AT**2
      if(IZT.eq.12) WDWID  = 10.0d0
      if(IZT.eq.13) WDWID  = 10.0d0
      if(IZT.eq.14) WDWID  = 10.0d0
      if(IZT.eq.20) WDWID  = 10.0d0
      if(IZT.eq.26) WDWID  = 12.0d0
      if(IZT.eq.28) WDWID  = 12.0d0
      if(IZT.eq.29) WDWID  = 12.0d0
      if(IZT.eq.39) WDWID  = 12.0d0
      if(IZT.eq.40) WDWID  = 12.0d0
      if(IZT.eq.41) WDWID  = 12.0d0
      if(IZT.eq.42) WDWID  = 12.0d0
      if(IZT.eq.50) WDWID  = 12.0d0
      if(IZT.eq.58) WDWID  = 12.0d0
      if(IZT.eq.73) WDWID  = 12.0d0
      if(IZT.eq.74) WDWID  = 12.0d0
      if(IZT.eq.79) WDWID  = 12.0d0
      if(IZT.eq.82) WDWID  = 13.0d0
      if(IZT.eq.83) WDWID  = 13.0d0
      if(IZT.eq.90) WDWID  = 13.0d0
      if(IZT.eq.92) WDWID  = 13.0d0

c     lambda for imaginary surface depth : lambda_D
c     -----------------------------------------------
      ALAWD = 0.0140d0

c     imaginary volume depth : W_V^{DISP}
c     -------------------------------------
      WCBW  = 17.0d0

c     imaginary volume depth : WID_V
c      ------------------------------------
      WCWID = 105.05d0
     +      - 14.13d0/(1.d0+exp((AT-100.21d0)/13.47d0))
      if(IZT.le.14) WCWID =  90.0d0
      if(IZT.eq.20) WCWID =  90.0d0
      if(IZT.eq.26) WCWID =  90.0d0
      if(IZT.eq.28) WCWID =  90.0d0
      if(IZT.eq.29) WCWID =  95.0d0
      if(IZT.eq.39) WCWID =  95.0d0
      if(IZT.eq.40) WCWID =  95.0d0
      if(IZT.eq.41) WCWID =  98.0d0
      if(IZT.eq.42) WCWID =  95.0d0
      if(IZT.eq.50) WCWID = 103.0d0
      if(IZT.eq.58) WCWID = 105.0d0
      if(IZT.eq.73) WCWID = 105.0d0
      if(IZT.eq.74) WCWID = 105.0d0
      if(IZT.eq.79) WCWID = 105.0d0
      if(IZT.eq.82) WCWID = 105.0d0
      if(IZT.eq.83) WCWID = 105.0d0
      if(IZT.eq.90) WCWID = 105.0d0
      if(IZT.eq.92) WCWID = 105.0d0


c     real, imaginary spin-orbit depth (the global K-D OMP)
c     ------------------------------------------------------
      VS    = 5.922d0 + 0.003d0 * AT
      ALASO = 0.005d0
      WSBW  = -3.1d0
      WSWID = 160.d0





c     radius parameter for real-volume, imaginary-surface and imajinary-volume
c     potential (no mass dipendence)
c     ---------------------------------------------------------------------------
      RR    = 1.210d0

c     diffuseness parameter for real-volume, imaginary-surface and imajinary-volume
c     --------------------------------------------------------------------------------
      ARO = 0.49d0 + 3.22d-3*AT - 2.07d-5*AT**2 + 4.47d-8*AT**3 
      if(IZT.eq.12) ARO   = 0.560d0
      if(IZT.eq.13) ARO   = 0.520d0
      if(IZT.eq.14) ARO   = 0.540d0
      if(IZT.eq.20) ARO   = 0.610d0
      if(IZT.eq.26) ARO   = 0.610d0
      if(IZT.eq.28) ARO   = 0.620d0
      if(IZT.eq.29) ARO   = 0.620d0
      if(IZT.eq.39) ARO   = 0.650d0
      if(IZT.eq.40) ARO   = 0.635d0
      if(IZT.eq.41) ARO   = 0.650d0
      if(IZT.eq.42) ARO   = 0.650d0
      if(IZT.eq.50) ARO   = 0.655d0
      if(IZT.eq.58) ARO   = 0.650d0
      if(IZT.eq.73) ARO   = 0.660d0
      if(IZT.eq.74) ARO   = 0.660d0
      if(IZT.eq.79) ARO   = 0.660d0
      if(IZT.eq.82) ARO   = 0.671d0
      if(IZT.eq.83) ARO   = 0.670d0
      if(IZT.eq.90) ARO   = 0.685d0
      if(IZT.eq.92) ARO   = 0.685d0

c     radius parameter for spin-orbit potential (the K-D global omp)
c     ----------------------------------------------------------------
      RS    = 1.18d0 - 0.65d0 * AT3

c     diffuseness parameter for spin-orbit potential (the K-D global omp)
c     ---------------------------------------------------------------------
      ASO   = 0.59d0

c     Coulomb radius parameter
c     -------------------------
      RZ    = 1.264d0

c     diffuseness parameter for the Coulomb potential
c     ------------------------------------------------
      AZ    = 0.341d0





c     Coulomb coefficient : C_{coul}
c     --------------------------------
      CCOUL = 0.9d0

c     isospin coefficient for real volume : C_{viso}
c     ------------------------------------------------
      CISO  = 24.3d0

c     isospin coefficient for imaginary surface : C_{wiso}
c     -----------------------------------------------------
      WCISO = 18.0d0
      return
      END



