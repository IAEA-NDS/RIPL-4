C 29/09/06                                                      ECIS06  KHCO-000
      SUBROUTINE KHCO(NCOLT,WV,IPI,WW,H,LO)                             KHCO-001
C INPUT:     NCOLT:   NUMBER OF NUCLEAR STATES (COUPLED OR NOT) WHEN    KHCO-002
C                     CALLED FROM COLF OR NUMBER OF CONTINUA FOR        KHCO-003
C                     COMPOUND NUCLEUS WHEN CALLED FROM CONU.           KHCO-004
C            WV(J,*): MASS OF THE PARTICLE FOR J=1,                     KHCO-005
C                     MASS OF THE TARGET FOR J=2,                       KHCO-006
C                     ENERGY IN THE CENTRE OF MASS IN MEV FOR J=3.      KHCO-007
C                     ENERGY IN THE LABORATORY SYSTEM IN MEV FOR J=12.  KHCO-008
C            IPI(J,I):PRODUCT OF CHARGES FOR J=4.                       KHCO-009
C            WW:      WV FOR THE GROUND STATE.                          KHCO-010
C            H:       INTEGRATION STEP SIZE IN FERMIS.                  KHCO-011
C            LO(I):   LOGICAL CONTROLS:                                 KHCO-012
C               LO(8)  =.TRUE. RELATIVISTIC KINEMATICS.                 KHCO-013
C               LO(93) =.TRUE. NO RECOIL CORRECTION FOR REACTIONS WITH  KHCO-014
C                              SMALL DIFFERENCES (LESS THAN .5) OF      KHCO-015
C                              TARGET MASSES.                           KHCO-016
C               LO(94) =.TRUE. NON RELATIVISTIC "REDUCED MASS" FOR DIRACKHCO-017
C                              EQUATION.                                KHCO-018
C               LO(95) =.TRUE. "REDUCED ENERGY" FOR NON COULOMB         KHCO-019
C                              INTERACTION IN RELATIVISTIC SCHROEDINGER KHCO-020
C                              EQUATION OR USE OF REST MASS IN DIRAC    KHCO-021
C                              EQUATION.                                KHCO-022
C               LO(96) =.TRUE. "REDUCED MASS" FOR COULOMB INTERACTION   KHCO-023
C                              IN RELATIVISTIC SCHROEDINGER EQUATION.   KHCO-024
C               LO(97) =.TRUE. SAME REDUCED MASS FOR ALL THE STATES     KHCO-025
C                              WITH SMALL DIFFERENCES (LESS THAN .5)    KHCO-026
C                              OF TARGET MASSES.                        KHCO-027
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    KHCO-028
C OUTPUT:    WV(J,*): WAVE NUMBER IN 1/FERMI FOR J=4,                   KHCO-029
C                     COULOMB PARAMETER FOR J=5,                        KHCO-030
C                     REDUCED MASS FOR J=6.                             KHCO-031
C                     REDUCED ENERGY FOR J=7.                           KHCO-032
C                     STEP SIZE FOR THIS LEVEL FOR J=8.                 KHCO-033
C                     SQUARE ROOT OF COEF. OF SCALAR POTENTIALS FOR J=9.KHCO-034
C                     SQUARE ROOT OF COEF. OF COUL. POTENTIALS FOR J109.KHCO-035
C                     K MULTIPLIED BY RATIO OF STEP SIZES FOR J=11.     KHCO-036
C                     REDUCED ENERGY TERM FOR J=12.                     KHCO-037
C                     SET TO 0. FOR J=14 TO 19.                         KHCO-038
C                                                                       KHCO-039
C FOR THE COMMON  /DCONS/ SEE CALC.                                     KHCO-040
C                                                                       KHCO-041
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     KHCO-042
C  CM:        ATOMIC MASS IN MEV.                                       KHCO-043
C  CHB:       PLANCK CONSTANT /(2*PI) IN MEV*FERMI.                     KHCO-044
C  CMB:       ATOMIC MASS CM DIVIDED BY H-BAR*C.                        KHCO-045
C  CCZ:       COULOMB ALPHA CONSTANT.                                   KHCO-046
C  CK:        H-BAR*C.                                                  KHCO-047
C   USED:     CM,CHB,CCZ,CK.                                            KHCO-048
C                                                                       KHCO-049
C***********************************************************************KHCO-050
      IMPLICIT REAL*8 (A-H,O-Z)                                         KHCO-051
      LOGICAL LO(150)                                                   KHCO-052
      DIMENSION IPI(11,*),WV(22,*),WW(22)                               KHCO-053
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            KHCO-054
      AMT1=DMAX1(WW(1),WW(2))                                           KHCO-055
      AMR=WW(3)/CM+WW(1)+WW(2)                                          KHCO-056
C WAVE NUMBER,COULOMB PARAMETER AND CALL TO COULOMB SUBROUTINES.        KHCO-057
      DO 6 I=1,NCOLT                                                    KHCO-058
      AMT2=DMAX1(WV(1,I),WV(2,I))                                       KHCO-059
      RC=1.D0                                                           KHCO-060
      IF ((.NOT.LO(93)).OR.(DABS(AMT1-AMT2).GT.0.5D0)) RC=AMT1/AMT2     KHCO-061
      HREC=H*RC                                                         KHCO-062
      IF ((.NOT.LO(97)).OR.(DABS(AMT1-AMT2).GT.0.5D0)) GO TO 1          KHCO-063
      AMP=WW(1)                                                         KHCO-064
      AMT=WW(2)                                                         KHCO-065
      GO TO 2                                                           KHCO-066
    1 AMP=WV(1,I)                                                       KHCO-067
      AMT=WV(2,I)                                                       KHCO-068
    2 IF (LO(8)) GO TO 3                                                KHCO-069
C NON RELATIVISTIC KINEMATICS: NO MASS CORRECTION.                      KHCO-070
      AMRE=AMP*AMT/(AMP+AMT)                                            KHCO-071
      AMRM=AMRE                                                         KHCO-072
      WSK2=CK*WV(3,I)*AMRM                                              KHCO-073
      GO TO 5                                                           KHCO-074
    3 WSK2=0.125D0*CK*WV(3,I)*(WV(3,I)/CM+2.D0*WV(1,I)+2.D0*AMT)*(WV(3,IKHCO-075
     1)/CM+2.D0*WV(1,I))*(WV(3,I)/CM+2.D0*AMT)/AMR**2                   KHCO-076
      WV(4,I)=DSQRT(DABS(WSK2))                                         KHCO-077
      AMRM=AMP*AMT/AMR                                                  KHCO-078
      IF (LO(109)) GO TO 4                                              KHCO-079
C RELATIVISTIC KINEMATICS FOR SCHROEDINGER EQUATION.                    KHCO-080
      AMRD=(AMR**4-(WV(1,I)**2-AMT**2)**2)/(4.D0*AMR**3)                KHCO-081
      AMRE=AMRD                                                         KHCO-082
      IF (LO(96)) AMRE=AMRM                                             KHCO-083
      IF (LO(95)) AMRM=AMRE                                             KHCO-084
      GO TO 5                                                           KHCO-085
    4 IF (LO(94)) AMRM=AMP*AMT/(AMP+AMT)                                KHCO-086
      IF (LO(95)) AMRM=AMP                                              KHCO-087
      AMRE=DSQRT(WSK2/CMB**2+AMRM**2)                                   KHCO-088
      IF (LO(96)) AMRE=DSQRT(WSK2/CMB**2+AMP**2)                        KHCO-089
    5 WV(4,I)=DSQRT(DABS(WSK2))                                         KHCO-090
      WV(5,I)=CM*CCZ*AMRE*IPI(4,I)/WV(4,I)/CHB**2                       KHCO-091
      WV(6,I)=CM*AMRM                                                   KHCO-092
      WV(7,I)=CM*AMRE                                                   KHCO-093
      WV(8,I)=HREC                                                      KHCO-094
      WV(9,I)=HREC*DSQRT(CK*AMRM)                                       KHCO-095
      WV(10,I)=HREC*DSQRT(CK*AMRE)                                      KHCO-096
      WV(11,I)=WV(4,I)*RC                                               KHCO-097
      WV(12,I)=HREC*HREC*WSK2                                           KHCO-098
      IF (DABS(WV(5,I)).GT.400.D0) WV(5,I)=400.D0*WV(5,I)/DABS(WV(5,I)) KHCO-099
    6 CONTINUE                                                          KHCO-100
      RETURN                                                            KHCO-101
      END                                                               KHCO-102
