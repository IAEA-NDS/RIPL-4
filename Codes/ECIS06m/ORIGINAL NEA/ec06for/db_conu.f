C 21/08/06                                                      ECIS06  CONU-000
      SUBROUTINE CONU(IX,IPI,WV,IPIM,WVM,NCI,XD,SCN,KXT,LO)             CONU-001
C DISCRETISATION OF CONTINUA FOR COMPOUND NUCLEUS.                      CONU-002
C INPUT:     IX:      0 FOR COUNT OF POINTS,                            CONU-003
C                     NOT 0 FOR COMPUTATION OF POINTS, WEIGHTS, ...     CONU-004
C            IPI:     INTEGER VALUES FOR THE DESCRIPTION OF CHANNELS    CONU-005
C                     (SEE CALX).                                       CONU-006
C            WV:      FLOATING VALUES FOR THE DESCRIPTION OF CHANNELS   CONU-007
C                     (SEE CALX).                                       CONU-008
C            SCN:     DESCRIPTIONS OF LEVEL DENSITIES:                  CONU-009
C                     1-SA  2-UX   3-TAU  4-SG   5-E0   6-EX   7-NZ     CONU-010
C                     (SEE INPUT DESCRIPTION).                          CONU-011
C            LO(I):   LOGICAL CONTROLS:                                 CONU-012
C               LO(8)  =.TRUE. RELATIVISTIC KINEMATICS.                 CONU-013
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     CONU-014
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    CONU-015
C OUTPUT:    IPIM:    IPI FOR CONTINUA OF COMPOUND NUCLEUS.             CONU-016
C            WVM:     WV FOR CONTINUA OF COMPOUND NUCLEUS.              CONU-017
C            NCI:     STARTING AND FINAL ADDRESSES FOR CONTINUA.        CONU-018
C            XD:      STEPS OF DISCRETISATION OF CONTINUA,              CONU-019
C                     ENERGY AND SPIN DEPENDENCE OF LEVEL DENSITIES     CONU-020
C            KXT:     NUMBER OF TRANSMISSION COEFFICIENTS               CONU-021
C                                                                       CONU-022
C FOR THE COMMON  /DCONS/ SEE CALC.                                     CONU-023
C FOR THE COMMON  /NCOMP/ SEE CALX.                                     CONU-024
C                                                                       CONU-025
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     CONU-026
C  CM:        ATOMIC MASS IN MEV.                                       CONU-027
C   USED:     CM.                                                       CONU-028
C                                                                       CONU-029
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     CONU-030
C  NCONT:     NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS.                 CONU-031
C  NCONS:     NUMBER OF LEVEL DENSITIES NEEDED.                         CONU-032
C  NIE:       NUMBER OF UNCOUPLED STATES ADDED FOR DISCRETISATION.      CONU-033
C  NCOLX:     TOTAL NUMBER OF LEVELS WITHOUT DISCRETISATION.            CONU-034
C  ACN1:      RATIO SIZE/STEP FOR DISCRETISATION OF A CONTINUUM.        CONU-035
C  ACN2:      MAXIMUM NUMBER OF STEPS BY MEV FOR A CONTINUUM.           CONU-036
C   DEFINED:  NIE.                                                      CONU-037
C   USED:     NCONT,NCONS,NIE,NCOLX,ACN1,ACN2.                          CONU-038
C***********************************************************************CONU-039
      IMPLICIT REAL*8 (A-H,O-Z)                                         CONU-040
      LOGICAL LO(150)                                                   CONU-041
      DIMENSION IPI(11,*),WV(22,*),IPIM(11,*),WVM(22,*),NCI(2,*),XD(3,*)CONU-042
     1,SCN(7,*)                                                         CONU-043
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            CONU-044
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQCONU-045
     1,ACN1,ACN2,AZ(18)                                                 CONU-046
      COMMON /INOUT/ MR,MW,MS                                           CONU-047
C DISCRETISATION OF CONTINUUM FOR COMPOUND NUCLEUS.                     CONU-048
      NIE=0                                                             CONU-049
      IF (IX.EQ.0) GO TO 1                                              CONU-050
      KXT=IPI(8,NCOLX+1-NCONT)                                          CONU-051
    1 DO 12 I=1,NCONT                                                   CONU-052
      J=I+NCONS-NCONT                                                   CONU-053
      SCNB=DEXP(2.D0*DSQRT(SCN(1,J)*SCN(2,J))-(SCN(6,J)-SCN(5,J))/SCN(3,CONU-054
     1J))*SCN(3,J)/(SCN(1,J)**.5D0*SCN(2,J)**1.5D0)                     CONU-055
      IV=I-NCONT+NCOLX                                                  CONU-056
      NCI(1,I)=NIE+1                                                    CONU-057
      IF (WV(3,IV).LE.0.D0) GO TO 11                                    CONU-058
      NE=0                                                              CONU-059
      E=WV(3,IV)                                                        CONU-060
      ACN3=1.D0/ACN2                                                    CONU-061
      IF (E*ACN2.LT.1.D0) ACN3=E/1.9D0                                  CONU-062
    2 IF (E.LE.0.D0) GO TO 11                                           CONU-063
      NIE=NIE+1                                                         CONU-064
      IF (NE.NE.0) GO TO 4                                              CONU-065
      IF (E.LT.ACN1*ACN3) GO TO 3                                       CONU-066
      ESTP=E/ACN1                                                       CONU-067
      GO TO 5                                                           CONU-068
    3 NE=IDINT(E/ACN3+1.D0)                                             CONU-069
      IF (NE.EQ.1) NE=2                                                 CONU-070
      ACN3=E/DFLOAT(NE)*1.000001D0                                      CONU-071
    4 ESTP=ACN3                                                         CONU-072
    5 IF (E.LT.ESTP) ESTP=E                                             CONU-073
      ECN=E-.5D0*ESTP                                                   CONU-074
      E=E-ESTP                                                          CONU-075
      IF (IX.EQ.0) GO TO 2                                              CONU-076
      DO 6 K=1,22                                                       CONU-077
    6 WVM(K,NIE)=WV(K,IV)                                               CONU-078
      DO 7 K=1,11                                                       CONU-079
    7 IPIM(K,NIE)=IPI(K,IV)                                             CONU-080
      WVM(3,NIE)=ECN                                                    CONU-081
      XD(2,NIE)=ESTP                                                    CONU-082
      IF (LO(8)) GO TO 8                                                CONU-083
      WVM(13,NIE)=WVM(3,NIE)*(WVM(1,NIE)+WVM(2,NIE))/WVM(2,NIE)         CONU-084
      GO TO 9                                                           CONU-085
    8 WVM(13,NIE)=WVM(3,NIE)*(WVM(3,NIE)/(2.D0*CM)+WVM(1,NIE)+WVM(2,NIE)CONU-086
     1)/WVM(2,NIE)                                                      CONU-087
    9 LY=4+IDINT(3.3D0*DSQRT(WVM(3,NIE)))                               CONU-088
      IPIM(10,NIE)=LY-1                                                 CONU-089
      IPIM(8,NIE)=KXT                                                   CONU-090
      KXT=KXT+LY*IPIM(2,NIE)                                            CONU-091
      IF (WV(3,1)-WVM(3,NIE).GT.SCN(6,J)) GO TO 10                      CONU-092
      XD(1,NIE)=XD(2,NIE)*DEXP((WV(3,1)-WVM(3,NIE)-SCN(5,J))/SCN(3,J))/SCONU-093
     1CN(3,J)                                                           CONU-094
      XD(3,NIE)=2.D0*SCN(4,J)**2                                        CONU-095
      GO TO 2                                                           CONU-096
   10 EXD=WV(3,1)-WVM(3,NIE)+SCN(2,J)-SCN(6,J)                          CONU-097
      EBYT=DSQRT(SCN(1,J)*EXD)                                          CONU-098
      XD(1,NIE)=XD(2,NIE)*DEXP(2.D0*EBYT)/(EBYT*EXD*SCNB)               CONU-099
      XD(3,NIE)=2.D0*SCN(4,J)**2*DSQRT(EXD/SCN(2,J))                    CONU-100
      GO TO 2                                                           CONU-101
   11 IF (IX.EQ.0) GO TO 12                                             CONU-102
      JY=NIE-NCI(1,I)+1                                                 CONU-103
      WRITE (MW,1000) I,JY                                              CONU-104
      NCI(2,I)=NIE                                                      CONU-105
   12 CONTINUE                                                          CONU-106
      IF (IX.NE.0) CALL KHCO(NIE,WVM,IPIM,WV,WV(8,1),LO)                CONU-107
      RETURN                                                            CONU-108
 1000 FORMAT (' ***** START OF',I4,'TH CONTINUUM *****',5X,I5,' DISCRETICONU-109
     1SATION POINTS.')                                                  CONU-110
      END                                                               CONU-111
