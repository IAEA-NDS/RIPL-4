C 31/08/06                                                      ECIS06  COLF-000
      SUBROUTINE COLF(NCOLT,NCOLL,IPI,WV,ISM,LMAX1,LMAX2,H,LM,NIV,LO,FG,COLF-001
     1XG,LMAX3,KXT,IEXP,Z)                                              COLF-002
C COULOMB FUNCTIONS AT THE MATCHING POINT RM=ISM*H.                     COLF-003
C INPUT:     NCOLT:   NUMBER OF NUCLEAR STATES (COUPLED OR NOT)         COLF-004
C                     AND CONTINUA FOR COMPOUND NUCLEUS.                COLF-005
C            NCOLL:   NUMBER OF COUPLED STATES.                         COLF-006
C            IPI(J,I):MULTIPLICITY OF PARTICLE AND TARGET FOR J=1 AND 2,COLF-007
C                     PRODUCT OF CHARGES FOR J=4,                       COLF-008
C                     NUMBER OF ANGULAR MOMENTA FOR I=1 AND J=10.       COLF-009
C            WV(J,*): MASS OF PARTICLE AND TARGET FOR J=1,2,            COLF-010
C                     CENTRE OF MASS ENERGY IN MEV FOR J=3.             COLF-011
C            ISM:     NUMBER OF POINTS FOR INTEGRATION.                 COLF-012
C            LMAX1:   MAXIMUM NUMBER OF COULOMB FUNCTIONS.              COLF-013
C            LMAX2:   MAXIMUM NUMBER OF COULOMB PHASE SHIFTS.           COLF-014
C            H:       INTEGRATION STEP SIZE IN FERMIS.                  COLF-015
C            LM:      LENGTH OF WORKING SPACE Z.                        COLF-016
C            NIV:     ADDRESSES OF COULOMB INTEGRALS IN NIV(*,*,3).     COLF-017
C            LO(I):   LOGICAL CONTROLS:                                 COLF-018
C               LO(8)  =.TRUE. RELATIVISTIC KINEMATICS.                 COLF-019
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     COLF-020
C               LO(108)=.TRUE. DIAGONAL COULOMB CORRECTIONS ARE NEEDED. COLF-021
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    COLF-022
C               LO(133)=.TRUE. STORE SCALAR AND COULOMB POTENTIAL       COLF-023
C                              INDEPENDENTLY.                           COLF-024
C OUTPUT:    IPI(J,I):STARTING ADDRESS OF PENETRABILITIES FOR I=NCOLL+1 COLF-025
C                     TO NCOLT AND J=8.                                 COLF-026
C                     NUMBER OF ANGULAR MOMENTA FOR J=10.               COLF-027
C            WV(J,*): OUTPUT OF SUBROUTINE KHCO FOR J=4 TO 11.          COLF-028
C            FG:      IN FG(I,J,K), COULOMB WAVE FUNCTIONS F, FP, G, GP COLF-029
C                     RESPECTIVELY FOR J=1 TO 4, FOR ANGULAR MOMENTUM   COLF-030
C                     I-1 AND CHANNEL K. IF COULOMB CORRECTIONS ARE     COLF-031
C                     REQUESTED, INTEGRALS OF F*F, F*G, G*F AND G*G     COLF-032
C                     DIVIDED BY R**2 FROM MATCHING RADIUS TO INFINITY  COLF-033
C                     WITH THE SAME ANGULAR MOMENTUM I-1 RESPECTIVELY   COLF-034
C                     FOR J+1 TO 4 AND K GIVEN BY NIV(*,*,3).           COLF-035
C            XG:      IN XG(I,K), SIGMA(I-1)-SIGMA(0) FOR ANGULAR       COLF-036
C                     MOMENTUM I-1 AND CHANNEL K. IF COULOMB CORRECTIONSCOLF-037
C                     ARE REQUESTED, INTEGRALS OF F*F DIVIDED BY R**2   COLF-038
C                     FROM ZERO TO INFINITY AND K GIVEN BY NIV(*,*,3).  COLF-039
C            LMAX3:   EFFECTIVE MAXIMUM NUMBER OF COULOMB FUNCTIONS.    COLF-040
C            KXT:     NUMBER OF PENETRABILITIES FOR UNCOUPLED STATES.   COLF-041
C WORKING AREA:                                                         COLF-042
C            IEXP:    POWERS OF 10 FOR LARGE VALUES OF FUNCTIONS        COLF-043
C                     (MULTIPLES OF 10**15).                            COLF-044
C            Z:       WORKING AREA FOR SUBROUTINE CORI.                 COLF-045
C                                                                       COLF-046
C FOR THE COMMON  /DCONS/ SEE CALC.                                     COLF-047
C                                                                       COLF-048
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     COLF-049
C  CM:        ATOMIC MASS IN MEV.                                       COLF-050
C  CHB:       PLANCK CONSTANT /(2*PI) IN MEV*FERMI.                     COLF-051
C  CCZ:       COULOMB ALPHA CONSTANT.                                   COLF-052
C  CK:        H-BAR*C.                                                  COLF-053
C  XZ:        CONVERSION FACTOR TO MILLIBARNS.                          COLF-054
C   DEFINED:  XZ.                                                       COLF-055
C   USED:     CM,CHB,CCZ,CK.                                            COLF-056
C                                                                       COLF-057
C***********************************************************************COLF-058
      IMPLICIT REAL*8 (A-H,O-Z)                                         COLF-059
      LOGICAL LO(150)                                                   COLF-060
      DIMENSION IPI(11,*),WV(22,*),FG(LMAX1,4,*),XG(LMAX2,*),IEXP(*),Z(*COLF-061
     1),NIV(NCOLL,NCOLL,3)                                              COLF-062
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            COLF-063
      COMMON /INOUT/ MR,MW,MS                                           COLF-064
      KXT=0                                                             COLF-065
      LO(133)=.FALSE.                                                   COLF-066
      WRITE (MW,1000)                                                   COLF-067
C CHECK OF THE LENGTH OF WORKING FIELD FOR IEXP.                        COLF-068
      IF (LM.LT.LMAX1) CALL MEMO('COLF',LM,LMAX1)                       COLF-069
      LMAX3=LMAX1                                                       COLF-070
C WAVE NUMBER,COULOMB PARAMETER AND CALL TO COULOMB SUBROUTINES.        COLF-071
      CALL KHCO(NCOLT,WV,IPI,WV,H,LO)                                   COLF-072
      XZ=10.D0/DFLOAT(IPI(2,1)*IPI(3,1))/WV(4,1)**2                     COLF-073
      RM=ISM*H                                                          COLF-074
      DO 6 I=1,NCOLT                                                    COLF-075
      LO(133)=LO(133).OR.((WV(5,I).NE.0.D0).AND.(WV(9,I).NE.WV(10,I)))  COLF-076
      IF (I.GT.NCOLL) GO TO 5                                           COLF-077
      WRITE (MW,1001) I,(WV(J,I),J=4,8)                                 COLF-078
      RAU=RM*WV(11,I)                                                   COLF-079
      ETA=WV(5,I)                                                       COLF-080
      IF (WV(3,I).GT.0.D0) GO TO 1                                      COLF-081
      CALL COCL(ETA,RAU,LMAX1-1,FG(1,1,I),FG(1,2,I),FG(1,3,I),FG(1,4,I),COLF-082
     1XG(1,I))                                                          COLF-083
      GO TO 2                                                           COLF-084
    1 LX=LMAX1                                                          COLF-085
      CALL FCOU(LX-1,ETA,RAU,FG(1,1,I),FG(1,2,I),FG(1,3,I),FG(1,4,I),IEXCOLF-086
     1P,XG(1,I))                                                        COLF-087
    2 IPI(10,I)=IPI(10,1)                                               COLF-088
      IF (WV(3,I).LT.0.D0) GO TO 6                                      COLF-089
C COMPUTATION OF COULOMB PHASE-SHIFTS.                                  COLF-090
      LMX1=LMAX3                                                        COLF-091
C CORRECTION OF LARGE VALUE AND SEARCH FOR MAXIMUM EFFECTIVE NUMBER     COLF-092
C OF COULOMB FUNCTIONS.                                                 COLF-093
      DO 3 J=1,LMX1                                                     COLF-094
      IF (IEXP(J).EQ.0) GO TO 3                                         COLF-095
      IF (LO(108).OR.IEXP(J).GT.15) LMAX3=MIN0(LMAX3,J)                 COLF-096
      FG(J,1,I)=FG(J,1,I)*1.D-15                                        COLF-097
      FG(J,2,I)=FG(J,2,I)*1.D-15                                        COLF-098
      FG(J,3,I)=FG(J,3,I)*1.D15                                         COLF-099
      FG(J,4,I)=FG(J,4,I)*1.D15                                         COLF-100
    3 CONTINUE                                                          COLF-101
      IF (LMX1.NE.LMAX3) WRITE (MW,1002) LMX1,LMAX3                     COLF-102
      IF (.NOT.LO(108)) IPI(10,I)=MIN0(LMAX3-1,IPI(10,1))               COLF-103
C COMPUTATION OF COULOMB INTEGRALS.                                     COLF-104
      XI=XG(1,I)                                                        COLF-105
      SI=WV(11,I)                                                       COLF-106
      DO 4 J=1,I                                                        COLF-107
      K=NIV(I,J,3)                                                      COLF-108
      IF (K.EQ.0) GO TO 4                                               COLF-109
      XJ=XG(1,J)                                                        COLF-110
      EJ=WV(5,J)                                                        COLF-111
      SJ=WV(11,J)                                                       COLF-112
      CALL CORI(ETA,EJ,RM*SI,RM*SJ,XG(1,K),XI,XJ,FG(1,1,I),FG(1,1,J),Z,LCOLF-113
     1MAX1,LMAX2,LM,LMAX3,FG(1,1,K))                                    COLF-114
    4 CONTINUE                                                          COLF-115
      GO TO 6                                                           COLF-116
    5 LY=0                                                              COLF-117
      IF (WV(3,I).GT.0.D0) LY=MIN0(IDINT(4.D0+3.3D0*DSQRT(WV(3,I))),LMAXCOLF-118
     11)                                                                COLF-119
      IPI(10,I)=LY-1                                                    COLF-120
      IPI(8,I)=KXT                                                      COLF-121
      KXT=KXT+LY*IPI(2,I)                                               COLF-122
    6 CONTINUE                                                          COLF-123
      DO 8 I=1,NCOLL                                                    COLF-124
      XG(1,I)=0.D0                                                      COLF-125
      DO 7 J=2,LMAX2                                                    COLF-126
    7 XG(J,I)=XG(J-1,I)+DATAN2(WV(5,I),DFLOAT(J-1))                     COLF-127
    8 CONTINUE                                                          COLF-128
      LO(133)=LO(133).AND.(.NOT.LO(109))                                COLF-129
      RETURN                                                            COLF-130
 1000 FORMAT (/' LEVEL       WAVE NUMBER  COULOMB PARAMETER   REDUCED MACOLF-131
     1SS       REDUCED ENERGY     STEP SIZE')                           COLF-132
 1001 FORMAT (1X,I5,5F18.10)                                            COLF-133
 1002 FORMAT (' NUMBER OF FINITE COULOMB INTEGRALS REDUCED FROM',I6,' TOCOLF-134
     1',I6)                                                             COLF-135
      END                                                               COLF-136
