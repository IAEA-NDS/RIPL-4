C 29/10/06                                                      ECIS06  ROTZ-000
      SUBROUTINE ROTZ(VR,XGN,WV,R,Q,VAL,J,IV,B,P,ID1,IQ1,LO)            ROTZ-001
C DEFORMATION OF THE SCHROEDINGER EQUATION WHICH IS EQUIVALENT TO THE   ROTZ-002
C DIRAC EQUATION ONLY FOR ELASTIC SCATTERING:                           ROTZ-003
C     NO ASYMMETRIC ROTATIONAL MODEL,                                   ROTZ-004
C     NO DIFFUSENESS OF COULOMB CHARGE AND NO COULOMB DEFORMATION,      ROTZ-005
C     NO VIBRATIONAL BAND IN ROTATIONAL MODEL,                          ROTZ-006
C     ONLY FIRST ORDER VIBRATIONAL MODEL.                               ROTZ-007
C INPUT:     VR(I,J): POTENTIAL FOR I=1, -1RST DERIVATIVE FOR I=5       ROTZ-008
C                     AND SECOND DERIVATIVE FOR I=6.                    ROTZ-009
C            XGN:     ABSCISSAE FOR LEGENDRE INTEGRATION.               ROTZ-010
C            WV:      MASS IN WV(6), ENERGY IN WV(7).                   ROTZ-011
C            R:       RADIUS.                                           ROTZ-012
C            Q(I,J,K):DERIVATIVES DR/DTHETA FOR K=2 AND 3.              ROTZ-013
C            VAL:     OPTICAL MODEL.                                    ROTZ-014
C            J:       POINT OF LEGENDRE ANGULAR INTEGRAL.               ROTZ-015
C            IV:      NUMBER OF PHONONS PLUS ONE.                       ROTZ-016
C            B:       DEFORMATIONS.                                     ROTZ-017
C            P(I,J):  PREVIOUS RESULTS FOR I=17 TO 26 AND WEIGHTS       ROTZ-018
C                     FOR I=28 AND UP IN THE ROTATIONAL MODEL.          ROTZ-019
C            ID1:     FIRST DIMENSION OF P.                             ROTZ-020
C            IQ1:     NUMBER OF TRANSITION FORM FACTORS PLUS 3.         ROTZ-021
C            LO(I):   LOGICAL CONTROLS:                                 ROTZ-022
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).ROTZ-023
C OUTPUT:    P(I,J)   FOR I=17 TO 26, POTENTIAL FOR J=1 AND             ROTZ-024
C                     TRANSITION FORM FACTORS FOR J=4 TO IQ1.           ROTZ-025
C                                                                       ROTZ-026
C FOR THE COMMON  /DCONS/ SEE CALC.                                     ROTZ-027
C                                                                       ROTZ-028
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     ROTZ-029
C  CCZ:       COULOMB ALPHA CONSTANT.                                   ROTZ-030
C   USED:     CCZ.                                                      ROTZ-031
C                                                                       ROTZ-032
C***********************************************************************ROTZ-033
      IMPLICIT REAL*8 (A-H,O-Z)                                         ROTZ-034
      LOGICAL LO(150)                                                   ROTZ-035
      DIMENSION VR(6,10),XGN(10),WV(22),Q(8,12,*),VAL(32),IV(*),B(10,*),ROTZ-036
     1P(ID1,*)                                                          ROTZ-037
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            ROTZ-038
      COMMON /INOUT/ MR,MW,MS                                           ROTZ-039
      XP=CHB**2/(2.D0*WV(6))                                            ROTZ-040
      DO 2 I=7,8                                                        ROTZ-041
      IF (R.GT.VAL(4*I-2)) GO TO 1                                      ROTZ-042
      VR(1,I)=CCZ*VAL(4*I-3)*(1.5D0-.5D0*(R*R)/VAL(4*I-2)**2)/VAL(4*I-2)ROTZ-043
      VR(6,I)=-CCZ*VAL(4*I-3)/VAL(4*I-2)**3                             ROTZ-044
      VR(5,I)=-VR(6,I)*R                                                ROTZ-045
      GO TO 2                                                           ROTZ-046
    1 VR(1,I)=CCZ*VAL(3*I-2)/R                                          ROTZ-047
      VR(5,I)=VR(1,I)/R                                                 ROTZ-048
      VR(6,I)=2.D0*VR(5,I)/R                                            ROTZ-049
    2 CONTINUE                                                          ROTZ-050
C COMPUTATION OF D(R) AND INVERSE                                       ROTZ-051
      C1=WV(7)+WV(6)+VR(1,1)-VR(1,3)-VR(1,7)                            ROTZ-052
      C2=VR(1,2)-VR(1,4)                                                ROTZ-053
      DD=C1**2+C2**2                                                    ROTZ-054
      DR=C1/DD                                                          ROTZ-055
      DI=-C2/DD                                                         ROTZ-056
      C3=WV(7)-WV(6)-VR(1,1)-VR(1,3)-VR(1,7)                            ROTZ-057
      C4=-VR(1,2)-VR(1,4)                                               ROTZ-058
      C5=(C1*C3-C2*C4-WV(7)**2+WV(6)**2)/(2.D0*WV(6))                   ROTZ-059
      C6=(C1*C4+C2*C3)/(2.D0*WV(6))                                     ROTZ-060
      D1=VR(5,3)-VR(5,1)+VR(5,7)                                        ROTZ-061
      D2=VR(5,4)-VR(5,2)                                                ROTZ-062
      CR=D1*DR-D2*DI+(VR(5,5)+VR(5,8))/WV(6)                            ROTZ-063
      CI=D2*DR+D1*DI+VR(5,6)/WV(6)                                      ROTZ-064
      CC=VR(1,7)*WV(7)/WV(6)                                            ROTZ-065
      BR=VR(6,1)-VR(6,7)-VR(6,3)                                        ROTZ-066
      BI=VR(6,2)-VR(6,4)                                                ROTZ-067
C SQUARE OF GRADIENT  AND LAPLACIAN     RADIAL TERM.                    ROTZ-068
      ER=BR*DR-BI*DI+(VR(5,5)+VR(5,8))/WV(6)*(CR+CR-(VR(5,5)+VR(5,8))/WVROTZ-069
     1(8))-VR(5,6)/WV(6)*(CI+CI-VR(5,6)/WV(6))-(VR(6,5)+VR(6,8))/WV(6)+2ROTZ-070
     2.D0*CR/R                                                          ROTZ-071
      EI=BR*DI+BI*DR+(VR(5,5)+VR(5,8))/WV(6)*(CI+CI-VR(5,6)/WV(6))+VR(5,ROTZ-072
     16)/WV(6)*(CR+CR-(VR(5,5)+VR(5,8))/WV(6))-VR(6,6)/WV(6)+2.D0*CI/R  ROTZ-073
      C5=C5-(.75D0*(CR**2-CI**2)-.5D0*ER)*XP                            ROTZ-074
      C6=C6-(1.5D0*CR*CI-.5D0*EI)*XP                                    ROTZ-075
      IF (.NOT.LO(1)) GO TO 5                                           ROTZ-076
C ROTATIONAL MODEL.                                                     ROTZ-077
      C=1.D0-XGN(J)**2                                                  ROTZ-078
      DO 3 I=1,8                                                        ROTZ-079
      VR(4,I)=(VR(6,I)*Q(I,J,2)**2*C+VR(5,I)*Q(I,J,3))/R**2             ROTZ-080
    3 VR(3,I)=VR(5,I)*Q(I,J,2)/R                                        ROTZ-081
C SQUARE OF GRADIENT IS D/DR.D/DR + 1/R**2 D/DTH.D/DTH.                 ROTZ-082
      C3=VR(3,1)-VR(3,3)-VR(3,7)                                        ROTZ-083
      C4=VR(3,2)-VR(3,4)                                                ROTZ-084
      C7=C3*DR-C4*DI                                                    ROTZ-085
      C8=C3*DI+C4*DR                                                    ROTZ-086
      C3=(VR(3,5)+VR(3,8))/WV(6)                                        ROTZ-087
      C4=VR(3,6)/WV(6)                                                  ROTZ-088
      C5=C5-(.75D0*(C7**2-C8**2)+.5D0*(C7*C3-C8*C4)+.75D0*(C3**2-C4**2))ROTZ-089
     1*C*XP                                                             ROTZ-090
      C6=C6-(1.5D0*C7*C8+.5D0*(C7*C4+C8*C3)+1.5D0*C3*C4)*C*XP           ROTZ-091
C LAPLACIAN IS D**2/DR**2 + 2/R D/DR + 1/R**2 (D**2/DTH**2+D/DTH/TG ).  ROTZ-092
      C3=VR(4,1)-VR(4,3)-VR(4,7)                                        ROTZ-093
      C4=VR(4,2)-VR(4,4)                                                ROTZ-094
      C5=C5+.5D0*(C3*DR-C4*DI+(VR(4,5)+VR(4,8))/WV(6))*XP               ROTZ-095
      C6=C6+.5D0*(C4*DR+C3*DI+VR(4,6)/WV(6))*XP                         ROTZ-096
C C5 + I C6   CENTRAL TERM,                                             ROTZ-097
C CR + I CI   SPIN-ORBIT TERM,                                          ROTZ-098
C DR + I DI   1/DENOMINATOR.                                            ROTZ-099
      P(17,1)=P(17,1)+(C5+CC)*P(28,J)                                   ROTZ-100
      P(18,1)=P(18,1)+C6*P(28,J)                                        ROTZ-101
      P(21,1)=P(21,1)+.5D0*XP*(CR-VR(5,8)/WV(6))*P(28,J)/R              ROTZ-102
      P(22,1)=P(22,1)+.5D0*XP*CI*P(28,J)/R                              ROTZ-103
      P(23,1)=P(23,1)+CC*P(28,J)                                        ROTZ-104
      P(24,1)=P(24,1)+VR(5,8)/R*P(28,J)                                 ROTZ-105
      IF (IQ1.LT.4) RETURN                                              ROTZ-106
C LOGARITHM OF D(R) C9 + I C0.                                          ROTZ-107
      C9=.5D0*DLOG(C1**2+C2**2)-(VR(1,5)+VR(1,8))/WV(6)                 ROTZ-108
      C0=DATAN2(C2,C1)-VR(1,6)/WV(6)                                    ROTZ-109
      DO 4 I=4,IQ1                                                      ROTZ-110
      IF (IV(I).GT.1) GO TO 7                                           ROTZ-111
      P(17,I)=P(17,I)+(C5+CC)*P(I+27,J)                                 ROTZ-112
      P(18,I)=P(18,I)+C6*P(I+27,J)                                      ROTZ-113
      P(21,I)=P(21,I)+.5D0*XP*(CR-VR(5,8)/WV(6))*P(I+27,J)/R            ROTZ-114
      P(22,I)=P(22,I)+.5D0*XP*CI*P(I+27,J)/R                            ROTZ-115
      P(23,I)=P(23,1)+CC*P(I+27,J)                                      ROTZ-116
      P(24,I)=P(24,1)+VR(5,8)/R*P(I+27,J)                               ROTZ-117
      P(25,I)=P(25,I)-.5D0*XP*C9*P(I+27,J)/R**2                         ROTZ-118
    4 P(26,I)=P(26,I)-.5D0*XP*C0*P(I+27,J)/R**2                         ROTZ-119
      RETURN                                                            ROTZ-120
C COMPUTATION OF ZEROTH ORDER.                                          ROTZ-121
    5 P(17,1)=C5+CC                                                     ROTZ-122
      P(18,1)=C6                                                        ROTZ-123
      P(21,1)=.5D0*(CR-VR(5,8)/WV(6))*XP/R                              ROTZ-124
      P(22,1)=.5D0*CI*XP/R                                              ROTZ-125
      P(23,1)=CC                                                        ROTZ-126
      P(24,1)=VR(5,8)/R                                                 ROTZ-127
      IF (IQ1.LT.4) RETURN                                              ROTZ-128
      DO 6 I=4,IQ1                                                      ROTZ-129
C COMPUTATION OF FIRST ORDER.                                           ROTZ-130
      IF (IV(I).NE.2) GO TO 8                                           ROTZ-131
      C1R=VR(5,1)*B(1,I)-VR(5,3)*B(3,I)                                 ROTZ-132
      C2R=VR(5,2)*B(2,I)-VR(5,4)*B(4,I)                                 ROTZ-133
      DRR=-(DR*DR-DI*DI)*C1R+2.D0*DR*DI*C2R                             ROTZ-134
      DIR=-(DR*DR-DI*DI)*C2R-2.D0*DR*DI*C1R                             ROTZ-135
      C3R=-VR(5,1)*B(1,I)-VR(5,3)*B(3,I)                                ROTZ-136
      C4R=-VR(5,2)*B(2,I)-VR(5,4)*B(4,I)                                ROTZ-137
      C5R=(C1R*C3-C2R*C4+C1*C3R-C2*C4R)/(2.D0*WV(6))                    ROTZ-138
      C6R=(C1R*C4+C2R*C3+C1*C4R+C2*C3R)/(2.D0*WV(6))                    ROTZ-139
      D1R=VR(6,3)*B(3,I)-VR(6,1)*B(1,I)                                 ROTZ-140
      D2R=VR(6,4)*B(4,I)-VR(6,2)*B(2,I)                                 ROTZ-141
      CRR=D1R*DR-D2R*DI+D1*DRR-D2*DIR+VR(6,5)*B(5,I)/WV(6)              ROTZ-142
      CIR=D2R*DR+D1R*DI+D2*DRR+D1*DIR+VR(6,6)*B(6,I)/WV(6)              ROTZ-143
      AR=BR*DRR-BI*DIR+(VR(5,5)+VR(5,8))/WV(6)*(CRR+CRR-VR(6,5)*B(5,I)/WROTZ-144
     1V(8))-VR(5,6)/WV(6)*(CIR+CIR-VR(6,6)*B(6,I)/WV(6))+VR(6,5)*B(5,I)/ROTZ-145
     2WV(6)*(CR+CR-(VR(5,5)+VR(5,8))/WV(6))-VR(6,6)*B(6,I)/WV(6)*(CI+CI-ROTZ-146
     3VR(5,6)/WV(6))+2.D0*CRR/R                                         ROTZ-147
      AI=BR*DIR+BI*DRR+(VR(5,5)+VR(5,8))/WV(6)*(CIR+CIR-VR(6,6)*B(6,I)/WROTZ-148
     1V(8))+VR(5,6)/WV(6)*(CRR+CRR-VR(6,5)*B(5,I)/WV(6))+VR(6,5)*B(5,I)/ROTZ-149
     2WV(6)*(CI+CI-VR(5,6)/WV(6))+VR(6,6)*B(6,I)/WV(6)*(CR+CR-VR(6,5)/WVROTZ-150
     3(8))+2.D0*CIR/R                                                   ROTZ-151
      P(17,I)=C5R-(1.5D0*(CR*CRR-CI*CIR)-.5D0*AR)*XP                    ROTZ-152
      P(18,I)=C6R-(1.5D0*(CRR*CI+CR*CIR)-.5D0*AI)*XP                    ROTZ-153
      P(21,I)=.5D0*CRR*XP/R                                             ROTZ-154
      P(22,I)=.5D0*CIR*XP/R                                             ROTZ-155
      P(23,I)=0.D0                                                      ROTZ-156
      P(24,I)=0.D0                                                      ROTZ-157
      P(25,I)=-.5D0*XP*(C1R*DR-C2R*DI-VR(5,5)*B(5,I)/WV(6))/R**2        ROTZ-158
    6 P(26,I)=-.5D0*XP*(C1R*DI+C2R*DR-VR(5,6)*B(6,I)/WV(6))/R**2        ROTZ-159
      RETURN                                                            ROTZ-160
    7 WRITE (MW,1000)                                                   ROTZ-161
      GO TO 9                                                           ROTZ-162
    8 WRITE (MW,1001)                                                   ROTZ-163
    9 WRITE (MW,1002)                                                   ROTZ-164
      STOP                                                              ROTZ-165
 1000 FORMAT (' NO VIBRATIONAL BAND IN ROTATIONAL MODEL.')              ROTZ-166
 1001 FORMAT (' ONLY FIRST ORDER VIBRATIONAL MODEL.')                   ROTZ-167
 1002 FORMAT (//' IN ROTZ  ...  STOP  ...')                             ROTZ-168
      END                                                               ROTZ-169
