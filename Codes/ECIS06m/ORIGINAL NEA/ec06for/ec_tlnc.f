C 28/02/07                                                      ECIS06  TLNC-000
      SUBROUTINE TLNC(HH,IPI,WV,TL,ISM,X,FF,LM,NN,DD,IL,V,VCO,LO)       TLNC-001
C TRANSMISSION COEFFICIENTS OF UNCOUPLED STATES FOR COMPOUND NUCLEUS.   TLNC-002
C INPUT:     HH:      STEP SIZE FOR THE GROUND STATE.                   TLNC-003
C            IPI,WV:  IPI AND WV FOR THIS STATE (SEE CALX).             TLNC-004
C            ISM:     NUMBER OF STEPS.                                  TLNC-005
C            LM:      NUMBER OF COULOMB FUNCTIONS NEEDED.               TLNC-006
C            IL:      LENGTH OF WORKING SPACE.                          TLNC-007
C            V:       POTENTIALS REAL AND IMAGINARY.                    TLNC-008
C            VCO:     STRENGTH OF LONG RANGE COULOMB CORRECTION.        TLNC-009
C            LO(I):   LOGICAL CONTROLS:                                 TLNC-010
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   TLNC-011
C               LO(26) =.TRUE. INTEGRATION STABILISED FOR LONG RANGE    TLNC-012
C                              CONSTANT POTENTIAL.                      TLNC-013
C               LO(27) =.TRUE. NUMEROV'S METHOD FOR SINGLE EQUATIONS.   TLNC-014
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     TLNC-015
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. TLNC-016
C               LO(129)=.TRUE. REAL SPIN-ORBIT OR DIRAC EQUATION.       TLNC-017
C               LO(130)=.TRUE. IMAGINARY SPIN-ORBIT OR DIRAC EQUATION.  TLNC-018
C               LO(133)=.TRUE. STORE SCALAR AND COULOMB POTENTIAL       TLNC-019
C                              INDEPENDENTLY.                           TLNC-020
C OUTPUT:    TL:      TRANSMISSION COEFFICIENTS OF UNCOUPLED LEVELS.    TLNC-021
C WORKING AREAS:                                                        TLNC-022
C            X:       FOR THE INTEGRATION.                              TLNC-023
C            FF:      FOR COULOMB FUNCTIONS AND CORRECTIONS.            TLNC-024
C            NN:      NEEDED BY SUBROUTINE FCOU.                        TLNC-025
C            DD:      NEEDED BY SUBROUTINE CORI.                        TLNC-026
C***********************************************************************TLNC-027
      IMPLICIT REAL*8 (A-H,O-Z)                                         TLNC-028
      LOGICAL LO(150)                                                   TLNC-029
      DIMENSION IPI(11),WV(22),TL(*),X(2,*),FF(LM,*),NN(*),DD(*),V(ISM,*TLNC-030
     1),VCO(2),FAM(4),B(4),G(4),AV(5)                                   TLNC-031
      IF (WV(3).LE.0.D0) RETURN                                         TLNC-032
      J=0                                                               TLNC-033
      CALL FCOU(IPI(10),WV(5),ISM*WV(8)*WV(11),FF,FF(1,2),FF(1,3),FF(1,4TLNC-034
     1),NN,FF(1,5))                                                     TLNC-035
      JC=0                                                              TLNC-036
      IF (.NOT.(LO(103).OR.LO(44))) GO TO 1                             TLNC-037
      JC=5                                                              TLNC-038
      RM=ISM*WV(8)*WV(11)                                               TLNC-039
      CALL CORI(WV(5),WV(5),RM,RM,DD,FF(1,5),FF(1,5),FF,FF,DD(200),LM,LMTLNC-040
     1,IL,LM,FF(1,6))                                                   TLNC-041
    1 L=IPI(10)+1                                                       TLNC-042
      IF (L.LE.0) RETURN                                                TLNC-043
      IJ=IPI(2)                                                         TLNC-044
      VPR=WV(9)**2                                                      TLNC-045
      VPC=WV(10)**2                                                     TLNC-046
      DO 24 LJ=1,L                                                      TLNC-047
      DO 23 JL=1,IJ                                                     TLNC-048
      JJ=2*(LJ+JL)-IJ-3                                                 TLNC-049
      J=J+1                                                             TLNC-050
      TL(J)=0.D0                                                        TLNC-051
      IF (JJ.LT.IJ+1-2*LJ) GO TO 23                                     TLNC-052
      L1=LJ-1                                                           TLNC-053
      CLL=DFLOAT(LJ*(LJ-1))                                             TLNC-054
      CLS=.25D0*DFLOAT(JJ*(JJ+2)-IJ*IJ+1)-CLL                           TLNC-055
C VALUES OF LONG RANGE TAILS OF CENTRAL POTENTIALS.                     TLNC-056
      F2=VCO(1)**2                                                      TLNC-057
      F3=VCO(2)*CLS                                                     TLNC-058
      IF (WV(3).GT.0.D0.AND.JC.GT.0) GO TO 2                            TLNC-059
      F2=0.D0                                                           TLNC-060
      F3=0.D0                                                           TLNC-061
C INTEGRATION REGION - SET UP OF POTENTIAL IN FIVE POINTS FOR           TLNC-062
C TRANSFORMATION OF MATCHING VALUES.                                    TLNC-063
    2 B1=HH*HH/48.D0                                                    TLNC-064
      C1=DFLOAT(ISM-1)*HH                                               TLNC-065
      A1=WV(11)**2                                                      TLNC-066
      IF (WV(3).LT.0.D0) A1=-A1                                         TLNC-067
      DO 3 K=1,5                                                        TLNC-068
      AV(K)=B1*(2.D0*WV(11)*WV(5)/C1-A1+(CLL-F2-F3/C1)/C1**2)           TLNC-069
    3 C1=C1+0.5D0*HH                                                    TLNC-070
C COMPUTATION OF COULOMB CORRECTIONS.                                   TLNC-071
      DO 4 K=1,4                                                        TLNC-072
    4 B(K)=FF(L1+1,K)                                                   TLNC-073
      IF (JC.LE.0) GO TO 8                                              TLNC-074
      IF (WV(5).NE.0.D0) GO TO 5                                        TLNC-075
      F2=F3*WV(11)                                                      TLNC-076
      F3=0.D0                                                           TLNC-077
    5 DO 6 K=1,4                                                        TLNC-078
    6 G(K)=-FF(L1+1,K+JC)*F2                                            TLNC-079
      IF (F3.EQ.0.D0) GO TO 7                                           TLNC-080
      B1=2.D0*WV(5)*DFLOAT(L1*(L1+1))                                   TLNC-081
      B2=DFLOAT(L1+1)**2+WV(5)**2                                       TLNC-082
      C1=-(DFLOAT(2*L1+1)*B2+2*WV(5)**2)/B1                             TLNC-083
      C2=DFLOAT(2*L1+3)*B2/B1                                           TLNC-084
      A1=DFLOAT(ISM)*WV(11)*WV(8)                                       TLNC-085
      D1=(B2+DFLOAT(L1+1)*WV(5)/A1)/A1/B1                               TLNC-086
      D2=-WV(5)*DSQRT(B2)/B1/A1                                         TLNC-087
      A1=B2/B1/A1                                                       TLNC-088
      A3=F3*WV(11)                                                      TLNC-089
      G(1)=G(1)-A3*(C1*FF(L1+1,1+JC)+C2*FF(L1+2,1+JC)-D1*FF(L1+1,1)**2-DTLNC-090
     12*2.D0*FF(L1+1,1)*FF(L1+2,1)-A1*FF(L1+2,1)**2)                    TLNC-091
      G(2)=G(2)-A3*(C1*FF(L1+1,2+JC)+C2*FF(L1+2,2+JC)-D1*FF(L1+1,1)*FF(LTLNC-092
     11+1,3)-D2*(FF(L1+1,1)*FF(L1+2,3)+FF(L1+2,1)*FF(L1+1,3))-A1*FF(L1+2TLNC-093
     2,1)*FF(L1+2,3))                                                   TLNC-094
      G(4)=G(4)-A3*(C1*FF(L1+1,4+JC)+C2*FF(L1+2,4+JC)-D1*FF(L1+1,3)**2-DTLNC-095
     12*2.D0*FF(L1+1,3)*FF(L1+2,3)-A1*FF(L1+2,3)**2)                    TLNC-096
    7 A4=1.D0+(G(1)*G(4)-G(2)**2)                                       TLNC-097
      G(3)=B(1)                                                         TLNC-098
      B(1)=(B(1)*(1.D0-G(2))+G(1)*B(3))/A4                              TLNC-099
      B(3)=(-G(3)*G(4)+(1.D0+G(2))*B(3))/A4                             TLNC-100
      G(3)=B(2)                                                         TLNC-101
      B(2)=(B(2)*(1.D0-G(2))+G(1)*B(4))/A4                              TLNC-102
      B(4)=(-G(3)*G(4)+(1.D0+G(2))*B(4))/A4                             TLNC-103
    8 A1=(1.D0-AV(2))/(2.D0+10.D0*AV(2))                                TLNC-104
      B1=(1.D0-AV(4))/(2.D0+10.D0*AV(4))                                TLNC-105
      A2=A1*(1.D0-AV(1))/(1.D0-4.D0*AV(1))                              TLNC-106
      B2=B1*(1.D0-AV(5))/(1.D0-4.D0*AV(5))                              TLNC-107
      C1=(2.D0+10.D0*AV(3))-(1.D0-AV(3))*(A1+B1)                        TLNC-108
      A1=(16.D0-144.D0*AV(2))/(2.D0+10.D0*AV(2))                        TLNC-109
      B1=(16.D0-144.D0*AV(4))/(2.D0+10.D0*AV(4))                        TLNC-110
      C2=(7.D0+A1*(1.D0-AV(1)))/(1.D0-4.D0*AV(1))                       TLNC-111
      D2=(7.D0+B1*(1.D0-AV(5)))/(1.D0-4.D0*AV(5))                       TLNC-112
      D1=(B1-A1)*(1.D0-AV(3))                                           TLNC-113
      A1=A2*D2+B2*C2                                                    TLNC-114
      B1=(C1*D2+D1*B2)/A1                                               TLNC-115
      B2=30.D0*HH*B2*WV(11)/A1                                          TLNC-116
      FAM(1)=B1*B(1)-B2*B(2)                                            TLNC-117
      FAM(3)=B1*B(3)-B2*B(4)                                            TLNC-118
      B1=(C2*C1-A2*D1)/A1                                               TLNC-119
      A2=-30.D0*HH*A2*WV(11)/A1                                         TLNC-120
      FAM(2)=B1*B(1)-A2*B(2)                                            TLNC-121
      FAM(4)=B1*B(3)-A2*B(4)                                            TLNC-122
      CSO=2.D0*VPR*CLS                                                  TLNC-123
C COMPUTATION OF THE REGULAR SOLUTION.                                  TLNC-124
      DO 9 IS=1,ISM                                                     TLNC-125
      X(1,IS+2)=WV(12)-CLL/DFLOAT(IS)**2+VPR*V(IS,1)                    TLNC-126
    9 X(2,IS+2)=VPR*V(IS,2)                                             TLNC-127
      IF ((.NOT.LO(101)).OR.(CSO.EQ.0.D0)) GO TO 12                     TLNC-128
      DO 10 IS=1,ISM                                                    TLNC-129
   10 X(1,IS+2)=X(1,IS+2)+CSO*V(IS,3)                                   TLNC-130
      IF (.NOT.LO(102)) GO TO 12                                        TLNC-131
      DO 11 IS=1,ISM                                                    TLNC-132
   11 X(2,IS+2)=X(2,IS+2)+CSO*V(IS,4)                                   TLNC-133
   12 IF (.NOT.LO(133)) GO TO 16                                        TLNC-134
      IF (WV(5).EQ.0.D0) GO TO 14                                       TLNC-135
      DO 13 IS=1,ISM                                                    TLNC-136
   13 X(1,IS+2)=X(1,IS+2)+VPC*V(IS,5)                                   TLNC-137
   14 IF (.NOT.LO(103)) GO TO 16                                        TLNC-138
      VPD=2.D0*CLS*VPC                                                  TLNC-139
      IF (VPD*VCO(2).EQ.0.D0) GO TO 16                                  TLNC-140
      DO 15 IS=1,ISM                                                    TLNC-141
   15 X(1,IS+2)=X(1,IS+2)+VPD*V(IS,6)                                   TLNC-142
   16 CONTINUE                                                          TLNC-143
      IF (LO(27)) GO TO 18                                              TLNC-144
C MODIFIED NUMEROV METHOD.                                              TLNC-145
      DO 17 IS=1,ISM                                                    TLNC-146
      A1=X(1,IS+2)**2-X(2,IS+2)**2                                      TLNC-147
      IF (LO(26)) A1=A1-X(1,IS+2)**3/30.D0                              TLNC-148
      X(2,IS+2)=X(2,IS+2)*(1.D0-X(1,IS+2)/6.D0)                         TLNC-149
   17 X(1,IS+2)=X(1,IS+2)-A1/12.D0                                      TLNC-150
      GO TO 20                                                          TLNC-151
C NUMEROV METHOD.                                                       TLNC-152
   18 DO 19 IS=1,ISM                                                    TLNC-153
      A2=(12.D0+X(1,IS+2))**2+X(2,IS+2)**2                              TLNC-154
      A1=12.D0*(X(1,IS+2)*(12.D0+X(1,IS+2))+X(2,IS+2)**2)/A2            TLNC-155
      IF (LO(26)) A1=A1+X(1,IS+2)**3/240.D0                             TLNC-156
      X(1,IS+2)=A1                                                      TLNC-157
   19 X(2,IS+2)=144.D0*X(2,IS+2)/A2                                     TLNC-158
   20 X(1,1)=0.D0                                                       TLNC-159
      X(2,1)=0.D0                                                       TLNC-160
      X(1,2)=1.D-15                                                     TLNC-161
      X(2,2)=0.D0                                                       TLNC-162
      DO 22 IS=1,ISM                                                    TLNC-163
      B1=X(1,IS+1)*X(1,IS+2)-X(2,IS+1)*X(2,IS+2)                        TLNC-164
      B2=X(2,IS+1)*X(1,IS+2)+X(1,IS+1)*X(2,IS+2)                        TLNC-165
      X(1,IS+2)=X(1,IS+1)+X(1,IS+1)-X(1,IS)-B1                          TLNC-166
      X(2,IS+2)=X(2,IS+1)+X(2,IS+1)-X(2,IS)-B2                          TLNC-167
      IF (DABS(X(1,IS+2)).LT.1.D15) GO TO 22                            TLNC-168
C RENORMALISATION OF LARGE FUNCTION VALUES.                             TLNC-169
      JS=2*IS+4                                                         TLNC-170
      DO 21 I=3,JS                                                      TLNC-171
   21 X(I,1)=X(I,1)*1.D-30                                              TLNC-172
   22 CONTINUE                                                          TLNC-173
C END OF INTEGRATION.                                                   TLNC-174
C MATCHING.                                                             TLNC-175
      A1=X(1,ISM)*FAM(4)-FAM(3)*X(1,ISM+2)                              TLNC-176
      A2=X(2,ISM)*FAM(4)-FAM(3)*X(2,ISM+2)                              TLNC-177
      B1=X(1,ISM)*FAM(2)-FAM(1)*X(1,ISM+2)                              TLNC-178
      B2=X(2,ISM)*FAM(2)-FAM(1)*X(2,ISM+2)                              TLNC-179
      A2=A2+B1                                                          TLNC-180
      A1=A1-B2                                                          TLNC-181
      A3=A1*A1+A2*A2                                                    TLNC-182
      A1=-A1/A3                                                         TLNC-183
      A2=A2/A3                                                          TLNC-184
      CR=A1*B1-A2*B2                                                    TLNC-185
      CI=B1*A2+B2*A1                                                    TLNC-186
      TL(J)=DMAX1(0.D0,4.D0*(CI-CR**2-CI**2))                           TLNC-187
   23 CONTINUE                                                          TLNC-188
   24 CONTINUE                                                          TLNC-189
      RETURN                                                            TLNC-190
      END                                                               TLNC-191
