C 04/07/06                                                      ECIS06  VIBM-000
      SUBROUTINE VIBM(NIV,IQ,T,IPI,NCOLL,IT,IPH,NVAR,VAR,NBETA,FAC,IDT,LVIBM-001
     1O)                                                                VIBM-002
C NUCLEAR REDUCED MATRIX ELEMENTS FOR THE HARMONIC VIBRATIONAL MODEL.   VIBM-003
C INPUT:     IPI(J,*):PARITY OF NUCLEAR STATES FOR J=1,                 VIBM-004
C                     MULTIPLICITY OF THE PARTICLE FOR J=2,             VIBM-005
C                     MULTIPLICITY OF THE TARGET FOR J=3.               VIBM-006
C            NCOLL:   NUMBER OF COUPLED STATES                          VIBM-007
C            IPH(I,J):DESCRIPTION OF VIBRATIONAL MODEL FOR LEVEL J:     VIBM-008
C                     NUMBER OF PHONONS (0, 1, 2 PHONONS OR 3 FOR       VIBM-009
C                     MIXTURE OF 1 AND 2 PHONONS STATES) FOR I=1.       VIBM-010
C                     ADDRESS OF THE DESCRIPTION OF TWO PHONONS AND     VIBM-011
C                     MIXED STATES WHICH ARE IN THE NVAR,VAR FOR I=2.   VIBM-012
C            NVAR,VAR:EQUIVALENT BY CALL, 1 AND 2 PHONONS MIXING        VIBM-013
C                     COEFFICIENTS.                                     VIBM-014
C            NBETA:   QUANTUM NUMBERS OF DEFORMATIONS IN NBETA(17,*),   VIBM-015
C                     0 IN NBETA(18,*) TO BE USED IN SECOND ORDER       VIBM-016
C                     MONOPOLE CORRECTION, ANYTHING NOT TO BE USED/     VIBM-017
C            NBT1:    NUMBER OF PHONONS.                                VIBM-018
C            FAC:     TABLE AND NUMBER OF LOGARITHMS OF FACTORIALS.     VIBM-019
C            IDT:     LENGTH AVAILABLE IN THIS SUBROUTINE.              VIBM-020
C            LO(I):   LOGICAL CONTROLS:                                 VIBM-021
C               LO(2)  =.TRUE. SECOND ORDER VIBRATIONAL OR CONSTRAINED  VIBM-022
C                              ASYMMETRIC ROTATIONAL MODEL.             VIBM-023
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      VIBM-024
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   VIBM-025
C               LO(117)=.TRUE. FOR ALL CALCULATIONS EXCEPT THE FIRST.   VIBM-026
C OUTPUT:    NIV:     NUMBER AND REFERENCES OF INTERACTION FORM FACTORS VIBM-027
C                     BETWEEN EACH CHANNELS. NIV(I1,I2,K): FIRST I OF   VIBM-028
C                     T(3,I) FOR THE PAIR OF NUCLEAR STATES I1,I2 FOR   VIBM-029
C                     K=1 AND LAST ONE FOR K=2.                         VIBM-030
C            IQ,T:    EQUIVALENT BY CALL OF IQ(6,I) AND T(3,I).         VIBM-031
C            IT:      NUMBER OF NUCLEAR MATRIX ELEMENTS IN IQ,T. FOR A  VIBM-032
C                     GIVEN VALUE OF I:                                 VIBM-033
C                     IQ(1,I): REFERENCE TO TABLE OF FORM FACTORS:      VIBM-034
C                              I FOR BETA(I) AND I+J*(NBT1+1) FOR       VIBM-035
C                              BETA(I)*BETA(J) WITH J LARGER THAN I IN  VIBM-036
C                              THE HARMONIC VIBRATIONAL MODEL (ORDER OF VIBM-037
C                              DERIVATIVE IN THE ANHARMONIC VIBRATIONAL VIBM-038
C                              MODEL)                                   VIBM-039
C                     IQ(2,I): REFERENCE TO TABLE OF ANGULAR MOMENTA,   VIBM-040
C                     IQ(3,I): ADDRESS OF THE ASSOCIATED SPIN-ORBIT FORMVIBM-041
C                              FACTOR OR 0,                             VIBM-042
C                     IQ(4,I): UNUSED,                                  VIBM-043
C                     T(3,I):  REDUCED NUCLEAR MATRIX ELEMENT MULTIPLIEDVIBM-044
C                              BY (-)**(L/2) WHERE L IS THE TRANSFERRED VIBM-045
C                              ANGULAR MOMENTUM                         VIBM-046
C                                                                       VIBM-047
C***********************************************************************VIBM-048
C  THE DEFORMATIONS BETA AND FACTORS 1/SQRT(4*PI) ARE NOT INCLUDED IN   VIBM-049
C THE MATRIX ELEMENTS WHICH ARE COMPUTED HERE. THE FULL EXPRESSIONS ARE:VIBM-050
C (0||Q2||0) = SUM ON BETA**2/(4*PI) WITH IQ=0,                         VIBM-051
C (0||Q1||I) = (-)**I BETA(I)/SQRT(4*PI) WITH IQ=I,                     VIBM-052
C (I||Q1||0) SAME VALUE WITHOUT (-)**I,                                 VIBM-053
C (IP||Q2||I) = (-)**I BETA(I)*BETA(IP)*DJCG(I,IP,0,0|IQ,0)/(2*PI)      VIBM-054
C    PLUS SUM ON BETA**2*SQRT(2*I+1)/(4*PI) WITH IQ=0 WHEN I=IP,        VIBM-055
C (0||Q2||L1,L2,I) = (-)**I BETA(L1)*BETA(L2) DJCG(L1,L2,0,0|I,0)/      VIBM-056
C    (2*PI*SQRT(1+DELTA(L1,L2)))   WITH IQ=I,                           VIBM-057
C (L1,L2,I||Q2||0) SAME VALUE WITHOUT (-)**I,                           VIBM-058
C (IP||Q1||L1,L2,I) = (-)**IQ BETA(*)*(DELTA(IQ,L1)*DELTA(IP,L2)+(-)**  VIBM-059
C    (IP+I+IQ)*DELTA(IQ,L2)*DELTA(IP,L1)) * SQRT((2*I+1)/((2*IQ+1)*     VIBM-060
C    (1+DELTA(L1,L2))*SQRT(4*PI)),                                      VIBM-061
C (L1,L2,I||Q1||IP) SAME EXPRESSION BUT WITH (-)**IP+I+IQ  IN FRONT,    VIBM-062
C (L3,L4,IP||Q2||L1,L2,I)  SUM ON BETA**2*(DELTA(L1,L3)*DELTA(L2,L4)+   VIBM-063
C    (-)**(L1+L2-I)*DELTA(L1,L4)*DELTA(L2,L3)) * SQRT(2*IP+1)/(4*PI*SQRTVIBM-064
C    ((1+DELTA(L1,L2))*(1+DELTA(L3,L4))) WITH IQ=0,WHEN I=IP            VIBM-065
C    PLUS,WHEN TWO PHONONS ARE IDENTICAL,SUM WITH ALL THE POSSIBLE VA-  VIBM-066
C    LUES FOR IQ OF BETA(L5)*BETA(L6)*DJCG(L5,L6,0,0|IQ,0)*DJ6J(L6,IP,  VIBM-067
C    L7,I,L5,IQ) * SQRT((2*I+1)*(2*IP+1))/(2*PI*SQRT((1+DELTA(L1,L2)*(1+VIBM-068
C    DELTA(L3,L4))   WHERE  L7 IS THE COMMON PHONON ,L5 AND L6 THE OTHERVIBM-069
C    PHONONS IN I AND IP, MULTIPLIED BY THE PHASE:                      VIBM-070
C      *(-)**(IP-L1)   IF L1=L3,        *(-)**(L3)       IF L1=L4,      VIBM-071
C      *(-)**(IP+L1-I) IF L2=L3,      *(-)**(L1+L3+L4-I) IF L2=L4.      VIBM-072
C (J||L=IQ||I) IS MULTIPLIED BY A PHASE (-)**((IQ+IPI(I)-IPI(J))/2)     VIBM-073
C AND THE FACTORS BETA/SQRT(4*PI) ARE SHIFTED ON THE FORM FACTORS.      VIBM-074
C                                                                       VIBM-075
C FOR THE COMMON  /COUPL/ SEE CALC.                                     VIBM-076
C                                                                       VIBM-077
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     VIBM-078
C  NBT1:      NUMBER OF PHONONS (VIBRATIONS).                           VIBM-079
C  NFA:       NUMBER OF LOGARITHMS OF FACTORIALS.                       VIBM-080
C  NVA:       NUMBER OF NUCLEAR PARAMETERS.                             VIBM-081
C   USED:     NBT1,NFA,NVA.                                             VIBM-082
C                                                                       VIBM-083
C***********************************************************************VIBM-084
      IMPLICIT REAL*8 (A-H,O-Z)                                         VIBM-085
      LOGICAL LIB(4),LO(150)                                            VIBM-086
      DIMENSION NIV(NCOLL,NCOLL,*),IQ(6,*),T(3,*),IPI(11,*),IPH(2,*),NVAVIBM-087
     1R(2,*),VAR(*),NBETA(18,*),FAC(*),IA(6,4),B(4),AA(2,2)             VIBM-088
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   VIBM-089
      COMMON /INOUT/ MR,MW,MS                                           VIBM-090
      IT=1                                                              VIBM-091
      NSP=0                                                             VIBM-092
      IF (LO(13).OR.LO(19)) NSP=10                                      VIBM-093
      DO 39 I1=1,NCOLL                                                  VIBM-094
      AA(1,1)=1.D0                                                      VIBM-095
      AA(1,2)=0.D0                                                      VIBM-096
      IF (IPH(1,I1).LE.2) GO TO 1                                       VIBM-097
      JVAR=IPH(2,I1)+1                                                  VIBM-098
      IVAR=NVAR(2,JVAR)                                                 VIBM-099
      AY=1.74532925D-02*VAR(IVAR)                                       VIBM-100
      AA(1,1)=DCOS(AY)                                                  VIBM-101
      AA(1,2)=DSIN(AY)                                                  VIBM-102
      IF (.NOT.LO(117)) WRITE (MW,1000) I1,VAR(IVAR),AA(1,1),AA(1,2)    VIBM-103
    1 AA(2,1)=AA(1,1)                                                   VIBM-104
      AA(2,2)=AA(1,2)                                                   VIBM-105
      DO 38 I2=I1,NCOLL                                                 VIBM-106
      IF (I1.EQ.I2) GO TO 2                                             VIBM-107
      AA(2,1)=1.D0                                                      VIBM-108
      AA(2,2)=0.D0                                                      VIBM-109
      IF (IPH(1,I2).LE.2) GO TO 2                                       VIBM-110
      IVAR=IPH(2,I2)+1                                                  VIBM-111
      JVAR=NVAR(2,IVAR)                                                 VIBM-112
      IF (JVAR.GT.NVA) GO TO 40                                         VIBM-113
      AY=1.74532925D-02*VAR(JVAR)                                       VIBM-114
      AA(2,1)=DCOS(AY)                                                  VIBM-115
      AA(2,2)=DSIN(AY)                                                  VIBM-116
    2 NIV(I1,I2,1)=IT                                                   VIBM-117
      L1=IPH(1,I1)                                                      VIBM-118
      AX=AA(1,1)                                                        VIBM-119
      IF (L1.GT.2) L1=1                                                 VIBM-120
      IF (DABS(AX).LT.1.D-6) GO TO 27                                   VIBM-121
    3 L2=IPH(1,I2)                                                      VIBM-122
      IF (L2.GT.2) L2=1                                                 VIBM-123
      AY=AX*AA(2,1)                                                     VIBM-124
      IF (DABS(AY).LT.1.D-6) GO TO 26                                   VIBM-125
    4 I=L1+L2+1                                                         VIBM-126
      IF (I.EQ.3.AND.L1.NE.L2) I=6                                      VIBM-127
      IF (L1.GT.L2) GO TO 5                                             VIBM-128
      J1=I1                                                             VIBM-129
      J2=I2                                                             VIBM-130
      GO TO 6                                                           VIBM-131
    5 J1=I2                                                             VIBM-132
      J2=I1                                                             VIBM-133
C  TRANSPOSITION.                                                       VIBM-134
      AY=AY*DFLOAT(1-MOD(IPI(3,I1)+IPI(3,I2)+2*(IPI(1,I1)+IPI(1,I2)+1),4VIBM-135
     1))                                                                VIBM-136
    6 GO TO ( 7 , 9 , 10 , 14 , 18 , 24 ) , I                           VIBM-137
C  (0||Q||0).                                                           VIBM-138
    7 IF ((IPI(3,J1).NE.1).OR.(IPI(3,J2).NE.1)) GO TO 41                VIBM-139
      IF (.NOT.LO(2)) GO TO 26                                          VIBM-140
      DO 8 L=1,NBT1                                                     VIBM-141
      IF (NBETA(18,L).NE.0) GO TO 8                                     VIBM-142
      IF (2*IT.GT.IDT) CALL MEMO('VIBM',IDT,2*IT)                       VIBM-143
      IQ(1,IT)=L*(NBT1+2)                                               VIBM-144
      IQ(2,IT)=0                                                        VIBM-145
      IQ(3,IT)=NSP                                                      VIBM-146
      T(3,IT)=AY                                                        VIBM-147
      IT=IT+1                                                           VIBM-148
    8 CONTINUE                                                          VIBM-149
      GO TO 26                                                          VIBM-150
C  (IP||Q||0)  WITH IP=J2.                                              VIBM-151
    9 IF (2*IT.GT.IDT) CALL MEMO('VIBM',IDT,2*IT)                       VIBM-152
      N2=IPH(2,J2)                                                      VIBM-153
      IF (IPH(1,J2).GT.2) N2=NVAR(1,N2+1)                               VIBM-154
      IF ((IPI(3,J2).NE.(2*NBETA(17,N2)+1)).OR.(IPI(3,J1).NE.1)) GO     VIBM-155
     1TO 41                                                             VIBM-156
      IQ(1,IT)=N2                                                       VIBM-157
      IQ(2,IT)=NBETA(17,N2)                                             VIBM-158
      IQ(3,IT)=NSP                                                      VIBM-159
      T(3,IT)=AY                                                        VIBM-160
      IF (MOD(IABS(IQ(2,IT)+IPI(1,J1)-IPI(1,J2)),4).NE.0) T(3,IT)=-T(3,IVIBM-161
     1T)                                                                VIBM-162
      GO TO 25                                                          VIBM-163
C  (IP||Q||I).                                                          VIBM-164
   10 IF (.NOT.LO(2)) GO TO 26                                          VIBM-165
      N1=IPH(2,J1)                                                      VIBM-166
      IF (IPH(1,J1).GT.2) N1=NVAR(1,N1+1)                               VIBM-167
      N2=IPH(2,J2)                                                      VIBM-168
      IF (IPH(1,J2).GT.2) N2=NVAR(1,N2+1)                               VIBM-169
      IF ((IPI(3,J2).NE.(2*NBETA(17,N2)+1)).OR.(IPI(3,J1).NE.(2*NBETA(17VIBM-170
     1,N1)+1))) GO TO 41                                                VIBM-171
      IF (N1.NE.N2) GO TO 12                                            VIBM-172
      AQ=DSQRT(DFLOAT(2*NBETA(17,N1)+1))                                VIBM-173
      DO 11 L=1,NBT1                                                    VIBM-174
      IF (NBETA(18,L).NE.0) GO TO 11                                    VIBM-175
      IF (2*IT.GT.IDT) CALL MEMO('VIBM',IDT,2*IT)                       VIBM-176
      IQ(1,IT)=L*(NBT1+2)                                               VIBM-177
      IQ(2,IT)=0                                                        VIBM-178
      IQ(3,IT)=NSP                                                      VIBM-179
      T(3,IT)=AQ*AY                                                     VIBM-180
      IT=IT+1                                                           VIBM-181
   11 CONTINUE                                                          VIBM-182
C FACTOR 2 FOR NON IDENTICAL PHONONS ADDED ON THE 10/03/81.             VIBM-183
   12 K1=IABS(NBETA(17,N2)-NBETA(17,N1))+1                              VIBM-184
      K2=NBETA(17,N2)+NBETA(17,N1)+1                                    VIBM-185
      FS=DFLOAT(2*(1-2*MOD(NBETA(17,N1)+IABS(IPI(1,J1)-IPI(1,J2)+K1-1)/2VIBM-186
     1,2)))                                                             VIBM-187
      DO 13 K=K1,K2,2                                                   VIBM-188
      IF (2*IT.GT.IDT) CALL MEMO('VIBM',IDT,2*IT)                       VIBM-189
      J=K-1                                                             VIBM-190
      AQ=FS*DJCG(IPI(3,J1)-1,IPI(3,J2)-1,2*J,0,0,FAC,NFA)               VIBM-191
      IQ(1,IT)=MAX0(N1,N2)*(NBT1+1)+MIN0(N1,N2)                         VIBM-192
      IQ(2,IT)=J                                                        VIBM-193
      IQ(3,IT)=NSP                                                      VIBM-194
      T(3,IT)=AQ*AY                                                     VIBM-195
      FS=-FS                                                            VIBM-196
   13 IT=IT+1                                                           VIBM-197
      GO TO 26                                                          VIBM-198
C  (L1,L2,IP||Q||I) WITH I=J1 AND IP=J2.                                VIBM-199
   14 I=IPH(2,J2)                                                       VIBM-200
      LB1=NVAR(1,I)                                                     VIBM-201
      LB2=NVAR(2,I)                                                     VIBM-202
      N1=IPH(2,J1)                                                      VIBM-203
      IF (IPH(1,J1).GT.2) N1=NVAR(1,N1+1)                               VIBM-204
      IF (IPI(3,J1).NE.(2*NBETA(17,N1)+1)) GO TO 41                     VIBM-205
      LIB(1)=LB1.EQ.N1                                                  VIBM-206
      LIB(2)=LB2.EQ.N1                                                  VIBM-207
      IF (2*IT.GT.IDT) CALL MEMO('VIBM',IDT,2*IT)                       VIBM-208
      IF (LIB(1).AND.LIB(2)) GO TO 16                                   VIBM-209
      IF (LIB(2)) GO TO 15                                              VIBM-210
      IF (.NOT.LIB(1)) GO TO 26                                         VIBM-211
      IQ(1,IT)=LB2                                                      VIBM-212
      IQ(2,IT)=NBETA(17,LB2)                                            VIBM-213
      IQ(3,IT)=NSP                                                      VIBM-214
      T(3,IT)=AY*DFLOAT(1-2*MOD(IABS(IPI(1,J1)-IPI(1,J2)+IQ(2,IT))/2,2))VIBM-215
      GO TO 17                                                          VIBM-216
   15 IQ(1,IT)=LB1                                                      VIBM-217
      IQ(2,IT)=NBETA(17,LB1)                                            VIBM-218
      IQ(3,IT)=NSP                                                      VIBM-219
      T(3,IT)=AY*DFLOAT(1-MOD(IABS(IPI(3,J1)+IPI(3,J2)-2+IPI(1,J1)-IPI(1VIBM-220
     1,J2)-IQ(2,IT)),4))                                                VIBM-221
      GO TO 17                                                          VIBM-222
   16 IF (MOD(IPI(3,J2),4).NE.1) GO TO 26                               VIBM-223
      IQ(1,IT)=N1                                                       VIBM-224
      IQ(2,IT)=NBETA(17,N1)                                             VIBM-225
      IQ(3,IT)=NSP                                                      VIBM-226
      T(3,IT)=DSQRT(2.D0)*AY*DFLOAT(1-MOD(IABS(IPI(1,J1)-IPI(1,J2)+IQ(2,VIBM-227
     1IT)),4))                                                          VIBM-228
   17 T(3,IT)=T(3,IT)*DSQRT(DFLOAT(IPI(3,J2))/DFLOAT(2*IQ(2,IT)+1))     VIBM-229
      GO TO 25                                                          VIBM-230
C  (L3,L4,IP||Q||L1,L2,I).                                              VIBM-231
   18 IF (.NOT.LO(2)) GO TO 26                                          VIBM-232
      I=IPH(2,I1)                                                       VIBM-233
      LB1=NVAR(1,I)                                                     VIBM-234
      LB2=NVAR(2,I)                                                     VIBM-235
      I=IPH(2,I2)                                                       VIBM-236
      LB3=NVAR(1,I)                                                     VIBM-237
      LB4=NVAR(2,I)                                                     VIBM-238
      LIB(1)=LB1.NE.LB3                                                 VIBM-239
      LIB(2)=LB2.NE.LB4                                                 VIBM-240
      LIB(3)=LB1.NE.LB4                                                 VIBM-241
      LIB(4)=LB2.NE.LB3                                                 VIBM-242
      IF (LIB(1).AND.LIB(2).AND.LIB(3).AND.LIB(4)) GO TO 26             VIBM-243
      JA1=(IPI(3,J1)-1)/2                                               VIBM-244
      JA2=(IPI(3,J2)-1)/2                                               VIBM-245
      IA(1,1)=NBETA(17,LB2)                                             VIBM-246
      IA(2,1)=NBETA(17,LB4)                                             VIBM-247
      IA(1,2)=NBETA(17,LB1)                                             VIBM-248
      IA(2,2)=NBETA(17,LB3)                                             VIBM-249
      IA(1,3)=IA(1,1)                                                   VIBM-250
      IA(2,3)=IA(2,2)                                                   VIBM-251
      IA(1,4)=IA(1,2)                                                   VIBM-252
      IA(2,4)=IA(2,1)                                                   VIBM-253
      IA(3,1)=IA(1,2)                                                   VIBM-254
      IA(3,2)=IA(1,1)                                                   VIBM-255
      IA(3,3)=IA(1,2)                                                   VIBM-256
      IA(3,4)=IA(1,1)                                                   VIBM-257
      IA(6,1)=MAX0(LB2,LB4)*(NBT1+1)+MIN0(LB2,LB4)                      VIBM-258
      IA(6,2)=MAX0(LB1,LB3)*(NBT1+1)+MIN0(LB1,LB3)                      VIBM-259
      IA(6,3)=MAX0(LB2,LB3)*(NBT1+1)+MIN0(LB2,LB3)                      VIBM-260
      IA(6,4)=MAX0(LB1,LB4)*(NBT1+1)+MIN0(LB1,LB4)                      VIBM-261
      IMIN=1000                                                         VIBM-262
      IMAX=0                                                            VIBM-263
      DO 19 I=1,4                                                       VIBM-264
      IF (LIB(I)) GO TO 19                                              VIBM-265
      IA(4,I)=IABS(IA(1,I)-IA(2,I))                                     VIBM-266
      IA(5,I)=IA(1,I)+IA(2,I)                                           VIBM-267
      IF (IA(4,I).LT.IMIN) IMIN=IA(4,I)                                 VIBM-268
      IF (IA(5,I).GT.IMAX) IMAX=IA(5,I)                                 VIBM-269
   19 CONTINUE                                                          VIBM-270
      B(1)=DFLOAT(1-2*MOD(JA2+IA(1,2),2))                               VIBM-271
      B(2)=DFLOAT(1-2*MOD(IA(1,2)+IA(2,2)+IA(2,1)+JA1,2))               VIBM-272
      B(3)=DFLOAT(1-2*MOD(IA(2,2),2))                                   VIBM-273
      B(4)=DFLOAT(1-2*MOD(JA2+JA1+IA(1,2),2))                           VIBM-274
      T0=DSQRT(DFLOAT((2*JA1+1)*(2*JA2+1)))*2.D0                        VIBM-275
      IF (LB1.EQ.LB2) T0=DSQRT(0.5D0)*T0                                VIBM-276
      IF (LB3.EQ.LB4) T0=DSQRT(0.5D0)*T0                                VIBM-277
      IF (JA1.NE.JA2) GO TO 21                                          VIBM-278
      TKQ=0.D0                                                          VIBM-279
      IF ((LB1.EQ.LB3).AND.(LB2.EQ.LB4)) TKQ=1.D0                       VIBM-280
      IF ((LB1.EQ.LB4).AND.(LB2.EQ.LB3)) TKQ=TKQ+DFLOAT(1-2*MOD(JA1+IA(1VIBM-281
     1,1)+IA(1,2),2))                                                   VIBM-282
      IF (TKQ.EQ.0.D0) GO TO 21                                         VIBM-283
      IF (LB1.EQ.LB2) TKQ=0.5D0*TKQ                                     VIBM-284
      DO 20 L=1,NBT1                                                    VIBM-285
      IF (NBETA(18,L).NE.0) GO TO 20                                    VIBM-286
      IF (2*IT.GT.IDT) CALL MEMO('VIBM',IDT,2*IT)                       VIBM-287
      IQ(1,IT)=L*(NBT1+2)                                               VIBM-288
      IQ(2,IT)=0                                                        VIBM-289
      IQ(3,IT)=NSP                                                      VIBM-290
      T(3,IT)=TKQ*DSQRT(2.D0*JA1+1.D0)*AY                               VIBM-291
      IT=IT+1                                                           VIBM-292
   20 CONTINUE                                                          VIBM-293
   21 IMIN=MAX0(IMIN,IABS(JA1-JA2))+1                                   VIBM-294
      IMAX=MIN0(IMAX,JA1+JA2)+1                                         VIBM-295
      DO 23 M=IMIN,IMAX                                                 VIBM-296
      J=M-1                                                             VIBM-297
      DO 22 K=1,4                                                       VIBM-298
      IF (LIB(K).OR.(J.LT.IA(4,K)).OR.(J.GT.IA(5,K)).OR.(MOD(J+IA(5,K),2VIBM-299
     1).NE.0)) GO TO 22                                                 VIBM-300
      IF (2*IT.GT.IDT) CALL MEMO('VIBM',IDT,2*IT)                       VIBM-301
      T3=B(K)*DJ6J(2*IA(1,K),2*IA(2,K),2*J,2*JA2,2*JA1,2*IA(3,K),FAC,NFAVIBM-302
     1)*DJCG(2*IA(1,K),2*IA(2,K),2*J,0,0,FAC,NFA)*DFLOAT(1-MOD(IABS(IPI(VIBM-303
     21,J1)-IPI(1,J2)+J),4))                                            VIBM-304
      IF (DABS(T3).LT.1.D-6) GO TO 22                                   VIBM-305
      IQ(1,IT)=IA(6,K)                                                  VIBM-306
      IQ(2,IT)=J                                                        VIBM-307
      IQ(3,IT)=NSP                                                      VIBM-308
      T(3,IT)=T3*T0*AY                                                  VIBM-309
      IT=IT+1                                                           VIBM-310
   22 CONTINUE                                                          VIBM-311
   23 CONTINUE                                                          VIBM-312
      GO TO 26                                                          VIBM-313
C  (L1,L2,IP||Q||0) WITH IP=J2.                                         VIBM-314
   24 IF (.NOT.LO(2)) GO TO 26                                          VIBM-315
      K3=(IPI(3,J2)-1)/2                                                VIBM-316
      IQ(2,IT)=K3                                                       VIBM-317
      I=IPH(2,J2)                                                       VIBM-318
      K1=NVAR(1,I)                                                      VIBM-319
      K2=NVAR(2,I)                                                      VIBM-320
      IQ(1,IT)=MAX0(K1,K2)*(NBT1+1)+MIN0(K1,K2)                         VIBM-321
      IQ(3,IT)=NSP                                                      VIBM-322
      W1=2.D0                                                           VIBM-323
      IF (K1.EQ.K2) W1=DSQRT(W1)                                        VIBM-324
      T(3,IT)=AY*DJCG(2*NBETA(17,K1),2*NBETA(17,K2),IPI(3,J2)-1,0,0,FAC,VIBM-325
     1NFA)*DFLOAT(1-MOD(IABS(IPI(1,J1)-IPI(1,J2)+K3),4))*W1             VIBM-326
   25 IT=IT+1                                                           VIBM-327
   26 IF ((AA(2,2).EQ.0.D0).OR.(L2.EQ.2)) GO TO 27                      VIBM-328
      L2=2                                                              VIBM-329
      AY=AX*AA(2,2)                                                     VIBM-330
      GO TO 4                                                           VIBM-331
   27 IF ((AA(1,2).EQ.0.D0).OR.(L1.EQ.2)) GO TO 28                      VIBM-332
      L1=2                                                              VIBM-333
      AX=AA(1,2)                                                        VIBM-334
      GO TO 3                                                           VIBM-335
   28 ITI=NIV(I1,I2,1)                                                  VIBM-336
      IF (ITI.EQ.IT) GO TO 37                                           VIBM-337
      ITF=IT-1                                                          VIBM-338
      IF (ITF.EQ.ITI) GO TO 32                                          VIBM-339
      IT1=ITI+1                                                         VIBM-340
      IT=ITI                                                            VIBM-341
      DO 31 I=IT1,ITF                                                   VIBM-342
      DO 29 J=ITI,IT                                                    VIBM-343
      IF ((IQ(1,I).NE.IQ(1,J)).OR.(IQ(2,I).NE.IQ(2,J)).OR.(IQ(3,I).NE.IQVIBM-344
     1(3,J))) GO TO 29                                                  VIBM-345
      T(3,J)=T(3,J)+T(3,I)                                              VIBM-346
      GO TO 31                                                          VIBM-347
   29 CONTINUE                                                          VIBM-348
      IT=IT+1                                                           VIBM-349
      DO 30 K=1,3                                                       VIBM-350
   30 IQ(K,IT)=IQ(K,I)                                                  VIBM-351
      T(3,IT)=T(3,I)                                                    VIBM-352
   31 CONTINUE                                                          VIBM-353
      IT=IT+1                                                           VIBM-354
   32 ITF=IT-1                                                          VIBM-355
      DO 36 I=ITI,ITF                                                   VIBM-356
      IF (IT.LE.I) GO TO 37                                             VIBM-357
   33 IF (DABS(T(3,I)).GT.1.D-12) GO TO 36                              VIBM-358
      IT1=I+1                                                           VIBM-359
      IT=IT-1                                                           VIBM-360
      IF (IT1.GT.IT) GO TO 37                                           VIBM-361
      DO 35 J=IT1,IT                                                    VIBM-362
      DO 34 K=1,3                                                       VIBM-363
   34 IQ(K,J-1)=IQ(K,J)                                                 VIBM-364
   35 T(3,J-1)=T(3,J)                                                   VIBM-365
      GO TO 33                                                          VIBM-366
   36 CONTINUE                                                          VIBM-367
   37 NIV(I1,I2,2)=IT-1                                                 VIBM-368
   38 CONTINUE                                                          VIBM-369
   39 CONTINUE                                                          VIBM-370
      IT=IT-1                                                           VIBM-371
      RETURN                                                            VIBM-372
   40 WRITE (MW,1001) JVAR,NVA                                          VIBM-373
      GO TO 42                                                          VIBM-374
   41 WRITE (MW,1002) I1,I2                                             VIBM-375
   42 WRITE (MW,1003)                                                   VIBM-376
      STOP                                                              VIBM-377
 1000 FORMAT (' STATE',I4,F15.5,' DEGREES      AMPLITUDES =',F15.7,'  1 VIBM-378
     1PHONON AND',F15.7,'  2 PHONONS.')                                 VIBM-379
 1001 FORMAT (' NUMBER OF VARIABLES USED:',I5,5X,'EXCEEDS NUMBER OF VARIVIBM-380
     1ABLES READ:',I6)                                                  VIBM-381
 1002 FORMAT (' INCORRECT DESCRIPTION OF LEVEL',I3,'  OR',I3)           VIBM-382
 1003 FORMAT (//' IN VIBM  ...  STOP  ...')                             VIBM-383
      END                                                               VIBM-384
