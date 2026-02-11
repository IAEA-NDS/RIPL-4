C 07/03/06                                                      ECIS06  STBF-000
      SUBROUTINE STBF(V,IVX,ISM,VAL,IVM,NX,IDX,XX,IV,HH,IEX,LO)         STBF-001
C COMPUTES BOUND STATE FORM FACTORS.                                    STBF-002
C INPUT:     IVX:     TABLE OF QUANTUM NUMBERS.                         STBF-003
C            ISM:     NUMBER OF INTEGRATION STEPS.                      STBF-004
C            VAL(J):  OSCILLATOR PARAMETER OR BINDING ENERGY FOR J=1,   STBF-005
C                     TOTAL MASS FOR J=2,                               STBF-006
C                     MASS OF BOUND PARTICLE FOR J=3,                   STBF-007
C                     PRODUCT OF CHARGES FOR J=4,                       STBF-008
C                     REAL POTENTIAL FOR UNBOUND STATES OR STARTING     STBF-009
C                     VALUE FOR THE SEARCH ON BOUND STATE FOR J=5,      STBF-010
C                     REDUCED RADIUS OF REAL POTENTIAL FOR J=6,         STBF-011
C                     DIFFUSENESS OF RADIUS OF REAL POTENTIAL FOR J=7,  STBF-012
C                     DEPTH OF SPIN-ORBIT POTENTIAL FOR J=8,            STBF-013
C                     REDUCED RADIUS OF SPIN-ORBIT POTENTIAL FOR J=9,   STBF-014
C                     DIFFUSENESS OF SPIN-ORBIT POTENTIAL FOR J=10,     STBF-015
C                     REDUCED RADIUS OF COULOMB POTENTIAL FOR J=11.     STBF-016
C            IVM:     STEP FACTOR FOR WOODS-SAXON WAVE FUNCTIONS.       STBF-017
C            NX:      LENGTH OF WORKING SPACE.                          STBF-018
C            IV:      IV=8 FOR SOLUTION IN A WOODS-SAXON POTENTIAL.     STBF-019
C            HH:      STEP SIZE FOR INTEGRATION.                        STBF-020
C            LO(I):   LOGICAL CONTROLS:                                 STBF-021
C               LO(51) =.TRUE. OUTPUT OF POTENTIALS.                    STBF-022
C               LO(47) =.TRUE. NO RECOIL CORRECTION FOR BOUND STATES.   STBF-023
C OUTPUT:    V(ISM):  BOUND STATE WAVE FUNCTION.                        STBF-024
C            IDX:     LENGTH OF WORKING FIELD USED.                     STBF-025
C WORKING AREA:                                                         STBF-026
C            XX:      TO COMPUTE BOUND FUNCTIONS.                       STBF-027
C            IEX:     FOR COULOMB FUNCTION, IN EQUIVALENCE WITH XX.     STBF-028
C                                                                       STBF-029
C FOR THE COMMON  /DCONS/ SEE CALC.                                     STBF-030
C                                                                       STBF-031
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     STBF-032
C  CCZ:       COULOMB ALPHA CONSTANT.                                   STBF-033
C  CK:        H-BAR*C.                                                  STBF-034
C   USED:     CCZ,CK.                                                   STBF-035
C                                                                       STBF-036
C***********************************************************************STBF-037
      IMPLICIT REAL*8 (A-H,O-Z)                                         STBF-038
      LOGICAL LO(150)                                                   STBF-039
      DIMENSION V(*),IVX(4),VAL(11),XX(*),IEX(*)                        STBF-040
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            STBF-041
      COMMON /INOUT/ MR,MW,MS                                           STBF-042
      L=IVX(2)                                                          STBF-043
      H=HH                                                              STBF-044
      IF ((.NOT.LO(47)).AND.(VAL(3).NE.VAL(2))) H=H*VAL(2)/(VAL(2)-VAL(3STBF-045
     1))                                                                STBF-046
      IF (IV.EQ.8) GO TO 7                                              STBF-047
      N=IVX(1)                                                          STBF-048
C COMPUTATION OF THE NORMALISATION OF LAGUERRE POLYNOMIALS              STBF-049
      A1=2.256758334191D0                                               STBF-050
      A3=L                                                              STBF-051
      IF (L.EQ.0) GO TO 2                                               STBF-052
      DO 1 I=1,L                                                        STBF-053
    1 A1=A1/(DFLOAT(I)+.5D0)                                            STBF-054
    2 IF (N.EQ.0) GO TO 4                                               STBF-055
      DO 3 I=1,N                                                        STBF-056
      A2=DFLOAT(I)                                                      STBF-057
    3 A1=A1*(A2+A3+.5D0)/A2                                             STBF-058
    4 A1=DSQRT(A1)                                                      STBF-059
C COMPUTATION OF LAGUERRE POLYNOMIALS                                   STBF-060
      RR=0.D0                                                           STBF-061
      DO 6 IS=1,ISM                                                     STBF-062
      RR=RR+H                                                           STBF-063
      X1=(VAL(1)*RR)**2                                                 STBF-064
      X3=A1*DEXP(.5D0*(A3*DLOG(X1)-X1))                                 STBF-065
      IF (N.EQ.0) GO TO 6                                               STBF-066
      X2=X3                                                             STBF-067
      DO  5 I=1,N                                                       STBF-068
      A2=DFLOAT(I)                                                      STBF-069
      X2=2.D0*X2*DFLOAT(I-N-1)*X1/(A2*(2.D0*(A3+A2)+1.D0))              STBF-070
    5 X3=X3+X2                                                          STBF-071
    6 V(IS)=X3                                                          STBF-072
      RETURN                                                            STBF-073
C BOUND STATE IN WOODS-SAXON POTENTIAL                                  STBF-074
C CONTROL OF STORAGE AND STEP SIZE                                      STBF-075
    7 IREB=IVM                                                          STBF-076
      IF (IREB.EQ.0) IREB=4                                             STBF-077
      KN=IREB*ISM                                                       STBF-078
      IDX=MAX0(IDX,3*KN+3)                                              STBF-079
      IF (IDX.GT.NX) CALL MEMO('STBF',NX,IDX)                           STBF-080
      A4=IREB                                                           STBF-081
      A4=H/A4                                                           STBF-082
      Y3=(VAL(2)-VAL(3))**0.33333333333333D0                            STBF-083
      A2=Y3*VAL(6)                                                      STBF-084
      KR=IDINT(A2/A4)                                                   STBF-085
      Y2=VAL(1)                                                         STBF-086
      IF (Y2.EQ.0.D0) Y2=0.01D0                                         STBF-087
      JL=IVX(4)                                                         STBF-088
      IF (IVX(3).EQ.0) JL=2*L                                           STBF-089
      IF (LO(51)) WRITE (MW,1000) IVX,KN,IREB,KR,Y2,(VAL(I),I=2,11)     STBF-090
      A1=A4*A4                                                          STBF-091
      Y1=VAL(3)*(VAL(2)-VAL(3))/VAL(2)                                  STBF-092
      Y5=CK*Y1*A1                                                       STBF-093
      A5=Y2*Y5                                                          STBF-094
      RR=DSQRT(DABS(A5))/A4                                             STBF-095
      X5=0.5D0*CK*CCZ*Y1*VAL(4)/RR                                      STBF-096
      IF (VAL(1).GT.0.D0) GO TO 8                                       STBF-097
C MATCHING OF UNBOUND FUNCTIONS                                         STBF-098
      IDX=MAX0(IDX,6*L+6)                                               STBF-099
      IF (IDX.GT.NX) CALL MEMO('STBF',NX,IDX)                           STBF-100
      CALL FCOU(L,X5,RR*A4*KN,XX(L+2),XX(2*L+3),XX(3*L+4),XX(4*L+5),IEX,STBF-101
     1XX(5*L+6))                                                        STBF-102
      Z1=XX(2*L+2)                                                      STBF-103
      Z2=XX(4*L+4)                                                      STBF-104
      CALL FCOU(L,X5,RR*A4*(KN-2),XX(L+2),XX(2*L+3),XX(3*L+4),XX(4*L+5),STBF-105
     1IEX,XX(5*L+6))                                                    STBF-106
      Z3=XX(2*L+2)                                                      STBF-107
      Z4=XX(4*L+4)                                                      STBF-108
    8 X1=Y3*VAL(11)                                                     STBF-109
      A3=Y3*VAL(9)                                                      STBF-110
      Y3=2.D0*VAL(8)*Y5                                                 STBF-111
      Y4=IVX(3)*(IVX(3)+2)-JL*(JL+2)+4*L*(L+1)                          STBF-112
      Y4=.25D0*Y4                                                       STBF-113
      A6=DFLOAT(L*(L+1))                                                STBF-114
      A6=A6*A1                                                          STBF-115
      X4=X5*2.D0*A1*RR                                                  STBF-116
      X3=0.D0                                                           STBF-117
      K1=KN+1                                                           STBF-118
      K2=KN+2                                                           STBF-119
C COMPUTATION OF OPTICAL POTENTIALS                                     STBF-120
      DO 10 I=1,K1                                                      STBF-121
      X3=X3+A4                                                          STBF-122
      X2=DEXP((X3-A3)/VAL(10))                                          STBF-123
      XX(K2+I)=A6/(X3*X3)+Y4*Y3*X2/((1.D0+X2)**2*X3*VAL(10))+A5         STBF-124
      XX(2*K2+I)=-1.D0/(1.D0+DEXP((X3-A2)/VAL(7)))                      STBF-125
      IF (X3.GT.X1) GO TO 9                                             STBF-126
      XX(K2+I)=XX(K2+I)+X4*(1.5D0-0.5D0*X3*X3/(X1*X1))/X1               STBF-127
      GO TO 10                                                          STBF-128
    9 XX(K2+I)=XX(K2+I)+X4/X3                                           STBF-129
   10 CONTINUE                                                          STBF-130
      IF (VAL(1).LT.0.D0) GO TO 15                                      STBF-131
C MATCHING CONDITIONS                                                   STBF-132
      X1=(X3-A4)*RR                                                     STBF-133
      X3=X3*RR                                                          STBF-134
      A6=A6/A1                                                          STBF-135
      N=IDINT(5.D0*(X3-X1)+1.D0)                                        STBF-136
      Y3=N                                                              STBF-137
      Y3=(X3-X1)/Y3                                                     STBF-138
      IG=IDINT(1.D0/Y3)                                                 STBF-139
      JG=MIN0(100,IG)                                                   STBF-140
      A5=Y3*Y3                                                          STBF-141
      A3=1.D0                                                           STBF-142
      DO 12 I=1,20                                                      STBF-143
      M=2*I*JG                                                          STBF-144
      A4=X3+Y3*DFLOAT(M)                                                STBF-145
      Y1=DEXP(Y3)                                                       STBF-146
      Y4=0.5D0/Y3+0.5D0                                                 STBF-147
      Y2=(A6/(A4*A4)+2.D0*X5/A4+1.D0)*A5                                STBF-148
      DO 11 J=1,M                                                       STBF-149
      Y1=(2.D0+Y2/(1.D0-Y2/12.D0))-1.D0/Y1                              STBF-150
      A2=Y2                                                             STBF-151
      A4=A4-Y3                                                          STBF-152
      Y2=(A6/(A4*A4)+2.D0*X5/A4+1.D0)*A5                                STBF-153
      A1=Y1*(1.D0-A2/12.D0)/(1.D0-Y2/12.D0)                             STBF-154
   11 Y4=Y4/(A1*A1)+1.D0                                                STBF-155
      IF (DABS(Y1-A3).LT.0.1D-4*DABS(Y1-1.D0)) GO TO 13                 STBF-156
   12 A3=Y1                                                             STBF-157
      IF (LO(51)) WRITE (MW,1001)                                       STBF-158
   13 Y4=(Y4-0.5D0)*Y3                                                  STBF-159
      IF (LO(51)) WRITE (MW,1002) I,M                                   STBF-160
      X2=1.D0/Y1                                                        STBF-161
      X3=1.D0                                                           STBF-162
      X4=Y2                                                             STBF-163
      DO 14 I=1,N                                                       STBF-164
      X1=X2                                                             STBF-165
      X2=X3                                                             STBF-166
      X3=(2.D0+X4/(1.D0-X4/12.D0))*X2-X1                                STBF-167
      A4=A4-Y3                                                          STBF-168
   14 X4=(A6/(A4*A4)+2.D0*X5/A4+1.D0)*A5                                STBF-169
      A3=X3*(1.D0-Y2/12.D0)/(1.D0-X4/12.D0)                             STBF-170
      GO TO 16                                                          STBF-171
   15 KR=K1                                                             STBF-172
   16 KS=IDINT(2.D0+DSQRT(A6/12.D0))                                    STBF-173
C STARTING VALUES                                                       STBF-174
      DO 17 I=1,KS                                                      STBF-175
   17 XX(I)=0.D0                                                        STBF-176
      XX(KS)=1.D0                                                       STBF-177
      IF (L.EQ.1) XX(KS-1)=-.2D0                                        STBF-178
      A4=VAL(5)                                                         STBF-179
      IK=2*IVX(1)                                                       STBF-180
      IST=0                                                             STBF-181
      X1=0.D0                                                           STBF-182
      Y1=-1.D0                                                          STBF-183
      IN1=0                                                             STBF-184
      IN2=0                                                             STBF-185
C SEARCH FOR THE EIGENVALUE                                             STBF-186
   18 X3=X1+A4                                                          STBF-187
   19 IST=IST+1                                                         STBF-188
      IN3=0                                                             STBF-189
      X4=X3*Y5                                                          STBF-190
C UPWARDS INTEGRATION                                                   STBF-191
      DO 20 I=KS,KR                                                     STBF-192
      X5=X4*XX(2*K2+I-1)+XX(K2+I-1)                                     STBF-193
      XX(I+1)=(2.D0+X5/(1.D0-X5/12.D0))*XX(I)-XX(I-1)                   STBF-194
      IF (XX(I+1)*XX(I).LT.0.D0) IN3=IN3+2                              STBF-195
   20 CONTINUE                                                          STBF-196
      IF (VAL(1).LT.0.D0) GO TO 29                                      STBF-197
      A1=XX(KR)                                                         STBF-198
      A2=XX(KR+1)                                                       STBF-199
      XX(K2)=1.D0                                                       STBF-200
      XX(K1)=A3*(1.D0-(X4*XX(3*KN+4)+XX(2*KN+2))/12.D0)/(1.D0-(X4*XX(3*KSTBF-201
     1N+5)+XX(2*KN+3))/12.D0)                                           STBF-202
C BACKWARDS INTEGRATION                                                 STBF-203
      DO 21 I=KR,KN                                                     STBF-204
      J=KN+KR-I                                                         STBF-205
      X5=X4*XX(2*K2+J)+XX(K2+J)                                         STBF-206
      IF (XX(J+1)*XX(J+2).LT.0.D0) IN3=IN3+2                            STBF-207
   21 XX(J)=(2.D0+X5/(1.D0-X5/12.D0))*XX(J+1)-XX(J+2)                   STBF-208
      Y3=A1/A2-XX(KR)/XX(KR+1)                                          STBF-209
      IF (Y3.GT.0.D0) IN3=IN3+1                                         STBF-210
      IF (LO(51)) WRITE (MW,1003) IST,IN3,X3,Y3                         STBF-211
      IF (DABS(Y3).LE.1.D-10.OR.X2.EQ.X3.OR.X1.EQ.X3) GO TO 27          STBF-212
      IF (IN2.NE.0) GO TO 24                                            STBF-213
      IF (IN3.GT.IK) GO TO 22                                           STBF-214
      Y1=Y3                                                             STBF-215
      X1=X3                                                             STBF-216
      IN1=IN3                                                           STBF-217
      A4=X1                                                             STBF-218
      GO TO 18                                                          STBF-219
   22 IF (LO(51)) WRITE (MW,1004) IK,IN1,IN3                            STBF-220
      X2=X3                                                             STBF-221
      IN2=IN3                                                           STBF-222
      Y2=Y3                                                             STBF-223
C INTERPOLATION OF SOLUTION                                             STBF-224
   23 X3=X1+0.5D0*(X2-X1)                                               STBF-225
      IF (IN2.EQ.IN1+1) X3=(Y1*X2-Y2*X1)/(Y1-Y2)                        STBF-226
      GO TO 19                                                          STBF-227
   24 X4=(Y1*Y2*X3*(X2-X1)+Y1*Y3*X2*(X1-X3)+Y2*Y3*X1*(X3-X2))/(Y1*Y2*(X2STBF-228
     1-X1)+Y1*Y3*(X1-X3)+Y2*Y3*(X3-X2))                                 STBF-229
      IF (IN2.GT.IN1+2) X4=1.D20                                        STBF-230
      IF (IN3.GT.IK) GO TO 25                                           STBF-231
      X1=X3                                                             STBF-232
      Y1=Y3                                                             STBF-233
      IN1=IN3                                                           STBF-234
      GO TO 26                                                          STBF-235
   25 X2=X3                                                             STBF-236
      Y2=Y3                                                             STBF-237
      IN2=IN3                                                           STBF-238
   26 IF ((X1-X4)*(X4-X2).LT.0.D0) GO TO 23                             STBF-239
      X3=X4                                                             STBF-240
      GO TO 19                                                          STBF-241
   27 Y3=A1/XX(KR)                                                      STBF-242
      DO 28 I=KR,K2                                                     STBF-243
   28 XX(I)=Y3*XX(I)                                                    STBF-244
   29 DO 30 I=1,KN                                                      STBF-245
   30 XX(I)=(XX(I)+10.D0*XX(I+1)+XX(I+2))/12.D0                         STBF-246
      IF (VAL(1).GT.0.D0) GO TO 31                                      STBF-247
      X1=XX(KN)*Z4-XX(KN-2)*Z2                                          STBF-248
      X2=XX(KN-2)*Z1-XX(KN)*Z3                                          STBF-249
      X3=Z1*Z4-Z2*Z3                                                    STBF-250
      Y2=DATAN2(X2,X1)                                                  STBF-251
      Y1=(X1**2+X2**2)/X3**2                                            STBF-252
      X3=VAL(5)                                                         STBF-253
      GO TO 33                                                          STBF-254
   31 Y1=0.D0                                                           STBF-255
      DO 32 I=1,KN                                                      STBF-256
   32 Y1=Y1+XX(I)**2                                                    STBF-257
      Y1=Y1-0.5D0*XX(KN)**2                                             STBF-258
      Y2=Y4*XX(KN)**2/RR                                                STBF-259
      Y3=IREB                                                           STBF-260
      Y1=Y1*H/Y3+Y2                                                     STBF-261
      Y2=100.D0*Y2/Y1                                                   STBF-262
   33 Y1=1.D0/DSQRT(Y1)                                                 STBF-263
      IS=0                                                              STBF-264
      RR=0.D0                                                           STBF-265
      DO 34 I=IREB,KN,IREB                                              STBF-266
      IS=IS+1                                                           STBF-267
      RR=RR+H                                                           STBF-268
   34 V(IS)=Y1*XX(I)/RR                                                 STBF-269
      IF (LO(51)) WRITE (MW,1005) X3,Y2,(IS,V(IS),IS=1,ISM)             STBF-270
      RETURN                                                            STBF-271
 1000 FORMAT (//' WOODS-SAXON POTENTIAL EIGENFUNCTION WITH N =',I2,3X,'LSTBF-272
     1 =',I3,3X,'2*S =',I2,3X,'2*J =',I3,3X,I4,' STEPS (DIVIDED BY',I3,'STBF-273
     2) MATCHING AT THE',I4,'TH'/' **** BINDING ENERGY',F12.6,' MEV ****STBF-274
     3',4X,'TOTAL MASS',F12.6,4X,'PARTICLE MASS',F12.6,4X,'PRODUCT OF CHSTBF-275
     4ARGES',F8.2/' SEARCH ON DEPTH OF REAL POTENTIAL FROM',F12.6,' WITHSTBF-276
     5 REDUCED RADIUS',F10.6,' FERMI AND DIFFUSENESS',F10.6,' FERMI'/' SSTBF-277
     6PIN-ORBIT POTENTIAL DEPTH:',F12.6,' MEV RADIUS:',F10.6,' FERMI DIFSTBF-278
     7FUSENESS:',F9.6,' FERMI COULOMB RADIUS:',F10.6,' FERMI')          STBF-279
 1001 FORMAT (10X,'NO CONVERGENCE ON MATCHING VALUES')                  STBF-280
 1002 FORMAT (10X,'MATCHING VALUE OBTAINED FOR I=',I3,' WITH',I6,' POINTSTBF-281
     1S')                                                               STBF-282
 1003 FORMAT (2X,I3,5X,'2*N (+0/1) =',I3,5X,'V =',D20.10,5X,'D =',D20.10STBF-283
     1)                                                                 STBF-284
 1004 FORMAT (' INTERPOLATION FOR 2*N =',I4,5X,'BETWEEN',I4,' AND',I4)  STBF-285
 1005 FORMAT (//' DEPTH USED',D20.10,' MEV',15X,'TAIL PERCENTAGE OR PHASSTBF-286
     1E-SHIFT',F10.5//(6(I6,D15.6)))                                    STBF-287
      END                                                               STBF-288
