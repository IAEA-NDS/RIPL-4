C 07/03/07                                                      ECIS06  SCAT-000
      SUBROUTINE SCAT(FR,MF,JMAX,KMAX,THETA,FCN,IPJ,IPK,NCO,COE,IPI,AB,ASCAT-001
     1,B,EX,WV,NCJ,JMIN,LO)                                             SCAT-002
C COMPUTATION OF OBSERVABLES.                                           SCAT-003
C INPUT:     FR:      SCATTERING COEFFICIENTS IN THE HELICITY FORMALISM.SCAT-004
C            MF:      TABLES OF HELICITY,DESCRIPTION OF OBSERVABLES,    SCAT-005
C                     .... ETC.  SEE DEPH AND OBSE.                     SCAT-006
C            JMAX:    DIMENSION FOR FR.                                 SCAT-007
C            KMAX:    DIMENSION FOR FCN.                                SCAT-008
C            THETA:   SCATTERING ANGLE.                                 SCAT-009
C            FCN:     COMPOUND NUCLEUS COEFFICIENTS.                    SCAT-010
C            IPJ:     NUMBER OF J VALUES USED FOR SCATTERING.           SCAT-011
C            IPK:     NUMBER OF L VALUES USED FOR COMPOUND NUCLEUS.     SCAT-012
C            NCO,COE: LOOPS AND COEFFICIENTS FOR OBSERVABLES (SEE OBSE).SCAT-013
C            IPI:     MULTIPLICITIES FOR PARTICLE AND TARGET NUCLEUS.   SCAT-014
C            WV:      MASSES, ENERGIES, ETC..    SEE CALX.              SCAT-015
C            NCJ:     NUMBER OF FACTORISATIONS OF 1/(1-X*COS(THETA)).   SCAT-016
C            JMIN:    TWICE MINIMUM OF THE TOTAL SPIN.                  SCAT-017
C            LO:      LOGICAL CONTROLS:                                 SCAT-018
C               LO(8)  =.TRUE. RELATIVISTIC KINEMATICS.                 SCAT-019
C               LO(18) =.TRUE. PROJECTILE-TARGET ANTISYMMETRISATION.    SCAT-020
C               LO(41) =.TRUE. FACTORISATION OF 1/(1-COS(THETA)).       SCAT-021
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             SCAT-022
C OUTPUT:    EX:      DIFFERENTIAL CROSS-SECTION FOLLOWED BY OBSERVABLESSCAT-023
C WORKING AREAS:                                                        SCAT-024
C            AB:      AMPLITUDES IN THE C. M. OR THE LAB. SYSTEM.       SCAT-025
C            A:       FOR INDEPENDENT AMPLITUDES.                       SCAT-026
C            B:       TO STORE THE ROTATION MATRIX ELEMENTS.            SCAT-027
C                                                                       SCAT-028
C THE COMMON /RESCT/ IS USED IN RESU AND SCAT.                          SCAT-029
C                                                                       SCAT-030
C FOR THE COMMON  /DCONS/ SEE CALC.                                     SCAT-031
C                                                                       SCAT-032
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     SCAT-033
C  CMB:       ATOMIC MASS CM DIVIDED BY H-BAR*C.                        SCAT-034
C  XZ:        CONVERSION FACTOR TO MILLIBARNS.                          SCAT-035
C   USED:     CMB,XZ.                                                   SCAT-036
C                                                                       SCAT-037
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /RESCT/:                     SCAT-038
C  EXCN:      COMPOUND-NUCLEUS CROSS-SECTION.                           SCAT-039
C  XA,XB:     CONSTANTS FOR CHANGE TO LABORATORY SYSTEM.                SCAT-040
C  JN1,JN2:   FIRST AND LAST AMPLITUDE IN:THE TABLE.                    SCAT-041
C  JN2:       NUMBER OF INDEPENDENT AMPLITUDES.                         SCAT-042
C  NTT:       TOTAL NUMBER OF AMPLITUDES FOR THIS STATE AND A GIVEN J.  SCAT-043
C  MTT:       SIZE OF THE SCATTERING MATRIX FOR THIS STATE AND ANGLE.   SCAT-044
C  MT1:       MULTIPLICITY OF THE OUTGOING PARTICLE.                    SCAT-045
C  MT2:       MULTIPLICITY OF THE RESIDUAL TARGET.                      SCAT-046
C  MT3:       MULTIPLICITY OF THE INCIDENT PARTICLE.                    SCAT-047
C  MT4:       MULTIPLICITY OF THE INITIAL TARGET.                       SCAT-048
C  NOUT:      NUCLEAR STATE CONSIDERED.                                 SCAT-049
C  NIX,NFX:   FIRST AND LAST OBSERVABLE IN THE TABLE.                   SCAT-050
C  KLT:       1 FOR EXPERIMENTAL DATA.                                  SCAT-051
C             2 FOR EQUIDISTANT ANGLES.                                 SCAT-052
C             3 FOR PURE COMPOUND NUCLEUS.                              SCAT-053
C  LKT:       LOGICAL=.TRUE. FOR CENTRE OF MASS SYSTEM.                 SCAT-054
C                                                                       SCAT-055
C***********************************************************************SCAT-056
      IMPLICIT REAL*8 (A-H,O-Z)                                         SCAT-057
      LOGICAL LO(150),LTT(6),LXY,LKT                                    SCAT-058
      DIMENSION FR(2,JMAX,*),MF(10,*),FCN(KMAX,*),NCO(20,*),COE(*),IPI(1SCAT-059
     12,*),AB(2,*),A(2,*),B(*),EX(*),WV(22,*),MO(18)                    SCAT-060
      EQUIVALENCE (MO(1),MI1),(MO(2),MP1),(MO(3),N1),(MO(4),L1),(MO(5),MSCAT-061
     1I2),(MO(6),MP2),(MO(7),N2),(MO(8),L2),(MO(9),MI3),(MO(10),MP3),(MOSCAT-062
     2(11),N3),(MO(12),L3),(MO(13),MI4),(MO(14),MP4),(MO(15),N4),(MO(16)SCAT-063
     3,L4),(MO(17),JQ),(MO(18),JP)                                      SCAT-064
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            SCAT-065
      COMMON /RESCT/ EXCN,XA,XB,JN1,JN2,NTT,MTT,MT1,MT2,MT3,MT4,NOUT,NIXSCAT-066
     1,NFX,KLT,LKT                                                      SCAT-067
      THET=1.74532925D-02*THETA                                         SCAT-068
      IF (KLT.EQ.3) GO TO 56                                            SCAT-069
C FOR EXPERIMENTAL DATA OR EQUIDISTANT ANGLES.                          SCAT-070
      W6=DCOS(THET)                                                     SCAT-071
      W7=DSIN(THET)                                                     SCAT-072
      IF (LKT) GO TO 1                                                  SCAT-073
      XC=(XA**2-XB**2)*W7**2+W6**2                                      SCAT-074
      IF (XC.GT.0.D0) XC=DSQRT(XC)                                      SCAT-075
      A3=XC*W6-XA*XB*W7**2                                              SCAT-076
      W7=W7*(XB*W6+XC*XA)                                               SCAT-077
      W6=A3                                                             SCAT-078
    1 IF (.NOT.LO(81)) GO TO 3                                          SCAT-079
C COMPOUND NUCLEUS CONTRIBUTION.                                        SCAT-080
      EXCN=FCN(1,NOUT)                                                  SCAT-081
      IF (IPK.EQ.1) GO TO 3                                             SCAT-082
      U1=0.D0                                                           SCAT-083
      U2=1.D0                                                           SCAT-084
      DO 2 J=2,IPK                                                      SCAT-085
      U1=(DFLOAT(4*J-7)*U2*W6-DFLOAT(2*J-4)*U1)/DFLOAT(2*J-3)           SCAT-086
      U2=(DFLOAT(4*J-5)*U1*W6-DFLOAT(2*J-3)*U2)/DFLOAT(2*J-2)           SCAT-087
    2 EXCN=EXCN+U2*FCN(J,NOUT)*DFLOAT(4*J-3)                            SCAT-088
C DIRECT INTERACTION CONTRIBUTION.                                      SCAT-089
    3 W3=1.D0                                                           SCAT-090
      X2=DSQRT(.5D0*(1.D0+W6))                                          SCAT-091
      X3=DSQRT(.5D0*(1.D0-W6))                                          SCAT-092
      DO 4 I=1,NTT                                                      SCAT-093
      A(1,I)=0.D0                                                       SCAT-094
    4 A(2,I)=0.D0                                                       SCAT-095
      IF (NOUT.NE.1) GO TO 7                                            SCAT-096
C COMPUTATION OF COULOMB AMPLITUDES FOR THE ELASTIC CHANNEL.            SCAT-097
      IF (WV(5,1).EQ.0.D0.OR.X3.LT.1.D-20) GO TO 7                      SCAT-098
      W2=2.D0*WV(5,1)*DLOG(X3)                                          SCAT-099
      W3=-0.5D0*WV(5,1)/(X3*X3)                                         SCAT-100
      W4=W3*DCOS(W2)                                                    SCAT-101
      W5=-W3*DSIN(W2)                                                   SCAT-102
      W3=W3**2                                                          SCAT-103
      IF (.NOT.LO(18)) GO TO 5                                          SCAT-104
C SYMMETRISATION BETWEEN PROJECTILE AND TARGET.                         SCAT-105
      FS=DFLOAT(2*MOD(IPI(2,1),2)-1)                                    SCAT-106
      Y2=2.D0*WV(5,1)*DLOG(X2)                                          SCAT-107
      Y1=-0.5D0*WV(5,1)/(X2*X2)                                         SCAT-108
      Y4=Y1*DCOS(Y2)                                                    SCAT-109
      Y5=-Y1*DSIN(Y2)                                                   SCAT-110
      W3=W4**2+W5**2+Y4**2+Y5**2+2.D0*FS*(W4*Y4+W5*Y5)/DFLOAT(IPI(2,1)) SCAT-111
C TRANSFORMATION OF THE COULOMB AMPLITUDES TO THE HELICITY FORMALISM.   SCAT-112
    5 DO 6 I=1,NTT                                                      SCAT-113
      M5=IPI(2,1)-1                                                     SCAT-114
      M6=IPI(3,1)-1                                                     SCAT-115
      M1=2*MF(1,I)-1-IPI(2,1)                                           SCAT-116
      M2=IPI(3,1)-2*MF(2,I)+1                                           SCAT-117
      M3=2*MF(3,I)-1-IPI(2,1)                                           SCAT-118
      M4=IPI(3,1)-2*MF(4,I)+1                                           SCAT-119
      CALL EMRO(M5,M3,M1,X2,X3,B,1)                                     SCAT-120
      CALL EMRO(M6,M4,M2,X2,X3,B(2),1)                                  SCAT-121
      W1=B(1)*B(2)                                                      SCAT-122
      A(1,I)=W4*W1                                                      SCAT-123
      A(2,I)=W5*W1                                                      SCAT-124
      IF (.NOT.LO(18)) GO TO 6                                          SCAT-125
      CALL EMRO(M5,M3,M2,X2,X3,B,1)                                     SCAT-126
      CALL EMRO(M6,M4,M1,X2,X3,B(2),1)                                  SCAT-127
      Y1=B(1)*B(2)*FS                                                   SCAT-128
      A(1,I)=A(1,I)+Y4*Y1                                               SCAT-129
      A(2,I)=A(2,I)+Y5*Y1                                               SCAT-130
    6 CONTINUE                                                          SCAT-131
      W3=W3/WV(4,1)**2                                                  SCAT-132
C COMPUTATION OF THE NUCLEAR AMPLITUDES.                                SCAT-133
    7 DO 13 I=JN1,JN2                                                   SCAT-134
      IJ=I-JN1+1                                                        SCAT-135
      W5=1.D0                                                           SCAT-136
      JN=IPJ                                                            SCAT-137
      IF (.NOT.LO(41)) GO TO 9                                          SCAT-138
C FACTORISATION OF 1/(1-X*COS(THETA)) WITH X AFTER THE SCATTERING       SCAT-139
C COEFFICIENT.                                                          SCAT-140
      DO 8 J=1,NCJ                                                      SCAT-141
      JN=JN+1                                                           SCAT-142
    8 W5=W5/(1.D0-FR(1,JN,I)*W6)                                        SCAT-143
    9 IF (MF(6,I).EQ.99999) GO TO 10                                    SCAT-144
      CALL EMRO(JMIN,MF(5,I),MF(6,I),X2,X3,B,JN)                        SCAT-145
      W1=W5                                                             SCAT-146
      GO TO 11                                                          SCAT-147
C ROTATION MATRIX ELEMENTS ARE THE SAME AS FOR THE LAST AMPLITUDE.      SCAT-148
   10 W1=W5*DFLOAT(MF(10,I))                                            SCAT-149
   11 DO 12 J=1,IPJ                                                     SCAT-150
      A(1,IJ)=A(1,IJ)+FR(1,J,I)*W1*B(J)                                 SCAT-151
   12 A(2,IJ)=A(2,IJ)+FR(2,J,I)*W1*B(J)                                 SCAT-152
C CONSTRUCTION OF THE TOTAL AMPLITUDE MATRIX.                           SCAT-153
      K1=MF(7,I)                                                        SCAT-154
      AB(1,K1)=A(1,IJ)                                                  SCAT-155
      AB(2,K1)=A(2,IJ)                                                  SCAT-156
      IF (MF(8,I).EQ.0) GO TO 13                                        SCAT-157
      K1=MF(8,I)                                                        SCAT-158
      W1=DFLOAT(MF(9,I))                                                SCAT-159
      AB(1,K1)=W1*A(1,IJ)                                               SCAT-160
      AB(2,K1)=W1*A(2,IJ)                                               SCAT-161
   13 CONTINUE                                                          SCAT-162
C CROSS-SECTION.                                                        SCAT-163
      DZ=0.D0                                                           SCAT-164
      IF (LO(81)) DZ=EXCN/XZ                                            SCAT-165
      DO 14 I1=1,MTT                                                    SCAT-166
   14 DZ=DZ+AB(1,I1)*AB(1,I1)+AB(2,I1)*AB(2,I1)                         SCAT-167
      EX(1)=DZ*XZ                                                       SCAT-168
      LX1=0                                                             SCAT-169
      LX2=0                                                             SCAT-170
      JEX=1                                                             SCAT-171
C LOOP ON THE OBSERVABLES.                                              SCAT-172
      DO 54 IEX=NIX,NFX                                                 SCAT-173
      JEX=JEX+1                                                         SCAT-174
      IF (MF(2,IEX).GT.1) GO TO 21                                      SCAT-175
      IF (MF(2,IEX)) 17 , 15 , 16                                       SCAT-176
C CROSS SECTION.                                                        SCAT-177
   15 EX(JEX)=EX(1)                                                     SCAT-178
      GO TO 54                                                          SCAT-179
C CROSS SECTION DIVIDED BY RUTHERFORD'S CROSS-SECTION.                  SCAT-180
   16 EX(JEX)=EX(1)/(10.D0*W3)                                          SCAT-181
      GO TO 54                                                          SCAT-182
   17 EX(JEX)=0.D0                                                      SCAT-183
      IF (MF(2,IEX).EQ.-3) GO TO 19                                     SCAT-184
      MTP=MT1*MT2*MT4                                                   SCAT-185
C VECTOR ANALYSING POWER FOR SPIN 1/2 OR 1.                             SCAT-186
      DO 18 I1=1,MTP                                                    SCAT-187
      I2=I1+MTP                                                         SCAT-188
   18 EX(JEX)=EX(JEX)+AB(2,I2)*AB(1,I1)-AB(1,I2)*AB(2,I1)               SCAT-189
      EX(JEX)=2.D0*EX(JEX)*XZ/EX(1)                                     SCAT-190
      IF (MT3.EQ.3) EX(JEX)=1.22474487D0*EX(JEX)                        SCAT-191
      GO TO 54                                                          SCAT-192
C VECTOR POLARISATION FOR SPIN 1/2 OR 1.                                SCAT-193
   19 DO 20 I1=1,MTT,MT1                                                SCAT-194
   20 EX(JEX)=EX(JEX)+AB(1,I1+1)*AB(2,I1)-AB(2,I1+1)*AB(1,I1)           SCAT-195
      EX(JEX)=2.D0*EX(JEX)*XZ/EX(1)                                     SCAT-196
      IF (MT1.EQ.3) EX(JEX)=1.22474487D0*EX(JEX)                        SCAT-197
      GO TO 54                                                          SCAT-198
C ALL THE OTHER OBSERVABLES.                                            SCAT-199
   21 EX(JEX)=0.D0                                                      SCAT-200
      K1=MF(3,IEX)                                                      SCAT-201
      K2=MF(4,IEX)                                                      SCAT-202
C LOOP ON THE COMPONENTS OF THE DESCRIPTION OF THE OBSERVABLE.          SCAT-203
      DO 53 II=K1,K2                                                    SCAT-204
      IX=1                                                              SCAT-205
      IXY=0                                                             SCAT-206
      IXX=IXY+MTT                                                       SCAT-207
      DO 22 I=1,18                                                      SCAT-208
   22 MO(I)=NCO(I,II)                                                   SCAT-209
      DO 23 IJ=1,6                                                      SCAT-210
      LTT(IJ)=MOD(JQ,2).EQ.1                                            SCAT-211
   23 JQ=JQ/2                                                           SCAT-212
C COMPUTATION OF THE ANGLE FOR A CHANGE OF FRAME.                       SCAT-213
      IF (MOD(JP,1000).EQ.0) GO TO 46                                   SCAT-214
      LXY=.TRUE.                                                        SCAT-215
      IF (LX1.EQ.MOD(JP,1000)) GO TO 45                                 SCAT-216
      LX1=MOD(JP,1000)                                                  SCAT-217
      LX2=0                                                             SCAT-218
      N=MT1                                                             SCAT-219
      IF (LX1.EQ.1) GO TO 26                                            SCAT-220
      GO TO 25                                                          SCAT-221
   24 IF (JP/1000.EQ.0) GO TO 46                                        SCAT-222
      LXY=.FALSE.                                                       SCAT-223
      IF (JP/1000.EQ.LX2) GO TO 45                                      SCAT-224
      LX2=JP/1000                                                       SCAT-225
      N=MT2                                                             SCAT-226
      IF (LX2.EQ.1) GO TO 27                                            SCAT-227
   25 IF (DABS(W7).LT.1.D-5.OR.N.EQ.1) GO TO 43                         SCAT-228
      W1=W6                                                             SCAT-229
      W2=-W7                                                            SCAT-230
      GO TO 31                                                          SCAT-231
C TRANSFORMATION TO THE LABORATORY SYSTEM FOR THE PARTICLE.             SCAT-232
   26 W1=WV(1,NOUT)                                                     SCAT-233
      GO TO 28                                                          SCAT-234
C TRANSFORMATION TO THE LABORATORY SYSTEM FOR THE TARGET.               SCAT-235
   27 W1=-WV(2,NOUT)                                                    SCAT-236
   28 IF (N.EQ.1) GO TO 43                                              SCAT-237
      IF (.NOT.LO(08)) GO TO 29                                         SCAT-238
      W1=WV(4,NOUT)/(CMB*W1)                                            SCAT-239
      W2=CMB*WV(2,1)/WV(4,1)                                            SCAT-240
      W3=DATAN(-W7/(W6*DSQRT(W1*W1+1.D0)+W1*DSQRT(W2*W2+1.D0)))         SCAT-241
      GO TO 30                                                          SCAT-242
   29 W3=DATAN(-W7/(WV(4,NOUT)*WV(2,1)/(WV(4,1)*W1)+W6))                SCAT-243
   30 IF (DABS(W3).LT.1.D-5.OR.N.EQ.1) GO TO 43                         SCAT-244
      W1=DCOS(W3)                                                       SCAT-245
      W2=DSIN(W3)                                                       SCAT-246
   31 DO 32 I=1,MTT                                                     SCAT-247
      AB(1,I+IXX)=0.D0                                                  SCAT-248
   32 AB(2,I+IXX)=0.D0                                                  SCAT-249
      W3=DSQRT(.5D0*(1.D0-W1))                                          SCAT-250
      W5=1.D0                                                           SCAT-251
      DO 33 I=2,N                                                       SCAT-252
   33 W5=-W5*W3                                                         SCAT-253
      IF (DABS(W5).LT.1.D-30) GO TO 43                                  SCAT-254
      W4=0.D0                                                           SCAT-255
      W8=.5D0*DFLOAT(N-1)                                               SCAT-256
      MTY=MT3*MT4                                                       SCAT-257
      X3=-W8                                                            SCAT-258
C TRANSFORMATION OF THE AMPLITUDE MATRIX.                               SCAT-259
      DO 42 I=1,N                                                       SCAT-260
      IF (I.EQ.1) GO TO 34                                              SCAT-261
      W3=W4                                                             SCAT-262
      IF (I.NE.2) W3=W3*DSQRT(DFLOAT((I-2)*(N+2-I)))                    SCAT-263
      W4=W5                                                             SCAT-264
      W5=(2.D0*(X3*W1-W8)*W4/W2-W3)/DSQRT(DFLOAT((I-1)*(N+1-I)))        SCAT-265
      X3=X3+1.D0                                                        SCAT-266
   34 X2=0.D0                                                           SCAT-267
      X5=W5                                                             SCAT-268
      X4=W8                                                             SCAT-269
      DO 41 J=I,N                                                       SCAT-270
      IF (J.EQ.I) GO TO 35                                              SCAT-271
      X1=X2                                                             SCAT-272
      IF (J.NE.I+1) X1=X1*DSQRT(DFLOAT((J-I-1)*(N+1-J+I)))              SCAT-273
      X2=X5                                                             SCAT-274
      X5=(2.D0*(X3-X4*W1)*X2/W2-X1)/DSQRT(DFLOAT((J-I)*(N+I-J)))        SCAT-275
      X4=X4-1.D0                                                        SCAT-276
   35 X=X5                                                              SCAT-277
      L=1+J-I                                                           SCAT-278
      M=N+1-I                                                           SCAT-279
      DO 40 K=1,2                                                       SCAT-280
      DO 39 I1=1,MTY                                                    SCAT-281
      IF (LXY) GO TO 37                                                 SCAT-282
      IXI=IXX+MT1*((I1-1)*MT2+L-1)                                      SCAT-283
      IXJ=IXY+MT1*((I1-1)*MT2+M-1)                                      SCAT-284
      DO 36 I4=1,MT1                                                    SCAT-285
      AB(1,I4+IXI)=AB(1,I4+IXI)+X*AB(1,I4+IXJ)                          SCAT-286
   36 AB(2,I4+IXI)=AB(2,I4+IXI)+X*AB(2,I4+IXJ)                          SCAT-287
      GO TO 39                                                          SCAT-288
   37 IXI=IXX+MT1*((I1-1)*MT2-1)+L                                      SCAT-289
      IXJ=IXY+MT1*((I1-1)*MT2-1)+M                                      SCAT-290
      DO 38 I4=MT1,MT2*MT1,MT1                                          SCAT-291
      AB(1,IXI+I4)=AB(1,IXI+I4)+X*AB(1,IXJ+I4)                          SCAT-292
   38 AB(2,IXI+I4)=AB(2,IXI+I4)+X*AB(2,IXJ+I4)                          SCAT-293
   39 CONTINUE                                                          SCAT-294
      IF (J.EQ.N) GO TO 41                                              SCAT-295
      ML=L                                                              SCAT-296
      L=M                                                               SCAT-297
      M=ML                                                              SCAT-298
      IF (MOD(J+N,2).NE.0) X=-X                                         SCAT-299
   40 CONTINUE                                                          SCAT-300
   41 CONTINUE                                                          SCAT-301
   42 CONTINUE                                                          SCAT-302
      GO TO 45                                                          SCAT-303
   43 DO 44 I=1,MTT                                                     SCAT-304
      AB(1,I+IXX)=AB(1,I+IXY)                                           SCAT-305
   44 AB(2,I+IXX)=AB(2,I+IXY)                                           SCAT-306
   45 IX=IX+1                                                           SCAT-307
      IXY=IXX                                                           SCAT-308
      IXX=IXY+MTT                                                       SCAT-309
      IF (LXY) GO TO 24                                                 SCAT-310
   46 IF (LTT(5)) GO TO 52                                              SCAT-311
      X=0.D0                                                            SCAT-312
C THE FOUR DO LOOPS.                                                    SCAT-313
      DO 51 I1=MI1,MP1                                                  SCAT-314
      J1=I1-N1                                                          SCAT-315
      X1=COE(10*II)                                                     SCAT-316
      IF (LTT(1)) X1=COE(10*II)*COE(I1+L1)                              SCAT-317
      DO 50 I2=MI2,MP2                                                  SCAT-318
      J2=I2-N2                                                          SCAT-319
      X2=X1                                                             SCAT-320
      IF (LTT(2)) X2=X1*COE(L2-I2)                                      SCAT-321
      DO 49 I3=MI3,MP3                                                  SCAT-322
      J3=I3-N3                                                          SCAT-323
      X3=X2                                                             SCAT-324
      IF (LTT(3)) X3=X2*COE(I3+L3)                                      SCAT-325
      DO 48 I4=MI4,MP4                                                  SCAT-326
      J4=I4-N4                                                          SCAT-327
      X4=X3                                                             SCAT-328
      IIJ=MT1*(I4+MT2*(J2-1+MT4*(J1-1+MT3*(IX-1)))-1)                   SCAT-329
      III=MT1*(J4+MT2*(I2-1+MT4*(I1-1+MT3*(IX-1)))-1)                   SCAT-330
      IF (LTT(4)) X4=X3*COE(L4-I4)                                      SCAT-331
      IF (LTT(6)) GO TO 47                                              SCAT-332
      X=X+(AB(1,I3+IIJ)*AB(1,J3+III)+AB(2,I3+IIJ)*AB(2,J3+III))*X4      SCAT-333
      GO TO 48                                                          SCAT-334
   47 X=X+(AB(2,I3+IIJ)*AB(1,J3+III)-AB(1,I3+IIJ)*AB(2,J3+III))*X4      SCAT-335
   48 CONTINUE                                                          SCAT-336
   49 CONTINUE                                                          SCAT-337
   50 CONTINUE                                                          SCAT-338
   51 CONTINUE                                                          SCAT-339
      GO TO 53                                                          SCAT-340
   52 X=COE(10*II)*DZ                                                   SCAT-341
   53 EX(JEX)=EX(JEX)+X                                                 SCAT-342
      EX(JEX)=EX(JEX)/DZ                                                SCAT-343
   54 CONTINUE                                                          SCAT-344
      IF (LKT) GO TO 55                                                 SCAT-345
      XC=DSQRT((XA*W6+XB)**2+W7**2)**3/DABS(XA+XB*W6)                   SCAT-346
      EX(1)=EX(1)*XC                                                    SCAT-347
   55 IF (KLT.EQ.1.OR.LKT) RETURN                                       SCAT-348
      EX(2)=EX(2)*XC                                                    SCAT-349
      IF (LO(81)) EXCN=EXCN*XC                                          SCAT-350
      RETURN                                                            SCAT-351
C PURE COMPOUND NUCLEUS FOR LEVELS JN1 TO JN2.                          SCAT-352
   56 IS=JN1-1                                                          SCAT-353
      DO 63 IV=JN1,JN2                                                  SCAT-354
      EX(IV-IS)=FCN(1,IV)                                               SCAT-355
      IF (IPK.EQ.1) GO TO 62                                            SCAT-356
      IF (LKT) GO TO 59                                                 SCAT-357
      IF (LO(8)) GO TO 57                                               SCAT-358
      XA=1.D0                                                           SCAT-359
      XB=WV(4,1)*WV(1,IV)/(WV(4,IV)*WV(2,1))                            SCAT-360
      GO TO 58                                                          SCAT-361
   57 XA=DSQRT(1.D0+(WV(4,1)/(CMB*WV(2,1)))**2)                         SCAT-362
      XB=DSQRT(WV(1,IV)**2+(WV(4,IV)/CMB)**2)*WV(4,1)/(WV(4,IV)*WV(2,1))SCAT-363
   58 W6=DCOS(THET)                                                     SCAT-364
      W7=DSIN(THET)                                                     SCAT-365
      XC=(XA**2-XB**2)*W7**2+W6**2                                      SCAT-366
      IF (XC.GT.0.D0) XC=DSQRT(XC)                                      SCAT-367
      W6=XC*W6-XA*XB*W7**2                                              SCAT-368
      GO TO 60                                                          SCAT-369
   59 W6=DCOS(THET)                                                     SCAT-370
   60 U1=0.D0                                                           SCAT-371
      U2=1.D0                                                           SCAT-372
      DO 61 J=2,IPK                                                     SCAT-373
      U1=(DFLOAT(4*J-7)*U2*W6-DFLOAT(2*J-4)*U1)/DFLOAT(2*J-3)           SCAT-374
      U2=(DFLOAT(4*J-5)*U1*W6-DFLOAT(2*J-3)*U2)/DFLOAT(2*J-2)           SCAT-375
   61 EX(IV-IS)=EX(IV-IS)+U2*FCN(J,IV)*DFLOAT(4*J-3)                    SCAT-376
   62 IF (.NOT.LKT) EX(IV-IS)=EX(IV-IS)*DSQRT((XA*W6+XB)**2+1.D0-W6**2)*SCAT-377
     1*3/DABS(XA+XB*W6)                                                 SCAT-378
   63 CONTINUE                                                          SCAT-379
      RETURN                                                            SCAT-380
      END                                                               SCAT-381
