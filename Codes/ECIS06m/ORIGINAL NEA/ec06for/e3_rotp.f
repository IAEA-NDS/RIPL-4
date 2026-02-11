C 07/03/07                                                      ECIS06  ROTP-000
      SUBROUTINE ROTP(BETA,NBETA,IP,IVY,VAL,ID1,ISM,WV,INVZ,IQ1,INY,PGN,ROTP-001
     1XGN,LO,V,VA,IV,PP,B,P)                                            ROTP-002
C INPUT:     BETA:    VIBRATIONAL DEFORMATIONS. FROM BETA(1,NBT1+1),    ROTP-003
C                     STATIC DEFORMATIONS (ROTATIONAL MODEL).           ROTP-004
C            NBETA:   QUANTUM NUMBERS OF DEFORMATIONS (IN EQUIVALENCE   ROTP-005
C                     BY CALL WITH BETA).                               ROTP-006
C            IP:      POTENTIAL NUMBER (SCHROEDINGER EQUATION), ADDRESS ROTP-007
C                     TO TEMPORARILY STORAGE OF DIRAC POTENTIALS IN V.  ROTP-008
C            IVY:     DESCRIPTION OF FORM FACTORS                       ROTP-009
C                     (SEE THIRD PART OF IQ IN REDM).                   ROTP-010
C            VAL:     OPTICAL PARAMETERS..                              ROTP-011
C            ID1:     FIRST DIMENSION OF THE WORKING SPACE P.           ROTP-012
C            ISM:     NUMBER OF STEPS FOR INTEGRATION.                  ROTP-013
C            WV:      INFORMATIONS FOR ROTZ, STEP SIZE IN WV(8).        ROTP-014
C            INVZ:    NUMBER OF TRANSITION FORM FACTORS.                ROTP-015
C            IQ1:     INY + NUMBER OF TRANSITION FORM FACTORS.          ROTP-016
C            INY:     1 + NUMBER OF DERIVATIVES OF CENTRAL              ROTP-017
C                     POTENTIALS ( SCHROEDINGER: 1, DIRAC: 2)           ROTP-018
C            PGN,XGN: WEIGHTS, ABSCISSAE OF LEGENDRE INTEGRAL.          ROTP-019
C            LO(I):   LOGICAL CONTROLS:                                 ROTP-020
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).ROTP-021
C               LO(3)  =.TRUE. ANHARMONIC VIBRATIONAL OR ASYMMETRIC     ROTP-022
C                              ROTATIONAL MODEL.                        ROTP-023
C               LO(6)  =.TRUE. USE DEFORMATION LENGTHS.                 ROTP-024
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    ROTP-025
C               LO(9)  =.TRUE. SYMMETRISED WOODS-SAXON FORM FACTORS WHENROTP-026
C                              THE RADIUS IS NEGATIVE.                  ROTP-027
C               LO(11) =.TRUE. DEFORMED COULOMB POTENTIAL.              ROTP-028
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            ROTP-029
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      ROTP-030
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. ROTP-031
C               LO(17) =.TRUE. FOLDING MODEL.                           ROTP-032
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   ROTP-033
C               LO(99) =.TRUE. SCHROEDINGER EQUIVALENT TO DIRAC         ROTP-034
C                              EQUATION.                                ROTP-035
C               LO(100)=.TRUE. DIRAC EQUATION.                          ROTP-036
C               LO(101)=.TRUE. THERE IS A REAL SPIN-ORBIT POTENTIAL.    ROTP-037
C               LO(102)=.TRUE. THERE IS AN IMAGINARY SPIN-ORBIT         ROTP-038
C                              POTENTIAL.                               ROTP-039
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. ROTP-040
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    ROTP-041
C OUTPUT:    V:       POTENTIALS AND TRANSITION FORM FACTORS,           ROTP-042
C                     ONLY POTENTIALS AND DERIVATIVES FOR DIRAC.        ROTP-043
C            VA:      TRANSITION FORM FACTORS FOR DIRAC EQUATION.       ROTP-044
C WORKING AREAS:                                                        ROTP-045
C            IV:      ORDER OF DERIVATIVE OF THE IQ1 FORM FACTORS.      ROTP-046
C            PP:      INTERMEDIATE RESULTS, WEIGHTS, .... IN ROTD.      ROTP-047
C            B:       FOR DEFORMATIONS.                                 ROTP-048
C            P:       WEIGHTS COMPUTED IN ROTD.                         ROTP-049
C                                                                       ROTP-050
C FOR THE COMMONS /COUPL/ AND /DCONS/ SEE CALC.                         ROTP-051
C FOR THE COMMON  /POTE1/ SEE REDM.                                     ROTP-052
C                                                                       ROTP-053
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     ROTP-054
C  IQM:       MAXIMUM L-VALUE OF DEFORMATION IN ROTATIONAL MODEL.       ROTP-055
C  NBT1:      NUMBER OF PHONONS (VIBRATIONS).                           ROTP-056
C   USED:     IQM,NBT1.                                                 ROTP-057
C                                                                       ROTP-058
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     ROTP-059
C  CCZ:       COULOMB ALPHA CONSTANT.                                   ROTP-060
C   USED:     CCZ.                                                      ROTP-061
C                                                                       ROTP-062
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE1/:                     ROTP-063
C  ITX(16):   STARTING ADDRESS OF DIFFERENT FORM FACTORS (SEE REDM).    ROTP-064
C  IMAX:      MAXIMUM ANGULAR MOMENTUM.                                 ROTP-065
C  INLS:      NUMBER OF SPIN-ORBIT FORM FACTORS NOT TAKING INTO ACCOUNT ROTP-066
C             MULTIPLICATION BY 2.                                      ROTP-067
C  INVD:      IDEM FOR COULOMB SPIN-ORBIT.                              ROTP-068
C   USED:     ITX,IMAX,INLS,INVD.                                       ROTP-069
C                                                                       ROTP-070
C***********************************************************************ROTP-071
      IMPLICIT REAL*8 (A-H,O-Z)                                         ROTP-072
      LOGICAL LO(150),LQ(8,5)                                           ROTP-073
      DIMENSION BETA(9,*),NBETA(18,*),V(ISM,*),VA(ISM,11,*),IVY(7,*),IV(ROTP-074
     1*),PP(*),B(10,*),P(ID1,*),VAL(4,9),PGN(10),XGN(10),IZ(8,2),EP(8),ZROTP-075
     2Z(2),DD(2),Q(8,36),VR(6,10),LDL(8),SRD(4),WV(22),PW(8)            ROTP-076
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   ROTP-077
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            ROTP-078
      COMMON /POTE1/ ITX(16),IMAX,INTC,INLS,INVC,INVD,ITXM              ROTP-079
      COMMON /INOUT/ MR,MW,MS                                           ROTP-080
      DATA LDL /7,12,7,12,13,14,11,19/                                  ROTP-081
      DATA SRD /1.D0,1.D0,2.D0,6.D0/                                    ROTP-082
C CHECKS OF DIFFUSENESSES.                                              ROTP-083
      DO 2 I=1,8                                                        ROTP-084
      J=LDL(I)                                                          ROTP-085
      LQ(I,5)=.NOT.LO(J)                                                ROTP-086
      PW(I)=1.D0+VAL(4,I)                                               ROTP-087
      IF (J.EQ.7) LQ(I,5)=LO(J)                                         ROTP-088
      LQ(I,1)=VAL(1,I).EQ.0.D0                                          ROTP-089
      LQ(I,2)=(I.LE.6).OR.(VAL(1,I)*VAL(3,I).NE.0.D0)                   ROTP-090
      LQ(I,3)=LQ(I,1).OR.(.NOT.LQ(I,2))                                 ROTP-091
      LQ(I,2)=LQ(I,2).OR.LQ(I,1)                                        ROTP-092
      IF (LQ(I,1)) GO TO 2                                              ROTP-093
      LQ(I,4)=(.NOT.LO(9)).OR.(VAL(2,I).GE.0.D0)                        ROTP-094
      IF (.NOT.LQ(I,2)) GO TO 1                                         ROTP-095
      IF (VAL(3,I).GE.0.02D0*WV(8)) GO TO 1                             ROTP-096
      WRITE (MW,1000) VAL(3,I),I,IP                                     ROTP-097
      VAL(3,I)=DMAX1(-VAL(3,I),0.02D0*WV(8))                            ROTP-098
    1 IF (VAL(2,I).GE.WV(8).OR.LO(9)) GO TO 2                           ROTP-099
      WRITE (MW,1001) VAL(2,I),I,IP                                     ROTP-100
      VAL(2,I)=DMAX1(-VAL(2,I),WV(8))                                   ROTP-101
    2 CONTINUE                                                          ROTP-102
      LQ(1,2)=LQ(7,2).AND.LQ(8,2)                                       ROTP-103
      IDG=0                                                             ROTP-104
      INX=INY+1                                                         ROTP-105
      SR=1.D0                                                           ROTP-106
      DO 11 I=1,IQ1                                                     ROTP-107
      IV(I)=I+3                                                         ROTP-108
      IF (I.LE.INY) GO TO 8                                             ROTP-109
      K=MOD(IVY(1,I-INY),1000)                                          ROTP-110
      IF (K.EQ.0.OR.LO(3)) GO TO 6                                      ROTP-111
      IF (K.GT.NBT1) GO TO 4                                            ROTP-112
      IV(I)=2                                                           ROTP-113
      IDG=MAX0(IDG,1)                                                   ROTP-114
      DO 3 J=1,8                                                        ROTP-115
      IF (.NOT.LO(6)) SR=DABS(VAL(2,J))                                 ROTP-116
    3 B(J,I)=0.282095D0*BETA(J,K)*SR                                    ROTP-117
      GO TO 10                                                          ROTP-118
    4 K1=MOD(K,NBT1+1)                                                  ROTP-119
      K2=K/(NBT1+1)                                                     ROTP-120
      IDG=MAX0(IDG,2)                                                   ROTP-121
      IV(I)=3                                                           ROTP-122
      DO 5 J=1,8                                                        ROTP-123
      IF (.NOT.LO(6)) SR=DABS(VAL(2,J))                                 ROTP-124
    5 B(J,I)=0.0397887D0*BETA(J,K1)*BETA(J,K2)*SR*SR                    ROTP-125
      GO TO 10                                                          ROTP-126
    6 IV(I)=K+1                                                         ROTP-127
      IDG=MAX0(IDG,K)                                                   ROTP-128
      IF (LO(1)) GO TO 8                                                ROTP-129
      DO 7 J=1,8                                                        ROTP-130
      IF (.NOT.LO(6)) SR=DABS(VAL(2,J))                                 ROTP-131
    7 B(J,I)=SR**K*0.282095D0*BETA(J,K+1)/SRD(K+1)                      ROTP-132
      GO TO 10                                                          ROTP-133
    8 DO 9 J=1,8                                                        ROTP-134
    9 B(J,I)=1.D0                                                       ROTP-135
   10 B(9,I)=B(5,I)                                                     ROTP-136
   11 B(10,I)=B(6,I)                                                    ROTP-137
      IV(1)=1                                                           ROTP-138
C SET UP OF FORM FACTOR COMPUTATION.                                    ROTP-139
      IDS=1                                                             ROTP-140
      IDR=1                                                             ROTP-141
      IF (LO(17).OR.LO(99)) IDS=0                                       ROTP-142
      IF (LO(109)) IDR=0                                                ROTP-143
      IDT=IDR+IDS                                                       ROTP-144
      DO 12 I=1,8                                                       ROTP-145
      IDZ=IDG                                                           ROTP-146
      IF (LQ(I,5)) IDZ=0                                                ROTP-147
      IZ(I,2)=IDZ+1                                                     ROTP-148
      IF (LQ(I,3)) GO TO 12                                             ROTP-149
      EP(I)=DEXP(WV(8)/VAL(3,I))                                        ROTP-150
      IF ((I.EQ.3).OR.(I.EQ.4)) IDZ=IDZ+IDR                             ROTP-151
      IF ((I.EQ.5).OR.(I.EQ.6).OR.(I.EQ.8)) IDZ=IDZ+IDS                 ROTP-152
      IF (LO(109)) IDZ=MAX0(IDZ,2)                                      ROTP-153
   12 IZ(I,1)=IDZ                                                       ROTP-154
C LQ(I,1) IS .TRUE. IF THE FORM FACTOR IS NOT USED.                     ROTP-155
C LQ(I,2) IS .TRUE. FOR COULOMB FORM FACTOR WITH DIFFUSENESS.           ROTP-156
C LQ(I,3) IS .FALSE. FOR ANY WOODS-SAXON FORM FACTORS.                  ROTP-157
C LQ(I,4) IS .FALSE. FOR SYMMETRISED WOODS-SAXON FORM FACTORS.          ROTP-158
C LQ(I,5) IS .FALSE. FOR DEFORMED FORM FACTORS.                         ROTP-159
      CALL ROTD(NBETA,BETA(1,NBT1+1),IVY,IV,PP,P,Q,VAL,ID1,INVZ,IMAX,IQ1ROTP-160
     1,IQM,INY,IDT,LQ,DD,PGN,XGN,LO)                                    ROTP-161
      ZZ(1)=CCZ*VAL(1,7)                                                ROTP-162
      ZZ(2)=CCZ*VAL(1,8)                                                ROTP-163
      IF (DD(1).NE.0.D0) ZZ(1)=ZZ(1)/DD(1)                              ROTP-164
      IF (DD(2).NE.0.D0) ZZ(2)=ZZ(2)/DD(2)                              ROTP-165
      DO 14 J=1,10                                                      ROTP-166
      DO 13 I=1,6                                                       ROTP-167
   13 VR(I,J)=0.D0                                                      ROTP-168
   14 CONTINUE                                                          ROTP-169
      R=0.D0                                                            ROTP-170
      DO 52 IS=1,ISM                                                    ROTP-171
      R=R+WV(8)                                                         ROTP-172
      DO 16 I=1,IQ1                                                     ROTP-173
      DO 15 J=17,27                                                     ROTP-174
   15 P(J,I)=0.D0                                                       ROTP-175
   16 CONTINUE                                                          ROTP-176
C INTEGRATION LOOP.                                                     ROTP-177
      DO 43 J=1,IDT                                                     ROTP-178
      DO 26 I=1,8                                                       ROTP-179
      IF (LQ(I,3)) GO TO 26                                             ROTP-180
      Q(I,J)=Q(I,J)+WV(8)                                               ROTP-181
      IF (P(I,J).NE.0.D0) GO TO 17                                      ROTP-182
      IF (Q(I,J)+50.D0*VAL(3,I).GT.0.D0) P(I,J)=DEXP(Q(I,J)/VAL(3,I))   ROTP-183
      GO TO 18                                                          ROTP-184
   17 IF (P(I,J).LT.1.D15) P(I,J)=P(I,J)*EP(I)                          ROTP-185
   18 IF ((.NOT.LQ(I,4)).AND.P(I+8,J).GT.1.D-15) P(I+8,J)=P(I+8,J)/EP(I)ROTP-186
      IDV=IZ(I,1)                                                       ROTP-187
      CALL WOSA(VR(1,I),VAL(3,I),VAL(4,I),P(I,J),Q(I,J),IDV+2,LQ(I,4))  ROTP-188
      DO 19 K=1,IDV+2                                                   ROTP-189
   19 VR(K,I)=VR(K,I)*VAL(1,I)                                          ROTP-190
      IF ((I.LT.3).OR.(I.GE.7).OR.(IDT.EQ.0)) GO TO 25                  ROTP-191
      IF (I.GT.4) GO TO 21                                              ROTP-192
      IF (IDR.EQ.0) GO TO 25                                            ROTP-193
      DO 20 K=1,IDV                                                     ROTP-194
   20 VR(K,I)=4.D0*VR(K+1,I)*VAL(3,I)                                   ROTP-195
      GO TO 25                                                          ROTP-196
   21 IF (IDS.EQ.0) GO TO 25                                            ROTP-197
      IF (LO(100)) GO TO 23                                             ROTP-198
      DO 22 K=1,IDV                                                     ROTP-199
      VR(K,I+4)=VR(K,I)/R**2                                            ROTP-200
   22 VR(K,I)=VR(K+1,I)/R                                               ROTP-201
      GO TO 25                                                          ROTP-202
   23 DO 24 K=1,IDV                                                     ROTP-203
   24 VR(K,I+4)=VR(K+1,I)                                               ROTP-204
   25 IF (.NOT.LO(109)) GO TO 26                                        ROTP-205
      VR(5,I)=VR(2,I)                                                   ROTP-206
      VR(6,I)=VR(3,I)                                                   ROTP-207
   26 CONTINUE                                                          ROTP-208
      IF (.NOT.LO(99)) GO TO 27                                         ROTP-209
      CALL ROTZ(VR,XGN,WV,R,Q,VAL,J,IV,B,P,ID1,IQ1,LO)                  ROTP-210
      GO TO 43                                                          ROTP-211
   27 DO 29 L=1,IQ1                                                     ROTP-212
      K=IV(L)                                                           ROTP-213
      DO 28 M=1,10                                                      ROTP-214
   28 P(M+16,L)=P(M+16,L)+VR(K,M)*P(L+27,J)*B(M,L)                      ROTP-215
   29 CONTINUE                                                          ROTP-216
      IF (LQ(1,2)) GO TO 43                                             ROTP-217
C DEFORMED COULOMB POTENTIAL.                                           ROTP-218
      DO 42 I=7,8                                                       ROTP-219
      IF (LQ(I,2)) GO TO 42                                             ROTP-220
      C=P(I,J)/R                                                        ROTP-221
      IF (I.EQ.8.AND.(.NOT.(LO(17).OR.LO(100)))) GO TO 31               ROTP-222
      IF (R.LT.P(I,J)) GO TO 30                                         ROTP-223
      P(I+16,1)=P(I+16,1)+ZZ(I-6)*(P(I,J)**2)*C*P(28,J)                 ROTP-224
      IF (.NOT.LO(100)) GO TO 32                                        ROTP-225
      P(I+16,2)=P(I+16,2)+ZZ(I-6)*P(I,J)*C**2*P(28,J)                   ROTP-226
      P(I+16,3)=P(I+16,3)+ZZ(I-6)*2.D0*C**3*P(28,J)                     ROTP-227
      GO TO 32                                                          ROTP-228
   30 P(I+16,1)=P(I+16,1)+(0.5D0*P(I,J)*P(I,J)-R*R/6.D0)*ZZ(I-6)*P(28,J)ROTP-229
     1*3.D0                                                             ROTP-230
      IF (.NOT.LO(100)) GO TO 32                                        ROTP-231
      P(I+16,2)=P(I+16,2)+R*ZZ(I-6)*P(28,J)                             ROTP-232
      P(I+16,3)=P(I+16,3)-ZZ(I-6)*P(28,J)                               ROTP-233
      GO TO 32                                                          ROTP-234
   31 C1=ZZ(I-6)*P(28,J)                                                ROTP-235
      IF (R.GT.P(I,J)) C1=C1*C**3                                       ROTP-236
      P(24,1)=P(24,1)+C1                                                ROTP-237
   32 IF (INVZ.EQ.0) GO TO 42                                           ROTP-238
      IF (.NOT.LO(8*I-45)) GO TO 42                                     ROTP-239
      DO 41 K=INX,IQ1                                                   ROTP-240
      L=IVY(7,K-INY)+1                                                  ROTP-241
      IF ((IVY(I-3,K-INY).EQ.0).OR.(L.EQ.1)) GO TO 41                   ROTP-242
      D=L                                                               ROTP-243
      N=IV(K)                                                           ROTP-244
      IF (R.LT.P(I,J)) GO TO 33                                         ROTP-245
      C1=(P(I,J)**2)*(C**L)*3.D0/((D+2.D0)*(2.D0*D-1.D0))               ROTP-246
      IF (N.GT.1) C1=C1*(D+2.D0)/P(I,J)                                 ROTP-247
      IF (N.GT.2) C1=C1*(D+1.D0)/P(I,J)                                 ROTP-248
      IF (N.GT.3) C1=C1*D/P(I,J)                                        ROTP-249
      GO TO 35                                                          ROTP-250
   33 IF (L.NE.3) GO TO 34                                              ROTP-251
      IF (N.EQ.1) C1=R*R*(0.2D0+DLOG(C))*0.6D0                          ROTP-252
      IF (N.GE.2) C1=0.6D0*R/C                                          ROTP-253
      IF (N.GE.3) C1=-C1/P(I,J)                                         ROTP-254
      IF (N.GE.4) C1=-2.D0*C1/P(I,J)                                    ROTP-255
      GO TO 35                                                          ROTP-256
   34 IF (N.EQ.1) C1=R*R*(1.D0/(D+2.D0)-1.D0/(C**(L-3)*(2.D0*D-1.D0)))*3ROTP-257
     1.D0/(D-3.D0)                                                      ROTP-258
      IF (N.GE.2) C1=R/C**(L-2)*3.D0/(2.D0*D-1.D0)                      ROTP-259
      IF (N.GE.3) C1=-C1*(D-2.D0)/P(I,J)                                ROTP-260
      IF (N.EQ.4) C1=-C1*(D-1.D0)/P(I,J)                                ROTP-261
   35 IF (I.EQ.7.OR.LO(17)) GO TO 40                                    ROTP-262
      IF (R.LT.P(I,J)) GO TO 36                                         ROTP-263
      C2=-D*C1/R                                                        ROTP-264
      GO TO 38                                                          ROTP-265
   36 IF (L.NE.3) GO TO 37                                              ROTP-266
      IF (N.EQ.1) C2=-1.2D0*(0.3D0*R-DLOG(C))*R                         ROTP-267
      IF (N.GE.2) C2=2.D0*C1/R                                          ROTP-268
      GO TO 38                                                          ROTP-269
   37 IF (N.EQ.1) C2=(2.D0/(D+2.D0)-(D-1.D0)/(C**(L-3)*(2.D0*D-1.D0)))*3ROTP-270
     1.D0/(D-3.D0)*R                                                    ROTP-271
      IF (N.NE.1)  C2=(D-1.D0)*C1/R                                     ROTP-272
   38 IF (LO(100)) GO TO 39                                             ROTP-273
      P(27,K)=P(27,K)-ZZ(2)*C1*P(K+27,J)*B(I,K)/R**2                    ROTP-274
      P(24,K)=P(24,K)-ZZ(2)*C2*P(K+27,J)*B(I,K)/R                       ROTP-275
      GO TO 41                                                          ROTP-276
   39 P(27,K)=P(27,K)-ZZ(2)*C2*P(K+27,J)*B(I,K)                         ROTP-277
      P(24,K)=P(24,K)+ZZ(2)*C1*P(K+27,J)*B(I,K)                         ROTP-278
      GO TO 41                                                          ROTP-279
   40 P(I+16,K)=P(I+16,K)+ZZ(I-6)*C1*P(K+27,J)*B(I,K)                   ROTP-280
   41 CONTINUE                                                          ROTP-281
   42 CONTINUE                                                          ROTP-282
   43 CONTINUE                                                          ROTP-283
C STORAGE OF FORM FACTORS.                                              ROTP-284
      IF (LO(100)) GO TO 47                                             ROTP-285
      V(IS,IP+ITX(1))=P(17,1)                                           ROTP-286
      V(IS,IP+ITX(2))=P(18,1)                                           ROTP-287
      V(IS,IP+ITX(3))=P(19,1)                                           ROTP-288
      V(IS,IP+ITX(4))=P(20,1)                                           ROTP-289
      IF (LO(101)) V(IS,IP+ITX(5))=P(21,1)                              ROTP-290
      IF (LO(102)) V(IS,IP+ITX(6))=P(22,1)                              ROTP-291
      V(IS,IP+ITX(7))=P(23,1)                                           ROTP-292
      IF (LO(103)) V(IS,IP+ITX(8))=P(24,1)                              ROTP-293
      IF (INVZ.LE.0) GO TO 52                                           ROTP-294
      DO 46 J=1,INVZ                                                    ROTP-295
      K=J+INY                                                           ROTP-296
      V(IS,J+ITX(9))=P(17,K)                                            ROTP-297
      V(IS,J+ITX(11))=P(19,K)                                           ROTP-298
      IF (.NOT.LO(12)) GO TO 44                                         ROTP-299
      V(IS,J+ITX(10))=P(18,K)                                           ROTP-300
      V(IS,J+ITX(12))=P(20,K)                                           ROTP-301
   44 IJ=IVY(4,J)                                                       ROTP-302
      IF (IJ.NE.0) V(IS,IJ+ITX(15))=P(23,K)                             ROTP-303
      IJ=IVY(3,J)                                                       ROTP-304
      IF (IJ.EQ.0) GO TO 45                                             ROTP-305
      V(IS,IJ+ITX(13))=P(21,K)                                          ROTP-306
      V(IS,IJ+ITX(13)+INLS)=-P(25,K)                                    ROTP-307
      IF (.NOT.LO(14)) GO TO 45                                         ROTP-308
      V(IS,IJ+ITX(14))=P(22,K)                                          ROTP-309
      V(IS,IJ+ITX(14)+INLS)=-P(26,K)                                    ROTP-310
   45 IJ=IVY(5,J)                                                       ROTP-311
      IF (IJ.EQ.0) GO TO 46                                             ROTP-312
      V(IS,IJ+ITX(16))=P(24,K)                                          ROTP-313
      V(IS,IJ+ITX(16)+INVD)=P(27,K)                                     ROTP-314
   46 CONTINUE                                                          ROTP-315
      GO TO 52                                                          ROTP-316
   47 K=IP                                                              ROTP-317
      DO 49 I=1,3                                                       ROTP-318
      DO 48 J=17,24                                                     ROTP-319
      V(IS,K)=P(J,I)                                                    ROTP-320
   48 K=K+1                                                             ROTP-321
   49 CONTINUE                                                          ROTP-322
      IF (INVZ.EQ.0) GO TO 52                                           ROTP-323
      DO 51 K=1,INVZ                                                    ROTP-324
      DO 50 J=1,11                                                      ROTP-325
   50 VA(IS,J,K)=P(J+16,K+3)                                            ROTP-326
   51 CONTINUE                                                          ROTP-327
   52 CONTINUE                                                          ROTP-328
C COULOMB POTENTIALS WITH DIFFUSE CHARGE DISTRIBUTION.                  ROTP-329
      IF (LQ(7,3).AND.LQ(8,3)) RETURN                                   ROTP-330
      DO 58 I=7,8                                                       ROTP-331
      IF (LQ(I,3)) GO TO 58                                             ROTP-332
      I1=IP+ITX(I)                                                      ROTP-333
      IF (I.EQ.8.AND.(.NOT.(LO(17).OR.LO(100)))) GO TO 53               ROTP-334
      IF (LO(100)) I1=IP+I-1                                            ROTP-335
      CALL COPO(V(1,I1),V(1,I1),P,ISM,WV(8),0,VAL(1,I),VAL(1,9),CCZ,ZT,.ROTP-336
     1FALSE.,.FALSE.)                                                   ROTP-337
      IF (.NOT.LO(100)) GO TO 54                                        ROTP-338
      CALL DERI(V(1,I1+8),V(1,I1),WV(8),ISM,.TRUE.)                     ROTP-339
      CALL DERI(V(1,I1+16),V(1,I1+8),WV(8),ISM,.TRUE.)                  ROTP-340
      GO TO 54                                                          ROTP-341
   53 CALL COPO(P,V(1,I1),P,ISM,WV(8),0,VAL(1,I),VAL(1,9),CCZ,ZT,.FALSE.ROTP-342
     1,.FALSE.)                                                         ROTP-343
      CALL DERI(V(1,I1),P,WV(8),ISM,.FALSE.)                            ROTP-344
   54 IF (INVZ.EQ.0.OR.(.NOT.LO(8*I-45))) GO TO 58                      ROTP-345
      DO 57 J=1,INVZ                                                    ROTP-346
      N=IVY(I-3,J)                                                      ROTP-347
      L=IVY(7,J)                                                        ROTP-348
      IF ((N.EQ.0).OR.(L.EQ.0)) GO TO 57                                ROTP-349
      I1=N+ITX(I+8)                                                     ROTP-350
      IF (LO(100)) I1=ITX(5)+11*J+I                                     ROTP-351
      I2=I1                                                             ROTP-352
      IF (.NOT.(LO(100).OR.(I.NE.8).OR.LO(17))) I2=I2+INVD              ROTP-353
      CALL COPO(V(1,I2),V(1,I1),P,ISM,WV(8),L,VAL(1,I),VAL(1,9),CCZ,ZT,.ROTP-354
     1FALSE.,.TRUE.)                                                    ROTP-355
      IF ((I.EQ.7).OR.LO(17)) GO TO 57                                  ROTP-356
      IF (.NOT.LO(100)) GO TO 55                                        ROTP-357
      CALL DERI(V(1,I1+3),V(1,I1),WV(8),ISM,.TRUE.)                     ROTP-358
      GO TO 57                                                          ROTP-359
   55 CALL DERI(V(1,I1),V(1,I2),WV(8),ISM,.FALSE.)                      ROTP-360
      RR=0.D0                                                           ROTP-361
      DO 56 IS=1,ISM                                                    ROTP-362
      RR=RR+WV(8)                                                       ROTP-363
   56 V(IS,I2)=-V(IS,I2)/RR**2                                          ROTP-364
   57 CONTINUE                                                          ROTP-365
   58 CONTINUE                                                          ROTP-366
      RETURN                                                            ROTP-367
 1000 FORMAT (' TOO SMALL DIFFUSENESS =',1P,D15.6,' FOR THE POTENTIAL',2ROTP-368
     1I4,'   CHANGED INTO ITS MINIMUM VALUE.')                          ROTP-369
 1001 FORMAT (' TOO SMALL COULOMB RADIUS =',1P,D15.6,' FOR THE POTENTIALROTP-370
     1',2I4,'   CHANGED INTO ITS MINIMUM VALUE.')                       ROTP-371
      END                                                               ROTP-372
