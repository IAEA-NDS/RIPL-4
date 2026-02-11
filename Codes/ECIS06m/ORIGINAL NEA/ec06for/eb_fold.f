C 20/08/06                                                      ECIS06  FOLD-000
      SUBROUTINE FOLD(V1,V2,VAL,NFO,IP,ISM,IST,IVY,INVZ,FR,PGN,XGN,WV,IZFOLD-001
     1Z,LO)                                                             FOLD-002
C V1 ARE THE POTENTIALS AND FORM FACTORS, UNFOLDED AS INPUT, FOLDED AS  FOLD-003
C OUTPUT. V2 ARE WORKING FIELDS IN WHICH 0. ARE STORED BEFORE THE CALL. FOLD-004
C INPUT:     V1:      POTENTIALS OR FORM FACTORS TO BE FOLDED.          FOLD-005
C            VAL(I,J):FOLDING PARAMETERS FOR THE FOLDING SETS J.        FOLD-006
C            NFO:     NUMBER OF SETS OF FOLDING PARAMETERS              FOLD-007
C            IP:      POTENTIAL TO BE FOLDED IF LO(7)=.FALSE.           FOLD-008
C            ISM:     NUMBER OF POINTS.                                 FOLD-009
C            IST:     MAXIMUM NUMBER OF STEPS FOR FOLDING FUNCTIONS.    FOLD-010
C            IVY:     TABLE OF FORM FACTORS (SEE REDM).                 FOLD-011
C            INVZ:    NUMBER OF TRANSITION FORM FACTORS TO FOLD.        FOLD-012
C            PGN:     WEIGHTS OF GAUSS-LEGENDRE INTEGRATION.            FOLD-013
C            XGN:   : ABSCISSAE OF GAUSS-LEGENDRE INTEGRATION.          FOLD-014
C            WV(J):   STEP SIZE FOR J=8.                                FOLD-015
C            LO(I):   LOGICAL CONTROLS:                                 FOLD-016
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    FOLD-017
C               LO(11) =.TRUE. DEFORMED COULOMB POTENTIAL.              FOLD-018
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            FOLD-019
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      FOLD-020
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. FOLD-021
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   FOLD-022
C               LO(100)=.TRUE. DIRAC EQUATION.                          FOLD-023
C               LO(101)=.TRUE. THERE IS A REAL SPIN-ORBIT POTENTIAL.    FOLD-024
C               LO(102)=.TRUE. THERE IS AN IMAGINARY SPIN-ORBIT         FOLD-025
C                              POTENTIAL.                               FOLD-026
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. FOLD-027
C OUTPUT:    V1:     :POTENTIALS OR FORM FACTORS FOLDED.                FOLD-028
C WORKING AREAS:                                                        FOLD-029
C            V2:      AREA INITIALISED TO 0.                            FOLD-030
C            FR(I,*): GAUSSIAN OR SAXON FOLDING FUNCTIONS FOR J=1,IMT.  FOLD-031
C                     FOR YUKAWA OR HULTHEN FOLDING: H-FUNCTIONS FOR    FOLD-032
C                     J=1,IMT, J-FUNCTIONS FOR J=IMT+1,2*IMT;           FOLD-033
C                     THESE FUNCTIONS ARE MULTIPLIED BY                 FOLD-034
C                     R*R*EXP(+\-R/VA). INTEGRALS WITH H AND J-FUNCTIONSFOLD-035
C                     FOR J=2*IMT+1,2*IMT+2.                            FOLD-036
C            IZZ(I,*):ADDRESS OF FOLDING PARAMETERS BELOW FOR I=1,      FOLD-037
C                     TYPE OF FORM FACTOR FROM 1 TO 16 FOR I=2,         FOLD-038
C                     ANGULAR MOMENTUM FOR I=3,                         FOLD-039
C                     ADDRESS OF STEP SIZE FOR I=4.                     FOLD-040
C                    (THIS PART IS AN INPUT WITH EXTERNAL FORM FACTORS) FOLD-041
C                                                                       FOLD-042
C FOR THE COMMON  /CONVE/ SEE CALX.                                     FOLD-043
C FOR THE COMMON  /POTE1/ SEE REDM.                                     FOLD-044
C                                                                       FOLD-045
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     FOLD-046
C  ACONV:     CONVERGENCE CRITERION FOR POTENTIAL AND FUNCTION.         FOLD-047
C   USED:     ACONV.                                                    FOLD-048
C                                                                       FOLD-049
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE1/:                     FOLD-050
C  ITX(16):   STARTING ADDRESS OF DIFFERENT FORM FACTORS (SEE REDM).    FOLD-051
C  INTC:      NUMBER OF FORM FACTORS WITHOUT DEFORMED SPIN-ORBIT        FOLD-052
C             INCLUDING CORRECTION TERMS.                               FOLD-053
C  INLS:      NUMBER OF SPIN-ORBIT FORM FACTORS NOT TAKING INTO ACCOUNT FOLD-054
C             MULTIPLICATION BY 2.                                      FOLD-055
C  INVD:      IDEM FOR COULOMB SPIN-ORBIT.                              FOLD-056
C  ITXM:      TOTAL NUMBER OF FORM FACTORS.                             FOLD-057
C   USED:     ITX,INTC,INLS,INVD,ITXM.                                  FOLD-058
C                                                                       FOLD-059
C FOLDING CONVENTIONS : VAL(*,1-3) REAL POTENTIAL, IMAGINARY, COULOMB.  FOLD-060
C IF VAL(1,*)=0  NO FOLDING.                                            FOLD-061
C IF VAL(3,*)=0 GAUSSIAN FORM FACTOR WITH RANGE VAL(2,*).               FOLD-062
C IF VAL(2,*)=0 HULTHEN FORM FACTOR WITH RANGES VAL(1,*) AND VAL(3,*).  FOLD-063
C IF VAL(2,*)=VAL(3,*)=0 YUKAWA FORM FACTOR WITH RANGE VAL(1,*).        FOLD-064
C ALL OTHER CASES  SAXON FORM FACTOR WITH RADIUS VAL(2,*) AND DIFFUSE-  FOLD-065
C NESS VAL(3,*) - ALL FORM FACTORS ARE NORMALISED - VAL(1,*) IS USED    FOLD-066
C ONLY FOR HULTHEN FORM FACTOR.                                         FOLD-067
C THE DIFFUSENESS OF A SAXON FORM FACTOR AND THE RANGES OF A HULTHEN OR FOLD-068
C YUKAWA FORM FACTOR ARE POSITIVE BY TAKING THE ABSOLUTE VALUE.         FOLD-069
C                                                                       FOLD-070
C***********************************************************************FOLD-071
      IMPLICIT REAL*8 (A-H,O-Z)                                         FOLD-072
      LOGICAL LO(150)                                                   FOLD-073
      DIMENSION V1(ISM,*),V2(ISM,*),VAL(3,3),IVY(7,*),FR(IST,*),PGN(10),FOLD-074
     1XGN(10),WV(22,*),IZZ(4,*),MM(3)                                   FOLD-075
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       FOLD-076
      COMMON /INOUT/ MR,MW,MS                                           FOLD-077
      COMMON /POTE1/ ITX(16),IMAX,INTC,INLS,INVC,INVD,ITXM              FOLD-078
      DATA PIM2 /6.2831853070D0/                                        FOLD-079
      AC=-DLOG(ACONV)                                                   FOLD-080
      JSL=ISM+1                                                         FOLD-081
      ITXN=ITXM                                                         FOLD-082
      IF (LO(100)) ITXN=ITXN-ITX(2)                                     FOLD-083
      IF (LO(7)) GO TO 9                                                FOLD-084
C INITIALISATION OF THE TABLE IZZ WHEN POTENTIALS ARE NOT EXTERNAL.     FOLD-085
      MM(1)=1                                                           FOLD-086
      MM(2)=2                                                           FOLD-087
      MM(3)=3                                                           FOLD-088
      IF ((VAL(1,2).EQ.VAL(1,1)).AND.(VAL(2,2).EQ.VAL(2,1)).AND.(VAL(3,2FOLD-089
     1).EQ.VAL(3,1))) MM(2)=1                                           FOLD-090
      IF ((VAL(1,3).EQ.VAL(1,2)).AND.(VAL(2,3).EQ.VAL(2,2)).AND.(VAL(3,3FOLD-091
     1).EQ.VAL(3,2))) MM(3)=2                                           FOLD-092
      IF ((VAL(1,3).EQ.VAL(1,1)).AND.(VAL(2,3).EQ.VAL(2,1)).AND.(VAL(3,3FOLD-093
     1).EQ.VAL(3,1))) MM(3)=1                                           FOLD-094
      NFOX=MAX0(MM(1),MM(2),MM(3))                                      FOLD-095
      DO 2 J=1,ITXN                                                     FOLD-096
      DO 1 I=1,3                                                        FOLD-097
    1 IZZ(I,J)=0                                                        FOLD-098
    2 IZZ(4,J)=IP                                                       FOLD-099
      DO 3 I=1,8                                                        FOLD-100
      IF ((I.EQ.5).AND.(.NOT.LO(101))) GO TO 3                          FOLD-101
      IF ((I.EQ.6).AND.(.NOT.LO(102))) GO TO 3                          FOLD-102
      IF ((I.EQ.8).AND.(.NOT.LO(103))) GO TO 3                          FOLD-103
      M=MOD(I+1,2)+1                                                    FOLD-104
      IF (I.GT.6) M=3                                                   FOLD-105
      N=I                                                               FOLD-106
      IF (.NOT.LO(100)) N=1+ITX(I)                                      FOLD-107
      IZZ(1,N)=-MM(M)                                                   FOLD-108
      IZZ(2,N)=I                                                        FOLD-109
    3 CONTINUE                                                          FOLD-110
      IF (INVZ.EQ.0) GO TO 9                                            FOLD-111
      IJ=ITX(3)-ITX(7)                                                  FOLD-112
      DO 8 J=1,INTC                                                     FOLD-113
      DO 7 I=1,8                                                        FOLD-114
      M=MOD(I+1,2)+1                                                    FOLD-115
      N=I+IJ                                                            FOLD-116
      IF (I.GT.6) GO TO 5                                               FOLD-117
      IF (I.GT.4) GO TO 4                                               FOLD-118
      IF ((M.EQ.2).AND.(.NOT.LO(12))) GO TO 7                           FOLD-119
      IF (.NOT.LO(100)) N=J+ITX(I+8)                                    FOLD-120
      GO TO 6                                                           FOLD-121
    4 IF (.NOT.LO(8+I)) GO TO 7                                         FOLD-122
      IF (.NOT.LO(100)) N=IVY(3,J)+ITX(I+8)                             FOLD-123
      GO TO 6                                                           FOLD-124
    5 IF (.NOT.LO(8*I-45)) GO TO 7                                      FOLD-125
      IF (.NOT.LO(100)) N=ITX(I+8)+IVY(I-3,J)                           FOLD-126
      M=3                                                               FOLD-127
    6 IZZ(1,N)=-MM(M)                                                   FOLD-128
      IZZ(2,N)=I+8                                                      FOLD-129
      IZZ(3,N)=IVY(7,J)                                                 FOLD-130
    7 CONTINUE                                                          FOLD-131
    8 IJ=IJ+11                                                          FOLD-132
C CLASSIFICATION OF STEP SIZES.                                         FOLD-133
    9 DO 11 I=2,ITXN                                                    FOLD-134
      IF (IZZ(1,I).EQ.0) GO TO 11                                       FOLD-135
      DO 10 J=2,I                                                       FOLD-136
      IF (IZZ(1,J-1).NE.IZZ(1,I)) GO TO 10                              FOLD-137
      IF (IZZ(4,J-1).EQ.IZZ(4,I)) GO TO 10                              FOLD-138
      I1=IZZ(4,I)                                                       FOLD-139
      J1=IZZ(4,J-1)                                                     FOLD-140
      IF (WV(8,I1).NE.WV(8,J1)) GO TO 10                                FOLD-141
      IZZ(4,I)=J1                                                       FOLD-142
      GO TO 11                                                          FOLD-143
   10 CONTINUE                                                          FOLD-144
   11 CONTINUE                                                          FOLD-145
      NFOT=-NFO                                                         FOLD-146
C FOR EXTERNAL POTENTIALS, ADD 1 FOR DUMMY FOLDING.                     FOLD-147
      IF (LO(7)) NFOX=NFO+1                                             FOLD-148
C LOOP ON FOLDINGS.                                                     FOLD-149
      DO 64 II=1,NFOX                                                   FOLD-150
      IF (II.GT.NFO) GO TO 54                                           FOLD-151
      IC=0                                                              FOLD-152
      DO 12 J=1,3                                                       FOLD-153
      IF (VAL(J,II).NE.0.D0) IC=IC+2**(3-J)                             FOLD-154
      IF (VAL(J,II).GE.0.D0) GO TO 12                                   FOLD-155
      WRITE (MW,1000) VAL(J,II),J,II                                    FOLD-156
      VAL(J,II)=-VAL(J,II)                                              FOLD-157
   12 CONTINUE                                                          FOLD-158
      IF (IC.LT.4) GO TO 64                                             FOLD-159
      IH=0                                                              FOLD-160
   13 IMT=1                                                             FOLD-161
      JH=0                                                              FOLD-162
      KH=0                                                              FOLD-163
C SEPARATION OF FOLDINGS WITH RESPECT TO THE STEP SIZE.                 FOLD-164
      DO 14 J=1,ITXN                                                    FOLD-165
      IF (IZZ(1,J)+II.NE.0) GO TO 14                                    FOLD-166
      IF (IZZ(4,J).LT.IH) GO TO 14                                      FOLD-167
      KH=KH+1                                                           FOLD-168
      IF ((JH.NE.0).AND.(IZZ(4,J).GT.IH)) GO TO 14                      FOLD-169
      JH=JH+1                                                           FOLD-170
      IH=IZZ(4,J)                                                       FOLD-171
      IZZ(1,J)=NFOT-1                                                   FOLD-172
      IMT=MAX0(IMT,IZZ(3,J)+1)                                          FOLD-173
   14 CONTINUE                                                          FOLD-174
      IF (JH.EQ.0) GO TO 64                                             FOLD-175
      NFOT=NFOT-1                                                       FOLD-176
      HH=WV(8,IH)                                                       FOLD-177
C COMPUTATION OF THE STRENGTH OF THE FOLDING FUNCTIONS.                 FOLD-178
      CN=1.D0                                                           FOLD-179
      IF (IC.GT.5) GO TO 15                                             FOLD-180
C YUKAWA FOLDING.                                                       FOLD-181
      IF (VAL(3,II).EQ.VAL(1,II)) VAL(3,II)=.99D0*VAL(1,II)             FOLD-182
      CN=HH/(2.D0*PIM2*(VAL(1,II)**2-VAL(3,II)**2))                     FOLD-183
      JST=ISM+5                                                         FOLD-184
      GO TO 18                                                          FOLD-185
   15 IF (IC.EQ.7) GO TO 16                                             FOLD-186
C GAUSSIAN FOLDING.                                                     FOLD-187
      ISY=1+IDINT(VAL(2,II)*DSQRT(AC)/HH)                               FOLD-188
      CN=HH/(0.5D0*PIM2*DSQRT(0.5D0*PIM2)*VAL(2,II)**3)                 FOLD-189
      GO TO 18                                                          FOLD-190
C SAXON FOLDING.                                                        FOLD-191
   16 IS=1+IDINT((VAL(2,II)+10.D0*VAL(3,II))/HH)                        FOLD-192
      ISY=1+IDINT((VAL(2,II)+2.D0*VAL(3,II)*AC)/HH)                     FOLD-193
      CN=0.D0                                                           FOLD-194
      DO 17 K=1,IS                                                      FOLD-195
      RR=HH*K                                                           FOLD-196
   17 CN=CN+RR*RR/(1.D0+DEXP((RR-VAL(2,II))/VAL(3,II)))                 FOLD-197
      CN=.5D0/(PIM2*CN)                                                 FOLD-198
C GAUSSIAN OR SAXON FOLDING.                                            FOLD-199
   18 AN=DFLOAT(IMT-1)                                                  FOLD-200
      IF (IC.LT.6) GO TO 36                                             FOLD-201
C DO LOOP ON THE POINTS OF POTENTIALS.                                  FOLD-202
      FS=0.D0                                                           FOLD-203
      IS=0                                                              FOLD-204
   19 IS=IS+1                                                           FOLD-205
      JST=IS+ISY                                                        FOLD-206
      JSM=MIN0(ISM,JST)                                                 FOLD-207
      R=VAL(2,II)                                                       FOLD-208
      A=VAL(3,II)                                                       FOLD-209
      X1=HH*DFLOAT(IS)                                                  FOLD-210
      IF (JST.LE.ISM) GO TO 20                                          FOLD-211
      JSZ=MAX0(ISM,1+IS+2*IDINT((R+2.D0*A)/HH))                         FOLD-212
      JST=MIN0(JST,JSZ+4)                                               FOLD-213
      ISZ=JSZ-1                                                         FOLD-214
      ISX=JST-JSZ                                                       FOLD-215
C LOOP ON THE POINTS OF THE FOLDING FUNCTION FOR WHICH THE SYMMETRY     FOLD-216
C BETWEEN THE TWO ARGUMENTS IS USED.                                    FOLD-217
   20 JI=IS-1                                                           FOLD-218
   21 JI=JI+1                                                           FOLD-219
      A3=HH*JI                                                          FOLD-220
      A4=2.D0*X1*A3                                                     FOLD-221
      Y=X1*X1+A3*A3                                                     FOLD-222
      IF (A.EQ.0.D0) GO TO 25                                           FOLD-223
C SAXON FOLDING FUNCTION  BY A 20-POINTS GAUSSIAN INTEGRATION           FOLD-224
C THE SAXON FORM FACTOR ITSELF IS THE VARIABLE.                         FOLD-225
      DO 22 K=1,IMT                                                     FOLD-226
   22 FR(JI,K)=0.D0                                                     FOLD-227
      R1=DMIN1(70.D0,(DABS(X1-A3)-R)/A)                                 FOLD-228
      R2=DMIN1(70.D0,(X1+A3-R)/A)                                       FOLD-229
      E1=1.D0/(1.D0+DEXP(R1))                                           FOLD-230
      E2=1.D0/(1.D0+DEXP(R2))                                           FOLD-231
      EM=0.5D0*(E1+E2)                                                  FOLD-232
      ES=0.5D0*(E1-E2)                                                  FOLD-233
      P1=ES*2.D0*A*PIM2/A4                                              FOLD-234
      DO 24 IJ=1,10                                                     FOLD-235
      PP=P1*PGN(IJ)                                                     FOLD-236
      E1=1.D0/(EM+ES*XGN(IJ))-1.D0                                      FOLD-237
      E2=1.D0/(EM-ES*XGN(IJ))-1.D0                                      FOLD-238
      U1=R+A*DLOG(E1)                                                   FOLD-239
      U2=R+A*DLOG(E2)                                                   FOLD-240
      R1=PP*U1*(1.D0+E1)/E1                                             FOLD-241
      R2=PP*U2*(1.D0+E2)/E2                                             FOLD-242
      FR(JI,1)=FR(JI,1)+R1+R2                                           FOLD-243
      IF (IMT.EQ.1) GO TO 24                                            FOLD-244
      E1=(Y-U1*U1)/A4                                                   FOLD-245
      E2=(Y-U2*U2)/A4                                                   FOLD-246
      R3=R1*E1                                                          FOLD-247
      R4=R2*E2                                                          FOLD-248
      FR(JI,2)=FR(JI,2)+R3+R4                                           FOLD-249
      IF (IMT.EQ.2) GO TO 24                                            FOLD-250
C RECURRENCE FOR LEGENDRE POLYNOMIALS.                                  FOLD-251
      DO 23 K=3,IMT                                                     FOLD-252
      R5=R1                                                             FOLD-253
      R6=R2                                                             FOLD-254
      R1=R3                                                             FOLD-255
      R2=R4                                                             FOLD-256
      R3=(DFLOAT(2*K-3)*R1*E1-R5*DFLOAT(K-2))/DFLOAT(K-1)               FOLD-257
      R4=(DFLOAT(2*K-3)*R2*E2-R6*DFLOAT(K-2))/DFLOAT(K-1)               FOLD-258
   23 FR(JI,K)=FR(JI,K)+R3+R4                                           FOLD-259
   24 CONTINUE                                                          FOLD-260
      GO TO 29                                                          FOLD-261
C GAUSSIAN FOLDING FUNCTION.                                            FOLD-262
   25 A4=A4/(R*R)                                                       FOLD-263
      Y=Y/(R*R)                                                         FOLD-264
      B1=DEXP(A4-Y)                                                     FOLD-265
      B2=0.D0                                                           FOLD-266
      IF (A4+Y.LT.50.) B2=DEXP(-A4-Y)                                   FOLD-267
      DN=1.D0/A4                                                        FOLD-268
   26 FR(JI,1)=PIM2*(B1-B2)*DN                                          FOLD-269
      IF (IMT.EQ.1) GO TO 29                                            FOLD-270
C DOWNWARDS RECURRENCE FOR SMALL ARGUMENTS.                             FOLD-271
      Q=DMAX1(DSQRT(10.5D0*A4)-0.5D0,AN)                                FOLD-272
      K=IDINT(Q+3.D0+21.D0*A4/(Q+Q+1.D0))                               FOLD-273
      A1=0.D0                                                           FOLD-274
   27 A1=A4/(2.D0*K+1.D0+A4*A1)                                         FOLD-275
      IF (K.LT.IMT) FR(JI,K+1)=A1                                       FOLD-276
      K=K-1                                                             FOLD-277
      IF (K.GT.0) GO TO 27                                              FOLD-278
      DO 28 K=2,IMT                                                     FOLD-279
   28 FR(JI,K)=FR(JI,K)*FR(JI,K-1)                                      FOLD-280
   29 IF (IC.LT.6) GO TO 39                                             FOLD-281
      IF (JI.LT.JST) GO TO 21                                           FOLD-282
      IK=1                                                              FOLD-283
      GO TO 46                                                          FOLD-284
   30 IK=2                                                              FOLD-285
      IF (JST.LE.ISM) GO TO 46                                          FOLD-286
C FORM THE COULOMB POTENTIALS,ASYMPTOTIC CORRECTION.                    FOLD-287
C THE FORM FACTORS ARE ASSUMED TO DECREASE AS (R)**(-L-1).              FOLD-288
      A1=DFLOAT(ISM)                                                    FOLD-289
      DO 32 I=JSL,JST                                                   FOLD-290
      Y=A1/DFLOAT(I)                                                    FOLD-291
      A4=DFLOAT(ISM*I)*HH**2                                            FOLD-292
      DO 31 J=1,IMT                                                     FOLD-293
      FR(I,J)=FR(I,J)*A4                                                FOLD-294
   31 A4=A4*Y                                                           FOLD-295
   32 CONTINUE                                                          FOLD-296
C CORRECTION OF THE LAST VALUE BY SUM BETWEEN ISM AND JSZ AND A PADE    FOLD-297
C APPROXIMATION LIMITED TO FOUR TERMS (BETWEEN JSZ+1 AND JST).          FOLD-298
      DO 35 I=1,IMT                                                     FOLD-299
      A1=0.D0                                                           FOLD-300
      A2=0.D0                                                           FOLD-301
      A3=0.D0                                                           FOLD-302
      A4=0.D0                                                           FOLD-303
      Y=0.D0                                                            FOLD-304
      IF (JSL.GT.ISZ) GO TO 34                                          FOLD-305
      DO 33 J=JSL,ISZ                                                   FOLD-306
   33 Y=Y+FR(J,I)                                                       FOLD-307
   34 IF (ISX.EQ.0.OR.FR(JSZ,I).EQ.0.D0) GO TO 35                       FOLD-308
      A1=-FR(JSZ+1,I)/FR(JSZ,I)                                         FOLD-309
      IF (ISX.EQ.1.OR.FR(JSZ+1,I).EQ.0.D0) GO TO 35                     FOLD-310
      B1=FR(JSZ+2,I)/FR(JSZ+1,I)                                        FOLD-311
      A2=-A1-B1                                                         FOLD-312
      IF (ISX.EQ.2.OR.A2.EQ.0.D0) GO TO 35                              FOLD-313
      B2=FR(JSZ+3,I)/FR(JSZ+1,I)                                        FOLD-314
      C1=-(B1*A1+B2)/A2                                                 FOLD-315
      A3=B1-C1                                                          FOLD-316
      IF (ISX.EQ.3.OR.A3.EQ.0.D0) GO TO 35                              FOLD-317
      A4=C1-(B2+(B2*A1+FR(JSZ+4,I)/FR(JSZ+1,I))/A2)/A3                  FOLD-318
   35 FR(JSL,I)=Y+FR(JSZ,I)/(1.D0+A1/(1.D0+A2/(1.D0+A3/(1.D0+A4))))     FOLD-319
      GO TO 46                                                          FOLD-320
C LOOP ON THE TWO YUKAWA FORM FACTORS.                                  FOLD-321
   36 FS=1.D0                                                           FOLD-322
      VA=VAL(1,II)                                                      FOLD-323
   37 IF (VA.EQ.0.D0) GO TO 53                                          FOLD-324
C COMPUTATION OF BESSEL FUNCTIONS MULTIPLIED BY R*R AND EXP(R/VA).      FOLD-325
      JI=0                                                              FOLD-326
      B1=1.D0                                                           FOLD-327
      B2=1.D0                                                           FOLD-328
      B4=DEXP(-HH/VA)                                                   FOLD-329
   38 JI=JI+1                                                           FOLD-330
      A2=HH*DFLOAT(JI)                                                  FOLD-331
      A4=A2/VA                                                          FOLD-332
      B2=B2*B4**2                                                       FOLD-333
      FR(JI,IMT+1)=FS*A2                                                FOLD-334
      DN=A2*VA                                                          FOLD-335
      GO TO 26                                                          FOLD-336
C UPWARDS RECURRENCE FOR THE IRREGULAR FUNCTION.                        FOLD-337
   39 FR(JI,IMT+2)=FR(JI,IMT+1)*(1.D0+1.D0/A4)                          FOLD-338
      IF (IMT.EQ.2) GO TO 41                                            FOLD-339
      DO 40 K=3,IMT                                                     FOLD-340
   40 FR(JI,K+IMT)=FR(JI,K-2+IMT)+FR(JI,K-1+IMT)*(2*K-3)/A4             FOLD-341
   41 IF (JI.LT.JST) GO TO 38                                           FOLD-342
C CORRECTIONS FOR SINGULAR FIRST DERIVATIVE (YUKAWA FORM FACTOR).       FOLD-343
      VR=-FS*HH*PIM2/6.D0                                               FOLD-344
      IK=1                                                              FOLD-345
      GO TO 46                                                          FOLD-346
C CORRECTION FOR COULOMB POTENTIALS.                                    FOLD-347
   42 A6=DFLOAT(ISM)                                                    FOLD-348
      DO 44 I=JSL,JST                                                   FOLD-349
      Y=A6/DFLOAT(I)                                                    FOLD-350
      A1=Y                                                              FOLD-351
      DO 43 J=1,IMT                                                     FOLD-352
      FR(I,J+IMT)=FR(I,J+IMT)*A1                                        FOLD-353
   43 A1=A1*Y                                                           FOLD-354
   44 CONTINUE                                                          FOLD-355
C ESTIMATION OF THE LAST VALUE BY A PADE OF FOUR TERMS.                 FOLD-356
      DO 45 I=1,IMT                                                     FOLD-357
      A6=-FR(JSL+1,I+IMT)/FR(JSL,I+IMT)                                 FOLD-358
      IF (I.EQ.1) GO TO 45                                              FOLD-359
      B2=FR(JSL+2,I+IMT)/FR(JSL+1,I+IMT)                                FOLD-360
      A2=-A6-B2                                                         FOLD-361
      B3=FR(JSL+3,I+IMT)/FR(JSL+1,I+IMT)                                FOLD-362
      C1=-(B2*A6+B3)/A2                                                 FOLD-363
      A4=B2-C1                                                          FOLD-364
      A1=C1-(B3+(B3*A6+FR(JSL+4,I+IMT)/FR(JSL+1,I+IMT))/A2)/A4          FOLD-365
      A6=A6/(1.D0+B4*A2/(1.D0+B4*A4/(1.D0+B4*A1)))                      FOLD-366
   45 FR(JSL,I+IMT)=FR(JSL,I+IMT)/(1.D0+B4*A6)                          FOLD-367
      IK=2                                                              FOLD-368
C FOLDING OF THE POTENTIALS.                                            FOLD-369
   46 DO 51 K=1,ITXN                                                    FOLD-370
      IF (IZZ(1,K).NE.NFOT) GO TO 51                                    FOLD-371
      IF ((IK.EQ.1).AND.MOD(IZZ(2,K),8).GT.6) GO TO 51                  FOLD-372
      IF ((IK.EQ.2).AND.MOD(IZZ(2,K),8).LE.6) GO TO 51                  FOLD-373
      N=IZZ(3,K)+1                                                      FOLD-374
      IF (IC.LT.6) GO TO 48                                             FOLD-375
      DO 47 JS=IS,JSM                                                   FOLD-376
      B1=(HH*DFLOAT(JS))**2                                             FOLD-377
      IF (JS.EQ.IS) GO TO 47                                            FOLD-378
      V2(JS,K)=V2(JS,K)+V1(IS,K)*FR(JS,N)*X1**2                         FOLD-379
   47 V2(IS,K)=V2(IS,K)+V1(JS,K)*FR(JS,N)*B1                            FOLD-380
      IF ((IK.EQ.2).AND.(JST.GT.ISM)) V2(IS,K)=V2(IS,K)+V1(ISM,K)*FR(JSLFOLD-381
     1,N)                                                               FOLD-382
      GO TO 51                                                          FOLD-383
   48 B2=0.D0                                                           FOLD-384
      B3=FR(JSL,N+IMT)*V1(ISM,K)                                        FOLD-385
      DO 49 IS=1,ISM                                                    FOLD-386
      JS=JSL-IS                                                         FOLD-387
      FR(JS,2*IMT+1)=B3*B4                                              FOLD-388
      B2=B2*B4+FR(IS,N)*V1(IS,K)                                        FOLD-389
      B3=B3*B4+FR(JS,N+IMT)*V1(JS,K)                                    FOLD-390
   49 FR(IS,2*IMT+2)=B2                                                 FOLD-391
      DO 50 IS=1,ISM                                                    FOLD-392
      RR=(IS*HH)**2                                                     FOLD-393
   50 V2(IS,K)=V2(IS,K)+(FR(IS,N)*FR(IS,2*IMT+1)+FR(IS,N+IMT)*FR(IS,2*IMFOLD-394
     1T+2))/RR+V1(IS,K)*VR                                              FOLD-395
   51 CONTINUE                                                          FOLD-396
      IF (IC.LT.6) GO TO 52                                             FOLD-397
      IF (IK.EQ.1) GO TO 30                                             FOLD-398
      IF (IS.LT.ISM) GO TO 19                                           FOLD-399
   52 IF (IK.EQ.1) GO TO 42                                             FOLD-400
   53 FS=-FS                                                            FOLD-401
      IF (FS.GE.0.D0) GO TO 55                                          FOLD-402
      VA=VAL(3,II)                                                      FOLD-403
      GO TO 37                                                          FOLD-404
C TRANSFER OF POTENTIALS FROM WORKING SPACE.                            FOLD-405
   54 NFOT=1                                                            FOLD-406
   55 DO 63 J=1,ITXN                                                    FOLD-407
      IF (IZZ(1,J).NE.NFOT) GO TO 63                                    FOLD-408
      IF (NFOT.EQ.1) GO TO 57                                           FOLD-409
      DO 56 IS=1,ISM                                                    FOLD-410
   56 V1(IS,J)=CN*V2(IS,J)                                              FOLD-411
      GO TO 58                                                          FOLD-412
C STEP SIZE FOR DUMMY FOLDING.                                          FOLD-413
   57 IH=IZZ(4,J)                                                       FOLD-414
      HH=WV(8,IH)                                                       FOLD-415
   58 IF (.NOT.LO(100).OR.(IZZ(2,J).GT.8)) GO TO 59                     FOLD-416
      CALL DERI(V1(1,J+8),V1(1,J),HH,ISM,.TRUE.)                        FOLD-417
      CALL DERI(V1(1,J+16),V1(1,J+8),HH,ISM,.TRUE.)                     FOLD-418
      GO TO 63                                                          FOLD-419
   59 IX=MOD(IZZ(2,J)-1,8)+1                                            FOLD-420
      IF ((IX.LT.5).OR.(IX.EQ.7)) GO TO 63                              FOLD-421
      IML=0                                                             FOLD-422
      INL=0                                                             FOLD-423
      IF (LO(100)) GO TO 60                                             FOLD-424
      IF (IZZ(2,J).LE.8) GO TO 61                                       FOLD-425
      INL=INLS                                                          FOLD-426
      IF (IX.EQ.8) INL=INVD                                             FOLD-427
      GO TO 61                                                          FOLD-428
   60 IF (IZZ(2,J).LE.8) GO TO 63                                       FOLD-429
      IML=4                                                             FOLD-430
      IF (IX.EQ.8) IML=3                                                FOLD-431
   61 CALL DERI(V2(1,J),V1(1,J),HH,ISM,LO(100))                         FOLD-432
      DO 62 IS=1,ISM                                                    FOLD-433
      IF (INL.EQ.0) GO TO 62                                            FOLD-434
      AY=DFLOAT(IS)*HH                                                  FOLD-435
      V1(IS,J+INL)=-V1(IS,J)/AY**2                                      FOLD-436
   62 V1(IS,J+IML)=V2(IS,J)                                             FOLD-437
   63 CONTINUE                                                          FOLD-438
      IF (KH.NE.JH.AND.NFOT.NE.1) GO TO 13                              FOLD-439
   64 CONTINUE                                                          FOLD-440
      RETURN                                                            FOLD-441
 1000 FORMAT (' NEGATIVE VALUE',D15.5,4X,'OF FOLDING PARAMETER VAL(',I3,FOLD-442
     1',',I3,')')                                                       FOLD-443
      END                                                               FOLD-444
