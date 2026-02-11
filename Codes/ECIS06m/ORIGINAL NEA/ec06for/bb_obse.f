C 21/02/07                                                      ECIS06  OBSE-000
      SUBROUTINE OBSE(MF,CMF,NT,NCOLR,IPI,JTN,LO,AM,JCAL,ITX,ITZ,BM,CBM,OBSE-001
     1NX,DM,NY)                                                         OBSE-002
C COMPUTES FOR EACH OBSERVABLE ALL THE INDICATIONS FOR DO LOOPS AND THE OBSE-003
C GEOMETRICAL COEFFICIENTS NEEDED IN SUBROUTINE SCAT.                   OBSE-004
C INPUT:     MF,CMF:  EQUIVALENT BY CALL, INFORMATIONS COMING FROM      OBSE-005
C                     SUBROUTINES DEPH AND LECD.                        OBSE-006
C            NT:      NUMBER OF ROWS OF MF.                             OBSE-007
C            NCOLR :  NUMBER OF EXPERIMENTAL ANGULAR DISTRIBUTIONS.     OBSE-008
C            IPI(J,*):MULTIPLICITY OF PARTICLE AND TARGET FOR J+2 AND 3.OBSE-009
C            JTN:     LENGTH AVAILABLE IN AM.                           OBSE-010
C            LO(I):   LOGICAL CONTROLS:                                 OBSE-011
C               LO(31) =.TRUE. INPUT OF EXPERIMENTAL DATA AND CHI2      OBSE-012
C                              CALCULATION.                             OBSE-013
C               LO(92) =.TRUE. NON STANDARD OBSERVABLES AT EQUIDISTANT  OBSE-014
C                              ANGLES.                                  OBSE-015
C               LO(126)=.TRUE. SOME OBSERVABLES IN LABORATORY SYSTEM.   OBSE-016
C OUTPUT:    MF,CMF:  SECOND PART DESCRIBED IN SUBROUTINE DEPH UPDATED. OBSE-017
C            AM:      CHARACTER*4 FOR TRANSFER, SEE DESCRIPTION BELOW.  OBSE-018
C            JCAL:    LENGTH OF AM.                                     OBSE-019
C            LO(I):   LOGICAL CONTROLS: LO(126) IS DEFINED HERE.        OBSE-020
C WORKING AREAS:                                                        OBSE-021
C            ITX:     EQUIVALENT WITH AM BY CALL, TO STORE QUANTUM      OBSE-022
C                     NUMBERS OF OBSERVABLES.                           OBSE-023
C            ITZ:     EQUIVALENT WITH AM BY CALL, TO STORE LEGENDS OF   OBSE-024
C                     OBSERVABLES.                                      OBSE-025
C            BM(9,*): IN EQUIVALENCE WITH AM SHIFTED OF 4 TIMES THE     OBSE-026
C                     NUMBER OF STANDARD AND NON STANDARD OBSERVABLES.  OBSE-027
C                     USED FOR COMPUTATION OF CLEBSCH-GORDAN            OBSE-028
C                     COEFFICIENTS BY RECURRENCE.                       OBSE-029
C            CBM:     IN EQUIVALENCE BY CALL WITH BM. USED FOR TRANSFER OBSE-030
C                     OF DATA TO AM.                                    OBSE-031
C            NX(18,*):IN EQUIVALENCE BY CALL WITH BM. USED FOR          OBSE-032
C                     MANIPULATION OF QUANTUM NUMBERS.                  OBSE-033
C            NY(*):   IN EQUIVALENCE BY CALL WITH NX.                   OBSE-034
C                                                                       OBSE-035
C***********************************************************************OBSE-036
C  ALL THE OBSERVABLES (AS GIVEN FOR SOME OF THEM IN SUBROUTINE DEPH FOROBSE-037
C     THE DESCRIPTION OF THE SECOND PART OF TABLE MF) ARE DEFINED AS A  OBSE-038
C     SUM OF OPERATORS "A(L1,M1,L2,M2 L1,M1,L2,M2)" WITH SOME NUMERICAL OBSE-039
C     COEFFICIENT. WITH SCATTERING AMPLITUDE B(I1,I2,I3,I4) FOR A GIVEN OBSE-040
C     ANGLE, THE CONTRIBUTION OF "A" IS GIVEN BY A QUADRUPLE SUM ON     OBSE-041
C     I1, I2, I3 AND I4 FROM IM1 TO IP1, IM2 TO IP2, IM3 TO IP3, IM4 TO OBSE-042
C     IP4 RESPECTIVELY OF:                                              OBSE-043
C     B(I1,I2,I3,I4)(COMPLEX CONJUGATE)*B(I1-NM1,I2-NM2,I3-NM3,I4-NM4)* OBSE-044
C     AM(IL1+I1)*AM(IL2-I2)*AM(IL3+I3)*AM(IL4-I4)                       OBSE-045
C     WITH AM(IL1+I1), AM(IL2-I2),... CLEBSCH-GORDAN COEFFICIENTS STOREDOBSE-046
C     IN THE SECOND PART OF TABLE AM ASSUMED TO BE DOUBLE PRECISION.    OBSE-047
C        DUE TO FREQUENT REDUCTION OF THESE CLEBSCH-GORDAN COEFFICIENTS OBSE-048
C     TO UNITY, THERE IS A LOGICAL FOR EACH DO LOOP WHICH ARE .TRUE. IF OBSE-049
C     THERE IS A TENSOR OPERATOR FOR THE RELATED DO LOOP,.FALSE. FOR    OBSE-050
C     IDENTITY OPERATOR; A FIFTH LOGICAL IS .TRUE. FOR NO TENSOR        OBSE-051
C     OPERATOR AT ALL AND A SIXTH IS .TRUE. FOR PURE IMAGINARY          OBSE-052
C     COEFFICIENT. THESE SIX LOGICALS ARE STORED ON A BINARY BASIS      OBSE-053
C     IN "JQ".                                                          OBSE-054
C        INDICATION FOR CHANGE OF FRAME, N1 FOR THE PARTICLE AND N2 THE OBSE-055
C     TARGET:                                                           OBSE-056
C            1 FOR LABORATORY SYSTEM,                                   OBSE-057
C            2 FOR AXIS ALONG THE INCIDENT DIRECTION,                   OBSE-058
C     ARE STORED UNDER THE FORM "JT"=N1+1000*N2.                        OBSE-059
C                                                                       OBSE-060
C  **** THE FIRST PART OF THE AM TABLE IS SETS OF 20 INTEGERS, EACH OF  OBSE-061
C     THEM RELATED TO A SINGLE "A" - THE 16 FIRST ARE INDICATIONS FOR   OBSE-062
C     THE 4 DO LOOPS: IM1, IP1, NM1, IL1, IM2, IP2, NM2, IL2, IM3, IP3, OBSE-063
C     NM3, IL3, IM4, IP4, NM4 AND IL4. THERE ARE FOLLOWED BY "JQ" AND   OBSE-064
C     "JT" AND THE COEFFICIENT OF THIS "A" IN DOUBLE PRECISION.         OBSE-065
C  **** THE SECOND PART ARE TENSOR OPERATOR MATRIX ELEMENTS WHICH ARE   OBSE-066
C     QUITE PURELY CLEBSCH-GORDAN COEFFICIENTS.                         OBSE-067
C                                                                       OBSE-068
C***********************************************************************OBSE-069
      IMPLICIT REAL*8 (A-H,O-Z)                                         OBSE-070
      LOGICAL LT1,LT2,LZ(6),LO(150)                                     OBSE-071
      DIMENSION MF(10,*),IPI(11,*),ITX(7,*),BM(9,*),IY(2,20),MT(4),NX(18OBSE-072
     1,*),MX(8,24),NZ(2),CX(24),DM(*),NY(*)                             OBSE-073
      CHARACTER*4 CMF(10,*),AM(*),ITZ(7,* ),CBM(*),IZ(5,20)             OBSE-074
      COMMON /INOUT/ MR,MW,MS                                           OBSE-075
      DATA IY /4*1,2*2,2*3,2*4,2*5,2*6,7,8,9,10,2*11,2*12,2*13,14,16,2*1OBSE-076
     17,18,19,20,21,2*22,2*23,2*24,2*1/                                 OBSE-077
      DATA CX /6*1.D0,2*-1.D0,1.D0,-1.D0,1.D0,2*-1.4142136573791504D0,3*OBSE-078
     1.5D0,1.4142136999999999D0,2*-1.D0,1.D0,-1.D0,1.D0,2*-1.41421365737OBSE-079
     291504D0/                                                          OBSE-080
      DATA MX,JTS /8*0,2*1,10*0,2*1,2*0,2,7*0,2,1,6*0,2*2,6*0,2*1,2*0,2*OBSE-081
     11,2*0,2*1,2*0,1,-1,2*0,2*1,2*0,2*1,2*0,2*1,2*0,1,-1,2*0,1,3*0,1,3*OBSE-082
     20,2*1,2*0,1,3*0,1,3*0,2*1,10*0,2*1,2*0,2*1,2*0,2*1,2*0,1,-1,4*0,2*OBSE-083
     31,4*0,4*1,4*0,3*1,-1,4*0,4*1,4*0,3*1,-1,4*0,1,0,1,5*0,3*1,5*0,1,0,OBSE-084
     42*1,5*0/                                                          OBSE-085
      DATA IZ /'   C','ROSS','-SEC','TION','    ','   C','. S.','/RUT','OBSE-086
     1HER.','    ','   A','SYM.',' OR ','IT11',2*'    ','VECT','. PO','LOBSE-087
     2AR.',4*'    ',' T20',4*'    ',' T21',4*'    ',' T22',3*'    ','KYYOBSE-088
     3 ','OR D',3*'    ','KXX ','OR R',2*'    ','   K','ZZ O','R A''',3*OBSE-089
     4'    ','KZX ','OR A',2*'    ','   K','XZ O','R R''',2*'    ','   SOBSE-090
     5','PIN-','FLIP',2*'    ','TARG','ET A','SYM.',4*'    ',' AYY',4*' OBSE-091
     6   ',' AXX',4*'    ' ,' AZZ',4*'    ',' AXZ',4*'    ',' AZX',2*'  OBSE-092
     7  ',' TOT','AL C','.-S.','    '/                                  OBSE-093
      IF (JTN.LT.270) CALL MEMO('OBSE',JTN,270)                         OBSE-094
C TRANSFER OF STANDARD DESCRIPTIONS.                                    OBSE-095
      DO 2 I=1,24                                                       OBSE-096
      DO 1 J=1,8                                                        OBSE-097
    1 NX(J+6,I)=MX(J,I)                                                 OBSE-098
    2 BM(3,I)=CX(I)                                                     OBSE-099
      DO 4 I=1,20                                                       OBSE-100
      ITX(1,I)=IY(1,I)                                                  OBSE-101
      ITX(2,I)=IY(2,I)                                                  OBSE-102
      DO 3 J=3,7                                                        OBSE-103
    3 ITZ(J,I)=IZ(J-2,I)                                                OBSE-104
    4 CONTINUE                                                          OBSE-105
      LO(126)=.FALSE.                                                   OBSE-106
      KIT=20                                                            OBSE-107
      KNX=24                                                            OBSE-108
      JCAL=0                                                            OBSE-109
      MT(1)=IPI(2,1)                                                    OBSE-110
      MT(2)=IPI(3,1)                                                    OBSE-111
C VERIFICATION OF INDICATIONS ALREADY IN MF(1,I).                       OBSE-112
      DO 5 I=1,NT                                                       OBSE-113
      IF (MF(2,I).LE.18) GO TO 5                                        OBSE-114
      IF (MF(2,I).EQ.19.AND.NT-I.GE.NCOLR) GO TO 83                     OBSE-115
      IF (MF(2,I).GT.19) GO TO 84                                       OBSE-116
    5 MF(6,I)=0                                                         OBSE-117
      IF (.NOT.(LO(92).OR.LO(31))) GO TO 56                             OBSE-118
C SEARCH FOR NON STANDARD OBSERVABLES.                                  OBSE-119
    6 DO 7 L=1,NT                                                       OBSE-120
      IF (MF(2,L).LT.0) GO TO 8                                         OBSE-121
    7 CONTINUE                                                          OBSE-122
      GO TO 56                                                          OBSE-123
    8 KZ=KIT                                                            OBSE-124
      KIT=KIT+1                                                         OBSE-125
C INPUT OF THE DESCRIPTION OF NON STANDARD OBSERVABLES.                 OBSE-126
      LT1=.FALSE.                                                       OBSE-127
      LT2=.FALSE.                                                       OBSE-128
      READ (MR,1000) LT1,LT2,LT3,KX,K,(ITZ(J,KIT),J=3,7)                OBSE-129
      K1=KNX+1                                                          OBSE-130
      KNX=KNX+K                                                         OBSE-131
      IF (9*KNX.GT.JTN) CALL MEMO('OBSE',JTN,9*KNX)                     OBSE-132
      READ (MR,1001) ((NX(J,I),J=7,14),I=K1,KNX)                        OBSE-133
      READ (MR,1002) (BM(3,I),I=K1,KNX)                                 OBSE-134
      IF (LT1.OR.LT2) READ (MR,1002) (BM(1,I),I=K1,KNX)                 OBSE-135
    9 K2=K1                                                             OBSE-136
      WRITE (MW,1003) KX,(ITZ(J,KIT),J=3,7),K                           OBSE-137
      IF (LT1) WRITE (MW,1004)                                          OBSE-138
      IF (LT3.EQ.1) WRITE (MW,1005)                                     OBSE-139
      IF (LT3.EQ.2) WRITE (MW,1006)                                     OBSE-140
      IF (.NOT.LT2) GO TO 22                                            OBSE-141
C OUTPUT OF THE DESCRIPTION NOT COMPLETELY IN TENSOR NOTATION           OBSE-142
C IF LT2=.TRUE. , OBSERVABLES CAN BE DEFINED AS MATRIX ELEMENTS         OBSE-143
C (MI,MF). IN THIS CASE,THE VALUES READ ARE (MI-S-1,MF-S-1) INSTEAD OF  OBSE-144
C THE QUANTUM NUMBERS OF TENSOR OPERATOR.                               OBSE-145
      WRITE (MW,1007)                                                   OBSE-146
   10 K3=MIN0(K2+5,KNX)                                                 OBSE-147
      WRITE (MW,1008) ((NX(I,J),I=11,14),J=K2,K3)                       OBSE-148
      WRITE (MW,1009) (BM(3,I),I=K2,K3)                                 OBSE-149
      WRITE (MW,1010) ((NX(I,J),I=7,10),J=K2,K3)                        OBSE-150
      WRITE (MW,1011) (BM(1,I),I=K2,K3)                                 OBSE-151
      K2=K3+1                                                           OBSE-152
      IF (K2.LE.KNX) GO TO 10                                           OBSE-153
      IV=MF(1,L)                                                        OBSE-154
      MT(3)=IPI(2,IV)                                                   OBSE-155
      MT(4)=IPI(3,IV)                                                   OBSE-156
      KLT=1                                                             OBSE-157
C SEARCH FOR NON TENSOR DESCRIPTION - PSEUDO-LOOP FOR IJ.               OBSE-158
      IJ=1                                                              OBSE-159
   11 K2=K1                                                             OBSE-160
      K6=KNX                                                            OBSE-161
   12 IF (K2.GT.K6) GO TO 14                                            OBSE-162
      DO 13 K=K2,K6                                                     OBSE-163
      IF (NX(2*IJ+5,K).LT.0) GO TO 15                                   OBSE-164
   13 CONTINUE                                                          OBSE-165
   14 IF (K6.NE.KNX) GO TO 27                                           OBSE-166
      GO TO 20                                                          OBSE-167
C RECURRENCE COMPUTATION OF (-)**(S-MF)*<S S MI -MF|L M>/SQRT(2*S+1)    OBSE-168
C STARTING WITH THE MINIMUM VALUE OF L.                                 OBSE-169
   15 M1=MT(IJ)+NX(2*IJ+5,K)                                            OBSE-170
      M2=MT(IJ)+NX(2*IJ+6,K)                                            OBSE-171
      IF (M1.LT.0.OR.M2.LT.0) GO TO 85                                  OBSE-172
      M4=M1-M2                                                          OBSE-173
      M3=IABS(M4)+1                                                     OBSE-174
      M5=MT(IJ)                                                         OBSE-175
      E2=0.D0                                                           OBSE-176
      E3=1.D0                                                           OBSE-177
      C1=DFLOAT(M1+M2-MT(IJ)+1)                                         OBSE-178
      D2=0.D0                                                           OBSE-179
      C2=0.D0                                                           OBSE-180
      DO 18 M=M3,M5                                                     OBSE-181
      IF (M.EQ.M3) GO TO 16                                             OBSE-182
      E1=E2                                                             OBSE-183
      E2=E3                                                             OBSE-184
      D1=D2                                                             OBSE-185
      F1=DFLOAT(((M-1)**2-M4**2)*(MT(IJ)**2-(M-1)**2))                  OBSE-186
      F2=DFLOAT((2*M-1)*(2*M-3))                                        OBSE-187
      D2=DSQRT(F1/F2)                                                   OBSE-188
      E3=(C1*E2-D1*E1)/D2                                               OBSE-189
   16 KNX=KNX+1                                                         OBSE-190
      DO 17 JI=7,14                                                     OBSE-191
   17 NX(JI,KNX)=NX(JI,K)                                               OBSE-192
      NX(2*IJ+5,KNX)=M-1                                                OBSE-193
      NX(2*IJ+6,KNX)=M4                                                 OBSE-194
      C2=C2+E3*E3                                                       OBSE-195
   18 BM(1,KNX)=E3                                                      OBSE-196
C NORMALISATION AND SIGN GIVEN BY COEFFICIENT FOR MAXIMUM VALUE OF L.   OBSE-197
      C2=C2*DFLOAT(MT(IJ))                                              OBSE-198
      C2=1.D0/DSQRT(C2)                                                 OBSE-199
      IF (E3.LT.0.D0) C2=-C2                                            OBSE-200
      IF (MOD(NX(2*IJ+6,K),2).EQ.0) C2=-C2                              OBSE-201
      DO 19 M=M3,M5                                                     OBSE-202
      JI=KNX+M-M5                                                       OBSE-203
      BM(3,JI)=BM(3,K)*BM(1,JI)*C2                                      OBSE-204
   19 BM(1,JI)=BM(1,K)*BM(1,JI)*C2                                      OBSE-205
      BM(3,K)=0.D0                                                      OBSE-206
      BM(1,K)=0.D0                                                      OBSE-207
      K2=K+1                                                            OBSE-208
      GO TO 12                                                          OBSE-209
   20 IJ=IJ+1                                                           OBSE-210
      IF (IJ.LE.4) GO TO 11                                             OBSE-211
      K2=K1                                                             OBSE-212
      LT2=.FALSE.                                                       OBSE-213
      IF (LT1) GO TO 22                                                 OBSE-214
C IF AXIS OF QUANTIFICATION NOT IN VERTICAL PLANE, TAKE IMAGINARY       OBSE-215
C AMPLITUDE FOR PURE IMAGINARY TENSORS.                                 OBSE-216
      DO 21 K=K1,KNX                                                    OBSE-217
      M1=NX(7,K)+NX(9,K)+NX(11,K)+NX(13,K)                              OBSE-218
      IF (MOD(M1,2).NE.0) BM(3,K)=BM(1,K)                               OBSE-219
   21 CONTINUE                                                          OBSE-220
   22 IF (LT1) GO TO 35                                                 OBSE-221
C THE FIRST NON ZERO MAGNETIC QUANTUM NUMBER MUST BE POSITIVE.          OBSE-222
      DO 26 K4=K1,KNX                                                   OBSE-223
      BM(1,K4)=0.D0                                                     OBSE-224
      DO 23 J=8,14,2                                                    OBSE-225
      IF (NX(J,K4)) 24 , 23 , 26                                        OBSE-226
   23 CONTINUE                                                          OBSE-227
      GO TO 26                                                          OBSE-228
   24 IK=0                                                              OBSE-229
      DO 25 J=8,14,2                                                    OBSE-230
      IK=IK+NX(J-1,K4)+NX(J,K4)                                         OBSE-231
   25 NX(J,K4)=-NX(J,K4)                                                OBSE-232
      IF (2*(IK/2).NE.IK) BM(3,K4)=-BM(3,K4)                            OBSE-233
   26 CONTINUE                                                          OBSE-234
      KLT=2                                                             OBSE-235
C REDUCTION OF THE DESCRIPTION:                                         OBSE-236
C FOR CHANGE INTO TENSORS IF KLT=1, FOR HERE IF KLT=2,                  OBSE-237
C FOR CHANGE TO AXIS IN THE REACTION PLANE IF KLT=3.                    OBSE-238
   27 IF (K1.GT.KNX) GO TO 86                                           OBSE-239
      DO 30 K4=K1,KNX                                                   OBSE-240
      IF (DABS(BM(3,K4))+DABS(BM(1,K4)).LT..1D-6) GO TO 32              OBSE-241
      K5=K4-1                                                           OBSE-242
      IF (K1.GT.K5) GO TO 30                                            OBSE-243
      DO 29 J=K1,K5                                                     OBSE-244
      DO 28 I=7,14                                                      OBSE-245
      IF (NX(I,J).NE.NX(I,K4)) GO TO 29                                 OBSE-246
   28 CONTINUE                                                          OBSE-247
      GO TO 31                                                          OBSE-248
   29 CONTINUE                                                          OBSE-249
   30 CONTINUE                                                          OBSE-250
      GO TO ( 20 , 35 , 50 ),KLT                                        OBSE-251
   31 BM(3,J)=BM(3,J)+BM(3,K4)                                          OBSE-252
      BM(1,J)=BM(1,J)+BM(1,K4)                                          OBSE-253
   32 KNX=KNX-1                                                         OBSE-254
      IF (KNX.LT.K4) GO TO 27                                           OBSE-255
      DO 34 K=K4,KNX                                                    OBSE-256
      DO 33 J=7,14                                                      OBSE-257
   33 NX(J,K)=NX(J,K+1)                                                 OBSE-258
      BM(1,K)=BM(1,K+1)                                                 OBSE-259
   34 BM(3,K)=BM(3,K+1)                                                 OBSE-260
      GO TO 27                                                          OBSE-261
C OUTPUT OF THE DESCRIPTION.                                            OBSE-262
   35 K3=MIN0(K2+5,KNX)                                                 OBSE-263
      WRITE (MW,1008) ((NX(I,J),I=11,14),J=K2,K3)                       OBSE-264
      WRITE (MW,1012) (BM(3,I),I=K2,K3)                                 OBSE-265
      WRITE (MW,1010) ((NX(I,J),I=7,10),J=K2,K3)                        OBSE-266
      IF (LT1) WRITE (MW,1011) (BM(1,I),I=K2,K3)                        OBSE-267
      K2=K3+1                                                           OBSE-268
      IF (K2.LE.KNX) GO TO 35                                           OBSE-269
      DO 37 K=K1,KNX                                                    OBSE-270
      DO 36 I=8,14,2                                                    OBSE-271
      IF (IABS(NX(I,K)).GT.NX(I-1,K).OR.NX(I-1,K).LT.0) GO TO 87        OBSE-272
   36 CONTINUE                                                          OBSE-273
   37 CONTINUE                                                          OBSE-274
      IF (.NOT.LT1) GO TO 54                                            OBSE-275
C CHANGE FROM VERTICAL AXIS OF QUANTIFICATION TO HELICITY DESCRIPTION   OBSE-276
C BY THE ROTATION R(PI/2,PI/2,PI/2).                                    OBSE-277
      DO 39 I=K1,KNX                                                    OBSE-278
      M1=IABS(NX(8,I))+IABS(NX(10,I))+IABS(NX(12,I))+IABS(NX(14,I))     OBSE-279
      M2=M1/2                                                           OBSE-280
      IF (2*M2.NE.M1) GO TO 88                                          OBSE-281
      IF (2*(M2/2).EQ.M2) GO TO 38                                      OBSE-282
      BM(3,I)=-BM(3,I)                                                  OBSE-283
      BM(1,I)=-BM(1,I)                                                  OBSE-284
   38 IF (M1.EQ.0) BM(1,I)=0.D0                                         OBSE-285
   39 CONTINUE                                                          OBSE-286
      KLT=3                                                             OBSE-287
C PSEUDO-LOOP ON IJ TO 50.                                              OBSE-288
      IJ=7                                                              OBSE-289
   40 J1=KNX                                                            OBSE-290
      DO 47 K3=K1,KNX                                                   OBSE-291
      K=NX(IJ,K3)                                                       OBSE-292
      N1=K+NX(IJ+1,K3)+1                                                OBSE-293
      N=2*K+1                                                           OBSE-294
      E3=1.D0                                                           OBSE-295
      IF (K.EQ.0) GO TO 42                                              OBSE-296
C ROTATION MATRIX ELEMENTS FOR PI/2.                                    OBSE-297
      DO 41 I=1,K                                                       OBSE-298
   41 E3=E3*.5D0                                                        OBSE-299
   42 C1=0.D0                                                           OBSE-300
      E2=0.D0                                                           OBSE-301
      FJ=DFLOAT(K)                                                      OBSE-302
      FS=-FJ                                                            OBSE-303
      DO 43 I=1,N1                                                      OBSE-304
      IF (I.EQ.1) GO TO 43                                              OBSE-305
      C2=C1                                                             OBSE-306
      C1=DSQRT(DFLOAT((I-1)*(1+N-I)))                                   OBSE-307
      E1=E2                                                             OBSE-308
      E2=E3                                                             OBSE-309
      E3=(2.D0*FJ*E2-E1*C2)/C1                                          OBSE-310
      FS=FS+1.D0                                                        OBSE-311
   43 F2=0.D0                                                           OBSE-312
      F3=E3                                                             OBSE-313
      D1=0.D0                                                           OBSE-314
      DO 46 J=1,N                                                       OBSE-315
      IF (J.EQ.1) GO TO 44                                              OBSE-316
      D2=D1                                                             OBSE-317
      D1=DSQRT(DFLOAT((J-1)*(1+N-J)))                                   OBSE-318
      F1=F2                                                             OBSE-319
      F2=F3                                                             OBSE-320
      F3=(2.D0*FS*F2-F1*D2)/D1                                          OBSE-321
   44 J1=J1+1                                                           OBSE-322
      IF (9*J1.GT.JTN) CALL MEMO('OBSE',JTN,9*J1)                       OBSE-323
      DO 45 L=7,14                                                      OBSE-324
   45 NX(L,J1)=NX(L,K3)                                                 OBSE-325
      BM(3,J1)=F3*BM(3,K3)                                              OBSE-326
      BM(1,J1)=F3*BM(1,K3)                                              OBSE-327
   46 NX(IJ+1,J1)=J-1-K                                                 OBSE-328
   47 CONTINUE                                                          OBSE-329
C REDUCTION OF THE DESCRIPTION.                                         OBSE-330
      J2=K1-1                                                           OBSE-331
      J3=KNX+1                                                          OBSE-332
      DO 49 J4=J3,J1                                                    OBSE-333
      J2=J2+1                                                           OBSE-334
      DO 48 L=7,14                                                      OBSE-335
   48 NX(L,J2)=NX(L,J4)                                                 OBSE-336
      BM(3,J2)=BM(3,J4)                                                 OBSE-337
   49 BM(1,J2)=BM(1,J4)                                                 OBSE-338
      KNX=J2                                                            OBSE-339
      GO TO 27                                                          OBSE-340
   50 IJ=IJ+2                                                           OBSE-341
      IF (IJ.LE.13) GO TO 40                                            OBSE-342
      DO 53 J1=K1,KNX                                                   OBSE-343
      M1=NX(7,J1)+NX(9,J1)+NX(11,J1)+NX(13,J1)                          OBSE-344
      M2=NX(8,J1)+NX(10,J1)+NX(12,J1)+NX(14,J1)+4*M1                    OBSE-345
      M3=M2/2                                                           OBSE-346
      IF (2*(M1/2).EQ.M1) GO TO 51                                      OBSE-347
      IF (2*M3.NE.M2) GO TO 52                                          OBSE-348
      BM(3,J1)=BM(1,J1)                                                 OBSE-349
      GO TO 52                                                          OBSE-350
   51 IF (2*M3.EQ.M2) GO TO 52                                          OBSE-351
      BM(3,J1)=-BM(1,J1)                                                OBSE-352
   52 IF (2*(M3/2).NE.M3) BM(3,J1)=-BM(3,J1)                            OBSE-353
   53 CONTINUE                                                          OBSE-354
      LT1=.FALSE.                                                       OBSE-355
      WRITE (MW,1013)                                                   OBSE-356
      K=KNX+1-K1                                                        OBSE-357
      GO TO 9                                                           OBSE-358
C STORAGE OF THE DESCRIPTION.                                           OBSE-359
   54 L1=0                                                              OBSE-360
      ITX(1,KIT)=K1                                                     OBSE-361
      ITX(2,KIT)=KNX                                                    OBSE-362
      DO 55 I=1,NT                                                      OBSE-363
      IF (MF(2,I).NE.-KX) GO TO 55                                      OBSE-364
      MF(2,I)=KZ                                                        OBSE-365
      MF(6,I)=LT3                                                       OBSE-366
      L1=L1+1                                                           OBSE-367
   55 CONTINUE                                                          OBSE-368
      IF (L1.EQ.0) WRITE (MW,1014) KX                                   OBSE-369
      GO TO 6                                                           OBSE-370
C COMPUTATION OF ALL THE INDICATIONS NEEDED FOR THE OBSERVABLES.        OBSE-371
C (BEGINNING AND END OF DO LOOPS,GEOMETRICAL COEFFICIENTS ...)          OBSE-372
   56 DO 57 I=1,NT                                                      OBSE-373
      I1=MF(2,I)+1                                                      OBSE-374
      MF(3,I)=ITX(1,I1)                                                 OBSE-375
   57 MF(4,I)=ITX(2,I1)                                                 OBSE-376
      LT1=.FALSE.                                                       OBSE-377
   58 LT1=.NOT.LT1                                                      OBSE-378
C LT1=.TRUE.  FIRST PASSAGE: NUMBER OF INFORMATIONS NEEDED FOR DO LOOPS OBSE-379
C LT1=.FALSE. SECOND ONE: COMPUTATION OF GEOMETRICAL COEFFICIENTS WHICH OBSE-380
C ARE STORED AFTER THE INDICATIONS FOR DO LOOPS.                        OBSE-381
      ICAL=0                                                            OBSE-382
      NCAL=0                                                            OBSE-383
      DO 78 I=1,NT                                                      OBSE-384
      IV=MF(1,I)                                                        OBSE-385
      MT(3)=IPI(2,IV)                                                   OBSE-386
      MT(4)=IPI(3,IV)                                                   OBSE-387
      IF (MF(2,I).LE.0.OR.MF(2,I).EQ.19) GO TO 78                       OBSE-388
      IF (MF(2,I).NE.1) GO TO 59                                        OBSE-389
      IF (IV.NE.1) GO TO 89                                             OBSE-390
      GO TO 78                                                          OBSE-391
   59 IF ((MF(2,I).EQ.2.AND.MT(1).LE.3).OR.(MF(2,I).EQ.3.AND.MT(3).LE.3)OBSE-392
     1) GO TO 62                                                        OBSE-393
      IF (I.EQ.1) GO TO 63                                              OBSE-394
      I1=I-1                                                            OBSE-395
      DO 60 J=1,I1                                                      OBSE-396
      IF (MF(1,I).EQ.MF(1,J).AND.MF(2,I).EQ.MF(2,J)) GO TO 61           OBSE-397
   60 CONTINUE                                                          OBSE-398
      GO TO 63                                                          OBSE-399
   61 MF(3,I)=MF(3,J)                                                   OBSE-400
      MF(4,I)=MF(4,J)                                                   OBSE-401
      GO TO 78                                                          OBSE-402
   62 MF(2,I)=-MF(2,I)                                                  OBSE-403
      IF ((MF(2,I).EQ.-2.AND.MT(1).EQ.0).OR.(MF(2,I).EQ.-3.AND.MT(3).EQ.OBSE-404
     10)) GO TO 90                                                      OBSE-405
      GO TO 78                                                          OBSE-406
   63 I1=MF(3,I)                                                        OBSE-407
      I2=MF(4,I)                                                        OBSE-408
      DO 77 KI=I1,I2                                                    OBSE-409
      KT=0                                                              OBSE-410
      DO 73 L=1,4                                                       OBSE-411
      K=NX(2*L+5,KI)                                                    OBSE-412
      NM=NX(2*L+6,KI)                                                   OBSE-413
      KT=KT+K                                                           OBSE-414
      IF (K.GT.MT(L)) GO TO 91                                          OBSE-415
      IF (ICAL.EQ.0) GO TO 65                                           OBSE-416
      DO 64 N=1,ICAL                                                    OBSE-417
      IF (MT(L).EQ.NX(15,N).AND.K.EQ.NX(16,N).AND.NM.EQ.NX(17,N)) GO TO OBSE-418
     1 71                                                               OBSE-419
   64 CONTINUE                                                          OBSE-420
   65 ICAL=ICAL+1                                                       OBSE-421
      IF (9*ICAL.GT.JTN) CALL MEMO('OBSE',JTN,9*ICAL)                   OBSE-422
      NX(15,ICAL)=MT(L)                                                 OBSE-423
      NX(16,ICAL)=K                                                     OBSE-424
      NX(17,ICAL)=NM                                                    OBSE-425
      NX(18,ICAL)=JCAL                                                  OBSE-426
      IF (LT1) GO TO 73                                                 OBSE-427
      IL=MT(L)                                                          OBSE-428
      NX(1,ICAL)=MAX0(1,1+NM)                                           OBSE-429
      NX(2,ICAL)=MIN0(IL,IL+NM)                                         OBSE-430
      NX(3,ICAL)=NM                                                     OBSE-431
      NX(4,ICAL)=JCAL-NX(1,ICAL)+1                                      OBSE-432
      IF (NX(16,ICAL).EQ.0) GO TO 70                                    OBSE-433
C RECURRENCE COMPUTATION OF GEOMETRICAL COEFFICIENTS.                   OBSE-434
      F3=DSQRT(DFLOAT(2*K+1))                                           OBSE-435
      DO 66 J=1,K                                                       OBSE-436
   66 F3=-F3*DSQRT(DFLOAT(IL-J)/DFLOAT(IL+J))                           OBSE-437
      JNM=NM                                                            OBSE-438
      INM=IABS(JNM)                                                     OBSE-439
      IF (JNM.EQ.0) GO TO 68                                            OBSE-440
      IF (INM.NE.JNM.AND.2*(INM/2).NE.INM) F3=-F3                       OBSE-441
      DO 67 J=1,INM                                                     OBSE-442
   67 F3=F3*DSQRT(DFLOAT((K+J)*(K-J+1))/DFLOAT(J*(IL-J)))               OBSE-443
   68 JCAL=JCAL+1                                                       OBSE-444
      IS=NX(2,ICAL)-NX(1,ICAL)                                          OBSE-445
      IF (JCAL+IS.GT.JTS) CALL MEMO('OBSE',JTS,JCAL+IS)                 OBSE-446
      DM(JCAL+KX)=F3                                                    OBSE-447
      IF (IS.LT.1) GO TO 70                                             OBSE-448
      F2=0.D0                                                           OBSE-449
      D1=K*(K+1)-(IL+1)*(INM-1)                                         OBSE-450
      C2=0.D0                                                           OBSE-451
      DO 69 J=1,IS                                                      OBSE-452
      C1=C2                                                             OBSE-453
      C2=DSQRT(DFLOAT(J*(IL-J)*(J+INM)*(IL-J-INM)))                     OBSE-454
      D1=D1+DFLOAT(2*(2*J-2-IL+INM))                                    OBSE-455
      F1=F2                                                             OBSE-456
      F2=F3                                                             OBSE-457
      F3=-(D1*F2+C1*F1)/C2                                              OBSE-458
      JCAL=JCAL+1                                                       OBSE-459
      DM(JCAL+KX)=F3                                                    OBSE-460
   69 CONTINUE                                                          OBSE-461
   70 N=ICAL                                                            OBSE-462
   71 IF (LT1) GO TO 73                                                 OBSE-463
      DO 72 MA=1,4                                                      OBSE-464
      MC=20*NCAL+4*L+MA-4                                               OBSE-465
   72 NY(MC+2*KX)=NX(MA,N)                                              OBSE-466
      LZ(L)=NX(16,N).NE.0                                               OBSE-467
   73 CONTINUE                                                          OBSE-468
      IF (LT1) GO TO 76                                                 OBSE-469
      LZ(5)=KT.EQ.0                                                     OBSE-470
      LZ(6)=2*(KT/2).NE.KT                                              OBSE-471
      NZ(1)=0                                                           OBSE-472
      NZ(2)=0                                                           OBSE-473
      IF (LZ(4)) NZ(2)=MF(6,I)                                          OBSE-474
      IF (LZ(3).OR.LZ(4)) NZ(1)=MF(6,I)                                 OBSE-475
      LO(126)=LO(126).OR.NZ(1).NE.0                                     OBSE-476
      IA1=0                                                             OBSE-477
      IA2=1                                                             OBSE-478
      DO 74 LI=1,6                                                      OBSE-479
      IF (LZ(LI)) IA1=IA1+IA2                                           OBSE-480
   74 IA2=2*IA2                                                         OBSE-481
      NY(20*NCAL+17+2*KX)=IA1                                           OBSE-482
      NY(20*NCAL+18+2*KX)=NZ(1)+1000*NZ(2)                              OBSE-483
      DM(10*NCAL+10+KX)=BM(3,KI)                                        OBSE-484
      DO 75 LI=2,4,2                                                    OBSE-485
      IF (.NOT.LZ(LI)) GO TO 75                                         OBSE-486
      IM=NY(20*NCAL+4*LI-3+2*KX)                                        OBSE-487
      NY(20*NCAL+4*LI-3+2*KX)=MT(LI)-NY(20*NCAL+4*LI-2+2*KX)+1          OBSE-488
      NY(20*NCAL+4*LI-2+2*KX)=MT(LI)-IM+1                               OBSE-489
      NY(20*NCAL+4*LI-1+2*KX)=-NY(20*NCAL+4*LI-1+2*KX)                  OBSE-490
      NY(20*NCAL+4*LI+2*KX)=NY(20*NCAL+4*LI+2*KX)+IM+NY(20*NCAL+4*LI-2+2OBSE-491
     1*KX)                                                              OBSE-492
   75 CONTINUE                                                          OBSE-493
   76 NCAL=NCAL+1                                                       OBSE-494
   77 CONTINUE                                                          OBSE-495
      IF (LT1) GO TO 78                                                 OBSE-496
      MF(4,I)=NCAL                                                      OBSE-497
      MF(3,I)=NCAL+I1-I2                                                OBSE-498
   78 CONTINUE                                                          OBSE-499
      IF (.NOT.LT1) GO TO 79                                            OBSE-500
      KX=9*MAX0(KNX,ICAL)                                               OBSE-501
      JCAL=10*NCAL                                                      OBSE-502
      JTS=JTN-18*(KX-1)                                                 OBSE-503
      IF (JCAL.GT.JTS) CALL MEMO('OBSE',JTS,JCAL)                       OBSE-504
      GO TO 58                                                          OBSE-505
C STORAGE OF LEGENDS AND COPY OF RESULTS.                               OBSE-506
   79 DO 81 I=1,NT                                                      OBSE-507
      I5=IABS(MF(2,I))+1                                                OBSE-508
      DO 80 J=3,7                                                       OBSE-509
   80 CMF(J+3,I)=ITZ(J,I5)                                              OBSE-510
   81 CONTINUE                                                          OBSE-511
      JN=2*JCAL                                                         OBSE-512
      DO 82 I=1,JN                                                      OBSE-513
   82 AM(I)=CBM(I+2*KX)                                                 OBSE-514
      RETURN                                                            OBSE-515
   83 WRITE (MW,1015) I,MF(2,I)                                         OBSE-516
      GO TO 92                                                          OBSE-517
   84 WRITE (MW,1016) I,MF(2,I)                                         OBSE-518
      GO TO 92                                                          OBSE-519
   85 WRITE (MW,1017) NX(2*IJ+5,K),NX(2*IJ+6,K),MT(IJ)                  OBSE-520
      GO TO 92                                                          OBSE-521
   86 WRITE (MW,1018)                                                   OBSE-522
      GO TO 92                                                          OBSE-523
   87 WRITE (MW,1019) NX(2*I-1,K1),NX(2*I,K1)                           OBSE-524
      GO TO 92                                                          OBSE-525
   88 WRITE (MW,1020)                                                   OBSE-526
      GO TO 92                                                          OBSE-527
   89 WRITE (MW,1021) MF(1,I)                                           OBSE-528
      GO TO 92                                                          OBSE-529
   90 WRITE (MW,1022) MF(1,I)                                           OBSE-530
      GO TO 92                                                          OBSE-531
   91 WRITE (MW,1023) MT(L),MF(2,I),L,K                                 OBSE-532
   92 WRITE (MW,1024)                                                   OBSE-533
      STOP                                                              OBSE-534
 1000 FORMAT (2L1,I1,I2,I5,5A4)                                         OBSE-535
 1001 FORMAT (8I5)                                                      OBSE-536
 1002 FORMAT (6F10.5)                                                   OBSE-537
 1003 FORMAT (/' OBSERVABLE',I3,' LABELLED',5X,5A4,10X,I3,' COMPONENTS')OBSE-538
 1004 FORMAT (' DEFINED WITH AN AXIS PERPENDICULAR TO THE REACTION PLANEOBSE-539
     1')                                                                OBSE-540
 1005 FORMAT (' DEFINED IN THE LABORATORY SYSTEM')                      OBSE-541
 1006 FORMAT (' DEFINED WITH RESPECT TO THE INCIDENT BEAM')             OBSE-542
 1007 FORMAT (' NOT COMPLETELY DEFINED BY TENSORS')                     OBSE-543
 1008 FORMAT (/6(11X,4I2,2X))                                           OBSE-544
 1009 FORMAT (6(' +',F8.4,'*M',9X))                                     OBSE-545
 1010 FORMAT (6(11X,4I2,2X))                                            OBSE-546
 1011 FORMAT (' IMAGINARY PARTS'/6(1X,F8.5,12X))                        OBSE-547
 1012 FORMAT (6(' +',F8.4,'*A',9X))                                     OBSE-548
 1013 FORMAT (//' AFTER TRANSFORMATION:'/)                              OBSE-549
 1014 FORMAT (' THE OBSERVABLE READ WITH NUMBER',I4,' IS NOT USED.')    OBSE-550
 1015 FORMAT (' THE',I4,'TH OBSERVABLE, OF KIND',I2,' MUST BE FOR EXPERIOBSE-551
     1MENTAL DATA.')                                                    OBSE-552
 1016 FORMAT (' THE',I4,'TH OBSERVABLE, OF KIND',I2,' IS NOT DEFINED.') OBSE-553
 1017 FORMAT (' NON TENSOR INDICATIONS',I4,' AND',I4,' INCORRECT FOR (2*OBSE-554
     1S+1) =',I3)                                                       OBSE-555
 1018 FORMAT (' ZERO OBSERVABLE.')                                      OBSE-556
 1019 FORMAT (/' TOO LARGE MAGNETIC QUANTUM NUMBER OR NEGATIVE MULTIPOLAOBSE-557
     1RITY',2I6)                                                        OBSE-558
 1020 FORMAT (' THE SUM OF MAGNETIC QUANTUM NUMBERS IS ODD FOR ONE COMPOOBSE-559
     1NENT.')                                                           OBSE-560
 1021 FORMAT (' NO CROSS SECTION DIVIDED BY RUTHERFORD''S FOR THE INELASOBSE-561
     1TIC CHANNEL',I3)                                                  OBSE-562
 1022 FORMAT (' NO POLARISATION FOR A ZERO SPIN IN THE CHANNEL',I3)     OBSE-563
 1023 FORMAT (5X,I5,' IS A TOO SMALL SPIN IN CHANNEL',I4,' AND PARTICLE'OBSE-564
     1,I4,' FOR A POLARISATION OF TENSOR ORDER',I4)                     OBSE-565
 1024 FORMAT (//' IN OBSE  ...  STOP  ...')                             OBSE-566
      END                                                               OBSE-567
