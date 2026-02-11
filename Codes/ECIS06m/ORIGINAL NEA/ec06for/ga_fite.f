C 28/06/06                                                      ECIS06  FITE-000
      SUBROUTINE FITE(KE,M,N,F,X,E,W,NW,IZ,WA)                          FITE-001
C CHI-SQUARE MINIMISING SUBROUTINE WRITTEN BY G. SCHWEIMER IN KARLSRUHE.FITE-002
C SIMPLIFIED FOR THE USE OF ECIS (NO GRADIENT).                         FITE-003
C SOLVES THE NONLINEAR LEAST SQUARES PROBLEM,                           FITE-004
C USING A LEAST SQUARES INTERPOLATION BETWEEN VARIABLES AND FUNCTIONS.  FITE-005
C CALLED SUBROUTINES: FIT2 /.TRUE./ (LINEAR LEAST SQUARES PROBLEM).     FITE-006
C                     FIT2 /.FALSE./(INVERSION OF A(TRANSPOSED)*A).     FITE-007
C                     FIT1(ONE DIMENSIONAL MINIMUM SEARCH).             FITE-008
C INPUT:     KE:      0 FOR THE FIRST CALL, 1 FOR THE FOLLOWING ONES.   FITE-009
C            M:       NUMBER OF FUNCTIONS, M GE N.                      FITE-010
C            N:       NUMBER OF VARIABLES, N GE 1.                      FITE-011
C            F:       FUNCTION VALUES AT THE POINT X.                   FITE-012
C            X:       STARTING VALUES OF THE VARIABLES.                 FITE-013
C            E:       ABSOLUTE SEARCH ACCURACIES FOR THE VARIABLES.     FITE-014
C            IZ(I)    NUMBER OF POINTS TO BE REMEMBERED FOR I=1         FITE-015
C                     (AT LEAST N+1).                                   FITE-016
C                     MAXIMUM NUMBER OF FUNCTION EVALUATIONS FOR I=2.   FITE-017
C OUTPUT:    KE:      ERROR CODE:                                       FITE-018
C                     KE=0: WITHOUT ERRORS.                             FITE-019
C                     KE=2: USER INTERRUPT.                             FITE-020
C                     KE=3: MAXIMUM NUMBER OF FUNCTION EVALUATIONS.     FITE-021
C                     KE=4: ROUNDING ERRORS.                            FITE-022
C                     KE=5: THE FUNCTIONS DO NOT DEPEND ON X(IZ(4)).    FITE-023
C                     KE=6: USELESS VARIABLES IN THE PREPARATORY CALLS, FITE-024
C                           THE LABELS OF THE VARIABLES ARE IZ(3),IZ(4).FITE-025
C                     KE=7: M LT N.                                     FITE-026
C            W(I):    FOR I<N, STANDARD ERROR OF THE VARIABLE I (FOR A  FITE-027
C                     VARIANCE AT BEST FIT EQUAL TO DEGREE OF FREEDOM). FITE-028
C                     FOR I+N, ERROR ENHANCEMENT (MANY/SINGLE VARIABLE  FITE-029
C                     RESULT).                                          FITE-030
C                     FOR N+N+I+(J*(J-1))/2, ERROR CORRELATION BETWEEN  FITE-031
C                     X(I) AND X(J), I<J.                               FITE-032
C            IZ(I)    NUMBER OF FUNCTION EVALUATIONS FOR I=3.           FITE-033
C                     FOR I=4, NUMBER OF DEGREES OF FREEDOM IF KE=0,2,3 FITE-034
C                     OR 4.  PLACE OF USELESS VARIABLE IF KE=5 OR 6.    FITE-035
C                     EVENTUALLY INDICATIONS FOR ONE DIMENSIONAL SEARCH FITE-036
C                     IF KE=1.                                          FITE-037
C WORKING AREAS:                                                        FITE-038
C            IZ(I):   FOR I=5 TO 4+IZ(1), LABELS OF THE QUADRATIC SUMS. FITE-039
C            W(I):    OF LENGTH MAX(14+N+K*(M+N+1),(N*(N+5))/2)         FITE-040
C                     EVENTUALLY, INDICATIONS FOR ONE DIMENSIONAL       FITE-041
C                     SEARCH. FROM MAX(15,N+1),K SETS OF (FUNCTIONS +   FITE-042
C                     VARIABLES + CHI2). SEE ABOVE THE ADDRESSES IN W   FITE-043
C                     AFTER THE SEARCH                                  FITE-044
C            NW:      EQUIVALENT BY CALL WITH W.                        FITE-045
C            WA(I):   OF LENGTH (M+1)*(N+1)+(N*(N+1))/2. FROM I=1,      FITE-046
C                     GRADIENT OR MATRIX A FOLLOWED BY CENTRAL FUNCTIONSFITE-047
C                     FROM I=1+M*(N+1), CENTRAL VARIABLES,              FITE-048
C                     FROM I=(N+1)*(M+1), MATRIX D.                     FITE-049
C                                                                       FITE-050
C THE WORKING FIELDS IZ AND W CONTAIN ALL INFORMATION TO CONTINUE       FITE-051
C THE SEARCH. THIS ALLOWS A SEARCH WITHIN ANOTHER SEARCH JUST CHANGING  FITE-052
C THE WORKING FIELDS.                                                   FITE-053
C                                                                       FITE-054
C FOR THE COMMON  /DCHI2/  SEE CALX.                                    FITE-055
C                                                                       FITE-056
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCHI2/:                     FITE-057
C  CHI2:      CHI-SQUARE COMPUTED IN SUBROUTINE RESU.                   FITE-058
C  YY(1):     STEP SIZE IN THE SEARCH.                                  FITE-059
C  YY(2):     HALF OF THE SUCCESS MULTIPLICATIVE FACTOR OF THE STEP.    FITE-060
C             THE STEP SIZE IS DIVIDED BY 4 TIMES ITS SQUARE AFTER      FITE-061
C             UNSUCCESS.                                                FITE-062
C  YY(3):     CHANGED X VALUE IN PREPARATORY CALLS,                     FITE-063
C             RATIO STEP/DISTANCE OF MINIMUM AFTER PREPARATORY CALLS    FITE-064
C             OR 0. FOR RANDOM PREDICTION,                              FITE-065
C             ERROR RENORMALISATION FACTOR WHEN KE=0,2 OR 3.(SQUARE ROOTFITE-066
C             OF VARIANCE AT BEST FIT DIVIDED BY DEGREE OF FREEDOM).    FITE-067
C   DEFINED:  CHI2,YY.                                                  FITE-068
C   USED:     CHI2,YY.                                                  FITE-069
C                                                                       FITE-070
C***********************************************************************FITE-071
      IMPLICIT REAL*8 (A-H,O-Z)                                         FITE-072
      DIMENSION F(*),X(*),E(*),W(*),NW(*),IZ(*),WA(*)                   FITE-073
      COMMON /DCHI2/ CHI2,CHI2M,YY(3)                                   FITE-074
      DATA EPS /1.D-5/                                                  FITE-075
      JD=(M+1)*(N+1)-1                                                  FITE-076
      JS=14+N                                                           FITE-077
      LM=M+N+1                                                          FITE-078
      K=IZ(1)                                                           FITE-079
      IF (KE.NE.0) GO TO 2                                              FITE-080
      IF (M.LT.N) GO TO 39                                              FITE-081
      IZ(3)=1                                                           FITE-082
      IZ(4)=0                                                           FITE-083
      DO 1 L=1,K                                                        FITE-084
      IZ(L+4)=1+K-L                                                     FITE-085
    1 W(JS+LM*L)=7.D35                                                  FITE-086
      IF (YY(2).LT.1.D0) YY(2)=1.D0                                     FITE-087
      KE=1                                                              FITE-088
    2 JM=JS+LM*IZ(5)-LM                                                 FITE-089
      J3=M*N                                                            FITE-090
      KV=K                                                              FITE-091
      IF (CHI2.LE.0.D0) GO TO 38                                        FITE-092
C ROW OF MATRIX S TO BE REPLACED BY NEW VALUES                          FITE-093
      L=IZ(K+4)                                                         FITE-094
      IF (W(JS+LM*L).EQ.7.D35) KV=L-1                                   FITE-095
      DO 3 I=1,K                                                        FITE-096
      J1=JS+LM*IZ(I+4)                                                  FITE-097
      IF (CHI2.LT.W(J1)) GO TO 4                                        FITE-098
    3 CONTINUE                                                          FITE-099
C ONE DIMENSIONAL SEARCH IS NECESSARY                                   FITE-100
      GO TO 30                                                          FITE-101
    4 IF (I.GT.MAX0(N+1,KV)) GO TO 30                                   FITE-102
C VECTOR OF LABELS OF THE QUADRATIC SUMS                                FITE-103
      IF (KV.LT.K) KV=KV+1                                              FITE-104
      I1=K+4                                                            FITE-105
      I2=K-I                                                            FITE-106
      IF ((IZ(3).GT.N+1).AND.(I.NE.1)) YY(1)=YY(1)/(4.D0*YY(2)**3)      FITE-107
      IF (I2.EQ.0) GO TO 6                                              FITE-108
      DO 5 J=1,I2                                                       FITE-109
      I1=I1-1                                                           FITE-110
    5 IZ(I1+1)=IZ(I1)                                                   FITE-111
      IZ(I1)=L                                                          FITE-112
      JM=JS+LM*IZ(5)-LM                                                 FITE-113
C NEW ROW                                                               FITE-114
    6 J1=JS+LM*(L-1)                                                    FITE-115
      DO 7 I=1,M                                                        FITE-116
      J1=J1+1                                                           FITE-117
    7 W(J1)=F(I)                                                        FITE-118
      DO 8 I=1,N                                                        FITE-119
      J1=J1+1                                                           FITE-120
    8 W(J1)=X(I)                                                        FITE-121
      W(J1+1)=CHI2                                                      FITE-122
      IF (IZ(3).GE.IZ(2)) GO TO 43                                      FITE-123
      IF (N.EQ.1) GO TO 35                                              FITE-124
      IF (IZ(3).GT.N+1) GO TO 13                                        FITE-125
C PREPARATORY FUNCTION EVALUATIONS                                      FITE-126
      MF=IZ(3)                                                          FITE-127
      IF (MF.EQ.1) GO TO 12                                             FITE-128
C SIGNIFICANCE OF THE NEW VARIABLE                                      FITE-129
      X(MF-1)=YY(3)                                                     FITE-130
      S=0.D0                                                            FITE-131
      DO 9 I=1,M                                                        FITE-132
      T=F(I)-W(JS+I)                                                    FITE-133
    9 S=S+T*T                                                           FITE-134
      J=2                                                               FITE-135
      IF (S.LT.EPS*EPS*W(JS+LM)) GO TO 41                               FITE-136
      YY(3)=S                                                           FITE-137
      W(MF-1)=DSQRT(YY(3))                                              FITE-138
      IF (MF.LE.2) GO TO 12                                             FITE-139
C INDEPENDENCE OF THE NEW VARIABLE                                      FITE-140
      DO 11 J=3,MF                                                      FITE-141
      I2=JS+LM*(J-2)                                                    FITE-142
      S=0.D0                                                            FITE-143
      DO 10 I=1,M                                                       FITE-144
   10 S=S+(W(I2+I)-W(JS+I))*(F(I)-W(JS+I))                              FITE-145
      IF (DABS(W(MF-1)*W(J-2)-DABS(S)).LT.EPS*DABS(S)) GO TO 40         FITE-146
   11 CONTINUE                                                          FITE-147
   12 IF (MF.EQ.N+1) GO TO 13                                           FITE-148
      YY(3)=X(MF)                                                       FITE-149
      X(MF)=X(MF)+YY(1)*E(MF)                                           FITE-150
      GO TO 64                                                          FITE-151
C END OF PREPARATORY FUNCTION EVALUATIONS                               FITE-152
C SUM OF INVERSES OF THE QUADRATIC SUMS                                 FITE-153
   13 S=0.D0                                                            FITE-154
      DO 14 L=1,KV                                                      FITE-155
      T=W(JS+LM*L)                                                      FITE-156
   14 S=S+1.D0/(T*T)                                                    FITE-157
      WJA=1.D0/S                                                        FITE-158
C CENTRE OF THE VARIABLES AND FUNCTIONS                                 FITE-159
      I1=M+N                                                            FITE-160
      DO 16 I=1,I1                                                      FITE-161
      J1=JS                                                             FITE-162
      S=0.D0                                                            FITE-163
      DO 15 L=1,KV                                                      FITE-164
      T=W(J1+LM)                                                        FITE-165
      S=S+W(J1+I)/(T*T)                                                 FITE-166
   15 J1=J1+LM                                                          FITE-167
   16 WA(J3+I)=S*WJA                                                    FITE-168
      IF (KE.NE.1) GO TO 47                                             FITE-169
C THE LINEAR APPROXIMATION MEANS THAT THE DIFFERENCE WITH CENTRAL VALUE FITE-170
C F(I,K)-F(I) = SUM F(I,J)*(X(J,K)-X(J)) WHERE F(I),X(J) ARE MEAN VALUESFITE-171
C AND F(I,J) IS THE FIRST DERIVATIVE WITH RESPECT TO X(J).              FITE-172
C USING THE MATRIX D(I,J) = WEIGHTED SUM ON (X(I,K)-X(I))*(X(J,K)-X(J)) FITE-173
C AND WRITING THE CHANGE DX = D*Y, THE BEST FIT IS OBTAINED FOR         FITE-174
C ||F-A*Y||=MIN(Y) WITH A(I,J)=WEIGHTED SUM(F(I,K)-F(I)*(X(J,K)-X(J))   FITE-175
C MATRIX A                                                              FITE-176
   17 J1=0                                                              FITE-177
      DO 19 I=1,N                                                       FITE-178
      U=WA(J3+M+I)                                                      FITE-179
      DO 19 J=1,M                                                       FITE-180
      J1=J1+1                                                           FITE-181
      J2=JS                                                             FITE-182
      S=0.D0                                                            FITE-183
      T=WA(J3+J)                                                        FITE-184
      DO 18 L=1,KV                                                      FITE-185
      V=W(J2+LM)                                                        FITE-186
      S=S+(W(J2+J)-T)*(W(J2+M+I)-U)/(V*V)                               FITE-187
   18 J2=J2+LM                                                          FITE-188
   19 WA(J1)=S*WJA                                                      FITE-189
      IF (KE.NE.1) GO TO 50                                             FITE-190
C LINEAR LEAST SQUARES PROBLEM                                          FITE-191
      CALL FIT2(M,N,WA,X,NW,IR,.TRUE.)                                  FITE-192
      IF (IR) 42 , 20 , 28                                              FITE-193
C MATRIX D                                                              FITE-194
   20 J1=JD                                                             FITE-195
      JA=J3+M                                                           FITE-196
      DO 23 I=1,N                                                       FITE-197
      T=WA(JA+I)                                                        FITE-198
      DO 22 J=1,I                                                       FITE-199
      J1=J1+1                                                           FITE-200
      J2=JS+M                                                           FITE-201
      S=0.D0                                                            FITE-202
      U=WA(JA+J)                                                        FITE-203
      DO 21 L=1,KV                                                      FITE-204
      V=W(J2+N+1)                                                       FITE-205
      S=S+(W(J2+I)-T)*(W(J2+J)-U)/(V*V)                                 FITE-206
   21 J2=J2+LM                                                          FITE-207
   22 WA(J1)=S*WJA                                                      FITE-208
   23 CONTINUE                                                          FITE-209
C NEW VARIABLES                                                         FITE-210
      DO 25 I=1,N                                                       FITE-211
      I2=1                                                              FITE-212
      J1=JD+(I*I-I)/2                                                   FITE-213
      S=0.D0                                                            FITE-214
      DO 24 J=1,N                                                       FITE-215
      J1=J1+I2                                                          FITE-216
      IF (J.GE.I) I2=J                                                  FITE-217
   24 S=S+WA(J1)*X(J)                                                   FITE-218
   25 WA(I)=WA(JA+I)-S                                                  FITE-219
C TEST OF CONVERGENCE                                                   FITE-220
      A=0.D0                                                            FITE-221
      DO 26 I=1,N                                                       FITE-222
      W(I)=WA(I)-W(JM+M+I)                                              FITE-223
   26 A=DMAX1(A,DABS(W(I)/E(I)))                                        FITE-224
      YY(1)=YY(1)*YY(2)                                                 FITE-225
      IF (A.LT.1.D0.OR.YY(1).LT.1.D0) GO TO 38                          FITE-226
      YY(3)=1.D0                                                        FITE-227
C STEP SIZE LIMITATION                                                  FITE-228
      IF (A.GT.2.D0*YY(1)) YY(3)=2.D0*YY(1)/A                           FITE-229
      DO 27 I=1,N                                                       FITE-230
   27 X(I)=W(JM+M+I)+YY(3)*W(I)                                         FITE-231
      IZ(4)=0                                                           FITE-232
      YY(1)=A*YY(3)                                                     FITE-233
      GO TO 64                                                          FITE-234
C RANDOM PREDICTION                                                     FITE-235
   28 DO 29 I=1,N                                                       FITE-236
   29 X(I)=W(JM+M+I)+YY(1)*E(I)*DFLOAT(MOD(IABS(NW(JM+I)),200)-100)/100.FITE-237
     1D0                                                                FITE-238
      YY(3)=0.D0                                                        FITE-239
      GO TO 64                                                          FITE-240
C ONE DIMENSIONAL SEARCH                                                FITE-241
   30 IF (N.EQ.1) GO TO 36                                              FITE-242
      IF (IZ(3).GE.IZ(2)) GO TO 43                                      FITE-243
      IF (IZ(4).EQ.2) GO TO 32                                          FITE-244
      IZ(4)=2                                                           FITE-245
      DO 31 I=1,N                                                       FITE-246
   31 W(I+14)=X(I)-W(JM+M+I)                                            FITE-247
      NW(1)=3                                                           FITE-248
      NW(2)=20                                                          FITE-249
      W(4)=0.5D0                                                        FITE-250
      W(7)=0.D0                                                         FITE-251
      W(8)=0.D0                                                         FITE-252
      W(9)=0.D0                                                         FITE-253
      W(10)=1.D0                                                        FITE-254
      W(12)=W(JM+LM)                                                    FITE-255
      W(13)=CHI2                                                        FITE-256
      GO TO 33                                                          FITE-257
   32 W(5)=CHI2                                                         FITE-258
      CALL FIT1(KE,NW,W(4))                                             FITE-259
   33 DO 34 I=1,N                                                       FITE-260
   34 X(I)=W(JM+M+I)+W(4)*W(I+14)                                       FITE-261
      IF (KE.EQ.3) KE=2                                                 FITE-262
      IF (KE.EQ.2) GO TO 43                                             FITE-263
      KE=1                                                              FITE-264
      YY(3)=W(4)                                                        FITE-265
      GO TO 64                                                          FITE-266
C ONLY ONE VARIABLE X                                                   FITE-267
   35 IF (IZ(3).GT.1) GO TO 36                                          FITE-268
      KE=0                                                              FITE-269
      W(6)=YY(1)*E(1)                                                   FITE-270
      W(7)=E(1)                                                         FITE-271
      W(8)=0.D0                                                         FITE-272
   36 NW(2)=IZ(2)                                                       FITE-273
      W(4)=X(1)                                                         FITE-274
      W(5)=CHI2                                                         FITE-275
      CALL FIT1(KE,NW,W(4))                                             FITE-276
      IZ(4)=2                                                           FITE-277
      X(1)=W(4)                                                         FITE-278
      IF (KE.EQ.1) GO TO 64                                             FITE-279
      IF (KE.GT.0) KE=KE+1                                              FITE-280
      YY(3)=0.D0                                                        FITE-281
      W(1)=0.D0                                                         FITE-282
      DO 37 J=1,M                                                       FITE-283
   37 F(J)=W(JM+I)                                                      FITE-284
      CHI2=W(JM+LM)                                                     FITE-285
      X(1)=W(JM+LM-1)                                                   FITE-286
      IF (NW(2).NE.0) GO TO 63                                          FITE-287
      W(1)=DSQRT(DABS((W(9)-W(11))/((W(12)-W(13))/(W(9)-W(10))-(W(13)-W(FITE-288
     114))/(W(10)-W(11)))))                                             FITE-289
      W(2)=1.D0                                                         FITE-290
      W(3)=1.D0                                                         FITE-291
      GO TO 60                                                          FITE-292
C END OF SEARCH                                                         FITE-293
   38 KE=0                                                              FITE-294
      IF (CHI2.EQ.0.D0.OR.IZ(2).LT.0) GO TO 64                          FITE-295
      GO TO 44                                                          FITE-296
C ERROR CODE DEFINITION                                                 FITE-297
   39 KE=KE+1                                                           FITE-298
   40 KE=KE+1                                                           FITE-299
   41 KE=KE+1                                                           FITE-300
   42 KE=KE+1                                                           FITE-301
   43 KE=KE+1                                                           FITE-302
      KE=KE+1                                                           FITE-303
C RESTORE OPTIMUM VALUES TO X AND F                                     FITE-304
   44 DO 45 I=1,M                                                       FITE-305
   45 F(I)=W(JM+I)                                                      FITE-306
      DO 46 I=1,N                                                       FITE-307
      X(I)=W(JM+M+I)                                                    FITE-308
   46 W(I)=0.D0                                                         FITE-309
      YY(3)=0.D0                                                        FITE-310
      CHI2=W(JM+LM)                                                     FITE-311
      IF (KE*(KE-3).NE.0.OR.(KE.EQ.3.AND.((YY(3).EQ.0.D0.AND.IZ(3).LE.N)FITE-312
     1))) GO TO 63                                                      FITE-313
C COMPUTATION OF THE ERRORS OF THE VARIABLES - RESTORE MATRIX G         FITE-314
      KV=MIN0(K,IZ(3))                                                  FITE-315
      GO TO 13                                                          FITE-316
C INVERSE OF MATRIX D                                                   FITE-317
   47 T=DSQRT(WJA)                                                      FITE-318
      J1=0                                                              FITE-319
      DO 49 I=1,N                                                       FITE-320
      S=WA(J3+M+I)                                                      FITE-321
      J2=JS+I-LM+M                                                      FITE-322
      DO 48 L=1,KV                                                      FITE-323
      J1=J1+1                                                           FITE-324
   48 WA(J1)=T*(W(J2+L*LM)-S)/W(JS+L*LM)                                FITE-325
   49 CONTINUE                                                          FITE-326
      CALL FIT2(KV,N,WA,WA(JD+1),NW,IR,.FALSE.)                         FITE-327
      IF (IR) 63 , 17 , 63                                              FITE-328
C MATRIX G = A*INVERSE OF D                                             FITE-329
   50 DO 54 L=1,M                                                       FITE-330
      J1=L-M                                                            FITE-331
      DO 52 I=1,N                                                       FITE-332
      I1=JD+(I*I-I)/2                                                   FITE-333
      I2=1                                                              FITE-334
      S=0.D0                                                            FITE-335
      DO 51 J=1,N                                                       FITE-336
      I1=I1+I2                                                          FITE-337
      IF (J.GE.I) I2=J                                                  FITE-338
   51 S=S+WA(I1)*WA(J1+J*M)                                             FITE-339
   52 W(I)=S                                                            FITE-340
      DO 53 J=1,N                                                       FITE-341
   53 WA(J1+J*M)=W(J)                                                   FITE-342
   54 CONTINUE                                                          FITE-343
C DIAGONAL ELEMENTS OF G(T)*G                                           FITE-344
      J1=0                                                              FITE-345
      DO 56 I=1,N                                                       FITE-346
      S=0.D0                                                            FITE-347
      DO 55 L=1,M                                                       FITE-348
      J1=J1+1                                                           FITE-349
   55 S=S+WA(J1)*WA(J1)                                                 FITE-350
   56 W(N+I)=DSQRT(S)                                                   FITE-351
C STANDARD ERRORS AND ERROR CORRELATIONS                                FITE-352
      CALL FIT2(M,N,WA,W(2*N+1),NW,IR,.FALSE.)                          FITE-353
      IF (IR.NE.0) GO TO 63                                             FITE-354
      DO 57 I=1,N                                                       FITE-355
      JDI=2*N+(I*I+I)/2                                                 FITE-356
      W(I)=DSQRT(W(JDI))                                                FITE-357
   57 W(N+I)=W(I)*W(N+I)                                                FITE-358
      J1=2*N                                                            FITE-359
      DO 59 I=1,N                                                       FITE-360
      DO 58 J=1,I                                                       FITE-361
      J1=J1+1                                                           FITE-362
   58 W(J1)=W(J1)/(W(I)*W(J))                                           FITE-363
   59 CONTINUE                                                          FITE-364
C ERROR RENORMALISATION FACTOR                                          FITE-365
   60 S=0.D0                                                            FITE-366
      DO 61 I=1,M                                                       FITE-367
   61 S=S+F(I)                                                          FITE-368
      YY(3)=DSQRT(DABS(CHI2-S*S/M)/MAX0(M-N-1,1))                       FITE-369
      DO 62 I=1,N                                                       FITE-370
   62 W(I)=W(I)*YY(3)                                                   FITE-371
   63 IZ(4)=M-N-1                                                       FITE-372
      IF ((KE-5)*(KE-6).NE.0) GO TO 64                                  FITE-373
      IZ(3)=J-2                                                         FITE-374
      IZ(4)=MF-1                                                        FITE-375
   64 IF (KE.EQ.1) IZ(3)=IZ(3)+1                                        FITE-376
      RETURN                                                            FITE-377
      END                                                               FITE-378
