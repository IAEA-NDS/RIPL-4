C 07/03/07                                                      ECIS06  COR0-000
      SUBROUTINE COR0(EI,EF,SSI,SSF,ZR,ZI,Z,LM3,RM,RS,FS,FI,FF,EX,CC)   COR0-001
C    COMPUTATION OF THE INTEGRALS FROM RS TO INFINITY OF PRODUCTS OF    COR0-002
C REGULAR OR IRREGULAR COULOMB FUNCTIONS OF L=0 AND L=1 DIVIDED BY R**2 COR0-003
C INPUT:     EI,EF:   COULOMB PARAMETERS.                               COR0-004
C            SSI,SSF: COULOMB PHASE-SHIFTS FOR L=0.                     COR0-005
C            LM3:     DIMENSION OF WORKING AREA.                        COR0-006
C            RM:      ACTUAL MATCHING RADIUS FOR INTEGRALS.             COR0-007
C            RS:      NEEDED MATCHING RADIUS FOR INTEGRALS.             COR0-008
C            FS:      SQUARE ROOT OF RATIO OF WAVE NUMBERS.             COR0-009
C            FI,FF:   SQUARES OF COULOMB PARAMETERS.                    COR0-010
C            EX:      DIFFERENCE OF COULOMB PARAMETERS.                 COR0-011
C OUTPUT:    CC:      INTEGRALS FROM RS TO INFINITY.                    COR0-012
C            ZR,ZI:   REAL, IMAGINARY PART OF LOG(GAMMA(0.,EI-EF))      COR0-013
C WORKING AREA:                                                         COR0-014
C            Z:       FOR THE COMPUTATION OF F AND G.                   COR0-015
C***********************************************************************COR0-016
      IMPLICIT REAL*8 (A-H,O-Z)                                         COR0-017
      DIMENSION Z(4,*),X1(4),X2(4),X3(4),X4(4),X5(4),C(2,2,2),Y(4,4),IE(COR0-018
     12),ABSC(20),P(20),B(7),CC(4,4),YY(8,2)                            COR0-019
      EQUIVALENCE (Y,YY)                                                COR0-020
      DATA PI /3.1415926535897932D0/                                    COR0-021
      DATA B /.8333333333333333D-01,-.2777777777777778D-2,.7936507936507COR0-022
     19365D-3,-.5952380952380952D-3,.84175084175D-3,-.19175269D-2,.641D-COR0-023
     202/                                                               COR0-024
      DATA ABSC( 1),P( 1) / 0.88114514472040D-03, 0.22606385492666D-02/ COR0-025
      DATA ABSC( 2),P( 2) / 0.46368806502715D-02, 0.52491422655764D-02/ COR0-026
      DATA ABSC( 3),P( 3) / 0.11370025008113D-01, 0.82105291909539D-02/ COR0-027
      DATA ABSC( 4),P( 4) / 0.21041590393104D-01, 0.11122924597084D-01/ COR0-028
      DATA ABSC( 5),P( 5) / 0.33593595860662D-01, 0.13968503490012D-01/ COR0-029
      DATA ABSC( 6),P( 6) / 0.48950596515563D-01, 0.16730097641274D-01/ COR0-030
      DATA ABSC( 7),P( 7) / 0.67020248393870D-01, 0.19391083987236D-01/ COR0-031
      DATA ABSC( 8),P( 8) / 0.87693884583344D-01, 0.21935454092836D-01/ COR0-032
      DATA ABSC( 9),P( 9) / 0.11084717428674D+00, 0.24347903817536D-01/ COR0-033
      DATA ABSC(10),P(10) / 0.13634087240504D+00, 0.26613923491968D-01/ COR0-034
      DATA ABSC(11),P(11) / 0.16402165769291D+00, 0.28719884549696D-01/ COR0-035
      DATA ABSC(12),P(12) / 0.19372305516601D+00, 0.30653121246465D-01/ COR0-036
      DATA ABSC(13),P(13) / 0.22526643745244D+00, 0.32402006728300D-01/ COR0-037
      DATA ABSC(14),P(14) / 0.25846209915691D+00, 0.33956022907617D-01/ COR0-038
      DATA ABSC(15),P(15) / 0.29311039781420D+00, 0.35305823695643D-01/ COR0-039
      DATA ABSC(16),P(16) / 0.32900295458712D+00, 0.36443291197902D-01/ COR0-040
      DATA ABSC(17),P(17) / 0.36592390749637D+00, 0.37361584528984D-01/ COR0-041
      DATA ABSC(18),P(18) / 0.40365120964931D+00, 0.38055180950313D-01/ COR0-042
      DATA ABSC(19),P(19) / 0.44195796466237D+00, 0.38519909082124D-01/ COR0-043
      DATA ABSC(20),P(20) / 0.48061379124697D+00, 0.38752973989212D-01/ COR0-044
C COMPUTATION OF THE INTEGRALS FROM RM TO INFINITY FOR L=0,1 BY         COR0-045
C RAWITSCHER'S METHOD. (COMP. PHYSICS COMM., VOL.11,183,1976)           COR0-046
C EXPANSIONS OF COULOMB WAVE FUNCTIONS ARE Z(N,2*M-1)+I*Z(N,2*M)=G+I*F  COR0-047
C   N=1: L=0,E=EI     N=2: L=0,E=EF     N=3: L=1,E=EI     N=4: L=1,E=EF COR0-048
      A1=RM/FS                                                          COR0-049
      B1=RM*FS                                                          COR0-050
      X1(1)=A1-EI*DLOG(2.D0*A1)+SSI                                     COR0-051
      X1(3)=B1-EF*DLOG(2.D0*B1)+SSF                                     COR0-052
      X1(2)=X1(1)+DATAN(EI)-.5D0*PI                                     COR0-053
      X1(4)=X1(3)+DATAN(EF)-.5D0*PI                                     COR0-054
C FIRST TERM OF ASYMPTOTIC EXPANSION OF H(+/-) FOR L=0,1.               COR0-055
      DO 1 I=1,4                                                        COR0-056
      X5(I)=0.D0                                                        COR0-057
      Z(I,1)=DCOS(X1(I))                                                COR0-058
    1 Z(I,2)=DSIN(X1(I))                                                COR0-059
      X1(1)=-A1-B1                                                      COR0-060
      X2(1)=EI+EF                                                       COR0-061
      X1(2)=B1-A1                                                       COR0-062
      X2(2)=-EX                                                         COR0-063
      IF (DABS(X2(2)).LT.1.D-8) X2(2)=0.D0                              COR0-064
      IF (X2(2).EQ.0.D0) GO TO 4                                        COR0-065
C COMPUTATION OF LOG(GAMMA(0.,EI-EF)).                                  COR0-066
      B2=-DATAN(EX/11.D0)                                               COR0-067
      B3=121.D0+EX**2                                                   COR0-068
      A2=.5D0*DLOG(B3)                                                  COR0-069
      ZR=.91893853320467274D0+10.5D0*A2+EX*B2-11.D0-DLOG(DABS(EX))      COR0-070
      ZI=10.5D0*B2+EX-EX*A2+0.5D0*PI                                    COR0-071
      IF (EX.LT.0.D0) ZI=ZI-PI                                          COR0-072
      A2=11.D0/B3                                                       COR0-073
      B2=EX/B3                                                          COR0-074
      A3=A2**2-B2**2                                                    COR0-075
      B3=2.D0*A2*B2                                                     COR0-076
      DO 2 I=1,7                                                        COR0-077
      ZR=ZR+B(I)*A2                                                     COR0-078
      ZI=ZI+B(I)*B2                                                     COR0-079
      C1=A2*A3-B2*B3                                                    COR0-080
      B2=A2*B3+B2*A3                                                    COR0-081
    2 A2=C1                                                             COR0-082
      DO 3 I=1,10                                                       COR0-083
      A2=DFLOAT(I)                                                      COR0-084
      ZR=ZR-.5D0*DLOG(A2**2+EX**2)                                      COR0-085
    3 ZI=ZI+DATAN(EX/A2)                                                COR0-086
    4 IF (DABS(X1(2)).GT.2.D0) GO TO 9                                  COR0-087
C COMPUTATION OF THE FIRST CONFLUENT HYPERGEOMETRIC FUNCTION.           COR0-088
      IF (X1(2).EQ.0.D0) GO TO 8                                        COR0-089
      IF (DABS(X2(2)).LT.1.D-3) GO TO 5                                 COR0-090
      A2=DEXP(-DSIGN(0.5D0*PI,X1(2))*X2(2)+ZR)                          COR0-091
      B2=DSIGN(ZI,X2(2))+X2(2)*DLOG(DABS(X1(2)))                        COR0-092
      X5(2)=A2*DCOS(B2)                                                 COR0-093
      X5(4)=A2*DSIN(B2)-1.D0/X2(2)                                      COR0-094
      GO TO 6                                                           COR0-095
C EXPANSION FOR SMALL EX.                                               COR0-096
    5 A2=DLOG(DABS(X1(2)))                                              COR0-097
      B2=PI*.5D0                                                        COR0-098
      IF (X1(2).LT.0.D0) B2=-B2                                         COR0-099
      A3=A2*(1.D0+EX*(B2-EX*((A2**2-3.D0*B2**2)/6.D0+EX*B2*(A2**2-B2**2)COR0-100
     1/6.D0)))                                                          COR0-101
      B3=B2-EX*((A2**2-B2**2)/2.D0+EX*(B2*(3.D0*A2**2-B2**2)/6.D0-EX*(A2COR0-102
     1**4-6.D0*A2**2*B2**2+B2**4)/24.D0))                               COR0-103
      A2=.5772156649015329D0-.400685634386331D0*EX**2                   COR0-104
      B2=-EX*(.822467033424113D0-0.270580808427784D0*EX**2)             COR0-105
      A4=A2*(1.D0+EX*(B2-EX*((A2**2-3.D0*B2**2)/6.D0+EX*B2*(A2**2-B2**2)COR0-106
     1/6.D0)))                                                          COR0-107
      B4=B2-EX*((A2**2-B2**2)/2.D0+EX*(B2*(3.D0*A2**2-B2**2)/6.D0-EX*(A2COR0-108
     1**4-6.D0*A2**2*B2**2+B2**4)/24.D0))                               COR0-109
      X5(2)=-A3-A4-EX*(A3*B4+B3*A4)                                     COR0-110
      X5(4)=-B3-B4+EX*(A3*A4-B3*B4)                                     COR0-111
    6 A2=X5(2)                                                          COR0-112
      B2=X5(4)                                                          COR0-113
      A3=1.D0                                                           COR0-114
      B3=0.D0                                                           COR0-115
      DO 7 J=1,1000                                                     COR0-116
      B4=J                                                              COR0-117
      C1=-X1(2)*(B4*B3-EX*A3)/(B4**2+EX**2)                             COR0-118
      B3=X1(2)*(B4*A3+EX*B3)/(B4**2+EX**2)                              COR0-119
      A3=C1                                                             COR0-120
      C2=(-X1(2)*B2+A3)/B4                                              COR0-121
      B2=(X1(2)*A2+B3)/B4                                               COR0-122
      A2=C2                                                             COR0-123
      X5(2)=X5(2)+A2                                                    COR0-124
      X5(4)=X5(4)+B2                                                    COR0-125
      IF (DABS(A2)+DABS(B2).LT.1.D-12*(DABS(X5(2))+DABS(X5(4)))) GO TO 9COR0-126
    7 CONTINUE                                                          COR0-127
      GO TO 9                                                           COR0-128
    8 X5(4)=-.5D0*PI                                                    COR0-129
      X5(2)=0.D0                                                        COR0-130
    9 DO 11 K=1,4                                                       COR0-131
      DO 10 J=1,4                                                       COR0-132
   10 Y(J,K)=0.D0                                                       COR0-133
   11 CONTINUE                                                          COR0-134
      N=IDINT(DMIN1(A1+DSQRT(A1**2-FI),B1+DSQRT(B1**2-FF))+2.D0)        COR0-135
      A3=0.D0                                                           COR0-136
      B3=0.D0                                                           COR0-137
C LOOP OF THE ASYMPTOTIC EXPANSION.                                     COR0-138
      DO 27 I=1,N                                                       COR0-139
      IF (4*I.GT.LM3) CALL MEMO('COR0',LM3,4*I)                         COR0-140
      A4=DFLOAT(I)-2                                                    COR0-141
      IF (I.EQ.1) GO TO 19                                              COR0-142
      IF (I.EQ.2) GO TO 13                                              COR0-143
      B4=1.D0-.5D0/A4                                                   COR0-144
      X3(1)=EI*B4/A1                                                    COR0-145
      X3(2)=X3(1)                                                       COR0-146
      X3(3)=EF*B4/B1                                                    COR0-147
      X3(4)=X3(3)                                                       COR0-148
      X4(1)=(FI-A4*(A4-1.D0))/(2.D0*A4*A1)                              COR0-149
      X4(2)=(FI-(A4+1.D0)*(A4-2.D0))/(2.D0*A4*A1)                       COR0-150
      X4(3)=(FF-A4*(A4-1.D0))/(2.D0*A4*B1)                              COR0-151
      X4(4)=(FF-(A4+1.D0)*(A4-2.D0))/(2.D0*A4*B1)                       COR0-152
C NEW TERM OF ASYMPTOTIC EXPANSION.                                     COR0-153
      DO 12 J=1,4                                                       COR0-154
      Z(J,2*I-3)=Z(J,2*I-5)*X3(J)-Z(J,2*I-4)*X4(J)                      COR0-155
   12 Z(J,2*I-2)=Z(J,2*I-5)*X4(J)+Z(J,2*I-4)*X3(J)                      COR0-156
   13 DO 16 L=1,2                                                       COR0-157
      DO 15 K=1,2                                                       COR0-158
      DO 14 J=1,2                                                       COR0-159
   14 C(J,K,L)=0.D0                                                     COR0-160
   15 CONTINUE                                                          COR0-161
   16 CONTINUE                                                          COR0-162
C PRODUCT OF THE ASYMPTOTIC EXPANSIONS IN C(1,N,M)+I*C(2,N,M)           COR0-163
C  N=1: H0(+)(EI)*H0(+)(EF)     N=2: H1(+)(EI)*H1(+)(EF)                COR0-164
C  M=1: HL(+)(EI)*HL(+)(EF)     M=2: HL(+)(EI)*HL(+)(EF).               COR0-165
      I1=I-1                                                            COR0-166
      DO 18 J=1,I1                                                      COR0-167
      M=I-J                                                             COR0-168
      DO 17 L=1,2                                                       COR0-169
      C(1,L,1)=C(1,L,1)+Z(L,2*J-1)*Z(L+2,2*M-1)-Z(L,2*J)*Z(L+2,2*M)     COR0-170
      C(2,L,1)=C(2,L,1)+Z(L,2*J-1)*Z(L+2,2*M)+Z(L,2*J)*Z(L+2,2*M-1)     COR0-171
      C(1,L,2)=C(1,L,2)+Z(L,2*J-1)*Z(L+2,2*M-1)+Z(L,2*J)*Z(L+2,2*M)     COR0-172
   17 C(2,L,2)=C(2,L,2)-Z(L,2*J-1)*Z(L+2,2*M)+Z(L,2*J)*Z(L+2,2*M-1)     COR0-173
   18 CONTINUE                                                          COR0-174
   19 B4=A4+2.D0                                                        COR0-175
      A5=0.D0                                                           COR0-176
      B5=0.D0                                                           COR0-177
C INTEGRATION FROM RS TO INFINITY.                                      COR0-178
      DO 26 M=1,2                                                       COR0-179
C TRANSFER OF PREVIOUS INTEGRALS.                                       COR0-180
      X4(M)=X5(M)*RM                                                    COR0-181
      X4(M+2)=X5(M+2)*RM                                                COR0-182
      IF (DABS(X1(M)).LE.2.D0.AND.M.EQ.2) GO TO 23                      COR0-183
C PADE METHOD FOR OMEGA FUNCTION.                                       COR0-184
C INSTEAD OF FORMULA (44), WE COMPUTE I*K*R*OMEGA IN TERMS OF I*K*R.    COR0-185
      A9=B4**2+(X2(M)+X1(M))**2                                         COR0-186
      A6=B4/A9                                                          COR0-187
      B6=-(X2(M)+X1(M))/A9                                              COR0-188
      A7=-B6*X1(M)                                                      COR0-189
      B7=A6*X1(M)                                                       COR0-190
      A8=A6                                                             COR0-191
      B8=B6+1.D0/X1(M)                                                  COR0-192
      B9=1.D0                                                           COR0-193
C STEED'S ALGORITHM.                                                    COR0-194
      DO 20 J=1,1000                                                    COR0-195
      C2=B7+X1(M)/B9                                                    COR0-196
      A9=(A7**2+C2**2)*B9                                               COR0-197
      C1=A7/A9                                                          COR0-198
      C2=-C2/A9                                                         COR0-199
      C3=-(C2*X1(M)+1.D0)*A8-C1*X1(M)*B8                                COR0-200
      C4=C1*X1(M)*A8-B8*(C2*X1(M)+1.D0)                                 COR0-201
      A7=1.D0+C1*(B4+B9)-C2*X2(M)                                       COR0-202
      B7=C1*X2(M)+(B4+B9)*C2                                            COR0-203
      A9=A7**2+B7**2                                                    COR0-204
      A7=A7/A9                                                          COR0-205
      B7=-B7/A9                                                         COR0-206
      A8=(A7-1.D0)*C3-B7*C4                                             COR0-207
      B8=(A7-1.D0)*C4+B7*C3                                             COR0-208
      B9=B9+1.D0                                                        COR0-209
      IF (DABS(C3)+DABS(C4)+DABS(A8)+DABS(B8).LT.1.D-12*(DABS(B6)+DABS(ACOR0-210
     16))) GO TO 21                                                     COR0-211
      A6=A6+A8+C3                                                       COR0-212
   20 B6=B6+B8+C4                                                       COR0-213
   21 X5(M)=0.D0                                                        COR0-214
      X5(M+2)=0.D0                                                      COR0-215
      J=J+1                                                             COR0-216
C DIRECT COMPUTATION OF PADE APPROXIMANT FOR MORE PRECISION.            COR0-217
      DO 22 K=1,J                                                       COR0-218
      C1=1.D0+B9*X5(M)                                                  COR0-219
      C2=B9*X5(M+2)                                                     COR0-220
      A9=C1**2+C2**2                                                    COR0-221
      B9=B9-1.D0                                                        COR0-222
      A7=((B4+B9)*C1+X2(M)*C2)/A9                                       COR0-223
      B7=X1(M)+(C1*X2(M)-C2*(B4+B9))/A9                                 COR0-224
      A9=A7**2+B7**2                                                    COR0-225
      X5(M)=A7/A9                                                       COR0-226
   22 X5(M+2)=-B7/A9                                                    COR0-227
      GO TO 24                                                          COR0-228
C TAYLOR EXPANSION OF CONFLUENT HYPERGEOMETRIC FUNCTION.                COR0-229
   23 IF (I.EQ.1) GO TO 26                                              COR0-230
      A7=1.D0+X1(M)*X5(M+2)                                             COR0-231
      B7=-X1(M)*X5(M)                                                   COR0-232
      A9=((B4-1.D0)**2+X2(M)**2)                                        COR0-233
      X5(M)=(A7*(B4-1.D0)+X2(M)*B7)/A9                                  COR0-234
      X5(M+2)=(B7*(B4-1.D0)-X2(M)*A7)/A9                                COR0-235
   24 IF (I.EQ.1) GO TO 26                                              COR0-236
C STORAGE OF INTEGRALS TO A FACTOR R IN Y(J,M)                          COR0-237
C  J=1 AND J=2  REAL AND IMAGINARY PARTS OF INTEGRAL OF H0(EI)*H0(EF)   COR0-238
C  J=3 AND J=4  REAL AND IMAGINARY PARTS OF INTEGRAL OF H1(EI)*H1(EF)   COR0-239
C  M=1 AND M=2 FOR H*H   M=3 AND M=4 FOR H*H*                           COR0-240
C  M=1 AND M=3 FOR 1/R   M=2 AND M=4 FOR 1/R**2                         COR0-241
      DO 25 J=1,2                                                       COR0-242
      A6=C(1,J,M)*X4(M)-C(2,J,M)*X4(M+2)                                COR0-243
      B6=C(1,J,M)*X4(M+2)+C(2,J,M)*X4(M)                                COR0-244
      A5=DMAX1(A5,A6**2+B6**2)                                          COR0-245
      Y(2*J-1,2*M-1)=Y(2*J-1,2*M-1)+A6                                  COR0-246
      Y(2*J,2*M-1)=Y(2*J,2*M-1)+B6                                      COR0-247
      A3=DMAX1(A3,Y(2*J-1,2*M-1)**2+Y(2*J,2*M-1)**2)                    COR0-248
      A6=C(1,J,M)*X5(M)-C(2,J,M)*X5(M+2)                                COR0-249
      B6=C(1,J,M)*X5(M+2)+C(2,J,M)*X5(M)                                COR0-250
      B5=DMAX1(B5,A6**2+B6**2)                                          COR0-251
      Y(2*J-1,2*M)=Y(2*J-1,2*M)+A6                                      COR0-252
      Y(2*J,2*M)=Y(2*J,2*M)+B6                                          COR0-253
   25 B3=DMAX1(B3,Y(2*J-1,2*M)**2+Y(2*J,2*M)**2)                        COR0-254
   26 CONTINUE                                                          COR0-255
      IF (I.NE.1.AND.A5.LT.A3*1.D-30.AND.B5.LT.B3*1.D-30) GO TO 28      COR0-256
   27 CONTINUE                                                          COR0-257
   28 A1=2.D0*RM                                                        COR0-258
C TRANSFORMATION FROM H(+/-) TO F AND G.                                COR0-259
      DO 29 I=1,4                                                       COR0-260
      CC(1,I)=(YY(2*I-1,2)-YY(2*I-1,1))/A1                              COR0-261
      CC(4,I)=(YY(2*I-1,2)+YY(2*I-1,1))/A1                              COR0-262
      CC(3,I)=(YY(2*I,1)+YY(2*I,2))/A1                                  COR0-263
   29 CC(2,I)=(YY(2*I,1)-YY(2*I,2))/A1                                  COR0-264
      IF (RM.EQ.RS) GO TO 35                                            COR0-265
C IF RM IS NOT RS, 40 POINTS GAUSS INTEGRATION FROM RS TO RM            COR0-266
C EACH GAUSS INTEGRATION IS FOR VARIATION OF R LESS THAN 20.            COR0-267
      IL=1+IDINT(DABS(RM-RS)*.05D0)                                     COR0-268
      A1=DFLOAT(IL)                                                     COR0-269
      A1=(RM-RS)/A1                                                     COR0-270
      A2=RM                                                             COR0-271
      DO 34 IT=1,IL                                                     COR0-272
      A3=A2-A1                                                          COR0-273
      DO 33 II=1,40                                                     COR0-274
      I=MIN0(II,41-II)                                                  COR0-275
      A4=ABSC(I)                                                        COR0-276
      IF (I.NE.II) A4=1.D0-A4                                           COR0-277
      A5=A3+A1*A4                                                       COR0-278
      CALL FCOU(1,EI,A5/FS,Y,Y(1,3),Y(3,1),Y(3,3),IE,X1)                COR0-279
      CALL FCOU(1,EF,A5*FS,Y(1,2),Y(1,4),Y(3,2),Y(3,4),IE,X1)           COR0-280
      DO 32 M=1,2                                                       COR0-281
      DO 31 J=1,2                                                       COR0-282
      DO 30 L=1,2                                                       COR0-283
      A6=P(I)*Y(2*J+M-2,1)*Y(2*L+M-2,2)/A5*A1                           COR0-284
      CC(2*L+J-2,M)=CC(2*L+J-2,M)+A6                                    COR0-285
   30 CC(2*L+J-2,M+2)=CC(2*L+J-2,M+2)+A6/A5                             COR0-286
   31 CONTINUE                                                          COR0-287
   32 CONTINUE                                                          COR0-288
   33 CONTINUE                                                          COR0-289
   34 A2=A3                                                             COR0-290
   35 RETURN                                                            COR0-291
      END                                                               COR0-292
