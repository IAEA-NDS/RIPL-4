C 27/03/07                                                      ECIS06  LCSP-000
      SUBROUTINE LCSP(F,FCN,JMAX,KMAX,IPI,NCOLL,NCOLS,MF,WV,JMIN,IPJ,IPKLCSP-001
     1,A,AX,JD,LO)                                                      LCSP-002
C COEFFICIENTS OF LEGENDRE POLYNOMIALS FOR CROSS-SECTIONS               LCSP-003
C SCATTERING COEFFICIENTS IN THE HELICITY REPRESENTATION                LCSP-004
C INPUT:     F:       S-MATRIX.                                         LCSP-005
C            FCN:     TRANSMISSION COEFFICIENTS.                        LCSP-006
C            JMAX:    MAXIMUM NUMBER OF CHANNEL SPINS, DIMENSION FOR F. LCSP-007
C            KMAX:    MAXIMUM NUMBER L FOR COMPOUND NUCLEUS, DIMENSION  LCSP-008
C                     FOR FCN.                                          LCSP-009
C            IPI(I,*):MULTIPLICITIES OF TARGET/PARTICLE FOR I=2/3,      LCSP-010
C                     FIRST/LAST CHANNEL NUMBER FOR I=6/7 (SEE DEPH).   LCSP-011
C            NCOLL:   NUMBER OF COUPLED LEVELS.                         LCSP-012
C            NCOLS:   NUMBER OF LEVELS WITH ANGULAR DISTRIBUTION.       LCSP-013
C            MF:      HELICITY NUMBERS (SEE DEPH).                      LCSP-014
C            WV:      WAVE NUMBER AND COULOMB PARAMETER (SEE COLF).     LCSP-015
C            JMIN:    TWICE MINIMUM CHANNEL SPIN.                       LCSP-016
C            IPJ:     NUMBER OF CHANNEL SPIN.                           LCSP-017
C            IPK:     NUMBER OF CHANNEL SPIN FOR COMPOUND NUCLEUS.      LCSP-018
C            JD:      DIMENSION OF WORKING FIELD AX.                    LCSP-019
C            LO(I):   LOGICAL CONTROLS:                                 LCSP-020
C               LO(18) =.TRUE. PROJECTILE-TARGET ANTISYMMETRISATION.    LCSP-021
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             LCSP-022
C WORKING AREAS:                                                        LCSP-023
C            A(7,*):  AMPLITUDE MULTIPLIED BY A LEGENDRE POLYNOMIAL     LCSP-024
C                     ALTERNATIVELY IN A(1-2,*) AND A(3-4,*),           LCSP-025
C                     COEFFICIENTS OF RECURRENCE IN A(5-7,*).           LCSP-026
C            AX:      1 TO JML  LEGENDRE COEFFICIENTS OF CROSS-SECTION: LCSP-027
C                     JML+1 TO JMT    FOR COMPOUND NUCLEUS,             LCSP-028
C                     JMT+1 TO JMX    PRODUCT SPIN ROTATION MATRICES.   LCSP-029
C                                                                       LCSP-030
C FOR THE COMMON  /ANGUL/ SEE LECT.                                     LCSP-031
C FOR THE COMMON  /DCONS/ SEE CALC.                                     LCSP-032
C                                                                       LCSP-033
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     LCSP-034
C  XZ:        CONVERSION FACTOR TO MILLIBARNS.                          LCSP-035
C   USED:     XZ.                                                       LCSP-036
C                                                                       LCSP-037
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ANGUL/:                     LCSP-038
C  NL(1):     POWER OF (1-COS(THETA)) FOR THE EXPANSION IN LEGENDRE     LCSP-039
C             POLYNOMIALS OF THE INTERFERENCE BETWEEN COULOMB AND       LCSP-040
C             NUCLEAR ELASTIC SCATTERING. POWER OF (1-COS(THETA)**2)    LCSP-041
C             IF LO(18) IS .TRUE..                                      LCSP-042
C  NL(2):     NUMBER OF LEGENDRE POLYNOMIALS FOR ELASTIC SCATTERING OF  LCSP-043
C             CHARGED PARTICLES.                                        LCSP-044
C  NL(3):     NUMBER OF LEGENDRE POLYNOMIALS FOR ELASTIC SCATTERING OF  LCSP-045
C             UNCHARGED PARTICLES, INELASTIC SCATTERING AND COMPOUND    LCSP-046
C             NUCLEUS.                                                  LCSP-047
C   USED:     NL.                                                       LCSP-048
C                                                                       LCSP-049
C***********************************************************************LCSP-050
      IMPLICIT REAL*8 (A-H,O-Z)                                         LCSP-051
      LOGICAL LO(150),LS                                                LCSP-052
      DIMENSION F(2,JMAX,*),FCN(KMAX,*),IPI(11,*),MF(10,*),WV(22,*),A(7,LCSP-053
     1*),AX(*),AT(2),NDL(2)                                             LCSP-054
      COMMON /ANGUL/ THETA1,THETA2,DTHETA,DTHE,NCJ,NL(3),JMM(2)         LCSP-055
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            LCSP-056
      COMMON /INOUT/ MR,MW,MS                                           LCSP-057
C DETERMINATION OF THE NUMBER OF LEGENDRE POLYNOMIALS.                  LCSP-058
      NDL(1)=NL(2)                                                      LCSP-059
      NDL(2)=NL(3)                                                      LCSP-060
      IF (NDL(1).EQ.0) NDL(1)=(6*IPJ-3+3*JMIN)/2                        LCSP-061
      IF (NDL(2).EQ.0) NDL(2)=2*IPJ-1+JMIN                              LCSP-062
      WRITE (65,1000) WV(1,1),WV(13,1),WV(2,1),IPI(4,1),NCOLS           LCSP-063
      WRITE (MW,1001)                                                   LCSP-064
C LOOP ON LEVELS.                                                       LCSP-065
      DO 36 JK=1,NCOLS                                                  LCSP-066
      IF (WV(3,JK).LT.0.D0) GO TO 35                                    LCSP-067
      IF (JK.GT.NCOLL) GO TO 32                                         LCSP-068
      JK1=IPI(6,JK)                                                     LCSP-069
      JK2=IPI(7,JK)                                                     LCSP-070
      JML=NDL(2)                                                        LCSP-071
      LS=JK.EQ.1.AND.WV(5,1).NE.0.D0                                    LCSP-072
      IF (LS) JML=NDL(1)                                                LCSP-073
      JMT=JML                                                           LCSP-074
      IF (JMT.GT.JD) CALL MEMO('LCSP',ID,JMT)                           LCSP-075
      DO 1 J=1,JMT                                                      LCSP-076
    1 AX(J)=0.D0                                                        LCSP-077
C LOOP ON AMPLITUDES.                                                   LCSP-078
      DO 29 I=JK1,JK2                                                   LCSP-079
      IF (MF(6,I).EQ.99999) GO TO 3                                     LCSP-080
      M1=MF(5,I)                                                        LCSP-081
      M2=MF(6,I)                                                        LCSP-082
      M3=((IABS(M1+M2)+IABS(M1-M2))/2-JMIN)/2+1                         LCSP-083
      IF (JK.GT.NCOLL) GO TO 3                                          LCSP-084
      A1=0.25D0*DFLOAT(M1*M2)                                           LCSP-085
      MJ=2*M3+JMIN-2                                                    LCSP-086
      A2=0.5D0*DFLOAT(MJ)                                               LCSP-087
      X1=DFLOAT((MJ+M2)/2+1-M3)                                         LCSP-088
      X2=DFLOAT((MJ-M2)/2+1-M3)                                         LCSP-089
      X3=DFLOAT((MJ+M1)/2+1-M3)                                         LCSP-090
      X4=DFLOAT((MJ-M1)/2+1-M3)                                         LCSP-091
      A(7,M3)=0.D0                                                      LCSP-092
C COEFFICIENTS OF THE RECURRENCE.                                       LCSP-093
      JM1=JML+IPJ                                                       LCSP-094
      JM=JM1+1                                                          LCSP-095
      IF (M3.GT.JM1) GO TO 3                                            LCSP-096
      DO 2 J=M3,JM1                                                     LCSP-097
      XJ=DFLOAT(J)                                                      LCSP-098
      A2=A2+1.D0                                                        LCSP-099
      A(6,J)=0.D0                                                       LCSP-100
      IF (A1.NE.0.D0) A(6,J)=A1/(A2*A2-A2)                              LCSP-101
      A(5,J)=DSQRT((X1+XJ)*(X2+XJ)*(X3+XJ)*(X4+XJ))/(A2*(2.D0*A2+1.D0)) LCSP-102
    2 A(7,J+1)=A(5,J)*(A2+A2+1.D0)/(A2+A2-1.D0)                         LCSP-103
    3 IF (M3.GT.JM) GO TO 29                                            LCSP-104
      DO 4 J=M3,JM                                                      LCSP-105
      A(1,J)=0.D0                                                       LCSP-106
      A(2,J)=0.D0                                                       LCSP-107
      A(3,J)=0.D0                                                       LCSP-108
      IF (J.GT.IPJ) GO TO 4                                             LCSP-109
      A(1,J)=F(1,J,I)                                                   LCSP-110
      A(2,J)=F(2,J,I)                                                   LCSP-111
    4 A(4,J)=0.D0                                                       LCSP-112
      IF (.NOT.LS) GO TO 24                                             LCSP-113
C COULOMB AMPLITUDE.                                                    LCSP-114
      NS1=IPI(2,1)-1                                                    LCSP-115
      NS2=IPI(3,1)-1                                                    LCSP-116
      NSM=IABS(NS1-NS2)-JMIN                                            LCSP-117
      NSP=NS1+NS2-JMIN                                                  LCSP-118
      XS3=DFLOAT(NS1+NS2+2)                                             LCSP-119
      XS4=DFLOAT(NS1-NS2)                                               LCSP-120
      NST=NSP/2+1                                                       LCSP-121
      JMX=JMT+NST                                                       LCSP-122
      JMY=JMX+2                                                         LCSP-123
      IF (LO(18)) JMY=JMY+2*NL(1)                                       LCSP-124
      IF (JMY.GT.JD) CALL MEMO('LCSP',ID,JMY)                           LCSP-125
      DO 5 J=1,NST                                                      LCSP-126
    5 AX(JMT+J)=1.D0                                                    LCSP-127
C REDUCED ROTATION MATRIX ELEMENT.                                      LCSP-128
      DO 9 I1=1,2                                                       LCSP-129
      MS1=2*MF(2*I1-1,I)-NS1-2                                          LCSP-130
      MS2=NS2-2*MF(2*I1,I)+2                                            LCSP-131
      MSP=IABS(MS1+MS2)-JMIN                                            LCSP-132
      NSQ=MAX0(NSM,MSP)/2+1                                             LCSP-133
      XS=MS1+MS2                                                        LCSP-134
      AT(I1)=1.D0                                                       LCSP-135
      IF (NSQ.EQ.1) GO TO 7                                             LCSP-136
      DO 6 J=2,NSQ                                                      LCSP-137
    6 AX(JMT+J-1)=0.D0                                                  LCSP-138
    7 IF (NSQ.EQ.NST) GO TO 9                                           LCSP-139
      NS=NSQ+1                                                          LCSP-140
      A2=0.D0                                                           LCSP-141
      A1=1.D0                                                           LCSP-142
      B2=0.D0                                                           LCSP-143
      DO 8 K=NS,NST                                                     LCSP-144
      A3=A2                                                             LCSP-145
      A2=A1                                                             LCSP-146
      B1=B2                                                             LCSP-147
      XJ=DFLOAT(2*K+JMIN-2)                                             LCSP-148
      B2=DSQRT((XJ**2-XS**2)*(XS3**2-XJ**2)*(XJ**2-XS4**2)/(XJ**2-1.D0))LCSP-149
     1/(2.D0*XJ)                                                        LCSP-150
      B3=DFLOAT(MS1-MS2)                                                LCSP-151
      IF (XJ.NE.2.D0) B3=B3-XS*XS3*XS4/(XJ*(XJ-2.D0))                   LCSP-152
      A1=(A2*B3-A3*B1)/B2                                               LCSP-153
      AX(JMT+K)=AX(JMT+K)*A1                                            LCSP-154
    8 AT(I1)=AT(I1)+A1**2                                               LCSP-155
    9 CONTINUE                                                          LCSP-156
      IJX=NL(1)                                                         LCSP-157
      A3=1.D0                                                           LCSP-158
      JP=IPJ                                                            LCSP-159
C MULTIPLICATION BY (1-COS), (1-COS**2) OR THEIR SQUARES.               LCSP-160
   10 DO 12 IJ=1,IJX                                                    LCSP-161
      A1=0.D0                                                           LCSP-162
      B1=0.D0                                                           LCSP-163
      DO 11 J=M3,JP                                                     LCSP-164
      A2=A1                                                             LCSP-165
      B2=B1                                                             LCSP-166
      A1=A(1,J)                                                         LCSP-167
      B1=A(2,J)                                                         LCSP-168
      A(1,J)=A1-A3*(A(7,J)*A2+A(6,J)*A1+A(5,J)*A(1,J+1))                LCSP-169
   11 A(2,J)=B1-A3*(A(7,J)*B2+A(6,J)*B1+A(5,J)*A(2,J+1))                LCSP-170
      JP=JP+1                                                           LCSP-171
      A(1,JP)=-A3*A1*A(6,JP)                                            LCSP-172
   12 A(2,JP)=-A3*B1*A(6,JP)                                            LCSP-173
      A3=-A3                                                            LCSP-174
      IF (LO(18).AND.A3.LT.0.D0) GO TO 10                               LCSP-175
C INITIALISATION OF COULOMB PHASE SHIFTS.                               LCSP-176
      C1=DSQRT(AT(1)*AT(2))                                             LCSP-177
      IF (AX(JMX).LT.0.D0) C1=-C1                                       LCSP-178
      C2=DFLOAT(IJX)                                                    LCSP-179
      C3=2.D0**(IJX-1)                                                  LCSP-180
      DO 23 IK=1,2                                                      LCSP-181
      NX=0                                                              LCSP-182
      NST1=NST                                                          LCSP-183
      DO 13 J=1,NST                                                     LCSP-184
      A1=2.D0*AX(JMT+J)/C1*DFLOAT(1-MOD(J+IK,2))                        LCSP-185
      IF (A1.NE.0.D0) NX=NX+1                                           LCSP-186
   13 A(3,J)=A1                                                         LCSP-187
      IF (NX.EQ.0) GO TO 23                                             LCSP-188
      AX(JMX+1)=-C2*C3*WV(5,1)/(C2*C2+WV(5,1)**2)                       LCSP-189
      AX(JMX+2)=-C3*WV(5,1)**2/(C2*C2+WV(5,1)**2)                       LCSP-190
      IF (.NOT.LO(18)) GO TO 15                                         LCSP-191
      AX(JMX+1)=4.D0*C3*AX(JMX+1)                                       LCSP-192
      AX(JMX+2)=4.D0*C3*AX(JMX+2)                                       LCSP-193
      DO 14 J=1,IJX                                                     LCSP-194
      A2=DFLOAT(J)*(DFLOAT(J+IJX)**2+WV(5,1)**2)                        LCSP-195
      A1=-DFLOAT(IJX-J+1)*(DFLOAT(IJX+J)*DFLOAT(IJX+J-1)+WV(5,1)**2)/A2 LCSP-196
      B1=DFLOAT(IJX-J+1)*WV(5,1)/A2                                     LCSP-197
      AX(JMX+2*J+1)=AX(JMX+2*J-1)*A1-AX(JMX+2*J)*B1                     LCSP-198
   14 AX(JMX+2*J+2)=AX(JMX+2*J-1)*B1+AX(JMX+2*J)*A1                     LCSP-199
   15 I1=3                                                              LCSP-200
      DO 21 L=1,JM                                                      LCSP-201
      B1=2*L-1                                                          LCSP-202
C COULOMB PHASE SHIFT.                                                  LCSP-203
      A1=AX(JMX+1)                                                      LCSP-204
      A2=AX(JMX+2)                                                      LCSP-205
      IF (.NOT.LO(18)) GO TO 17                                         LCSP-206
      DO 16 J=1,IJX                                                     LCSP-207
      A1=A1+AX(JMX+2*J+1)                                               LCSP-208
   16 A2=A2+AX(JMX+2*J+2)                                               LCSP-209
   17 DO 18 J=1,NST1                                                    LCSP-210
      IF (LO(18).AND.MOD(L+IK,2).EQ.1) GO TO 18                         LCSP-211
      A(1,J)=A(1,J)+A1*A(I1,J)*B1                                       LCSP-212
      A(2,J)=A(2,J)+A2*A(I1,J)*B1                                       LCSP-213
   18 CONTINUE                                                          LCSP-214
      I2=7-I1                                                           LCSP-215
      IF (L.EQ.JM) GO TO 21                                             LCSP-216
C PRODUCT OF ROTATION MATRIX ELEMENTS BY NEXT PL.                       LCSP-217
      B2=-(B1-1.D0)/(B1+1.D0)                                           LCSP-218
      B3=2.D0*B1/(B1+1.D0)                                              LCSP-219
      A1=0.D0                                                           LCSP-220
      DO 19 J=M3,NST1                                                   LCSP-221
      A2=A1                                                             LCSP-222
      A1=A(I1,J)                                                        LCSP-223
      A(I2,J)=B2*A(I2,J)+(A(7,J)*A2+A(6,J)*A1+A(5,J)*A(I1,J+1))*B3      LCSP-224
      IF (DABS(A(I2,J)).LT.1.D-10) A(I2,J)=0.D0                         LCSP-225
   19 CONTINUE                                                          LCSP-226
      IF (NST1.NE.JM) A(I2,NST1+1)=A(7,NST1+1)*A1*B3                    LCSP-227
      NST1=MIN0(NST1+1,JM)                                              LCSP-228
      A1=(B1+1.D0)*WV(5,1)                                              LCSP-229
      B2=2.D0*C2                                                        LCSP-230
      A3=.25D0*(B1+1.D0+B2)**2+WV(5,1)**2                               LCSP-231
      A2=.25D0*(B1+1.D0-B2)*(B1+1.D0+B2)-WV(5,1)**2                     LCSP-232
      B3=(AX(JMX+1)*A2-AX(JMX+2)*A1)/A3                                 LCSP-233
      AX(JMX+2)=(AX(JMX+1)*A1+AX(JMX+2)*A2)/A3                          LCSP-234
      AX(JMX+1)=B3                                                      LCSP-235
      IF (.NOT.LO(18)) GO TO 21                                         LCSP-236
      DO 20 J=1,IJX                                                     LCSP-237
      B2=B2+2.D0                                                        LCSP-238
      A3=.25D0*(B1+1.D0+B2)**2+WV(5,1)**2                               LCSP-239
      A2=.25D0*(B1+1.D0-B2)*(B1+1.D0+B2)-WV(5,1)**2                     LCSP-240
      B3=(AX(JMX+2*J+1)*A2-AX(JMX+2*J+2)*A1)/A3                         LCSP-241
      AX(JMX+2*J+2)=(AX(JMX+2*J+1)*A1+AX(JMX+2*J+2)*A2)/A3              LCSP-242
   20 AX(JMX+2*J+1)=B3                                                  LCSP-243
   21 I1=I2                                                             LCSP-244
      DO 22 J=M3,JM                                                     LCSP-245
      A(3,J)=0.D0                                                       LCSP-246
   22 A(4,J)=0.D0                                                       LCSP-247
   23 CONTINUE                                                          LCSP-248
   24 B3=JMIN+2*M3-1                                                    LCSP-249
      I1=1                                                              LCSP-250
      I2=3                                                              LCSP-251
C COEFFICIENT OF OTHER PL.                                              LCSP-252
      DO 28 L=1,JML                                                     LCSP-253
      A3=0.D0                                                           LCSP-254
      C1=B3                                                             LCSP-255
      IF (LO(18).AND.MOD(L,2).NE.1) GO TO 26                            LCSP-256
      DO 25 J=M3,IPJ                                                    LCSP-257
      A3=A3+(F(1,J,I)*A(I1,J)+F(2,J,I)*A(I1+1,J))/C1                    LCSP-258
   25 C1=C1+2.D0                                                        LCSP-259
      IF (MF(8,I).NE.0.) A3=A3+A3                                       LCSP-260
      AX(L)=AX(L)+A3*XZ                                                 LCSP-261
   26 IF (L.EQ.JML) GO TO 28                                            LCSP-262
      C1=DFLOAT(L)                                                      LCSP-263
      C2=-(C1-1.D0)/C1                                                  LCSP-264
      C3=(2.D0*C1-1.D0)/C1                                              LCSP-265
      I1=I2                                                             LCSP-266
      I2=4-I1                                                           LCSP-267
      A1=0.D0                                                           LCSP-268
      B1=0.D0                                                           LCSP-269
      DO 27 J=M3,JM1                                                    LCSP-270
      A2=A1                                                             LCSP-271
      B2=B1                                                             LCSP-272
      A1=A(I2,J)                                                        LCSP-273
      B1=A(I2+1,J)                                                      LCSP-274
      A(I1,J)=C2*A(I1,J)+(A(7,J)*A2+A(6,J)*A1+A(5,J)*A(I2,J+1))*C3      LCSP-275
   27 A(I1+1,J)=C2*A(I1+1,J)+(A(7,J)*B2+A(6,J)*B1+A(5,J)*A(I2+1,J+1))*C3LCSP-276
      A(I1,JM)=A(7,JM)*A1*C3                                            LCSP-277
      A(I1+1,JM)=A(7,JM)*B1*C3                                          LCSP-278
   28 CONTINUE                                                          LCSP-279
   29 CONTINUE                                                          LCSP-280
      WRITE (65,1002) JK,JML                                            LCSP-281
      IF (LS) WRITE (MW,1003) JK,IJX,JML                                LCSP-282
      IF (.NOT.LS) WRITE (MW,1004) JK,JML                               LCSP-283
      DO 30 LL=1,JML,5                                                  LCSP-284
      L=LL-1                                                            LCSP-285
      LM=MIN0(JML,LL+4)                                                 LCSP-286
   30 WRITE (MW,1005) L,(AX(M),M=LL,LM)                                 LCSP-287
      DO 31 LL=1,JML                                                    LCSP-288
      L=LL-1                                                            LCSP-289
      IF (LS) WRITE (65,1006) JK,L,AX(LL),IJX                           LCSP-290
      IF (.NOT.LS) WRITE (65,1007) JK,L,AX(LL)                          LCSP-291
   31 CONTINUE                                                          LCSP-292
   32 IF (.NOT.LO(81)) GO TO 36                                         LCSP-293
      WRITE (65,1008) JK,IPK                                            LCSP-294
      WRITE (MW,1009) JK,IPK                                            LCSP-295
      DO 33 LL=1,IPK,5                                                  LCSP-296
      L=2*LL-2                                                          LCSP-297
      LM=MIN0(IPK,LL+4)                                                 LCSP-298
   33 WRITE (MW,1005) L,(FCN(M,JK),M=LL,LM)                             LCSP-299
      DO 34 LL=1,IPK                                                    LCSP-300
      L=2*LL-2                                                          LCSP-301
   34 WRITE (65,1010) JK,L,FCN(LL,JK)                                   LCSP-302
      GO TO 36                                                          LCSP-303
   35 WRITE (65,1011) JK                                                LCSP-304
   36 CONTINUE                                                          LCSP-305
      RETURN                                                            LCSP-306
 1000 FORMAT ('<LEGENDRE>',F10.2,1P,D20.8,0P,F10.2,2I5)                 LCSP-307
 1001 FORMAT ('1 COEFFICIENTS OF LEGENDRE POLYNOMIALS DESCRIBING THE CROLCSP-308
     1SS-SECTIONS.')                                                    LCSP-309
 1002 FORMAT (2I5,' COUPLED LEVEL, NUMBER OF VALUES')                   LCSP-310
 1003 FORMAT (/' DIRECT INTERACTION FOR LEVEL',I2,' MULTIPLIED BY (1-COSLCSP-311
     1(THETA))**(',I2,' )',I4,' COEFFICIENTS'/4X,'L',8X,'C(L)',14X,'C(L+LCSP-312
     21)',13X,'C(L+2)',13X,'C(L+3)',13X,'C(L+4)')                       LCSP-313
 1004 FORMAT (/' DIRECT INTERACTION FOR LEVEL',I2,',',I4,' COEFFICIENTS'LCSP-314
     1/4X,'L',8X,'C(L)',14X,'C(L+1)',13X,'C(L+2)',13X,'C(L+3)',13X,'C(L+LCSP-315
     24)')                                                              LCSP-316
 1005 FORMAT (I5,1P,5D19.10)                                            LCSP-317
 1006 FORMAT (2I5,1P,D20.10,5X,'*(1-COS(THETA)**(',I2,' )')             LCSP-318
 1007 FORMAT (2I5,1P,D20.10)                                            LCSP-319
 1008 FORMAT (2I5,' UNCOUPLED LEVEL, NUMBER OF VALUES')                 LCSP-320
 1009 FORMAT (/' COMPOUND NUCLEUS FOR LEVEL',I2,',',I4,' COEFFICIENTS (OLCSP-321
     1NLY EVEN ONES ARE GIVEN)'/4X,'L',8X,'C(L)',14X,'C(L+2)',13X,'C(L+4LCSP-322
     2)',13X,'C(L+6)',13X,'C(L+8)')                                     LCSP-323
 1010 FORMAT (2I5,1P,D20.10,5X,'COMPOUND NUCLEUS')                      LCSP-324
 1011 FORMAT (I5,4X,'0 CLOSED CHANNEL')                                 LCSP-325 
      END                                                               LCSP-326
