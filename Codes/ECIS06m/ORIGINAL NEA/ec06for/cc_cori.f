C 28/02/07                                                      ECIS06  CORI-000
      SUBROUTINE CORI(EI,EF,SI,SF,T,SSI,SSF,FGI,FGF,Z,LM1,LM2,LM3,LM,W) CORI-001
C    COMPUTATION OF THE INTEGRALS FROM RS TO INFINITY OF PRODUCTS OF    CORI-002
C REGULAR OR IRREGULAR COULOMB FUNCTIONS OF SAME L-VALUE DIVIDED BY R**2CORI-003
C FOR L=0 TO L=LM-1 AND OF REGULAR FUNCTIONS OF SAME L DIVIDED BY R**2  CORI-004
C L=0 TO L=LM2-1.CALLS CORZ IF THERE IS NO CHARGE.                      CORI-005
C INPUT:     EI,EF:   COULOMB PARAMETERS.                               CORI-006
C            SI,SF:   WAVE NUMBERS MULTIPLIED BY MATCHING RADIUS.       CORI-007
C            SSI,SSF: COULOMB PHASE-SHIFTS FOR L=0.                     CORI-008
C            FGI,FGF: REGULAR AND IRREGULAR COULOMB FUNCTIONS.          CORI-009
C            LM1:     MAXIMUM NUMBER OF INTEGRALS FROM RS TO INF.       CORI-010
C            LM2:     NUMBER OF INTEGRALS FROM 0 TO INFINITY            CORI-011
C            LM3:     DIMENSION OF WORKING AREA.                        CORI-012
C            LM:      ACTUAL NUMBER OF INTEGRALS FROM RS TO INF.        CORI-013
C OUTPUT:    T:       INTEGRALS FROM 0 TO INFINITY FOR PRODUCTS OF      CORI-014
C                           REGULAR FUNCTIONS.                          CORI-015
C            W:       INTEGRALS FROM RS TO INFINITY.                    CORI-016
C WORKING AREA:                                                         CORI-017
C            Z:       FOR THE COMPUTATION OF HYPERGEOMETRIC F3.         CORI-018
C                                                                       CORI-019
C THE SUBROUTINES CORI, CORO, CORZ AND CORA FOR COULOMB CORRECTIONS     CORI-020
C ARE IDENTICAL IN THE CODES ECIS AND DWBA.                             CORI-021
C THE PRESENT VERSION IS OF FEBRUARY 2006. FOR DWBA, THE NO-CHARGE      CORI-022
C PART AND ODD INTEGRALS COULD BE SUPPRESSED.                           CORI-023
C***********************************************************************CORI-024
      IMPLICIT REAL*8 (A-H,O-Z)                                         CORI-025
      LOGICAL LL(4)                                                     CORI-026
      DIMENSION FGI(LM1,*),FGF(LM1,*),Z(4,*),W(LM1,*),X1(4),X2(4),X3(4),CORI-027
     1X4(4),T(*),Y(4,4),AZ(2),CC(4,4)                                   CORI-028
      COMMON /INOUT/ MR,MW,MS                                           CORI-029
      DATA PI /3.1415926535897932D0/                                    CORI-030
      IF ((EI.EQ.0.D0).AND.(EF.EQ.0.D0)) GO TO 38                       CORI-031
C COMPUTATION OF SOME CONSTANTS.                                        CORI-032
      IF (SI.GT.2.2D0*SF) GO TO 39                                      CORI-033
      RS=DSQRT(SI*SF)                                                   CORI-034
      FS=DSQRT(SF/SI)                                                   CORI-035
      FI=EI**2                                                          CORI-036
      FF=EF**2                                                          CORI-037
      EIF=EI*EF                                                         CORI-038
      EX=EF-EI                                                          CORI-039
      SIF=SI/SF+SF/SI                                                   CORI-040
      CX=FI*SI**2+FF*SF**2                                              CORI-041
      DX=FI*SI**2-FF*SF**2                                              CORI-042
      IF (DABS(DX).GT.1.D-10) GO TO 1                                   CORI-043
      DT=1.D0                                                           CORI-044
      DX=0.D0                                                           CORI-045
      DY=EIF                                                            CORI-046
      GO TO 2                                                           CORI-047
    1 DT=SI**2-SF**2                                                    CORI-048
      DY=(FI*SI**4-FF*SF**4)/(SI*SF)                                    CORI-049
      IF (DABS(DT).LT.1.D-10) GO TO 2                                   CORI-050
      DX=DX/DT                                                          CORI-051
      DY=DY/DT                                                          CORI-052
      DT=1.D0                                                           CORI-053
    2 RM=DMAX1((5.D0*DABS(EI)+22.5D0)*FS,(5.D0*DABS(EF)+22.5D0)/FS,.18D0CORI-054
     1*EIF)/3.D0                                                        CORI-055
      IF (RM.LT.RS) RM=RS                                               CORI-056
C COMPUTATION OF INTEGRALS FROM RS TO INFINITY FOR L=0 AND L=1.         CORI-057
      CALL COR0(EI,EF,SSI,SSF,ZR,ZI,Z,LM3,RM,RS,FS,FI,FF,EX,CC)         CORI-058
C COMPUTATION OF INTEGRALS FROM 0 TO INFINITY FOR L=0 AND L=1.          CORI-059
      IF (EI.NE.EF.OR.SI.NE.SF) GO TO 4                                 CORI-060
C WHEN EI=EF=E   I(L,L,2)=(PI-PI*COTH(PI*E)-1/E+SUM ON N FROM 0 TO L OF CORI-061
C                2*E/(N**2+E**2))/2/(2*L+1).                            CORI-062
      A2=PI*.5D0                                                        CORI-063
      B2=DEXP(-2.D0*PI*EI)                                              CORI-064
      IF (EI.NE.0.D0) A2=A2-.5D0*PI*(B2+1.D0)/(1.D0-B2)+.5D0/EI         CORI-065
      B3=0.D0                                                           CORI-066
      DO 3 I=1,LM2                                                      CORI-067
      B2=B3                                                             CORI-068
      B3=B3+1.D0                                                        CORI-069
      T(I)=A2/(B2+B3)                                                   CORI-070
    3 A2=A2+EI/(B3**2+FI)                                               CORI-071
      LL(1)=.TRUE.                                                      CORI-072
      LL(2)=.TRUE.                                                      CORI-073
      GO TO 13                                                          CORI-074
C COMPUTATION OF I(0,0,1) AND I(1,1,1)                                  CORI-075
C  I(L,L,1) = 2*(SI*SF)**(1/2)/(SI+SF)**2 * EXP( PI*SIGN(SF-SI)*EX/2) * CORI-076
C REAL PART OF (|SF-SI|/(SF+SI))**(I*EX)*GAMMA(-I*EX)*EXP(I*(XEF-XEI))* CORI-077
C   2F1(L+1-I*EI,L+1+I*EF;1+I*EF-I*EI;((SF-SI)/(SF+SI))**2)             CORI-078
C  WHERE EX=EF-EI.                                                      CORI-079
    4 A2=EX**2                                                          CORI-080
      B2=((SI-SF)/(SI+SF))**2                                           CORI-081
      A3=2.D0*DEXP(ZR-PI*.5D0*EX)/(SI+SF)**2*RS                         CORI-082
      IF (SF.GT.SI) A3=A3*DEXP(PI*EX)                                   CORI-083
      A4=ZI+SSF-SSI+.5D0*DLOG(B2)*EX                                    CORI-084
      A5=A3*DCOS(A4)                                                    CORI-085
      A6=A3*DSIN(A4)                                                    CORI-086
      A3=A5                                                             CORI-087
      A4=A6                                                             CORI-088
      A7=0.D0                                                           CORI-089
      DO 7 L=1,2                                                        CORI-090
      AZ(L)=A3                                                          CORI-091
      A8=0.D0                                                           CORI-092
      DO 5 N=1,500                                                      CORI-093
      A8=A8+1.D0                                                        CORI-094
      B3=(A7+A8)**2+EIF                                                 CORI-095
      B4=(A7+A8)*EX                                                     CORI-096
      B7=B2/(A8*(A8**2+A2))                                             CORI-097
      B5=(B3*A8+B4*EX)*B7                                               CORI-098
      B6=(A8*B4-EX*B3)*B7                                               CORI-099
      A9=A3*B5-A4*B6                                                    CORI-100
      A4=A3*B6+A4*B5                                                    CORI-101
      A3=A9                                                             CORI-102
      AZ(L)=AZ(L)+A3                                                    CORI-103
      IF (DABS(A3)+DABS(A4).LT.DABS(AZ(L))*1.D-16) GO TO 6              CORI-104
    5 CONTINUE                                                          CORI-105
      LI=L-1                                                            CORI-106
      WRITE (MW,1000) LI,LI,EI,EF,AZ(L),A3,A4                           CORI-107
    6 B3=1.D0+EIF                                                       CORI-108
      B5=(1.D0-B2)/DSQRT(B3*B3+EX*EX)                                   CORI-109
      A7=1.D0                                                           CORI-110
      A3=B5*(A5*B3-A6*EX)                                               CORI-111
    7 A4=B5*(A6*B3+A5*EX)                                               CORI-112
C COMPUTATION OF I(0,0,2) AND I(1,1,2) . ONLY THE FIRST IS NEEDED       CORI-113
C   FOR BACKWARD RECURRENCE:                                            CORI-114
C I(L,L,2) = PI/(2*SH(PI*EX)*( (Y/SI-1/(Y*SF))/(2*L+1)+((SI-SF)/(SF*SI))CORI-115
C  *EXP(-PI*EX/2) * REAL PART OF EXP(I*(XEI-XEF))*((SI-SF)/2)**(-I*EX)  CORI-116
C  *SF**(I*EF)*SI**(-I*EI)*F3(-L+I*EI,-L-I*EF,L+1+I*EI,L+1-I*EF,2-I*EX; CORI-117
C  (SI-SF)/(2*SI),(SF-SI)/(2*SF))                                       CORI-118
C   WHERE  Y=EXP(PI*EX/2)*(SI/SF)**L*|GAMMA(L+1+I*EF)/GAMMA(L+1+I*EI)|. CORI-119
      IF (EI.GT.EF) GO TO 8                                             CORI-120
      EEF=EF                                                            CORI-121
      EEI=EI                                                            CORI-122
      SXF=SF                                                            CORI-123
      SXI=SI                                                            CORI-124
      SSS=SSI-SSF-ZI                                                    CORI-125
      GO TO 9                                                           CORI-126
    8 EEF=EI                                                            CORI-127
      EEI=EF                                                            CORI-128
      SXF=SI                                                            CORI-129
      SXI=SF                                                            CORI-130
      SSS=SSF-SSI+ZI                                                    CORI-131
    9 EEX=EEF-EEI                                                       CORI-132
      SFI=SXI-SXF                                                       CORI-133
      A2=EEX/(2.D0*EEF)                                                 CORI-134
      A3=-A2*SFI/(2.D0*SXF)                                             CORI-135
      A4=A2*SFI/(2.D0*SXI)                                              CORI-136
      LL(1)=DFLOAT(LM2)*DABS(EX).GT.3.D0*DABS(EI+EF)                    CORI-137
      LL(2)=.NOT.LL(1)                                                  CORI-138
      B3=DEXP(-ZR)/EEX/DSQRT(1.D0+EX**2)                                CORI-139
      IF (SFI.LT.0.D0) B3=B3*DEXP(EEX*PI)                               CORI-140
      B4=SSS+DATAN(EEX)-.5D0*PI+EEF*DLOG(SXF)-EEI*DLOG(SXI)-EEX*DLOG(.5DCORI-141
     10*DABS(SFI))                                                      CORI-142
      A5=B3*DCOS(B4)                                                    CORI-143
      A6=B3*DSIN(B4)                                                    CORI-144
      B5=A5                                                             CORI-145
      B6=A6                                                             CORI-146
      B2I=-2.D0*PI+4.D0*PI**2*EEI                                       CORI-147
      B2F=-2.D0*PI+4.D0*PI**2*EEF                                       CORI-148
      IF (DABS(EEI).GT.1.D-8) B2I=(DEXP(-2.D0*PI*EEI)-1.D0)/EEI         CORI-149
      IF (DABS(EEF).GT.1.D-8) B2F=(DEXP(-2.D0*PI*EEF)-1.D0)/EEF         CORI-150
      B2=DSQRT(B2I/B2F)                                                 CORI-151
      A7=DEXP(-.5D0*PI*EEX)                                             CORI-152
      A8=PI*A7**2/(1.D0-A7**4)*RS                                       CORI-153
      A7=-A7*SFI/(SF*SI)                                                CORI-154
      Z(1,1)=1.D0                                                       CORI-155
      Z(2,1)=0.D0                                                       CORI-156
      Z(3,1)=1.D0                                                       CORI-157
      Z(4,1)=0.D0                                                       CORI-158
      DO 12 L=1,2                                                       CORI-159
      LI=L-1                                                            CORI-160
      T(L)=B5                                                           CORI-161
      N=0                                                               CORI-162
   10 N=N+1                                                             CORI-163
      IF (2*N.GT.LM3) CALL MEMO('CORI',LM3,2*N)                         CORI-164
      B4=DFLOAT(N+1)**2+EX**2                                           CORI-165
      B9=(B5*DFLOAT(N+1)-EEX*B6)/(B4*A2)                                CORI-166
      B6=(B5*EEX+B6*DFLOAT(N+1))/(B4*A2)                                CORI-167
      B5=B9                                                             CORI-168
      B3=DFLOAT((L+N-1)*(N-L))-EEI**2                                   CORI-169
      B4=EEI*DFLOAT(2*N-1)                                              CORI-170
      Z(1,N+1)=(Z(1,N)*B3-Z(2,N)*B4)*A4/DFLOAT(N)                       CORI-171
      Z(2,N+1)=(Z(1,N)*B4+Z(2,N)*B3)*A4/DFLOAT(N)                       CORI-172
      B3=DFLOAT((L+N-1)*(N-L))-EEF**2                                   CORI-173
      B4=-EEF*DFLOAT(2*N-1)                                             CORI-174
      Z(3,N+1)=(Z(3,N)*B3-Z(4,N)*B4)*A3/DFLOAT(N)                       CORI-175
      Z(4,N+1)=(Z(3,N)*B4+Z(4,N)*B3)*A3/DFLOAT(N)                       CORI-176
      B7=Z(1,N+1)                                                       CORI-177
      B8=Z(2,N+1)                                                       CORI-178
      DO 11 J=1,N                                                       CORI-179
      M=N+2-J                                                           CORI-180
      B7=B7+Z(1,J)*Z(3,M)-Z(2,J)*Z(4,M)                                 CORI-181
   11 B8=B8+Z(2,J)*Z(3,M)+Z(1,J)*Z(4,M)                                 CORI-182
      IF (DABS(B7)+DABS(B8).GT.1.D30) GO TO 40                          CORI-183
      B7=B7*B5-B8*B6                                                    CORI-184
      T(L)=T(L)+B7                                                      CORI-185
      IF (DABS(B7).GT.DABS(T(L))*1.D-16) GO TO 10                       CORI-186
      T(L)=A8*(A7*T(L)+(B2/SXI-1.D0/(B2*SXF))/DFLOAT(2*L-1))            CORI-187
      IF (LL(1)) GO TO 13                                               CORI-188
      B2=B2*SXF/SXI*DSQRT((1.D0+EEF**2)/(1.D0+EEI**2))                  CORI-189
      B4=DSQRT((1.D0+EEI**2)*(1.D0+EEF**2))                             CORI-190
      B3=(1.D0+EIF)/B4                                                  CORI-191
      B4=-EEX/B4                                                        CORI-192
      B5=B3*A5-B4*A6                                                    CORI-193
   12 B6=B3*A6+B4*A5                                                    CORI-194
C COMP. OF THE INTEGRALS FROM 0 TO INFINITY.                            CORI-195
C UPWARDS RECURRENCES FOR THE OTHER INTEGRALS:  STARTING VALUES:        CORI-196
   13 DO 14 J=1,4                                                       CORI-197
      W(1,J)=CC(J,3)                                                    CORI-198
      W(2,J)=CC(J,4)                                                    CORI-199
      X2(J)=CC(J,1)/RS                                                  CORI-200
   14 X3(J)=CC(J,2)/RS                                                  CORI-201
      RM=RS**2                                                          CORI-202
      IM=1                                                              CORI-203
      LX=LM                                                             CORI-204
      IF (.NOT.LL(1)) LX=MAX0(LX,LM2)                                   CORI-205
      LL(3)=LL(1)                                                       CORI-206
      I=2                                                               CORI-207
      A1=DSQRT(1.D0+FI)                                                 CORI-208
      B1=DSQRT(1.D0+FF)                                                 CORI-209
      IF (LL(1)) GO TO 15                                               CORI-210
      C4=AZ(1)                                                          CORI-211
      C5=AZ(2)                                                          CORI-212
   15 I=I+1                                                             CORI-213
      A2=A1                                                             CORI-214
      B2=B1                                                             CORI-215
      A=DFLOAT(I-1)                                                     CORI-216
      A1=DSQRT(A**2+FI)                                                 CORI-217
      B1=DSQRT(A**2+FF)                                                 CORI-218
      LL(4)=I.GT.LM                                                     CORI-219
   16 A9=(A-1.D0)**2*DT+DX                                              CORI-220
      A3=(2.D0*A-3.D0)*(A**2*DT+DX)*A2*B2                               CORI-221
      A4=-(2.D0*A-1.D0)*(DX*CX/RM+DY*((A-1.D0)**2+A**2)+DT*(A*A-A)**2*SICORI-222
     1F)                                                                CORI-223
      A5=(2.D0*A+1.D0)*A9*A1*B1                                         CORI-224
      C1=DT*(2.D0*A-1.D0)*(EI*SI+EF*SF)*A2*B2/(A-1.D0)                  CORI-225
      C2=DT*(A-.5D0)*(B2**2*((A-1.D0)*EF*SF-(A+1.D0)*EI*SI)*SF/SI+A2**2*CORI-226
     1((A-1.D0)*EI*SI-(A+1.D0)*EF*SF)*SI/SF)/(A-1.D0)                   CORI-227
      A6=A*A2*B2                                                        CORI-228
      A7=-(2.D0*A-1.D0)*(EIF+(A-1.D0)*A*SIF*.5D0)                       CORI-229
      A8=(A-1.D0)*A1*B1                                                 CORI-230
      IF (LL(4)) GO TO 17                                               CORI-231
      B6=(A-.5D0)/(RM*RS*(A-1.D0)**2)                                   CORI-232
      B3=-2.D0*RM*DX*B6*A2*B2                                           CORI-233
      B4=B6*A9*SF*(2.D0*(A*A-A)*A+(A+1.D0)*EI*SI-(A-1.D0)*EF*SF)*B2     CORI-234
      B5=B6*A9*SI*(2.D0*(A*A-A)*A+(A+1.D0)*EF*SF-(A-1.D0)*EI*SI)*A2     CORI-235
      B6=-B6*(A9*(2.D0*(A*A-A)**2*(2.D0*A-1.D0)+(EI*SI+EF*SF)*(A-1.D0)*(CORI-236
     1A+1.D0)*(2.D0*A-1.D0)+2.D0*EIF*RM*A)-(A-1.D0)*(CX*DX+(A-1.D0)*((A+CORI-237
     21.D0)*CX*DT-2.D0*RM*DY)))                                         CORI-238
      B9=(A-.5D0)*A/RS                                                  CORI-239
      B7=B9*B2/SI                                                       CORI-240
      B8=B9*A2/SF                                                       CORI-241
      B9=B9*(-(2.D0*A-1.D0)*(A-1.D0)-(EI*SI+EF*SF))/RM                  CORI-242
   17 IF (IM.NE.1) GO TO 26                                             CORI-243
      IF (LL(3)) GO TO 18                                               CORI-244
      T(I)=-(A3*T(I-2)+A4*T(I-1)+C1*C4+C2*C5)/A5                        CORI-245
      C3=C4                                                             CORI-246
      C4=C5                                                             CORI-247
      C5=-(A6*C3+A7*C4)/A8                                              CORI-248
   18 IF (LL(4)) GO TO 23                                               CORI-249
      DO 20 J=1,4                                                       CORI-250
      J1=MOD(J+1,2)+I-2                                                 CORI-251
      J2=(J+1)/2+I-3                                                    CORI-252
      DO 19 K=1,4                                                       CORI-253
      K1=MOD(K+1,2)*2+1                                                 CORI-254
      K2=2*((K+1)/2)-1                                                  CORI-255
   19 Y(K,J)=FGI(J1,K1)*FGF(J2,K2)                                      CORI-256
   20 CONTINUE                                                          CORI-257
C RECURRENCE FOR (0,2;L).                                               CORI-258
      DO 21 J=1,4                                                       CORI-259
      X1(J)=X2(J)                                                       CORI-260
      X2(J)=X3(J)                                                       CORI-261
   21 W(I,J)=-(B3*Y(J,1)+B4*Y(J,2)+B5*Y(J,3)+B6*Y(J,4)+A3*W(I-2,J)+A4*W(CORI-262
     1I-1,J)+C1*X1(J)+C2*X2(J))/A5                                      CORI-263
      DO 22 J=1,4                                                       CORI-264
   22 X3(J)=-(B7*Y(J,2)+B8*Y(J,3)+B9*Y(J,4)+A6*X1(J)+A7*X2(J))/A8       CORI-265
   23 IF (I.LT.LX) GO TO 15                                             CORI-266
      LN=MAX0(LM,IDINT(RS))+50                                          CORI-267
      LX=LN                                                             CORI-268
      IF (.NOT.LL(2)) LX=MAX0(LX,LM2+6*MIN0(IDINT((SI+SF)/DABS(SI-SF)),5CORI-269
     100))                                                              CORI-270
      LL(3)=LL(2)                                                       CORI-271
      IM=2                                                              CORI-272
      X4(1)=1.D0                                                        CORI-273
      X4(2)=0.D0                                                        CORI-274
      X4(3)=0.D0                                                        CORI-275
      X4(4)=0.D0                                                        CORI-276
      A1=LX-1                                                           CORI-277
      A2=DSQRT(A1**2+FI)                                                CORI-278
      B2=DSQRT(A1**2+FF)                                                CORI-279
      DO 24 J=1,4                                                       CORI-280
      X2(J)=0.D0                                                        CORI-281
   24 X3(J)=0.D0                                                        CORI-282
      X3(1)=1.D-20                                                      CORI-283
      X3(2)=1.D-20                                                      CORI-284
   25 LX=LX-1                                                           CORI-285
      A=LX                                                              CORI-286
      A1=A2                                                             CORI-287
      B1=B2                                                             CORI-288
      A3=DFLOAT(LX-1)                                                   CORI-289
      LL(4)=LX.GT.LN                                                    CORI-290
      A2=DSQRT(A3**2+FI)                                                CORI-291
      B2=DSQRT(A3**2+FF)                                                CORI-292
      GO TO 16                                                          CORI-293
   26 IF (LL(3)) GO TO 27                                               CORI-294
      X4(1)=-A6/(A7+A8*X4(1))                                           CORI-295
      X4(3)=X4(1)*X4(2)                                                 CORI-296
      X4(2)=-(A4*X4(3)+A5*X4(4)*X4(1)+C1+C2*X4(1))/A3                   CORI-297
      X4(4)=X4(3)                                                       CORI-298
      IF (LX.LE.LM2) T(LX)=X4(3)/X4(2)                                  CORI-299
   27 IF (LL(4)) GO TO 33                                               CORI-300
C DOWNWARD RECURRENCE FOR F(EI) AND F(EF).                              CORI-301
      DO 28 J=1,4                                                       CORI-302
      X1(J)=X2(J)                                                       CORI-303
   28 X2(J)=X3(J)                                                       CORI-304
      X3(1)=((2.D0*A-1.D0)*(EI+(A-1.D0)*A/SI)*X2(1)-(A-1.D0)*A1*X1(1))/(CORI-305
     1A*A2)                                                             CORI-306
      X3(2)=((2.D0*A-1.D0)*(EF+(A-1.D0)*A/SF)*X2(2)-(A-1.D0)*B1*X1(2))/(CORI-307
     1A*B2)                                                             CORI-308
C DOWNWARDS RECURRENCE FOR (1,2;L).                                     CORI-309
      X3(4)=-(A7*X2(4)+A8*X1(4)+B7*X2(1)*X3(2)+B8*X3(1)*X2(2)+B9*X2(1)*XCORI-310
     12(2))/A6                                                          CORI-311
C DOWNWARD RECURRENCE FOR (0,2;L).                                      CORI-312
      X3(3)=-(A4*X2(3)+A5*X1(3)+C1*X3(4)+C2*X2(4)+B3*X3(1)*X3(2)+B4*X2(1CORI-313
     1)*X3(2)+B5*X3(1)*X2(2)+B6*X2(1)*X2(2))/A3                         CORI-314
      I=LX-1                                                            CORI-315
      IF (I.LE.LM) W(I,1)=X3(3)                                         CORI-316
      IF (DABS(X3(1)).LT.1.D10) GO TO 29                                CORI-317
      X3(1)=X3(1)*1.D-20                                                CORI-318
      X2(1)=X2(1)*1.D-20                                                CORI-319
      GO TO 30                                                          CORI-320
   29 IF (DABS(X3(2)).LT.1.D10) GO TO 33                                CORI-321
      X2(2)=X2(2)*1.D-20                                                CORI-322
      X3(2)=X3(2)*1.D-20                                                CORI-323
   30 DO 31 J=3,4                                                       CORI-324
      X3(J)=X3(J)*1.D-20                                                CORI-325
   31 X2(J)=X2(J)*1.D-20                                                CORI-326
      IF (I.GT.LM) GO TO 33                                             CORI-327
      DO 32 J=I,LM                                                      CORI-328
   32 W(J,1)=W(J,1)*1.D-20                                              CORI-329
   33 IF (LX.GT.2) GO TO 25                                             CORI-330
      IF (LL(2)) GO TO 36                                               CORI-331
C COMPUTATION OF THE MIXTURE OF THE DECREASING SOLUTION OF HOMOGENEOUS  CORI-332
C EQUATION  (S2/S1)**L * |GAMMA(L+1+I*E2)/GAMMA(L+1+I*E1)|/(2*L+1)      CORI-333
C WHERE S2 IS THE SMALLER OF THE TWO VALUES SI AND SF.                  CORI-334
      A1=X4(2)*AZ(1)                                                    CORI-335
      A2=T(1)-A1                                                        CORI-336
      A3=0.D0                                                           CORI-337
      DO 35 I=2,LM2                                                     CORI-338
      A3=A3+1.D0                                                        CORI-339
      A1=T(I)*A1                                                        CORI-340
      A4=SF*DSQRT((A3**2+FF)/(A3**2+FI))/SI                             CORI-341
      IF (SI.GT.SF) GO TO 34                                            CORI-342
      A2=A2/A4                                                          CORI-343
      GO TO 35                                                          CORI-344
   34 A2=A2*A4                                                          CORI-345
   35 T(I)=A1+A2/(2.D0*A3+1.D0)                                         CORI-346
C THE VALUES OBTAINED ARE UNNORMALISED INTEGRALS FROM RS TO 0.          CORI-347
C NORMALISATION AND TRANSFORMATION TO INTEGRAL FROM RS TO INFINITY.     CORI-348
   36 A1=FGI(1,1)*FGF(1,1)/(X3(1)*X3(2))                                CORI-349
      DO 37 J=1,LM                                                      CORI-350
   37 W(J,1)=A1*W(J,1)+T(J)                                             CORI-351
      IF (SI.NE.SF) WRITE (MW,1001) CC(1,3),W(1,1),CC(1,4),W(2,1)       CORI-352
      IF ((DABS(CC(1,3)-W(1,1))/(DABS(CC(1,3))+DABS(W(1,1))).GT.1.D-4).OCORI-353
     1R.(DABS(CC(1,4)-W(2,1))/(DABS(CC(1,4))+DABS(W(2,1))).GT.1.D-4)) GOCORI-354
     2 TO 41                                                            CORI-355
      RETURN                                                            CORI-356
   38 CALL CORZ(SI,SF,T,FGI,FGF,LM1,LM2,LM,W)                           CORI-357
      RETURN                                                            CORI-358
   39 WRITE (MW,1002) SF,SI                                             CORI-359
      GO TO 42                                                          CORI-360
   40 WRITE (MW,1003) N,LI,LI,EI,EF,T(L),B5,B6                          CORI-361
      GO TO 42                                                          CORI-362
   41 WRITE (MW,1004)                                                   CORI-363
      IF (SI.EQ.SF) WRITE (MW,1001) CC(1,3),W(1,1),CC(1,4),W(2,1)       CORI-364
   42 WRITE (MW,1005)                                                   CORI-365
      STOP                                                              CORI-366
 1000 FORMAT (' NO CONVERGENCE WITH 500 TERMS FOR I(',I1,',',I1,') WITH CORI-367
     1EI =',F15.6,31X,'AND EF =',F15.6,5X,'IN CORI'/5X,'VALUE =',D15.8,5CORI-368
     2X,'LAST TERM =',2D16.8/' ... CONTINUE ...')                       CORI-369
 1001 FORMAT (' INTEGRALS WITH REGULAR FUNCTIONS: (L+1)       DIRECT:   CORI-370
     1   BACKWARDS RECURRENCE:'/37X,'1',2D20.10/37X,'2',2D20.10)        CORI-371
 1002 FORMAT (' TOO LARGE RATIO SF/SI =',F15.6,1H/,F15.6,5X,'IN CORI.') CORI-372
 1003 FORMAT (' THE',I4,' TERM IS TOO LARGE IN THE COMPUTATION OF I(',I1CORI-373
     1,',',I1,',2)    WITH EI =',F15.6,5X,'AND EF =',F15.6,5X,'IN CORH:'CORI-374
     2/5X,'VALUE =',D15.8,5X,'LAST TERM =',2D16.8)                      CORI-375
 1004 FORMAT (' NO COINCIDENCE BETWEEN DIRECT AND BACKWARDS RECURRENCES CORI-376
     1FOR INTEGRALS WITH REGULAR FUNCTIONS.')                           CORI-377
 1005 FORMAT (//' IN CORI  ...  STOP  ...')                             CORI-378
      END                                                               CORI-379
