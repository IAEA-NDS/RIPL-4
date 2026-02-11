C 02/03/07                                                      ECIS06  SCHE-000
      SUBROUTINE SCHE(F,JMAX,KMAX,IPI,MT1,MT2,MF,TX,BM,MC,FA,XG,LMAX1,WVSCHE-001
     1,KAB,KBA,KCB,JMIN,IPJ,IPK,FN,NCOLL,NCOLS,NCT,FGN,AM,JIT,JTI,NLT,IDSCHE-002
     21,LO)                                                             SCHE-003
C SCATTERING COEFFICIENTS IN THE HELICITY REPRESENTATION.               SCHE-004
C INPUT:     F:       S-MATRIX.                                         SCHE-005
C            JMAX:    MAXIMUM NUMBER OF CHANNEL SPINS, DIMENSION FOR F. SCHE-006
C            KMAX:    MAXIMUM NUMBER OF COMPOUND NUCLEUS, DIMENSION FOR SCHE-007
C                     FN.                                               SCHE-008
C            IPI(I,*):PARITY OF CHANNEL FOR I=1,                        SCHE-009
C                     MULTIPLICITIES OF PARTICLE FOR I=2,               SCHE-010
C                     MULTIPLICITIES OF TARGET FOR I=3,                 SCHE-011
C                     PRODUCT OF CHARGES FOR I=4,                       SCHE-012
C                     FIRST/LAST CHANNEL NUMBER FOR I=6, 7 (SEE DEPH),  SCHE-013
C                     MAXIMUM ANGULAR MOMENTUM FOR I=11.                SCHE-014
C            MT1,MT2: MAXIMUM 2*SPIN+1 FOR PART. AND TARGET.            SCHE-015
C            MF:      HELICITY NUMBERS (SEE DEPH).                      SCHE-016
C            TX:      TOTAL REACTION CROSS SECTION IN MILLIBARNS        SCHE-017
C                     FOLLOWED BY THE TOTAL CROSS SECTION FOR EACH      SCHE-018
C                     LEVEL, THE COMPOUND NUCLEUS CROSS SECTIONS,       SCHE-019
C                     THE FISSION AND THE GAMMA CROSS SECTIONS.         SCHE-020
C            XG:      COULOMB PHASE-SHIFTS.                             SCHE-021
C            LMAX1:   DIMENSION FOR XG.                                 SCHE-022
C            WV:      WAVE NUMBER AND COULOMB PARAMETER (SEE COLF).     SCHE-023
C            KAB:     DIMENSION FOR FA.                                 SCHE-024
C            KBA:     NUMBER OF INDEPENDENT AMPLITUDES.                 SCHE-025
C            KCB:     DIMENSION FOR MC.                                 SCHE-026
C            JMIN:    TWICE MINIMUM CHANNEL SPIN.                       SCHE-027
C            IPJ:     NUMBER OF THE CHANNEL SPIN.                       SCHE-028
C            IPK:     NUMBER OF L VALUES FOR COMPOUND NUCLEUS.          SCHE-029
C            FN:      COMPOUND NUCLEUS CONTRIBUTION.                    SCHE-030
C            NCOLL:   NUMBER OF COUPLED LEVELS.                         SCHE-031
C            NCOLS:   NUMBER OF LEVELS WITH ANGULAR DISTRIBUTION.       SCHE-032
C            NCT:     NUMBER OF EQUATIONS AND SOLUTIONS FOR EACH PARITY.SCHE-033
C            FGN:     COEFFICIENTS OF LEGENDRE POLYNOMIALS FOR COMPOUND SCHE-034
C                     NUCLEUS.                                          SCHE-035
C            JIT:     NUMBER OF DIFFERENT RATES OF INTERPOLATION.       SCHE-036
C            JTI:     LIMITS AND STEPS OF INTERPOLATION.                SCHE-037
C            NLT:     MEMORIES NEEDED FOR LEGENDRE POLYNOMIALS.         SCHE-038
C            ID1:     LENGTH AVAILABLE FOR BM.                          SCHE-039
C            LO(I):   LOGICAL CONTROLS:                                 SCHE-040
C               LO(18) =.TRUE. PROJECTILE-TARGET ANTISYMMETRISATION.    SCHE-041
C               LO(41) =.TRUE. FACTORISATION OF 1/(1-COS(THETA)).       SCHE-042
C               LO(43) =.TRUE. INTERPOLATION ON TOTAL SPIN.             SCHE-043
C               LO(56) =.TRUE. OUTPUT OF S-MATRIX ELEMENTS.             SCHE-044
C               LO(60) =.TRUE. S-MATRIX ELEMENTS WRITTEN ON FILE 60.    SCHE-045
C               LO(65) =.TRUE. PRINT COEFFICIENTS OF THE EXPANSION IN   SCHE-046
C                              LEGENDRE POLYNOMIALS ON FILE 65.         SCHE-047
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             SCHE-048
C               LO(82) =.TRUE. OLD SIMPLIFIED COMPOUND NUCLEUS.         SCHE-049
C               LO(116)=.TRUE. NO OUTPUT.                               SCHE-050
C               LO(119)=.TRUE. RESULTS WITH THE LAST CALCULATION.       SCHE-051
C               LO(120)=.TRUE. OUTPUT AND LAST CALCULATION BEST ONE.    SCHE-052
C               LO(122)=.TRUE. IDENTICAL PARTICLES WITHOUT SPIN.        SCHE-053
C               LO(132)=.TRUE. STORE FISSION AND GAMMA TRANSMISSION     SCHE-054
C                              COEFFICIENTS FOR INTERPOLATION.          SCHE-055
C OUTPUT:    F:       HELICITY SCATTERING COEFFICIENTS.                 SCHE-056
C            TX,FGN:  ADDITION OF INTERPOLATED VALUES.                  SCHE-057
C WORKING AREAS:                                                        SCHE-058
C            BM:      FOR FACTORS 1/(-X*COS(THETA)), 3J COEFFICIENTS,...SCHE-059
C            MC:      NUCLEAR STATE NUMBERS AND ANGULAR MOMENTA.        SCHE-060
C            FA:      FOR STORAGE OF S-MATRIX FOR A GIVEN TOTAL SPIN.   SCHE-061
C            AM:      FOR PRODUCTS OF C.G. COEFFICIENTS.                SCHE-062
C                                                                       SCHE-063
C FOR THE COMMON  /ANGUL/ SEE LECT.                                     SCHE-064
C FOR THE COMMON  /DCONS/ SEE CALC.                                     SCHE-065
C FOR THE COMMON  /NCOMP/ SEE CALX.                                     SCHE-066
C                                                                       SCHE-067
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     SCHE-068
C  XZ:        CONVERSION FACTOR TO MILLIBARNS.                          SCHE-069
C   USED:     XZ.                                                       SCHE-070
C                                                                       SCHE-071
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ANGUL/:                     SCHE-072
C  NCJ:       NUMBER OF FACTORISATIONS OF 1/(1-COS(THETA)) IN AMPLITUDE.SCHE-073
C   DEFINED:  NCJ.                                                      SCHE-074
C   USED:     NCJ.                                                      SCHE-075
C                                                                       SCHE-076
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     SCHE-077
C  NSP(1):    NUMBER OF UNCOUPLED LEVELS FOR COMPOUND NUCLEUS           SCHE-078
C             INCLUDING DISCRETISATION OF CONTINUUM.                    SCHE-079
C  NSP(2):    NUMBER OF THESE LEVELS WITH ANGULAR DISTRIBUTION.         SCHE-080
C  NSP(3):    NUMBER OF THESE LEVELS WITHOUT ANGULAR DISTRIBUTION.      SCHE-081
C   USED:     NSP.                                                      SCHE-082
C                                                                       SCHE-083
C***********************************************************************SCHE-084
      IMPLICIT REAL*8 (A-H,O-Z)                                         SCHE-085
      LOGICAL LO(150)                                                   SCHE-086
      DIMENSION F(2,JMAX,*),IPI(11,*),MF(10,*),TX(*),BM(*),MC(KCB,2,*),FSCHE-087
     1A(2,KAB,*),XG(LMAX1,*),WV(22,*),FN(KMAX,*),NCT(6),FGN(KMAX,*),AM(MSCHE-088
     2T1,MT2,*),JTI(2,JIT)                                              SCHE-089
      CHARACTER*1 IP(2),IPP                                             SCHE-090
      COMMON /ANGUL/ THETA1,THETA2,DTHETA,DTHE,NCJ,NL(3),JMM(2)         SCHE-091
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            SCHE-092
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQSCHE-093
     1,ACN1,ACN2,AZ(18)                                                 SCHE-094
      COMMON /INOUT/ MR,MW,MS                                           SCHE-095
      DATA IP,PI,K2,K3,K4 /'+','-',3.1415926535897932D0,3*0/            SCHE-096
      IF (LO(120)) GO TO 60                                             SCHE-097
      IF (LO(119)) RETURN                                               SCHE-098
C TABLES OF QUANTUM NUMBERS.                                            SCHE-099
      IPD=1                                                             SCHE-100
      IF (LO(122)) IPD=2                                                SCHE-101
      IPZ=2*IPD                                                         SCHE-102
      IJA=JMIN-IPZ                                                      SCHE-103
      IAJ=IJA                                                           SCHE-104
      DO 4 JPI=1,2                                                      SCHE-105
      NC=0                                                              SCHE-106
      DO 3 I=1,NCOLS                                                    SCHE-107
      NJ1=IJA-IPI(3,I)+1                                                SCHE-108
      NJ=IPI(3,I)                                                       SCHE-109
      DO 2 J=1,NJ                                                       SCHE-110
      L1=NJ1-IPI(2,I)+1                                                 SCHE-111
      NM=IPI(2,I)                                                       SCHE-112
      DO 1 K=1,NM                                                       SCHE-113
      IF (MOD(L1+2*IPI(1,I)+2*(JPI+IPD),4).EQ.0) GO TO 1                SCHE-114
      NC=NC+1                                                           SCHE-115
      MC(NC,JPI,1)=I                                                    SCHE-116
      MC(NC,JPI,2)=L1                                                   SCHE-117
      MC(NC,JPI,3)=NJ1                                                  SCHE-118
    1 L1=L1+2                                                           SCHE-119
    2 NJ1=NJ1+2                                                         SCHE-120
    3 CONTINUE                                                          SCHE-121
    4 CONTINUE                                                          SCHE-122
      NXY=0                                                             SCHE-123
      IF (LO(81).AND.(.NOT.LO(82))) NXY=NSP(3)+2                        SCHE-124
      IF (.NOT.LO(43)) GO TO 38                                         SCHE-125
C INTERPOLATION FOR THE S-MATRIX.                                       SCHE-126
      RZ=4.D0*PI*XZ*DFLOAT(IPD)                                         SCHE-127
      I1=0                                                              SCHE-128
      I4=0                                                              SCHE-129
C LOOP ON PARITIES.                                                     SCHE-130
      DO 17 J=1,2                                                       SCHE-131
      I2=I1+NCT(J)                                                      SCHE-132
      I3=NCT(J+2)                                                       SCHE-133
      I1=I1+1                                                           SCHE-134
      IF (I1.GT.I2) GO TO 17                                            SCHE-135
C LOOP ON SOLUTIONS.                                                    SCHE-136
      DO 16 K=1,I3                                                      SCHE-137
      JI4=(3+MAX0(IPI(3,1)-1-IAJ-MC(K,J,3),IPI(2,1)-1-MC(K,J,2)-MC(K,J,3SCHE-138
     1),-2*MC(K,J,2),-2*MC(K,J,3)))/4-IPD+1                             SCHE-139
      I6=0                                                              SCHE-140
C LOOP ON S MATRIX ELEMENTS.                                            SCHE-141
      DO 15 I=I1,I2                                                     SCHE-142
      I4=I4+1                                                           SCHE-143
      I6=I6+1                                                           SCHE-144
      IV=MC(I6,J,1)                                                     SCHE-145
      JI3=(3+MAX0(IPI(3,IV)-1-IAJ-MC(I6,J,3),IPI(2,IV)-1-MC(I6,J,2)-MC(ISCHE-146
     16,J,3),-2*MC(I6,J,2),-2*MC(I6,J,3)))/4-IPD+1                      SCHE-147
      JI1=MAX0(JI3,JI4)                                                 SCHE-148
      JT1=MIN0(IPJ,IPI(10,1)-MC(K,J,2)/2,IPI(10,IV)-MC(IV,J,2)/2)       SCHE-149
      N1=0                                                              SCHE-150
      N3=1                                                              SCHE-151
    5 IF (N3.LT.JI1) GO TO 6                                            SCHE-152
      N1=N1+1                                                           SCHE-153
      K1=K2                                                             SCHE-154
      K2=K3                                                             SCHE-155
      K3=K4                                                             SCHE-156
      K4=N3                                                             SCHE-157
      IF (N1.GT.3) GO TO 9                                              SCHE-158
    6 DO 7 L=1,JIT                                                      SCHE-159
      IF (N3.LE.JTI(1,L)) GO TO 8                                       SCHE-160
    7 N3=N3+JTI(2,L)*IPD                                                SCHE-161
      L=JIT                                                             SCHE-162
    8 N3=N3+IPD                                                         SCHE-163
      IF (N3.LE.JT1) GO TO 5                                            SCHE-164
      IF (N3.GT.JTI(1,L)+JTI(2,L)*IPD+IPD) GO TO 14                     SCHE-165
      GO TO 69                                                          SCHE-166
    9 JI2=K3-IPD                                                        SCHE-167
      IF (N3.EQ.JT1) JI2=K4-IPD                                         SCHE-168
      IF (JI1.GT.JI2) GO TO 14                                          SCHE-169
      M1=(K1-K2)*(K1-K3)*(K1-K4)                                        SCHE-170
      M2=(K2-K1)*(K2-K3)*(K2-K4)                                        SCHE-171
      M3=(K3-K1)*(K3-K2)*(K3-K4)                                        SCHE-172
      M4=(K4-K1)*(K4-K2)*(K4-K3)                                        SCHE-173
      N2=1                                                              SCHE-174
      DO 13 N4=JI1,JI2,IPD                                              SCHE-175
      IF (N4.EQ.K1.OR.N4.EQ.K2.OR.N4.EQ.K3.OR.N4.EQ.K4) GO TO 13        SCHE-176
      SZ=RZ*DFLOAT(2*N4+JMIN-1)                                         SCHE-177
      IF (N2.NE.1) GO TO 10                                             SCHE-178
      N2=3                                                              SCHE-179
      A2=F(1,K1,I4)**2+F(2,K1,I4)**2                                    SCHE-180
      B2=F(1,K2,I4)**2+F(2,K2,I4)**2                                    SCHE-181
      C2=F(1,K3,I4)**2+F(2,K3,I4)**2                                    SCHE-182
      D2=F(1,K4,I4)**2+F(2,K4,I4)**2                                    SCHE-183
      IF (A2*B2*C2*D2.EQ.0.D0) GO TO 10                                 SCHE-184
      A1=DATAN2(F(2,K1,I4),F(1,K1,I4))                                  SCHE-185
      B1=DATAN2(F(2,K2,I4),F(1,K2,I4))                                  SCHE-186
      C1=DATAN2(F(2,K3,I4),F(1,K3,I4))                                  SCHE-187
      D1=DATAN2(F(2,K4,I4),F(1,K4,I4))                                  SCHE-188
      IF ((A1-B1)*(B1-C1).LT.0.D0.OR.(B1-C1)*(C1-D1).LT.0.D0) GO TO 10  SCHE-189
      N2=2                                                              SCHE-190
      A2=DLOG(A2)                                                       SCHE-191
      B2=DLOG(B2)                                                       SCHE-192
      C2=DLOG(C2)                                                       SCHE-193
      D2=DLOG(D2)                                                       SCHE-194
   10 A3=DFLOAT((N4-K2)*(N4-K3)*(N4-K4))/DFLOAT(M1)                     SCHE-195
      B3=DFLOAT((N4-K1)*(N4-K3)*(N4-K4))/DFLOAT(M2)                     SCHE-196
      C3=DFLOAT((N4-K1)*(N4-K2)*(N4-K4))/DFLOAT(M3)                     SCHE-197
      D3=DFLOAT((N4-K1)*(N4-K2)*(N4-K3))/DFLOAT(M4)                     SCHE-198
      IF (N2.EQ.3) GO TO 11                                             SCHE-199
      A4=A1*A3+B1*B3+C1*C3+D1*D3                                        SCHE-200
      B4=DEXP(0.5D0*(A2*A3+B2*B3+C2*C3+D2*D3))                          SCHE-201
      F(1,N4,I4)=B4*DCOS(A4)                                            SCHE-202
      F(2,N4,I4)=B4*DSIN(A4)                                            SCHE-203
      GO TO 12                                                          SCHE-204
   11 F(1,N4,I4)=F(1,K1,I4)*A3+F(1,K2,I4)*B3+F(1,K3,I4)*C3+F(1,K4,I4)*D3SCHE-205
      F(2,N4,I4)=F(2,K1,I4)*A3+F(2,K2,I4)*B3+F(2,K3,I4)*C3+F(2,K4,I4)*D3SCHE-206
   12 TX(IV+1)=TX(IV+1)+(F(1,N4,I4)**2+F(2,N4,I4)**2)*SZ                SCHE-207
      IF (I6.EQ.K) TX(1)=TX(1)+F(2,N4,I4)*SZ                            SCHE-208
   13 CONTINUE                                                          SCHE-209
   14 JI1=JI2+IPZ                                                       SCHE-210
      IF (JI1.LT.JT1) GO TO 6                                           SCHE-211
   15 CONTINUE                                                          SCHE-212
   16 CONTINUE                                                          SCHE-213
   17 I1=I2                                                             SCHE-214
      IF (.NOT.(LO(81).AND.LO(132))) GO TO 38                           SCHE-215
C INTERPOLATION FOR THE COMPOUND NUCLEUS.                               SCHE-216
      I1=0                                                              SCHE-217
      I4=0                                                              SCHE-218
      DO 37 J=1,3                                                       SCHE-219
      IF (J.LE.2) GO TO 18                                              SCHE-220
      IF (NXY.EQ.0) GO TO 37                                            SCHE-221
      I2=I1                                                             SCHE-222
      I3=NXY                                                            SCHE-223
      GO TO 19                                                          SCHE-224
   18 I2=I1+NCT(J+4)                                                    SCHE-225
      I3=NCT(J+2)                                                       SCHE-226
      I1=I1+1                                                           SCHE-227
   19 IF (I1.GT.I2) GO TO 37                                            SCHE-228
      I6=0                                                              SCHE-229
C LOOP ON SOLUTIONS.                                                    SCHE-230
      DO 36 K=1,I3                                                      SCHE-231
      IF (J.EQ.3) GO TO 20                                              SCHE-232
      JI4=(3+MAX0(IPI(3,1)-1-IAJ-MC(K,J,3),IPI(2,1)-1-MC(K,J,2)-MC(K,J,3SCHE-233
     1),-2*MC(K,J,2),-2*MC(K,J,3)))/4-IPD+1                             SCHE-234
   20 I6=0                                                              SCHE-235
C LOOP ON S MATRIX ELEMENTS.                                            SCHE-236
      DO 35 I=I1,I2                                                     SCHE-237
      I4=I4+1                                                           SCHE-238
      I6=I6+1                                                           SCHE-239
      JI1=1                                                             SCHE-240
      IF (J.EQ.3) GO TO 21                                              SCHE-241
      IV=MC(I6,J,1)                                                     SCHE-242
      IVL=IV+NCOLL+1                                                    SCHE-243
      JI3=(3+MAX0(IPI(3,IV)-1-IAJ-MC(I6,J,3),IPI(2,IV)-1-MC(I6,J,2)-MC(ISCHE-244
     16,J,3),-2*MC(I6,J,2),-2*MC(I6,J,3)))/4-IPD+1                      SCHE-245
      JI1=MAX0(JI3,JI4)                                                 SCHE-246
      GO TO 22                                                          SCHE-247
   21 IVL=K+NCOLL+NCOLS+1                                               SCHE-248
   22 NT=0                                                              SCHE-249
      JPJ=0                                                             SCHE-250
      DO 23 L=JI1,IPK                                                   SCHE-251
      IF (DABS(FN(L,I4)).LT.1.D-15) GO TO 23                            SCHE-252
      NT=NT+1                                                           SCHE-253
      JPJ=L                                                             SCHE-254
   23 CONTINUE                                                          SCHE-255
      NT=MIN0(NT,4)-1                                                   SCHE-256
      IF (NT.LT.1) GO TO 35                                             SCHE-257
      N1=0                                                              SCHE-258
      N3=1                                                              SCHE-259
   24 IF (N3.LT.JI1) GO TO 25                                           SCHE-260
      N1=N1+1                                                           SCHE-261
      K1=K2                                                             SCHE-262
      K2=K3                                                             SCHE-263
      K3=K4                                                             SCHE-264
      K4=N3                                                             SCHE-265
      IF (N1.GT.NT) GO TO 28                                            SCHE-266
   25 DO 26 L=1,JIT                                                     SCHE-267
      IF (N3.LE.JTI(1,L)) GO TO 27                                      SCHE-268
   26 N3=N3+JTI(2,L)*IPD                                                SCHE-269
      L=JIT                                                             SCHE-270
   27 N3=N3+IPD                                                         SCHE-271
      IF (N3.LE.JPJ) GO TO 24                                           SCHE-272
      IF (N3.GT.JTI(1,L)+JTI(2,L)*IPD+IPD) GO TO 34                     SCHE-273
   28 JI2=K3-IPD                                                        SCHE-274
      IF (N3.EQ.JPJ) JI2=K4-IPD                                         SCHE-275
      IF (JI1.GT.JI2) GO TO 34                                          SCHE-276
      DO 33 N4=JI1,JI2,IPD                                              SCHE-277
      IF (N4.EQ.K1.OR.N4.EQ.K2.OR.N4.EQ.K3.OR.N4.EQ.K4) GO TO 33        SCHE-278
      SZ=RZ*DFLOAT(2*N4+JMIN-1)                                         SCHE-279
      IF (NT.LT.3) GO TO 29                                             SCHE-280
      A3=DFLOAT((N4-K2)*(N4-K3)*(N4-K4))/DFLOAT((K1-K2)*(K1-K3)*(K1-K4))SCHE-281
      B3=DFLOAT((N4-K1)*(N4-K3)*(N4-K4))/DFLOAT((K2-K1)*(K2-K3)*(K2-K4))SCHE-282
      C3=DFLOAT((N4-K1)*(N4-K2)*(N4-K4))/DFLOAT((K3-K1)*(K3-K2)*(K3-K4))SCHE-283
      D3=DFLOAT((N4-K1)*(N4-K2)*(N4-K3))/DFLOAT((K4-K1)*(K4-K2)*(K4-K3))SCHE-284
      FN(N4,I4)=DABS(FN(K1,I4)*A3+FN(K2,I4)*B3+FN(K3,I4)*C3+FN(K4,I4)*D3SCHE-285
     1)                                                                 SCHE-286
      GO TO 31                                                          SCHE-287
   29 IF (NT.EQ.1) GO TO 30                                             SCHE-288
      B3=DFLOAT((N4-K3)*(N4-K4))/DFLOAT((K2-K3)*(K2-K4))                SCHE-289
      C3=DFLOAT((N4-K2)*(N4-K4))/DFLOAT((K3-K2)*(K3-K4))                SCHE-290
      D3=DFLOAT((N4-K2)*(N4-K3))/DFLOAT((K4-K2)*(K4-K3))                SCHE-291
      FN(N4,I4)=DABS(FN(K2,I4)*B3+FN(K3,I4)*C3+FN(K4,I4)*D3)            SCHE-292
      GO TO 31                                                          SCHE-293
   30 C3=DFLOAT(N4-K4)/DFLOAT(K3-K4)                                    SCHE-294
      D3=DFLOAT(N4-K3)/DFLOAT(K4-K3)                                    SCHE-295
      FN(N4,I4)=DABS(FN(K3,I4)*C3+FN(K4,I4)*D3)                         SCHE-296
   31 TX(IVL)=TX(IVL)+FN(N4,I4)*SZ                                      SCHE-297
      IF ((J.EQ.3).OR.(FN(N4,I4).EQ.0.D0)) GO TO 33                     SCHE-298
      CALL COCN(MC(K,J,2)+IPZ*N4,MC(K,J,2)+IPZ*N4,MC(K,J,3)+IPZ*N4,MC(K,SCHE-299
     1J,3)+IPZ*N4,IPI(3,1)-1,IPI(2,1)-1,IJA+IPZ*N4,N4,BM,BM(N4+1),ID1-N4SCHE-300
     2)                                                                 SCHE-301
      CALL COCN(MC(I6,J,2)+IPZ*N4,MC(I6,J,2)+IPZ*N4,MC(I6,J,3)+IPZ*N4,MCSCHE-302
     1(I6,J,3)+IPZ*N4,IPI(3,IV)-1,IPI(2,IV)-1,IJA+IPZ*N4,N4,BM(N4+1),BM(SCHE-303
     22*N4+1),ID1-2*N4)                                                 SCHE-304
      IF ((IV.EQ.1).AND.(IPI(3,1).GT.1).AND.(K.NE.I6)) CALL COCN(MC(K,J,SCHE-305
     12)+IPZ*N4,MC(I6,J,2)+IPZ*N4,MC(K,J,3)+IPZ*N4,MC(I6,J,3)+IPZ*N4,IPISCHE-306
     2(3,IV)-1,IPI(2,IV)-1,IJA+IPZ*N4,N4,BM(2*N4+1),BM(3*N4+1),ID1-3*N4)SCHE-307
      DO 32 LL=1,N4                                                     SCHE-308
      Z=BM(LL)*BM(LL+N4)                                                SCHE-309
      IF ((IV.EQ.1).AND.(IPI(3,1).GT.1).AND.(K.NE.I6)) Z=BM(LL+2*N4)**2+SCHE-310
     1Z                                                                 SCHE-311
   32 FGN(LL,IV)=FGN(LL,IV)+Z*FN(N4,I4)*XZ                              SCHE-312
   33 CONTINUE                                                          SCHE-313
   34 JI1=JI2+IPZ                                                       SCHE-314
      IF (JI1.LT.JPJ) GO TO 25                                          SCHE-315
   35 CONTINUE                                                          SCHE-316
   36 CONTINUE                                                          SCHE-317
   37 I1=I2                                                             SCHE-318
   38 IK=1                                                              SCHE-319
      NPT=0                                                             SCHE-320
      DO 57 IJ=1,IPJ,IPD                                                SCHE-321
C TRANSFER OF S MATRIX FOR A GIVEN ANGULAR MOMENTUM.                    SCHE-322
      I1=0                                                              SCHE-323
      I4=0                                                              SCHE-324
      DO 41 J=1,2                                                       SCHE-325
      I2=I1+NCT(J+2)                                                    SCHE-326
      I1=I1+1                                                           SCHE-327
      IF (I1.GT.I2) GO TO 41                                            SCHE-328
      I3=NCT(J)                                                         SCHE-329
      DO 40 I=I1,I2                                                     SCHE-330
      DO 39 K=1,I3                                                      SCHE-331
      I4=I4+1                                                           SCHE-332
      FA(1,K,I)=F(1,IJ,I4)                                              SCHE-333
   39 FA(2,K,I)=F(2,IJ,I4)                                              SCHE-334
   40 CONTINUE                                                          SCHE-335
   41 I1=I2                                                             SCHE-336
      DO 42 I=1,KBA                                                     SCHE-337
      F(1,IJ,I)=0.D0                                                    SCHE-338
   42 F(2,IJ,I)=0.D0                                                    SCHE-339
      J1=0                                                              SCHE-340
      IAJ=IAJ+IPZ                                                       SCHE-341
      AAJ=DFLOAT(IAJ+1)                                                 SCHE-342
      DO 54 JI=1,2                                                      SCHE-343
      J2=J1+NCT(JI+2)                                                   SCHE-344
      J1=J1+1                                                           SCHE-345
      IF (J1.GT.J2) GO TO 54                                            SCHE-346
      NC1=0                                                             SCHE-347
      NC2=0                                                             SCHE-348
      NC=NCT(JI)                                                        SCHE-349
C GEOMETRIC COEFFICIENT FOR THE TRANSFORMATION TO HELICITY COEFFICIENTS.SCHE-350
      DO 48 IC=1,NC                                                     SCHE-351
      MC(IC,JI,2)=MC(IC,JI,2)+IPZ                                       SCHE-352
      MC(IC,JI,3)=MC(IC,JI,3)+IPZ                                       SCHE-353
      IV=MC(IC,JI,1)                                                    SCHE-354
      NI=IPI(2,IV)                                                      SCHE-355
      MI=IPI(3,IV)                                                      SCHE-356
      MC(IC,1,4)=MIN0(MC(IC,JI,3)-IABS(IAJ+1-MI),MC(IC,JI,2)-IABS(MC(IC,SCHE-357
     1JI,3)+1-NI))                                                      SCHE-358
      IF (MC(IC,JI,2).GT.2*IPI(10,IV)) MC(IC,1,4)=-1                    SCHE-359
      IF (MC(IC,1,4).LT.0) GO TO 48                                     SCHE-360
      IF (IV.EQ.1) NC1=NC1+1                                            SCHE-361
      NC2=NC2+1                                                         SCHE-362
      A1=0.D0                                                           SCHE-363
      YM=DFLOAT(MC(IC,JI,3)-NI+1)                                       SCHE-364
      XB1=0.5D0*DFLOAT(NI**2+(MC(IC,JI,3)-MC(IC,JI,2))*(MC(IC,JI,2)+MC(ISCHE-365
     1C,JI,3)+2)-1)                                                     SCHE-366
      C3=0.D0                                                           SCHE-367
C COUPLING FOR PARTICLE HELICITY.                                       SCHE-368
      DO 43 I1=1,NI                                                     SCHE-369
      BM(I1)=0.D0                                                       SCHE-370
      IF (IABS(2*I1-NI-1).GT.MC(IC,JI,3)) GO TO 43                      SCHE-371
      N3=(MC(IC,JI,3)+2*I1-NI-1)*(MC(IC,JI,3)-2*I1+NI+3)*(I1-1)         SCHE-372
      IF (N3.EQ.0) BM(I1)=DFLOAT(2*MOD(I1,2)-1)                         SCHE-373
      IF (N3.LE.0) GO TO 43                                             SCHE-374
      C2=C3                                                             SCHE-375
      C3=DSQRT(DFLOAT((NI-I1+1)*(I1-1))*(YM+DFLOAT(2*I1-2))*(YM+DFLOAT(2SCHE-376
     1*NI-2*I1+2)))                                                     SCHE-377
      BM(I1)=(XB1-DFLOAT(2*I1-NI-3)**2)*BM(I1-1)/C3                     SCHE-378
      IF (I1.GE.3) BM(I1)=BM(I1)-C2*BM(I1-2)/C3                         SCHE-379
   43 A1=A1+BM(I1)**2                                                   SCHE-380
      DO 47 I1=1,NI                                                     SCHE-381
      A2=0.D0                                                           SCHE-382
      IF (DABS(BM(I1)).LT.1.D-10) GO TO 45                              SCHE-383
      IA=2*I1-NI-1                                                      SCHE-384
      XB1=0.5D0*DFLOAT(MI**2+(IAJ-MC(IC,JI,3))*(IAJ+MC(IC,JI,3)+2)-1)   SCHE-385
      C3=0.D0                                                           SCHE-386
C COUPLING FOR TARGET HELICITY.                                         SCHE-387
      DO 44 I2=1,MI                                                     SCHE-388
      BM(NI+I2)=0.D0                                                    SCHE-389
      IF (IABS(2*I2-MI-1-IA).GT.IAJ) GO TO 44                           SCHE-390
      N3=(IAJ+IA-2*I2+MI+3)*(IAJ-IA+2*I2-MI-1)*(I2-1)                   SCHE-391
      IF (N3.EQ.0) BM(NI+I2)=DFLOAT(2*MOD(I2,2)-1)                      SCHE-392
      IF (N3.LE.0) GO TO 44                                             SCHE-393
      C2=C3                                                             SCHE-394
      IB=2*I2-MI-3                                                      SCHE-395
      C3=DSQRT(DFLOAT((MI-I2+1)*(I2-1))*(AAJ**2-DFLOAT(2*I1-2*I2+MI-NI+1SCHE-396
     1)**2))                                                            SCHE-397
      BM(NI+I2)=((XB1-DFLOAT(IB*(IB-IA)))*BM(NI+I2-1)-C2*BM(NI+I2-2))/C3SCHE-398
   44 A2=A2+BM(NI+I2)**2                                                SCHE-399
      IF (A2*A1.NE.0.D0) A2=DFLOAT(MOD(1+MC(IC,JI,3)-IAJ+MI,4)-1)*DSQRT(SCHE-400
     1AAJ/(A1*A2))                                                      SCHE-401
   45 DO 46 I2=1,MI                                                     SCHE-402
   46 AM(I1,I2,IC)=BM(I1)*BM(NI+I2)*A2                                  SCHE-403
   47 CONTINUE                                                          SCHE-404
   48 CONTINUE                                                          SCHE-405
      IF (NC1.EQ.0) GO TO 54                                            SCHE-406
      NCIN=NCT(JI+2)                                                    SCHE-407
C TRANSFORMATION.                                                       SCHE-408
      IF (IK.NE.IJ) GO TO 49                                            SCHE-409
      IF (.NOT.(LO(56).OR.LO(60))) GO TO 49                             SCHE-410
      BJ=.5D0*DFLOAT(IAJ)                                               SCHE-411
      JIJ=1+MOD(IJ+JI,2)                                                SCHE-412
      NC1=NC1*NC2                                                       SCHE-413
      IF (LO(56)) WRITE (MW,1000) BJ,IP(JIJ)                            SCHE-414
      NPT=NPT+1                                                         SCHE-415
      IF (LO(60)) WRITE (99,1001) BJ,IP(JIJ),NC2,NC1                    SCHE-416
   49 NC1=0                                                             SCHE-417
      DO 53 IC=1,NCIN                                                   SCHE-418
      IF (MC(IC,1,4).LT.0) GO TO 53                                     SCHE-419
      ICX=IC                                                            SCHE-420
      NC1=NC1+1                                                         SCHE-421
      IF (JI.EQ.2) ICX=ICX+NCT(3)                                       SCHE-422
      LCI=MC(IC,JI,2)/2+1                                               SCHE-423
      NC2=0                                                             SCHE-424
      DO 52 ICP=1,NC                                                    SCHE-425
      IF (MC(ICP,1,4).LT.0) GO TO 52                                    SCHE-426
      LCP=MC(ICP,JI,2)/2+1                                              SCHE-427
      IV=MC(ICP,JI,1)                                                   SCHE-428
      NC2=NC2+1                                                         SCHE-429
      C1=XG(LCI,1)+XG(LCP,IV)                                           SCHE-430
      IF (IK.NE.IJ) GO TO 50                                            SCHE-431
      IF (.NOT.(LO(56).OR.LO(60))) GO TO 50                             SCHE-432
      B1=-2.D0*FA(2,ICP,ICX)                                            SCHE-433
      IF (IC.EQ.ICP) B1=B1+1.D0                                         SCHE-434
      B2=2.D0*FA(1,ICP,ICX)                                             SCHE-435
      B3=DSQRT(B1**2+B2**2)                                             SCHE-436
      D1=0.D0                                                           SCHE-437
      IF (B3.NE.0.D0) D1=DATAN2(B2,B1)                                  SCHE-438
      D2=DMOD(D1+C1+PI,2.D0*PI)-PI                                      SCHE-439
      LC=LCP-1                                                          SCHE-440
      BJ=0.5D0*DFLOAT(MC(ICP,JI,3))                                     SCHE-441
      IF (LO(56)) WRITE (MW,1002) NC1,NC2,IV,LC,BJ,B1,B2,B3,D1,D2       SCHE-442
      IF (LO(60)) WRITE (99,1003) NC1,NC2,IV,LC,BJ,B1,B2,B3             SCHE-443
C MULTIPLICATION BY THE COULOMB PHASE.                                  SCHE-444
   50 A1=DCOS(C1)                                                       SCHE-445
      A2=DSIN(C1)                                                       SCHE-446
      C2=FA(1,ICP,ICX)*A1-FA(2,ICP,ICX)*A2                              SCHE-447
      C3=FA(1,ICP,ICX)*A2+FA(2,ICP,ICX)*A1                              SCHE-448
      I1=IPI(6,IV)                                                      SCHE-449
      I2=IPI(7,IV)                                                      SCHE-450
C  HELICITY SCATTERING COEFFICIENTS.                                    SCHE-451
      DO 51 ID=I1,I2                                                    SCHE-452
      MF1=MF(1,ID)                                                      SCHE-453
      MF2=MF(2,ID)                                                      SCHE-454
      MF3=MF(3,ID)                                                      SCHE-455
      MF4=MF(4,ID)                                                      SCHE-456
      C1=AM(MF1,MF2,ICP)*AM(MF3,MF4,IC)                                 SCHE-457
      F(1,IJ,ID)=F(1,IJ,ID)+C2*C1                                       SCHE-458
   51 F(2,IJ,ID)=F(2,IJ,ID)+C3*C1                                       SCHE-459
   52 CONTINUE                                                          SCHE-460
   53 CONTINUE                                                          SCHE-461
   54 CONTINUE                                                          SCHE-462
      IF ((.NOT.LO(43)).OR.IK.NE.IJ) GO TO 56                           SCHE-463
      DO 55 L=1,JIT                                                     SCHE-464
      IF (IK.LE.JTI(1,L)) GO TO 56                                      SCHE-465
   55 IK=IK+JTI(2,L)*IPD                                                SCHE-466
   56 IK=IK+IPD                                                         SCHE-467
   57 CONTINUE                                                          SCHE-468
      IF (.NOT.LO(60)) GO TO 60                                         SCHE-469
      WRITE (60,1004) WV(1,1),WV(13,1),WV(2,1),IPI(4,1),NPT             SCHE-470
      REWIND 99                                                         SCHE-471
      DO 59 I=1,NPT                                                     SCHE-472
      READ (99,1001) U1,IPP,K1,K2                                       SCHE-473
      WRITE (60,1001) U1,IPP,K1,K2                                      SCHE-474
      DO 58 K=1,K2                                                      SCHE-475
      READ (99,1003) K1,K2,K3,K4,BJ,B1,B2,B3                            SCHE-476
   58 WRITE (60,1005) K1,K2,K3,K4,BJ,B1,B2,B3                           SCHE-477
   59 CONTINUE                                                          SCHE-478
      CLOSE (99,STATUS='DELETE')                                        SCHE-479
   60 IF (.NOT.LO(65)) GO TO 61                                         SCHE-480
      NSA=7*(NLT+2*IPJ+1)+1                                             SCHE-481
      IF (NSA.GT.ID1) CALL MEMO('SCHE',ID1,NSA)                         SCHE-482
      CALL LCSP(F,FGN,JMAX,KMAX,IPI,NCOLL,NCOLS,MF,WV,JMIN,IPJ,IPJ,BM,BMSCHE-483
     1(NSA),ID1-NSA,LO)                                                 SCHE-484
   61 IF (LO(120).OR.(.NOT.LO(41))) RETURN                              SCHE-485
C ELIMINATION OF FACTORS 1/(1-X*COS(THETA)).                            SCHE-486
      NCJ=MIN0(NCJ,IPJ-1)                                               SCHE-487
      IF (NCJ.EQ.0) RETURN                                              SCHE-488
      DO 68 NJ=1,NCJ                                                    SCHE-489
      IPJ=IPJ-1                                                         SCHE-490
C LOOP ON THE INDEPENDENT AMPLITUDES.                                   SCHE-491
      DO 67 K=1,KBA                                                     SCHE-492
      IF (MF(6,K).EQ.99999) GO TO 62                                    SCHE-493
      M1=MF(5,K)                                                        SCHE-494
      M2=MF(6,K)                                                        SCHE-495
      M3=((IABS(M1+M2)+IABS(M1-M2))/2-JMIN)/2+1                         SCHE-496
   62 B3=0.25D0*DFLOAT(M1*M2)                                           SCHE-497
      IF (M3.GT.IPJ) GO TO 70                                           SCHE-498
      D1=0.D0                                                           SCHE-499
      D2=0.D0                                                           SCHE-500
      MJ=2*M3+JMIN-2                                                    SCHE-501
      C3=0.5D0*DFLOAT(MJ)                                               SCHE-502
      F1=DFLOAT((MJ+M2)/2+1-M3)                                         SCHE-503
      F2=DFLOAT((MJ-M2)/2+1-M3)                                         SCHE-504
      F3=DFLOAT((MJ+M1)/2+1-M3)                                         SCHE-505
      F4=DFLOAT((MJ-M1)/2+1-M3)                                         SCHE-506
      C1=0.D0                                                           SCHE-507
      B1=0.D0                                                           SCHE-508
      A1=0.D0                                                           SCHE-509
      B4=0.D0                                                           SCHE-510
C CALC. OF X WHICH MINIMISES THE DIFFERENCES WITH A WEIGHT (J+1)**2     SCHE-511
C FOR THE 5 LAST ONE.                                                   SCHE-512
      NPJ=IPJ-5                                                         SCHE-513
      DO 63 I=M3,IPJ                                                    SCHE-514
      FI=I                                                              SCHE-515
      A2=A1                                                             SCHE-516
      B2=B1                                                             SCHE-517
      A1=F(1,I,K)                                                       SCHE-518
      B1=F(2,I,K)                                                       SCHE-519
      C3=C3+1.D0                                                        SCHE-520
      C2=C1                                                             SCHE-521
      C1=DSQRT((F1+FI)*(F2+FI)*(F3+FI)*(F4+FI))/(C3*(2.D0*C3+1.D0))     SCHE-522
      IF (B3.NE.0.D0) B4=B3/(C3*C3-C3)                                  SCHE-523
      BM(2*I-1)=A1*B4+A2*C2+F(1,I+1,K)*C1                               SCHE-524
      BM(2*I)=B1*B4+B2*C2+F(2,I+1,K)*C1                                 SCHE-525
      IF (I.LE.NPJ) GO TO 63                                            SCHE-526
      D1=D1+C3*C3*(F(1,I,K)*BM(2*I-1)+F(2,I,K)*BM(2*I))                 SCHE-527
      D2=D2+C3*C3*(BM(2*I-1)**2+BM(2*I)**2)                             SCHE-528
   63 C1=C1*(C3+C3+1.D0)/(C3+C3-1.D0)                                   SCHE-529
      A3=.9999999D0                                                     SCHE-530
      IF (LO(18)) GO TO 64                                              SCHE-531
      IF (D2.NE.0.D0) A3=D1/D2                                          SCHE-532
      A4=A3                                                             SCHE-533
C X IS FIXED BETWEEN +1 AND -1                                          SCHE-534
      IF (A4.GT..9999999D0) A4=.9999999D0                               SCHE-535
      IF (A4.LT.-.9999999D0) A4=-.9999999D0                             SCHE-536
C CALCULATION OF THE NEW SCATTERING COEFFICIENTS.                       SCHE-537
      GO TO 65                                                          SCHE-538
   64 A4=A3*DFLOAT(2*MOD(NJ,2)-1)                                       SCHE-539
   65 DO 66 I=M3,IPJ                                                    SCHE-540
      F(1,I,K)=F(1,I,K)-A4*BM(2*I-1)                                    SCHE-541
   66 F(2,I,K)=F(2,I,K)-A4*BM(2*I)                                      SCHE-542
      IF (.NOT.LO(116)) WRITE (MW,1006) K,A4,A3,F(1,IPJ,K),F(2,IPJ,K),F(SCHE-543
     11,IPJ+1,K),F(2,IPJ+1,K)                                           SCHE-544
      F(1,IPJ+1,K)=A4                                                   SCHE-545
   67 CONTINUE                                                          SCHE-546
   68 CONTINUE                                                          SCHE-547
      RETURN                                                            SCHE-548
   69 WRITE (MW,1007) N1                                                SCHE-549
      GO TO 71                                                          SCHE-550
   70 WRITE (MW,1008) NJ                                                SCHE-551
   71 WRITE (MW,1009)                                                   SCHE-552
      STOP                                                              SCHE-553
 1000 FORMAT (//' CHANNEL SPIN AND PARITY =',F7.1,A1//'  IC ICP N    L  SCHE-554
     1  J',18X,'S MATRIX',20X,'|S|',7X,'PHASE /WITH COUL.')             SCHE-555
 1001 FORMAT (1X,F9.1,4X,A1,1X,I4,1X,I4)                                SCHE-556
 1002 FORMAT (1X,3I3,I5,F7.1,4X,1P,2D15.7,' I',4X,0P,3F11.8)            SCHE-557
 1003 FORMAT (1X,3(I2,1X),I3,1X,F5.1,1X,2(1P,D15.7,0P,1X),5X,F11.8)     SCHE-558
 1004 FORMAT ('<S-MATRIX>',F10.2,1P,D20.8,0P,F10.2,2I5)                 SCHE-559
 1005 FORMAT (1X,3(I2,1X),I3,1X,F5.1,1X,2(1P,D15.7,0P,1X),'I',4X,F11.8) SCHE-560
 1006 FORMAT (' AMPLITUDE =',I3,D15.7,' (',D15.7,')  NEW',2D15.7,3X,'OLDSCHE-561
     1',2D15.7)                                                         SCHE-562
 1007 FORMAT (5X,I2,' AMPLITUDES INSUFFICIENT TO INTERPOLATE.')         SCHE-563
 1008 FORMAT (' AMPLITUDES INSUFFICIENT TO FACTORISE (1-COS)',I2,' TIMESSCHE-564
     1')                                                                SCHE-565
 1009 FORMAT (//' IN SCHE  ...  STOP  ...')                             SCHE-566
      END                                                               SCHE-567
