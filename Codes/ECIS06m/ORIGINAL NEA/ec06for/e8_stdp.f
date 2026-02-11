C 07/03/07                                                      ECIS06  STDP-000
      SUBROUTINE STDP(V,IVY,ISM,VAL,NVAL,NX,IDX,X,WV,PGN,XGN,NPP,IZZ,P,LSTDP-001
     1O)                                                                STDP-002
C COMPUTES FORM FACTORS INDEPENDENTLY OF MODELS.                        STDP-003
C INPUT:     IVY:     TABLE OF FORM FACTORS (SEE REDM).                 STDP-004
C            ISM:     NUMBER OF INTEGRATION STEPS.                      STDP-005
C            VAL,NVAL:FOR OPTICAL MODEL PARAMETERS TO USE HERE.         STDP-006
C            NX:      LENGTH OF WORKING SPACE.                          STDP-007
C            WV:      STEP SIZE IN WV(8,*).                             STDP-008
C            PGN:     WEIGHTS OF GAUSS-LEGENDRE INTEGRATION             STDP-009
C            XGN:   : ABSCISSAE OF GAUSS-LEGENDRE INTEGRATION.          STDP-010
C            NPP:     NUMBER OF OPTICAL POTENTIALS.                     STDP-011
C            ACONV:   VALUE BELOW WHICH THE FOLDING IS NEGLECTED.       STDP-012
C            LO(I):   LOGICAL CONTROLS:                                 STDP-013
C               LO(6)  =.TRUE. USE DEFORMATION LENGTHS.                 STDP-014
C               LO(9)  =.TRUE. SYMMETRISED WOODS-SAXON FORM FACTORS WHENSTDP-015
C                              THE RADIUS IS NEGATIVE.                  STDP-016
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            STDP-017
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. STDP-018
C               LO(100)=.TRUE. DIRAC EQUATION.                          STDP-019
C OUTPUT:    V:       ELASTIC AND INELASTIC FORM FACTORS IN THE SEQUENCESTDP-020
C                     CENTRAL-REAL, SPIN-ORBIT-REAL, TRANSITION REAL,   STDP-021
C                     TRANSITION SPIN-ORBIT-REAL, IMAGINARY POTENTIALS  STDP-022
C                     FOLLOWED BY COULOMB, COULOMB TRANSITION POTENTIALSSTDP-023
C            IDX:     LENGTH OF WORKING FIELD USED.                     STDP-024
C WORKING AREAS:                                                        STDP-025
C            X,P:     TO COMPUTE BOUND FUNCTIONS, FOLD COULOMB POTENTIALSTDP-026
C                     AND COMPUTE ROTATIONAL FORM FACTORS AND BESSEL    STDP-027
C                     FUNCTIONS, IN EQUIVALENCE BY CALL.                STDP-028
C            IZZ:     FOR FOLDING.                                      STDP-029
C                                                                       STDP-030
C FOR THE COMMON  /DCONS/ SEE CALC.                                     STDP-031
C FOR THE COMMON  /POTE1/ SEE REDM.                                     STDP-032
C                                                                       STDP-033
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     STDP-034
C  CCZ:       COULOMB ALPHA CONSTANT.                                   STDP-035
C   USED:     CCZ.                                                      STDP-036
C                                                                       STDP-037
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE1/:                     STDP-038
C  ITX(16):   STARTING ADDRESS OF DIFFERENT FORM FACTORS (SEE REDM).    STDP-039
C  IMAX:      MAXIMUM ANGULAR MOMENTUM.                                 STDP-040
C  INTC:      NUMBER OF FORM FACTORS WITHOUT DEFORMED SPIN-ORBIT        STDP-041
C             INCLUDING CORRECTION TERMS.                               STDP-042
C  INLS:      NUMBER OF SPIN-ORBIT FORM FACTORS NOT TAKING INTO ACCOUNT STDP-043
C             MULTIPLICATION BY 2.                                      STDP-044
C  INVD:      IDEM FOR COULOMB SPIN-ORBIT.                              STDP-045
C  ITXM:      TOTAL NUMBER OF FORM FACTORS.                             STDP-046
C   USED:     ITX,IMAX,INTC,INLS,INVD,ITXM.                             STDP-047
C                                                                       STDP-048
C***********************************************************************STDP-049
      IMPLICIT REAL*8 (A-H,O-Z)                                         STDP-050
      LOGICAL LO(150),LT(9)                                             STDP-051
      DIMENSION V(ISM,*),IVY(7,*),VAL(*),NVAL(2,*),X(3,*),WV(22,*),PGN(1STDP-052
     10),XGN(10),IZZ(4,*),P(*),NIJ(3),ITZ(10),ZB(77),ZB1(40),ZB2(37),Y(3STDP-053
     2),VR(5),CL(8)                                                     STDP-054
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       STDP-055
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            STDP-056
      COMMON /INOUT/ MR,MW,MS                                           STDP-057
      COMMON /POTE1/ ITX(16),IMAX,INTC,INLS,INVC,INVD,ITXM              STDP-058
      EQUIVALENCE (ZB1,ZB),(ZB2,ZB(41))                                 STDP-059
      DATA NIJ /1,24,11/                                                STDP-060
      DATA PIS2 /1.5707963267949D0/                                     STDP-061
      DATA ITZ /3,3,3,3,4,5,3,3,5,5/                                    STDP-062
      DATA ZB1 /4.49340945790906D0,7.72525183693771D0,5.76345919689455D0STDP-063
     1,9.09501133047635D0,12.3229409705666D0,6.98793200050052D0,10.41711STDP-064
     285473794D0,13.6980231532492D0,16.9236212852138D0,8.18256145257124DSTDP-065
     30,11.7049071545704D0,15.0396647076165D0,18.3012559595420D0,21.5254STDP-066
     4177333999D0,9.35581211104275D0,12.9665301727743D0,16.3547096393505STDP-067
     5D0,19.6531521018212D0,22.9045506479037D0,26.1277501372255D0,10.512STDP-068
     68354080940D0,14.2073924588425D0,17.6479748701659D0,20.983463068944STDP-069
     78D0,24.2627680423970D0,27.5078683649043D0,30.7303807316466D0,11.65STDP-070
     870321925164D0,15.4312892102684D0,18.9229991985461D0,22.29534801913STDP-071
     908D0,25.6028559538106D0,28.8703733470427D0,32.1111962396826D0,35.3STDP-072
     A331941827165D0,12.7907817119721D0,16.6410028815122D0,20.1824707649STDP-073
     B492D0,23.5912748179830D0,26.9270407788180D0/                      STDP-074
      DATA ZB2 /30.2172627093614D0,33.4768008195015D0,36.7145291272447D0STDP-075
     1,39.9361278108677D0,13.9158226105049D0,17.8386431992053D0,21.42848STDP-076
     269721154D0,24.8732139238751D0,28.2371343599681D0,31.5501883818318DSTDP-077
     30,34.8286965376857D0,38.0824790873276D0,41.3178646902445D0,44.5391STDP-078
     4446334095D0,15.0334693037434D0,19.0258535361278D0,22.6627206581361STDP-079
     5D0,26.1427676433791D0,29.5346341078439D0,32.8705345976875D0,36.168STDP-080
     61571359112D0,39.4382144800081D0,42.6876512846611D0,45.921201763835STDP-081
     76D0,49.1422214247461D0,16.1447429423013D0,20.2039426328117D0,23.88STDP-082
     865307559684D0,27.4012592588663D0,30.8207940864510D0,34.17947466648STDP-083
     932D0,37.4962736357858D0,40.7827470981251D0,44.0464252109438D0,47.2STDP-084
     A924656052694D0,50.5245397255712D0,53.7453428657930D0/             STDP-085
C MEANING OF INTERNAL LOGICAL LT:                                       STDP-086
C LT(1) SPIN-ORBIT TRANSITION (SCHROEDINGER),                           STDP-087
C LT(2) SPIN-ORBIT POTENTIAL DERIVATIVE,                                STDP-088
C LT(3) FORM FACTOR GIVEN BY POINTS OR BOUND STATE,                     STDP-089
C LT(4) COMPUTATION OF A COULOMB POTENTIAL FROM THE CHARGE DENSITY,     STDP-090
C LT(5) ZERO DIFFUSENESS COULOMB,                                       STDP-091
C LT(6) SURFACE POTENTIAL,                                              STDP-092
C LT(7) FIRST OR SECOND PASSAGE IN WOODS-SAXON/BESSEL LOOP,             STDP-093
C LT(8) WOODS-SAXON FORM FACTOR,                                        STDP-094
C LT(9) NOT SYMMETRIC.                                                  STDP-095
      IDX=0                                                             STDP-096
      HM=1000.D0                                                        STDP-097
      NFO=0                                                             STDP-098
      N=ITXM                                                            STDP-099
      IF (LO(100)) N=N-ITX(7)                                           STDP-100
      DO 2 I=1,N                                                        STDP-101
      DO 1 J=1,4                                                        STDP-102
    1 IZZ(J,I)=0                                                        STDP-103
    2 CONTINUE                                                          STDP-104
      IF (INTC.EQ.0) GO TO 8                                            STDP-105
      NTT=24*NPP                                                        STDP-106
      DO 7 I=1,INTC                                                     STDP-107
      IF (LO(100)) GO TO 5                                              STDP-108
      DO 3 J=9,12                                                       STDP-109
      IF ((.NOT.LO(12)).AND.(MOD(J,2).EQ.0)) GO TO 3                    STDP-110
      K=ITX(J)+I                                                        STDP-111
      IZZ(3,K)=IVY(7,I)                                                 STDP-112
    3 CONTINUE                                                          STDP-113
      K=IVY(3,I)                                                        STDP-114
      IF (K.EQ.0) GO TO 4                                               STDP-115
      IZZ(3,K+ITX(13))=IVY(7,I)                                         STDP-116
      IF (LO(14)) IZZ(3,K+ITX(10))=IVY(7,I)                             STDP-117
    4 K=IVY(4,I)                                                        STDP-118
      IF (K.NE.0) IZZ(3,K+ITX(15))=IVY(7,I)                             STDP-119
      K=IVY(5,I)                                                        STDP-120
      IF (K.NE.0) IZZ(3,K+ITX(16))=IVY(7,I)                             STDP-121
      GO TO 7                                                           STDP-122
    5 DO 6 J=1,11                                                       STDP-123
    6 IZZ(3,NTT+J)=IVY(7,I)                                             STDP-124
    7 NTT=NTT+11                                                        STDP-125
    8 NMA=NVAL(1,1)                                                     STDP-126
    9 IF (NMA.GE.NVAL(1,2)) GO TO 97                                    STDP-127
      I1=NVAL(1,NMA)                                                    STDP-128
      IV=NVAL(2,NMA)                                                    STDP-129
      ITV=MOD(I1-1,8)+1                                                 STDP-130
      J1=(I1-1)/8                                                       STDP-131
      ITT=ITV                                                           STDP-132
      LT(1)=(ITV-5)*(ITV-6)*(ITV-8).EQ.0                                STDP-133
      INL=0                                                             STDP-134
      IF (J1.GT.NPP) ITT=ITV+8                                          STDP-135
      IF (LO(100)) GO TO 10                                             STDP-136
      JI=1                                                              STDP-137
      IJ=1                                                              STDP-138
      IF (J1.GT.NPP) J1=J1-NPP                                          STDP-139
      L1=J1+ITX(ITT)                                                    STDP-140
      IF ((ITT.EQ.13).OR.(ITT.EQ.14)) INL=INLS                          STDP-141
      IF (ITT.EQ.16) INL=INVD                                           STDP-142
      IF (INL.NE.0) IJ=2                                                STDP-143
      GO TO 13                                                          STDP-144
   10 IF (J1.LE.NPP) GO TO 11                                           STDP-145
      JI=3                                                              STDP-146
      IF (LT(1)) INL=4                                                  STDP-147
      IF (ITT.EQ.16) INL=3                                              STDP-148
      L1=24*NPP+11*(J1-NPP-1)+ITV                                       STDP-149
      GO TO 12                                                          STDP-150
   11 INL=8                                                             STDP-151
      JI=2                                                              STDP-152
      L1=ITV+24*(J1-1)                                                  STDP-153
   12 IJ=1+INL/3                                                        STDP-154
   13 IF (IV.NE.16) GO TO 16                                            STDP-155
C FORM FACTOR COPIED.                                                   STDP-156
      L3=0                                                              STDP-157
      L2=L1+NIJ(JI)*(NVAL(2,NMA+1)-I1)/8                                STDP-158
      DO 15 J=1,IJ                                                      STDP-159
      DO 14 IS=1,ISM                                                    STDP-160
   14 V(IS,L1+L3)=V(IS,L2+L3)*VAL(NMA+3)                                STDP-161
   15 L3=L3+INL                                                         STDP-162
      IZZ(1,L1)=-NVAL(1,NMA+1)                                          STDP-163
      IZZ(2,L1)=IZZ(2,L2)                                               STDP-164
      IZZ(4,L1)=IZZ(4,L2)                                               STDP-165
      GO TO 96                                                          STDP-166
   16 J=MAX0(IV,1)                                                      STDP-167
      NMB=NMA+ITZ(J)                                                    STDP-168
      LT(7)=.FALSE.                                                     STDP-169
      IF (IV.GE.7) GO TO 19                                             STDP-170
      IF (VAL(NMB).NE.0.D0) GO TO 19                                    STDP-171
C ZERO FORM FACTORS AND GO TO NEXT.                                     STDP-172
      L2=L1                                                             STDP-173
      DO 18 J=1,IJ                                                      STDP-174
      DO 17 IS=1,ISM                                                    STDP-175
   17 V(IS,L2)=0.D0                                                     STDP-176
   18 L2=L2+INL                                                         STDP-177
      GO TO 96                                                          STDP-178
   19 LT(1)=LT(1).AND.(.NOT.LO(100))                                    STDP-179
      LT(2)=LT(1).AND.(ITV.EQ.ITT).AND.(NVAL(1,NMA+1).EQ.0)             STDP-180
      LT(3)=.FALSE.                                                     STDP-181
      LT(4)=ITV.GT.6                                                    STDP-182
      LT(8)=(IV.GT.0).AND.(IV.LT.7)                                     STDP-183
      K=IABS(NVAL(1,NMA+2))                                             STDP-184
      HH=WV(8,K)                                                        STDP-185
      H=HH                                                              STDP-186
      HM=DMIN1(H,HM)                                                    STDP-187
      IZZ(1,L1)=-NVAL(1,NMA+1)                                          STDP-188
      IZZ(2,L1)=ITT                                                     STDP-189
      IZZ(4,L1)=K                                                       STDP-190
      IF (NVAL(1,NMA+1).EQ.0) GO TO 20                                  STDP-191
      IJ=1                                                              STDP-192
      NFO=MAX0(NFO,NVAL(1,NMA+1))                                       STDP-193
   20 L=IZZ(3,L1)                                                       STDP-194
      JI=IJ                                                             STDP-195
      IF (LT(2)) JI=JI+1                                                STDP-196
      IF ((IV.GT.0).AND.(IV.LT.9)) GO TO 28                             STDP-197
      AP=VAL(NMB)                                                       STDP-198
      IF (IV.GE.9) GO TO 33                                             STDP-199
      LT(2)=LT(2).AND.(NVAL(2,NMA+1).NE.0)                              STDP-200
      LT(4)=LT(4).AND.((NVAL(2,NMA+1).NE.0))                            STDP-201
C FORM FACTOR INTERPOLATED.                                             STDP-202
      IR=0                                                              STDP-203
      DO 25 IS=1,ISM                                                    STDP-204
      X0=IS*H*VAL(NMB+1)                                                STDP-205
   21 IF (X0.LT.VAL(NMB+IR+6)) GO TO 22                                 STDP-206
      IF (IR.GE.-2*IV-8) GO TO 22                                       STDP-207
      IR=IR+1                                                           STDP-208
      GO TO 21                                                          STDP-209
   22 V(IS,L1)=0.D0                                                     STDP-210
      DO 24 K=2,8,2                                                     STDP-211
      X1=1.D0                                                           STDP-212
      DO 23 J=2,8,2                                                     STDP-213
      IF (K.EQ.J) GO TO 23                                              STDP-214
      X1=X1*(X0-VAL(NMB+IR+J))/(VAL(NMB+IR+K)-VAL(NMB+IR+J))            STDP-215
   23 CONTINUE                                                          STDP-216
   24 V(IS,L1)=V(IS,L1)+X1*VAL(NMB+IR+K+1)                              STDP-217
   25 V(IS,L1)=VAL(NMB)*V(IS,L1)                                        STDP-218
      IF (IJ.EQ.1) GO TO 31                                             STDP-219
   26 L2=L1                                                             STDP-220
      DO 27 J=2,IJ                                                      STDP-221
      CALL DERI(V(1,L2+INL),V(1,L2),H,ISM,.TRUE.)                       STDP-222
   27 L2=L2+INL                                                         STDP-223
      IF (IV.GT.0) GO TO 94                                             STDP-224
      GO TO 31                                                          STDP-225
   28 IF (LT(8)) GO TO 37                                               STDP-226
C BOUND STATE FORM FACTOR.                                              STDP-227
      K=NVAL(1,NMB)                                                     STDP-228
      IVM=NVAL(2,NMB-1)                                                 STDP-229
      IVX=IV                                                            STDP-230
      JV=NVAL(1,NMB+1)/2+1                                              STDP-231
      NMC=NMB+3                                                         STDP-232
      IF (K.EQ.1) GO TO 29                                              STDP-233
      JV=JV+ISM                                                         STDP-234
      IF (JV.GT.NX) CALL MEMO('STDF',NX,JV)                             STDP-235
      IDX=IDX-2*JV                                                      STDP-236
      CALL STBF(P,NVAL(2,NMB),ISM,VAL(NMC+3),IVM,NX-JV,IDX,P(JV+1),IVX,HSTDP-237
     1H,IZZ(1,ITXM+1),LO)                                               STDP-238
      NMC=NMC+8*IV-50                                                   STDP-239
      IF (K.EQ.3) IVX=IVX-1                                             STDP-240
      NMB=NMB+2                                                         STDP-241
   29 CALL STBF(V(1,L1),NVAL(2,NMB),ISM,VAL(NMC),IVM,NX-JV,IDX,P(JV+1),ISTDP-242
     1VX,HH,IZZ(1,ITXM+1),LO)                                           STDP-243
      IF (K.EQ.1) GO TO 96                                              STDP-244
      IDX=IDX+JV                                                        STDP-245
      DO 30 IS=1,ISM                                                    STDP-246
   30 V(IS,L1)=V(IS,L1)*P(IS)*VAL(NMB+3)                                STDP-247
      GO TO 96                                                          STDP-248
   31 LT(3)=.TRUE.                                                      STDP-249
      IF (.NOT.LT(2).OR.LT(4)) GO TO 32                                 STDP-250
      IDX=MAX0(IDX,ISM)                                                 STDP-251
      IF (IDX.GT.2*NX) CALL MEMO('STDF',NX,(IDX+1)/2)                   STDP-252
      CALL DERI(P,V(1,L1),H,ISM,.TRUE.)                                 STDP-253
   32 IF (LT(2).OR.(NVAL(2,NMA+1).NE.0)) GO TO 52                       STDP-254
      IF (IJ.EQ.1) GO TO 96                                             STDP-255
      GO TO 94                                                          STDP-256
C BESSEL EXPANSION.                                                     STDP-257
   33 NMC=NMB+1                                                         STDP-258
      L2X=NVAL(2,NMA+2)                                                 STDP-259
      LL=NVAL(1,NMA+3)                                                  STDP-260
      IF (IV.EQ.9) GO TO 34                                             STDP-261
      JL=NVAL(2,NMA+3)                                                  STDP-262
      LJ=JI+JL-1                                                        STDP-263
      IF (VAL(NMC).EQ.0.D0) VAL(NMC)=1.D0                               STDP-264
      GO TO 52                                                          STDP-265
   34 LJ=LL+JI+NVAL(2,NMA+3)                                            STDP-266
      IDX=MAX0(IDX,3*MAX0(L2X,LJ))                                      STDP-267
      IF (IDX.GT.2*NX) CALL MEMO('STDF',NX,(IDX+1)/2)                   STDP-268
C COMPUTATION OF ZEROS OF BESSEL FUNCTIONS - THE L+1 ZEROS FOR L=1 TO   STDP-269
C L=11 ARE IN DATA ZB. THE OTHERS ARE COMPUTED WITH MC MAHON'S FORMULA  STDP-270
C PAGE 371, HANDBOOK OF MATH. FUNCTIONS, ABRAMOVITZ AND STEGUN.         STDP-271
      IF (VAL(NMC).EQ.0.D0) VAL(NMC)=ISM*H                              STDP-272
      DO 36 II=1,L2X                                                    STDP-273
      IF ((LL.NE.0).AND.((LL.LT.12).AND.(II.LE.LL+1))) GO TO 35         STDP-274
      X(1,II)=(2*II+LL)*PIS2                                            STDP-275
      IF (LL.EQ.0) GO TO 36                                             STDP-276
      A1=2.D0*X(1,II)                                                   STDP-277
      A2=DFLOAT(LL*(LL+1))                                              STDP-278
      X(1,II)=X(1,II)-A2*(1.D0+(7.D0*A2-6.D0+((166.D0*A2-408.D0)*A2+360.STDP-279
     1D0+(((6949.D0*A2-33252.D0)*A2+81180.D0)*A2-75600.D0)/(7.D0*A1**2))STDP-280
     2/(5.D0*A1**2))/(3.D0*A1**2))/A1                                   STDP-281
      GO TO 36                                                          STDP-282
   35 K=(LL*(LL+1))/2+II-1                                              STDP-283
      X(1,II)=ZB(K)                                                     STDP-284
   36 X(1,II)=X(1,II)/VAL(NMC)                                          STDP-285
      LM=LL+1                                                           STDP-286
      LT(6)=LT(2).AND.(NVAL(2,NMA+1).EQ.0)                              STDP-287
      GO TO 52                                                          STDP-288
C WOODS-SAXON AND ITS DERIVATIVES TO SOME POWER.                        STDP-289
   37 NMC=NMB+4                                                         STDP-290
      LT(5)=((VAL(NMB+2).EQ.0.D0).AND.LT(4))                            STDP-291
      LT(4)=((VAL(NMB+2).NE.0.D0).AND.LT(4))                            STDP-292
      LT(6)=(((ITV.EQ.3).OR.(ITV.EQ.4)).AND.(.NOT.LO(100)))             STDP-293
      LT(9)=(.NOT.LO(9)).OR.(VAL(NMB+1).GE.0.D0).OR.LT(5)               STDP-294
      IF (ITV.GT.6) NMC=NMC+1                                           STDP-295
      IF (LT(5)) GO TO 38                                               STDP-296
      SEP=DEXP(HH/VAL(NMB+2))                                           STDP-297
      IF (VAL(NMB+2).GT.0.02D0*H)  GO TO 39                             STDP-298
      WRITE (MW,1000) VAL(NMB+2),NMA,I                                  STDP-299
      VAL(NMB+2)=DMAX1(-VAL(NMB+2),0.02D0*H)                            STDP-300
      GO TO 39                                                          STDP-301
   38 IF (VAL(NMB+1).GE.H) GO TO 39                                     STDP-302
      WRITE (MW,1001) VAL(NMB+1),NMA,I                                  STDP-303
      VAL(NMB+1)=DMAX1(-VAL(NMB+1),H)                                   STDP-304
   39 IF (IV.LE.4) GO TO 49                                             STDP-305
C INITIALISATION OF DO LOOPS FOR DEFORMED POTENTIALS.                   STDP-306
      IQM=NVAL(2,NMA+2)                                                 STDP-307
      IX=20                                                             STDP-308
      IQ=MAX0(L,IQM)                                                    STDP-309
      IF (IV.EQ.6) IQ=MAX0(IQ,NVAL(2,NMA+2))                            STDP-310
      IDX=MAX0(IDX,3*IQ+120)                                            STDP-311
      IF (IDX.GT.2*NX) CALL MEMO('STDF',NX,(IDX+1)/2)                   STDP-312
      A2=0.D0                                                           STDP-313
      A1=1.D0                                                           STDP-314
      IF (LO(6)) A1=DABS(VAL(NMB+1))                                    STDP-315
      DO 40 I=1,IQM                                                     STDP-316
   40 X(1,40+I)=VAL(NMC+I-1)*DSQRT(DFLOAT(2*I+1)/(8.D0*PIS2))/A1        STDP-317
      DO 48 II=1,20                                                     STDP-318
      I=1+MOD(II-1,10)                                                  STDP-319
      X(2,22)=XGN(I)                                                    STDP-320
      IF (I.NE.II) X(2,22)=-X(2,22)                                     STDP-321
      X(2,21)=1.D0                                                      STDP-322
      DO 41 J=2,IQ                                                      STDP-323
   41 X(2,J+21)=(DFLOAT(2*J-1)*X(2,22)*X(2,J+20)+DFLOAT(1-J)*X(2,J+19))/STDP-324
     1DFLOAT(J)                                                         STDP-325
      X(3,II+20)=0.5D0*PGN(I)                                           STDP-326
      X(2,II)=0.D0                                                      STDP-327
      A5=X(2,L+21)                                                      STDP-328
      IF (IV.EQ.5) GO TO 46                                             STDP-329
      KL=NVAL(2,NMA+2)                                                  STDP-330
      IF (NVAL(1,NMA+3).NE.0) GO TO 42                                  STDP-331
      A5=A5*X(2,KL+21)*DSQRT(DFLOAT(2*KL+1))                            STDP-332
      GO TO 46                                                          STDP-333
C COMPUTATION OF Y(L,KK) * Y(KL,KK).                                    STDP-334
   42 A5=DSQRT(DFLOAT(2*KL+1))                                          STDP-335
      KK=NVAL(1,NMA+3)                                                  STDP-336
      DO 43 N=1,KK                                                      STDP-337
   43 A5=A5*(1.D0-XGN(I)**2)*DFLOAT(2*N-1)**2/DSQRT(DFLOAT((KL+N)*(L+N)*STDP-338
     1(KL-N+1)*(L-N+1)))                                                STDP-339
      KZ=L-KK                                                           STDP-340
      DO 45 J=1,2                                                       STDP-341
      IF (KZ.LE.0) GO TO 45                                             STDP-342
      A4=0.D0                                                           STDP-343
      DO 44 K=1,KZ                                                      STDP-344
      A3=A4                                                             STDP-345
      A4=A5                                                             STDP-346
   44 A5=A3+(A4*X(2,22)-A3)*(2.D0*(KK+K)-1.D0)/K                        STDP-347
   45 KZ=KL-KK                                                          STDP-348
   46 X(3,II)=X(3,II+20)*A5*DSQRT(2.D0*L+1.D0)                          STDP-349
      RR=1.D0                                                           STDP-350
      DO 47 K=1,IQM                                                     STDP-351
   47 RR=RR+X(1,40+K)*X(2,K+21)                                         STDP-352
      X(2,II)=-RR*DABS(VAL(NMB+1))                                      STDP-353
      X(1,II)=0.D0                                                      STDP-354
      X(1,20+II)=1.D-16                                                 STDP-355
      IF (.NOT.LT(5)) GO TO 48                                          STDP-356
      X(2,II)=DABS(X(2,II))                                             STDP-357
      A2=A2+X(2,II)**3*X(3,20+II)                                       STDP-358
   48 CONTINUE                                                          STDP-359
      GO TO 50                                                          STDP-360
C INITIALISATION OF DO LOOPS FOR NOT DEFORMED POTENTIALS.               STDP-361
   49 IQM=0                                                             STDP-362
      IX=1                                                              STDP-363
      IDX=MAX0(IDX,33)                                                  STDP-364
      IF (IDX.GT.2*NX) CALL MEMO('STDF',NX,(IDX+1)/2)                   STDP-365
      X(1,1)=0.D0                                                       STDP-366
      X(2,1)=-DABS(VAL(NMB+1))                                          STDP-367
      X(3,1)=1.D0                                                       STDP-368
      IF (.NOT.LT(5)) GO TO 50                                          STDP-369
      X(2,1)=-X(2,1)                                                    STDP-370
      A2=X(2,1)**3                                                      STDP-371
   50 IV1=1+MOD(IV-1,4)                                                 STDP-372
      AP=VAL(NMB)                                                       STDP-373
      AN=VAL(NMB+1)                                                     STDP-374
      IF (LO(6)) AN=AN/DABS(VAL(NMB+1))                                 STDP-375
      IF (.NOT.LT(9)) AN=-AN                                            STDP-376
      IF (IV1.GT.1) AP=AP*AN/DSQRT(8.D0*PIS2)                           STDP-377
      IF (IV1.GT.2) AP=AP*AN*0.5D0                                      STDP-378
      IF (IV1.GT.3) AP=AP*AN/3.D0                                       STDP-379
      IF (VAL(NMB+3).GT.-1.D0) GO TO 51                                 STDP-380
      WRITE (MW,1002) VAL(NMB+3),NMA                                    STDP-381
      VAL(NMB+3)=-.8D0                                                  STDP-382
   51 A1=VAL(NMB+3)                                                     STDP-383
      A5=DFLOAT(L+1)                                                    STDP-384
      IV2=IV1                                                           STDP-385
      IV3=IV2-1                                                         STDP-386
      IV2=IV3+JI                                                        STDP-387
      IF (LT(6)) IV2=IV2+1                                              STDP-388
      IF (ITV.GT.6) AP=AP*CCZ                                           STDP-389
      IF (LT(5)) AP=AP/A2                                               STDP-390
   52 IF (LT(4)) JI=1                                                   STDP-391
      RR=0.D0                                                           STDP-392
      AN=0.D0                                                           STDP-393
      A6=0.D0                                                           STDP-394
      DO 87 IS=1,ISM                                                    STDP-395
      RR=RR+HH                                                          STDP-396
      IF (LT(3)) GO TO 86                                               STDP-397
      DO 53 I=1,3                                                       STDP-398
   53 Y(I)=0.D0                                                         STDP-399
      IF (IV.GE.9) GO TO 67                                             STDP-400
C INTEGRATION LOOP.                                                     STDP-401
      DO 66 I=1,IX                                                      STDP-402
      IF (LT(5)) GO TO 57                                               STDP-403
      X(2,I)=X(2,I)+HH                                                  STDP-404
      IF (X(1,I).NE.0.D0) GO TO 54                                      STDP-405
      IF (X(2,I)+50.D0*VAL(NMB+2).GT.0.D0) X(1,I)=DEXP(X(2,I)/VAL(NMB+2)STDP-406
     1)                                                                 STDP-407
      GO TO 55                                                          STDP-408
   54 IF (X(1,I).LT.1.D15) X(1,I)=X(1,I)*SEP                            STDP-409
   55 CALL WOSA(VR,VAL(NMB+2),A1,X(1,I),X(2,I),IV2,LT(9))               STDP-410
      IF (LT(6)) VR(IV1)=VR(IV2)*4.D0*VAL(NMB+2)                        STDP-411
      IF (LT(4)) A6=A6+VR(1)*RR**2*X(3,I+20)*(1.D0+VAL(NMB+4)*RR**2)    STDP-412
      DO 56 J=1,JI                                                      STDP-413
   56 Y(J)=Y(J)+VR(IV3+J)*X(3,I)                                        STDP-414
      GO TO 66                                                          STDP-415
C DEFORMED COULOMB POTENTIAL.                                           STDP-416
   57 A4=X(2,I)/RR                                                      STDP-417
      IF (ITV.NE.ITT) GO TO 59                                          STDP-418
      IF (A4.GT.1.D0) GO TO 58                                          STDP-419
      Y(1)=Y(1)+(X(2,I)**2)*A4*X(3,I)                                   STDP-420
      Y(2)=Y(2)+X(2,I)*A4**2*X(3,I)                                     STDP-421
      Y(3)=Y(3)+2.D0*A4**3*X(3,I)                                       STDP-422
      GO TO 66                                                          STDP-423
   58 Y(1)=Y(1)+(0.5D0*X(2,I)*X(2,I)-RR*RR/6.D0)*X(3,I)*3.D0            STDP-424
      Y(2)=Y(2)+RR*X(3,I)                                               STDP-425
      Y(3)=Y(3)-X(3,I)                                                  STDP-426
      GO TO 66                                                          STDP-427
   59 IF (A4.GT.1.D0) GO TO 60                                          STDP-428
      A3=(X(2,I)**2)*(A4**(L+1))*3.D0/((A5+2.D0)*(2.D0*A5-1.D0))        STDP-429
      IF (IV1.GT.1) A3=A3*(A5+2.D0)/X(2,I)                              STDP-430
      IF (IV1.GT.2) A3=A3*(A5+1.D0)/X(2,I)                              STDP-431
      IF (IV1.GT.3) A3=A3*A5/X(2,I)                                     STDP-432
      GO TO 62                                                          STDP-433
   60 IF (L.NE.2) GO TO 61                                              STDP-434
      IF (IV1.EQ.1) A3=RR*RR*(0.2D0+DLOG(A4))*0.6D0                     STDP-435
      IF (IV1.GE.2) A3=0.6D0*RR/A4                                      STDP-436
      IF (IV1.GE.3) A3=-A3/X(2,I)                                       STDP-437
      IF (IV1.GE.4) A3=-2.D0*A3/X(2,I)                                  STDP-438
      GO TO 62                                                          STDP-439
   61 IF (IV1.EQ.1) A3=RR*RR*(1.D0/(A5+2.D0)-1.D0/(A4**(L-2)*(2.D0*A5-1.STDP-440
     1D0)))*3.D0/(A5-3.D0)                                              STDP-441
      IF (IV1.GE.2) A3=RR/A4**(L-1)*3.D0/(2.D0*A5-1.D0)                 STDP-442
      IF (IV1.GE.3) A3=-A3*(A5-2.D0)/X(2,I)                             STDP-443
      IF (IV1.EQ.4) A3=-A3*(A5-1.D0)/X(2,I)                             STDP-444
   62 Y(1)=Y(1)+A3*X(3,I)                                               STDP-445
      IF (JI.EQ.1) GO TO 66                                             STDP-446
      IF (A4.GT.1.D0) GO TO 63                                          STDP-447
      A3=-A5*A3/RR                                                      STDP-448
      GO TO 65                                                          STDP-449
   63 IF (L.NE.2) GO TO 64                                              STDP-450
      IF (IV2.EQ.1) A3=-1.2D0*(0.3D0*RR-DLOG(A4))*RR                    STDP-451
      IF (IV2.GE.2) A3=2.D0*A3/RR                                       STDP-452
      GO TO 65                                                          STDP-453
   64 IF (IV2.EQ.1) A3=(2.D0/(A5+2.D0)-(A5-1.D0)/(A4**(L-2)*(2.D0*A5-1.DSTDP-454
     10)))*3.D0/(A5-3.D0)*RR                                            STDP-455
      IF (IV2.NE.1) A3=(A5-1.D0)*A3/RR                                  STDP-456
   65 Y(2)=Y(2)-A3*X(3,I)                                               STDP-457
   66 CONTINUE                                                          STDP-458
      GO TO 84                                                          STDP-459
C COMPUTATION OF BESSEL FUNCTIONS.                                      STDP-460
   67 IF (IV.GT.9) GO TO 76                                             STDP-461
      IF (RR.GT.VAL(NMC)) GO TO 87                                      STDP-462
      DO 75 II=1,L2X                                                    STDP-463
      A1=RR*X(1,II)                                                     STDP-464
      X(2,1)=DSIN(A1)/A1                                                STDP-465
      IF (LJ.EQ.1) GO TO 72                                             STDP-466
      K=IDINT(1.D0+A1)                                                  STDP-467
      IF (K.LT.LJ) GO TO 69                                             STDP-468
      X(2,2)=(X(2,1)-DCOS(A1))/A1                                       STDP-469
      IF (LJ.EQ.2) GO TO 72                                             STDP-470
      DO 68 J=3,LJ                                                      STDP-471
   68 X(2,J)=(2*J-3)*X(2,J-1)/A1-X(2,J-2)                               STDP-472
      GO TO 72                                                          STDP-473
   69 A3=LJ                                                             STDP-474
      A4=DMAX1(DSQRT(10.5D0*A1)-0.5D0,A3)                               STDP-475
      K=IDINT(A4+3.D0+21.D0*A1/(A4+A4+1.D0))                            STDP-476
      A2=0.D0                                                           STDP-477
   70 A2=A1/(2.D0*K+1.D0-A2*A1)                                         STDP-478
      IF (K.LE.LJ) X(2,K+1)=A2                                          STDP-479
      K=K-1                                                             STDP-480
      IF (K.GE.1) GO TO 70                                              STDP-481
      DO 71 K=2,LJ                                                      STDP-482
   71 X(2,K)=X(2,K)*X(2,K-1)                                            STDP-483
C COMPUTATION OF DERIVATIVES OF BESSEL FUNCTIONS (- DERIVATIVE).        STDP-484
   72 JL=LJ-1                                                           STDP-485
      DO 74 J=LM,LJ                                                     STDP-486
      KK=JI+J-LJ                                                        STDP-487
      IF (KK.GE.1) Y(KK)=Y(KK)+VAL(NMC+II)*X(2,LM)                      STDP-488
      IF (KK.EQ.JI) GO TO 75                                            STDP-489
      A3=0.D0                                                           STDP-490
      A5=0.D0                                                           STDP-491
      DO 73 K=1,JL                                                      STDP-492
      A2=(A5+1.D0)*X(2,K+1)-A5*A3                                       STDP-493
      A3=X(2,K)                                                         STDP-494
      X(2,K)=X(1,II)*A2/(2.D0*A5+1.D0)                                  STDP-495
   73 A5=A5+1.D0                                                        STDP-496
   74 JL=JL-1                                                           STDP-497
   75 CONTINUE                                                          STDP-498
      GO TO 84                                                          STDP-499
C COMPUTATION OF LAGUERRE POLYNOMIALS X**LL L(2X**2) DEXP(-X**2).       STDP-500
   76 DO 83 J=1,L2X,2                                                   STDP-501
      IF (J.NE.1) GO TO 77                                              STDP-502
      A1=RR/VAL(NMC)                                                    STDP-503
      A2=A1*A1                                                          STDP-504
      CL(1)=A1**LL*DEXP(-0.5D0*A2)                                      STDP-505
      CL(2)=(LL+1.5D0-A2)*CL(1)                                         STDP-506
      GO TO 78                                                          STDP-507
   77 CL(1)=(CL(2)*(DFLOAT(LL+2*J)-2.5D0-A2)-CL(1)*(DFLOAT(LL+J)-1.5D0))STDP-508
     1/DFLOAT(J-1)                                                      STDP-509
      CL(2)=(CL(1)*(DFLOAT(LL+2*J)-0.5D0-A2)-CL(2)*(DFLOAT(LL+J)-0.5D0))STDP-510
     1/DFLOAT(J)                                                        STDP-511
   78 IF (JL.NE.0) GO TO 79                                             STDP-512
      Y(1)=Y(1)+CL(1)*VAL(NMC+J)                                        STDP-513
      IF (J.LT.L2X) Y(1)=Y(1)+CL(2)*VAL(NMC+J+1)                        STDP-514
   79 IF (LJ.EQ.0) GO TO 83                                             STDP-515
      DO 80 L=3,6                                                       STDP-516
   80 CL(L)=0.D0                                                        STDP-517
C -DERIVATIVE OF LAGUERRE POLYNOMIALS X**LL L(2X**2) DEXP(-X**2).       STDP-518
      DO 82 K=1,LJ                                                      STDP-519
      DO 81 L=1,6                                                       STDP-520
   81 CL(L+2)=CL(L)                                                     STDP-521
      CL(1)=((DFLOAT(LL+2*J+K)-A2)*CL(3)-DFLOAT(2*J)*CL(4)-DFLOAT(2*(K-1STDP-522
     1)*(K-2))*CL(7))/A1+DFLOAT(4*(K-1))*CL(5)                          STDP-523
      CL(2)=-((DFLOAT(LL+2*J+1-K)-A2)*CL(4)-DFLOAT(2*J+1+2*LL)*CL(3)-DFLSTDP-524
     1OAT(2*(K-1)*(K-2))*CL(8))/A1-DFLOAT(4*(K-1))*CL(6)                STDP-525
      IF (K.LT.JL) GO TO 83                                             STDP-526
      Y(K+1-JL)=Y(K+1-JL)+CL(1)*VAL(NMC+J)                              STDP-527
      IF (J.LT.L2X) Y(K+1-JL)=Y(K+1-JL)+CL(2)*VAL(NMC+J+1)              STDP-528
   82 CONTINUE                                                          STDP-529
   83 CONTINUE                                                          STDP-530
   84 L2=L1                                                             STDP-531
      DO 85 J=1,JI                                                      STDP-532
      V(IS,L2)=AP*Y(J)                                                  STDP-533
   85 L2=L2+INL                                                         STDP-534
      GO TO 87                                                          STDP-535
   86 Y(1)=V(IS,L1)/AP                                                  STDP-536
      IF (LT(2).AND.(.NOT.LT(4))) V(IS,L1)=P(IS)                        STDP-537
   87 AN=AN+AP*Y(1)*RR**(L+2)                                           STDP-538
      IF (LT(7)) GO TO 94                                               STDP-539
      IF (.NOT.LT(4)) GO TO 89                                          STDP-540
C FOLDING OF CHARGE DISTRIBUTION WITH COULOMB POTENTIAL.                STDP-541
      IDX=MAX0(IDX,5*ISM)                                               STDP-542
      IF (IDX.GT.2*NX) CALL MEMO('STDF',NX,(IDX+1)/2)                   STDP-543
      A7=0.D0                                                           STDP-544
      IF (LT(8)) A7=VAL(NMB+4)                                          STDP-545
      CALL COPO(V(1,L1),V(1,L1),X,ISM,HH,L,AP,A7,CCZ,A6,.TRUE.,.FALSE.) STDP-546
      IF (A6.EQ.0.D0) A6=DABS(V(ISM,L1)*DFLOAT(2*L+1)*(ISM*HH)**(L+1)/(CSTDP-547
     1CZ*VAL(NMB)))                                                     STDP-548
      DO 88 IS=1,ISM                                                    STDP-549
   88 V(IS,L1)=V(IS,L1)/A6                                              STDP-550
   89 IF (NVAL(2,NMA+1).EQ.0) GO TO 93                                  STDP-551
      AN=DABS(AN*HH/AP)                                                 STDP-552
      IF (ITV.GT.6) AN=DABS(V(ISM,L1)*DFLOAT(2*L+1)*(ISM*HH)**(L+1)/(CCZSTDP-553
     1*VAL(NMB)))                                                       STDP-554
      IF (LT(2).AND.(ITT.EQ.8).AND.(NVAL(1,NMA+1).EQ.0).AND.(VAL(NMB+2).STDP-555
     1EQ.0.D0).AND.LT(8)) AN=AN*HH*ISM                                  STDP-556
      IF (NVAL(2,NMA+1).GT.0) GO TO 90                                  STDP-557
      NVAL(2,NMA+1)=-NVAL(2,NMA+1)                                      STDP-558
      AP=VAL(NMB)                                                       STDP-559
      VAL(NMB)=AP*AN                                                    STDP-560
      WRITE (MW,1003) AP,VAL(NMB),ITT,J1                                STDP-561
      GO TO 93                                                          STDP-562
   90 L2=L1                                                             STDP-563
      DO 92 J=1,IJ                                                      STDP-564
      DO 91 IS=1,ISM                                                    STDP-565
   91 V(IS,L2)=V(IS,L2)/AN                                              STDP-566
   92 L2=L2+INL                                                         STDP-567
   93 IF (.NOT.LT(4)) GO TO 94                                          STDP-568
      LT(4)=.FALSE.                                                     STDP-569
      LT(7)=.TRUE.                                                      STDP-570
      IF (LT(2)) GO TO 31                                               STDP-571
      IF (IJ.NE.1) GO TO 26                                             STDP-572
   94 IF (.NOT.LT(1).OR.(NVAL(1,NMA+1).NE.0)) GO TO 96                  STDP-573
      RR=0.D0                                                           STDP-574
      DO 95 IS=1,ISM                                                    STDP-575
      RR=RR+H                                                           STDP-576
      A1=V(IS,L1+INL)                                                   STDP-577
      V(IS,L1+INL)=-V(IS,L1)/RR**2                                      STDP-578
   95 V(IS,L1)=A1/RR                                                    STDP-579
   96 NMA=NVAL(2,I1-6)+1                                                STDP-580
      GO TO 9                                                           STDP-581
   97 IF (NFO.EQ.0) RETURN                                              STDP-582
C FOLDING.                                                              STDP-583
      NNF=ISM*ITXM+1                                                    STDP-584
      IF (LO(100)) NNF=NNF-ISM*ITX(7)                                   STDP-585
      DO 98 I=1,NNF                                                     STDP-586
   98 P(I)=0.D0                                                         STDP-587
      N1=NVAL(2,1)                                                      STDP-588
      IST=0                                                             STDP-589
      DO 99 N=1,NFO                                                     STDP-590
      IS=5+ISM+IDINT(2.D0*(DABS(VAL(N1+3*N-2))+2.D0*DABS(VAL(N1+3*N-1)))STDP-591
     1/HM)                                                              STDP-592
      IF (VAL(N1+3*N-3)*VAL(N1+3*N-2).EQ.0.D0) IS=ISM+5                 STDP-593
   99 IST=MAX0(IS,IST)                                                  STDP-594
      IDY=NNF+2*IST*(IMAX+1)                                            STDP-595
      IDX=MAX0(IDY,IDX)                                                 STDP-596
      IF (IDX.GT.2*NX) CALL MEMO('STDF',NX,(IDX+1)/2)                   STDP-597
      CALL FOLD(V,P,VAL(N1),NFO,1,ISM,IST,IVY,INTC,P(NNF),PGN,XGN,WV,IZZSTDP-598
     1,LO)                                                              STDP-599
      RETURN                                                            STDP-600
 1000 FORMAT (' TOO SMALL DIFFUSENESS =',D15.8,'   CHANGED INTO MINIMUM STDP-601
     1VALUE IN STDP FOR NMA =',I5,' AND I =',I2)                        STDP-602
 1001 FORMAT (' TOO SMALL COULOMB RADIUS =',D15.8,'   CHANGED INTO MINIMSTDP-603
     1UM VALUE IN STDP FOR NMA =',I5,' AND I =',I2)                     STDP-604
 1002 FORMAT (' POWER 1+',D15.6,'  CHANGED TO .2 FOR NMA =',I4)         STDP-605
 1003 FORMAT (' STRENGTH',D15.6,' REPLACED BY THE INTEGRAL',D15.6,' FOR STDP-606
     1THE FORM FACTOR',I4,' OF THE POTENTIAL',I4)                       STDP-607
      END                                                               STDP-608
