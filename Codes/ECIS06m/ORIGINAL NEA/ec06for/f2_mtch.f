C 08/03/07                                                      ECIS06  MTCH-000
      SUBROUTINE MTCH(NVI,NCOLL,KAB,WV,MC,BG,NAT,AT,AA,ISM,LMAX2,NIV,IVZMTCH-001
     1,FG,LMAX1,LMAX3,KR,VCO,VDO,FAM,JD,LO)                             MTCH-002
C COMPUTATION OF COULOMB CORRECTIONS AS INTEGRALS FROM THE MATCHING     MTCH-003
C POINT IF LO(127)=.FALSE. OR FROM THE ORIGIN IF LO(127)=.TRUE..        MTCH-004
C INPUT:     NVI:     ADDRESSES OF COUPLINGS IN TABLE AT.               MTCH-005
C            NCOLL:   NUMBER OF NUCLEAR STATES.                         MTCH-006
C            KAB:     MAXIMUM NUMBER OF COUPLED CHANNELS.               MTCH-007
C            WV:      WAVE NUMBERS AND COULOMB PARAMETERS.              MTCH-008
C            MC:      NUCLEAR STATE NUMBER AND ANGULAR MOMENTUM.        MTCH-009
C            BG:      TABLE OF COULOMB INTEGRALS FROM 0 TO INFINITY     MTCH-010
C                     FOR COUPLED EQUATIONS OR WHEN LO(127)=.TRUE..     MTCH-011
C            NAT,AT:  TABLE OF COUPLING COEFFICIENTS.                   MTCH-012
C            ISM:     NUMBER OF INTEGRATION POINTS.                     MTCH-013
C            LMAX2:   FIRST DIMENSION OF TABLE BG.                      MTCH-014
C            NIV:     ADDRESS IN THE TABLE OF REDUCED MATRIX ELEMENT.   MTCH-015
C            IVZ:     TABLE FOR USE OF FORM FACTORS (SEE REDM 3RD PART).MTCH-016
C            FG:      COULOMB FUNCTIONS.                                MTCH-017
C            LMAX1:   FIRST DIMENSION OF TABLE FG.                      MTCH-018
C            LMAX3:   MAXIMUM NUMBER OF COULOMB INTEGRALS FROM THE      MTCH-019
C                     MATCHING POINT TO INFINITY.                       MTCH-020
C            VCO:     STRENGTH OF TAILS OF COULOMB POTENTIALS.          MTCH-021
C            VDO:     STRENGTH OF TAILS OF COULOMB TRANSITIONS.         MTCH-022
C            FAM(I,J):OUTPUT OF PREVIOUS CALL FROM WHICH J=1 TO 10 ARE  MTCH-023
C                     NEEDED IF UNCOUPLED SOLUTIONS ARE COPIED.         MTCH-024
C            JD:      FIRST DIMENSION OF NAT,ST.                        MTCH-025
C            LO(I):   LOGICAL CONTROLS:                                 MTCH-026
C               LO(11) =.TRUE. DEFORMED COULOMB POTENTIAL.              MTCH-027
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   MTCH-028
C               LO(29) =.TRUE. NO DIAGONAL TERMS IN SECOND MEMBER.      MTCH-029
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     MTCH-030
C               LO(57) =.TRUE. OUTPUT PHASE-SHIFTS AT EACH ITERATION.   MTCH-031
C               LO(58) =.TRUE. OUTPUT OF THE COEFFICIENTS OF EACH FORM  MTCH-032
C                              FACTOR FOR ALL SETS OF EQUATIONS.        MTCH-033
C               LO(100)=.TRUE. DIRAC EQUATION.                          MTCH-034
C               LO(125)=.TRUE. USUAL COUPLED EQUATIONS.                 MTCH-035
C               LO(127)=.TRUE. COULOMB CORRECTIONS IN ASYMPTOTIC REGION.MTCH-036
C OUTPUT:    AA:      COULOMB INTEGRALS FROM THE MATCHING POINT TO      MTCH-037
C                     INFINITY FOR ITERATIONS LO(125)=.FALSE..          MTCH-038
C            FAM(I,J):MATCHING VALUES AND WAVE NUMBER FOR J=1 TO 6,     MTCH-039
C                     COEFFICIENT OF CENTRAL POTENTIAL FOR J=7,         MTCH-040
C                     COEFFICIENT OF SPIN-ORBIT POTENTIAL FOR J=9,      MTCH-041
C                     ENERGY FOR J=8, CENTRIFUGAL POTENTIAL FOR J=10.   MTCH-042
C                     IF LO(133)=.TRUE.,COEFFICIENT OF COULOMB POTENTIALMTCH-043
C                     FOR J=11 AND OF SPIN-ORBIT COULOMB FOR J=12.      MTCH-044
C WORKING AREAS:                                                        MTCH-045
C            AA(1,1,I):FOR I=7,10 IF LO(127) OR LO(125)=.TRUE.          MTCH-046
C            KR:      WORKING AREA IN SUBROUTINE LINS.                  MTCH-047
C                                                                       MTCH-048
C LOCAL TABLES LA(3,11) AND BA(2,11) ARE SET FOR COULOMB CORRECTIONS    MTCH-049
C LIMITED BY MCM(1)=5 AND MCM(2)=4.                                     MTCH-050
C                                                                       MTCH-051
C FOR THE COMMON  /DCONS/ SEE CALC.                                     MTCH-052
C FOR THE COMMON  /NOEQU/ SEE QUAN.                                     MTCH-053
C FOR THE COMMON  /POTE2/ SEE REDM.                                     MTCH-054
C                                                                       MTCH-055
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     MTCH-056
C  CHB:       PLANCK CONSTANT /(2*PI) IN MEV*FERMI.                     MTCH-057
C   USED:     CHB.                                                      MTCH-058
C                                                                       MTCH-059
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NOEQU/:                     MTCH-060
C  NCXN:      NUMBER OF SOLUTIONS NEEDED.                               MTCH-061
C  NC:        NUMBER OF EQUATIONS FOR IDENTICAL PARTICLES.              MTCH-062
C   USED:     NCXN,NC.                                                  MTCH-063
C                                                                       MTCH-064
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     MTCH-065
C  ITY(5):    STARTING ADDRESS OF REAL CENTRAL TRANSITION.              MTCH-066
C        ITY(5)=0 IS USED FOR DIRAC EQUATIONS.                          MTCH-067
C  INTV:      NUMBER OF TRANSITION FORM FACTORS WITHOUT SPIN-ORBIT,     MTCH-068
C             TAKING INTO ACCOUNT DISPERSION.                           MTCH-069
C  INSL:      NUMBER OF SPIN-ORBIT FORM TRANSITION FACTORS NOT TAKING   MTCH-070
C                  INTO ACCOUNT MULTIPLICATION BY 2.                    MTCH-071
C   USED:     ITY(5),INTV,INSL.                                         MTCH-072
C                                                                       MTCH-073
C***********************************************************************MTCH-074
      IMPLICIT REAL*8 (A-H,O-Z)                                         MTCH-075
      LOGICAL LO(150),LV                                                MTCH-076
      DIMENSION NVI(KAB,KAB,3),WV(22,*),MC(KAB,6),BG(LMAX2,*),NAT(2*JD,*MTCH-077
     1),AT(JD,*),AA(KAB,KAB,*),NIV(NCOLL,NCOLL,3),IVZ(4,*),FG(LMAX1,4,*)MTCH-078
     2,VCO(2,*),VDO(2,*),FAM(KAB,12),B(4),C(4),G(4),AB(4,2),KR(*),AV(5),MTCH-079
     3LA(3,11),BA(2,11),N1(4),N2(4),N3(4),N4(4)                         MTCH-080
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            MTCH-081
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         MTCH-082
      COMMON /NOEQU/ NCXN,NIC,NCI,NC,NCIN,NIN,JPI,IPJ,R1(2),NAJ         MTCH-083
      COMMON /INOUT/ MR,MW,MS                                           MTCH-084
      DATA N1,N2 /1,1,3,3,1,3,1,3/                                      MTCH-085
      DATA N3,N4 /1,1,2,2,1,2,1,2/                                      MTCH-086
      LV=LO(44).AND.(LO(57).OR.LO(58))                                  MTCH-087
      LL4=4                                                             MTCH-088
      IF (LO(127)) LL4=1                                                MTCH-089
      IF (LV) WRITE (MW,1000)                                           MTCH-090
C LOOPS ON EQUATIONS.                                                   MTCH-091
      DO 55 IC=1,NC                                                     MTCH-092
      I1=MC(IC,1)                                                       MTCH-093
      K1=MC(IC,4)                                                       MTCH-094
      L1=MC(IC,2)                                                       MTCH-095
      IF (K1.GE.0.OR.LO(127)) GO TO 2                                   MTCH-096
C TRANSFER OF INFORMATIONS WHEN UNCOUPLED FUNCTIONS ARE NOT RECOMPUTED. MTCH-097
      K1=-K1                                                            MTCH-098
      IF (.NOT.LO(100)) MC(IC,6)=0                                      MTCH-099
      IF (K1.EQ.IC) GO TO 5                                             MTCH-100
      DO 1 I=1,10                                                       MTCH-101
    1 FAM(IC,I)=FAM(K1,I)                                               MTCH-102
      GO TO 5                                                           MTCH-103
C VALUES OF LONG RANGE TAILS OF CENTRAL POTENTIALS.                     MTCH-104
    2 F2=VCO(1,K1)**2                                                   MTCH-105
      F3=VCO(2,K1)*MC(IC,6)                                             MTCH-106
      JC=0                                                              MTCH-107
      IF (I1.LE.NCOLL) JC=NIV(I1,I1,3)                                  MTCH-108
      IF (WV(3,I1).GT.0.D0.AND.JC.GT.0) GO TO 3                         MTCH-109
      F2=0.D0                                                           MTCH-110
      F3=0.D0                                                           MTCH-111
    3 IF (LO(127)) GO TO 5                                              MTCH-112
C INTEGRATION REGION - SET UP OF POTENTIAL IN FIVE POINTS FOR           MTCH-113
C TRANSFORMATION OF MATCHING VALUES.                                    MTCH-114
      B1=WV(8,1)*WV(8,1)/48.D0                                          MTCH-115
      C1=(ISM-1)*WV(8,1)                                                MTCH-116
      A1=WV(11,I1)**2                                                   MTCH-117
      IF (WV(3,I1).LT.0.D0) A1=-A1                                      MTCH-118
      DO 4 I=1,5                                                        MTCH-119
      AV(I)=B1*(2.D0*WV(11,I1)*WV(5,I1)/C1-A1+(MC(IC,5)-F2-F3/C1)/C1**2)MTCH-120
    4 C1=C1+0.5D0*WV(8,1)                                               MTCH-121
C COMPUTATION OF COULOMB CORRECTIONS.                                   MTCH-122
    5 AZ=ISM*WV(8,I1)                                                   MTCH-123
      DO 47 IP=1,IC                                                     MTCH-124
      DO 7 K=1,2                                                        MTCH-125
      DO 6 I=1,4                                                        MTCH-126
    6 AB(I,K)=0.D0                                                      MTCH-127
    7 CONTINUE                                                          MTCH-128
      I2=MC(IP,1)                                                       MTCH-129
      ILM=1                                                             MTCH-130
      IF (.NOT.LO(44)) GO TO 44                                         MTCH-131
      L2=MC(IP,2)                                                       MTCH-132
      I3=NIV(I2,I1,3)                                                   MTCH-133
      IF (I3.EQ.0) GO TO 44                                             MTCH-134
      AY=DSQRT(WV(11,I1)*WV(11,I2))                                     MTCH-135
      AW=AY*AZ                                                          MTCH-136
C SCAN THE COUPLINGS BETWEEN EQUATIONS.                                 MTCH-137
      IS=0                                                              MTCH-138
      IF ((.NOT.LO(100)).AND.LO(19).AND.(IC.NE.IP)) ILM=2               MTCH-139
      K1=NVI(IP,IC,1)                                                   MTCH-140
      K2=NVI(IP,IC,2)                                                   MTCH-141
      IF ((.NOT.LO(11)).OR.(K1.GT.K2)) GO TO 16                         MTCH-142
      IML=1                                                             MTCH-143
C CENTRAL CONTRIBUTION.                                                 MTCH-144
      I=1                                                               MTCH-145
    8 IF (VDO(1,I).EQ.0.D0) GO TO 30                                    MTCH-146
      II=I+ITY(5)                                                       MTCH-147
      DO 9 K=K1,K2                                                      MTCH-148
      IF (NAT(1,K).EQ.II) GO TO 10                                      MTCH-149
    9 CONTINUE                                                          MTCH-150
      GO TO 30                                                          MTCH-151
   10 IF ((AT(2,K).EQ.0.D0).AND.(.NOT.LO(100))) GO TO 30                MTCH-152
      LA1=MC(IC,2)                                                      MTCH-153
      LA2=MC(IP,2)                                                      MTCH-154
      LA3=IVZ(4,I)                                                      MTCH-155
      ZT=AT(2,K)*AY**LA3*VDO(1,I)/WV(8,I1)**2                           MTCH-156
      IF (LO(100)) ZT=ZT*.5D0*DSQRT((1.D0+WV(6,I1)/WV(7,I1))*(1.D0+WV(6,MTCH-157
     1I2)/WV(7,I1)))                                                    MTCH-158
      IM=1                                                              MTCH-159
   11 IF (IM.LT.5) YT=ZT                                                MTCH-160
      IF ((YT.EQ.0.D0).AND.(ZT.EQ.0.D0)) GO TO 14                       MTCH-161
      IF (IS.EQ.0) GO TO 13                                             MTCH-162
      DO 12 J=1,IS                                                      MTCH-163
      IF ((LA(1,J).NE.LA1).OR.(LA(2,J).NE.LA2).OR.(LA(3,J).NE.LA3)) GO TMTCH-164
     1O 12                                                              MTCH-165
      BA(1,J)=BA(1,J)+ZT                                                MTCH-166
      BA(2,J)=BA(2,J)+YT                                                MTCH-167
      GO TO 14                                                          MTCH-168
   12 CONTINUE                                                          MTCH-169
   13 IS=IS+1                                                           MTCH-170
      LA(1,IS)=LA1                                                      MTCH-171
      LA(2,IS)=LA2                                                      MTCH-172
      LA(3,IS)=LA3                                                      MTCH-173
      BA(1,IS)=ZT                                                       MTCH-174
      BA(2,IS)=YT                                                       MTCH-175
   14 GO TO ( 15 , 30 , 20 , 30 , 27 , 28 , 29 , 30 ),IM                MTCH-176
   15 IF (.NOT.LO(100)) GO TO 30                                        MTCH-177
C VECTOR CONTRIBUTION IN DIRAC EQUATION.                                MTCH-178
      IM=2                                                              MTCH-179
      LA1=MC(IC,3)-MC(IC,2)                                             MTCH-180
      LA2=MC(IP,3)-MC(IP,2)                                             MTCH-181
      LA3=IVZ(4,I)                                                      MTCH-182
      ZT=AT(4,K)*AY**LA3*VDO(1,I)*CHB**2*WV(11,I1)*WV(11,I2)/(2.D0*WV(7,MTCH-183
     1I1)*(WV(6,I2)+WV(7,I1))*WV(8,I1)**2)                              MTCH-184
      IF ((LA1-MC(IC,2))*(LA2-MC(IP,2)).LT.0) ZT=-ZT                    MTCH-185
      GO TO 11                                                          MTCH-186
   16 IF (.NOT.LO(19)) GO TO 31                                         MTCH-187
      K3=NVI(IP,IC,3)                                                   MTCH-188
C TENSOR CONTRIBUTION IN DIRAC EQUATION.                                MTCH-189
      IF (.NOT.LO(100)) GO TO 21                                        MTCH-190
      IML=2                                                             MTCH-191
      K4=K2+1                                                           MTCH-192
      IF (K4.GT.K3) GO TO 31                                            MTCH-193
      I=1                                                               MTCH-194
   17 IF (VDO(2,I).EQ.0.D0) GO TO 30                                    MTCH-195
      DO 18 K=K4,K3                                                     MTCH-196
      IF (NAT(1,K).EQ.IVZ(3,I)+INTV) GO TO 19                           MTCH-197
   18 CONTINUE                                                          MTCH-198
      GO TO 30                                                          MTCH-199
   19 IM=3                                                              MTCH-200
      LA1=MC(IC,3)-MC(IC,2)                                             MTCH-201
      LA2=MC(IP,2)                                                      MTCH-202
      LA3=IVZ(4,I)+1                                                    MTCH-203
      ZT=.5D0*(AT(4,K)+DFLOAT(LA3)*AT(2,K))*AY**LA3*VDO(2,I)*CHB*WV(11,IMTCH-204
     11)/(WV(8,I1)**2*WV(7,I1))                                         MTCH-205
      IF ((LA1-MC(IC,2)).GT.0) ZT=-ZT                                   MTCH-206
      GO TO 11                                                          MTCH-207
   20 IM=4                                                              MTCH-208
      LA1=MC(IC,2)                                                      MTCH-209
      LA2=MC(IP,3)-LA2                                                  MTCH-210
      ZT=-.5D0*(AT(4,K)-DFLOAT(LA3)*AT(2,K))*AY**LA3*VDO(2,I)*CHB*WV(11,MTCH-211
     1I2)/(WV(8,I2)**2*WV(7,I1))                                        MTCH-212
      IF ((LA2-MC(IP,2)).GT.0) ZT=-ZT                                   MTCH-213
      GO TO 11                                                          MTCH-214
   21 I=1                                                               MTCH-215
      IML=3                                                             MTCH-216
C SPIN-ORBIT CONTRIBUTION IN SCHROEDINGER EQUATION.                     MTCH-217
      KP1=NVI(IC,IP,1)                                                  MTCH-218
      KP2=NVI(IC,IP,2)                                                  MTCH-219
      KP3=NVI(IC,IP,3)                                                  MTCH-220
   22 IF (VDO(2,I).EQ.0.D0) GO TO 30                                    MTCH-221
      II=IVZ(3,I)+ITY(5)+INSL+INTV                                      MTCH-222
      AX1=0.D0                                                          MTCH-223
      AX2=0.D0                                                          MTCH-224
      AX3=0.D0                                                          MTCH-225
      IF (K1.GT.K3) GO TO 24                                            MTCH-226
      DO 23 K=K1,K3                                                     MTCH-227
      IF (NAT(1,K)+INSL.EQ.II) AX1=AT(2,K)                              MTCH-228
      IF ((K.LE.K2).AND.(NAT(1,K).EQ.II)) AX2=AT(2,K)                   MTCH-229
      IF ((K.GT.K2).AND.(NAT(1,K).EQ.II)) AX3=AT(2,K)                   MTCH-230
   23 CONTINUE                                                          MTCH-231
   24 AY1=0.D0                                                          MTCH-232
      AY2=0.D0                                                          MTCH-233
      AY3=0.D0                                                          MTCH-234
      IF (KP1.GT.KP3) GO TO 26                                          MTCH-235
      DO 25 K=KP1,KP3                                                   MTCH-236
      IF (NAT(1,K)+INSL.EQ.II) AY1=AT(2,K)                              MTCH-237
      IF ((K.LE.KP2).AND.(NAT(1,K).EQ.II)) AY2=AT(2,K)                  MTCH-238
      IF ((K.GT.KP2).AND.(NAT(1,K).EQ.II)) AY3=AT(2,K)                  MTCH-239
   25 CONTINUE                                                          MTCH-240
   26 IF ((K1.GT.K3).AND.(KP1.GT.KP3)) GO TO 30                         MTCH-241
      LA1=MC(IC,2)                                                      MTCH-242
      LA2=MC(IP,2)                                                      MTCH-243
      LA3=IVZ(4,I)+2                                                    MTCH-244
      IF (LA3.GT.LA1+LA2) GO TO 30                                      MTCH-245
      XT=AY**LA3*VDO(2,I)/WV(8,I1)**2                                   MTCH-246
      ZT=(AX2-DFLOAT(LA3-1)*AX1+DFLOAT(LA1+1)*AX3)*XT                   MTCH-247
      YT=(AY2-DFLOAT(LA3-1)*AY1+DFLOAT(LA2+1)*AY3)*XT                   MTCH-248
      IM=5                                                              MTCH-249
      GO TO 11                                                          MTCH-250
   27 LA3=LA3-1                                                         MTCH-251
      XT=XT/AY                                                          MTCH-252
      ZT=AX3*WV(11,I1)*WV(5,I1)/DFLOAT(LA1+1)*XT                        MTCH-253
      YT=AY3*WV(11,I2)*WV(5,I2)/DFLOAT(LA2+1)*XT                        MTCH-254
      IM=6                                                              MTCH-255
      GO TO 11                                                          MTCH-256
   28 LA2=LA2+1                                                         MTCH-257
      YT=-AY3*DSQRT(1.D0+(WV(5,I2)/DFLOAT(LA2))**2)*WV(11,I2)*XT        MTCH-258
      ZT=0.D0                                                           MTCH-259
      IM=7                                                              MTCH-260
      GO TO 11                                                          MTCH-261
   29 LA1=LA1+1                                                         MTCH-262
      LA2=LA2-1                                                         MTCH-263
      ZT=-AX3*DSQRT(1.D0+(WV(5,I1)/DFLOAT(LA1))**2)*WV(11,I1)*XT        MTCH-264
      YT=0.D0                                                           MTCH-265
      IM=8                                                              MTCH-266
      GO TO 11                                                          MTCH-267
   30 I=I+1                                                             MTCH-268
      IF (I.LE.INTV) GO TO ( 8 , 17 , 22 ),IML                          MTCH-269
      IF (IML.EQ.1) GO TO 16                                            MTCH-270
   31 IF (IS.EQ.0) GO TO 44                                             MTCH-271
      DO 43 IK=1,IS                                                     MTCH-272
      IF (DABS(BA(1,IK))+DABS(BA(2,IK)).LT.1.D-8) GO TO 43              MTCH-273
      LI=LA(1,IK)                                                       MTCH-274
      LF=LA(2,IK)                                                       MTCH-275
      LL=LA(3,IK)                                                       MTCH-276
      EI=WV(5,I1)                                                       MTCH-277
      EF=WV(5,I2)                                                       MTCH-278
      XI=WV(11,I1)*AZ                                                   MTCH-279
      XF=WV(11,I2)*AZ                                                   MTCH-280
      L3=(LI+LF-LL+3)/2                                                 MTCH-281
      IF (L3.LE.0) GO TO 43                                             MTCH-282
      CALL CORA(LL,LI,LF,EI,EF,XI,XF,B,C,LO(127))                       MTCH-283
C ORDER IN C   (LI,LF),(LI,LF+1),(LI+1,LF),(LI+1,LF+1).                 MTCH-284
C ORDER IN FG   F(EI)*F(EF),G(EI)*F(EF),F(EI)*G(EF),G(EI)*G(EF).        MTCH-285
      IF (LO(127)) GO TO 40                                             MTCH-286
      IF (L3+3.GT.LMAX3) GO TO 64                                       MTCH-287
C INTEGRATION REGION - INTEGRALS FROM MATCHING POINT TO INFINITY.       MTCH-288
C COMPUTATION OF THE INTEGRALS FROM THE MATCHING POINT USING B AND C.   MTCH-289
      DO 33 I=1,4                                                       MTCH-290
      G(I)=0.D0                                                         MTCH-291
      DO 32 J=1,4                                                       MTCH-292
   32 G(I)=G(I)+B(J)*FG(L3+J-1,I,I3)-C(J)*FG(L3+N3(J),N2(I),I1)*FG(L3+N4MTCH-293
     1(J),N1(I),I2)                                                     MTCH-294
   33 CONTINUE                                                          MTCH-295
      DO 35 I=1,ILM                                                     MTCH-296
      DO 34 J=1,4                                                       MTCH-297
   34 AB(J,I)=AB(J,I)+BA(I,IK)*G(J)                                     MTCH-298
   35 CONTINUE                                                          MTCH-299
      IF (LO(125).OR.(LO(29).AND.(IP.EQ.IC))) GO TO 38                  MTCH-300
      C0=-DFLOAT(LL+1)/AZ+6.D0/WV(8,I1)                                 MTCH-301
      DO 37 I=1,ILM                                                     MTCH-302
      YT=BA(I,IK)*WV(8,I1)**2*AY/(12.D0*AW**(LL+1))                     MTCH-303
C WITH THE GREEN'S FUNCTIONS METHOD,WE MUST ADD H**2*FP(R)/12           MTCH-304
C WHICH IS VRE(ISM-1)*(KI*FP(I)*G(F)+KF*F(I)*GP(F)-(LL+1)/R)/12.        MTCH-305
C FINITE STEP CORRECTIONS.                                              MTCH-306
      DO 36 J=1,4                                                       MTCH-307
   36 AB(J,I)=AB(J,I)-YT*(WV(11,I1)*FG(LI+1,N2(J)+1,I1)*FG(LF+1,N1(J),I2MTCH-308
     1)+WV(11,I2)*FG(LI+1,N2(J),I1)*FG(LF+1,N1(J)+1,I2)+C0*FG(LI+1,N2(J)MTCH-309
     2,I1)*FG(LF+1,N1(J),I2))                                           MTCH-310
   37 CONTINUE                                                          MTCH-311
      GO TO 43                                                          MTCH-312
C CORRECTION OF THE POTENTIAL IN FIVE POINTS FOR MATCHING VALUES.       MTCH-313
   38 C1=DFLOAT(ISM-1)*WV(8,I1)*AY                                      MTCH-314
      DO 39 I=1,5                                                       MTCH-315
      AV(I)=AV(I)+BA(1,IK)*AY*WV(8,I1)**2/C1**(LL+1)/48.D0              MTCH-316
   39 C1=C1+0.5D0*WV(8,I1)*AY                                           MTCH-317
      GO TO 43                                                          MTCH-318
C COMPUTATION OF INTEGRALS FROM 0 TO INFINITY USING B.                  MTCH-319
   40 AX=0.D0                                                           MTCH-320
      IF (L3+3.GT.LMAX2) GO TO 65                                       MTCH-321
      DO 41 I=1,4                                                       MTCH-322
   41 AX=AX+B(I)*BG(L3+I-1,I3)                                          MTCH-323
      DO 42 I=1,ILM                                                     MTCH-324
   42 AB(1,I)=AB(1,I)+BA(I,IK)*AX                                       MTCH-325
   43 CONTINUE                                                          MTCH-326
C LIMITATION FOR TOO LARGE INTEGRAL OF IRREGULAR FUNCTIONS.             MTCH-327
   44 IF (DABS(AB(4,1)).GT.WV(11,I1)) AB(4,1)=AB(1,1)                   MTCH-328
      IF ((ILM.EQ.2).AND.(DABS(AB(4,2)).GT.WV(11,I2))) AB(4,2)=AB(1,2)  MTCH-329
      IF (.NOT.LV) GO TO 45                                             MTCH-330
      WRITE (MW,1001) IC,IP,L1,L2,(AB(I,1),I=1,LL4)                     MTCH-331
      IF (ILM.EQ.2) WRITE (MW,1001) IP,IC,L2,L1,(AB(I,2),I=1,LL4)       MTCH-332
C BUILD UP OF MATRIX OF CORRECTIONS.                                    MTCH-333
   45 DO 46 L=1,LL4                                                     MTCH-334
      IF (IC.NE.IP) AA(IP,IC,L)=AB(L,ILM)                               MTCH-335
   46 AA(IC,IP,L)=AB(L,1)                                               MTCH-336
      IF (LO(127).OR.(IC.EQ.IP)) GO TO 47                               MTCH-337
      AA(IP,IC,2)=AB(3,ILM)                                             MTCH-338
      AA(IP,IC,3)=AB(2,ILM)                                             MTCH-339
   47 CONTINUE                                                          MTCH-340
      IF (.NOT.LO(127)) GO TO 48                                        MTCH-341
      IF (JC.LE.0) GO TO 55                                             MTCH-342
      IF (WV(5,I1).EQ.0.D0) F2=F3*WV(11,I1)                             MTCH-343
      AA(IC,IC,1)=AA(IC,IC,1)-F2*BG(L1+1,JC)*WV(11,I1)                  MTCH-344
      IF ((WV(5,I1).EQ.0.D0).OR.(F3.EQ.0.D0)) GO TO 55                  MTCH-345
      B1=2.D0*WV(5,I1)*DFLOAT(L1*(L1+1))                                MTCH-346
      B2=DFLOAT(L1+1)**2+WV(5,I1)**2                                    MTCH-347
      C1=-(DFLOAT(2*L1+1)*B2+2*WV(5,I1)**2)/B1                          MTCH-348
      C2=DFLOAT(2*L1+3)*B2/B1                                           MTCH-349
      AA(IC,IC,1)=AA(IC,IC,1)-F3*(C1*BG(L1+1,JC)+C2*BG(L1+2,JC))*WV(11,IMTCH-350
     11)**2                                                             MTCH-351
      GO TO 55                                                          MTCH-352
   48 IF (MC(IC,4).LT.0) GO TO 55                                       MTCH-353
C MATCHING VALUES.                                                      MTCH-354
      DO 49 K=1,4                                                       MTCH-355
   49 B(K)=FG(L1+1,K,I1)                                                MTCH-356
      IF (JC.LE.0) GO TO 53                                             MTCH-357
      IF ((.NOT.LO(29)).AND.(F2.EQ.0.D0.AND.F3.EQ.0.D0)) GO TO 53       MTCH-358
      IF (WV(5,I1).NE.0.D0) GO TO 50                                    MTCH-359
      F2=F3*WV(11,I1)                                                   MTCH-360
      F3=0.D0                                                           MTCH-361
   50 DO 51 K=1,4                                                       MTCH-362
   51 G(K)=-FG(L1+1,K,JC)*F2+AB(K,1)/WV(11,I1)                          MTCH-363
      IF (F3.EQ.0.D0) GO TO 52                                          MTCH-364
      B1=2.D0*WV(5,I1)*DFLOAT(L1*(L1+1))                                MTCH-365
      B2=DFLOAT(L1+1)**2+WV(5,I1)**2                                    MTCH-366
      C1=-(DFLOAT(2*L1+1)*B2+2*WV(5,I1)**2)/B1                          MTCH-367
      C2=DFLOAT(2*L1+3)*B2/B1                                           MTCH-368
      A1=DFLOAT(ISM)*WV(8,I1)*WV(11,I1)                                 MTCH-369
      D1=(B2+DFLOAT(L1+1)*WV(5,I1)/A1)/A1/B1                            MTCH-370
      D2=-WV(5,I1)*DSQRT(B2)/B1/A1                                      MTCH-371
      A1=B2/B1/A1                                                       MTCH-372
      A3=F3*WV(11,I1)                                                   MTCH-373
      G(1)=G(1)-A3*(C1*FG(L1+1,1,JC)+C2*FG(L1+2,1,JC)-D1*FG(L1+1,1,I1)**MTCH-374
     12-D2*2.D0*FG(L1+1,1,I1)*FG(L1+2,1,I1)-A1*FG(L1+2,1,I1)**2)        MTCH-375
      G(2)=G(2)-A3*(C1*FG(L1+1,2,JC)+C2*FG(L1+2,2,JC)-D1*FG(L1+1,1,I1)*FMTCH-376
     1G(L1+1,3,I1)-D2*(FG(L1+1,1,I1)*FG(L1+2,3,I1)+FG(L1+2,1,I1)*FG(L1+1MTCH-377
     2,3,I1))-A1*FG(L1+2,1,I1)*FG(L1+2,3,I1))                           MTCH-378
      G(4)=G(4)-A3*(C1*FG(L1+1,4,JC)+C2*FG(L1+2,4,JC)-D1*FG(L1+1,3,I1)**MTCH-379
     12-D2*2.D0*FG(L1+1,3,I1)*FG(L1+2,3,I1)-A1*FG(L1+2,3,I1)**2)        MTCH-380
   52 A4=1.D0+(G(1)*G(4)-G(2)**2)                                       MTCH-381
      G(3)=B(1)                                                         MTCH-382
      B(1)=(B(1)*(1.D0-G(2))+G(1)*B(3))/A4                              MTCH-383
      B(3)=(-G(3)*G(4)+(1.D0+G(2))*B(3))/A4                             MTCH-384
      G(3)=B(2)                                                         MTCH-385
      B(2)=(B(2)*(1.D0-G(2))+G(1)*B(4))/A4                              MTCH-386
      B(4)=(-G(3)*G(4)+(1.D0+G(2))*B(4))/A4                             MTCH-387
   53 A1=(1.D0-AV(2))/(2.D0+10.D0*AV(2))                                MTCH-388
      B1=(1.D0-AV(4))/(2.D0+10.D0*AV(4))                                MTCH-389
      A2=A1*(1.D0-AV(1))/(1.D0-4.D0*AV(1))                              MTCH-390
      B2=B1*(1.D0-AV(5))/(1.D0-4.D0*AV(5))                              MTCH-391
      C1=(2.D0+10.D0*AV(3))-(1.D0-AV(3))*(A1+B1)                        MTCH-392
      A1=(16.D0-144.D0*AV(2))/(2.D0+10.D0*AV(2))                        MTCH-393
      B1=(16.D0-144.D0*AV(4))/(2.D0+10.D0*AV(4))                        MTCH-394
      C2=(7.D0+A1*(1.D0-AV(1)))/(1.D0-4.D0*AV(1))                       MTCH-395
      D2=(7.D0+B1*(1.D0-AV(5)))/(1.D0-4.D0*AV(5))                       MTCH-396
      D1=(B1-A1)*(1.D0-AV(3))                                           MTCH-397
      A1=A2*D2+B2*C2                                                    MTCH-398
      B1=(C1*D2+D1*B2)/A1                                               MTCH-399
      B2=30.D0*WV(8,1)*B2*WV(11,I1)/A1                                  MTCH-400
      FAM(IC,1)=B1*B(1)-B2*B(2)                                         MTCH-401
      FAM(IC,3)=B1*B(3)-B2*B(4)                                         MTCH-402
      B1=(C2*C1-A2*D1)/A1                                               MTCH-403
      A2=-30.D0*WV(8,1)*A2*WV(11,I1)/A1                                 MTCH-404
      FAM(IC,2)=B1*B(1)-A2*B(2)                                         MTCH-405
      FAM(IC,4)=B1*B(3)-A2*B(4)                                         MTCH-406
      FAM(IC,5)=WV(11,I1)                                               MTCH-407
      IF (LO(100)) FAM(IC,5)=FAM(IC,5)*CHB/(WV(7,I1)+WV(6,I1))          MTCH-408
      BT=FAM(IC,2)*FAM(IC,3)-FAM(IC,1)*FAM(IC,4)                        MTCH-409
      IF (BT.EQ.0.D0) BT=1.D0                                           MTCH-410
      DO 54 I=1,4                                                       MTCH-411
   54 FAM(IC,I)=FAM(IC,I)/BT                                            MTCH-412
      FAM(IC,6)=FAM(IC,5)/BT                                            MTCH-413
      IF (LO(100)) FAM(IC,5)=FAM(IC,5)*WV(8,I1)                         MTCH-414
      FAM(IC,7)=WV(9,I1)**2                                             MTCH-415
      FAM(IC,9)=2.D0*FAM(IC,7)*DFLOAT(MC(IC,6))                         MTCH-416
      FAM(IC,8)=WV(12,I1)                                               MTCH-417
      FAM(IC,10)=DFLOAT(MC(IC,5))                                       MTCH-418
      FAM(IC,11)=WV(10,I1)**2                                           MTCH-419
      FAM(IC,12)=2.D0*FAM(IC,11)*DFLOAT(MC(IC,6))                       MTCH-420
      IF (.NOT.LO(103)) FAM(IC,12)=0.D0                                 MTCH-421
      IF (WV(5,I1).EQ.0.D0) FAM(IC,11)=0.D0                             MTCH-422
   55 CONTINUE                                                          MTCH-423
      DO 58 IC=1,NC                                                     MTCH-424
      I1=MC(IC,1)                                                       MTCH-425
      DO 57 L=1,LL4                                                     MTCH-426
      DO 56 IP=1,NC                                                     MTCH-427
   56 AA(IC,IP,L)=AA(IC,IP,L)/WV(11,I1)                                 MTCH-428
   57 CONTINUE                                                          MTCH-429
   58 CONTINUE                                                          MTCH-430
      IF (.NOT.LO(127)) RETURN                                          MTCH-431
      IF (.NOT.LO(125)) GO TO 61                                        MTCH-432
C TRANSPOSITION IF COUPLED EQUATIONS ARE USED.                          MTCH-433
      DO 60 IC=1,NC                                                     MTCH-434
      DO 59 IP=1,IC                                                     MTCH-435
      AX=AA(IP,IC,1)                                                    MTCH-436
      AA(IP,IC,1)=AA(IC,IP,1)                                           MTCH-437
   59 AA(IC,IP,1)=AX                                                    MTCH-438
   60 CONTINUE                                                          MTCH-439
C LINEAR SYSTEM FOR COMPUTATION OF THE S-MATRIX FROM THE K-MATRIX.      MTCH-440
   61 DO 63 IC=1,NC                                                     MTCH-441
      DO 62 IP=1,NC                                                     MTCH-442
      AA(IC,IP,2)=0.D0                                                  MTCH-443
      AA(IC,IP,4)=0.D0                                                  MTCH-444
   62 AA(IC,IP,5)=-AA(IC,IP,1)                                          MTCH-445
   63 AA(IC,IC,4)=1.D0                                                  MTCH-446
      CALL LINS(AA(1,1,4),KAB,AA,KAB,AA(1,1,5),KAB,AA(1,1,2),KAB,NC,NCXNMTCH-447
     1,KR,IER)                                                          MTCH-448
      RETURN                                                            MTCH-449
   64 WRITE (MW,1002) L3,LI,LF,LMAX3                                    MTCH-450
      RETURN                                                            MTCH-451
   65 WRITE (MW,1003) L3,LMAX2                                          MTCH-452
      RETURN                                                            MTCH-453
 1000 FORMAT (/' CHANNELS  L-VALUES     COULOMB INTEGRALS: F*F,G*F,F*G AMTCH-454
     1ND G*G:')                                                         MTCH-455
 1001 FORMAT (2X,2I3,2X,2I4,3X,4D18.10)                                 MTCH-456
 1002 FORMAT (' STARTING VALUES',I3,' FOR INTEGRALS AND',2I4,' FOR COULOMTCH-457
     1MB FUNCTIONS TOO LARGE . LIMITATION AT',I4)                       MTCH-458
 1003 FORMAT (' STARTING VALUE',I4,' TOO LARGE FOR COULOMB INTEGRALS WHIMTCH-459
     1CH ARE COMPUTED UP TO',I5)                                        MTCH-460
      END                                                               MTCH-461
