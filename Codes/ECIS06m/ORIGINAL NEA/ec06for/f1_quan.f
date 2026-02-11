C 20/08/06                                                      ECIS06  QUAN-000
      SUBROUTINE QUAN(NCOLL,WV,IPI,NIV,IQ,TQ,IVQ,IVZ,MC,NAT,AT,IT,NVI,KAQUAN-001
     1B,KBC,AA,FAC,NFA,NMAX,JD,LO)                                      QUAN-002
C INPUT:     NCOLL:   NUMBER OF COUPLED NUCLEAR STATES.                 QUAN-003
C            WV:      DESCRIPTION OF THE CHANNELS. SEE CALX.            QUAN-004
C            IPI:     PROJECTILE AND TARGET MULTIPLICITY, PARITY        QUAN-005
C                     AND MAXIMUM ANGULAR MOMENTUM.  SEE CALX.          QUAN-006
C            NIV:     FIRST/LAST ADDRESS IN THE TABLE OF                QUAN-007
C                     REDUCED NUCLEAR MATRIX ELEMENTS.  SEE REDM.       QUAN-008
C            IQ,TQ:   TABLE OF REDUCED NUCLEAR MATRIX ELEMENTS.         QUAN-009
C            IVQ:     TABLE OF MULTIPOLES.  SEE REDM.                   QUAN-010
C            IVZ:     TABLE OF FORM FACTORS.  SEE REDM.                 QUAN-011
C            KAB:     MAXIMUM NUMBER OF COUPLED CHANNELS.               QUAN-012
C            KBC:     MAXIMUM NUMBER OF SOLUTIONS.                      QUAN-013
C            FAC:     TABLE OF LOGARITHMS OF FACTORIALS.                QUAN-014
C            NFA:     LENGTH OF FAC.                                    QUAN-015
C            NMAX:    AVAILABLE LENGTH OF AT/NAT(LESS 100 FOR 9J).      QUAN-016
C            JD:      FIRST DIMENSION OF TABLES NAT AND AT.             QUAN-017
C            LO(I):   LOGICAL CONTROLS:                                 QUAN-018
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            QUAN-019
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      QUAN-020
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. QUAN-021
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   QUAN-022
C               LO(58) =.TRUE. OUTPUT OF THE COEFFICIENTS OF EACH FORM  QUAN-023
C                              FACTOR FOR ALL SETS OF EQUATIONS.        QUAN-024
C               LO(100)=.TRUE. DIRAC EQUATION.                          QUAN-025
C               LO(101)=.TRUE. THERE IS A REAL SPIN-ORBIT POTENTIAL.    QUAN-026
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. QUAN-027
C               LO(110)=.TRUE. DERIVATIVES ARE NEEDED.                  QUAN-028
C               LO(123)=.TRUE. IDENTICAL PARTICLES WITH SPIN.           QUAN-029
C               LO(127)=.TRUE. COULOMB CORRECTIONS IN ASYMPTOTIC REGION.QUAN-030
C               LO(128)=.TRUE. NO COPY OF UNCOUPLED FUNCTIONS AND       QUAN-031
C                              PHASE-SHIFTS.                            QUAN-032
C OUTPUT:    MC:      NUCLEAR STATE INDICATIONS. FOR IC=1,NC:           QUAN-033
C                     MC(IC,1): NUCLEAR STATE NUMBERS.                  QUAN-034
C                     MC(IC,2): ORBITAL MOMENTUM LC.                    QUAN-035
C                     MC(IC,3): TWICE ORBITAL SPIN.                     QUAN-036
C                     MC(IC,4): REFERENCE TO POTENTIAL OR COPY.         QUAN-037
C                     MC(IC,5): LC*(LC+1).                              QUAN-038
C                     MC(IC,6): EIGENVALUE OF L.S.                      QUAN-039
C            NAT,AT:  IN EQUIVALENCE BY CALL. TABLE OF COUPLING         QUAN-040
C                     COEFFICIENTS. FOR I=1,IT:                         QUAN-041
C                     NAT(1,I): ADDRESS OF THE REAL FORM FACTOR,        QUAN-042
C                     NAT(2,I): ADDRESS OF THE IMAGINARY FORM FACTOR,   QUAN-043
C                     AT(2,I):  GEOMETRICAL COEFFICIENT,                QUAN-044
C                     AND ONLY FOR DIRAC EQUATION:                      QUAN-045
C                     AT(3,I):  IDEM FOR SMALL COMPONENTS,              QUAN-046
C            IT:      NUMBER OF NON ZERO COUPLING COEFFICIENTS.         QUAN-047
C            NVI:     ADDRESS IN THE TABLE AT(I), ANALOGOUS TO THE NIV  QUAN-048
C                     ADDRESSES: NVI(I1,I2,1) TO NVI(I1,I2,2)           QUAN-049
C                     FOR THE NON DERIVATIVES COUPLINGS, NVI(I1,I2,2)+1 QUAN-050
C                     TO NVI(I1,I2,3)  FOR THE DERIVATIVE COUPLINGS.    QUAN-051
C                     FOR SMALL COMPONENTS OF DIRAC EQUATION, IDEM WITH QUAN-052
C                     NVI(I1,I2,4) TO NVI(I1,I2,6).                     QUAN-053
C            AA:      COEFFICIENTS OF THE SYMMETRISATION FOR IDENTICAL  QUAN-054
C                     PARTICLE AND TARGET, MC RESULTS ARE IN MC(*,*+7), QUAN-055
C                     NAT AND AT RESULTS ARE AFTER THE USUAL ONES, NVI  QUAN-056
C                     RESULTS ARE IN NVI(*,*,*+3), NC AND NCIN ARE      QUAN-057
C                     DIFFERENT FROM NIC AND NCI.                       QUAN-058
C            LO(I):   LOGICAL CONTROLS: DEFINED HERE: LO(110).          QUAN-059
C WORKING AREA:                                                         QUAN-060
C            AT(I,J): FOR J>IT, WORKING AREA USED IN DJ9J.              QUAN-061
C                                                                       QUAN-062
C THE COMMON /NOEQU/ IS USED IN CAL1, QUAN, MTCH AND SCAM.              QUAN-063
C                                                                       QUAN-064
C FOR THE COMMON  /NCOMP/ SEE CALX.                                     QUAN-065
C FOR THE COMMON  /POTE2/ SEE REDM.                                     QUAN-066
C                                                                       QUAN-067
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     QUAN-068
C  AZ(6):     DEFORMED SPIN-ORBIT PARAMETERS. (SEE INPUT DESCRIPTION    QUAN-069
C             AND COMMENTS IN THIS SUBROUTINE).                         QUAN-070
C                                                                       QUAN-071
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NOEQU/:                     QUAN-072
C  NCXN:      NUMBER OF SOLUTIONS NEEDED.                               QUAN-073
C  NIC:       NUMBER OF EQUATIONS AT THE CHANNEL SPIN.                  QUAN-074
C  NCI:       NUMBER OF SOLUTIONS AT THE CHANNEL SPIN.                  QUAN-075
C  NC:        NUMBER OF EQUATIONS FOR IDENTICAL PARTICLES.              QUAN-076
C  NCIN:      NUMBER OF SOLUTIONS FOR IDENTICAL PARTICLES.              QUAN-077
C  NIN:       NUMBER OF COUPLING POTENTIALS.                            QUAN-078
C  JPI:       PARITY 0 OR 1.                                            QUAN-079
C  IPJ:       VALUE OF J+1 OR J+1/2 WHERE J IS THE CHANNEL SPIN.        QUAN-080
C  R1(2):     MAXIMUM OF SCATTERING AND COMPOUND COEFFICIENT.           QUAN-081
C  NAJ:       TWICE THE CHANNEL SPIN J.                                 QUAN-082
C   DEFINED:  NIC,NCI,NC,NCIN,NIN.                                      QUAN-083
C   USED:     NIC,NCI,NC,NCIN,NIN,JPI,NAJ.                              QUAN-084
C   NOT USED: NCXN,IPJ,R1.                                              QUAN-085
C                                                                       QUAN-086
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     QUAN-087
C  ITY(5):    STARTING ADDRESS OF REAL CENTRAL TRANSITION.              QUAN-088
C             ITY(5)=0 IS USED FOR DIRAC EQUATIONS.                     QUAN-089
C  INTV:      NUMBER OF TRANSITION FORM FACTORS WITHOUT SPIN-ORBIT,     QUAN-090
C             TAKING INTO ACCOUNT DISPERSION.                           QUAN-091
C  INSL:      NUMBER OF SPIN-ORBIT FORM TRANSITION FACTORS NOT TAKING   QUAN-092
C             INTO ACCOUNT MULTIPLICATION BY 2.                         QUAN-093
C   USED:     ITY(5),INTV,INSL.                                         QUAN-094
C                                                                       QUAN-095
C***********************************************************************QUAN-096
      IMPLICIT REAL*8 (A-H,O-Z)                                         QUAN-097
      LOGICAL LO(150),LLO,LLC,LLP                                       QUAN-098
      DIMENSION WV(22,*),IPI(11,*),NIV(NCOLL,NCOLL,3),IQ(6,*),TQ(3,*),IVQUAN-099
     1Q(3,*),IVZ(4,*),MC(KAB,12),NAT(2*JD,*),AT(JD,*),NVI(KAB,KAB,6),AA(QUAN-100
     2KBC,*),FAC(*)                                                     QUAN-101
      CHARACTER*1 IP(2)                                                 QUAN-102
      COMMON /INOUT/ MR,MW,MS                                           QUAN-103
      COMMON /NCOMP/ NSP(12),ACN(2),AZ(6),BZ(12)                        QUAN-104
      COMMON /NOEQU/ NCXN,NIC,NCI,NC,NCIN,NIN,JPI,IPJ,R1(2),NAJ         QUAN-105
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         QUAN-106
      DATA IP /'+','-'/                                                 QUAN-107
      LLC=LO(133).AND.LO(11)                                            QUAN-108
      LLP=LO(133).AND.LO(19)                                            QUAN-109
C FIND QUANTUM NUMBER OF COUPLED CHANNELS.                              QUAN-110
      NC=0                                                              QUAN-111
      AJ=.5D0*DFLOAT(NAJ)                                               QUAN-112
      DO 6 I=1,NCOLL                                                    QUAN-113
      NJ1=NAJ-IPI(3,I)+1                                                QUAN-114
      NJ2=IABS(NJ1)                                                     QUAN-115
      NJ=IPI(3,I)                                                       QUAN-116
      NL=IPI(2,I)                                                       QUAN-117
      DO 5 J=1,NJ                                                       QUAN-118
      L1=(NAJ-NJ-NL)/2+J                                                QUAN-119
      L2=IABS(L1)                                                       QUAN-120
      DO 4 K=1,NL                                                       QUAN-121
      IF (MOD(L1+IPI(1,I)+JPI,2).NE.0) GO TO 3                          QUAN-122
      IF ((L1.LT.L2).OR.(NJ1.LT.NJ2).OR.(L1.GT.IPI(10,I))) GO TO 3      QUAN-123
      NC=NC+1                                                           QUAN-124
      MC(NC,4)=IPI(11,I)                                                QUAN-125
      IF (LO(128)) GO TO 2                                              QUAN-126
      IF (LO(100)) MC(NC,4)=I                                           QUAN-127
      DO 1 M=NC,KAB                                                     QUAN-128
      IF (MC(M,1).NE.I) GO TO 1                                         QUAN-129
      IF (MC(M,2).NE.L1) GO TO 1                                        QUAN-130
      IF (MC(M,3).NE.NJ1) GO TO 1                                       QUAN-131
      MC(NC,4)=-M                                                       QUAN-132
      GO TO 2                                                           QUAN-133
    1 CONTINUE                                                          QUAN-134
    2 MC(NC,1)=I                                                        QUAN-135
      MC(NC,2)=L1                                                       QUAN-136
      MC(NC,3)=NJ1                                                      QUAN-137
      MC(NC,5)=MC(NC,2)*(MC(NC,2)+1)                                    QUAN-138
      MC(NC,6)=((NJ1-IPI(2,I)+1)*(NJ1+IPI(2,I)+1))/4-MC(NC,5)           QUAN-139
    3 L1=L1+1                                                           QUAN-140
    4 CONTINUE                                                          QUAN-141
    5 NJ1=NJ1+2                                                         QUAN-142
      IF (I.NE.1) GO TO 6                                               QUAN-143
      NCIN=NC                                                           QUAN-144
      IF (NCIN.EQ.0) RETURN                                             QUAN-145
    6 CONTINUE                                                          QUAN-146
      IF (.NOT.LO(58)) GO TO 8                                          QUAN-147
      WRITE (MW,1000) AJ,IP(JPI+1),NC,NCIN                              QUAN-148
      DO 7 J=1,NC                                                       QUAN-149
      SJ=0.5D0*DFLOAT(MC(J,3))                                          QUAN-150
    7 WRITE (MW,1001) J,MC(J,1),MC(J,2),SJ,MC(J,4),MC(J,5),MC(J,6)      QUAN-151
C COMPUTATION OF COUPLING COEFFICIENTS FOR THE FORM FACTORS             QUAN-152
C WHEN THE SPIN-ORBIT IS DEFORMED,THERE IS NO SYMMETRY, THE TOTAL       QUAN-153
C TABLE IS CALCULATED. WITH NO SPIN-ORBIT DEFORMATION,ONLY ONE HALF IS  QUAN-154
C CALCULATED. LO(110) RETURNS .TRUE. IF THEY ARE DERIVATIVE COUPLINGS.  QUAN-155
    8 LO(110)=.FALSE.                                                   QUAN-156
      IT=0                                                              QUAN-157
      NIN=0                                                             QUAN-158
      INVZ=INTV+ITY(5)                                                  QUAN-159
      DO 41 I1=1,NC                                                     QUAN-160
      J1=MC(I1,1)                                                       QUAN-161
      K1=IPI(11,J1)                                                     QUAN-162
      IF (LO(100)) K1=J1                                                QUAN-163
      IF (LO(127)) MC(I1,4)=K1                                          QUAN-164
      DO 40 I2=1,I1                                                     QUAN-165
      J2=MC(I2,1)                                                       QUAN-166
      A2=WV(9,J1)*WV(9,J2)                                              QUAN-167
      A3=WV(10,J1)*WV(10,J2)                                            QUAN-168
      IT1=IT+1                                                          QUAN-169
      NVI(I2,I1,1)=IT1                                                  QUAN-170
      K1=NIV(J2,J1,1)                                                   QUAN-171
      K2=NIV(J2,J1,2)                                                   QUAN-172
      IF (K1.GT.K2) GO TO 32                                            QUAN-173
      IF (2*(K2-K1+3+IT).GT.NMAX) CALL MEMO('QUAN',NMAX,2*(K2-K1+3+IT)) QUAN-174
      IL=2*MC(I2,2)                                                     QUAN-175
      ILP=2*MC(I1,2)                                                    QUAN-176
      ISJ=MC(I2,3)                                                      QUAN-177
      ISJP=MC(I1,3)                                                     QUAN-178
      IAI=IPI(3,J2)-1                                                   QUAN-179
      IAIP=IPI(3,J1)-1                                                  QUAN-180
      IS=IPI(2,J2)-1                                                    QUAN-181
      ISP=IPI(2,J1)-1                                                   QUAN-182
      DO 22 K=K1,K2                                                     QUAN-183
      NI=IQ(2,K)                                                        QUAN-184
C COEFFICIENT OF A CENTRAL MULTIPOLE WITH ANGULAR MOMENTUM IQ=IVQ(1,NI) QUAN-185
C TRANSFER OF SPIN IS=IVQ(2,NI)/2 AND OF TOTAL MOMENTUM IJ=IVQ(3,NI)/2  QUAN-186
C NOTING THE SPIN S,THE ORBITAL ANGULAR MOMENTUM L,THE TOTAL SPIN OF THEQUAN-187
C PARTICLE J,THE SPIN OF THE TARGET AI,THE EIGENVALUE OF L.S G......    QUAN-188
C AND "K"=SQRT(2*K+1)   THE MOST GENERAL EXPRESSION OF (I2||IQ||I1) IS  QUAN-189
C  (-)**(AJ+AI2+1+J1+(L2+IQ-L1)/2)*"L1"*"L2"*"J1"*"J2"*"IQ"*"IJ"*       QUAN-190
C CGS(IQ,L2,L1)*C6J(J1,IJ,J2,AI2,AJ,AI1)*C9J(L2,L1,IQ,S2,S1,IS,J2,J1,IJ)QUAN-191
C IN A MACROSCOPIC MODEL, IS=0 AND IJ=IQ, THE EXPRESSIONS USED ARE      QUAN-192
C  (-)**(AJ+AI2+1-S+(L1+L2+IQ)/2)*"L1"*"L2"*"J1"*"J2"*"IQ"*CGS(IQ,L2,L1)QUAN-193
C    *C6J(J1,IQ,J2,AI2,AJ,AI1)*C6J(L1,IQ,L2,J2,S,J1)                    QUAN-194
C FOR S=1/2 (-)**(AJ+AI2-J2+(L1-L2+IQ)/2)*"IQ"*CGS(IQ,J2,J1)            QUAN-195
C FOR S=0   (-)**(AJ+AI2+1+(L1+L2+IQ)/2)*"L1"*"L2"*"IQ"*CGS(IQ,L2,L1)*  QUAN-196
C           C6J(L1,IQ,L2,AI2,AJ,AI1).                                   QUAN-197
      IIQ=2*IVQ(1,NI)                                                   QUAN-198
      DO 16 IX=2,JD                                                     QUAN-199
      CX=0.D0                                                           QUAN-200
      CMATEL=0.D0                                                       QUAN-201
      IF ((IIQ.GT.IL+ILP).OR.(IIQ.LT.IABS(IL-ILP))) GO TO 13            QUAN-202
      FX=1.D0                                                           QUAN-203
      KK=MC(I1,2)+MC(I2,2)+IVQ(1,NI)                                    QUAN-204
      IIS=IVQ(2,NI)                                                     QUAN-205
      IF (IIS.LT.0) GO TO 11                                            QUAN-206
      IIJ=IVQ(3,NI)                                                     QUAN-207
      CMATEL=DJ6J(ISJ,IIJ,ISJP,IAIP,NAJ,IAI,FAC,NFA)                    QUAN-208
      A1=DFLOAT(IIJ+1)                                                  QUAN-209
      IF (IS.GT.1.OR.IIS.NE.0) GO TO 9                                  QUAN-210
      CMATEL=CMATEL*DCGS(IIJ,ISJP,ISJ,FAC,NFA)                          QUAN-211
      KK=KK+NAJ+IAIP                                                    QUAN-212
      IF (IS.EQ.1) KK=KK+2+IL-ISJP                                      QUAN-213
      IF (IS.NE.1) A1=A1*DFLOAT((ISJ+1)*(ISJP+1))                       QUAN-214
      GO TO 12                                                          QUAN-215
    9 A1=A1*DFLOAT((IL+1)*(ILP+1))                                      QUAN-216
      A1=A1*DFLOAT((ISJ+1)*(ISJP+1))                                    QUAN-217
      KK=(NAJ+IAIP-IS)+KK                                               QUAN-218
      CMATEL=CMATEL*DCGS(IIQ,ILP,IL,FAC,NFA)                            QUAN-219
      IF (IIS.NE.0) GO TO 10                                            QUAN-220
      FX=DJ6J(IL,ISJ,IS,ISJP,ILP,IIJ,FAC,NFA)                           QUAN-221
      GO TO 12                                                          QUAN-222
   10 CMATEL=CMATEL*DJ9J(IL,ILP,IIQ,IS,ISP,IIS,ISJ,ISJP,IIJ,AT(1,IT+1),1QUAN-223
     100)                                                               QUAN-224
      A1=A1*DFLOAT((IIQ+1)*(IIS+1))                                     QUAN-225
      KK=KK+ISJP+ILP+IS                                                 QUAN-226
      GO TO 12                                                          QUAN-227
C MAGNETIC COULOMB EXCITATION OF THE PARTICLE:                          QUAN-228
   11 IIJ=IIQ-2                                                         QUAN-229
      IF (IIQ.GT.IL+ILP+2.OR.IIQ.LT.IABS(IL-ILP)+2) GO TO 13            QUAN-230
      CMATEL=DJ6J(ISJ,ISJP,IIJ,ILP,IL,IS,FAC,NFA)*DJ6J(ISJ,ISJP,IIQ-2,IAQUAN-231
     1IP,IAI,NAJ,FAC,NFA)*DCGS(IIQ,ILP,IL,FAC,NFA)                      QUAN-232
      KK=KK+NAJ+IAIP-IS-2                                               QUAN-233
      A1=DFLOAT((IIQ+1)**3)*DFLOAT((ISJ+1)*(ISJP+1))*DFLOAT((IL+1)*(ILP+QUAN-234
     11))*DFLOAT((IL+ILP+IIQ+2)*(IL+ILP-IIJ))*(IL-ILP+IIQ)*(ILP-IL+IIQ)/QUAN-235
     2DFLOAT((IIQ-1)*IIJ)**2/16.D0                                      QUAN-236
   12 IF (MOD(KK,4).EQ.0) CMATEL=-CMATEL                                QUAN-237
      CMATEL=TQ(3,K)*CMATEL*DSQRT(A1)                                   QUAN-238
      CX=CMATEL*FX                                                      QUAN-239
   13 IF (IX.EQ.3) GO TO 14                                             QUAN-240
      IT=IT+1                                                           QUAN-241
      NAT(1,IT)=IQ(1,K)                                                 QUAN-242
      NAT(2,IT)=0                                                       QUAN-243
   14 AT(IX,IT)=CX*A2                                                   QUAN-244
      IF (.NOT.LO(100)) GO TO 17                                        QUAN-245
      IF (IX.EQ.3) GO TO 22                                             QUAN-246
      IF (IIQ.NE.IIJ) GO TO 15                                          QUAN-247
      AT(3,IT)=AT(2,IT)                                                 QUAN-248
      IF (IIS.NE.0) AT(3,IT)=-AT(3,IT)                                  QUAN-249
      GO TO 22                                                          QUAN-250
   15 IL=2*ISJ-IL                                                       QUAN-251
   16 ILP=2*ISJP-ILP                                                    QUAN-252
   17 NAT(1,IT)=NAT(1,IT)+ITY(5)                                        QUAN-253
      IF (LO(133)) AT(3,IT)=0.D0                                        QUAN-254
      IF (LLC) AT(3,IT)=CX*A3                                           QUAN-255
      IF (LO(12)) NAT(2,IT)=NAT(1,IT)                                   QUAN-256
      IF (IQ(3,K).EQ.0) GO TO 22                                        QUAN-257
C DEFORMED SPIN-ORBIT:                                                  QUAN-258
C FORM FACTOR (1/R)(D/DR)V(R)         COEFFICIENT:  (I2||IQ||I1)*G1     QUAN-259
C PARAMETRISATION:  (I2||IQ||I1)*(AZ(3)*G1+AZ(4)*G2+AZ(1)).             QUAN-260
      A=CX*(MC(I2,6)*AZ(3)+MC(I1,6)*AZ(4)+AZ(1))                        QUAN-261
      IF (DABS(A).LT.1.D-10) GO TO 18                                   QUAN-262
      IT=IT+1                                                           QUAN-263
      NAT(1,IT)=IQ(3,K)+INVZ                                            QUAN-264
      NAT(2,IT)=0                                                       QUAN-265
      IF (LO(14)) NAT(2,IT)=NAT(1,IT)                                   QUAN-266
      AT(2,IT)=2.D0*A*A2                                                QUAN-267
      IF (LO(133)) AT(3,IT)=0.D0                                        QUAN-268
      IF (LLP) AT(3,IT)=2.D0*A*A3                                       QUAN-269
C FORM FACTOR V(R)/R**2   COEFFICIENT: (I2||IQ||I1)*(S*(IQ*(IQ+1)-L1*(L1QUAN-270
C   +1)-L2*(L2+1)+G2*(1+G1)/(2*S))+XXX)                                 QUAN-271
C PARAMETRISATION:   (I2||IQ||I1)*AZ(5)*(S*(IQ*(IQ+1)*AZ(2)........     QUAN-272
C XXX IS -SQRT(F1*F2)*C6J(L1,J1,S-1,J2,L2,IQ)/(2*S*C6J(L1,J1,S,J2,L2,IQ)QUAN-273
C  WITH F=(L*L+L-(J-S)*(J-S+1))*((J+S)*(J+S+1)-L*L-L)                   QUAN-274
C  FOR S=1/2  XXX=0                                                     QUAN-275
C  FOR S=1  XXX=-4*L1*L2*(L1+1)*(L2+1)/IQ*IQ+IQ-L1*L1-L1-L2*L2-L2)      QUAN-276
C WHEN J1=L1 AND J2=L2 , ELSE  XXX=0.                                   QUAN-277
   18 AZI=DFLOAT(IVQ(1,NI)*(IVQ(1,NI)+1))                               QUAN-278
      ASPI=DFLOAT(IS)                                                   QUAN-279
      SMATEL=(0.5D0*DFLOAT(IPI(2,J1)-1)*(AZ(2)*AZI-DFLOAT(MC(I1,5)+MC(I2QUAN-280
     1,5)))+DFLOAT(MC(I1,6))*(1.D0+DFLOAT(MC(I2,6))/ASPI))*FX           QUAN-281
      IF (IS-2) 21 , 19 , 20                                            QUAN-282
   19 IF (IL.NE.ISJ.OR.ILP.NE.ISJP) GO TO 21                            QUAN-283
      A1=AZI-MC(I1,5)-MC(I2,5)                                          QUAN-284
      IF (A1.EQ.0.D0) GO TO 20                                          QUAN-285
      SMATEL=SMATEL-4.D0*DFLOAT(MC(I1,5)*MC(I2,5))*FX/A1                QUAN-286
      GO TO 21                                                          QUAN-287
   20 G1=DFLOAT((MC(I1,3)-IPI(2,J1))/2)                                 QUAN-288
      G2=DFLOAT((MC(I1,3)+IPI(2,J1))/2)                                 QUAN-289
      F1=(DFLOAT(MC(I1,5))-G1*(G1+1.D0))*(G2*(G2+1.D0)-DFLOAT(MC(I1,5)))QUAN-290
      G1=DFLOAT((MC(I2,3)-IPI(2,J1))/2)                                 QUAN-291
      G2=DFLOAT((MC(I2,3)+IPI(2,J1))/2)                                 QUAN-292
      F2=(DFLOAT(MC(I2,5))-G1*(G1+1.D0))*(G2*(G2+1.D0)-DFLOAT(MC(I2,5)))QUAN-293
      F3=F1*F2                                                          QUAN-294
      IF (F3.GT.0.D0) SMATEL=SMATEL-DSQRT(F3)*DJ6J(IL,ISJ,IS-2,ISJP,ILP,QUAN-295
     1IIQ,FAC,NFA)/ASPI                                                 QUAN-296
   21 A=CMATEL*AZ(5)*SMATEL                                             QUAN-297
      IF (DABS(A).LT.1.D-10) GO TO 22                                   QUAN-298
      IT=IT+1                                                           QUAN-299
      NAT(1,IT)=IQ(3,K)+INVZ+INSL                                       QUAN-300
      NAT(2,IT)=0                                                       QUAN-301
      IF (LO(14)) NAT(2,IT)=NAT(1,IT)                                   QUAN-302
      AT(2,IT)=2.D0*A*A2                                                QUAN-303
      IF (LO(133)) AT(3,IT)=0.D0                                        QUAN-304
      IF (LLP) AT(3,IT)=2.D0*A*A3                                       QUAN-305
   22 CONTINUE                                                          QUAN-306
      LLO=.FALSE.                                                       QUAN-307
   23 IF (IT-IT1) 32 , 27 , 24                                          QUAN-308
C  SUMMATION OF COEFFICIENTS RELATED TO THE SAME FORM FACTOR.           QUAN-309
   24 IT2=IT-1                                                          QUAN-310
      DO 26 I=IT1,IT2                                                   QUAN-311
      DO 25 J=I,IT2                                                     QUAN-312
      IF ((NAT(1,I).NE.NAT(1,J+1)).OR.(NAT(2,I).NE.NAT(2,J+1))) GO TO 25QUAN-313
      AT(2,I)=AT(2,I)+AT(2,J+1)                                         QUAN-314
      AT(2,J+1)=0.D0                                                    QUAN-315
      IF (JD.EQ.2) GO TO 25                                             QUAN-316
      AT(3,I)=AT(3,I)+AT(3,J+1)                                         QUAN-317
      AT(3,J+1)=0.D0                                                    QUAN-318
   25 CONTINUE                                                          QUAN-319
   26 CONTINUE                                                          QUAN-320
C ELIMINATION OF TOO SMALL COEFFICIENTS.                                QUAN-321
   27 IT2=IT                                                            QUAN-322
      IT=IT1-1                                                          QUAN-323
      NMR=0                                                             QUAN-324
      NMI=0                                                             QUAN-325
      DO 28 I=IT1,IT2                                                   QUAN-326
      CX=DABS(AT(2,I))                                                  QUAN-327
      IF (JD.EQ.3) CX=CX+DABS(AT(3,I))                                  QUAN-328
      IF (CX.LT.1.D-10) GO TO 28                                        QUAN-329
      IT=IT+1                                                           QUAN-330
      NAT(1,IT)=NAT(1,I)                                                QUAN-331
      NAT(2,IT)=NAT(2,I)                                                QUAN-332
      AT(2,IT)=AT(2,I)                                                  QUAN-333
      IF (JD.EQ.3) AT(3,IT)=AT(3,I)                                     QUAN-334
      NMR=NMR+1                                                         QUAN-335
      IF (NAT(2,IT).NE.0) NMI=NMI+1                                     QUAN-336
   28 CONTINUE                                                          QUAN-337
      IF (NMR.NE.0) NIN=NIN+1                                           QUAN-338
      IF (NMI.NE.0) NIN=NIN+1                                           QUAN-339
      IF (LLO) GO TO 36                                                 QUAN-340
      IF (.NOT.(LO(13).OR.LO(19))) GO TO 32                             QUAN-341
      NVI(I2,I1,2)=IT                                                   QUAN-342
C COEFFICIENTS OF THE DERIVATIVE COUPLING                               QUAN-343
C FORM FACTOR  V(R)/R**2  COEFFICIENT: (I2||IQ||I1)*(G1-G2)             QUAN-344
C PARAMETRISATION: (I2||IQ||I1)*(G1-G2)*AZ(6).                          QUAN-345
      IT2=IT                                                            QUAN-346
      IF (LO(100)) GO TO 30                                             QUAN-347
      DO 29 I=IT1,IT2                                                   QUAN-348
      IF (NAT(1,I).GT.INVZ) GO TO 29                                    QUAN-349
      IJ=NAT(1,I)-ITY(5)                                                QUAN-350
      IF (IVZ(3,IJ).EQ.0) GO TO 29                                      QUAN-351
      A=DFLOAT(MC(I2,6)-MC(I1,6))*AZ(6)*2.D0                            QUAN-352
      IF (A.EQ.0.D0) GO TO 29                                           QUAN-353
      LO(110)=.TRUE.                                                    QUAN-354
      IT=IT+1                                                           QUAN-355
      NAT(1,IT)=IVZ(3,IJ)+INVZ+INSL                                     QUAN-356
      NAT(2,IT)=0                                                       QUAN-357
      IF (LO(14)) NAT(2,IT)=NAT(1,IT)                                   QUAN-358
      AT(2,IT)=AT(2,I)*A                                                QUAN-359
      IF (JD.EQ.3) AT(3,IT)=0.D0                                        QUAN-360
      IF (LLP) AT(3,IT)=AT(3,I)*A                                       QUAN-361
   29 CONTINUE                                                          QUAN-362
      GO TO 33                                                          QUAN-363
   30 DO 31 I=IT1,IT2                                                   QUAN-364
      IJ=NAT(1,I)                                                       QUAN-365
      IF (IVZ(3,IJ).EQ.0) GO TO 31                                      QUAN-366
      LO(110)=.TRUE.                                                    QUAN-367
      IT=IT+1                                                           QUAN-368
      NAT(1,IT)=IVZ(3,IJ)+INTV                                          QUAN-369
      NAT(2,IT)=0                                                       QUAN-370
      AT(2,IT)=AT(2,I)                                                  QUAN-371
      AT(3,IT)=AT(2,IT)*DFLOAT(MC(I2,6)-MC(I1,6))                       QUAN-372
   31 CONTINUE                                                          QUAN-373
      IF (IT.GT.IT2) NIN=NIN+1                                          QUAN-374
      GO TO 33                                                          QUAN-375
   32 NVI(I2,I1,2)=IT                                                   QUAN-376
   33 NVI(I2,I1,3)=IT                                                   QUAN-377
      IF (LO(100).OR.(NVI(I2,I1,3).EQ.NVI(I2,I1,2))) GO TO 38           QUAN-378
C COPY OF THE COEFFICIENTS AND CORRECTIONS IN ORDER TO OBTAIN           QUAN-379
C AN HERMITIAN INTERACTION.                                             QUAN-380
      K1=NVI(I2,I1,1)                                                   QUAN-381
      K2=NVI(I2,I1,2)                                                   QUAN-382
      IT1=IT+1                                                          QUAN-383
      NVI(I1,I2,1)=IT1                                                  QUAN-384
      DO 34 K=K1,K2                                                     QUAN-385
      IT=IT+1                                                           QUAN-386
      NAT(1,IT)=NAT(1,K)                                                QUAN-387
      NAT(2,IT)=NAT(2,K)                                                QUAN-388
      IF (JD.EQ.3) AT(3,IT)=AT(3,K)                                     QUAN-389
   34 AT(2,IT)=AT(2,K)                                                  QUAN-390
      K1=K2+1                                                           QUAN-391
      K2=NVI(I2,I1,3)                                                   QUAN-392
      DO 35 K=K1,K2                                                     QUAN-393
      IT=IT+1                                                           QUAN-394
      NAT(1,IT)=NAT(1,K)-INSL                                           QUAN-395
      NAT(2,IT)=0                                                       QUAN-396
      IF (NAT(2,K).NE.0) NAT(2,IT)=NAT(2,K)-INSL                        QUAN-397
      AT(2,IT)=-AT(2,K)                                                 QUAN-398
      IF (JD.EQ.3) AT(3,IT)=-AT(3,K)                                    QUAN-399
      IT=IT+1                                                           QUAN-400
      NAT(1,IT)=NAT(1,K)                                                QUAN-401
      NAT(2,IT)=NAT(2,K)                                                QUAN-402
      IF (JD.EQ.3) AT(3,IT)=AT(3,K)                                     QUAN-403
   35 AT(2,IT)=AT(2,K)                                                  QUAN-404
      LLO=.TRUE.                                                        QUAN-405
      GO TO 23                                                          QUAN-406
   36 NVI(I1,I2,2)=IT                                                   QUAN-407
      DO 37 K=K1,K2                                                     QUAN-408
      IT=IT+1                                                           QUAN-409
      NAT(1,IT)=NAT(1,K)                                                QUAN-410
      NAT(2,IT)=NAT(2,K)                                                QUAN-411
      IF (JD.EQ.3) AT(3,IT)=-AT(3,K)                                    QUAN-412
   37 AT(2,IT)=-AT(2,K)                                                 QUAN-413
      NVI(I1,I2,3)=IT                                                   QUAN-414
      NIN=NIN+2                                                         QUAN-415
      IF (LO(14)) NIN=NIN+2                                             QUAN-416
      GO TO 40                                                          QUAN-417
C SYMMETRISATION OF THE TABLE WHEN THERE IS NO DEFORMED SPIN-ORBIT.     QUAN-418
   38 DO 39 K=1,3                                                       QUAN-419
   39 NVI(I1,I2,K)=NVI(I2,I1,K)                                         QUAN-420
   40 CONTINUE                                                          QUAN-421
   41 CONTINUE                                                          QUAN-422
      IF (.NOT.LO(58)) GO TO 42                                         QUAN-423
C OUTPUT OF COUPLING COEFFICIENTS.                                      QUAN-424
      WRITE (MW,1002) ((J,I,(NVI(J,I,K),K=1,3),I=1,NC),J=1,NC)          QUAN-425
      IF (IT.EQ.0) GO TO 42                                             QUAN-426
      IF (JD.EQ.2) WRITE (MW,1003) (I,NAT(1,I),NAT(2,I),AT(2,I),I=1,IT) QUAN-427
      IF (JD.EQ.3) WRITE (MW,1004) (I,NAT(1,I),NAT(2,I),AT(2,I),AT(3,I),QUAN-428
     1I=1,IT)                                                           QUAN-429
   42 NCI=NCIN                                                          QUAN-430
      NIC=NC                                                            QUAN-431
      IF (.NOT.LO(123)) RETURN                                          QUAN-432
C FIND QUANTUM NUMBER OF COUPLED CHANNELS FOR IDENTICAL PARTICLES. NOT  QUAN-433
C USED FOR DIRAC FORMALISM (INCORRECT MEANING OF MC(*,10[3]) IN MTCH).  QUAN-434
      NCIN=0                                                            QUAN-435
      JA=NAJ/2                                                          QUAN-436
      ISI=IPI(2,1)-1                                                    QUAN-437
      NSM=ISI+1                                                         QUAN-438
      NSB=JPI+1                                                         QUAN-439
      DO 44 IS=NSB,NSM,2                                                QUAN-440
      LP=JA+IS                                                          QUAN-441
      LM=IABS(JA-IS+1)+1                                                QUAN-442
      DO 43 L=LM,LP                                                     QUAN-443
      IF (MOD(L+JPI,2).NE.1) GO TO 43                                   QUAN-444
      NCIN=NCIN+1                                                       QUAN-445
      MC(NCIN,10)=IPI(11,1)                                             QUAN-446
      MC(NCIN,7)=1                                                      QUAN-447
      MC(NCIN,8)=L-1                                                    QUAN-448
      MC(NCIN,9)=2*IS-2                                                 QUAN-449
      MC(NCIN,11)=MC(NCIN,8)*(MC(NCIN,8)+1)                             QUAN-450
      MC(NCIN,12)=(JA*(JA+1)-L*(L-1)-IS*(IS-1))/2                       QUAN-451
   43 CONTINUE                                                          QUAN-452
   44 CONTINUE                                                          QUAN-453
      NC=NCIN                                                           QUAN-454
      IF (LO(58)) WRITE (MW,1005) NCI,NCIN                              QUAN-455
      IF (NCIN.EQ.0) RETURN                                             QUAN-456
      IF (NIC.EQ.NCI) GO TO 46                                          QUAN-457
      N=NCI+1                                                           QUAN-458
      DO 45 I=N,NIC                                                     QUAN-459
      NC=NC+1                                                           QUAN-460
      MC(NC,7)=MC(I,1)                                                  QUAN-461
      MC(NC,8)=MC(I,2)                                                  QUAN-462
      MC(NC,9)=MC(I,3)                                                  QUAN-463
      MC(NC,10)=MC(I,4)                                                 QUAN-464
      MC(NC,11)=MC(I,5)                                                 QUAN-465
   45 MC(NC,12)=MC(I,6)                                                 QUAN-466
   46 IF (.NOT.LO(58)) GO TO 48                                         QUAN-467
      WRITE (MW,1000) AJ,IP(JPI+1),NC,NCIN                              QUAN-468
      DO 47 J=1,NC                                                      QUAN-469
      SJ=0.5D0*DFLOAT(MC(J,9))                                          QUAN-470
   47 WRITE (MW,1001) J,MC(J,7),MC(J,8),SJ,MC(J,10),MC(J,11),MC(J,12)   QUAN-471
C COMPUTATION OF TRANSFORMATION COEFFICIENTS.                           QUAN-472
   48 ITI=IT+1                                                          QUAN-473
      DO 50 J=1,NCIN                                                    QUAN-474
      DO 49 I=1,NCI                                                     QUAN-475
      AA(I,J)=0.D0                                                      QUAN-476
      IF (MC(J,8).NE.MC(I,2)) GO TO 49                                  QUAN-477
      IJ=MC(I,3)                                                        QUAN-478
      IS=MC(J,9)                                                        QUAN-479
      AA(I,J)=DFLOAT(1-MOD(2*MC(I,2)+NAJ+2*ISI,4))*DSQRT((IJ+1.D0)*(IS+1QUAN-480
     1.D0))*DJ6J(2*MC(I,2),ISI,IJ,ISI,NAJ,IS,FAC,NFA)                   QUAN-481
   49 CONTINUE                                                          QUAN-482
   50 CONTINUE                                                          QUAN-483
      LO(110)=.FALSE.                                                   QUAN-484
      NIN=0                                                             QUAN-485
      NI=NCI-NCIN                                                       QUAN-486
      IBB=2                                                             QUAN-487
      IF (.NOT.((LO(101).OR.LO(103)).AND.LO(100))) IBB=1                QUAN-488
      DO 70 I1=1,NC                                                     QUAN-489
      IF (I1.GT.NCIN) GO TO 51                                          QUAN-490
      J1=1                                                              QUAN-491
      K1=NCI                                                            QUAN-492
      GO TO 52                                                          QUAN-493
   51 J1=I1+NI                                                          QUAN-494
      K1=J1                                                             QUAN-495
   52 DO 69 I2=1,NC                                                     QUAN-496
      IF ((I2.GT.I1).AND.(IBB.EQ.1)) GO TO 70                           QUAN-497
      DO 67 IB=1,IBB                                                    QUAN-498
      IF (I2.GT.NCIN) GO TO 53                                          QUAN-499
      J2=1                                                              QUAN-500
      K2=NCI                                                            QUAN-501
      GO TO 54                                                          QUAN-502
   53 J2=I2+NI                                                          QUAN-503
      K2=J2                                                             QUAN-504
   54 IT1=IT+1                                                          QUAN-505
      IF (IB.EQ.1) NVI(I2,I1,4)=IT1                                     QUAN-506
      DO 58 L1=J1,K1                                                    QUAN-507
      A1=1.D0                                                           QUAN-508
      IF (L1.LE.NCI) A1=AA(L1,I1)                                       QUAN-509
      IF (A1.EQ.0.D0) GO TO 58                                          QUAN-510
      DO 57 L2=J2,K2                                                    QUAN-511
      A2=1.D0                                                           QUAN-512
      IF (L2.LE.NCI) A2=AA(L2,I2)                                       QUAN-513
      IF (A2.EQ.0.D0) GO TO 57                                          QUAN-514
      M1=NVI(L2,L1,IB)                                                  QUAN-515
      M2=NVI(L2,L1,IB+1)                                                QUAN-516
      IF (IB.EQ.2) M1=M1+1                                              QUAN-517
      IF (M1.GT.M2) GO TO 57                                            QUAN-518
      IF (2*(M2-M1+3+IT).GT.NMAX) CALL MEMO('QUAN',NMAX,2*(M2-M1+3+IT)) QUAN-519
      DO 56 M=M1,M2                                                     QUAN-520
      IT=IT+1                                                           QUAN-521
      NAT(1,IT)=NAT(1,M)                                                QUAN-522
      NAT(2,IT)=NAT(2,M)                                                QUAN-523
      DO 55 IX=2,JD                                                     QUAN-524
   55 AT(IX,IT)=AT(IX,M)*A1*A2                                          QUAN-525
   56 CONTINUE                                                          QUAN-526
   57 CONTINUE                                                          QUAN-527
   58 CONTINUE                                                          QUAN-528
      IF (IT-IT1) 66 , 63 , 59                                          QUAN-529
C  SUMMATION OF COEFFICIENTS RELATED TO THE SAME FORM FACTOR.           QUAN-530
   59 IT2=IT-1                                                          QUAN-531
      DO 62 I=IT1,IT2                                                   QUAN-532
      DO 61 J=I,IT2                                                     QUAN-533
      IF ((NAT(1,I).NE.NAT(1,J+1)).OR.(NAT(2,I).NE.NAT(2,J+1))) GO TO 61QUAN-534
      DO 60 IX=2,JD                                                     QUAN-535
      AT(IX,I)=AT(IX,I)+AT(IX,J+1)                                      QUAN-536
   60 AT(IX,J+1)=0.D0                                                   QUAN-537
   61 CONTINUE                                                          QUAN-538
   62 CONTINUE                                                          QUAN-539
C ELIMINATION OF TOO SMALL COEFFICIENTS.                                QUAN-540
   63 IT2=IT                                                            QUAN-541
      IT=IT1-1                                                          QUAN-542
      NMR=0                                                             QUAN-543
      NMI=0                                                             QUAN-544
      DO 65 I=IT1,IT2                                                   QUAN-545
      CX=DABS(AT(2,I))                                                  QUAN-546
      IF (LO(100)) CX=CX+DABS(AT(3,I))                                  QUAN-547
      IF (CX.LT.1.D-10) GO TO 65                                        QUAN-548
      IT=IT+1                                                           QUAN-549
      NAT(1,IT)=NAT(1,I)                                                QUAN-550
      NAT(2,IT)=NAT(2,I)                                                QUAN-551
      DO 64 IX=2,JD                                                     QUAN-552
   64 AT(IX,IT)=AT(IX,I)                                                QUAN-553
      IF (LO(100)) AT(3,IT)=AT(3,I)                                     QUAN-554
      NMR=NMR+1                                                         QUAN-555
      IF (NAT(2,IT).NE.0) NMI=NMI+1                                     QUAN-556
   65 CONTINUE                                                          QUAN-557
      IF (NMR.NE.0) NIN=NIN+1                                           QUAN-558
      IF (NMI.NE.0) NIN=NIN+1                                           QUAN-559
      IF ((IB.EQ.2).AND.(NMR.NE.0)) LO(110)=.TRUE.                      QUAN-560
   66 NVI(I2,I1,IB+4)=IT                                                QUAN-561
   67 CONTINUE                                                          QUAN-562
      IF (IBB.EQ.2) GO TO 70                                            QUAN-563
      NVI(I2,I1,6)=IT                                                   QUAN-564
C SYMMETRISATION OF THE TABLE WHEN THERE IS NO DEFORMED SPIN-ORBIT.     QUAN-565
      DO 68 K=4,6                                                       QUAN-566
   68 NVI(I1,I2,K)=NVI(I2,I1,K)                                         QUAN-567
   69 CONTINUE                                                          QUAN-568
   70 CONTINUE                                                          QUAN-569
      IF (.NOT.LO(58)) RETURN                                           QUAN-570
C OUTPUT OF COUPLING COEFFICIENTS.                                      QUAN-571
      WRITE (MW,1002) ((I,J,(NVI(I,J,K),K=4,6),I=1,NC),J=1,NC)          QUAN-572
      IF (IT.LT.ITI) RETURN                                             QUAN-573
      IF (.NOT.LO(100)) WRITE (MW,1003) (I,NAT(1,I),NAT(2,I),AT(2,I),I=IQUAN-574
     1TI,IT)                                                            QUAN-575
      IF (LO(100)) WRITE (MW,1004) (I,NAT(1,I),NAT(2,I),AT(2,I),AT(3,I),QUAN-576
     1I=ITI,IT)                                                         QUAN-577
      RETURN                                                            QUAN-578
 1000 FORMAT (/' CHANNEL SPIN AND PARITY',F6.1,A1,I11,' COUPLED CHANNELSQUAN-579
     1 AND',I3,' SOLUTIONS'//8X,' I',3X,' V',3X,' L',3X,' J',9X,' POT',5QUAN-580
     2X,' CL',6X,' CJ'/)                                                QUAN-581
 1001 FORMAT (5X,3I5,F6.1,I11,2I9)                                      QUAN-582
 1002 FORMAT (/' PAIRS OF CHANNELS N1 N2, AND COUPLING COEFFICIENT NUMBEQUAN-583
     1RS NVI(N1,N2,K),K=1,3'/(1X,6(I3,I3,',',3I4,';')))                 QUAN-584
 1003 FORMAT (//5X,'COEFFICIENTS'/(4(2X,3I3,1P,D15.6)))                 QUAN-585
 1004 FORMAT (//5X,'COEFFICIENTS'/(3(2X,3I3,1P,2D15.6)))                QUAN-586
 1005 FORMAT (/' NUMBER OF SOLUTIONS REDUCED FROM',I3,' TO',I3)         QUAN-587
      END                                                               QUAN-588
