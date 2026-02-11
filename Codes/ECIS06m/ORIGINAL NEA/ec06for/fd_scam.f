C 03/03/07                                                      ECIS06  SCAM-000
      SUBROUTINE SCAM(F,FCN,TX,GCN,JMAX,KMAX,MC,MD,FAR,FAI,WV,NCOLL,NCOLSCAM-001
     1S,KAB,KBA,IPI,IPIM,GAM,FISS,TL,NCJ,XD,P,Q,V,NCT,IQ,AA,KBC,IDT,LO) SCAM-002
C STORAGE OF SCATTERING COEFFICIENTS.                                   SCAM-003
C INPUT:     JMAX:    MAXIMUM NUMBER OF CHANNEL SPINS, DIMENSION FOR F. SCAM-004
C            KMAX:    MAXIMUM NUMBER OF L OF COMPOUND NUCLEUS,          SCAM-005
C                     DIMENSION FOR FCN, GCN.                           SCAM-006
C            MC:      NUCLEAR STATE NUMBERS, ANGULAR MOMENTA,  SEE QUAN.SCAM-007
C            MD:      EXTENDED MC FOR IDENTICAL PARTICLES WITH SPIN.    SCAM-008
C            FAR/FAI: REAL/IMAGINARY PART OF SCATTERING COEFFICIENTS.   SCAM-009
C            WV:      WAVE NUMBER AND COULOMB PARAMETER.  SEE COLF.     SCAM-010
C            NCOLL:   NUMBER OF COUPLED LEVELS.                         SCAM-011
C            NCOLS:   NUMBER OF LEVELS WITH ANGULAR DISTRIBUTION.       SCAM-012
C            KAB:     DIMENSION FOR MC.                                 SCAM-013
C            KBA:     NUMBER OF INDEPENDENT AMPLITUDES.                 SCAM-014
C            IPI:     PARITY, MULTIPLICITIES, ADDRESSES IN F (SEE CALX).SCAM-015
C            IPIM:    IPI FOR CONTINUA OF COMPOUND NUCLEUS.             SCAM-016
C            GAM:     GAMMA TRANSMISSION COEFFICIENTS READ.             SCAM-017
C            FISS:    FISSION COEFFICIENTS FOR COMPOUND NUCLEUS.        SCAM-018
C            TL:      TRANSMISSION COEFFICIENTS OF UNCOUPLED LEVELS.    SCAM-019
C            NCJ:     STARTING AND FINAL ADDRESSES FOR CONTINUA.        SCAM-020
C            XD:      ENERGY AND SPIN DEPENDENCE OF LEVEL DENSITIES.    SCAM-021
C            NCT:     NUMBER OF EQUATIONS AND SOLUTIONS FOR EACH PARITY.SCAM-022
                                                                        SCAM-023
C            IQ:      DIMENSION FOR V.                                  SCAM-024
C            AA:      COEFFICIENTS OF SYMMETRISATION.                   SCAM-025
C            KBC:     DIMENSION FOR AA.                                 SCAM-026
C            IDT:     SIZE AVAILABLE FOR THE Q IN DOUBLE PRECISION.     SCAM-027
C            LO(I):   LOGICAL CONTROLS:                                 SCAM-028
C               LO(18) =.TRUE. PROJECTILE-TARGET ANTISYMMETRISATION.    SCAM-029
C               LO(55) =.TRUE. OUTPUT OF C-MATRIX ELEMENTS AND OF       SCAM-030
C                              COMPOUND NUCLEUS INTERMEDIATE RESULTS.   SCAM-031
C               LO(63) =.TRUE. PENETRABILITIES WRITTEN ON FILE 63.      SCAM-032
C               LO(82) =.TRUE. OLD SIMPLIFIED COMPOUND NUCLEUS.         SCAM-033
C               LO(83) =.TRUE. NO ENGELBRETCH-WEIDENMULLER TRANSFORM.   SCAM-034
C               LO(85) =.TRUE. FISSION TRANSMISSION COEFFICIENTS.       SCAM-035
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      SCAM-036
C               LO(87) =.TRUE. NO WIDTH FLUCTUATIONS.                   SCAM-037
C               LO(123)=.TRUE. IDENTICAL PARTICLES WITH SPIN.           SCAM-038
C               LO(124)=.TRUE. COMPUTE TRANSMISSION COEFFICIENTS.       SCAM-039
C               LO(125)=.TRUE. USUAL COUPLED EQUATIONS.                 SCAM-040
C               LO(131)=.TRUE. TOTAL SPIN IS TOO LARGE FOR COMPOUND     SCAM-041
C                              NUCLEUS.                                 SCAM-042
C               LO(132)=.TRUE. STORE FISSION AND GAMMA TRANSMISSION     SCAM-043
C                              COEFFICIENTS FOR INTERPOLATION.          SCAM-044
C OUTPUT:    F:       SCATTERING COEFFICIENTS.                          SCAM-045
C            FCN:     COMPOUND NUCLEUS CONTRIBUTIONS.                   SCAM-046
C            TX:      TOTAL REACTION CROSS SECTION IN MB FOLLOWED BY    SCAM-047
C                     THE TOTAL DIRECT CROSS SECTIONS OF COUPLED        SCAM-048
C                     LEVELS, THE TOTAL COMPOUND CROSS SECTIONS OF      SCAM-049
C                     LEVELS, THE FISSION AND GAMMA CROSS SECTIONS.     SCAM-050
C            GCN:     COMPOUND NUCLEUS COEFFICIENTS OF LEGENDRE         SCAM-051
C                     POLYNOMIALS.                                      SCAM-052
C WORKING AREAS:                                                        SCAM-053
C            P:       FOR DIAGONALISATION OF THE S MATRIX.              SCAM-054
C            Q:       FOR ANGULAR DISTRIBUTION OF COMPOUND NUCLEUS.     SCAM-055
C            V(1,*):  LEVEL AND QUANTUM NUMBERS.                        SCAM-056
C            V(2,*):  WEIGH FOR CONTINUA.                               SCAM-057
C            V(3,*):  TRANSMISSION COEFFICIENTS.                        SCAM-058
C            V(4,*):  WIDTH FLUCTUATION PARAMETER.                      SCAM-059
C            V(5,*):  CONTRIBUTION TO COMPOUND NUCLEUS.                 SCAM-060
C                                                                       SCAM-061
C FOR THE COMMON  /DCONS/ SEE CALC.                                     SCAM-062
C FOR THE COMMON  /NCOMP/ SEE CALX.                                     SCAM-063
C FOR THE COMMON  /NOEQU/ SEE QUAN.                                     SCAM-064
C                                                                       SCAM-065
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     SCAM-066
C  XZ:        CONVERSION FACTOR TO MILLIBARNS.                          SCAM-067
C   USED:     XZ.                                                       SCAM-068
C                                                                       SCAM-069
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     SCAM-070
C  NFISS:     NUMBER OF FISSION TRANSMISSION COEFFICIENTS.              SCAM-071
C  NRD:       NUMBER OF GAMMA TRANSMISSION COEFFICIENTS.                SCAM-072
C  NCONT:     NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS.                 SCAM-073
C  NCOJ:      NUMBER OF SPINS OF THE TARGET FOR A CONTINUUM.            SCAM-074
C  NCOLX:     TOTAL NUMBER OF LEVELS WITHOUT DISCRETISATION.            SCAM-075
C  NDP:       ADDRESS OF FISSION AND GAMMA FINAL RESULTS.               SCAM-076
C  NDQ:       ADDRESS OF FISSION AND GAMMA INTERMEDIATE RESULTS.        SCAM-077
C  BZ:        HAUSER-FESHBACH AND MOLDAUER'S PARAMETERS DESCRIBED BELOW.SCAM-078
C   BZ(1):    SQUARE ROOT OF ELASTIC ENHANCEMENT.                       SCAM-079
C   BZ(2):    IF LO(82)=.TRUE., SPIN CUT-OFF PARAMETER.                 SCAM-080
C             IF LO(82)=.FALSE., PARTICLE DEGREES OF FREEDOM.           SCAM-081
C   BZ(3):    SQUARE ROOT OF LEVEL DENSITY PARAMETER. IF LO(82)=LO(87)= SCAM-082
C             .FALSE., PARAMETER BZ(3) IN MOLDAUER'S FORMULA OF INPUT   SCAM-083
C             DESCRIPTION.                                              SCAM-084
C   BZ(4):    IF LO(82)=LO(87)=.FALSE., PARAMETER BZ(4) IN SAME FORMULA.SCAM-085
C   BZ(5):    IF LO(82)=LO(87)=.FALSE., PARAMETER BZ(5) IN SAME FORMULA.SCAM-086
C  FNUG:      RADIATIVE DEGREES OF FREEDOM.                             SCAM-087
C  TG1:       DERIVED DATA FOR GAMMA IN COMPOUND NUCLEUS.               SCAM-088
C  SGSQ:      DERIVED DATA FOR GAMMA IN COMPOUND NUCLEUS.               SCAM-089
C   USED:     NFISS,NRD,NCONT,NCOJ,NCOLX,NDP,NDQ,BZ,FNUG,TG1,SGSQ       SCAM-090
C                                                                       SCAM-091
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NOEQU/:                     SCAM-092
C  NIC:       NUMBER OF EQUATIONS AT THE CHANNEL SPIN.                  SCAM-093
C  NCI:       NUMBER OF SOLUTIONS AT THE CHANNEL SPIN.                  SCAM-094
C  NC:        NUMBER OF EQUATIONS FOR IDENTICAL PARTICLES.              SCAM-095
C  NCIN:      NUMBER OF SOLUTIONS FOR IDENTICAL PARTICLES.              SCAM-096
C  JPI:       PARITY 0 OR 1.                                            SCAM-097
C  IPJ:       VALUE OF J+1 OR J+1/2 WHERE J IS THE CHANNEL SPIN.        SCAM-098
C  R1(2):     MAXIMUM OF SCATTERING AND COMPOUND COEFFICIENT.           SCAM-099
C  NAJ:       TWICE THE CHANNEL SPIN J.                                 SCAM-100
C   DEFINED:  R1.                                                       SCAM-101
C   USED:     NIC,NCI,NC,NCIN,JPI,IPJ,R1,NAJ.                           SCAM-102
C                                                                       SCAM-103
C***********************************************************************SCAM-104
      IMPLICIT REAL*8 (A-H,O-Z)                                         SCAM-105
      LOGICAL LO(150),LG(2)                                             SCAM-106
      DIMENSION F(2,JMAX,*),FCN(KMAX,*),TX(*),GCN(KMAX,*),MC(KAB,12),MD(SCAM-107
     1KAB,12),FAR(KAB,*),FAI(KAB,*),WV(22,*),IPI(11,*),IPIM(11,*),GAM(*)SCAM-108
     2,FISS(2,*),TL(*),NCJ(2,*),XD(3,*),P(NC,NC,4),Q(*),V(IQ,*),NCT(6),ASCAM-109
     3A(KBC,*),X(20),W(20),TP(20),SGF(2),TG(2),FNU(2),MCX(4,2),VCX(4,2) SCAM-110
      CHARACTER*1 IP(2)                                                 SCAM-111
      CHARACTER*8 AL(2)                                                 SCAM-112
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            SCAM-113
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQSCAM-114
     1,ACN(8),BZ(5),TG0,BN,FNUG,EGD,GGD,TG1,SGSQ                        SCAM-115
      COMMON /NOEQU/ NCXN,NIC,NCI,NC,NCIN,NIN,JPI,IPJ,R1(2),NAJ         SCAM-116
      COMMON /INOUT/ MR,MW,MS                                           SCAM-117
      DATA IP,AL,PI,NSY /'+','-',' FISSION','   GAMMA',3.141592653589793SCAM-118
     12D0,0/                                                            SCAM-119
      DATA X /7.0539889691988753D-02,3.7212681800161144D-01,9.1658210248SCAM-120
     1327356D-01,1.7073065310283439D+00,2.7491992553094321D+00,4.0489253SCAM-121
     2138508869D+00,5.6151749708616165D+00,7.4590174536710633D+00,9.5943SCAM-122
     3928695810968D+00,1.2038802546964316D+01,1.4814293442630740D+01,1.7SCAM-123
     4948895520519376D+01,2.1478788240285011D+01,2.5451702793186906D+01,SCAM-124
     52.9932554631700612D+01,3.5013434240479000D+01,4.0833057056728571D+SCAM-125
     601,4.7619994047346502D+01,5.5810795750063899D+01,6.652441652561575SCAM-126
     74D+01/                                                            SCAM-127
      DATA W /1.8108006241898926D-01,4.2255676787856397D-01,6.6690954670SCAM-128
     1184815D-01,9.1535237278307367D-01,1.1695397071955460D+00,1.4313549SCAM-129
     2859282060D+00,1.7029811379850227D+00,1.9870158907927472D+00,2.2866SCAM-130
     3357812534308D+00,2.6058347275538333D+00,2.9497837342139509D+00,3.3SCAM-131
     4253957820093196D+00,3.7422554705898109D+00,4.2142367102518804D+00,SCAM-132
     54.7625184614902093D+00,5.4217260442455743D+00,6.2540123569324213D+SCAM-133
     600,7.3873143890544346D+00,9.1513287309874796D+00,1.289338864593999SCAM-134
     77D+01/                                                            SCAM-135
      IK(I,L,J)=((J-NAJ+IPI(3,I)-1)/2*IPI(2,I)+(L-J+IPI(2,I)-1)/2)/2    SCAM-136
      NZ=MOD(IPJ+JPI+1,2)                                               SCAM-137
      NS=0                                                              SCAM-138
      IF (NZ.NE.0) NS=NCT(1)*NCT(3)                                     SCAM-139
      R1(1)=0.D0                                                        SCAM-140
      R1(2)=0.D0                                                        SCAM-141
      JC=0                                                              SCAM-142
      AJ=0.5D0*DFLOAT(NAJ)                                              SCAM-143
      DZ=1.D0                                                           SCAM-144
      IF (LO(18)) DZ=2.D0                                               SCAM-145
      RZ=PI*DZ*XZ*DFLOAT(NAJ+1)                                         SCAM-146
      IF (LO(55)) WRITE (MW,1000) AJ,IP(JPI+1),NC,NCIN                  SCAM-147
C COMPUTATION OF TRANSMISSION COEFFICIENTS.                             SCAM-148
      NXX=NCIN                                                          SCAM-149
      IF (LO(124)) NXX=NC                                               SCAM-150
      DO 5 IC=1,NXX                                                     SCAM-151
      IV=MC(IC,1)                                                       SCAM-152
      IF (WV(3,IV).LT.0.D0) GO TO 5                                     SCAM-153
      IF (IC.LE.NCIN) TX(1)=TX(1)+4.D0*FAI(IC,IC)*RZ*DZ                 SCAM-154
      IF (.NOT.LO(124)) GO TO 2                                         SCAM-155
      JC=JC+1                                                           SCAM-156
      V(1,JC)=131072.D0*(65536.D0*DFLOAT(IV)+DFLOAT(MC(IC,2)))+DFLOAT(MCSCAM-157
     1(IC,3))                                                           SCAM-158
      V(2,JC)=1.D0                                                      SCAM-159
      V(3,JC)=4.D0*FAI(IC,IC)                                           SCAM-160
      DO 1 N=4,IQ                                                       SCAM-161
    1 V(N,JC)=0.D0                                                      SCAM-162
    2 DO 4 ICP=1,NC                                                     SCAM-163
      IVQ=MC(ICP,1)                                                     SCAM-164
      A1=DSQRT(WV(11,IVQ)/WV(11,IV))                                    SCAM-165
      IF (LO(125)) A1=1.D0/A1                                           SCAM-166
      FAR(ICP,IC)=FAR(ICP,IC)*A1                                        SCAM-167
      FAI(ICP,IC)=FAI(ICP,IC)*A1                                        SCAM-168
      IF (WV(3,IVQ).LT.0.D0) GO TO 3                                    SCAM-169
      IF (LO(124)) V(3,JC)=V(3,JC)-4.D0*(FAR(ICP,IC)**2+FAI(ICP,IC)**2) SCAM-170
      IF (IC.GT.NCIN) GO TO 4                                           SCAM-171
C TEST OF CONVERGENCE.                                                  SCAM-172
      R1(1)=DMAX1(R1(1),DABS(FAR(ICP,IC))+DABS(FAI(ICP,IC)))            SCAM-173
      A1=FAR(ICP,IC)*FAR(ICP,IC)+FAI(ICP,IC)*FAI(ICP,IC)                SCAM-174
      IF (A1.GT.1.D0) WRITE (MW,1001) AJ,IP(JPI+1),NC,NCIN,IC,ICP,FAR(ICSCAM-175
     1P,IC),FAI(ICP,IC)                                                 SCAM-176
      TX(IVQ+1)=TX(IVQ+1)+4.D0*A1*RZ*DZ                                 SCAM-177
      FJ=0.5D0*DFLOAT(MC(ICP,3))                                        SCAM-178
      IF (.NOT.LO(55)) GO TO 4                                          SCAM-179
C PRINT OUT OF THE AMPLITUDES.                                          SCAM-180
      A1=DSQRT(A1)                                                      SCAM-181
      B1=0.D0                                                           SCAM-182
      IF (A1.NE.0.D0) B1=DATAN2(FAI(ICP,IC),FAR(ICP,IC))                SCAM-183
      WRITE (MW,1002) ICP,IC,IVQ,MC(ICP,2),FJ,FAR(ICP,IC),FAI(ICP,IC),A1SCAM-184
     1,B1                                                               SCAM-185
      GO TO 4                                                           SCAM-186
    3 IF ((IC.LE.NCIN).AND.LO(55)) WRITE (MW,1003) ICP,IC,IVQ,MC(ICP,2),SCAM-187
     1FJ,FAR(ICP,IC),FAI(ICP,IC)                                        SCAM-188
    4 CONTINUE                                                          SCAM-189
      IF (.NOT.LO(124)) GO TO 5                                         SCAM-190
      IF (V(3,JC).LT.0.D0) V(3,JC)=0.D0                                 SCAM-191
      IF (IC.LE.NCIN) R1(2)=DMAX1(R1(2),V(3,JC))                        SCAM-192
    5 CONTINUE                                                          SCAM-193
      IF ((.NOT.(LO(63).OR.LO(55))).OR.(JC.EQ.0)) GO TO 8               SCAM-194
C OUTPUT OF TRANSMISSION COEFFICIENTS FOR COUPLED CHANNELS.             SCAM-195
      IF (LO(55)) WRITE (MW,1004) AJ,IP(JPI+1)                          SCAM-196
      IC=1+(JC-1)/4                                                     SCAM-197
      IF (LO(63)) WRITE (99,1005) AJ,IP(JPI+1),JC                       SCAM-198
      DO 7 I=1,IC                                                       SCAM-199
      J1=4*(I-1)                                                        SCAM-200
      J2=MIN0(JC-J1,4)                                                  SCAM-201
      DO 6 J=1,J2                                                       SCAM-202
      N=IDINT(V(1,J+J1)/131072.D0+.01D0)                                SCAM-203
      VCX(J,1)=.5D0*DMOD(V(1,J+J1),131072.D0)                           SCAM-204
      VCX(J,2)=V(3,J+J1)                                                SCAM-205
      MCX(J,2)=MOD(N,65536)                                             SCAM-206
    6 MCX(J,1)=N/65536                                                  SCAM-207
      IF (LO(63)) WRITE (99,1006) (MCX(J,1),MCX(J,2),VCX(J,1),VCX(J,2),JSCAM-208
     1=1,J2)                                                            SCAM-209
      IF (LO(55)) WRITE (MW,1007) (MCX(J,1),MCX(J,2),VCX(J,1),VCX(J,2),JSCAM-210
     1=1,J2)                                                            SCAM-211
    7 CONTINUE                                                          SCAM-212
    8 IF (LO(131)) GO TO 62                                             SCAM-213
      NSS=NS                                                            SCAM-214
      IF (NZ.NE.0) NSS=NCT(5)*NCT(3)                                    SCAM-215
C FISSION COEFFICIENT.                                                  SCAM-216
      TG(1)=0.D0                                                        SCAM-217
      FNU(1)=.5D0                                                       SCAM-218
      IF (.NOT.LO(86)) GO TO 9                                          SCAM-219
      KN=2*IPJ+JPI-1                                                    SCAM-220
      IF (KN.GE.NFISS) GO TO 9                                          SCAM-221
      TG(1)=FISS(1,KN)                                                  SCAM-222
      FNU(1)=.5D0*FISS(2,KN)                                            SCAM-223
    9 LG(1)=(.NOT.LO(85)).OR.TG(1).EQ.0.D0                              SCAM-224
C GAMMA COEFFICIENT.                                                    SCAM-225
      TG(2)=0.D0                                                        SCAM-226
      FNU(2)=.5D0*FNUG                                                  SCAM-227
      IF (.NOT.LO(86)) GO TO 12                                         SCAM-228
      IF (NRD.NE.0) GO TO 11                                            SCAM-229
      A1=0.D0                                                           SCAM-230
      N1=NAJ-2                                                          SCAM-231
      N1=IABS(N1)                                                       SCAM-232
      N2=NAJ+2                                                          SCAM-233
      DO 10 J=N1,N2,2                                                   SCAM-234
      A2=-DFLOAT((J+1)*(J+1))/(4.D0*SGSQ)                               SCAM-235
   10 A1=A1+DEXP(A2)*DFLOAT(J+1)/SGSQ                                   SCAM-236
      TG(2)=6.283185307D0*TG1*A1                                        SCAM-237
      GO TO 12                                                          SCAM-238
   11 IF (IPJ.LE.NRD) TG(2)=GAM(IPJ)                                    SCAM-239
   12 LG(2)=(.NOT.LO(86)).OR.TG(2).EQ.0.D0                              SCAM-240
      NC1=JC                                                            SCAM-241
      NJC=JC                                                            SCAM-242
      IF (NCOLL.EQ.NCOLX) GO TO 24                                      SCAM-243
C TRANSMISSION COEFFICIENTS FOR UNCOUPLED STATES.                       SCAM-244
      NSP1=NCOLX-NCOLL-NCONT                                            SCAM-245
      IF (NSP1.LT.1) GO TO 17                                           SCAM-246
      DO 16 I=1,NSP1                                                    SCAM-247
      II=I+NCOLL                                                        SCAM-248
      IF (WV(3,II).LT.0.D0) GO TO 16                                    SCAM-249
      NM=IPI(3,II)                                                      SCAM-250
      MN=IPI(2,II)                                                      SCAM-251
      DO 15 J=1,NM                                                      SCAM-252
      NN=NAJ-NM-1+2*J                                                   SCAM-253
      IF (NN.LT.IABS(NM-1-NAJ)) GO TO 15                                SCAM-254
      DO 14 K=1,MN                                                      SCAM-255
      MM=NN+MN+1-2*K                                                    SCAM-256
      IF (MM.LT.IABS(MN-1-NN).OR.MM.GE.2*IPI(10,II)+2) GO TO 14         SCAM-257
      L=MM/2                                                            SCAM-258
      IF (MOD(L+JPI+IPI(1,II),2).NE.0) GO TO 14                         SCAM-259
      M=MN*L+K+IPI(8,II)                                                SCAM-260
      JC=JC+1                                                           SCAM-261
      V(1,JC)=131072.D0*(65536.D0*DFLOAT(II)+DFLOAT(L))+DFLOAT(NN)      SCAM-262
      V(2,JC)=1.D0                                                      SCAM-263
      V(3,JC)=TL(M)                                                     SCAM-264
      DO 13 N=4,IQ                                                      SCAM-265
   13 V(N,JC)=0.D0                                                      SCAM-266
   14 CONTINUE                                                          SCAM-267
   15 CONTINUE                                                          SCAM-268
      IF (II.LE.NCOLS) NC1=JC                                           SCAM-269
   16 CONTINUE                                                          SCAM-270
   17 IF (NCONT.EQ.0) GO TO 24                                          SCAM-271
C TRANSMISSION COEFFICIENTS FOR CONTINUA.                               SCAM-272
      DO 23 I=1,NCONT                                                   SCAM-273
      IJ=NCJ(1,I)                                                       SCAM-274
      JI=NCJ(2,I)                                                       SCAM-275
      IF (IJ.GT.JI) GO TO 23                                            SCAM-276
      MN=IPIM(2,IJ)                                                     SCAM-277
      IT=MOD(IPIM(3,IJ)+1,2)                                            SCAM-278
      DO 22 II=IJ,JI                                                    SCAM-279
      L=IPIM(10,II)+1                                                   SCAM-280
      IF (L.EQ.0) GO TO 22                                              SCAM-281
      M=IPIM(8,II)                                                      SCAM-282
      II1=II+NCOLL+NSP1                                                 SCAM-283
      DO 21 LJ=1,L                                                      SCAM-284
      DO 20 JL=1,MN                                                     SCAM-285
      M=M+1                                                             SCAM-286
      IF (TL(M).EQ.0.D0) GO TO 20                                       SCAM-287
      JJ=2*(LJ+JL)-MN-3                                                 SCAM-288
      IF (JJ.LT.0) GO TO 20                                             SCAM-289
      IKM=IABS(JJ-NAJ)/2+1                                              SCAM-290
      IKP=MIN0((JJ+NAJ)/2+1,NCOJ)                                       SCAM-291
      IF (IKM.GT.IKP) GO TO 20                                          SCAM-292
      AI=0.D0                                                           SCAM-293
      DO 18 KK=IKM,IKP                                                  SCAM-294
      AR=.5D0*DFLOAT(IT+2*KK-1)                                         SCAM-295
   18 AI=AI+AR*DEXP(-AR*AR/XD(3,II))/XD(3,II)                           SCAM-296
      JC=JC+1                                                           SCAM-297
      V(1,JC)=131072.D0*(65536.D0*DFLOAT(II1)+DFLOAT(LJ-1))+DFLOAT(JJ)  SCAM-298
      V(2,JC)=AI*XD(1,II)                                               SCAM-299
      V(3,JC)=TL(M)                                                     SCAM-300
      DO 19 N=4,IQ                                                      SCAM-301
   19 V(N,JC)=0.D0                                                      SCAM-302
   20 CONTINUE                                                          SCAM-303
   21 CONTINUE                                                          SCAM-304
   22 CONTINUE                                                          SCAM-305
   23 CONTINUE                                                          SCAM-306
C COMPOUND NUCLEUS.                                                     SCAM-307
   24 SGF(1)=0.D0                                                       SCAM-308
      SGF(2)=0.D0                                                       SCAM-309
      IF (LO(82)) GO TO 58                                              SCAM-310
      IF (LO(83)) GO TO 28                                              SCAM-311
C COMPUTATION OF SATCHLER P-MATRIX ("PR","PI").                         SCAM-312
      I=0                                                               SCAM-313
      DO 27 II=1,NC                                                     SCAM-314
      IF (WV(3,MC(II,1)).LT.0.D0) GO TO 27                              SCAM-315
      I=I+1                                                             SCAM-316
      J=0                                                               SCAM-317
      DO 26 JJ=1,NC                                                     SCAM-318
      IF (WV(3,MC(JJ,1)).LT.0.D0) GO TO 26                              SCAM-319
      J=J+1                                                             SCAM-320
      P(I,J,1)=2.D0*(FAI(I,J)+FAI(J,I))                                 SCAM-321
      P(I,J,2)=0.D0                                                     SCAM-322
      P(I,J,3)=0.D0                                                     SCAM-323
      P(I,J,4)=0.D0                                                     SCAM-324
      K=0                                                               SCAM-325
      DO 25 KK=1,NC                                                     SCAM-326
      IF (WV(3,MC(KK,1)).LT.0.D0) GO TO 25                              SCAM-327
      K=K+1                                                             SCAM-328
      P(I,J,1)=P(I,J,1)-4.D0*(FAR(I,K)*FAR(J,K)+FAI(I,K)*FAI(J,K))      SCAM-329
   25 P(I,J,2)=P(I,J,2)+4.D0*(FAR(I,K)*FAI(J,K)-FAI(I,K)*FAR(J,K))      SCAM-330
   26 CONTINUE                                                          SCAM-331
   27 P(I,I,3)=1.D0                                                     SCAM-332
      CALL DIAG(P,P(1,1,2),P(1,1,3),P(1,1,4),NC,NJC,1.D-12,A1,IERR)     SCAM-333
      IF (IERR.EQ.0) GO TO 28                                           SCAM-334
      WRITE (MW,1008)                                                   SCAM-335
      LO(83)=.TRUE.                                                     SCAM-336
   28 BIR=1.D-12+TG(1)+TG(2)                                            SCAM-337
      DO 29 IC=1,JC                                                     SCAM-338
      IF ((.NOT.LO(83)).AND.IC.LE.NJC) V(3,IC)=DMAX1(0.D0,P(IC,IC,1))   SCAM-339
   29 BIR=BIR+V(3,IC)*V(2,IC)                                           SCAM-340
C FLUCTUATION PARAMETER NU=2*FNU, (P.A.M.,N.P.A344(1980)185).           SCAM-341
C COMMON FACTOR "TP" OF WIDTH FLUCTUATION INTEGRAL                      SCAM-342
C (P.A.M.,PRC 11(1975)426).                                             SCAM-343
      IF (LO(87).OR.(BIR.LT.0.0001D0)) GO TO 39                         SCAM-344
      DO 30 M=1,20                                                      SCAM-345
   30 TP(M)=1.D0                                                        SCAM-346
      EFB=DEXP(-BZ(5)*BIR)                                              SCAM-347
      DO 36 IC=1,JC                                                     SCAM-348
      IF (BZ(2).NE.0.D0) GO TO 31                                       SCAM-349
      V(4,IC)=(1.D0+BZ(4)+(V(3,IC)**BZ(3)-BZ(4))*EFB)/2.D0              SCAM-350
      GO TO 32                                                          SCAM-351
   31 V(4,IC)=0.5D0*BZ(2)                                               SCAM-352
   32 A1=V(3,IC)/V(4,IC)/BIR                                            SCAM-353
      IF (A1.EQ.0.D0) GO TO 36                                          SCAM-354
      IF (A1.GT.1.D-9) GO TO 34                                         SCAM-355
      DO 33 M=1,20                                                      SCAM-356
   33 TP(M)=TP(M)*DEXP(X(M)*V(3,IC)*V(2,IC)/BIR*(1.D0-.5D0*X(M)*A1))    SCAM-357
      GO TO 36                                                          SCAM-358
   34 DO 35 M=1,20                                                      SCAM-359
   35 TP(M)=TP(M)*((BIR+X(M)*V(3,IC)/V(4,IC))/BIR)**(V(4,IC)*V(2,IC))   SCAM-360
   36 CONTINUE                                                          SCAM-361
      DO 38 I=1,2                                                       SCAM-362
      IF (LG(I)) GO TO 38                                               SCAM-363
      DO 37 M=1,20                                                      SCAM-364
   37 TP(M)=TP(M)*(1.D0+X(M)*TG(I)/(BIR*FNU(I)))**FNU(I)                SCAM-365
   38 CONTINUE                                                          SCAM-366
C STORAGE OF COMPOUND TERMS.                                            SCAM-367
   39 TQ=0.D0                                                           SCAM-368
      IF (.NOT.LO(55)) GO TO 40                                         SCAM-369
      WRITE (MW,1009) AJ,IP(JPI+1)                                      SCAM-370
      IF (.NOT.LO(83)) WRITE (MW,1010)                                  SCAM-371
   40 NCX=NJC                                                           SCAM-372
      IF (LO(83)) NCX=NCIN                                              SCAM-373
      DO 56 IC=1,NCX                                                    SCAM-374
      G=1.D0                                                            SCAM-375
      BRI=BIR                                                           SCAM-376
      IF (LO(87).AND.(IC.LE.NCIN)) BRI=BRI+V(3,IC)*BZ(3)                SCAM-377
      AR=1.D0                                                           SCAM-378
      IF (LO(83)) GO TO 42                                              SCAM-379
      AR=0.D0                                                           SCAM-380
      DO 41 IA=1,NCIN                                                   SCAM-381
   41 AR=AR+(P(IA,IC,3)**2+P(IA,IC,4)**2)                               SCAM-382
   42 DO 52 ICP=1,JC                                                    SCAM-383
      IV=IDINT(V(1,ICP)/131072.D0/65536.D0+.01D0)                       SCAM-384
      IF (LO(87).OR.(BIR.LT.0.0001D0)) GO TO 44                         SCAM-385
      G=0.D0                                                            SCAM-386
      DO 43 M=1,20                                                      SCAM-387
   43 G=G+W(M)/(TP(M)*(1.D0+X(M)*V(3,IC)/(V(4,IC)*BRI))*(1.D0+X(M)*V(3,ISCAM-388
     1CP)/(V(4,ICP)*BRI)))                                              SCAM-389
      IF (ICP.EQ.IC) G=G+G/V(4,IC)                                      SCAM-390
   44 TQ=V(3,IC)*V(3,ICP)*G/BRI                                         SCAM-391
      IF (LO(87).AND.(ICP.EQ.IC)) TQ=TQ+TQ*BZ(3)                        SCAM-392
      IF ((ICP.GT.NC1).OR.LO(83)) GO TO 50                              SCAM-393
C INVERSE E-W TRANSFORMATION (P.A.M.,PRC 12(1975)744).                  SCAM-394
      IF (ICP.GT.NJC) GO TO 48                                          SCAM-395
      DO 47 IA=1,NCIN                                                   SCAM-396
      DO 46 IB=1,NJC                                                    SCAM-397
      A1=P(IA,IC,3)*P(IB,ICP,3)-P(IA,IC,4)*P(IB,ICP,4)                  SCAM-398
      B1=P(IA,IC,3)*P(IB,ICP,4)+P(IA,IC,4)*P(IB,ICP,3)                  SCAM-399
      A2=A1                                                             SCAM-400
      B2=B1                                                             SCAM-401
      IF (ICP.EQ.IC) GO TO 45                                           SCAM-402
      A2=A1+P(IA,ICP,3)*P(IB,IC,3)-P(IA,ICP,4)*P(IB,IC,4)               SCAM-403
      B2=B1+P(IA,ICP,3)*P(IB,IC,4)+P(IA,ICP,4)*P(IB,IC,3)               SCAM-404
   45 V(IA+4,IB)=V(IA+4,IB)+(A1*A2+B1*B2)*TQ                            SCAM-405
   46 CONTINUE                                                          SCAM-406
   47 CONTINUE                                                          SCAM-407
      GO TO 51                                                          SCAM-408
   48 DO 49 IA=1,NCIN                                                   SCAM-409
   49 V(IA+4,ICP)=V(IA+4,ICP)+(P(IA,IC,3)**2+P(IA,IC,4)**2)*TQ          SCAM-410
      GO TO 51                                                          SCAM-411
   50 V(5,ICP)=V(5,ICP)+AR*TQ                                           SCAM-412
   51 IF ((IC.GT.NCIN).OR.(.NOT.LO(55))) GO TO 52                       SCAM-413
      GNU=2.D0*V(4,ICP)                                                 SCAM-414
      L=IDINT(DMOD(V(1,ICP)/131072.D0,65536.D0)+.01D0)                  SCAM-415
      FJ=.5D0*DMOD(V(1,ICP),131072.D0)                                  SCAM-416
      WRITE (MW,1011) IC,ICP,IV,L,FJ,V(3,ICP),TQ,GNU,G                  SCAM-417
   52 CONTINUE                                                          SCAM-418
      IF (LO(83).AND.IC.GT.NCIN) GO TO 56                               SCAM-419
      DO 55 I=1,2                                                       SCAM-420
      IF (LG(I)) GO TO 55                                               SCAM-421
      IF (LO(87)) GO TO 54                                              SCAM-422
      G=0.D0                                                            SCAM-423
      DO 53 M=1,20                                                      SCAM-424
   53 G=G+W(M)/(TP(M)*(1.D0+X(M)*V(3,IC)/(V(4,IC)*BRI))*(1.D0+X(M)*TG(I)SCAM-425
     1/(BRI*FNU(I))))                                                   SCAM-426
   54 A1=V(3,IC)*TG(I)*G/BRI                                            SCAM-427
      SGF(I)=SGF(I)+A1*AR                                               SCAM-428
      FN=2.D0*FNU(I)                                                    SCAM-429
      IF (LO(55).AND.(IC.LE.NCIN)) WRITE (MW,1012) AL(I),TG(I),SGF(I),FNSCAM-430
     1,G                                                                SCAM-431
   55 CONTINUE                                                          SCAM-432
   56 CONTINUE                                                          SCAM-433
      IF (LO(55)) WRITE (MW,1013) BIR                                   SCAM-434
      DO 57 I=1,2                                                       SCAM-435
      IF (LG(I)) GO TO 57                                               SCAM-436
      IF (LO(132)) FCN(IPJ,NDQ+I)=FCN(IPJ,NDQ+I)+.25D0*SGF(I)           SCAM-437
      TX(NDP+I)=TX(NDP+I)+RZ*SGF(I)                                     SCAM-438
   57 CONTINUE                                                          SCAM-439
      GO TO 62                                                          SCAM-440
C STORAGE FOR SIMPLIFIED COMPOUND NUCLEUS.                              SCAM-441
   58 BIR=BZ(3)**2*(2.D0*AJ+1.D0)*DEXP(-AJ*(AJ+1.D0)/(2.D0*BZ(2)**2))   SCAM-442
      DO 59 IC=1,NJC                                                    SCAM-443
   59 BIR=BIR+4.D0*V(3,IC)                                              SCAM-444
      IF (LO(55)) WRITE (MW,1014) BIR                                   SCAM-445
      DO 61 IC=1,NCIN                                                   SCAM-446
      NT=NS+IK(1,2*MC(IC,2),MC(IC,3))*NCT(NZ+1)                         SCAM-447
      DO 60 ICP=1,NJC                                                   SCAM-448
      TQ=V(3,IC)*V(3,ICP)/BIR                                           SCAM-449
      IF (IV.EQ.1) TQ=TQ*BZ(1)**2                                       SCAM-450
   60 V(IC+4,ICP)=V(IC+4,ICP)+4.D0*TQ                                   SCAM-451
   61 CONTINUE                                                          SCAM-452
   62 IF (LO(55).AND.LO(123)) WRITE (MW,1000) AJ,IP(JPI+1),NC,NCIN      SCAM-453
C STORAGE OF THE AMPLITUDES.                                            SCAM-454
      NI=NCI-NCIN                                                       SCAM-455
      DO 70 I1=1,NCI                                                    SCAM-456
      NT=NS+IK(1,2*MD(I1,2),MD(I1,3))*NCT(NZ+1)                         SCAM-457
      IF (LO(123)) GO TO 63                                             SCAM-458
      J1=I1                                                             SCAM-459
      K1=J1                                                             SCAM-460
      GO TO 64                                                          SCAM-461
   63 J1=1                                                              SCAM-462
      K1=NCIN                                                           SCAM-463
   64 IVS=0                                                             SCAM-464
      DO 69 I2=1,NIC                                                    SCAM-465
      IV=MD(I2,1)                                                       SCAM-466
      IF ((IV.NE.IVS).AND.(IV.NE.1)) NT=NT+(IPI(2,IVS)*IPI(3,IVS)+NSY)/2SCAM-467
      NSY=MOD(IPI(1,IV)+JPI+IPJ+(IPI(2,IV)+IPI(3,IV))/2+1,2)            SCAM-468
      IVS=IV                                                            SCAM-469
      ID=NT+IK(IV,2*MD(I2,2),MD(I2,3))+1                                SCAM-470
      IF (LO(123).AND.(I2.LE.NCI)) GO TO 65                             SCAM-471
      J2=I2-NI                                                          SCAM-472
      K2=J2                                                             SCAM-473
      GO TO 66                                                          SCAM-474
   65 J2=1                                                              SCAM-475
      K2=NCIN                                                           SCAM-476
   66 AR=0.D0                                                           SCAM-477
      AI=0.D0                                                           SCAM-478
      DO 68 L1=J1,K1                                                    SCAM-479
      A1=DZ                                                             SCAM-480
      IF (LO(123)) A1=A1*AA(I1,L1)                                      SCAM-481
      IF (A1.EQ.0.D0) GO TO 68                                          SCAM-482
      DO 67 L2=J2,K2                                                    SCAM-483
      A2=1.D0                                                           SCAM-484
      IF (LO(123).AND.(L2.LE.NCIN)) A2=AA(I2,L2)                        SCAM-485
      IF (A2.EQ.0.D0) GO TO 67                                          SCAM-486
      AR=AR+FAR(L2,L1)*A1*A2                                            SCAM-487
      AI=AI+FAI(L2,L1)*A1*A2                                            SCAM-488
   67 CONTINUE                                                          SCAM-489
   68 CONTINUE                                                          SCAM-490
      F(1,IPJ,ID)=AR                                                    SCAM-491
      F(2,IPJ,ID)=AI                                                    SCAM-492
      IF (.NOT.(LO(55).AND.LO(123))) GO TO 69                           SCAM-493
      A1=DSQRT(AR*AR+AI*AI)                                             SCAM-494
      A2=0.D0                                                           SCAM-495
      IF (A1.NE.0.D0) A2=DATAN2(AI,AR)                                  SCAM-496
      FJ=0.5D0*DFLOAT(MD(I2,3))                                         SCAM-497
      WRITE (MW,1002) I2,I1,IV,MD(I2,2),FJ,AR,AI,A1,A2                  SCAM-498
   69 CONTINUE                                                          SCAM-499
   70 CONTINUE                                                          SCAM-500
      IF (LO(131)) RETURN                                               SCAM-501
C COMPOUND NUCLEUS.                                                     SCAM-502
      JCX=JC+NI                                                         SCAM-503
      IIV=IPI(2,1)-1                                                    SCAM-504
      III=IPI(3,1)-1                                                    SCAM-505
      DO 81 IC=1,NCI                                                    SCAM-506
      LA=2*MD(IC,2)                                                     SCAM-507
      JA=MD(IC,3)                                                       SCAM-508
      NT=NSS+IK(IV,LA,JA)*NCT(NZ+5)                                     SCAM-509
      IVS=0                                                             SCAM-510
      CALL COCN(LA,LA,JA,JA,III,IIV,NAJ,IPJ,Q,Q(IPJ+1),IDT-2*IPJ)       SCAM-511
      IF (LO(123)) GO TO 71                                             SCAM-512
      J1=IC                                                             SCAM-513
      K1=J1                                                             SCAM-514
      GO TO 72                                                          SCAM-515
   71 J1=1                                                              SCAM-516
      K1=NCIN                                                           SCAM-517
   72 DO 80 ICP=1,JCX                                                   SCAM-518
      IF (LO(123).AND.(ICP.LE.NCI)) GO TO 73                            SCAM-519
      J2=ICP-NI                                                         SCAM-520
      K2=J2                                                             SCAM-521
      GO TO 74                                                          SCAM-522
   73 J2=1                                                              SCAM-523
      K2=NCIN                                                           SCAM-524
   74 AR=0.D0                                                           SCAM-525
      DO 76 L1=J1,K1                                                    SCAM-526
      A1=DZ                                                             SCAM-527
      IF (LO(123)) A1=A1*AA(I1,L1)                                      SCAM-528
      IF (A1.EQ.0.D0) GO TO 79                                          SCAM-529
      DO 75 L2=J2,K2                                                    SCAM-530
      A2=1.D0                                                           SCAM-531
      IF (LO(123).AND.(L2.GT.NCIN)) A2=AA(I2,L2)                        SCAM-532
      IF (A2.EQ.0.D0) GO TO 75                                          SCAM-533
      AR=AR+V(L1+4,L2)*A1*A2                                            SCAM-534
   75 CONTINUE                                                          SCAM-535
   76 CONTINUE                                                          SCAM-536
      IF (AR.EQ.0.D0) GO TO 81                                          SCAM-537
      IV=IDINT(V(1,ICP)/131072.D0/65536.D0+.01D0)                       SCAM-538
      IF ((IV.NE.IVS).AND.(IV.NE.1)) NT=NT+(IPI(2,IVS)*IPI(3,IVS)+NSY)/2SCAM-539
      NSY=MOD(IPI(1,IV)+JPI+IPJ+(IPI(2,IV)+IPI(3,IV))/2+1,2)            SCAM-540
      IVS=IV                                                            SCAM-541
      TX(NCOLL+IV+1)=TX(NCOLL+IV+1)+AR*V(2,ICP)*RZ                      SCAM-542
      IF (ICP.GT.NC1) GO TO 79                                          SCAM-543
      JB=DMOD(V(1,ICP),131072.D0)                                       SCAM-544
c     LB=IDINT(V(1,ICP)/512.D0)-2048*IV                                 SCAM-545
      N=IDINT(V(1,ICP)/131072.D0+.01D0)                                 SCAM-546
      LB=2*MOD(N,65536)                                                 SCAM-547
c     MCX(J,2)=MOD(N,65536)                                             SCAM-548
c   6 MCX(J,1)=N/65536                                                  SCAM-549
c     JB=IDINT(V(1,ICP))-1048576*IV-512*LB                              SCAM-550
      IF (.NOT.LO(132)) GO TO 77                                        SCAM-551
      ID=NT+IK(IV,LB,JB)+1                                              SCAM-552
      FCN(IPJ,ID)=FCN(IPJ,ID)+.25D0*AR                                  SCAM-553
   77 CALL COCN(LB,LB,JB,JB,IPI(3,IV)-1,IPI(2,IV)-1,NAJ,IPJ,Q(IPJ+1),Q(2SCAM-554
     1*IPJ+1),IDT-3*IPJ)                                                SCAM-555
      IF ((IV.EQ.1).AND.(III.NE.0).AND.(ICP.NE.IC)) CALL COCN(LA,LB,JA,JSCAM-556
     1B,IPI(3,IV)-1,IPI(2,IV)-1,NAJ,IPJ,Q(2*IPJ+1),Q(3*IPJ+1),IDT-4*IPJ)SCAM-557
      DO 78 LL=1,IPJ                                                    SCAM-558
      AI=Q(LL)*Q(LL+IPJ)                                                SCAM-559
      IF ((IV.EQ.1).AND.(III.NE.0).AND.(ICP.NE.IC)) AI=AI+Q(LL+2*IPJ)**2SCAM-560
   78 GCN(LL,IV)=GCN(LL,IV)+.25D0*AI*AR*XZ                              SCAM-561
      GO TO 80                                                          SCAM-562
   79 IF (.NOT.LO(132)) GO TO 80                                        SCAM-563
      J=IV+KBA-NCOLS                                                    SCAM-564
      FCN(IPJ,J)=FCN(IPJ,J)+0.25D0*AR*V(2,ICP)                          SCAM-565
   80 CONTINUE                                                          SCAM-566
   81 JCX=NC1+NI                                                        SCAM-567
      RETURN                                                            SCAM-568
 1000 FORMAT (//' CHANNEL SPIN AND PARITY =',F7.1,A1,I11,' COUPLED CHANNSCAM-569
     1ELS AND',I3,' SOLUTIONS'//'  IC  ICP N    L    J',19X,'C MATRIX',2SCAM-570
     20X,'|C|',6X,'PHASE')                                              SCAM-571
 1001 FORMAT (' C MATRIX LARGER THAN 1 FOR AJ IPI NC NCIN IC ICP =',F5.1SCAM-572
     1,A1,4I3,2X,2D12.4)                                                SCAM-573
 1002 FORMAT (1X,I3,I4,I3,I5,F7.1,4X,1P,2D15.7,' I',4X,0P,2F11.8)       SCAM-574
 1003 FORMAT (1X,I3,I4,I3,I5,F7.1,4X,1P,2D15.7,' I',4X,'CLOSED')        SCAM-575
 1004 FORMAT (/' TRANSMISSION COEFFICIENTS FOR CHANNEL SPIN AND PARITY =SCAM-576
     1',F7.1,A1/4(' LEVEL   L     J',16X))                              SCAM-577
 1005 FORMAT (1X,F9.1,4X,A1,1X,I4)                                      SCAM-578
 1006 FORMAT (1X,I2,I6,F9.1,2X,1P,D18.8,0P)                             SCAM-579
 1007 FORMAT (4(1X,I2,I4,F6.1,2X,1P,D14.7,0P,3X))                       SCAM-580
 1008 FORMAT (' ERROR IN EIGENSYSTEM.  E-W CORRECTION DISCONTINUED.')   SCAM-581
 1009 FORMAT (//' CHANNEL SPIN AND PARITY =',F7.1, A1//'  IC  ICP N    LSCAM-582
     1    J',9X,'TL',9X,'HF',11X,'NU',9X,'G')                           SCAM-583
 1010 FORMAT (' TL  IN TRANSFORMED CHANNEL SPACE:')                     SCAM-584
 1011 FORMAT (1X,I3,I4,I3,I5,F7.1,F12.6,D14.6,2F10.6)                   SCAM-585
 1012 FORMAT (15X,A8,F12.6,D14.6,F10.4,F10.6)                           SCAM-586
 1013 FORMAT (20X,'SUM',F12.6)                                          SCAM-587
 1014 FORMAT (15X,'H.-F.',D20.6)                                        SCAM-588
      END                                                               SCAM-589
