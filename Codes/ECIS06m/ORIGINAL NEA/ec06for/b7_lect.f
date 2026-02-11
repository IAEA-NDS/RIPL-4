C 01/10/06                                                      ECIS06  LECT-000
      SUBROUTINE LECT(NBET,IPI,IPP,WV,RM,IDT,LO,NBETA,BETA,VAL,FISS,GAM,LECT-001
     1ISM,SCN)                                                          LECT-002
C LECT READS ALL INPUT EXCEPT THE 5 FIRST DATA CARDS,THE LEVELS, THE    LECT-003
C EXPERIMENTAL DATA,THE SEARCH CONDITIONS,THE NUCLEAR REDUCED MATRIX    LECT-004
C ELEMENTS AND THE FORM FACTORS.                                        LECT-005
C INPUT:     NBET:    NUMBER OF PHONONS.                                LECT-006
C            IPI(J,*):MULTIPLICITY OF PARTICLE FOR J=2.                 LECT-007
C                     PRODUCT OF CHARGES FOR J=4.                       LECT-008
C            IPI(1,J):FIRST LEVEL USING POTENTIAL J.                    LECT-009
C            WV(J,*): MASS OF THE PARTICLE FOR J=1.                     LECT-010
C                     MASS OF THE TARGET FOR J=2.                       LECT-011
C                     ENERGY IN THE CENTRE OF MASS IN MEV FOR J=3.      LECT-012
C                     ENERGY IN THE LABORATORY SYSTEM IN MEV FOR J=13.  LECT-013
C            RM:      MATCHING RADIUS (VALUE READ).                     LECT-014
C            IDT:     LENGTH FREE FOR BETA.                             LECT-015
C            LO(I):   LOGICAL CONTROLS:                                 LECT-016
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).LECT-017
C               LO(3)  =.TRUE. ANHARMONIC VIBRATIONAL OR ASYMMETRIC     LECT-018
C                              ROTATIONAL MODEL.                        LECT-019
C               LO(4)  =.TRUE. PARAMETRISED SPIN-ORBIT DEFORMATION.     LECT-020
C               LO(5)  =.TRUE  DIFFERENT DEFORMATION FOR EACH POTENTIAL.LECT-021
C               LO(6)  =.TRUE. USE DEFORMATION LENGTHS.                 LECT-022
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    LECT-023
C               LO(8)  =.TRUE. RELATIVISTIC KINEMATICS.                 LECT-024
C               LO(11) =.TRUE. DEFORMED COULOMB POTENTIAL.              LECT-025
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            LECT-026
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      LECT-027
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. LECT-028
C               LO(16) =.TRUE. HEAVY-ION DEFINITION OF REDUCED RADII ANDLECT-029
C                              DEFORMATIONS.                            LECT-030
C               LO(17) =.TRUE. FOLDING MODEL.                           LECT-031
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   LECT-032
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     LECT-033
C               LO(46) =.TRUE. RESTRICTED COULOMB CORRECTIONS.          LECT-034
C               LO(66) =.TRUE. NO CALCULATION AT EQUIDISTANT ANGLES.    LECT-035
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             LECT-036
C               LO(82) =.TRUE. OLD SIMPLIFIED COMPOUND NUCLEUS.         LECT-037
C               LO(85) =.TRUE. FISSION TRANSMISSION COEFFICIENTS.       LECT-038
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      LECT-039
C               LO(99) =.TRUE. SCHROEDINGER EQUIVALENT TO DIRAC         LECT-040
C                              EQUATION.                                LECT-041
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    LECT-042
C OUTPUT:    NBET:    NUMBER OF PHONONS PLUS DEFORMATIONS.              LECT-043
C            RM:      MATCHING RADIUS (ACTUAL VALUE).                   LECT-044
C            NBETA:   QUANTUM NUMBERS OF DEFORMATIONS IN NBETA(J,*) FOR LECT-045
C                     J=17,18.                                          LECT-046
C            BETA:    IN EQUIVALENCE BY CALL WITH NBETA: NUCLEAR        LECT-047
C                     DEFORMATION FOR POTENTIALS IN BETA(J,*) FOR J=1,8.LECT-048
C            VAL:     (25 OPTICAL +9 FOLDING PARAMETERS)* NPP           LECT-049
C            FISS:    FISSION DATA FOR COMPOUND NUCLEUS                 LECT-050
C            GAM:     GAMMA DATA FOR COMPOUND NUCLEUS                   LECT-051
C            ISM:     NUMBER OF INTEGRATION STEPS                       LECT-052
C            SCN:     DESCRIPTIONS OF LEVEL DENSITIES                   LECT-053
C            LO:      LOGICAL CONTROLS:                                 LECT-054
C               LO(11) =.TRUE. DEFORMED COULOMB POTENTIAL.              LECT-055
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      LECT-056
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. LECT-057
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   LECT-058
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      LECT-059
C               LO(101)=.TRUE. THERE IS A REAL SPIN-ORBIT POTENTIAL.    LECT-060
C               LO(102)=.TRUE. THERE IS AN IMAGINARY SPIN-ORBIT         LECT-061
C                              POTENTIAL.                               LECT-062
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. LECT-063
C               LO(108)=.TRUE. DIAGONAL COULOMB CORRECTIONS ARE NEEDED. LECT-064
C               LO(129)=.TRUE. REAL SPIN-ORBIT OR DIRAC EQUATION.       LECT-065
C               LO(130)=.TRUE. IMAGINARY SPIN-ORBIT OR DIRAC EQUATION.  LECT-066
C                                                                       LECT-067
C THE COMMON /ANGUL/ IS USED IN CALX, LECT, SCHE, LCSP, RESU AND REST.  LECT-068
C                                                                       LECT-069
C FOR THE COMMON  /CONVE/, /COUPL/, /DCONS/ AND /NCOMP/ SEE CALC.       LECT-070
C                                                                       LECT-071
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ANGUL/:                     LECT-072
C  THETA1:    FIRST ANGLE FOR COMPUTATION AT EQUIDISTANT ANGLES.        LECT-073
C  DTHETA:    STEP FOR COMPUTATION AT EQUIDISTANT ANGLES.               LECT-074
C  THETA2:    LAST ANGLE FOR COMPUTATION AT EQUIDISTANT ANGLES.         LECT-075
C  DTHE:      AVERAGING ANGLE.                                          LECT-076
C  NCJ:       NUMBER OF FACTORISATIONS OF 1/(1-COS(THETA)) IN AMPLITUDE.LECT-077
C  NL(1):     POWER OF (1-COS(THETA)) FOR THE EXPANSION IN LEGENDRE     LECT-078
C             POLYNOMIALS OF THE INTERFERENCE BETWEEN COULOMB AND       LECT-079
C             NUCLEAR ELASTIC SCATTERING. POWER OF (1-COS(THETA)**2)    LECT-080
C             IF LO(18) IS .TRUE..                                      LECT-081
C  NL(2):     NUMBER OF LEGENDRE POLYNOMIALS FOR ELASTIC SCATTERING OF  LECT-082
C             CHARGED PARTICLES.                                        LECT-083
C  NL(3):     NUMBER OF LEGENDRE POLYNOMIALS FOR ELASTIC SCATTERING OF  LECT-084
C             UNCHARGED PARTICLES, INELASTIC SCATTERING AND COMPOUND    LECT-085
C             NUCLEUS.                                                  LECT-086
C  JMM(1):    NUMBER OF CHANNEL SPINS USED FOR MINIMUM CHI2.            LECT-087
C  JMM(2):    NUMBER OF L FOR COMPOUND NUCLEUS FOR MINIMUM CHI2.        LECT-088
C   DEFINED:  THETA1,DTHETA,THETA2,DTHE                                 LECT-089
C   USED:     THETA1,DTHETA,THETA2,DTHE                                 LECT-090
C   NOT USED: NCJ,NL,JMM.                                               LECT-091
C                                                                       LECT-092
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     LECT-093
C  H:         STEP SIZE FOR INTEGRATION.                                LECT-094
C  ACONV:     CONVERGENCE CRITERION FOR POTENTIAL AND FUNCTION.         LECT-095
C   DEFINED:  H,ACONV.                                                  LECT-096
C   USED:     H,ACONV.                                                  LECT-097
C                                                                       LECT-098
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     LECT-099
C  IQM:       MAXIMUM L-VALUE OF DEFORMATION IN ROTATIONAL MODEL.       LECT-100
C  IQMAX:     MAXIMUM L-EXPANSION IN ROTATIONAL MODEL.                  LECT-101
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             LECT-102
C  NSPIN:     TWICE THE K-VALUE OF THE ROTATIONAL BAND.                 LECT-103
C   DEFINED:  IQM,IQMAX,NSPIN.                                          LECT-104
C   USED:     NPP.                                                      LECT-105
C                                                                       LECT-106
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     LECT-107
C  CM:        ATOMIC MASS IN MEV.                                       LECT-108
C  CK:        H-BAR*C.                                                  LECT-109
C   USED:     CM,CK.                                                    LECT-110
C                                                                       LECT-111
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     LECT-112
C  NFISS:     NUMBER OF FISSION TRANSMISSION COEFFICIENTS.              LECT-113
C  NRD:       NUMBER OF GAMMA TRANSMISSION COEFFICIENTS.                LECT-114
C  NCONT:     NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS.                 LECT-115
C  NCONS:     NUMBER OF LEVEL DENSITIES NEEDED.                         LECT-116
C  NCOLX:     TOTAL NUMBER OF LEVELS WITHOUT DISCRETISATION.            LECT-117
C  AZ(6):     DEFORMED SPIN-ORBIT PARAMETERS. SEE ALSO COMMENT IN       LECT-118
C             INPUT DESCRIPTION AND SUBROUTINE QUAN.                    LECT-119
C  BZ(6)      HAUSER-FESHBACH AND MOLDAUER'S PARAMETERS DESCRIBED BELOW.LECT-120
C   BZ(1):    SQUARE ROOT OF ELASTIC ENHANCEMENT.                       LECT-121
C   BZ(2):    IF LO(82)=.TRUE., SPIN CUT-OFF PARAMETER.                 LECT-122
C             IF LO(82)=.FALSE., PARTICLE DEGREES OF FREEDOM.           LECT-123
C   BZ(3):    SQUARE ROOT OF LEVEL DENSITY PARAMETER. IF LO(82)=LO(87)= LECT-124
C             .FALSE., PARAMETER BZ(3) IN MOLDAUER'S FORMULA OF INPUT   LECT-125
C             DESCRIPTION.                                              LECT-126
C   BZ(4):    IF LO(82)=LO(87)=.FALSE., PARAMETER BZ(4) IN SAME FORMULA.LECT-127
C   BZ(5):    IF LO(82)=LO(87)=.FALSE., PARAMETER BZ(5) IN SAME FORMULA.LECT-128
C  TG0:       SLOW S-WAVE NEUTRON GAMMA WIDTH/SPACING FOR NORMALISATION.LECT-129
C  BN:        NEUTRON SEPARATION ENERGY.                                LECT-130
C  FNUG:      RADIATIVE DEGREES OF FREEDOM.                             LECT-131
C  EGD:       ENERGY OF THE GIANT DIPOLE RESONANCE.                     LECT-132
C  GGD:       RESONANCE WIDTH.                                          LECT-133
C   DEFINED:  AZ,BZ,TG0,BN,FNUG,EGD,GGD.                                LECT-134
C   USED:     NFISS,NRD,NCONT,NCONS,NCOLX,AZ,BZ,TG0,BN,FNUG,EGD,GGD.    LECT-135
C                                                                       LECT-136
C***********************************************************************LECT-137
      IMPLICIT REAL*8 (A-H,O-Z)                                         LECT-138
      DIMENSION IPI(11,*),IPP(34,*),WV(22,*),NBETA(18,*),BETA(9,*),VAL(4LECT-139
     12,*),FISS(2,*),GAM(*),SCN(7,*),RO(6)                              LECT-140
      CHARACTER*4 ALEG(6),BLEG(2)                                       LECT-141
      LOGICAL LO(150),LT(6)                                             LECT-142
      COMMON /ANGUL/ THETA1,THETA2,DTHETA,DTHE,NCJ,NL(3),JMM(2)         LECT-143
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       LECT-144
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   LECT-145
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            LECT-146
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQLECT-147
     1,ACN1,ACN2,AZ(6),BZ(5),TG0,BN,FNUG,EGD,GGD,TG1,SGSQ               LECT-148
      COMMON /INOUT/ MR,MW,MS                                           LECT-149
      DATA BLEG /'NOT ','    '/                                         LECT-150
      IQM=0                                                             LECT-151
      IF (NBET.EQ.0) GO TO 4                                            LECT-152
C INPUT OF PHONON DEFORMATIONS IN VIBRATIONAL AND ROTATIONAL MODELS     LECT-153
      DO 3 I=1,NBET                                                     LECT-154
      READ (MR,1000) NBETA(17,I),NBETA(18,I),BETA(1,I),K                LECT-155
      DO 1 J=2,8                                                        LECT-156
    1 BETA(J,I)=BETA(1,I)                                               LECT-157
      IF (LO(5)) READ (MR,1001) (BETA(J,I),J=2,8)                       LECT-158
      WRITE (MW,1002) I,NBETA(17,I),NBETA(18,I),(BETA(J,I),J=1,8)       LECT-159
      IF (.NOT.LO(16)) GO TO 3                                          LECT-160
      IF (K.EQ.0) K=1                                                   LECT-161
      K=MIN0(K,NCOLX)                                                   LECT-162
      AM3=WV(2,K)**.33333333333333D0                                    LECT-163
      BM3=WV(1,K)**.33333333333333D0                                    LECT-164
      CM3=AM3/(AM3+BM3)                                                 LECT-165
      DM3=1.D0                                                          LECT-166
      IF (LO(6)) DM3=CM3                                                LECT-167
      DO 2 J=1,6                                                        LECT-168
    2 BETA(J,I)=BETA(J,I)*CM3/DM3                                       LECT-169
      BETA(7,I)=BETA(7,I)*CM3**NBETA(17,I)/DM3                          LECT-170
      BETA(8,I)=BETA(8,I)*CM3**NBETA(17,I)/DM3                          LECT-171
      WRITE (MW,1003) K,(BETA(J,I),J=1,8)                               LECT-172
    3 CONTINUE                                                          LECT-173
    4 IF (LO(7)) GO TO 18                                               LECT-174
      IF (.NOT.LO(1)) GO TO 13                                          LECT-175
C INPUT OF DEFORMATIONS FOR THE ROTATIONAL MODEL                        LECT-176
      READ (MR,1000) IQM,IQMAX,ASPIN,K                                  LECT-177
      NSPIN=IDINT(2.D0*ASPIN+0.1D0)                                     LECT-178
      IQ=NBET+1                                                         LECT-179
      JQ=NBET+IQM/2                                                     LECT-180
      IF (.NOT.LO(3)) GO TO 5                                           LECT-181
      JQ=NBET+IQM                                                       LECT-182
      IF ((IQM.GT.35).OR.(IQMAX.GT.8)) GO TO 33                         LECT-183
    5 IF (9*JQ.GT.IDT) CALL MEMO('LECT',IDT,9*JQ)                       LECT-184
      READ (MR,1001) (BETA(1,I),I=IQ,JQ)                                LECT-185
      M=0                                                               LECT-186
      L=0                                                               LECT-187
      DO 9 I=IQ,JQ                                                      LECT-188
      IF (LO(3)) GO TO 6                                                LECT-189
      L=2*(I-NBET)                                                      LECT-190
      GO TO 7                                                           LECT-191
    6 M=M+2                                                             LECT-192
      IF (M.LE.L) GO TO 7                                               LECT-193
      L=L+2                                                             LECT-194
      M=0                                                               LECT-195
    7 NBETA(17,I)=L                                                     LECT-196
      NBETA(18,I)=M                                                     LECT-197
      DO 8 J=2,8                                                        LECT-198
    8 BETA(J,I)=BETA(1,I)                                               LECT-199
    9 CONTINUE                                                          LECT-200
      IF (LO(5)) READ (MR,1001) ((BETA(J,I),J=2,8),I=IQ,JQ)             LECT-201
      WRITE (MW,1004) IQMAX,ASPIN                                       LECT-202
      WRITE (MW,1005) (I,NBETA(17,I),NBETA(18,I),(BETA(J,I),J=1,8),I=IQ,LECT-203
     1JQ)                                                               LECT-204
      IF (.NOT.LO(16)) GO TO 12                                         LECT-205
      IF (K.EQ.0) K=1                                                   LECT-206
      K=MIN0(K,NCOLX)                                                   LECT-207
      AM3=WV(2,K)**.33333333333333D0                                    LECT-208
      BM3=WV(1,K)**.33333333333333D0                                    LECT-209
      CM3=AM3/(AM3+BM3)                                                 LECT-210
      DM3=1.D0                                                          LECT-211
      IF (LO(6)) DM3=CM3                                                LECT-212
      DO 11 I=IQ,JQ                                                     LECT-213
      IF (NBETA(18,I).NE.0) GO TO 11                                    LECT-214
      DO 10 K=1,6                                                       LECT-215
   10 BETA(K,I)=BETA(K,I)*CM3/DM3                                       LECT-216
      BETA(7,I)=BETA(7,I)*CM3**NBETA(17,I)/DM3                          LECT-217
      BETA(8,I)=BETA(8,I)*CM3**NBETA(17,I)/DM3                          LECT-218
   11 CONTINUE                                                          LECT-219
      WRITE (MW,1006) K                                                 LECT-220
      WRITE (MW,1005) (I,NBETA(17,I),NBETA(18,I),(BETA(J,I),J=1,8),I=IQ,LECT-221
     1JQ)                                                               LECT-222
   12 NBET=JQ                                                           LECT-223
      GO TO 18                                                          LECT-224
C ANHARMONIC VIBRATIONAL MODEL WITH DIFFERENT DEFORMATIONS.             LECT-225
   13 IF (.NOT.LO(3)) GO TO 18                                          LECT-226
      IF (36.GT.IDT) CALL MEMO('LECT',IDT,36)                           LECT-227
      WRITE (MW,1007)                                                   LECT-228
      DO 15 I=1,4                                                       LECT-229
      DO 14 J=1,8                                                       LECT-230
   14 BETA(J,I)=1.D0                                                    LECT-231
      NBETA(17,I)=0                                                     LECT-232
      NBETA(18,I)=0                                                     LECT-233
      IF (LO(5)) READ (MR,1001) (BETA(J,I),J=1,8)                       LECT-234
   15 WRITE (MW,1008) I,NBETA(17,I),(BETA(J,I),J=1,8)                   LECT-235
      NBET=4                                                            LECT-236
      IF (.NOT.LO(16)) GO TO 18                                         LECT-237
      READ (MR,1009) (NBETA(17,J),J=1,4),K                              LECT-238
      IF (K.EQ.0) K=1                                                   LECT-239
      K=MIN0(K,NCOLX)                                                   LECT-240
      AM3=WV(2,K)**.33333333333333D0                                    LECT-241
      BM3=WV(1,K)**.33333333333333D0                                    LECT-242
      CM3=AM3/(AM3+BM3)                                                 LECT-243
      DM3=1.D0                                                          LECT-244
      IF (LO(6)) DM3=CM3                                                LECT-245
      WRITE (MW,1006) K                                                 LECT-246
      DO 17 I=2,4                                                       LECT-247
      IF (NBETA(17,I).EQ.0) NBETA(17,I)=2                               LECT-248
      IF (NBETA(17,I).LT.0) NBETA(17,I)=0                               LECT-249
      DO 16 J=1,6                                                       LECT-250
   16 BETA(J,I)=BETA(J,I)*(CM3/DM3)**(I-1)                              LECT-251
      BETA(7,I)=BETA(7,I)*(CM3**NBETA(17,I)/DM3)**(I-1)                 LECT-252
      BETA(8,I)=BETA(8,I)*(CM3**NBETA(17,I)/DM3)**(I-1)                 LECT-253
   17 WRITE (MW,1010) I,NBETA(17,I),(BETA(J,I),J=1,8)                   LECT-254
C INPUT OF OPTICAL MODEL PARAMETERS.                                    LECT-255
   18 W2=10.D20                                                         LECT-256
      LO(101)=.FALSE.                                                   LECT-257
      LO(102)=.FALSE.                                                   LECT-258
      LO(103)=.FALSE.                                                   LECT-259
      LO(108)=.FALSE.                                                   LECT-260
      LT(4)=.TRUE.                                                      LECT-261
      LT(6)=.TRUE.                                                      LECT-262
      W1=0.D0                                                           LECT-263
      W3=WV(4,1)/(ACONV*WV(3,1))                                        LECT-264
      DO 21 IP=1,NPP                                                    LECT-265
      IJ=IPP(1,IP)                                                      LECT-266
      IF (IJ.NE.-1) GO TO 19                                            LECT-267
      IPP(1,IP)=1                                                       LECT-268
      WRITE (MW,1011) IP                                                LECT-269
      IJ=1                                                              LECT-270
   19 IF (LO(7)) GO TO 21                                               LECT-271
      AM3=WV(2,IJ)**.33333333333333D0                                   LECT-272
      IF (LO(16)) AM3=AM3+WV(1,IJ)**.33333333333333D0                   LECT-273
      WRITE (MW,1012) IP,AM3                                            LECT-274
      DO 20 I=1,6                                                       LECT-275
      ALEG(I)=BLEG(2)                                                   LECT-276
      READ (MR,1013) VAL(4*I-3,IP),RO(I),VAL(4*I-1,IP),VAL(4*I,IP)      LECT-277
      VAL(4*I-2,IP)=AM3*RO(I)                                           LECT-278
      IF (VAL(4*I-3,IP).EQ.0.D0) GO TO 20                               LECT-279
      W1=DMAX1(W1,VAL(4*I-2,IP)+DLOG(W3*DABS(VAL(4*I-3,IP)))*VAL(4*I-1,ILECT-280
     1P))                                                               LECT-281
      W2=DMIN1(W2,VAL(4*I-1,IP))                                        LECT-282
   20 CONTINUE                                                          LECT-283
      LT(1)=IPI(2,IJ).EQ.1.OR.VAL(17,IP).EQ.0.D0                        LECT-284
      LT(2)=IPI(2,IJ).EQ.1.OR.VAL(21,IP).EQ.0.D0                        LECT-285
      LO(101)=LO(101).OR.(.NOT.LT(1))                                   LECT-286
      LO(102)=LO(102).OR.(.NOT.LT(2))                                   LECT-287
C OUTPUT OF OPTICAL PARAMETERS.                                         LECT-288
      IF (.NOT.LO(12)) ALEG(2)=BLEG(1)                                  LECT-289
      IF (LT(1).OR.(.NOT.LO(13))) ALEG(5)=BLEG(1)                       LECT-290
      IF (LT(2).OR.(.NOT.LO(14))) ALEG(6)=BLEG(1)                       LECT-291
      ALEG(4)=ALEG(2)                                                   LECT-292
      WRITE (MW,1014) (ALEG(I),VAL(4*I-3,IP),VAL(4*I-2,IP),RO(I),VAL(4*ILECT-293
     1-1,IP),VAL(4*I,IP),I=1,6)                                         LECT-294
      READ (MR,1013) RO(1),VAL(27,IP),VAL(33,IP),VAL(28,IP)             LECT-295
      READ (MR,1013) VAL(29,IP),RO(2),VAL(31,IP),VAL(33,IP)             LECT-296
      IF (LO(99)) VAL(27,IP)=0.D0                                       LECT-297
      IF (LO(99)) VAL(31,IP)=0.D0                                       LECT-298
      VAL(25,IP)=IPI(4,IJ)                                              LECT-299
      VAL(26,IP)=AM3*RO(1)                                              LECT-300
      LT(3)=VAL(25,IP).EQ.0.D0                                          LECT-301
      LT(4)=LT(3).AND.LT(4)                                             LECT-302
      ALEG(1)=BLEG(2)                                                   LECT-303
      IF (LT(3).OR.(.NOT.LO(11))) ALEG(1)=BLEG(1)                       LECT-304
      VAL(30,IP)=AM3*RO(2)                                              LECT-305
      LT(5)=VAL(29,IP).EQ.0.OR.IPI(2,IJ).EQ.1                           LECT-306
      LO(103)=LO(103).OR.(.NOT.LT(5))                                   LECT-307
      LT(6)=LT(5).AND.LT(6)                                             LECT-308
      ALEG(2)=BLEG(2)                                                   LECT-309
      IF (LT(5).OR.(.NOT.LO(19))) ALEG(2)=BLEG(1)                       LECT-310
      WRITE (MW,1015) (ALEG(I),VAL(4*I+21,IP),VAL(4*I+22,IP),RO(I),VAL(4LECT-311
     1*I+23,IP),VAL(4*I+24,IP),I=1,2),VAL(25,IP)                        LECT-312
      LO(108)=LO(108).OR.VAL(29,IP).NE.0.D0.OR.(LO(109).AND.VAL(19,IP).NLECT-313
     1E.0.D0)                                                           LECT-314
      W1=DMAX1(W1,VAL(26,IP)+10.D0*VAL(27,IP),VAL(30,IP)+10.D0*VAL(31,IPLECT-315
     1))                                                                LECT-316
      IF (VAL(27,IP).NE.0.D0) W2=DMIN1(W2,VAL(27,IP))                   LECT-317
      IF (VAL(33,IP).NE.0.D0) W2=DMIN1(W2,VAL(33,IP))                   LECT-318
      IF (.NOT.LO(17)) GO TO 21                                         LECT-319
C INPUT OF FOLDING PARAMETERS.                                          LECT-320
      READ (MR,1016) (VAL(I,IP),I=34,42)                                LECT-321
      WRITE (MW,1017) (VAL(I,IP),I=34,42)                               LECT-322
   21 CONTINUE                                                          LECT-323
      IF (LO(99)) LO(102)=.TRUE.                                        LECT-324
      IF (.NOT.LO(7)) GO TO 23                                          LECT-325
      DO 22 IJ=1,NCOLX                                                  LECT-326
      LT(4)=LT(4).AND.IPI(4,IJ).EQ.0                                    LECT-327
      LT(6)=LT(6).AND.IPI(2,IJ).EQ.1                                    LECT-328
   22 LO(102)=LO(102).OR.IPI(2,IJ).NE.1                                 LECT-329
      W1=20.D0                                                          LECT-330
      LO(101)=LO(102)                                                   LECT-331
      LO(103)=LO(102).AND.(.NOT.LT(6))                                  LECT-332
      LO(108)=LO(109).OR.LO(103)                                        LECT-333
   23 LO(101)=LO(101).OR.LO(102).OR.LO(103)                             LECT-334
      LO(129)=LO(101).OR.LO(109)                                        LECT-335
      LO(130)=LO(102).OR.LO(109)                                        LECT-336
      LO(13)=LO(13).AND.LO(101)                                         LECT-337
      LO(14)=LO(14).AND.LO(102)                                         LECT-338
      LO(11)=LO(11).AND.(.NOT.LT(4))                                    LECT-339
      LO(19)=LO(19).AND.(.NOT.LT(6))                                    LECT-340
      IF (LO(19)) LO(11)=.TRUE.                                         LECT-341
      LO(108)=LO(108).OR.LO(44)                                         LECT-342
      IF ((.NOT.LO(44)).AND.LO(46)) LO(108)=.FALSE.                     LECT-343
C DEFAULT VALUE OF MATCHING RADIUS: MAX(R*A**1/3+10*A).                 LECT-344
C DEFAULT VALUE FOR INTEGRATION STEP: MIN(MIN(A)/2,1/(2*K)).            LECT-345
      IF (RM.LE.0.D0) RM=W1                                             LECT-346
      IF (H.LE.0.D0) H=DMIN1(W2/2.D0,0.5D0/WV(4,1))                     LECT-347
      ISM=IDINT(RM/H+0.5D0)                                             LECT-348
      H=DFLOAT(ISM)                                                     LECT-349
      H=RM/H                                                            LECT-350
      RM=H*ISM                                                          LECT-351
      WRITE (MW,1018) H,RM                                              LECT-352
      IF (LO(66)) GO TO 24                                              LECT-353
C INPUT OF LIMITS FOR CALCULATION AT EQUIDISTANT ANGLES.                LECT-354
      READ (MR,1001) THETA1,DTHETA,THETA2,DTHE                          LECT-355
      IF (DTHETA.EQ.0.D0) DTHETA=1.D0                                   LECT-356
      IF ((THETA2-THETA1)*DTHETA.LT.0.D0) DTHETA=-DTHETA                LECT-357
      WRITE (MW,1019) THETA1,DTHETA,THETA2,DTHE                         LECT-358
C INPUT OF DEFORMED SPIN-ORBIT PARAMETERS.                              LECT-359
   24 AZ(1)=0.D0                                                        LECT-360
      AZ(2)=1.D0                                                        LECT-361
      AZ(3)=1.D0                                                        LECT-362
      AZ(4)=0.D0                                                        LECT-363
      AZ(5)=1.D0                                                        LECT-364
      AZ(6)=1.D0                                                        LECT-365
      IF (.NOT.LO(4)) GO TO 25                                          LECT-366
      READ (MR,1001) AZ                                                 LECT-367
      WRITE (MW,1020) AZ                                                LECT-368
   25 IF (.NOT.LO(81)) GO TO 32                                         LECT-369
      READ (MR,1001) BZ                                                 LECT-370
      IF (BZ(1).EQ.0.D0) BZ(1)=1.4142D0                                 LECT-371
      IF (.NOT.LO(82)) GO TO 26                                         LECT-372
      IF (BZ(2).EQ.0.D0) BZ(2)=3.5D0                                    LECT-373
      IF (BZ(3).EQ.0.D0) BZ(3)=100.D0                                   LECT-374
      WRITE (MW,1021) (BZ(I),I=1,3)                                     LECT-375
      GO TO 32                                                          LECT-376
   26 IF (BZ(3).EQ.0.D0) BZ(3)=1.212D0                                  LECT-377
      IF (BZ(4).EQ.0.D0) BZ(4)=0.78D0                                   LECT-378
      IF (BZ(5).EQ.0.D0) BZ(5)=0.228D0                                  LECT-379
      WRITE (MW,1022) (BZ(I),I=1,5)                                     LECT-380
      IF (.NOT.LO(85)) GO TO 28                                         LECT-381
      WRITE (MW,1023) NFISS                                             LECT-382
      DO 27 I=1,NFISS                                                   LECT-383
      READ (MR,1001) FISS(1,I),FISS(2,I)                                LECT-384
      IF (FISS(2,I).LT.0.5D0) FISS(2,I)=0.D0                            LECT-385
   27 WRITE (MW,1024) FISS(1,I),FISS(2,I)                               LECT-386
   28 IF (.NOT.LO(86)) GO TO 30                                         LECT-387
      WRITE (MW,1025)                                                   LECT-388
      IF (NRD.GT.0) GO TO 29                                            LECT-389
      READ (MR,1001) TG0,BN,FNUG,EGD,GGD                                LECT-390
      READ (MR,1001) SCN(7,1),(SCN(I,1),I=1,6)                          LECT-391
      IF (BN.EQ.0.D0) BN=8.D0                                           LECT-392
      IF (FNUG.LT.1.D0) FNUG=20.D0                                      LECT-393
      IF (TG0.EQ.0.D0) LO(86)=.FALSE.                                   LECT-394
      WRITE (MW,1026) TG0,BN,FNUG,EGD,GGD                               LECT-395
      J=1                                                               LECT-396
      WRITE (MW,1027) J,(SCN(I,J),I=1,7)                                LECT-397
      NA=IDINT(WV(2,1)+WV(1,1)+.5D0)                                    LECT-398
      CALL LDEN(NA,SCN(1,1))                                            LECT-399
      GO TO 30                                                          LECT-400
   29 READ (MR,1001) (GAM(I),I=1,NRD)                                   LECT-401
      WRITE (MW,1028) (GAM(I),I=1,NRD)                                  LECT-402
   30 IF (NCONT.EQ.0) GO TO 32                                          LECT-403
      WRITE (MW,1029)                                                   LECT-404
      J1=1+NCONS-NCONT                                                  LECT-405
      DO 31 J=J1,NCONS                                                  LECT-406
      READ (MR,1001) SCN(7,J),(SCN(I,J),I=1,6)                          LECT-407
      WRITE (MW,1027) J,(SCN(I,J),I=1,7)                                LECT-408
      NA=IDINT(WV(2,NCOLX+J-NCONS)+.5D0)                                LECT-409
   31 CALL LDEN(NA,SCN(1,J))                                            LECT-410
   32 RETURN                                                            LECT-411
   33 WRITE (MW,1030) IQM,IQMAX                                         LECT-412
      STOP                                                              LECT-413
 1000 FORMAT (2I5,F10.5,I5)                                             LECT-414
 1001 FORMAT (7F10.5)                                                   LECT-415
 1002 FORMAT (10X,'PHONON',I3,5X,'L =',I3,5X,'K =',I3,5X,'BETA =',8F8.5)LECT-416
 1003 FORMAT (12X,'AFTER HEAVY ION CORRECTION I =',I2,'  BETA =',8F8.5) LECT-417
 1004 FORMAT (/10X,'MULTIPOLE EXPANSION UP TO IQMAX =',I3,20X,'K =',F7.2LECT-418
     1,' BAND.')                                                        LECT-419
 1005 FORMAT (/' BETA(I,J) FOR  L   K ',8X,'V',9X,'W',8X,'VS',8X,'WS',7XLECT-420
     1,'VSO',7X,'WSO',6X,'COUL S.O. COUL'/(5X,I5,5X,I2,2X,I2,2X,8F10.5))LECT-421
 1006 FORMAT (/10X,'AFTER HEAVY ION CORRECTION I =',I3/)                LECT-422
 1007 FORMAT (/10X,'RATIOS OF ANHARMONIC DEFORMATIONS'/25X,'V',9X,'W',8XLECT-423
     1,'VS',8X,'WS',7X,'VSO',7X,'WSO',6X,'COUL S.O. COUL')              LECT-424
 1008 FORMAT (' ORDER',I3,6X,I3,2X,8F10.5)                              LECT-425
 1009 FORMAT (14I5)                                                     LECT-426
 1010 FORMAT (' ORDER',I3,'  IQ =',I3,2X,8F10.5)                        LECT-427
 1011 FORMAT (/' ***** NO STATE FOR THE POTENTIAL',I3,' *** WE USE THE GLECT-428
     1ROUND STATE'/)                                                    LECT-429
 1012 FORMAT (//' OPTICAL POTENTIALS  **',I3,' **     REDUCED RADIUS MULLECT-430
     1TIPLIED BY  ',1P,D15.6/)                                          LECT-431
 1013 FORMAT (4F10.5)                                                   LECT-432
 1014 FORMAT (2X,A4,'DEFORMED  VOLUME/SCALAR REAL POTENTIAL'/' DEPTH',F1LECT-433
     12.6,' MEV   RADIUS',F10.6,' FERMI (REDUCED VALUE',F9.6,')   DIFFUSLECT-434
     2ENESS',F9.6,' FERMI AT THE POWER (1+',F9.6,')'//2X,A4,'DEFORMED  VLECT-435
     3OLUME/SCALAR IMAGINARY POTENTIAL'/' DEPTH',F12.6,' MEV   RADIUS',FLECT-436
     410.6,' FERMI (REDUCED VALUE',F9.6,')   DIFFUSENESS',F9.6,' FERMI ALECT-437
     5T THE POWER (1+',F9.6,')'//2X,A4,'DEFORMED  SURFACE/VECTOR REAL POLECT-438
     6TENTIAL'/' DEPTH',F12.6,' MEV   RADIUS',F10.6,' FERMI (REDUCED VALLECT-439
     7UE',F9.6,')   DIFFUSENESS',F9.6,' FERMI AT THE POWER (1+',F9.6,')'LECT-440
     8//2X,A4,'DEFORMED  SURFACE/VECTOR IMAGINARY POTENTIAL'/' DEPTH',F1LECT-441
     92.6,' MEV   RADIUS',F10.6,' FERMI (REDUCED VALUE',F9.6,')   DIFFUSLECT-442
     AENESS',F9.6,' FERMI AT THE POWER (1+',F9.6,')'//2X,A4,'DEFORMED  RLECT-443
     BEAL SPIN-ORBIT/TENSOR POTENTIAL'/' DEPTH',F12.6,' MEV   RADIUS',F1LECT-444
     C0.6,' FERMI (REDUCED VALUE',F9.6,')   DIFFUSENESS',F9.6,' FERMI ATLECT-445
     D THE POWER (1+',F9.6,')'//2X,A4,'DEFORMED  IMAGINARY SPIN-ORBIT/TELECT-446
     ENSOR POTENTIAL'/' DEPTH',F12.6,' MEV   RADIUS',F10.6,' FERMI (REDULECT-447
     FCED VALUE',F9.6,')   DIFFUSENESS',F9.6,' FERMI AT THE POWER (1+',FLECT-448
     H9.6,')')                                                          LECT-449
 1015 FORMAT (/2X,A4,'DEFORMED  COULOMB POTENTIAL'/'  CHARGE PRODUCT',F7LECT-450
     1.0,'  RADIUS',F10.6,' FERMI (REDUCED VALUE',F9.6,')   DIFFUSENESS'LECT-451
     2,F9.6,' FERMI AT THE POWER (1+',F9.6,')'//2X,A4,'DEFORMED  SPIN-ORLECT-452
     3BIT COULOMB POTENTIAL'/' DEPTH',F12.6,' MEV   RADIUS',F10.6,' FERMLECT-453
     4I (REDUCED VALUE',F9.6,')   DIFFUSENESS',F9.6,' FERMI AT THE POWERLECT-454
     5 (1+',F9.6,')'/14X,'THIRD CHARGE PARAMETER',F9.6)                 LECT-455
 1016 FORMAT (3F10.5)                                                   LECT-456
 1017 FORMAT (/' *** FOLDING MODEL ***'/' REAL PART      V =',F10.4,6X,'LECT-457
     1R =',F10.4,6X,'A =',F10.4/' IMAGINARY PART V =',F10.4,6X,'R =',F10LECT-458
     2.4,6X,'A =',F10.4/' COULOMB PART   V =',F10.4,6X,'R =',F10.4,6X,'ALECT-459
     3 =',F10.4/)                                                       LECT-460
 1018 FORMAT (/' INTEGRATION STEP SIZE =',F8.5,10X,'MATCHING RADIUS =',FLECT-461
     18.3,' FERMI')                                                     LECT-462
 1019 FORMAT (/' SCATTERING ANGLES FROM',F7.3,' IN STEPS OF',F7.3,' UP TLECT-463
     1O',F8.3,' DEGREES   AVERAGED WITH +/-',F8.3,' DEGREES.'/)         LECT-464
 1020 FORMAT (10X,'****** AZ ******',1P,D15.8)                          LECT-465
 1021 FORMAT (/' SQUARE ROOT OF ELASTIC ENHANCEMENT FACTOR',8X,1P,D15.8/LECT-466
     1' SPIN CUT-OFF PARAMETER',27X,D15.8/' SQUARE ROOT OF LEVEL DENSITYLECT-467
     2 PARAMETER',11X,D15.8)                                            LECT-468
 1022 FORMAT (/' SQUARE ROOT OF ELASTIC ENHANCEMENT FACTOR WITHOUT FLUCTLECT-469
     1UATIONS',7X,1P,D15.8/' PARTICLE WIDTH FLUCTUATION DEGREE OF FREEDOLECT-470
     2M',25X,D15.8/' PARAMETERS OF MOLDAUER''S FORMULA',22X,3D15.8)     LECT-471
 1023 FORMAT (/I5,' FISSION DATA:')                                     LECT-472
 1024 FORMAT (10X,1P,2D15.8)                                            LECT-473
 1025 FORMAT (/' GAMMA DATA:')                                          LECT-474
 1026 FORMAT (5X,'TG0:',1P,D13.6,6X,'BN:',D13.6,4X,'FNUG:',D13.6,5X,'EGDLECT-475
     1:',D13.6,5X,'GGD:',D13.6/' DENSITY OF STATES GIVEN BY')           LECT-476
 1027 FORMAT (' READ VALUES'/(1X,I3,'  SA:',1P,D13.6,6X,'UX:',D13.6,5X,'LECT-477
     1TAU:',D13.6,6X,'SG:',D13.6/28X,'E0:',D13.6,6X,'EX:',D13.6,7X,'Z:',LECT-478
     20P,F5.0))                                                         LECT-479
 1028 FORMAT (1P,8D15.8)                                                LECT-480
 1029 FORMAT (' DENSITY OF STATES FOR CONTINUUM GIVEN BY')              LECT-481
 1030 FORMAT (' IQM =',I3,' OR IQMAX =',I3,' ARE LARGER THAN THE MAXIMUMLECT-482
     1 VALUES 35 OR 8 OF THE ASYMMETRIC ROTATIONAL MODEL'///' IN LECT  .LECT-483
     2..  STOP  ...')                                                   LECT-484
      END                                                               LECT-485
