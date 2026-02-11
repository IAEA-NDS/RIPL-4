C 07/03/07                                                      ECIS06  CALX-000
      SUBROUTINE CALX(NW,CW,DW,LO)                                      CALX-001
C CALX AND THE SUBROUTINES CALLED BY IT READ ALL THE INPUT EXCEPT FOR   CALX-002
C REDUCED NUCLEAR MATRIX ELEMENTS AND EXTERNAL FORM FACTORS.            CALX-003
C INPUT:     NW:      WORKING AREA FOR INTEGERS.                        CALX-004
C            CW:      WORKING AREA FOR CHARACTERS.                      CALX-005
C            DW:      WORKING AREA FOR DOUBLE PRECISION.                CALX-006
C                     NW,CW,DW ARE IN EQUIVALENCE BY CALL.              CALX-007
C            LO(I):   LOGICAL CONTROLS:                                 CALX-008
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).CALX-009
C               LO(2)  =.TRUE. SECOND ORDER VIBRATIONAL OR CONSTRAINED  CALX-010
C                              ASYMMETRIC ROTATIONAL MODEL.             CALX-011
C               LO(3)  =.TRUE. ANHARMONIC VIBRATIONAL OR ASYMMETRIC     CALX-012
C                              ROTATIONAL MODEL.                        CALX-013
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    CALX-014
C               LO(10) =.TRUE. DISPERSION RELATIONS FOR POTENTIALS.     CALX-015
C               LO(11) =.TRUE. DEFORMED COULOMB POTENTIAL.              CALX-016
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            CALX-017
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. CALX-018
C               LO(18) =.TRUE. PROJECTILE-TARGET ANTISYMMETRISATION.    CALX-019
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   CALX-020
C               LO(20) =.TRUE. DISPERSION RELATIONS FOR TRANSITION      CALX-021
C                              FORM-FACTORS.                            CALX-022
C               LO(21) =.TRUE. USUAL COUPLED EQUATIONS.                 CALX-023
C               LO(22) =.TRUE. NO USE OF PADE APPROXIMANTS.             CALX-024
C               LO(28) =.TRUE. COMPUTATION UP TO J-CONVERGENCE-(F: STOP CALX-025
C                              WHEN ALL COUPLING TERMS NEGLIGIBLE).     CALX-026
C               LO(29) =.TRUE. NO DIAGONAL TERMS IN SECOND MEMBER.      CALX-027
C               LO(30) =.TRUE. PURE DWBA CALCULATION.                   CALX-028
C               LO(31) =.TRUE. INPUT OF EXPERIMENTAL DATA AND CHI2      CALX-029
C                              CALCULATION.                             CALX-030
C               LO(131)=.TRUE. TOTAL SPIN IS TOO LARGE FOR COMPOUND     CALX-031
C                              NUCLEUS.                                 CALX-032
C               LO(36) =.TRUE. RESTART A SEARCH FROM TAPE 8.            CALX-033
C               LO(32) =.TRUE. AUTOMATIC SEARCH ON SOME PARAMETERS.     CALX-034
C               LO(41) =.TRUE. FACTORISATION OF 1/(1-COS(THETA)).       CALX-035
C               LO(42) =.TRUE. SCHMIDT'S ORTHOGONALISATION OF SOLUTIONS CALX-036
C                              IN USUAL COUPLED EQUATIONS.              CALX-037
C               LO(43) =.TRUE. INTERPOLATION ON TOTAL SPIN.             CALX-038
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     CALX-039
C               LO(45) =.TRUE. LIMITATION ON THE ANGULAR MOMENTA FOR    CALX-040
C                              COULOMB CORRECTIONS.                     CALX-041
C               LO(63) =.TRUE. PENETRABILITIES WRITTEN ON FILE 63.      CALX-042
C               LO(65) =.TRUE. PRINT COEFFICIENTS OF THE EXPANSION IN   CALX-043
C                              LEGENDRE POLYNOMIALS ON FILE 65.         CALX-044
C               LO(66) =.TRUE. NO CALCULATION AT EQUIDISTANT ANGLES.    CALX-045
C               LO(71) =.TRUE. NO DETAILED OUTPUT OF LOGICAL CONTROLS.  CALX-046
C               LO(74) =.TRUE. OUTPUT OF TIME IN DIFFERENT STEPS.       CALX-047
C               LO(75) =.TRUE. NO COMPLETE OUTPUT AT THE FIRST RUN OF A CALX-048
C                              SEARCH.                                  CALX-049
C               LO(76) =.TRUE. LO(51) TO LO(65) ARE ALWAYS USED.        CALX-050
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             CALX-051
C               LO(82) =.TRUE. OLD SIMPLIFIED COMPOUND NUCLEUS.         CALX-052
C               LO(84) =.TRUE. UNCOUPLED LEVELS FOR COMPOUND NUCLEUS.   CALX-053
C               LO(85) =.TRUE. FISSION TRANSMISSION COEFFICIENTS.       CALX-054
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      CALX-055
C               LO(99) =.TRUE. SCHROEDINGER EQUIVALENT TO DIRAC         CALX-056
C                              EQUATION.                                CALX-057
C               LO(100)=.TRUE. DIRAC EQUATION.                          CALX-058
C               LO(108)=.TRUE. DIAGONAL COULOMB CORRECTIONS ARE NEEDED. CALX-059
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    CALX-060
C               LO(115)=.TRUE. FIRST COMPUTATION FOR THIS ENERGY.       CALX-061
C               LO(116)=.TRUE. NO OUTPUT.                               CALX-062
C               LO(123)=.TRUE. IDENTICAL PARTICLES WITH SPIN.           CALX-063
C               LO(124)=.TRUE. COMPUTE TRANSMISSION COEFFICIENTS.       CALX-064
C               LO(128)=.TRUE. NO COPY OF UNCOUPLED FUNCTIONS AND       CALX-065
C                              PHASE-SHIFTS.                            CALX-066
C                                                                       CALX-067
C USE OF SOME PARTS OF THE WORKING AREA (W,NW,DW):                      CALX-068
C  1:   IPI(11,NCOLX) INTEGER VALUES FOR THE DESCRIPTION OF CHANNELS.   CALX-069
C           FIRST INDEX:                                                CALX-070
C       1 - PARITY (0 FOR + AND 1 FOR -).                               CALX-071
C       2 - MULTIPLICITY OF INCIDENT PARTICLE.                          CALX-072
C       3 - MULTIPLICITY OF THE TARGET.                                 CALX-073
C       4 - PRODUCT OF CHARGES.                                         CALX-074
C       5 - INDEX OF POTENTIAL.                                         CALX-075
C       6,7,8,9 - BEGINNING AND END IN THE TWO PARTS OF TABLE MF-FM     CALX-076
C                           (SEE DEPH).                                 CALX-077
C      10 - MAXIMUM ANGULAR MOMENTUM.                                   CALX-078
C      11 - INDEX OF POTENTIAL ENERGY DEPENDENT BY DISPERSION RELATIONS.CALX-079
C  NIPP:    IPP(2,18,NPP)/PIP(18,NPP) DISPERSION PARAMETERS.            CALX-080
C     1,1 - FIRST LEVEL USING THIS POTENTIAL BEFORE INPUT OF PARAMETERS.CALX-081
C           REPLACED BY INDICATION TO USE LABORATORY ENERGY (1), CENTRE CALX-082
C           OF MASS ENERGY (-1) OR READ COEFFICIENTS (0).               CALX-083
C     2,1 - N2 - POWER FOR LARGE NEGATIVE ENERGY CORRECTIONS.           CALX-084
C     1,2 - NV - |NV| POWER FOR VOLUME POTENTIAL.                       CALX-085
C     2,2 - NS - |NS| POWER FOR SURFACE POTENTIAL.                      CALX-086
C             SUM OF TWO TERMS IF NV OR NS ARE NEGATIVE,                CALX-087
C             ABSENCE OF VOLUME OR SURFACE TERM IF NV OR NS ARE 0.      CALX-088
C     1,3 - NL - |NL| POWER FOR SPIN-ORBIT POTENTIAL.                   CALX-089
C       4 - ENERGY CORRESPONDING TO THE DEPTHS READ.                    CALX-090
C       5 - FERMI ENERGY EF.                                            CALX-091
C       6 - THRESHOLD ENERGY EP.                                        CALX-092
C       7 - LARGE ENERGY STARTING VALUE ABOVE FERMI ENERGY EA.          CALX-093
C       8 - EXPONENTIAL VARIATION OF REAL SPIN-ORBIT.                   CALX-094
C       9 - LINEAR VARIATION OF IMAGINARY SPIN-ORBIT.                   CALX-095
C      10 - BV PARAMETER FOR VOLUME POTENTIALS.                         CALX-096
C      11 - STRENGTH OF LARGE POSITIVE ENERGY TERM IN VOLUME POTENTIAL, CALX-097
C           SECOND BV' PARAMETER FOR VOLUME POTENTIALS IF NV<0.         CALX-098
C      12 - EXPONENTIAL DECREASE IN SQRT|E| FOR LARGE ENERGY TERMS OF   CALX-099
C           VOLUME POTENTIALS OR FRACTION IN THE FIRST TERM IF NV<0.    CALX-100
C      13 - BS PARAMETER FOR SURFACE POTENTIALS.                        CALX-101
C      14 - EXPONENTIAL DECREASE OF A SURFACE POTENTIAL OR SECOND BS'   CALX-102
C           PARAMETER FOR SURFACE POTENTIALS IF NS<0.                   CALX-103
C      15 - NON-LOCALITY RANGE PARAMETER OF SURFACE POTENTIAL OR        CALX-104
C           FRACTION IN THE FIRST TERM IF NS<0.                         CALX-105
C      16 - BL PARAMETER FOR SPIN-ORBIT POTENTIALS.                     CALX-106
C      17 - EXPONENTIAL VARIATION OF H.F. REAL VOLUME POTENTIAL.        CALX-107
C  NWV:     WV(20,NCOLX) IBM-DOUBLE PRECISION VALUES FOR THE CHANNELS.  CALX-108
C       1 - MASS OF THE PARTICLE.                                       CALX-109
C       2 - MASS OF THE TARGET.                                         CALX-110
C       3 - ENERGY IN THE CENTRE OF MASS IN MEV.                        CALX-111
C       4 - K WAVE NUMBER.                                              CALX-112
C       5 - COULOMB PARAMETER.                                          CALX-113
C       6 - REDUCED MASS.                                               CALX-114
C       7 - REDUCED ENERGY.                                             CALX-115
C       8 - STEP SIZE FOR THIS LEVEL.                                   CALX-116
C       9 - SQUARE ROOT OF COEFFICIENT OF SCALAR POTENTIALS.            CALX-117
C      10 - SQUARE ROOT OF COEFFICIENT OF COULOMB POTENTIALS.           CALX-118
C      11 - K WAVE NUMBER MULTIPLIED BY RATIO OF STEP SIZES.            CALX-119
C      12 - REDUCED ENERGY TERM.                                        CALX-120
C      13 - ENERGY IN THE LABORATORY SYSTEM IN MEV.                     CALX-121
C      14,15,16,17,18,19,20 - DISPERSION CORRECTIONS (SEE DISP).        CALX-122
C                                                                       CALX-123
C THE COMMON /CONVE/ IS USED IN CALC, CALX, LECT, CAL1, STDP, FOLD,     CALX-124
C                               INTI, INSH, INSI, INRH, INRI AND REST.  CALX-125
C THE COMMON /DCHI2/ IS USED IN CALC, CALX, RESU, VARI, EVAL, REST      CALX-126
C                               AND FITE.                               CALX-127
C THE COMMON /NCOMP/ IS USED IN CALC, CALX, LECT, CONU, GGDR, CAL1,     CALX-128
C                               QUAN, SCAM, SCHE, RESU, VARI, EVAL,     CALX-129
C                               AND REST.                               CALX-130
C                                                                       CALX-131
C FOR THE COMMON  /ANGUL/ SEE LECT.                                     CALX-132
C FOR THE COMMON  /ADDRE/, /COUPL/, /INTEG/ AND /TITRE/ SEE CALC.       CALX-133
C                                                                       CALX-134
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ADDRE/:                     CALX-135
C  NIPH:      NUMBER OF PHONONS FOR THE HARMONIC VIBRATIONAL MODEL.     CALX-136
C  NJIT:      DATA FOR INTERPOLATION ON TOTAL SPIN.                     CALX-137
C  NWV:       NON INTEGER VALUES FOR THE CHANNELS.                      CALX-138
C  NIPP:      DISPERSION PARAMETERS.                                    CALX-139
C  NPAR:      INDICATIONS FOR NUCLEAR PARAMETERS.                       CALX-140
C  NPAA:      VALUES OF NUCLEAR PARAMETERS.                             CALX-141
C  NSCN:      LEVEL DENSITY DESCRIPTION.                                CALX-142
C  NFIS:      FISSION DATA FOR COMPOUND NUCLEUS.                        CALX-143
C  NGAM:      GAMMA DATA FOR COMPOUND NUCLEUS.                          CALX-144
C  NNIV:      ADDRESSES OF REDUCED NUCLEAR MATRIX ELEMENTS AND COULOMB  CALX-145
C             INTEGRALS IN NIV(NCOLL,NCOLL,3).                          CALX-146
C  NPOT:      OPTICAL POTENTIAL PARAMETERS.                             CALX-147
C  NBETA:     DEFORMATION PARAMETERS.                                   CALX-148
C  NFM:       HELICITIES AND OBSERVABLES (SEE DEPH,LECD AND OBSE).      CALX-149
C  NTGX:      BEGINNING OF CHI2 AND NORMALISATION OF DATA.              CALX-150
C  NDONN:     EXPERIMENTAL DATA.                                        CALX-151
C  NNVI:      TABLE OF ADDRESSES OF COUPLING COEFFICIENTS.              CALX-152
C  NDE:       SEARCH ACCURACIES.                                        CALX-153
C  NISE:      INDEXES OF THE VARIABLE PARAMETERS IN SEARCH.             CALX-154
C  NRC:       PERMANENT WORKING FIELD FOR THE SEARCH.                   CALX-155
C  NIW:       INTEGER WORKING FIELD FOR THE SEARCH.                     CALX-156
C  NNWI:      SAME AS NNVI FOR SYMMETRISED EQUATIONS.                   CALX-157
C  NCC:       TABLE OF ENERGIES, L*(L+1) AND L*S  (SEE QUAN).           CALX-158
C  MCC:       SAME AS NCC FOR SYMMETRISED EQUATIONS.                    CALX-159
C  NXA:       TABLE OF COEFFICIENTS OF SYMMETRISATION.                  CALX-160
C  NAM1:      COMPUTATION OF OBSERVABLES (SEE OBSE).                    CALX-161
C  NFAC:      TABLE OF LOG OF FACTORIALS FOR GEOMETRIC COEFFICIENTS.    CALX-162
C  NFG:       COULOMB FUNCTIONS AND FINITE INTEGRALS.                   CALX-163
C  NXG:       COULOMB PHASES AND INFINITE INTEGRALS.                    CALX-164
C  NRES:      FUNCTIONS FOR THE SEARCH.                                 CALX-165
C  NXX:       VARIABLES FOR THE SEARCH.                                 CALX-166
C  NT:        TABLE OF REDUCED NUCLEAR MATRIX ELEMENTS.                 CALX-167
C   DEFINED:  NIPH,NJIT,NWV,NIPP,NPAR,NPAA,NSCN,NFIS,NGAM,NNIV,NPOT,    CALX-168
C             NBETA,NFM,NTGX,NDONN,NNVI,NDE,NISE,NRC,NIW,NNWI,NCC,MCC,  CALX-169
C             NXA,NAM1,NFAC,NFG,NXG,NRES,NXX,NT.                        CALX-170
C   USED:     NIPH,NJIT,NWV,NIPP,NPAR,NPAA,NSCN,NFIS,NGAM,NNIV,NPOT,    CALX-171
C             NBETA,NFM,NTGX,NDONN,NNVI,NDE,NISE,NRC,NIW,NNWI,NCC,MCC,  CALX-172
C             NXA,NAM1,NFAC,NFG,NXG,NRES,NXX,NT.                        CALX-173
C                                                                       CALX-174
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ANGUL/:                     CALX-175
C  THETA1:    FIRST ANGLE FOR COMPUTATION AT EQUIDISTANT ANGLES.        CALX-176
C  DTHETA:    STEP FOR COMPUTATION AT EQUIDISTANT ANGLES.               CALX-177
C  THETA2:    LAST ANGLE FOR COMPUTATION AT EQUIDISTANT ANGLES.         CALX-178
C  NCJ:       NUMBER OF FACTORISATIONS OF 1/(1-COS(THETA)) IN AMPLITUDE.CALX-179
C  NL(1):     POWER OF (1-COS(THETA)) FOR THE EXPANSION IN LEGENDRE     CALX-180
C             POLYNOMIALS OF THE INTERFERENCE BETWEEN COULOMB AND       CALX-181
C             NUCLEAR ELASTIC SCATTERING. POWER OF (1-COS(THETA)**2)    CALX-182
C             IF LO(18) IS .TRUE..                                      CALX-183
C  NL(2):     NUMBER OF LEGENDRE POLYNOMIALS FOR ELASTIC SCATTERING OF  CALX-184
C             CHARGED PARTICLES.                                        CALX-185
C  NL(3):     NUMBER OF LEGENDRE POLYNOMIALS FOR ELASTIC SCATTERING OF  CALX-186
C             UNCHARGED PARTICLES, INELASTIC SCATTERING AND COMPOUND    CALX-187
C             NUCLEUS.                                                  CALX-188
C   DEFINED:  NCJ,NL.                                                   CALX-189
C   USED:     THETA1,DTHETA,THETA2,NCJ,NL.                              CALX-190
C                                                                       CALX-191
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     CALX-192
C  H:         STEP SIZE FOR INTEGRATION.                                CALX-193
C  BJM:       CONVERGENCE COEFFICIENT OF IMAGINARY POTENTIAL.           CALX-194
C  EITER:     CONVERGENCE CRITERION FOR S-MATRIX.                       CALX-195
C  ACONV:     CONVERGENCE CRITERION FOR POTENTIAL AND FUNCTION.         CALX-196
C  CONJ:      CONVERGENCE CRITERION FOR THE SCATTERING AMPLITUDES.      CALX-197
C  HCONV:     ACONV*H*H CONVERGENCE CRITERION FOR SECOND MEMBERS.       CALX-198
C   DEFINED:  H,BJM,EITER,ACONV,CONJ,HCONV.                             CALX-199
C   USED:     BJM,EITER,ACONV,CONJ.                                     CALX-200
C                                                                       CALX-201
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     CALX-202
C  IQM:       MAXIMUM L-VALUE OF DEFORMATION IN ROTATIONAL MODEL.       CALX-203
C  IQMAX:     MAXIMUM L-EXPANSION IN ROTATIONAL MODEL.                  CALX-204
C  NBT1:      NUMBER OF PHONONS (VIBRATIONS).                           CALX-205
C  NFA:       NUMBER OF LOGARITHMS OF FACTORIALS.                       CALX-206
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             CALX-207
C  NSPIN:     TWICE THE K-VALUE OF THE ROTATIONAL BAND.                 CALX-208
C  NVA:       NUMBER OF NUCLEAR PARAMETERS.                             CALX-209
C   DEFINED:  IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA.                         CALX-210
C   USED:     NFA,NPP,NVA.                                              CALX-211
C                                                                       CALX-212
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCHI2/:                     CALX-213
C  CHI2:      CHI-SQUARE COMPUTED IN SUBROUTINE RESU.                   CALX-214
C  CHI2M:     MINIMUM CHI-SQUARE IN THE SEARCH.                         CALX-215
C  YY(1):     STEP SIZE IN THE SEARCH.                                  CALX-216
C  YY(2):     HALF OF THE SUCCESS MULTIPLICATIVE FACTOR OF THE STEP.    CALX-217
C  YY(3):     VARIOUS MEANINGS.  SEE FITE.                              CALX-218
C   DEFINED:  YY.                                                       CALX-219
C   NOT USED: CHI2,CHI2M.                                               CALX-220
C                                                                       CALX-221
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /INTEG/:                     CALX-222
C  IDMT:      TOTAL WORKING FIELD LENGTH AS SINGLE PRECISION.           CALX-223
C  NPLACE:    MAXIMUM WORKING FIELD USED.                               CALX-224
C  NCOLL:     NUMBER OF COUPLED CHANNELS.                               CALX-225
C  NJMAX:     MAXIMUM NUMBER OF J-VALUES.                               CALX-226
C  ITERM:     MAXIMUM NUMBER OF ITERATIONS.                             CALX-227
C  JDM:       MINIMUM NUMBER OF TOTAL SPIN.                             CALX-228
C  JIT:       NUMBER OF RATES OF INTERPOLATION ON TOTAL SPIN.           CALX-229
C  KMIN:      MINIMUM NUMBER OF J VALUES FOR COMPOUND NUCLEUS.          CALX-230
C  KMAX:      MAXIMUM NUMBER OF J VALUES FOR COMPOUND NUCLEUS.          CALX-231
C  NCOLS:     NUMBER OF CHANNELS WITH ANGULAR DISTRIBUTIONS.            CALX-232
C  NCOLT:     NUMBER OF CHANNELS INCLUDING UNCOUPLED STATES.            CALX-233
C  NBET:      NUMBER OF DIFFERENT DEFORMATIONS (VIBRATIONS+ROTATIONS).  CALX-234
C  LMX:       DIFFERENCE BETWEEN NUMBER OF J AND OF COULOMB FUNCTIONS.  CALX-235
C  LMAX1:     NUMBER OF L-VALUES FOR COULOMB FUNCTIONS.                 CALX-236
C  NLT:       MEMORIES NEEDED FOR LEGENDRE POLYNOMIALS.                 CALX-237
C  ISM:       NUMBER OF INTEGRATION STEPS.                              CALX-238
C  NJC:       MAXIMUM NUMBER OF OBSERVABLES AT EQUIDISTANT ANGLES.      CALX-239
C  JTX:       MAXIMUM NUMBER OF CALCULATED VALUES FOR A PLOT.           CALX-240
C  KCC:       NUMBER OF INDEPENDENT AMPLITUDES WITH UNCOUPLED STATES.   CALX-241
C  MS1:       LARGEST PARTICLE MULTIPLICITY.                            CALX-242
C  MS2:       LARGEST TARGET MULTIPLICITY.                              CALX-243
C  KBA:       NUMBER OF INDEPENDENT AMPLITUDES WITHOUT UNCOUPLED STATES.CALX-244
C  KAB:       MAXIMUM NUMBER OF EQUATIONS.                              CALX-245
C  KBC:       MAXIMUM NUMBER OF SOLUTIONS.                              CALX-246
C  JTH:       MAXIMUM NUMBER OF ANGLES FOR A PLOT.                      CALX-247
C  NCOLR:     NUMBER OF EXPERIMENTAL ANGULAR DISTRIBUTIONS.             CALX-248
C  NREC:      NUMBER OF VARIABLES IN SEARCH.                            CALX-249
C  NTOT:      NUMBER OF EXPERIMENTAL DATA.                              CALX-250
C  LMAX2:     NUMBER OF L VALUES FOR COULOMB PHASES.                    CALX-251
C  KE:        CONTROL OF SEARCH (SEE FITE).                             CALX-252
C  ITEMM:     MAXIMUM NUMBER OF ITERATIONS AT THE BEGINNING OF A RUN.   CALX-253
C  MCM(1):    MAXIMUM ANGULAR MOMENTUM FOR CENTRAL COULOMB CORRECTIONS. CALX-254
C  MCM(2):    MAXIMUM ANGULAR MOMENTUM FOR SPIN-ORBIT COULOMB CORRECT.  CALX-255
C  NCT(1):    NUMBER OF EQUATIONS FOR POSITIVE PARITY.                  CALX-256
C  NCT(2):    NUMBER OF EQUATIONS FOR NEGATIVE PARITY.                  CALX-257
C  NCT(3):    NUMBER OF SOLUTIONS FOR POSITIVE PARITY.                  CALX-258
C  NCT(4):    NUMBER OF SOLUTIONS FOR NEGATIVE PARITY.                  CALX-259
C  NCT(5):    NUMBER OF COMPOUND NUCLEUS COEFFICIENTS, POSITIVE PARITY. CALX-260
C  NCT(6):    NUMBER OF COMPOUND NUCLEUS COEFFICIENTS, NEGATIVE PARITY. CALX-261
C   DEFINED:  NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCOLS,NCOLT,   CALX-262
C             NBET,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,   CALX-263
C             JTH,NCOLR,NREC,NTOT,LMAX2,KE,ITEMM,MCM,NCT.               CALX-264
C   USED:     IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMAX,NCOLS,NCOLT,   CALX-265
C             NBET,LMX,LMAX1,JTX,KCC,KAB,KBC,JTH,NCOLR,NREC,NTOT,LMAX2. CALX-266
C                                                                       CALX-267
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     CALX-268
C  NSP(1):    NUMBER OF UNCOUPLED LEVELS FOR COMPOUND NUCLEUS           CALX-269
C             INCLUDING DISCRETISATION OF CONTINUUM.                    CALX-270
C  NSP(2):    NUMBER OF THESE LEVELS WITH ANGULAR DISTRIBUTION.         CALX-271
C  NSP(3):    NUMBER OF THESE LEVELS WITHOUT ANGULAR DISTRIBUTION.      CALX-272
C  NFISS:     NUMBER OF FISSION TRANSMISSION COEFFICIENTS.              CALX-273
C  NRD:       NUMBER OF GAMMA TRANSMISSION COEFFICIENTS.                CALX-274
C  NCONT:     NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS.                 CALX-275
C  NCOJ:      NUMBER OF SPINS OF THE TARGET FOR A CONTINUUM.            CALX-276
C  NCONS:     NUMBER OF LEVEL DENSITIES NEEDED.                         CALX-277
C  NIE:       NUMBER OF UNCOUPLED STATES ADDED FOR DISCRETISATION.      CALX-278
C  NCOLX:     TOTAL NUMBER OF LEVELS WITHOUT DISCRETISATION.            CALX-279
C  NDP:       ADDRESS OF FISSION AND GAMMA FINAL RESULTS.               CALX-280
C  NDQ:       ADDRESS OF FISSION AND GAMMA INTERMEDIATE RESULTS.        CALX-281
C  ACN1:      RATIO SIZE/STEP FOR DISCRETISATION OF A CONTINUUM.        CALX-282
C  ACN2:      MAXIMUM NUMBER OF STEPS BY MEV FOR A CONTINUUM.           CALX-283
C  AZ(6):     DEFORMED SPIN-ORBIT PARAMETERS. SEE ALSO COMMENT IN       CALX-284
C             INPUT DESCRIPTION AND SUBROUTINE QUAN.                    CALX-285
C  BZ(5):     HAUSER-FESHBACH AND MOLDAUER'S PARAMETERS DESCRIBED BELOW.CALX-286
C   BZ(1):    SQUARE ROOT OF ELASTIC ENHANCEMENT.                       CALX-287
C   BZ(2):    IF LO(82)=.TRUE., SPIN CUT-OFF PARAMETER.                 CALX-288
C             IF LO(82)=.FALSE., PARTICLE DEGREES OF FREEDOM.           CALX-289
C   BZ(3):    SQUARE ROOT OF LEVEL DENSITY PARAMETER. IF LO(82)=LO(87)= CALX-290
C             .FALSE., PARAMETER BZ(3) IN MOLDAUER'S FORMULA OF INPUT   CALX-291
C             DESCRIPTION.                                              CALX-292
C   BZ(4):    IF LO(82)=LO(87)=.FALSE., PARAMETER BZ(4) IN SAME FORMULA.CALX-293
C   BZ(5):    IF LO(82)=LO(87)=.FALSE., PARAMETER BZ(5) IN SAME FORMULA.CALX-294
C  TG0:       SLOW S-WAVE NEUTRON GAMMA WIDTH/SPACING FOR NORMALISATION.CALX-295
C  BN:        NEUTRON SEPARATION ENERGY.                                CALX-296
C  FNUG:      RADIATIVE DEGREES OF FREEDOM.                             CALX-297
C  EGD:       ENERGY OF THE GIANT DIPOLE RESONANCE.                     CALX-298
C  GGD:       RESONANCE WIDTH.                                          CALX-299
C  TG1:       DERIVED DATA FOR GAMMA IN COMPOUND NUCLEUS.               CALX-300
C  SGSQ:      DERIVED DATA FOR GAMMA IN COMPOUND NUCLEUS.               CALX-301
C   DEFINED:  NSP,NFISS,NRD,NCONT,NCOJ,NCONS,NCOLX,ACN1,ACN2.           CALX-302
C   USED:     NSP,NFISS,NRD,NCONT,NCOJ,NCONS,NCOLX,ACN1,ACN2.           CALX-303
C   NOT USED: NIE,NDP,NDQ,AZ,BZ,TG0,BN,FNUG,EGD,TGD,TG1,SGSQ.           CALX-304
C             AZ(6),BZ(5) AND FOLLOWING ARE SOMETIMES GROUPED IN AZ(18);CALX-305
C             THIS IS NECESSARY FOR AUTOMATIC SEARCH.                   CALX-306
C                                                                       CALX-307
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /TITRE/:                     CALX-308
C  TITLE(18): TITLE OF THE RUN PRINTED AS HEADING OF RESULTS.           CALX-309
C   DEFINED:  TITLE.                                                    CALX-310
C   USED:     TITLE.                                                    CALX-311
C                                                                       CALX-312
C***********************************************************************CALX-313
      IMPLICIT REAL*8 (A-H,O-Z)                                         CALX-314
      LOGICAL LO(150)                                                   CALX-315
      DIMENSION NW(2,*),NGR(2),NPR(2),DESC(3),DW(*)                     CALX-316
      CHARACTER*4 FIN,DESC,TITLE,CW(2,1)                                CALX-317
      COMMON /ADDRE/ NIPH,NJIT,NWV,NIPP,NPAR,NPAA,NSCN,NFIS,NGAM,NNIV,NPCALX-318
     1OT,NBETA,NFM,NTGX,NDONN,NNVI,NDE,NISE,NRC,NIW,NNWI,NCC,MCC,NXA,NAMCALX-319
     21,NFAC,NFG,NXG,NRES,NXX,NT,NIVQ,NIVY,NIVZ,NCOI,MIPI,NXD,MWV,NIXT,NCALX-320
     3TY,NSM,NTX,NRY,NRCO,NRDO,NVC1,NVC2,NNC,NCX                        CALX-321
      COMMON /ANGUL/ THETA1,THETA2,DTHETA,DTHE,NCJ,NL(3),JMM(2)         CALX-322
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       CALX-323
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   CALX-324
      COMMON /DCHI2/ CHI2,CHI2M,YY(3)                                   CALX-325
      COMMON /INOUT/ MR,MW,MS                                           CALX-326
      COMMON /INTEG/ IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCOCALX-327
     1LS,NCOLT,NBET,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,JTCALX-328
     2H,NCOLR,NREC,NTOT,LMAX2,KE,ITEMM,KXT,LMAX3,NRZ,NTZ,IPM,IPK,MCM(2),CALX-329
     3NCT(6)                                                            CALX-330
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQCALX-331
     1,ACN1,ACN2,AZ(6),BZ(5),TG0,BN,FNUG,EGD,GGD,TG1,SGSQ               CALX-332
      COMMON /TITRE/ TITLE(18)                                          CALX-333
      DATA FIN,DESC /'FIN ','DESC','RIPT','ION '/                       CALX-334
    1 READ (MR,1000) TITLE                                              CALX-335
      DO 2 J=1,3                                                        CALX-336
      IF (TITLE(J).NE.DESC(J)) GO TO 3                                  CALX-337
    2 CONTINUE                                                          CALX-338
      CALL INPA                                                         CALX-339
      CALL INPB                                                         CALX-340
      CALL INPC                                                         CALX-341
      GO TO 1                                                           CALX-342
    3 IF (TITLE(1).EQ.FIN) RETURN                                       CALX-343
      DO 4 I=1,100                                                      CALX-344
    4 LO(I)=.FALSE.                                                     CALX-345
      READ (MR,1001) (LO(I),I=1,50)                                     CALX-346
C IF LO(36)=.TRUE. RETURN TO RESTART A SEARCH SAVED ON TAPE 8.          CALX-347
      IF (LO(36)) RETURN                                                CALX-348
      READ (MR,1001) (LO(I),I=51,100)                                   CALX-349
      READ (MR,1002) NCOLL,NJMAX,ITERM,NPP,NCJ,NGR,NPR,LMZ,JDM,LML,JIT,MCALX-350
     1N                                                                 CALX-351
C ELIMINATION OF CONTRADICTIONS BETWEEN LOGICAL CONTROLS.               CALX-352
      IF (NCOLL.LE.1) LO(28)=.TRUE.                                     CALX-353
      IF (.NOT.LO(31)) LO(32)=.FALSE.                                   CALX-354
      IF (LO(14)) LO(13)=.TRUE.                                         CALX-355
      IF (LO(19)) LO(11)=.TRUE.                                         CALX-356
      IF (LO(30)) ITERM=1                                               CALX-357
      IF (LO(18)) LO(29)=.TRUE.                                         CALX-358
      IF (LO(100)) LO(29)=.FALSE.                                       CALX-359
      IF (LO(100).AND.LO(95)) LO(94)=.FALSE.                            CALX-360
      LO(109)=LO(99).OR.LO(100)                                         CALX-361
      IF (LO(109)) LO(8)=.TRUE.                                         CALX-362
      IF (LO(82)) LO(84)=.FALSE.                                        CALX-363
      IF (LO(82)) LO(85)=.FALSE.                                        CALX-364
      IF (LO(82)) LO(86)=.FALSE.                                        CALX-365
      IF (LO(82).OR.LO(84).OR.LO(85).OR.LO(86)) LO(81)=.TRUE.           CALX-366
      LO(44)=LO(44).AND.LO(11)                                          CALX-367
      IF (.NOT.LO(99)) GO TO 5                                          CALX-368
      LO(11)=.FALSE.                                                    CALX-369
      LO(19)=.FALSE.                                                    CALX-370
      IF (LO(1)) LO(3)=.FALSE.                                          CALX-371
      IF (.NOT.LO(1)) LO(2)=.FALSE.                                     CALX-372
      IF ((.NOT.LO(10)).OR.(.NOT.LO(12))) LO(20)=.FALSE.                CALX-373
    5 LO(124)=LO(81).OR.LO(63)                                          CALX-374
      WRITE (MW,1003) IDMT                                              CALX-375
C OUTPUT OF LOGICAL CONTROLS.                                           CALX-376
      IF (LO(71)) GO TO 6                                               CALX-377
      WRITE (MW,1004) (LO(I),LO(I+53),I=1,3),(LO(I+53),LO(I),I=4,5),LO(6CALX-378
     1),LO(59),(LO(I),I=7,8),(LO(I),LO(51+I),I=9,10)                    CALX-379
      WRITE (MW,1005) (LO(I),LO(I+51),I=11,13),LO(14),LO(15),LO(65),LO(1CALX-380
     16),(LO(I),LO(I+49),I=17,20),(LO(I),LO(I+50),I=21,22)              CALX-381
      WRITE (MW,1006) (LO(I),LO(I+50),I=23,25),LO(76),LO(26),LO(27),LO(7CALX-382
     17),LO(78),LO(28),LO(29),LO(81),LO(30),LO(82),LO(83),LO(31)        CALX-383
      WRITE (MW,1007) (LO(I),LO(I+52),I=32,34),LO(87),LO(35),LO(36),LO(9CALX-384
     11),(LO(I+51),LO(I),I=41,44),LO(45),LO(46),LO(96)                  CALX-385
      WRITE (MW,1008) LO(47),LO(97),(LO(I),LO(I+47),I=51,53)            CALX-386
      GO TO 7                                                           CALX-387
    6 WRITE (MW,1009) (LO(I),I=1,100)                                   CALX-388
    7 WRITE (MW,1010) TITLE                                             CALX-389
      IF (LO(109).AND.(.NOT.LO(8))) WRITE (MW,1011)                     CALX-390
      READ (MR,1012) H,RM,BJM,EITER,ACONV,CONJ                          CALX-391
C DEFECT VALUES OF NJMAX,ITERM,NPP,NCJ,NPR,NGR,EITER,ACONV,CONJ.        CALX-392
C FOR NPR AND NGR, SEE DEPH.                                            CALX-393
      IF (NJMAX.EQ.0) NJMAX=20                                          CALX-394
      IF (ITERM.EQ.0) ITERM=20                                          CALX-395
      IF (NPP.EQ.0) NPP=1                                               CALX-396
      IF (NPR(1).EQ.0) NPR(1)=1                                         CALX-397
      IF (NPR(2).EQ.0) NPR(2)=1                                         CALX-398
      IF (NGR(1).EQ.0) NGR(1)=2                                         CALX-399
      IF (NGR(2).EQ.0) NGR(2)=2                                         CALX-400
      IF (NCJ.LE.0) NCJ=1                                               CALX-401
      IF (JIT.EQ.0) JIT=1                                               CALX-402
      IF (.NOT.LO(43)) JIT=0                                            CALX-403
      IF (MN.LE.0) MN=1                                                 CALX-404
      NJMAX=MN*NJMAX                                                    CALX-405
      IF (EITER.EQ.0.D0) EITER=1.D-5                                    CALX-406
      IF (ACONV.EQ.0.D0) ACONV=1.D-5                                    CALX-407
      IF (CONJ.EQ.0.D0) CONJ=1.D-5                                      CALX-408
C OUTPUT OF TITLE, MASSES, ETC...                                       CALX-409
      WRITE (MW,1013) NJMAX,JDM,CONJ,NGR,NPR                            CALX-410
      IF (LO(41)) WRITE (MW,1014) NCJ                                   CALX-411
      IF (LML.NE.0) WRITE (MW,1015) LML                                 CALX-412
      IF (.NOT.LO(21)) WRITE (MW,1016) ITERM,EITER,ACONV                CALX-413
      LO(128)=LO(21).OR.LO(29).OR.BJM.NE.0.D0                           CALX-414
      IF (BJM.NE.0.D0) WRITE (MW,1017) BJM                              CALX-415
      IF (LO(21)) WRITE (MW,1018)                                       CALX-416
      IF (LO(21).AND.LO(42)) WRITE (MW,1019) ITERM                      CALX-417
C INITIALISATION OF A SEARCH OR A SINGLE RUN.                           CALX-418
      DO 8 I=111,120                                                    CALX-419
    8 LO(I)=.FALSE.                                                     CALX-420
      LO(115)=.TRUE.                                                    CALX-421
C LIMITATION ON ANGULAR MOMENTUM FOR COULOMB CORRECTIONS.               CALX-422
      MCM(1)=3                                                          CALX-423
      MCM(2)=2                                                          CALX-424
      IF (.NOT.LO(45)) GO TO 9                                          CALX-425
      READ (MR,1002) MC1,MC2                                            CALX-426
      IF (MC1.GT.0) MCM(1)=MIN0(MC1,5)                                  CALX-427
      IF (MC2.GT.0) MCM(2)=MIN0(MC2,4)                                  CALX-428
      IF (MC1.LT.0) MCM(1)=0                                            CALX-429
      IF (MC2.LT.0) MCM(2)=0                                            CALX-430
      WRITE (MW,1021) MCM                                               CALX-431
    9 NL2=3*NJMAX                                                       CALX-432
      NL3=2*NJMAX                                                       CALX-433
C LEGENDRE POLYNOMIALS DATA.                                            CALX-434
      IF (.NOT.LO(65)) GO TO 10                                         CALX-435
      READ (MR,1002) NL                                                 CALX-436
      IF (NL(1).EQ.0) NL(1)=2                                           CALX-437
      IF (NL(2).NE.0) NL2=NL(2)                                         CALX-438
      IF (NL(3).NE.0) NL3=NL(3)                                         CALX-439
      IF (LO(65)) WRITE (MW,1022) NL                                    CALX-440
C COMPOUND NUCLEUS DATA.                                                CALX-441
   10 NSP(1)=0                                                          CALX-442
      NSP(2)=0                                                          CALX-443
      NSP(3)=0                                                          CALX-444
      NFISS=0                                                           CALX-445
      NRD=0                                                             CALX-446
      NCONT=0                                                           CALX-447
      NCONS=0                                                           CALX-448
      IF (.NOT.(LO(84).OR.LO(85).OR.LO(86))) GO TO 11                   CALX-449
      READ (MR,1023) NSP(1),NSP(2),NFISS,NRD,NCONT,NCOJ,KMIN,KMAX,ACN1,ACALX-450
     1CN2                                                               CALX-451
      IF (.NOT.LO(84)) NSP(1)=0                                         CALX-452
      NSP(2)=MIN0(NSP(1),NSP(2))                                        CALX-453
      NSP(3)=NSP(1)-NSP(2)                                              CALX-454
      IF (NSP(3).LT.NCONT) GO TO 20                                     CALX-455
      IF (.NOT.LO(85)) NFISS=0                                          CALX-456
      IF (.NOT.LO(86)) NRD=0                                            CALX-457
      IF (NFISS.EQ.0) LO(85)=.FALSE.                                    CALX-458
      IF (LO(84).OR.LO(85).OR.LO(86)) WRITE (MW,1024) NSP,NFISS,NRD,NCONCALX-459
     1T                                                                 CALX-460
      NCONS=NCONT                                                       CALX-461
      IF (LO(86).AND.NRD.EQ.0) NCONS=NCONS+1                            CALX-462
      IF (NCONT.EQ.0) GO TO 11                                          CALX-463
      IF (NCOJ.LE.0) NCOJ=30                                            CALX-464
      IF (ACN1.LE.1.D0) ACN1=8.D0                                       CALX-465
      IF (ACN2.LE.1.D0) ACN2=8.D0                                       CALX-466
      WRITE (MW,1025) NCOJ,ACN1,ACN2                                    CALX-467
   11 NCOLX=NCOLL+NSP(1)                                                CALX-468
      NCOLS=NCOLL+NSP(2)                                                CALX-469
      NCOLT=NCOLX                                                       CALX-470
      NIPH=(11*NCOLX+1)/2+1                                             CALX-471
      NJIT=NIPH+NCOLL                                                   CALX-472
      IF (LO(7)) NJIT=NIPH                                              CALX-473
      NWV=NJIT+JIT                                                      CALX-474
      NIPP=NWV+22*NCOLX                                                 CALX-475
      NPAR=NIPP+17*NPP                                                  CALX-476
      NPLACE=NPAR                                                       CALX-477
      CALL MEMO('CALX',IDMT,NPLACE)                                     CALX-478
C INTERPOLATION DATA.                                                   CALX-479
      IF (.NOT.LO(43)) GO TO 13                                         CALX-480
      K=NJIT                                                            CALX-481
      READ (MR,1002) (NW(1,NJIT+I-1),NW(2,NJIT+I-1),I=1,JIT)            CALX-482
      WRITE (MW,1026) (NW(1,NJIT+I-1),NW(2,NJIT+I-1),I=1,JIT)           CALX-483
      M=0                                                               CALX-484
      L=-1                                                              CALX-485
      DO 12 I=1,JIT                                                     CALX-486
      IF (NW(1,NJIT+I-1).LT.L) GO TO 21                                 CALX-487
      L=NW(1,NJIT+I-1)                                                  CALX-488
      NW(2,NJIT+I-1)=NW(2,NJIT+I-1)-M                                   CALX-489
      M=M+NW(2,NJIT+I-1)                                                CALX-490
   12 K=K+2                                                             CALX-491
      WRITE (MW,1027) (NW(2,NJIT+I-1),I=1,JIT)                          CALX-492
C INPUT OF LEVEL DESCRIPTIONS.                                          CALX-493
   13 CALL LECL(NCOLX,NCOLL,NCONT,IDMT-NPAR,LO,NW,NW(1,NIPH),DW(NWV),NW(CALX-494
     11,NIPP),NW(1,NPAR),DW(NPAR),NA,NB,NIMAX,NBET)                     CALX-495
      NPAA=NPAR+NA                                                      CALX-496
      NSCN=NPAA+NVA+NB                                                  CALX-497
      IF (KMIN.EQ.0) KMIN=IABS(NW(2,1)-NW(1,2))/2+NW(1,1)               CALX-498
      IF (KMAX.EQ.0) KMAX=NJMAX                                         CALX-499
      LO(132)=LO(43).AND.((NW(1,NJIT).LT.KMAX).OR.(NW(2,NJIT).NE.0))    CALX-500
C INPUT OF POTENTIALS, DEFORMATIONS ......                              CALX-501
      NBT1=NBET                                                         CALX-502
      LMX=NIMAX/2+2                                                     CALX-503
      LMAX1=NJMAX+LMX                                                   CALX-504
      IF (LML.EQ.0) LML=LMAX1                                           CALX-505
      NW(2,5)=LML                                                       CALX-506
      NFIS=NSCN+7*NCONS                                                 CALX-507
      NGAM=NFIS+2*NFISS                                                 CALX-508
      NNIV=NGAM+NRD                                                     CALX-509
      NPOT=NNIV+(3*NCOLL*NCOLL+1)/2                                     CALX-510
      NBETA=NPOT+42*NPP                                                 CALX-511
      IF (LO(7)) NBETA=NPOT                                             CALX-512
      IF (NBETA+9*NBET.GT.IDMT) CALL MEMO('CALX',IDMT,NBETA+9*NBET)     CALX-513
      IF (DW(NWV+4).EQ.0.) NL2=NL3                                      CALX-514
      NLT=MAX0(NL2,NL3)                                                 CALX-515
      IF (.NOT.LO(65)) NLT=0                                            CALX-516
      CALL LECT(NBET,NW,NW(1,NIPP),DW(NWV),RM,IDMT-NBETA,LO,NW(1,NBETA),CALX-517
     1DW(NBETA),DW(NPOT),DW(NFIS),DW(NGAM),ISM,DW(NSCN))                CALX-518
      HCONV=ACONV*H*H                                                   CALX-519
      NFM=NBETA+9*NBET                                                  CALX-520
      IF (LO(74)) CALL HORA                                             CALX-521
C HELICITY AMPLITUDES AND OBSERVABLES.                                  CALX-522
      NJX=NCOLL                                                         CALX-523
      CALL DEPH(NCOLL,DW(NWV),NW,NGR,NPR,IDMT-NFM,LO,NW(1,NFM),KTGR,NW(1CALX-524
     1,NNIV),NJX,NJY)                                                   CALX-525
      JTH=0                                                             CALX-526
      IF (.NOT.LO(66)) JTH=IDINT((THETA2-THETA1)/DTHETA+1.5D0)          CALX-527
      JTX=JTX*JTH                                                       CALX-528
      NGRM=NFM+5*KCC                                                    CALX-529
      NNVI=NFM+5*KTGR                                                   CALX-530
      NTGX=1                                                            CALX-531
      NDONN=1                                                           CALX-532
      NREC=0                                                            CALX-533
      NTOT=0                                                            CALX-534
      IF (.NOT.LO(31)) GO TO 17                                         CALX-535
C NUMBER OF PARAMETERS IN SEARCH AND EXPERIMENTAL DATA.                 CALX-536
      READ (MR,1028) NCOLR,NREC,NFIT,NESSAI,YY(1),YY(2)                 CALX-537
C SEARCH CONDITIONS.                                                    CALX-538
C DEFECT VALUES FOR NESSAI (MAXIMUM NUMBER OF EVALUATIONS) AND ECH.     CALX-539
C NFIT IS THE NUMBER OF FUNCTIONS STORED FOR SEARCH BEYOND NREC+1.      CALX-540
      IF (NESSAI.EQ.0) NESSAI=100                                       CALX-541
      KFIT=NREC+NFIT+1                                                  CALX-542
      IF (YY(1).EQ.0.D0) YY(1)=20.D0                                    CALX-543
      IF (YY(2).LT.1.D0) YY(2)=1.D0                                     CALX-544
      YY(3)=0.D0                                                        CALX-545
      WRITE (MW,1029) NCOLR                                             CALX-546
      IF (NCOLR.LE.0) GO TO 15                                          CALX-547
      KTGR=KTGR+NCOLR                                                   CALX-548
      NTGX=NFM+5*KTGR                                                   CALX-549
      NDONN=NTGX+7*NCOLR                                                CALX-550
C EXPERIMENTAL DATA.                                                    CALX-551
      IF (NDONN.GT.IDMT) CALL MEMO('CALX',IDMT,NDONN)                   CALX-552
      CALL LECD(DW(NWV),NGR,KFIT,NESSAI,YY,NJY,IDMT-NDONN,LO,NW(1,NNVI),CALX-553
     1NW(1,NTGX),DW(NTGX),DW(NDONN),DW(NDONN),NW(1,NDONN),NMX)          CALX-554
      JTX=MAX0(JTH,JTX)                                                 CALX-555
      NNVI=NDONN+6*NTOT                                                 CALX-556
      IF (.NOT.LO(32)) GO TO 17                                         CALX-557
      IF (NREC.LE.0) GO TO 16                                           CALX-558
      NDE=NNVI                                                          CALX-559
      NISE=NDE+NREC                                                     CALX-560
      NRC=NISE+NMX/2                                                    CALX-561
      NIW=NRC+MAX0(14+NREC+KFIT*(NTOT+NREC+1),(NREC*(NREC+5))/2)        CALX-562
      IF (NIW.GT.IDMT) CALL MEMO('CALX',IDMT,NIW)                       CALX-563
      NW(1,NIW)=KFIT                                                    CALX-564
      NW(2,NIW)=NESSAI                                                  CALX-565
      NW(1,NIW+1)=1                                                     CALX-566
      NNVI=NIW+(KFIT+5)/2                                               CALX-567
      I1=51                                                             CALX-568
      IF (LO(76).OR.(.NOT.LO(75))) I1=59                                CALX-569
      DO 14 I=I1,65                                                     CALX-570
      LO(I+85)=LO(I)                                                    CALX-571
   14 LO(I)=.FALSE.                                                     CALX-572
      LO(116)=LO(75)                                                    CALX-573
      GO TO 17                                                          CALX-574
   15 WRITE (MW,1030)                                                   CALX-575
      LO(31)=.FALSE.                                                    CALX-576
   16 IF (LO(32)) WRITE (MW,1031)                                       CALX-577
      LO(32)=.FALSE.                                                    CALX-578
   17 IF (LO(74)) CALL HORA                                             CALX-579
      NNWI=NNVI                                                         CALX-580
      IF (LO(123)) NNWI=NNWI+(3*KAB*KAB+1)/2                            CALX-581
      NCC=NNWI+2*KAB*KAB                                                CALX-582
      MCC=NCC                                                           CALX-583
      IF (LO(123)) MCC=MCC+3*KAB                                        CALX-584
      NXA=MCC+3*KAB                                                     CALX-585
      NAM1=NXA                                                          CALX-586
      IF (LO(123)) NAM1=NAM1+KBC*KBC                                    CALX-587
      NAM2=NAM1+4*(NJY+16)                                              CALX-588
      IF (NAM2.GT.IDMT) CALL MEMO('CALX',IDMT,NAM2)                     CALX-589
C DO LOOPS AND CG-COEFFICIENTS FOR OBSERVABLES.                         CALX-590
      CALL OBSE(NW(1,NGRM),CW(1,NGRM),KTGR-KCC,NCOLR,NW,IDMT-NAM2,LO,CW(CALX-591
     11,NAM1),JCAL,NW(1,NAM1),CW(1,NAM1),DW(NAM2),CW(1,NAM2),NW(1,NAM2),CALX-592
     2DW(NAM2),NW(1,NAM2))                                              CALX-593
      NFA=4*LMAX1+10-2*NJMAX                                            CALX-594
      NFAC=NAM1+JCAL                                                    CALX-595
      NFG=NFAC+NFA+1                                                    CALX-596
      LMAX2=LMAX1                                                       CALX-597
      IF (.NOT.LO(108)) GO TO 18                                        CALX-598
      LMAX1=LMZ+LMX                                                     CALX-599
      IF (LMZ.EQ.0) LMAX1=LMAX2-NJMAX/2                                 CALX-600
      IF (LMZ.LT.0) LMAX1=LMAX2                                         CALX-601
      WRITE (MW,1032) LMAX2,LMAX1                                       CALX-602
   18 NXG=NFG+4*NJX*LMAX1                                               CALX-603
      NRES=NXG+NJX*LMAX2                                                CALX-604
      KE=0                                                              CALX-605
      NXX=NRES+NTOT                                                     CALX-606
      NT=NXX+NREC                                                       CALX-607
      ITEMM=ITERM                                                       CALX-608
      NPLACE=NT                                                         CALX-609
      IF (NPLACE.GT.IDMT) CALL MEMO('CALX',IDMT,NPLACE)                 CALX-610
      DW(NFAC)=0.D0                                                     CALX-611
      DO 19 I=1,NFA                                                     CALX-612
   19 DW(NFAC+I)=DW(NFAC+I-1)+DLOG(DFLOAT(I))                           CALX-613
      RETURN                                                            CALX-614
   20 WRITE (MW,1033) NSP(3),NCONT                                      CALX-615
      GO TO 22                                                          CALX-616
   21 WRITE (MW,1034) I,NW(K,1),L                                       CALX-617
   22 WRITE (MW,1035)                                                   CALX-618
      STOP                                                              CALX-619
 1000 FORMAT (18A4)                                                     CALX-620
 1001 FORMAT (50L1)                                                     CALX-621
 1002 FORMAT (14I5)                                                     CALX-622
 1003 FORMAT ('1',20X,'E. C. I. S. CODE FOR COUPLED CHANNELS CALCULATIONCALX-623
     1S ( 2006 ) WORKING FIELD LENGTH =',I10//24X,'TO OBTAIN THE INPUT DCALX-624
     2ESCRIPTION, PUNCH ''DESCRIPTION '' IN COLUMNS 1-12 OF THE FIRST CACALX-625
     3RD.'//)                                                           CALX-626
 1004 FORMAT ('  1',L2,' - ROTATIONAL MODEL-(.F.: VIBRATIONAL MODEL).',9CALX-627
     1X,'| 54',L2,' - OUTPUT OF THE LENGTH USED IN THE WORKING FIELD.'/'CALX-628
     2  2',L2,' - SECOND ORDER VIBRATIONAL MODEL-(.F.: FIRST ORDER)  | 5CALX-629
     35',L2,' - OUTPUT OF C-MATRIX ELEMENTS AND COMPOUND NUCLEUS'/8X,'ORCALX-630
     4 CONSTRAINED ASYMMETRIC ROTATIONAL MODEL.',8X,'|',8X,'INTERMEDIATECALX-631
     5 RESULTS.'/'  3',L2,' - ANHARMONIC VIBRATIONAL MODEL-(.F.: HARMONICALX-632
     6C) OR    | 56',L2,' - OUTPUT OF S-MATRIX ELEMENTS.'/8X,'ASYMMETRICCALX-633
     7 ROTATIONAL MODEL-(.F.: SYMMETRIC).',6X,'| 57',L2,' - OUTPUT OF PHCALX-634
     8ASE-SHIFTS AT EACH E.C.I.S. ITERATION.'/'  4',L2,' - PARAMETRISED CALX-635
     9SPIN-ORBIT DEFORMATION-(.F.:STANDARD) | 58',L2,' - OUTPUT OF THE CCALX-636
     AOEFFICIENTS OF EACH FORM FACTOR'/'  5',L2,' - DIFFERENT DEFORMATIOCALX-637
     BN FOR EACH POTENTIAL(.F.: SAME)|',8X,'FOR ALL SETS OF EQUATIONS.'/CALX-638
     C'  6',L2,' - USE OF DEFORMATION LENGTHS.',24X,'| 59',L2,' - TOTAL CALX-639
     DELASTIC REACTION CROSS-SECTIONS WRITTEN ON'/'  7',L2,' - NUCLEAR MCALX-640
     EATRIX ELEMENTS AND FORM FACTORS ON CARDS. |',8X,'FILE 58 AND TOTALCALX-641
     F INELASTIC REACTION'/'  8',L2,' - RELATIVISTIC KINEMATICS.',27X,'|CALX-642
     G',8X,'CROSS-SECTIONS ON FILE 59.'/'  9',L2,' - SYMMETRISED WOODS-SCALX-643
     HAXON FORM FOR NEGATIVE RADII.   | 60',L2,' - S-MATRIX ELEMENTS WRICALX-644
     ITTEN ON FILE 60.'/' 10',L2,' - DISPERSION RELATIONS FOR POTENTIALSCALX-645
     J.',15X,'|'/59X,'| 61',L2,' - REDUCED NUCLEAR MATRIX ELEMENTS WRITTCALX-646
     KEN ON FILE 61.')                                                  CALX-647
 1005 FORMAT (' 11',L2,' - DEFORMED COULOMB POTENTIAL.',24X,'| 62',L2,' CALX-648
     1- POTENTIALS WRITTEN ON FILE 62.'/' 12',L2,' - DEFORMED IMAGINARY CALX-649
     2POTENTIAL.',22X,'| 63',L2,' - PENETRABILITIES WRITTEN ON FILE 63.'CALX-650
     3/' 13',L2,' - DEFORMED REAL SPIN-ORBIT/TENSOR POTENTIAL. ',8X,'| 6CALX-651
     44',L2,' - RESULTS FOR EXPERIMENTAL DATA ON FILE 64'/' 14',L2,' - DCALX-652
     5EFORMED IMAGINARY SPIN-ORBIT/TENSOR POTENTIAL.    |',8X,'AND AT EQCALX-653
     6UIDISTANT ANGLES WRITTEN ON FILE 66.'/' 15',L2,' - REDUCED NUCLEARCALX-654
     7 MATRIX ELEMENTS READ FROM CARDS.   | 65',L2,' - LEGENDRE EXPANSIOCALX-655
     8N FOR CROSS-SECTIONS WRITTEN ON'/' 16',L2,' - HEAVY-ION DEFINITIONCALX-656
     9 OF RADII AND DEFORMATIONS.    |',8X,'FILE 65.'/' 17',L2,' - FOLDICALX-657
     ANG MODEL.',37X,'| 66',L2,' - NO CALCULATION AT EQUIDISTANT ANGLES.CALX-658
     B'/' 18',L2,' - PROJECTILE-TARGET ANTISYMMETRISATION.',14X,'| 67',LCALX-659
     C2,' - NO PLOT OF EXPERIMENTAL DATA.'/' 19',L2,' - DEFORMED COULOMBCALX-660
     D SPIN-ORBIT POTENTIAL.',13X,'| 68',L2,' - NO PLOT OF CROSS-SECTIONCALX-661
     ES AT EQUIDISTANT ANGLES.'/' 20',L2,' - DISPERSION RELATIONS FOR TRCALX-662
     FANSITION FORM-FACTORS.  | 69',L2,' - NO PLOT OF POLARISATIONS AT ECALX-663
     GQUIDISTANT ANGLES.'/59X,'|'/' 21',L2,' - USUAL COUPLED EQUATIONS-(CALX-664
     H.F.: ITERATIONS).',9X,'| 71',L2,' - NO DETAILED OUTPUT OF LOGICAL CALX-665
     ICONTROLS.'/' 22',L2,' - NO USE OF PADE FOR CONVERGENCE OF THE ITERCALX-666
     JATIONS.  | 72',L2,' - NO OUTPUT OF EXPERIMENTAL DATA AS THEY ARE RCALX-667
     KEAD.')                                                            CALX-668
 1006 FORMAT (' 23',L2,' - NO USE OF PADE WITHOUT CONVERGENCE.',16X,'| 7CALX-669
     13',L2,' - NO OUTPUT OF EXTERNAL POTENTIALS AS THEY ARE READ.'/' 24CALX-670
     2',L2,' - COUPLING POTENTIALS COMPUTED AT EACH ITERATION.    | 74',CALX-671
     3L2,' - OUTPUT OF VARIATIONS IN STORAGE.'/' 25',L2,' - COMPLETE CALCALX-672
     4CULATION UP TO THE END-(.F.: ONE ITE-  | 75',L2,' - NO COMPLETE OUCALX-673
     5TPUT AT THE FIRST RUN OF A SEARCH.'/8X,'RATION ONLY AS SOON AS TWOCALX-674
     6 ITERATIONS ARE ENOUGH). | 76',L2,' - LO(51) TO LO(65) ARE ALWAYS CALX-675
     7USED-(.F.: ONLY FOR'/' 26',L2,' - INTEGRATION STABILISED FOR LONG CALX-676
     8RANGE POTENTIAL.   |',8X,'COMPLETE OUTPUT)'/' 27',L2,' - NUMEROV''CALX-677
     9S METHOD FOR SINGLE EQUATIONS-(.F.: MODI-  | 77',L2,' - OUTPUT OF CALX-678
     ATIME DIFFERENCES DURING THE SEARCH.'/8X,'FIED NUMEROV''S METHOD).'CALX-679
     B,28X,'| 78',L2,' - OUTPUT OF DIFFERENCES BETWEEN EXPERIMENTAL'/' 2CALX-680
     C8',L2,' - COMPUTATION UP TO J-CONVERGENCE-(.F.: STOP WHEN    |',8XCALX-681
     D,'AND CALCULATED VALUES.'/8X,'ALL THE INHOMOGENEOUS TERMS ARE NEGLCALX-682
     EIGIBLE).',7X,'|'/' 29',L2,' - NO DIAGONAL TERM IN SECOND MEMBERS.'CALX-683
     F,16X,'| 81',L2,' - HAUSER-FESHBACH CORRECTIONS TO CROSS-SECTIONS.'CALX-684
     G/' 30',L2,' - PURE DWBA CALCULATION.',29X,'| 82',L2,' - SIMPLEST CCALX-685
     HOMPOUND NUCLEUS FORMALISM.'/59X,'| 83',L2,' - NO ENGELBRETCH-WEIDECALX-686
     INMULLER TRANSFORMATION IN'/' 31',L2,' - INPUT OF EXPERIMENTAL DATACALX-687
     J AND CHI2 CALCULATION.   |',8X,'IN COMPOUND NUCLEUS.')            CALX-688
 1007 FORMAT (' 32',L2,' - AUTOMATIC SEARCH ON SOME PARAMETERS.',15X,'| CALX-689
     184',L2,' - UNCOUPLED LEVELS FOR COMPOUND NUCLEUS.'/' 33',L2,' - SYCALX-690
     2MMETRISED CHI2 FOR CROSS-SECTIONS.',15X,'| 85',L2,' - FISSION DATACALX-691
     3 IN COMPOUND NUCLEUS.'/' 34',L2,' - NEXT CALCULATION CHANGING ENERCALX-692
     4GY AND/OR SOME',7X,'| 86',L2,' - GAMMA EMISSION DATA IN COMPOUND NCALX-693
     5UCLEUS.'/8X,'PARAMETERS.',40X,'| 87',L2,' - NO WIDTH FLUCTUATIONS.CALX-694
     6'/' 35',L2,' - SEARCH SAVED ON TAPE 8 IF ENDED OR LACK OF TIME.   CALX-695
     7|'/' 36',L2,' - RESTART A SEARCH FROM TAPE 8.',22X,'| 91',L2,' - ECALX-696
     8QUIDISTANT ANGLES IN THE LABORATORY SYSTEM.'/59X,'| 92',L2,' - NONCALX-697
     9 STANDARD OBSERVABLES AT EQUIDISTANT ANGLES.'/' 41',L2,' - FACTORICALX-698
     ASATION OF 1/(1-COS(THETA)) IN AMPLITUDES.   | 93',L2,' - NO RECOILCALX-699
     B CORRECTION FOR REACTIONS.'/' 42',L2,' - SCHMIDT''S ORTHOGONAL. INCALX-700
     C USUAL COUPLED EQUATIONS.  | 94',L2,' - NON RELATIVISTIC "REDUCED CALX-701
     DMASS" FOR DIRAC EQUATION.'/' 43',L2,' - INTERPOLATION ON TOTAL SPICALX-702
     EN.',23X,'| 95',L2,' - "REDUCED ENERGY" FOR NON COULOMB INTERACTIONCALX-703
     F IN'/' 44',L2,' - COULOMB CORRECTIONS.',31X,'|',8X,'RELATIVISTIC SCALX-704
     GCHROEDINGER EQUATION OR USE OF REST'/' 45',L2,' - ANGULAR MOMENTUMCALX-705
     H LIMIT FOR COULOMB CORRECTIONS.    |',8X,'MASS IN DIRAC EQUATION.'CALX-706
     I/' 46',L2,' - RESTRICTED COULOMB CORRECTIONS.',20X,'| 96',L2,' - "CALX-707
     JREDUCED MASS" FOR COULOMB INTERACTION IN')                        CALX-708
 1008 FORMAT (' 47',L2,' - NO RECOIL CORRECTION FOR BOUND STATES.',13X,'CALX-709
     1|',8X,'RELATIVISTIC SCHROEDINGER EQUATION.'/59X,'| 97',L2,' - SAMECALX-710
     2 REDUCED MASS FOR ALL THE STATES.'/' 51',L2,' - OUTPUT OF POTENTIACALX-711
     3LS.',30X,'| 98',L2,' - Q ADDED TO THE MASS OF RESIDUAL NUCLEUS.'/'CALX-712
     4 52',L2,' - OUTPUT OF REDUCED NUCLEAR MATRIX ELEMENTS.',9X,'| 99',CALX-713
     5L2,' - SCHROEDINGER EQUIVALENT TO DIRAC EQUATION.'/' 53',L2,' - OUCALX-714
     6TPUT OF THE NUMBER OF ITERATIONS.',16X,'|100',L2,' - COMPLETE DIRACALX-715
     7C EQUATION.'/'1')                                                 CALX-716
 1009 FORMAT (//' **** FIRST CONTROL CARD ****',2X,'1 ',9(' 1'),' 2 ',9(CALX-717
     1' 2'),' 3 ',9(' 3'),' 4 ',9(' 4'),' 5'/ 11X,5('  1 2 3 4 5 6 7 8 9CALX-718
     2 0')/11X,5(1X,10L2)//' *** SECOND CONTROL CARD ****',2X,'1 ',9(' 1CALX-719
     3'),' 2 ',9(' 2'),' 3 ',9(' 3'),' 4 ',9(' 4'),' 5'/11X,5('  1 2 3 4CALX-720
     4 5 6 7 8 9 0')/11X,5(1X,10L2)/)                                   CALX-721
 1010 FORMAT (10X,82('*')/10X,'*',80X,'*'/10X,'*',4X,18A4,4X,'*'/10X,'*'CALX-722
     1,80X,'*'/10X,82('*')/)                                            CALX-723
 1011 FORMAT (//' ******** YOU ARE USING DIRAC EQUATION WITHOUT RELATIVICALX-724
     1STIC KINEMATICS ********'//)                                      CALX-725
 1012 FORMAT (7F10.5)                                                   CALX-726
 1013 FORMAT (' MAXIMUM NUMBER OF J',I6,'  (MIN',I3,') STOP WHEN MAXIMUMCALX-727
     1 S-MATRIX ELEMENT IS LESS THAN',D18.4//' PLOT CONDITIONS FOR CROSSCALX-728
     2-SECTIONS',2I5,5X,'FOR POLARISATIONS',2I5/)                       CALX-729
 1014 FORMAT (4X,I3,' FACTORISATION OF 1/(1-COS(THETA)) IN THE AMPLITUDECALX-730
     1S.')                                                              CALX-731
 1015 FORMAT (' ANGULAR MOMENTA LIMITED TO',I5)                         CALX-732
 1016 FORMAT (10X,'ITERATION METHOD:  MAXIMUM NUMBER OF ITERATIONS',I10/CALX-733
     110X,'CONVERGENCE CRITERION:',1P,D15.4,' FOR S-MATRIX',D15.4,' FOR CALX-734
     2POTENTIALS AND FUNCTIONS.'/)                                      CALX-735
 1017 FORMAT (' IMAGINARY POTENTIAL INCREASED WITH A FACTOR',F10.5,' FORCALX-736
     1 BETTER CONVERGENCE.')                                            CALX-737
 1018 FORMAT (' USUAL COUPLED EQUATIONS.')                              CALX-738
 1019 FORMAT (' SCHMIDT''S ORTHOGONALISATION EVERY',I5,'  STEPS.')      CALX-739
 1021 FORMAT (' COULOMB CORRECTIONS LIMITED TO ANGULAR MOMENTA',I2,' FORCALX-740
     1 CENTRAL TERM AND',I2,' FOR SPIN-ORBIT TERM')                     CALX-741
 1022 FORMAT (/' INDICATIONS FOR EXPANSION OF CROSS-SECTIONS IN LEGENDRECALX-742
     1 POLYNOMIALS:',3I5)                                               CALX-743
 1023 FORMAT (8I5,2F10.5)                                               CALX-744
 1024 FORMAT (/' COMPOUND NUCLEUS INPUT:',I5,' UNCOUPLED STATES',I4,' WICALX-745
     1TH ANGULAR DISTRIBUTION AND',I4,' WITHOUT ANGULAR DISTRIBUTION'/24CALX-746
     2X,I5,' FISSION TRANSM. COEFF.'/24X,I5,' GAMMA TRANSM. COEFF.'/24X,CALX-747
     3I5,' CONTINUA')                                                   CALX-748
 1025 FORMAT (' MAXIMUM NUMBER OF SPIN VALUES IN THE CONTINUUM:',I5/' DICALX-749
     1SCRETISATION WITH:',2F10.5)                                       CALX-750
 1026 FORMAT (' INTERPOLATION OF S-MATRIX',3(2X,'FROM',I6,'  BY STEPS OFCALX-751
     1',I4,'+1')/(26X,3(2X,'FROM',I6,'  BY STEPS OF',I4,'+1')))         CALX-752
 1027 FORMAT (' UNCUMULATED INCREASES  ',3(27X,I4)/(24X,3(27X,I4)))     CALX-753
 1028 FORMAT (4I5,2F10.5)                                               CALX-754
 1029 FORMAT (//5X,I5,'  EXPERIMENTAL ANGULAR DISTRIBUTIONS.'/)         CALX-755
 1030 FORMAT (' THERE ARE NO EXPERIMENTAL DATA   ...  NO SEARCH.')      CALX-756
 1031 FORMAT (' NO PARAMETER IN SEARCH.')                               CALX-757
 1032 FORMAT (2X,I10,' COULOMB PHASES AND INDEFINITE INTEGRALS'/2X,I10,'CALX-758
     1 COULOMB FUNCTIONS AND FINITE INTEGRALS.')                        CALX-759
 1033 FORMAT (' NUMBER OF UNCOUPLED STATES WITHOUT ANGULAR DISTRIBUTION'CALX-760
     1,I5,' LESS THAN THE NUMBER OF CONTINUA:',I5)                      CALX-761
 1034 FORMAT (2X,I3,'TH LIMIT OF INTERPOLATION',I6,' SMALLER THAN PREVIOCALX-762
     1US ONE',I6)                                                       CALX-763
 1035 FORMAT (//' IN CALX  ...  STOP  ...')                             CALX-764
      END                                                               CALX-765
