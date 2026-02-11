C 23/02/07                                                      ECIS06  CALC-000
      SUBROUTINE CALC(NW,CW,DW,IDMX)                                    CALC-001
C MAIN SUBROUTINE OF THE PROGRAM.                                       CALC-002
C INPUT:     NW:      WORKING AREA FOR INTEGERS.                        CALC-003
C            CW:      WORKING AREA FOR CHARACTERS.                      CALC-004
C            DW:      WORKING AREA FOR DOUBLE PRECISION.                CALC-005
C            IDMX:    LENGTH OF DW.                                     CALC-006
C               NW,CW,DW ARE IN EQUIVALENCE BY CALL.                    CALC-007
C                                                                       CALC-008
C THE COMMON /ADDRE/ IS USED IN CALC, CALX, CAL1, VARI, EVAL AND REST.  CALC-009
C THE COMMON /COUPL/ IS USED IN CALC, CALX, LECL, LECT, REDM, VIBM,     CALC-010
C                               ROTM, ROAM, EXTP, LDIS, CAL1, POTE,     CALC-011
C                               ROTP, VARI, EVAL AND REST.              CALC-012
C THE COMMON /DCONS/ IS USED IN ECIS, CALC, LECL, LECT, COLF, KHCO,     CALC-013
C                               CONU, CAL1, POTE, ROTP, ROTZ, STDP,     CALC-014
C                               STBF, MTCH, SCAM, SCHE, LCSP, RESU,     CALC-015
C                               SCAT AND REST.                          CALC-016
C THE COMMON /INTEG/ IS USED IN CALC, CALX, DEPH, LECD, CAL1, VARI,     CALC-017
C                               EVAL AND REST.                          CALC-018
C THE COMMON /TITRE/ IS USED IN CALC, CALX, RESU, EVAL AND REST.        CALC-019
C                                                                       CALC-020
C FOR THE COMMON  /CONVE/, /DCHI2/ AND /NCOMP/ SEE CALX.                CALC-021
C FOR THE COMMONS /POTE1/ AND /POTE2/ SEE REDM.                         CALC-022
C                                                                       CALC-023
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ADDRE/:                     CALC-024
C  NIPH:      NUMBER OF PHONONS FOR THE HARMONIC VIBRATIONAL MODEL.     CALC-025
C  NJIT:      DATA FOR INTERPOLATION ON TOTAL SPIN.                     CALC-026
C  NWV:       NON INTEGER VALUES FOR THE CHANNELS.                      CALC-027
C  NIPP:      DISPERSION PARAMETERS.                                    CALC-028
C  NPAR:      INDICATIONS FOR NUCLEAR PARAMETERS.                       CALC-029
C  NPAA:      VALUES OF NUCLEAR PARAMETERS.                             CALC-030
C  NSCN:      LEVEL DENSITY DESCRIPTION.                                CALC-031
C  NFIS:      FISSION DATA FOR COMPOUND NUCLEUS.                        CALC-032
C  NGAM:      GAMMA DATA FOR COMPOUND NUCLEUS.                          CALC-033
C  NNIV:      ADDRESSES OF REDUCED NUCLEAR MATRIX ELEMENTS AND COULOMB  CALC-034
C             INTEGRALS IN NIV(NCOLL,NCOLL,3).                          CALC-035
C  NPOT:      OPTICAL POTENTIAL PARAMETERS.                             CALC-036
C  NBETA:     DEFORMATION PARAMETERS.                                   CALC-037
C  NFM:       HELICITIES AND OBSERVABLES (SEE DEPH,LECD AND OBSE).      CALC-038
C  NTGX:      BEGINNING OF CHI2 AND NORMALISATION OF DATA.              CALC-039
C  NDONN:     EXPERIMENTAL DATA.                                        CALC-040
C  NNVI:      TABLE OF ADDRESSES OF COUPLING COEFFICIENTS.              CALC-041
C  NDE:       SEARCH ACCURACIES.                                        CALC-042
C  NISE:      INDEXES OF THE VARIABLE PARAMETERS IN SEARCH.             CALC-043
C  NRC:       PERMANENT WORKING FIELD FOR THE SEARCH.                   CALC-044
C  NIW:       INTEGER WORKING FIELD FOR THE SEARCH.                     CALC-045
C  NNWI:      SAME AS NNVI FOR SYMMETRISED EQUATIONS.                   CALC-046
C  NCC:       TABLE OF ENERGIES, L*(L+1) AND L*S  (SEE QUAN).           CALC-047
C  MCC:       SAME AS NCC FOR SYMMETRISED EQUATIONS.                    CALC-048
C  NXA:       TABLE OF COEFFICIENTS OF SYMMETRISATION.                  CALC-049
C  NAM1:      COMPUTATION OF OBSERVABLES (SEE OBSE).                    CALC-050
C  NFAC:      TABLE OF LOG OF FACTORIALS FOR GEOMETRIC COEFFICIENTS.    CALC-051
C  NFG:       COULOMB FUNCTIONS AND FINITE INTEGRALS.                   CALC-052
C  NXG:       COULOMB PHASES AND INFINITE INTEGRALS.                    CALC-053
C  NRES:      FUNCTIONS FOR THE SEARCH.                                 CALC-054
C  NXX:       VARIABLES FOR THE SEARCH.                                 CALC-055
C  NT:        TABLE OF REDUCED NUCLEAR MATRIX ELEMENTS.                 CALC-056
C  NIVQ:      TABLE OF MULTIPOLES.                                      CALC-057
C  NIVY:      TABLE OF FORM FACTOR IDENTIFICATION IVY (FOR COMPUTATION).CALC-058
C  NIVZ:      TABLE OF FORM FACTOR IDENTIFICATION IVZ (FOR USE).        CALC-059
C  NCOI:      ADDRESS OF THE TABLE FOR DISCRETISATION OF CONTINUUM.     CALC-060
C  MIPI:      ADDRESS OF "IPIM" FOR THE DISCRETISATION OF CONTINUUM.    CALC-061
C  NXD:       ADDRESS OF WEIGHTS AND STEPS OF CONTINUUM.                CALC-062
C  MWV:       SAME AS NWV FOR THE DISCRETISATION OF CONTINUUM.          CALC-063
C  NIXT:      TRANSMISSION COEFFICIENTS OF UNCOUPLED STATES.            CALC-064
C  NTY:       TEMPORARY RESULTS OF COMPOUND NUCLEUS FOR INTERPOLATION.  CALC-065
C  NSM:       STANDARD AND HELICITY SCATTERING MATRIX ELEMENTS.         CALC-066
C  NTX:       PARTIAL REACTION CROSS SECTIONS AND COMPOUND NUCLEUS.     CALC-067
C  NRY:       COMPOUND NUCLEUS COEFFICIENTS.                            CALC-068
C  NRCO:      STRENGTH OF COULOMB CENTRAL POTENTIALS FOR CORRECTIONS.   CALC-069
C  NRDO:      STRENGTH OF COULOMB TRANSITION POTENTIALS FOR CORRECTIONS.CALC-070
C  NVC1:      REAL POTENTIALS.                                          CALC-071
C  NVC2:      IMAGINARY POTENTIALS.                                     CALC-072
C  NNC:       FIRST FREE ADDRESS AFTER COMPUTATION OF POTENTIALS.       CALC-073
C  NCX:       FIRST FREE ADDRESS FOR COMPUTATION OF POTENTIALS.         CALC-074
C   DEFINED:  NPOT,NT,NIVQ,NIVY,NIVZ,NCOI,MIPI,NXD,MWV,NIXT,NTY,NSM,NTX,CALC-075
C             NRY,,NRDO,NVC1,NVC2,NNC,NCX.                              CALC-076
C   USED:     NIPH,NWV,NIPP,NPAR,NPAA,NSCN,NNIV,NPOT,NBETA,NDE,NRC,NIW, CALC-077
C             NFAC,NFG,NXG,NRES,NXX,NT,NIVQ,NIVY,NIVZ,NCOI,MIPI,NXD,MWV,CALC-078
C             NIXT,NTY,NSM,NTX,NRY,NRCO,NRDO,NVC1,NNC,NCX.              CALC-079
C   NOT USED: NJIT,NFIS,NGAM,NFM,NTGX,NDONN,NNVI,NISE,NNWI,NCC,MCC,NXA, CALC-080
C             NAM1.                                                     CALC-081
C                                                                       CALC-082
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     CALC-083
C  H:         STEP SIZE FOR INTEGRATION.                                CALC-084
C   USED:     H                                                         CALC-085
C                                                                       CALC-086
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     CALC-087
C  IQM:       MAXIMUM L-VALUE OF DEFORMATION IN ROTATIONAL MODEL.       CALC-088
C  IQMAX:     MAXIMUM L-EXPANSION IN ROTATIONAL MODEL.                  CALC-089
C  NBT1:      NUMBER OF PHONONS (VIBRATIONS).                           CALC-090
C  NFA:       NUMBER OF LOGARITHMS OF FACTORIALS.                       CALC-091
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             CALC-092
C  NSPIN:     TWICE THE K-VALUE OF THE ROTATIONAL BAND.                 CALC-093
C  NVA:       NUMBER OF NUCLEAR PARAMETERS.                             CALC-094
C   USED:     IQMAX,NFA,NPP,NSPIN                                       CALC-095
C   NOT USED: IQM,NBT1,NVA                                              CALC-096
C                                                                       CALC-097
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCHI2/:                     CALC-098
C  CHI2M:     MINIMUM CHI-SQUARE IN THE SEARCH.                         CALC-099
C   DEFINED:  CHI2M.                                                    CALC-100
C                                                                       CALC-101
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     CALC-102
C  CM:        ATOMIC MASS IN MEV.                                       CALC-103
C  CHB:       PLANCK CONSTANT /(2*PI) IN MEV*FERMI.                     CALC-104
C  CZ:        ELECTRIC CONSTANT.                                        CALC-105
C  CMB:       ATOMIC MASS CM DIVIDED BY H-BAR*C.                        CALC-106
C  CCZ:       COULOMB ALPHA CONSTANT.                                   CALC-107
C  CK:        H-BAR*C.                                                  CALC-108
C  XZ:        CONVERSION FACTOR TO MILLIBARNS.                          CALC-109
C   DEFINED:  CMB,CCZ,CK.                                               CALC-110
C   USED:     CM,CHB,CZ.                                                CALC-111
C   NOT USED: XZ.                                                       CALC-112
C                                                                       CALC-113
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /INTEG/:                     CALC-114
C  IDMT:      TOTAL WORKING FIELD LENGTH AS SINGLE PRECISION.           CALC-115
C  NPLACE:    MAXIMUM WORKING FIELD USED.                               CALC-116
C  NCOLL:     NUMBER OF COUPLED CHANNELS.                               CALC-117
C  NJMAX:     MAXIMUM NUMBER OF J-VALUES.                               CALC-118
C  ITERM:     MAXIMUM NUMBER OF ITERATIONS.                             CALC-119
C  JDM:       MINIMUM NUMBER OF TOTAL SPIN.                             CALC-120
C  JIT:       NUMBER OF RATES OF INTERPOLATION ON TOTAL SPIN.           CALC-121
C  KMIN:      MINIMUM NUMBER OF J VALUES FOR COMPOUND NUCLEUS.          CALC-122
C  KMAX:      MAXIMUM NUMBER OF J VALUES FOR COMPOUND NUCLEUS.          CALC-123
C  NCOLS:     NUMBER OF CHANNELS WITH ANGULAR DISTRIBUTIONS.            CALC-124
C  NCOLT:     NUMBER OF CHANNELS INCLUDING UNCOUPLED STATES.            CALC-125
C  NBET:      NUMBER OF DIFFERENT DEFORMATIONS (VIBRATIONS+ROTATIONS).  CALC-126
C  LMX:       DIFFERENCE BETWEEN NUMBER OF J AND OF COULOMB FUNCTIONS.  CALC-127
C  LMAX1:     NUMBER OF L-VALUES FOR COULOMB FUNCTIONS.                 CALC-128
C  NLT:       MEMORIES NEEDED FOR LEGENDRE POLYNOMIALS.                 CALC-129
C  ISM:       NUMBER OF INTEGRATION STEPS.                              CALC-130
C  NJC:       MAXIMUM NUMBER OF OBSERVABLES AT EQUIDISTANT ANGLES.      CALC-131
C  JTX:       MAXIMUM NUMBER OF CALCULATED VALUES FOR A PLOT.           CALC-132
C  KCC:       NUMBER OF INDEPENDENT AMPLITUDES WITH UNCOUPLED STATES.   CALC-133
C  MS1:       LARGEST PARTICLE MULTIPLICITY.                            CALC-134
C  MS2:       LARGEST TARGET MULTIPLICITY.                              CALC-135
C  KBA:       NUMBER OF INDEPENDENT AMPLITUDES WITHOUT UNCOUPLED STATES.CALC-136
C  KAB:       MAXIMUM NUMBER OF EQUATIONS.                              CALC-137
C  KBC:       MAXIMUM NUMBER OF SOLUTIONS.                              CALC-138
C  JTH:       MAXIMUM NUMBER OF ANGLES FOR A PLOT.                      CALC-139
C  NCOLR:     NUMBER OF EXPERIMENTAL ANGULAR DISTRIBUTIONS.             CALC-140
C  NREC:      NUMBER OF VARIABLES IN SEARCH.                            CALC-141
C  NTOT:      NUMBER OF EXPERIMENTAL DATA.                              CALC-142
C  LMAX2:     NUMBER OF L VALUES FOR COULOMB PHASES.                    CALC-143
C  KE:        CONTROL OF SEARCH (SEE FITE).                             CALC-144
C  ITEMM:     MAXIMUM NUMBER OF ITERATIONS AT THE BEGINNING OF A RUN.   CALC-145
C  KXT:       NUMBER OF PENETRABILITIES FOR UNCOUPLED STATES.           CALC-146
C  LMAX3:     EFFECTIVE NUMBER OF COULOMB FUNCTIONS.                    CALC-147
C  NRZ:       NUMBER OF RESULTS TO SAVE FOR MINIMUM CHI2..              CALC-148
C  NTZ:       NUMBER OF MEMORIES TO INITIALISE TO ZERO FOR SCATTERING   CALC-149
C             MATRIX AND COMPOUND NUCLEUS RESULTS..                     CALC-150
C  IPM:       NUMBER PF J VALUES USED FOR SCATTERING MATRIX (IN CAL1).  CALC-151
C  IPK:       NUMBER PF J VALUES USED FOR COMPOUND NUCLEUS (IN CAL1).   CALC-152
C  MCM(1):    MAXIMUM ANGULAR MOMENTUM FOR CENTRAL COULOMB CORRECTIONS. CALC-153
C  MCM(2):    MAXIMUM ANGULAR MOMENTUM FOR SPIN-ORBIT COULOMB CORRECT.  CALC-154
C  NCT(1):    NUMBER OF EQUATIONS FOR POSITIVE PARITY.                  CALC-155
C  NCT(2):    NUMBER OF EQUATIONS FOR NEGATIVE PARITY.                  CALC-156
C  NCT(3):    NUMBER OF SOLUTIONS FOR POSITIVE PARITY.                  CALC-157
C  NCT(4):    NUMBER OF SOLUTIONS FOR NEGATIVE PARITY.                  CALC-158
C  NCT(5):    NUMBER OF COMPOUND NUCLEUS COEFFICIENTS, POSITIVE PARITY. CALC-159
C  NCT(6):    NUMBER OF COMPOUND NUCLEUS COEFFICIENTS, NEGATIVE PARITY. CALC-160
C   DEFINED:  IDMT,NPLACE,KE,NRZ,NTZ.                                   CALC-161
C   USED:     IDMT,NPLACE,NCOLL,NJMAX,KMAX,NCOLS,NCOLT,NBET,LMAX1,ISM,  CALC-162
C             KCC,KBA,NREC,NTOT,LMAX2,KE,KXT,LMAX3.                     CALC-163
C   NOT USED: ITERM,JDM,JIT,KMIN,LMX,NLT,NJC,JTX,MS1,MS2,KAB,KBC,JTH,   CALC-164
C             NCOLR,ITMM,IPM,IPK,MCM,NCT.                               CALC-165
C                                                                       CALC-166
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     CALC-167
C  NSP(1):    NUMBER OF UNCOUPLED LEVELS FOR COMPOUND NUCLEUS           CALC-168
C             INCLUDING DISCRETISATION OF CONTINUUM.                    CALC-169
C  NSP(2):    NUMBER OF THESE LEVELS WITH ANGULAR DISTRIBUTION.         CALC-170
C  NSP(3):    NUMBER OF THESE LEVELS WITHOUT ANGULAR DISTRIBUTION.      CALC-171
C  NIE:       NUMBER OF UNCOUPLED STATES ADDED FOR DISCRETISATION.      CALC-172
C  NCOLX:     TOTAL NUMBER OF LEVELS WITHOUT DISCRETISATION.            CALC-173
C  NDP:       ADDRESS OF FISSION AND GAMMA FINAL RESULTS.               CALC-174
C  NDQ:       ADDRESS OF FISSION AND GAMMA INTERMEDIATE RESULTS.        CALC-175
C   DEFINED:  NSP,NIE,NCOLX,NDP,NDQ.                                    CALC-176
C   USED:     NSP,NIE,NCOLX.                                            CALC-177
C                                                                       CALC-178
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE1/:                     CALC-179
C  ITX(16):   STARTING ADDRESS OF DIFFERENT FORM FACTORS.               CALC-180
C             FOR SCHROEDINGER EQUATION, ITX(I)+1 IS THE STARTING       CALC-181
C             ADDRESS OF THE FORM FACTOR READ IN EXTP WITH ITYP=I       CALC-182
C             (POTENTIAL FOR I=1 TO 8, TRANSITION FOR I=9 TO 16).       CALC-183
C             FOR DIRAC EQUATIONS, ITX(1)=0,                            CALC-184
C             ITX(2)+1=ADDRESS OFF FIRST TRANSITION FORM FACTOR,        CALC-185
C             ITX(7)=ADDRESS OF LAST TRANSITION FORM FACTOR,            CALC-186
C             ITX(3)=ADDRESS OF LAST TEMPORARY CENTRAL POTENTIAL,       CALC-187
C             ITX(4)=ITX(7)-24,ITX(5)=ITX(3)-11,ITX(6)=ITX(2)-4.        CALC-188
C             ALL ARE USED FOR SCHROEDINGER, THE FIRST 8 FOR DIRAC.     CALC-189
C  INTC:      NUMBER OF FORM FACTORS WITHOUT DEFORMED SPIN-ORBIT        CALC-190
C             INCLUDING CORRECTION TERMS.                               CALC-191
C  ITXM:      TOTAL NUMBER OF FORM FACTORS.                             CALC-192
C   USED:     ITX(1),ITX(2),INTC,ITXM.                                  CALC-193
C                                                                       CALC-194
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     CALC-195
C  ITY(2):    STARTING ADDRESS OF IMAGINARY CENTRAL POTENTIAL.          CALC-196
C        ITY(2)=14*NCOLL IS USED FOR DIRAC EQUATIONS.                   CALC-197
C  INTV:      SAME AS INVT, TAKING INTO ACCOUNT DISPERSION.             CALC-198
C  NPX:       NUMBER OF POTENTIALS TAKING INTO ACCOUNT DISPERSION.      CALC-199
C   USED:     ITY(2),INTV,NPX.                                          CALC-200
C                                                                       CALC-201
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /TITRE/:                     CALC-202
C  TITLE(18): TITLE OF THE RUN PRINTED AS HEADING OF RESULTS.           CALC-203
C   USED:     TITLE.                                                    CALC-204
C                                                                       CALC-205
C ******** MEANING OF THE LOGICAL CONTROLS LO ********                  CALC-206
C                                                                       CALC-207
C LO(135+I)=LO(50+I) FOR I=1,15 DURING SEARCH IF LO(51) TO LO(65) ARE   CALC-208
C USED ONLY FOR COMPLETE OUTPUT.                                        CALC-209
C THE VALUES OF THE FIRST 100 LO ARE READ IN CALX ON TWO DATA CARDS.    CALC-210
C THE FIRST DATA CARD IS FOR 1-50,THE SECOND FOR 51-100.                CALC-211
C ONLY,THE VALUES LISTED ON THE INPUT DESCRIPTION ARE USED.             CALC-212
C *** FOR THEIR MEANING, SEE THIS DESCRIPTION *****                     CALC-213
C EXCEPTION IF LO(36)=.TRUE.: ONLY THE FIRST CARD IS READ, THE CODE DOESCALC-214
C NOT TAKE THEM INTO ACCOUNT BUT CONTINUE A SEARCH SAVED ON TAPE MS.    CALC-215
C                                                                       CALC-216
C *** MEANING OF THE 100 FIRST LO USED HERE ********                    CALC-217
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    CALC-218
C               LO(10) =.TRUE. DISPERSION RELATIONS FOR POTENTIALS.     CALC-219
C               LO(32) =.TRUE. AUTOMATIC SEARCH ON SOME PARAMETERS.     CALC-220
C               LO(34) =.TRUE. NEXT CALCULATION CHANGING ENERGY AND/OR  CALC-221
C                              SOME PARAMETERS.                         CALC-222
C               LO(35) =.TRUE. SEARCH SAVED ON TAPE 8 IF CORRECTLY ENDEDCALC-223
C                              OR STOPPED BY THE NUMBER OF EVALUATIONS. CALC-224
C               LO(36) =.TRUE. RESTART A SEARCH FROM TAPE 8.            CALC-225
C               LO(52) =.TRUE. REDUCED NUCLEAR MATRIX ELEMENTS OUTPUT.  CALC-226
C               LO(54) =.TRUE. OUTPUT LENGTH USED IN THE WORKING FIELD. CALC-227
C               LO(61) =.TRUE. REDUCED NUCLEAR MATRIX ELEMENTS WRITTEN  CALC-228
C                              ON FILE 61.                              CALC-229
C               LO(75) =.TRUE. NO COMPLETE OUTPUT AT THE FIRST RUN OF A CALC-230
C                              SEARCH.                                  CALC-231
C               LO(76) =.TRUE. LO(51) TO LO(65) ARE ALWAYS USED.        CALC-232
C               LO(77) =.TRUE. OUTPUT OF TIME DIFFERENCES DURING THE    CALC-233
C                              SEARCH.                                  CALC-234
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             CALC-235
C               LO(100)=.TRUE. DIRAC EQUATION.                          CALC-236
C                                                                       CALC-237
C *** MEANING OF LO(I) FOR I GREATER THAN 100 ****                      CALC-238
C LO(101) IS TRUE IF THERE IS A REAL SPIN-ORBIT POTENTIAL.              CALC-239
C LO(102) IS TRUE IF THERE IS A IMAGINARY SPIN-ORBIT POTENTIAL.         CALC-240
C LO(103) IS TRUE IF THERE IS A COULOMB SPIN-ORBIT POTENTIAL.           CALC-241
C LO(104) IS TRUE IF CONVERGENCE IS OBTAINED IN THE ITERATION.          CALC-242
C LO(105) IS TRUE IF CONVERGENCE IS OBTAINED FOR THIS EQUATION.         CALC-243
C LO(106) IS TRUE WHEN THE ITERATION IS NOT THE LAST ONE PERMITTED.     CALC-244
C LO(107) IS TRUE IF ALL THE COUPLINGS HAVE TO BE CALCULATED BEFOREHAND.CALC-245
C LO(108) IS TRUE IF THE DIAGONAL COULOMB CORRECTIONS ARE NEEDED.       CALC-246
C LO(109) IS TRUE FOR DIRAC POTENTIALS.                                 CALC-247
C LO(110) IS TRUE IF DERIVATIVES ARE NEEDED.                            CALC-248
C LO(111) IS TRUE IF NUCLEAR PARAMETERS ARE CHANGED IN SEARCH.          CALC-249
C LO(112) IS TRUE IF SPIN-ORBIT OR COMPOUND NUCLEUS PARAMETERS ARE      CALC-250
C                 CHANGED IN SEARCH.                                    CALC-251
C LO(113) IS TRUE IF DISPERSION RELATION IS CHANGED IN SEARCH.          CALC-252
C LO(114) IS TRUE IF COMPOUND NUCLEUS CONTINUUM IS CHANGED IN SEARCH.   CALC-253
C LO(115) IS TRUE IF IT IS THE FIRST COMPUTATION FOR THIS ENERGY.       CALC-254
C LO(116) IS TRUE FOR NO OUTPUT.                                        CALC-255
C LO(117) IS TRUE FOR ALL THE CALCULATIONS EXCEPT THE FIRST.            CALC-256
C LO(118) IS TRUE FOR LAST RESULTS.                                     CALC-257
C LO(119) IS TRUE FOR RESULTS WITHOUT DOING THE CALCULATION AGAIN.      CALC-258
C LO(120) IS TRUE FOR OUTPUT AND LAST CALCULATION IS THE BEST ONE.      CALC-259
C LO(121) IS TRUE FOR OPTICAL MODEL WITHOUT COUPLING.                   CALC-260
C LO(122) IS TRUE FOR IDENTICAL PARTICLES WITHOUT SPIN.                 CALC-261
C LO(123) IS TRUE FOR IDENTICAL PARTICLES WITH SPIN.                    CALC-262
C LO(124) IS TRUE TO COMPUTE TRANSMISSION COEFFICIENTS.                 CALC-263
C LO(125) IS TRUE IN CAL1 FOR USUAL COUPLED EQUATIONS SUBROUTINES.      CALC-264
C LO(126) IS TRUE IF THERE ARE OBSERVABLES IN THE LABORATORY SYSTEM.    CALC-265
C LO(127) IS TRUE FOR COULOMB CORRECTIONS WITH PURE REGULAR FUNCTIONS.  CALC-266
C LO(128) IS TRUE FOR NO COPY OF UNCOUPLED FUNCTIONS AND PHASE-SHIFT.   CALC-267
C LO(129) IS TRUE FOR REAL SPIN-ORBIT OR DIRAC EQUATION.                CALC-268
C LO(130) IS TRUE FOR IMAGINARY SPIN-ORBIT OR DIRAC EQUATION.           CALC-269
C LO(131) IS TRUE IF THE TOTAL SPIN IS TOO LARGE FOR COMPOUND NUCLEUS.  CALC-270
C LO(132) IS TRUE STORE FISSION AND GAMMA TRANSMISSION COEFFICIENTS FOR CALC-271
C                 INTERPOLATION.                                        CALC-272
C LO(133) IS TRUE TO STORE SCALAR AND COULOMB POTENTIAL INDEPENDENTLY.  CALC-273
C LO(134) AND LO(135) ARE NOT USED.                                     CALC-274
C LO(111) TO LO(120) ARE INITIALISED TO .FALSE. IN CALX.                CALC-275
C LO(116) TO LO(120) ARE INITIALISED TO .FALSE. IN EVAL.                CALC-276
C                                                                       CALC-277
C***********************************************************************CALC-278
      IMPLICIT REAL*8 (A-H,O-Z)                                         CALC-279
      LOGICAL LO(150)                                                   CALC-280
      DIMENSION NW(2,IDMX),DW(IDMX)                                     CALC-281
      CHARACTER*4 FIN,TITLE,CW(2*IDMX)                                  CALC-282
      COMMON /ADDRE/ NIPH,NJIT,NWV,NIPP,NPAR,NPAA,NSCN,NFIS,NGAM,NNIV,NPCALC-283
     1OT,NBETA,NFM,NTGX,NDONN,NNVI,NDE,NISE,NRC,NIW,NNWI,NCC,MCC,NXA,NAMCALC-284
     21,NFAC,NFG,NXG,NRES,NXX,NT,NIVQ,NIVY,NIVZ,NCOI,MIPI,NXD,MWV,NIXT,NCALC-285
     3TY,NSM,NTX,NRY,NRCO,NRDO,NVC1,NVC2,NNC,NCX                        CALC-286
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       CALC-287
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   CALC-288
      COMMON /DCHI2/ CHI2,CHI2M,YY(3)                                   CALC-289
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            CALC-290
      COMMON /INOUT/ MR,MW,MS                                           CALC-291
      COMMON /INTEG/ IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCOCALC-292
     1LS,NCOLT,NBET,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,JTCALC-293
     2H,NCOLR,NREC,NTOT,LMAX2,KE,ITEMM,KXT,LMAX3,NRZ,NTZ,IPM,IPK,MCM(2),CALC-294
     3NCT(6)                                                            CALC-295
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQCALC-296
     1,ACN1,ACN2,AZ(18)                                                 CALC-297
      COMMON /POTE1/ ITX(16),IMAX,INTC,INLS,INVC,INVD,ITXM              CALC-298
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         CALC-299
      COMMON /TITRE/ TITLE(18)                                          CALC-300
      DATA FIN /'FIN '/                                                 CALC-301
      IDMT=IDMX                                                         CALC-302
      CMB=CM/CHB                                                        CALC-303
      CCZ=CHB/CZ                                                        CALC-304
      CK=2.D0*CM/CHB**2                                                 CALC-305
C MAIN INPUT.                                                           CALC-306
    1 CHI2M=1.D35                                                       CALC-307
      CALL CALX(NW,CW,DW,LO)                                            CALC-308
      IF (TITLE(1).EQ.FIN) RETURN                                       CALC-309
      NSP1D=NSP(1)                                                      CALC-310
      IF (LO(36)) GO TO 17                                              CALC-311
    2 NSP(1)=NSP1D                                                      CALC-312
      NSP(3)=NSP(1)-NSP(2)                                              CALC-313
      CALL COLF(NCOLX,NCOLL,NW,DW(NWV),ISM,LMAX1,LMAX2,H,IDMT-NT,NW(1,NNCALC-314
     1IV),LO,DW(NFG),DW(NXG),LMAX3,KXT,NW(1,NT),DW(NT))                 CALC-315
C COMPUTATION OF NUCLEAR MATRIX ELEMENTS AND SPACE FOR FORM FACTORS.    CALC-316
    3 IF (LO(61)) WRITE (61,1000) DW(NWV),DW(NWV+12),DW(NWV+1),NW(2,2),NCALC-317
     1COLL                                                              CALC-318
      CALL REDM(NW,NCOLL,NW(1,NBETA),DW(NBETA),NW(1,NIPH),NW(1,NPAR),DW(CALC-319
     1NPAA),DW(NFAC),IDMT-NT,LO,NW(1,NNIV),NW(1,NT),NW(1,NT),DW(NT),IT,ICALC-320
     2M)                                                                CALC-321
      NIVQ=NT+3*IT                                                      CALC-322
      NIVY=NIVQ+(3*IM+1)/2                                              CALC-323
      NIVZ=NIVY+(7*INTC+1)/2                                            CALC-324
      NCOI=NIVZ+2*INTV                                                  CALC-325
      IF ((.NOT.LO(7)).OR.LO(117)) GO TO 6                              CALC-326
      NPOT=NCOI                                                         CALC-327
      CALL EXTP(NCOLL,NCOLT,DW(NWV),NW(1,NNIV),NW(1,NT),NW(1,NIVY),NW(1,CALC-328
     1NIVQ),NW,DW(NFAC),DW(NPOT),NW(1,NPOT),NW(1,NPOT+2),NW(1,NIPP),DW(NCALC-329
     2IPP),IDMT-NPOT,LO)                                                CALC-330
C PERMUTATION OF POTENTIALS AND INFORMATIONS ON TRANSITIONS.            CALC-331
      N=NPOT-NT                                                         CALC-332
      M=NPOT+NW(2,NPOT+1)                                               CALC-333
      IF (N+M.GT.IDMT) CALL MEMO('CALC',IDMT,N+M)                       CALC-334
      NN=2*N                                                            CALC-335
      N1=2*NT-2                                                         CALC-336
      N2=2*M-2                                                          CALC-337
      DO 4 I=1,NN                                                       CALC-338
    4 CW(N2+I)=CW(N1+I)                                                 CALC-339
      NN=NN+2*NW(2,NPOT+1)                                              CALC-340
      N2=N2-2*NW(2,NPOT+1)                                              CALC-341
      DO 5 I=1,NN                                                       CALC-342
    5 CW(N1+I)=CW(N2+I)                                                 CALC-343
      NPOT=NT                                                           CALC-344
      NT=NPOT+NW(2,NPOT+1)                                              CALC-345
      NIVQ=NT+3*IT                                                      CALC-346
      NIVY=NIVQ+(3*IM+1)/2                                              CALC-347
      NIVZ=NIVY+(7*INTC+1)/2                                            CALC-348
    6 NCOI=NIVZ+2*INTV                                                  CALC-349
      MIPI=NCOI+2*NCONT                                                 CALC-350
      NXD=MIPI                                                          CALC-351
      MWV=NXD                                                           CALC-352
      NIE=0                                                             CALC-353
      NIXT=MWV+22*NIE                                                   CALC-354
    7 CALL DISP(NW,DW(NWV),NW(1,NIPP),DW(NIPP),NW(1,NPOT),DW(NPOT),NCOLXCALC-355
     1-NCONT,LO)                                                        CALC-356
      IF (NCONT.EQ.0) GO TO 10                                          CALC-357
      IF (.NOT.(LO(114).OR.LO(115))) GO TO 9                            CALC-358
      CALL CONU(0,NW,DW(NWV),NW(1,MIPI),DW(MWV),NW(1,NCOI),DW(NXD),DW(NSCALC-359
     1CN),KXT,LO)                                                       CALC-360
      NXD=MIPI+(11*NIE+1)/2                                             CALC-361
      MWV=NXD+3*NIE                                                     CALC-362
      NIXT=MWV+22*NIE                                                   CALC-363
      IF (NIXT.GT.IDMT) CALL MEMO('CALC',IDMT,NIXT)                     CALC-364
      NCOLT=NCOLX+NIE-NCONT                                             CALC-365
      NSP(1)=NSP(1)+NCOLT-NCOLX                                         CALC-366
      NSP(3)=NSP(3)+NCOLT-NCOLX                                         CALC-367
    8 CALL CONU(1,NW,DW(NWV),NW(1,MIPI),DW(MWV),NW(1,NCOI),DW(NXD),DW(NSCALC-368
     1CN),KXT,LO)                                                       CALC-369
    9 CALL DISP(NW(1,MIPI),DW(MWV),NW(1,NIPP),DW(NIPP),NW(1,NPOT),DW(NPOCALC-370
     1T),-NIE,LO)                                                       CALC-371
   10 NDP=2*NCOLL+NSP(1)+1                                              CALC-372
      NDQ=KBA+NSP(3)-NSP(1)                                             CALC-373
      NTY=NIXT+KXT                                                      CALC-374
      NSM=NTY                                                           CALC-375
      IF (LO(81).AND.LO(132)) NSM=NSM+KMAX*(KCC+2+NCOLT-NCOLS)          CALC-376
      NTX=NSM+2*NJMAX*KBA                                               CALC-377
      NRY=NTX+NCOLS+1                                                   CALC-378
      IF (LO(81)) NRY=NTX+NCOLL+NCOLT+3                                 CALC-379
      NRCO=NRY                                                          CALC-380
      IF (LO(81)) NRCO=NRCO+KMAX*NCOLS                                  CALC-381
      NRZ=NRCO-NSM                                                      CALC-382
      NTZ=NRCO-NTY                                                      CALC-383
      IF (LO(32)) NRCO=2*NRCO-NSM                                       CALC-384
      NRDO=NRCO+2*NPX                                                   CALC-385
      IF (LO(100)) NRDO=NRCO+2*NCOLT                                    CALC-386
      NVC1=NRDO+2*INTV                                                  CALC-387
      NVC2=NVC1+ITY(2)*ISM                                              CALC-388
      NNC=NVC1+ITX(1)*ISM                                               CALC-389
      IF (LO(100)) NNC=NVC1+ITX(7)*ISM                                  CALC-390
      NCX=NVC1+ITXM*ISM                                                 CALC-391
      NPLACE=MAX0(NPLACE,NCX,NNC+((NREC+1)*(2*NTOT+2+NREC)/2+1)/2)      CALC-392
   11 CALL GGDR(NW,DW(NWV),DW(NSCN),LO)                                 CALC-393
      IF (LO(113)) GO TO 12                                             CALC-394
      IF (LO(54)) WRITE (MW,1001) NPLACE                                CALC-395
      IF (NPLACE.GT.IDMT) CALL MEMO('CALC',IDMT,NPLACE)                 CALC-396
C MAIN COMPUTATION.                                                     CALC-397
   12 CALL CAL1(NW,CW,DW,LO)                                            CALC-398
      IF (LO(77)) CALL HORA                                             CALC-399
      IF (LO(118).OR.(.NOT.LO(32))) GO TO 18                            CALC-400
C AUTOMATIC SEARCH.                                                     CALC-401
      IF (NW(1,NIW+1).GT.1.OR.LO(76).OR.LO(75)) GO TO 14                CALC-402
C CHANGE OF CONTROLS IF FULL OUTPUT WAS NOT REQUESTED AT THE FIRST RUN. CALC-403
      DO 13 I=51,58                                                     CALC-404
      LO(I+85)=LO(I)                                                    CALC-405
   13 LO(I)=.FALSE.                                                     CALC-406
      LO(116)=.TRUE.                                                    CALC-407
C IDENTIFICATION OF VARIABLES.                                          CALC-408
   14 CALL VARI(1,NW,DW,LO)                                             CALC-409
C SAVE THE SEARCH ON TAPE MS.                                           CALC-410
      IF (LO(35).AND.NW(1,NIW+1).GE.NW(2,NIW)) CALL REST(0,NW,DW,IDMT,LOCALC-411
     1)                                                                 CALC-412
C HANDLING OF VARIABLES.                                                CALC-413
   15 CALL FITE(KE,NTOT,NREC,DW(NRES),DW(NXX),DW(NDE),DW(NRC),NW(1,NRC),CALC-414
     1NW(1,NIW),DW(NNC))                                                CALC-415
C TRANSFORMATION OF VARIABLES INTO PARAMETERS.                          CALC-416
   16 CALL VARI(0,NW,DW,LO)                                             CALC-417
      IF (KE.EQ.1) GO TO 19                                             CALC-418
      LO(111)=LO(52).OR.LO(61).OR.LO(111)                               CALC-419
      LO(112)=LO(112).OR.LO(10)                                         CALC-420
      GO TO 19                                                          CALC-421
C CONTINUATION OF A PREVIOUS SEARCH.                                    CALC-422
   17 CALL REST(1,NW,DW,IDMT,LO)                                        CALC-423
      IF (LO(77)) CALL HORA                                             CALC-424
      IF ((KE.NE.1).AND.(NW(1,NIW+1).NE.1)) GO TO 16                    CALC-425
      GO TO 15                                                          CALC-426
   18 N=IDMT-NPLACE                                                     CALC-427
      WRITE (MW,1002) NPLACE,N                                          CALC-428
      IF (LO(35).AND.LO(32).AND.(KE.EQ.0)) CALL REST(0,NW,DW,IDMT,LO)   CALC-429
      IF (.NOT.LO(34)) GO TO 1                                          CALC-430
      CALL EVAL(NW,DW,CM,LO)                                            CALC-431
      KE=0                                                              CALC-432
      IF (LO(115)) GO TO 2                                              CALC-433
   19 IF (LO(111)) GO TO 3                                              CALC-434
      IF (LO(113)) GO TO 7                                              CALC-435
      IF (LO(114)) GO TO 8                                              CALC-436
      IF (LO(112)) GO TO 11                                             CALC-437
      GO TO 12                                                          CALC-438
 1000 FORMAT ('<RED.MAT.>',F10.2,1P,D20.8,0P,F10.2,2I5)                 CALC-439
 1001 FORMAT (10X,'REQUIRED WORKING FIELD',I10)                         CALC-440
 1002 FORMAT (//' *** WORKSPACE USED IN THIS COMPUTATION',I10,'  ***',  CALC-441
     110X,I10,' MEMORIES NOT USED.')                                    CALC-442
      END                                                               CALC-443
