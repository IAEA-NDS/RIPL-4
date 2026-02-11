C 07/03/07                                                      ECIS06  CAL1-000
      SUBROUTINE CAL1(NW,CW,DW,LO)                                      CAL1-001
C THIS SUBROUTINE COMPUTES FORM FACTORS (POTE), LOOKS FOR QUANTUM       CAL1-002
C NUMBERS (QUAN), CHECKS FOR CONVERGENCE WITH RESPECT TO TOTAL SPIN,    CAL1-003
C COMPUTES HELICITY AMPLITUDES (SCHE) AND COMPUTES CHI2 (RESU).         CAL1-004
C INPUT:     NW:      WORKING AREA FOR INTEGERS.                        CAL1-005
C            CW:      WORKING AREA FOR CHARACTERS.                      CAL1-006
C            DW:      WORKING AREA FOR DOUBLE PRECISION.                CAL1-007
C                     NW,CW,DW ARE IN EQUIVALENCE BY CALL.              CAL1-008
C            LO(I):   LOGICAL CONTROLS:                                 CAL1-009
C               LO(21) =.TRUE. USUAL COUPLED EQUATIONS.                 CAL1-010
C               LO(23) =.TRUE. NO USE OF PADE AND SHIFT TO USUAL COUPLEDCAL1-011
C                              EQUATIONS WHEN THERE IS NO CONVERGENCE.  CAL1-012
C               LO(24) =.TRUE. COMPUTE THE COUPLINGS AT EACH ITERATION. CAL1-013
C               LO(25) =.TRUE. COMPLETE CALCULATION UP TO THE END - (F: CAL1-014
C                              ONE ITERATION IF TWO ARE ENOUGH).        CAL1-015
C               LO(28) =.TRUE. COMPUTATION UP TO J-CONVERGENCE-(F: STOP CAL1-016
C                              WHEN ALL COUPLING TERMS NEGLIGIBLE).     CAL1-017
C               LO(43) =.TRUE. INTERPOLATION ON TOTAL SPIN.             CAL1-018
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     CAL1-019
C               LO(51) =.TRUE. OUTPUT OF POTENTIALS.                    CAL1-020
C               LO(54) =.TRUE. OUTPUT LENGTH USED IN THE WORKING FIELD. CAL1-021
C               LO(63) =.TRUE. PENETRABILITIES WRITTEN ON FILE 63.      CAL1-022
C               LO(74) =.TRUE. OUTPUT OF TIME IN DIFFERENT STEPS.       CAL1-023
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             CAL1-024
C               LO(82) =.TRUE. OLD SIMPLIFIED COMPOUND NUCLEUS.         CAL1-025
C               LO(100)=.TRUE. DIRAC EQUATION.                          CAL1-026
C               LO(104)=.TRUE. CONVERGENCE IS OBTAINED IN THE ITERATION.CAL1-027
C               LO(107)=.TRUE. ALL THE COUPLINGS CALCULATED BEFOREHAND. CAL1-028
C               LO(108)=.TRUE. DIAGONAL COULOMB CORRECTIONS ARE NEEDED. CAL1-029
C               LO(110)=.TRUE. DERIVATIVES ARE NEEDED.                  CAL1-030
C               LO(115)=.TRUE. FIRST COMPUTATION FOR THIS ENERGY.       CAL1-031
C               LO(117)=.TRUE. FOR ALL CALCULATIONS EXCEPT THE FIRST.   CAL1-032
C               LO(118)=.TRUE. FOR LAST RESULTS.                        CAL1-033
C               LO(119)=.TRUE. RESULTS WITH THE LAST CALCULATION.       CAL1-034
C               LO(120)=.TRUE. OUTPUT AND LAST CALCULATION BEST ONE.    CAL1-035
C               LO(122)=.TRUE. IDENTICAL PARTICLES WITHOUT SPIN.        CAL1-036
C               LO(124)=.TRUE. COMPUTE TRANSMISSION COEFFICIENTS.       CAL1-037
C               LO(125)=.TRUE. USUAL COUPLED EQUATIONS.                 CAL1-038
C               LO(127)=.TRUE. COULOMB CORRECTIONS IN ASYMPTOTIC REGION.CAL1-039
C               LO(128)=.TRUE. NO COPY OF UNCOUPLED FUNCTIONS AND       CAL1-040
C                              PHASE-SHIFTS.                            CAL1-041
C               LO(131)=.TRUE. TOTAL SPIN IS TOO LARGE FOR COMPOUND     CAL1-042
C                              NUCLEUS.                                 CAL1-043
C OUTPUT:     LO:     LOGICAL CONTROLS. DEFINED HERE:                   CAL1-044
C                     LO(24), LO(107), LO(115), LO(117), LO(125),       CAL1-045
C                     LO(127), LO(128) AND LO(131).                     CAL1-046
C                                                                       CAL1-047
C FOR THE COMMON  /ADDRE/, /COUPL/, /DCONS/ AND /INTEG/ SEE CALC.       CAL1-048
C FOR THE COMMON  /CONVE/ SEE CALX.                                     CAL1-049
C FOR THE COMMON  /NOEQU/ SEE QUAN.                                     CAL1-050
C FOR THE COMMON  /POTE2/ SEE REDM.                                     CAL1-051
C                                                                       CAL1-052
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ADDRE/:                     CAL1-053
C  NJIT:      DATA FOR INTERPOLATION ON TOTAL SPIN.                     CAL1-054
C  NWV:       NON INTEGER VALUES FOR THE CHANNELS.                      CAL1-055
C  NIPP:      DISPERSION PARAMETERS.                                    CAL1-056
C  NFIS:      FISSION DATA FOR COMPOUND NUCLEUS.                        CAL1-057
C  NGAM:      GAMMA DATA FOR COMPOUND NUCLEUS.                          CAL1-058
C  NNIV:      ADDRESSES OF REDUCED NUCLEAR MATRIX ELEMENTS AND COULOMB  CAL1-059
C             INTEGRALS IN NIV(NCOLL,NCOLL,3).                          CAL1-060
C  NPOT:      OPTICAL POTENTIAL PARAMETERS.                             CAL1-061
C  NBETA:     DEFORMATION PARAMETERS.                                   CAL1-062
C  NFM:       HELICITIES AND OBSERVABLES (SEE DEPH,LECD AND OBSE).      CAL1-063
C  NTGX:      BEGINNING OF CHI2 AND NORMALISATION OF DATA.              CAL1-064
C  NDONN:     EXPERIMENTAL DATA.                                        CAL1-065
C  NNVI:      TABLE OF ADDRESSES OF COUPLING COEFFICIENTS.              CAL1-066
C  NNWI:      SAME AS NNVI FOR SYMMETRISED EQUATIONS.                   CAL1-067
C  NCC:       TABLE OF ENERGIES, L*(L+1) AND L*S  (SEE QUAN).           CAL1-068
C  MCC:       SAME AS NCC FOR SYMMETRISED EQUATIONS.                    CAL1-069
C  NXA:       TABLE OF COEFFICIENTS OF SYMMETRISATION.                  CAL1-070
C  NAM1:      COMPUTATION OF OBSERVABLES (SEE OBSE).                    CAL1-071
C  NFAC:      TABLE OF LOG OF FACTORIALS FOR GEOMETRIC COEFFICIENTS.    CAL1-072
C  NFG:       COULOMB FUNCTIONS AND FINITE INTEGRALS.                   CAL1-073
C  NXG:       COULOMB PHASES AND INFINITE INTEGRALS.                    CAL1-074
C  NRES:      FUNCTIONS FOR THE SEARCH.                                 CAL1-075
C  NT:        TABLE OF REDUCED NUCLEAR MATRIX ELEMENTS.                 CAL1-076
C  NIVQ:      TABLE OF MULTIPOLES.                                      CAL1-077
C  NIVY:      TABLE OF FORM FACTOR IDENTIFICATION IVY (FOR COMPUTATION).CAL1-078
C  NIVZ:      TABLE OF FORM FACTOR IDENTIFICATION IVZ (FOR USE).        CAL1-079
C  NCOI:      ADDRESS OF THE TABLE FOR DISCRETISATION OF CONTINUUM.     CAL1-080
C  MIPI:      ADDRESS OF "IPIM" FOR THE DISCRETISATION OF CONTINUUM.    CAL1-081
C  NXD:       ADDRESS OF WEIGHTS AND STEPS OF CONTINUUM.                CAL1-082
C  MWV:       SAME AS NWV FOR THE DISCRETISATION OF CONTINUUM.          CAL1-083
C  NIXT:      TRANSMISSION COEFFICIENTS OF UNCOUPLED STATES.            CAL1-084
C  NTY:       TEMPORARY RESULTS OF COMPOUND NUCLEUS FOR INTERPOLATION.  CAL1-085
C  NSM:       STANDARD AND HELICITY SCATTERING MATRIX ELEMENTS.         CAL1-086
C  NTX:       PARTIAL REACTION CROSS SECTIONS AND COMPOUND NUCLEUS.     CAL1-087
C  NRY:       COMPOUND NUCLEUS COEFFICIENTS.                            CAL1-088
C  NRCO:      STRENGTH OF COULOMB CENTRAL POTENTIALS FOR CORRECTIONS.   CAL1-089
C  NRDO:      STRENGTH OF COULOMB TRANSITION POTENTIALS FOR CORRECTIONS.CAL1-090
C  NVC1:      REAL POTENTIALS.                                          CAL1-091
C  NVC2:      IMAGINARY POTENTIALS.                                     CAL1-092
C  NNC:       FIRST FREE ADDRESS AFTER COMPUTATION OF POTENTIALS.       CAL1-093
C  NCX:       FIRST FREE ADDRESS FOR COMPUTATION OF POTENTIALS.         CAL1-094
C   USED:     NJIT,NWV,NIPP,NFIS,NGAM,NNIV,NPOT,NBETA,NFM,NTGX,NDONN,   CAL1-095
C             NNVI,NNWI,NCC,MCC,NXA,NAM1,NFAC,NFG,NXG,NRES,NT,NIVQ,NIVY,CAL1-096
C             NIVZ,NCOI,MIPI,NXD,MWV,NIXT,NTY,NSM,NTX,NRY,NRCO,NRDO,    CAL1-097
C             NVC1,NVC2,NNC,NCX.                                        CAL1-098
C                                                                       CAL1-099
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     CAL1-100
C  H:         STEP SIZE FOR INTEGRATION.                                CAL1-101
C  BJM:       CONVERGENCE COEFFICIENT OF IMAGINARY POTENTIAL.           CAL1-102
C  EITER:     CONVERGENCE CRITERION FOR S-MATRIX.                       CAL1-103
C  ACONV:     CONVERGENCE CRITERION FOR POTENTIAL AND FUNCTION.         CAL1-104
C  CONJ:      CONVERGENCE CRITERION FOR THE SCATTERING AMPLITUDES.      CAL1-105
C  HCONV:     ACONV*H*H CONVERGENCE CRITERION FOR SECOND MEMBERS.       CAL1-106
C   USED:     H,BJM,EITER,ACONV,CONJ.                                   CAL1-107
C                                                                       CAL1-108
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     CAL1-109
C  IQM:       MAXIMUM L-VALUE OF DEFORMATION IN ROTATIONAL MODEL.       CAL1-110
C  NBT1:      NUMBER OF PHONONS (VIBRATIONS).                           CAL1-111
C  NFA:       NUMBER OF LOGARITHMS OF FACTORIALS.                       CAL1-112
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             CAL1-113
C   USED:     IQM,NBT1,NFA,NPP.                                         CAL1-114
C                                                                       CAL1-115
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     CAL1-116
C  CHB:       PLANCK CONSTANT /(2*PI) IN MEV*FERMI.                     CAL1-117
C  CMB:       ATOMIC MASS CM DIVIDED BY H-BAR*C.                        CAL1-118
C  CCZ:       COULOMB ALPHA CONSTANT.                                   CAL1-119
C  XZ:        CONVERSION FACTOR TO MILLIBARNS.                          CAL1-120
C   USED:     CHB,CMB,CCZ,XZ.                                           CAL1-121
C                                                                       CAL1-122
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /INTEG/:                     CAL1-123
C  IDMT:      TOTAL WORKING FIELD LENGTH AS SINGLE PRECISION.           CAL1-124
C  NPLACE:    MAXIMUM WORKING FIELD USED.                               CAL1-125
C  NCOLL:     NUMBER OF COUPLED CHANNELS.                               CAL1-126
C  NJMAX:     MAXIMUM NUMBER OF J-VALUES.                               CAL1-127
C  ITERM:     MAXIMUM NUMBER OF ITERATIONS.                             CAL1-128
C  JDM:       MINIMUM NUMBER OF TOTAL SPIN.                             CAL1-129
C  JIT:       NUMBER OF RATES OF INTERPOLATION ON TOTAL SPIN.           CAL1-130
C  KMIN:      MINIMUM NUMBER OF J VALUES FOR COMPOUND NUCLEUS.          CAL1-131
C  KMAX:      MAXIMUM NUMBER OF J VALUES FOR COMPOUND NUCLEUS.          CAL1-132
C  NCOLS:     NUMBER OF CHANNELS WITH ANGULAR DISTRIBUTIONS.            CAL1-133
C  NCOLT:     NUMBER OF CHANNELS INCLUDING UNCOUPLED STATES.            CAL1-134
C  LMX:       DIFFERENCE BETWEEN NUMBER OF J AND OF COULOMB FUNCTIONS.  CAL1-135
C  LMAX1:     NUMBER OF L-VALUES FOR COULOMB FUNCTIONS.                 CAL1-136
C  NLT:       MEMORIES NEEDED FOR LEGENDRE POLYNOMIALS.                 CAL1-137
C  ISM:       NUMBER OF INTEGRATION STEPS.                              CAL1-138
C  NJC:       MAXIMUM NUMBER OF OBSERVABLES AT EQUIDISTANT ANGLES.      CAL1-139
C  JTX:       MAXIMUM NUMBER OF CALCULATED VALUES FOR A PLOT.           CAL1-140
C  KCC:       NUMBER OF INDEPENDENT AMPLITUDES WITH UNCOUPLED STATES.   CAL1-141
C  MS1:       LARGEST PARTICLE MULTIPLICITY.                            CAL1-142
C  MS2:       LARGEST TARGET MULTIPLICITY.                              CAL1-143
C  KBA:       NUMBER OF INDEPENDENT AMPLITUDES WITHOUT UNCOUPLED STATES.CAL1-144
C  KAB:       MAXIMUM NUMBER OF EQUATIONS.                              CAL1-145
C  KBC:       MAXIMUM NUMBER OF SOLUTIONS.                              CAL1-146
C  JTH:       MAXIMUM NUMBER OF ANGLES FOR A PLOT.                      CAL1-147
C  NCOLR:     NUMBER OF EXPERIMENTAL ANGULAR DISTRIBUTIONS.             CAL1-148
C  LMAX2:     NUMBER OF L VALUES FOR COULOMB PHASES.                    CAL1-149
C  ITEMM:     MAXIMUM NUMBER OF ITERATIONS AT THE BEGINNING OF A RUN.   CAL1-150
C  KXT:       NUMBER OF PENETRABILITIES FOR UNCOUPLED STATES.           CAL1-151
C  LMAX3:     EFFECTIVE NUMBER OF COULOMB FUNCTIONS.                    CAL1-152
C  NRZ:       NUMBER OF RESULTS TO SAVE FOR MINIMUM CHI2..              CAL1-153
C  NTZ:       NUMBER OF MEMORIES TO INITIALISE TO ZERO FOR SCATTERING   CAL1-154
C             MATRIX AND COMPOUND NUCLEUS RESULTS..                     CAL1-155
C  IPM:       NUMBER PF J VALUES USED FOR SCATTERING MATRIX (IN CAL1).  CAL1-156
C  IPK:       NUMBER PF J VALUES USED FOR COMPOUND NUCLEUS (IN CAL1).   CAL1-157
C  MCM(1):    MAXIMUM ANGULAR MOMENTUM FOR CENTRAL COULOMB CORRECTIONS. CAL1-158
C  MCM(2):    MAXIMUM ANGULAR MOMENTUM FOR SPIN-ORBIT COULOMB CORRECT.  CAL1-159
C  NCT(1):    NUMBER OF EQUATIONS FOR POSITIVE PARITY.                  CAL1-160
C  NCT(2):    NUMBER OF EQUATIONS FOR NEGATIVE PARITY.                  CAL1-161
C  NCT(3):    NUMBER OF SOLUTIONS FOR POSITIVE PARITY.                  CAL1-162
C  NCT(4):    NUMBER OF SOLUTIONS FOR NEGATIVE PARITY.                  CAL1-163
C  NCT(5):    NUMBER OF COMPOUND NUCLEUS COEFFICIENTS, POSITIVE PARITY. CAL1-164
C  NCT(6):    NUMBER OF COMPOUND NUCLEUS COEFFICIENTS, NEGATIVE PARITY. CAL1-165
C   DEFINED:  NPLACE,ITERM,IPM,IPK.                                     CAL1-166
C   USED:     IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCOLS,    CAL1-167
C             NCOLT,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,  CAL1-168
C             JTH,NCOLR,LMAX2,ITEMM,KXT,LMAX3,NRZ,NTZ,IPM,IPK,MCM,NCT.  CAL1-169
C                                                                       CAL1-170
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NOEQU/:                     CAL1-171
C  NCXN:      NUMBER OF SOLUTIONS NEEDED.                               CAL1-172
C  NC:        NUMBER OF EQUATIONS FOR IDENTICAL PARTICLES.              CAL1-173
C  NCIN:      NUMBER OF SOLUTIONS FOR IDENTICAL PARTICLES.              CAL1-174
C  NIN:       NUMBER OF COUPLING POTENTIALS.                            CAL1-175
C  JPI:       PARITY 0 OR 1.                                            CAL1-176
C  IPJ:       VALUE OF J+1 OR J+1/2 WHERE J IS THE CHANNEL SPIN.        CAL1-177
C  R1(2):     MAXIMUM OF SCATTERING AND COMPOUND COEFFICIENT.           CAL1-178
C  NAJ:       TWICE THE CHANNEL SPIN J.                                 CAL1-179
C   DEFINED:  NCXN,NIN,JPI,IPJ,NAJ.                                     CAL1-180
C   USED:     NCXN,NC,NCIN,NIN,R1.                                      CAL1-181
C                                                                       CAL1-182
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     CAL1-183
C  NPX:       NUMBER OF POTENTIALS TAKING INTO ACCOUNT DISPERSION.      CAL1-184
C   USED:     NPX.                                                      CAL1-185
C                                                                       CAL1-186
C***********************************************************************CAL1-187
      IMPLICIT REAL*8 (A-H,O-Z)                                         CAL1-188
      LOGICAL LO(150)                                                   CAL1-189
      CHARACTER*4 CW(2,*)                                               CAL1-190
      CHARACTER*1 JP                                                    CAL1-191
      DIMENSION NW(2,*),DW(*)                                           CAL1-192
      COMMON /ADDRE/ NIPH,NJIT,NWV,NIPP,NPAR,NPAA,NSCN,NFIS,NGAM,NNIV,NPCAL1-193
     1OT,NBETA,NFM,NTGX,NDONN,NNVI,NDE,NISE,NRC,NIW,NNWI,NCC,MCC,NXA,NAMCAL1-194
     21,NFAC,NFG,NXG,NRES,NXX,NT,NIVQ,NIVY,NIVZ,NCOI,MIPI,NXD,MWV,NIXT,NCAL1-195
     3TY,NSM,NTX,NRY,NRCO,NRDO,NVC1,NVC2,NNC,NCX                        CAL1-196
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       CAL1-197
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   CAL1-198
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            CAL1-199
      COMMON /INOUT/ MR,MW,MS                                           CAL1-200
      COMMON /INTEG/ IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCOCAL1-201
     1LS,NCOLT,NBET,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,JTCAL1-202
     2H,NCOLR,NREC,NTOT,LMAX2,KE,ITEMM,KXT,LMAX3,NRZ,NTZ,IPM,IPK,MCM(2),CAL1-203
     3NCT(6)                                                            CAL1-204
      COMMON /NCOMP/ NSP(5),NCONT,NCOJ(3),NCOLX,NDP,NDQ,ACN(20)         CAL1-205
      COMMON /NOEQU/ NCXN,NIC,NCI,NC,NCIN,NIN,JPI,IPJ,R1(2),NAJ         CAL1-206
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         CAL1-207
      LMD=2                                                             CAL1-208
      IF (LO(100).OR.LO(133)) LMD=3                                     CAL1-209
      JMIN=MOD(NW(2,1)+NW(1,2),2)                                       CAL1-210
      NFAM=NNC                                                          CAL1-211
      NPAD=NFAM+KAB*(2*KBC+10)                                          CAL1-212
      IF (LO(124)) NPAD=NFAM+KAB*(2*KAB+10)                             CAL1-213
      NR5=NPAD+2*ITEMM*(KAB+2)                                          CAL1-214
      IF (LO(21).OR.LO(22)) NR5=NPAD                                    CAL1-215
      NWW=NR5+2*KAB                                                     CAL1-216
      NXC=NWW+4*ISM*KAB                                                 CAL1-217
      IF (LO(120).AND.(.NOT.LO(51))) GO TO 2                            CAL1-218
      IF (LO(119)) GO TO 22                                             CAL1-219
      IF (LO(100)) NXC=NXC+4*ISM*KAB                                    CAL1-220
      IF (LO(21)) NXC=NR5                                               CAL1-221
      ID1=IDMT-NCX                                                      CAL1-222
C COMPUTATION OF POTENTIALS AND FORM FACTORS.                           CAL1-223
      N=NCX-1                                                           CAL1-224
      DO 1 I=NRCO,N                                                     CAL1-225
    1 DW(I)=0.D0                                                        CAL1-226
      CALL POTE(DW(NBETA),NW(1,NBETA),NW(1,NIVQ),NW(1,NIVY),NW(1,NIVZ),DCAL1-227
     1W(NPOT),NW(1,NPOT),DW(NWV),DW(MWV),ISM,NCOLL,NCOLX-NCONT,NCOLT,ID1CAL1-228
     2,NW,NW(1,MIPI),NW(1,NIPP),MCM,LO,DW(NVC1),DW(NRCO),DW(NRDO),DW(NIXCAL1-229
     3T),DW(NCX),NW(1,NCX))                                             CAL1-230
    2 IF (LO(120)) GO TO 22                                             CAL1-231
      NPLACE=MAX0(NPLACE,NXC+ID1)                                       CAL1-232
      JPI=0                                                             CAL1-233
      IP1=1                                                             CAL1-234
      IF (LO(122)) IP1=2                                                CAL1-235
      NPT=0                                                             CAL1-236
      N=NTY+NTZ-1                                                       CAL1-237
      DO 3 I=NTY,N                                                      CAL1-238
    3 DW(I)=0.D0                                                        CAL1-239
      IF ((.NOT.LO(81)).OR.LO(82)) GO TO 4                              CAL1-240
      NX=NNC+2*ISM                                                      CAL1-241
      N=NX+2*(ISM+2)                                                    CAL1-242
      IF (N.GT.IDMT) CALL MEMO('CAL1',IDMT,N)                           CAL1-243
      NPLACE=MAX0(NPLACE,N)                                             CAL1-244
C LOOP ON THE PARITIES  JPI IS THE PARITY (0 OR 1).                     CAL1-245
    4 IPJ=1                                                             CAL1-246
      AJM=-1.D0                                                         CAL1-247
      LO(127)=.FALSE.                                                   CAL1-248
      LO(131)=.NOT.LO(81)                                               CAL1-249
      ITERM=ITEMM                                                       CAL1-250
      KB=(KAB+1)/2                                                      CAL1-251
      DO 5 K=1,KB                                                       CAL1-252
      NW(1,NCC+k-1)=0                                                   CAL1-253
    5 NW(2,NCC+k-1)=0                                                   CAL1-254
C LOOP ON THE VALUES OF J.                                              CAL1-255
C SEARCH FOR QUANTUM NUMBERS AND COUPLING COEFFICIENTS.                 CAL1-256
    6 LO(125)=LO(21)                                                    CAL1-257
      NAJ=JMIN+2*IPJ-2                                                  CAL1-258
      AJ=0.5D0*DFLOAT(NAJ)                                              CAL1-259
      ITERR=0                                                           CAL1-260
    7 CALL QUAN(NCOLL,DW(NWV),NW,NW(1,NNIV),NW(1,NT),DW(NT),NW(1,NIVQ),NCAL1-261
     1W(1,NIVZ),NW(1,NCC),NW(1,NXC),DW(NXC),IU,NW(1,NNVI),KAB,KBC,DW(NXACAL1-262
     2),DW(NFAC),NFA,IDMT-NXC-200,LMD,LO)                               CAL1-263
      IF (NCIN.EQ.0) GO TO 15                                           CAL1-264
      AJM=AJ                                                            CAL1-265
      NX=NXC+IU*LMD                                                     CAL1-266
      IPX=IPJ                                                           CAL1-267
      LO(131)=LO(131).OR.IPJ.GT.KMAX                                    CAL1-268
      NX1=NX                                                            CAL1-269
      NCXN=NCIN                                                         CAL1-270
      IF (LO(124)) NCXN=NC                                              CAL1-271
      NWR=2*NX1+20*KAB*KAB-1                                            CAL1-272
      IF (NWR+10*KAB**2.GT.IDMT) CALL MEMO('CAL1',IDMT,NWR+10*KAB**2)   CAL1-273
      CALL MTCH(NW(1,NNWI),NCOLL,KAB,DW(NWV),NW(1,MCC),DW(NXG),NW(1,NXC)CAL1-274
     1,DW(NXC),DW(NX1),ISM,LMAX2,NW(1,NNIV),NW(1,NIVZ),DW(NFG),LMAX1,LMACAL1-275
     2X3,NW(1,NWR),DW(NRCO),DW(NRDO),DW(NFAM),LMD,LO)                   CAL1-276
      IF (LO(127)) GO TO 12                                             CAL1-277
      IF (LO(44)) NX=NX+4*KAB*KAB                                       CAL1-278
      IF (LO(74)) CALL HORA                                             CAL1-279
C NX FIRST ADDRESS WHICH CAN BE USED.                                   CAL1-280
      IF (LO(125)) GO TO 12                                             CAL1-281
C ITERATIONS.                                                           CAL1-282
      NFAR=NFAM+10*KAB                                                  CAL1-283
      NFAI=NFAR+NCXN*KAB                                                CAL1-284
      IF (LO(100)) GO TO 8                                              CAL1-285
      NWR=2*NX+4*(ISM+2)-1                                              CAL1-286
      NR4=NWR+2*NC*ISM                                                  CAL1-287
      IF (LO(110)) NR4=NWR+4*NC*ISM                                     CAL1-288
      N=NR4+2*ISM                                                       CAL1-289
      GO TO 9                                                           CAL1-290
    8 NWR=NX+2*(ISM+2)                                                  CAL1-291
      NR4=NWR+4*NC*ISM                                                  CAL1-292
      NIN=4*NIN                                                         CAL1-293
C  N IS THE FIRST FREE ADDRESS AND NIN THE NUMBER OF COUPLING POTENTIALSCAL1-294
C WHICH CAN BE STORED.                                                  CAL1-295
      N=NR4+4*ISM                                                       CAL1-296
    9 N3=N                                                              CAL1-297
      IF (N.GT.IDMT) CALL MEMO('CAL1',IDMT,N)                           CAL1-298
      LO(107)=(.NOT.LO(24)).AND.ITERM.GT.1                              CAL1-299
      IF (BJM.NE.0.D0) NIN=NIN+NC                                       CAL1-300
      IF (LO(107)) N3=N3+NIN*ISM                                        CAL1-301
      IF (LO(54)) WRITE (MW,1000) N3,AJM                                CAL1-302
      IF (N3.LE.IDMT) GO TO 10                                          CAL1-303
      LO(107)=.FALSE.                                                   CAL1-304
      LO(24)=.TRUE.                                                     CAL1-305
      N3=N                                                              CAL1-306
      WRITE (MW,1001)                                                   CAL1-307
   10 NPLACE=MAX0(NPLACE,N3)                                            CAL1-308
      IF (LO(100)) GO TO 11                                             CAL1-309
      CALL INTI(DW(NFAM),DW(NX),DW(NWW),DW(NWR),DW(NR4),DW(NPAD),NW(1,NRCAL1-310
     15),ISM,KAB,DW(N),ITERM,NC,DW(NVC1),NW(1,NNWI),NW(1,MCC),NW(1,NXC),CAL1-311
     2DW(NXC),LMD,DW(NX1),NCXN,NNI,ITERR,LO)                            CAL1-312
      IF (LO(104).OR.ITERM.EQ.1) GO TO 13                               CAL1-313
      IF (.NOT.LO(23)) GO TO 13                                         CAL1-314
      IF (LO(110)) WRITE (MW,1002) AJM,JPI                              CAL1-315
      LO(125)=.TRUE.                                                    CAL1-316
      LO(128)=.TRUE.                                                    CAL1-317
      GO TO 7                                                           CAL1-318
   11 CALL INTR(DW(NFAM),DW(NX),DW(NWW),DW(NWR),DW(NR4),DW(NPAD),NW(1,NRCAL1-319
     15),ISM,KAB,DW(N),ITERM,NC,DW(NVC1),DW(NVC2),NW(1,NNWI),NW(1,MCC),NCAL1-320
     2W(1,NXC),DW(NXC),DW(NX1),NCXN,NNI,ITERR,LO)                       CAL1-321
      GO TO 13                                                          CAL1-322
C USUAL COUPLED CHANNELS CALCULATION.                                   CAL1-323
C NWR IS THE FIRST ADDRESS TO STORE POTENTIALS AND NM THE MAXIMUM       CAL1-324
C NUMBER OF POINTS   FOR ARGUMENTS  SEE INCH.                           CAL1-325
   12 NWR=NX+6*KAB*KAB                                                  CAL1-326
      IF (LO(54)) WRITE (MW,1000) NWR,AJM                               CAL1-327
      NFAI=NX+KAB*KAB                                                   CAL1-328
      NFAR=NFAI+3*KAB*KAB                                               CAL1-329
      IF (LO(127)) GO TO 13                                             CAL1-330
      N=NWR+2*NC*NC*(ISM+1)                                             CAL1-331
      IF (N.GT.IDMT) CALL MEMO('CAL1',IDMT,N)                           CAL1-332
      NM=(IDMT-NWR)/(2*NC*NC)-1                                         CAL1-333
      NM=MIN0(NM,ISM)                                                   CAL1-334
      IF (LO(54).AND.NM.NE.ISM) WRITE (MW,1003) NM                      CAL1-335
      N3=NWR+2*(NM+1)*NC*NC                                             CAL1-336
      NPLACE=MAX0(NPLACE,N3,NFAR+3*KAB*KAB)                             CAL1-337
      CALL INCH(DW(NVC1),NW(1,MCC),NW(1,NXC),DW(NXC),LMD,NW(1,NNWI),DW(NCAL1-338
     1FAM),DW(NX1),DW(NX),DW(NWR),ISM,KAB,NC,NCXN,NM,ITERM,NW(1,NWR),LO)CAL1-339
   13 IF (LO(74)) CALL HORA                                             CAL1-340
C NFAR AND NFAI  ADDRESSES OF REAL AND IMAGINARY PARTS OF S-MATRIX.     CAL1-341
      NMX=NWR+(4+NCIN)*(NC+KXT)                                         CAL1-342
      N=NMX+MAX0(4*NC*NC,LMAX2+4*IPJ)                                   CAL1-343
      IF (N.GT.IDMT) CALL MEMO('CAL1',IDMT,N)                           CAL1-344
      NPT=NPT+1                                                         CAL1-345
      CALL SCAM(DW(NSM),DW(NTY),DW(NTX),DW(NRY),NJMAX,KMAX,NW(1,MCC),NW(CAL1-346
     11,NCC),DW(NFAR),DW(NFAI),DW(NWV),NCOLL,NCOLS,KAB,KCC,NW,NW(1,MIPI)CAL1-347
     2,DW(NGAM),DW(NFIS),DW(NIXT),NW(1,NCOI),DW(NXD),DW(NMX),DW(NMX),DW(CAL1-348
     3NWR),NCT,NCIN+4,DW(NXA),KBC,IDMT-NMX,LO)                          CAL1-349
      IF (LO(131)) GO TO 14                                             CAL1-350
      LO(131)=R1(2).LT.CONJ.AND.IPJ.GT.KMIN                             CAL1-351
      IF (.NOT.LO(131)) IPK=IPJ                                         CAL1-352
      IF (LO(131).AND.(LO(118).OR.LO(115))) WRITE (MW,1004) AJM,R1(2)   CAL1-353
   14 IF (LO(74)) CALL HORA                                             CAL1-354
      IF (LO(21)) GO TO 15                                              CAL1-355
C REDUCTION OF MAXIMUM NUMBER OF ITERATIONS IF 2 WHERE SUFFICIENT(ITERR)CAL1-356
      IF ((.NOT.LO(25)).AND.(ITERR.LE.2).AND.(IPJ.GT.JDM)) ITERM=1      CAL1-357
C CHECKS OF CONVERGENCE.                                                CAL1-358
      IF ((NNI.EQ.NC).AND.(IPJ.GT.JDM)) LO(127)=.TRUE.                  CAL1-359
      IF ((.NOT.LO(28)).AND.NNI.EQ.NC.AND.IPJ.GT.JDM+1) GO TO 18        CAL1-360
C INCREASE OF THE TOTAL ANGULAR MOMENTUM.                               CAL1-361
   15 IF (.NOT.LO(43)) GO TO 17                                         CAL1-362
      DO 16 I=1,JIT                                                     CAL1-363
      IF (IPJ.LE.NW(1,NJIT+I-1)) GO TO 17                               CAL1-364
   16 IPJ=IPJ+NW(2,NJIT+I-1)*IP1                                        CAL1-365
   17 IPJ=IPJ+IP1                                                       CAL1-366
      IF (IPJ.GT.LMAX3-LMX) LO(127)=.TRUE.                              CAL1-367
      IF (.NOT.LO(108)) LO(127)=.FALSE.                                 CAL1-368
      IF ((IPJ.LE.NJMAX).AND.((R1(1).GE.CONJ).OR.(IPJ.LE.JDM).OR.(NCIN.ECAL1-369
     1Q.0))) GO TO 6                                                    CAL1-370
C CHANGE OF PARITY.                                                     CAL1-371
   18 NAJ=JMIN+2*IPX-2                                                  CAL1-372
      IF (LO(118).OR.LO(115)) WRITE (MW,1005) AJM,R1(1)                 CAL1-373
      IF (JPI.EQ.1) GO TO 19                                            CAL1-374
      JPI=JPI+1                                                         CAL1-375
      IPY=IPX                                                           CAL1-376
      IPZ=IPK                                                           CAL1-377
      GO TO 4                                                           CAL1-378
   19 IPM=MAX0(IPX,IPY)                                                 CAL1-379
      IPK=MAX0(IPK,IPZ)                                                 CAL1-380
      IF (.NOT.LO(63)) GO TO 22                                         CAL1-381
      WRITE (63,1006) DW(NWV),DW(NWV+12),DW(NWV+1),NW(2,2),NPT          CAL1-382
      REWIND 99                                                         CAL1-383
      DO 21 I=1,NPT                                                     CAL1-384
      READ (99,1007) U1,JP,K3                                           CAL1-385
      WRITE (63,1007) U1,JP,K3                                          CAL1-386
      DO 20 K=1,K3                                                      CAL1-387
      READ (99,1008) K1,K2,U1,U2                                        CAL1-388
   20 WRITE (63,1008) K1,K2,U1,U2                                       CAL1-389
   21 CONTINUE                                                          CAL1-390
      CLOSE (99,STATUS='DELETE')                                        CAL1-391
   22 KCB=MAX0(NCT(5),NCT(6))                                           CAL1-392
      NDX=NXC+4*KCB                                                     CAL1-393
      NDY=NDX+JTH                                                       CAL1-394
      NDZ=NDY+JTH                                                       CAL1-395
      NMY=NDZ+JTX                                                       CAL1-396
      NGA=NMY+3*NJC                                                     CAL1-397
      NMA=NGA+4*KBC*KAB                                                 CAL1-398
      NMC=NMA+MS1*MS2*KCB                                               CAL1-399
      N=NMC+2*IPM                                                       CAL1-400
      IF (N.GT.IDMT) CALL MEMO('CAL1',IDMT,N)                           CAL1-401
      CALL SCHE(DW(NSM),NJMAX,KMAX,NW,MS1,MS2,NW(1,NFM),DW(NTX),DW(NMC),CAL1-402
     1NW(1,NXC),DW(NGA),DW(NXG),LMAX2,DW(NWV),KAB,KBA,KCB,JMIN,IPM,IPK,DCAL1-403
     2W(NTY),NCOLL,NCOLS,NCT,DW(NRY),DW(NMA),JIT,NW(1,NJIT),NLT,IDMT-NMCCAL1-404
     3,LO)                                                              CAL1-405
      ID1=IDMT-NGA                                                      CAL1-406
C FOR ARGUMENTS    SEE CALX AND RESU.                                   CAL1-407
      CALL RESU(NW,DW(NSM),DW(NTX),DW(NSM),NJMAX,KMAX,NCOLL,NCOLS,NW(1,NCAL1-408
     1FM),CW(1,NFM),NW(1,NTGX),DW(NTGX),IPM,IPK,DW(NDONN),NCOLR,NW(1,NAMCAL1-409
     21),DW(NAM1),DW(NWV),DW(MWV),DW(NRY),NW(1,NCOI),DW(NXD),JMIN,NRZ,NJCAL1-410
     3C,DW(NRES),DW(NGA),DW(NMY),DW(NDX),DW(NDY),DW(NDZ),ID1,LO)        CAL1-411
      NPLACE=MAX0(NPLACE,NGA+ID1)                                       CAL1-412
      IF (LO(74)) CALL HORA                                             CAL1-413
      LO(117)=.TRUE.                                                    CAL1-414
      LO(115)=.FALSE.                                                   CAL1-415
      RETURN                                                            CAL1-416
 1000 FORMAT (10X,'REQUIRED WORKING FIELD',I10,'  FOR J =',F7.1)        CAL1-417
 1001 FORMAT (/' WORKING FIELD TOO SMALL TO STORE ALL THE POTENTIALS'/1XCAL1-418
     1,52('*')//' THE 24TH CONTROL IS SET .TRUE.'//)                    CAL1-419
 1002 FORMAT (' WARNING: FOR J =',F6.1,' PARITY (-1)**',I1,' THE DERIVATCAL1-420
     1IVE TERMS ARE NEGLECTED.')                                        CAL1-421
 1003 FORMAT ('+',60X,'COMPUTATION OF POTENTIALS BY',I5,'  AT A TIME')  CAL1-422
 1004 FORMAT (' MAXIMUM J-VALUE =',F6.1,' FOR COMPOUND NUCLEUS',6X,'MAXICAL1-423
     1MUM COEFFICIENT AT THE END',D12.3)                                CAL1-424
 1005 FORMAT (' MAXIMUM J-VALUE =',F6.1,16X,'MAXIMUM SCATTERING COEFFICICAL1-425
     1ENT AT THE END',D12.3)                                            CAL1-426
 1006 FORMAT ('<TLJ     >',F10.2,1P,D20.8,0P,F10.2,2I5)                 CAL1-427
 1007 FORMAT (1X,F9.1,4X,A1,1X,I4)                                      CAL1-428
 1008 FORMAT (1X,I2,I6,F9.1,2X,1P,D18.8,0P)                             CAL1-429
      END                                                               CAL1-430
