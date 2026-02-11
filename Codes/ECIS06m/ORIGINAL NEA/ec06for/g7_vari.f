C 29/10/06                                                      ECIS06  VARI-000
      SUBROUTINE VARI(KF,NW,DW,LO)                                      VARI-001
C ASSIGNS THE VALUES OF THE VARIABLE PARAMETERS AND PRINTS THE OUTPUT   VARI-002
C DURING THE SEARCH.                                                    VARI-003
C INPUT:     KF:      1 BEFORE THE CALL TO THE SEARCH SUBROUTINES:      VARI-004
C                     AT THE FIRST CALL(KE=0) DEFINES THE VARIABLES     VARI-005
C                     IN SEARCH,                                        VARI-006
C                     AT ANY CALL, PRINT CHI2 AND VARIABLES.            VARI-007
C                     0 AFTER THE CALL TO THE SEARCH SUBROUTINES:       VARI-008
C                     AT ANY CALL SETS PARAMETERS TO THEIR VALUE FOR    VARI-009
C                     NEXT EVALUATION,                                  VARI-010
C                     AT THE LAST CALL AND IF KE IS NOT 1, PRINT        VARI-011
C                     ERRORS AND PARAMETERS.                            VARI-012
C            NW:      WORKING AREA FOR INTEGERS.                        VARI-013
C            DW:      WORKING AREA FOR DOUBLE PRECISION IN EQUIVALENCE  VARI-014
C                     BY CALL WITH NW.                                  VARI-015
C            LO:      LOGICAL CONTROLS:                                 VARI-016
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).VARI-017
C               LO(3)  =.TRUE. ANHARMONIC VIBRATIONAL OR ASYMMETRIC     VARI-018
C                              ROTATIONAL MODEL.                        VARI-019
C               LO(4)  =.TRUE. PARAMETRISED SPIN-ORBIT DEFORMATION.     VARI-020
C               LO(6)  =.TRUE. USE DEFORMATION LENGTHS.                 VARI-021
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    VARI-022
C               LO(16) =.TRUE. HEAVY-ION DEFINITION OF REDUCED RADII ANDVARI-023
C                              DEFORMATIONS.                            VARI-024
C               LO(17) =.TRUE. FOLDING MODEL.                           VARI-025
C               LO(41) =.TRUE. FACTORISATION OF 1/(1-COS(THETA)).       VARI-026
C               LO(53) =.TRUE. OUTPUT OF THE NUMBER OF ITERATIONS.      VARI-027
C               LO(55) =.TRUE. OUTPUT OF C-MATRIX ELEMENTS AND OF       VARI-028
C                              COMPOUND NUCLEUS INTERMEDIATE RESULTS.   VARI-029
C               LO(56) =.TRUE. OUTPUT OF S-MATRIX ELEMENTS.             VARI-030
C               LO(57) =.TRUE. OUTPUT PHASE-SHIFTS AT EACH ITERATION.   VARI-031
C               LO(58) =.TRUE. OUTPUT OF THE COEFFICIENTS OF EACH FORM  VARI-032
C                              FACTOR FOR ALL SETS OF EQUATIONS.        VARI-033
C               LO(60) =.TRUE. S-MATRIX ELEMENTS WRITTEN ON FILE 60.    VARI-034
C               LO(62) =.TRUE. POTENTIALS WRITTEN ON FILE 62.           VARI-035
C               LO(64) =.TRUE. PRINT RESULTS ON FILES 64 AND 66.        VARI-036
C               LO(76) =.TRUE. LO(51) TO LO(65) ARE ALWAYS USED.        VARI-037
C               LO(78) =.TRUE. OUTPUT OF DIFFERENCES BETWEEN            VARI-038
C                              EXPERIMENTAL AND CALCULATED VALUES.      VARI-039
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             VARI-040
C               LO(82) =.TRUE. OLD SIMPLIFIED COMPOUND NUCLEUS.         VARI-041
C               LO(87) =.TRUE. NO WIDTH FLUCTUATIONS.                   VARI-042
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      VARI-043
C               LO(111)=.TRUE. NUCLEAR PARAMETERS ARE CHANGED IN SEARCH.VARI-044
C               LO(112)=.TRUE. SPIN-ORBIT OR COMPOUND NUCLEUS PARAMETERSVARI-045
C                              ARE CHANGED IN SEARCH.                   VARI-046
C               LO(113)=.TRUE. DISPERSION RELATION IS CHANGED IN SEARCH.VARI-047
C               LO(114)=.TRUE. COMPOUND NUCLEUS CONTINUUM IS CHANGED IN VARI-048
C                              SEARCH.                                  VARI-049
C               LO(116)=.TRUE. NO OUTPUT.                               VARI-050
C               LO(118)=.TRUE. FOR LAST RESULTS.                        VARI-051
C               LO(119)=.TRUE. RESULTS WITH THE LAST CALCULATION.       VARI-052
C               LO(120)=.TRUE. OUTPUT AND LAST CALCULATION BEST ONE.    VARI-053
C                                                                       VARI-054
C FOR THE COMMON  /ADDRE/, /COUPL/ AND /INTEG/ SEE CALC.                VARI-055
C FOR THE COMMON  /DCHI2/ AND /NCOMP/ SEE CALX.                         VARI-056
C                                                                       VARI-057
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ADDRE/:                     VARI-058
C  NWV:       NON INTEGER VALUES FOR THE CHANNELS.                      VARI-059
C  NIPP:      DISPERSION PARAMETERS.                                    VARI-060
C  NPAA:      VALUES OF NUCLEAR PARAMETERS.                             VARI-061
C  NSCN:      LEVEL DENSITY DESCRIPTION.                                VARI-062
C  NFIS:      FISSION DATA FOR COMPOUND NUCLEUS.                        VARI-063
C  NGAM:      GAMMA DATA FOR COMPOUND NUCLEUS.                          VARI-064
C  NPOT:      OPTICAL POTENTIAL PARAMETERS.                             VARI-065
C  NBETA:     DEFORMATION PARAMETERS.                                   VARI-066
C  NTGX:      BEGINNING OF CHI2 AND NORMALISATION OF DATA.              VARI-067
C  NISE:      INDEXES OF THE VARIABLE PARAMETERS IN SEARCH.             VARI-068
C  NRC:       PERMANENT WORKING FIELD FOR THE SEARCH.                   VARI-069
C  NIW:       INTEGER WORKING FIELD FOR THE SEARCH.                     VARI-070
C  NRES:      FUNCTIONS FOR THE SEARCH.                                 VARI-071
C  NXX:       VARIABLES FOR THE SEARCH.                                 VARI-072
C  NT:        TABLE OF REDUCED NUCLEAR MATRIX ELEMENTS.                 VARI-073
C  NIVQ:      TABLE OF MULTIPOLES.                                      VARI-074
C  NIVY:      TABLE OF FORM FACTOR IDENTIFICATION IVY (FOR COMPUTATION).VARI-075
C   USED:     NWV,NIPP,NPAA,NSCN,NFIS,NGAM,NPOT,NBETA,NTGX,NISE,NRC,NIW,VARI-076
C             NRES,NXX,NT,NIVQ,NIVY.                                    VARI-077
C                                                                       VARI-078
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     VARI-079
C  IQM:       MAXIMUM L-VALUE OF DEFORMATION IN ROTATIONAL MODEL.       VARI-080
C  IQMAX:     MAXIMUM L-EXPANSION IN ROTATIONAL MODEL.                  VARI-081
C  NBT1:      NUMBER OF PHONONS (VIBRATIONS).                           VARI-082
C  NFA:       NUMBER OF LOGARITHMS OF FACTORIALS.                       VARI-083
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             VARI-084
C  NSPIN:     TWICE THE K-VALUE OF THE ROTATIONAL BAND.                 VARI-085
C  NVA:       NUMBER OF NUCLEAR PARAMETERS.                             VARI-086
C   USED:     NPP,NVA.                                                  VARI-087
C                                                                       VARI-088
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCHI2/:                     VARI-089
C  CHI2:      CHI-SQUARE COMPUTED IN SUBROUTINE RESU.                   VARI-090
C  CHI2M:     MINIMUM CHI-SQUARE IN THE SEARCH.                         VARI-091
C  YY(1):     STEP SIZE IN THE SEARCH.                                  VARI-092
C  YY(2):     HALF OF THE SUCCESS MULTIPLICATIVE FACTOR OF THE STEP.    VARI-093
C  YY(3):     VARIOUS MEANINGS.  SEE FITE.                              VARI-094
C   USED:     CHI2,CHI2M,YY.                                            VARI-095
C                                                                       VARI-096
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /INTEG/:                     VARI-097
C  NBET:      NUMBER OF DIFFERENT DEFORMATIONS (VIBRATIONS+ROTATIONS).  VARI-098
C  NCOLR:     NUMBER OF EXPERIMENTAL ANGULAR DISTRIBUTIONS.             VARI-099
C  NREC:      NUMBER OF VARIABLES IN SEARCH.                            VARI-100
C  NTOT:      NUMBER OF EXPERIMENTAL DATA.                              VARI-101
C  KE:        CONTROL OF SEARCH (SEE FITE).                             VARI-102
C   USED:     NBET,NCOLR,NREC,NTOT,KE.                                  VARI-103
C                                                                       VARI-104
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     VARI-105
C  NFISS:     NUMBER OF FISSION TRANSMISSION COEFFICIENTS.              VARI-106
C  NRD:       NUMBER OF GAMMA TRANSMISSION COEFFICIENTS.                VARI-107
C  NCONT:     NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS.                 VARI-108
C  NCONS:     NUMBER OF LEVEL DENSITIES NEEDED.                         VARI-109
C  AZ(I):     SPIN-ORBIT PARAMETRISATION FOR I-1 TO 6,                  VARI-110
C             HAUSER FESHBACH PARAMETERS FOR J=7 TO 11,                 VARI-111
C             GIANT DIPOLE RESONANCE FOR I=12 TO 16.                    VARI-112
C   DEFINED:  AZ.                                                       VARI-113
C   USED:     NFISS,NRD,NCONT,NCONS,AZ.                                 VARI-114
C                                                                       VARI-115
C PARAMETERS IN SEARCH ARE GIVEN BY NW(2-MOD(I,2),NISE+I/2)=INDEX(I)    VARI-116
C STARTING AT NW(1,NISE) FOR I=1 TO NREC:                               VARI-117
C   POSITIVE VALUE: PARAMETER DEFINED BY INDEX(I)                       VARI-118
C   NEGATIVE VALUE -L: PARAMETERS DEFINED BY INDEX(J) FOR J=L+1 TO      VARI-119
C     J=L+INDEX(L) ARE DEFINED AS THE SAME VARIABLE.                    VARI-120
C      1-1000      OPTICAL MODEL,FOLDING PARAMETERS.                    VARI-121
C   1001-2000      DEFORMATIONS FOR A GIVEN POTENTIAL:  - LX(1)=.TRUE.  VARI-122
C   2001-3000      DEFORMATIONS FOR A GIVEN MULTIPOLE:  - LX(1)=.TRUE.  VARI-123
C   3001-4000      INDIVIDUAL DEFORMATION:              - LX(1)=.TRUE.  VARI-124
C   4001-5000      NUCLEAR MODEL PARAMETER:             - LO(111)=.TRUE.VARI-125
C   5001-6000      NUCLEAR MATRIX ELEMENT:              - LX(2)=.TRUE.  VARI-126
C   6001-7000      SPIN-ORBIT AND COMPOUND NUCLEUS PARAMETRISATION:     VARI-127
C    6001-6006     SPIN-ORBIT PARAMETRISATION:          - LX(3)=.TRUE.  VARI-128
C    6007-6011     BZ1, BZ2, BZ3, BZ4, BZ5:     - LO(112)=LX(4)=.TRUE.  VARI-129
C    6012-6016     TGO, BN, FNUG, EGD, GGD:     - LO(112)=LX(5)=.TRUE.  VARI-130
C    6017-6...     SEE DESCRIPTION OF INPUT.                            VARI-131
C   FOR GAMMA      SA, UX, TAU, SG, E0, EX:     - LO(112)=LX(6)=.TRUE.  VARI-132
C   FOR CONTINUUM  SA, UX, TAU, SG, E0, EX:     - LO(114)=LX(7)=.TRUE.  VARI-133
C                  GAMMA TRANSMISSION FACTORS:          - LX(8)=.TRUE.  VARI-134
C                  FISSION TRANSMISSION COEFFICIENT,                    VARI-135
C                  DEGREE OF FREEDOM:                  - LX(9)=.TRUE.   VARI-136
C   7001-8000      DISPERSION RELATIONS PARAMETRISATION:- LO(113)=.TRUE.VARI-137
C  10001-99999     EXTERNAL OPTICAL MODEL (PARAMETERS ABOVE 1000)       VARI-138
C IT STOPS THE CALCULATION FOR AN INDEX OF PARAMETER OUT OF RANGE.      VARI-139
C***********************************************************************VARI-140
      IMPLICIT REAL*8 (A-H,O-Z)                                         VARI-141
      LOGICAL LO(150),LX(9)                                             VARI-142
      DIMENSION NW(2,*),DW(*),R0(8)                                     VARI-143
      COMMON /ADDRE/ NIPH,NJIT,NWV,NIPP,NPAR,NPAA,NSCN,NFIS,NGAM,NNIV,NPVARI-144
     1OT,NBETA,NFM,NTGX,NDONN,NNVI,NDE,NISE,NRC,NIW,NNWI,NCC,MCC,NXA,NAMVARI-145
     21,NFAC,NFG,NXG,NRES,NXX,NT,NIVQ,NIVY,NIVZ,NCOI,MIPI,NXD,MWV,NIXT,NVARI-146
     3TY,NSM,NTX,NRY,NRCO,NRDO,NVC1,NVC2,NNC,NCX                        VARI-147
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   VARI-148
      COMMON /DCHI2/ CHI2,CHI2M,YY(3)                                   VARI-149
      COMMON /INOUT/ MR,MW,MS                                           VARI-150
      COMMON /INTEG/ IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCOVARI-151
     1LS,NCOLT,NBET,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,JTVARI-152
     2H,NCOLR,NREC,NTOT,LMAX2,KE,ITEMM,KXT,LMAX3,NRZ,NTZ,IPM,IPK,MCM(2),VARI-153
     3NCT(6)                                                            VARI-154
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQVARI-155
     1,ACN1,ACN2,AZ(18)                                                 VARI-156
C TRANSFER BETWEEN VARIABLES AND PARAMETERS                             VARI-157
      IF (KF.EQ.1.AND.KE.NE.0) GO TO 36                                 VARI-158
      DO 1 I=111,114                                                    VARI-159
    1 LO(I)=.FALSE.                                                     VARI-160
      DO 2 I=1,9                                                        VARI-161
    2 LX(I)=.FALSE.                                                     VARI-162
      DO 35 I=1,NREC                                                    VARI-163
      J=NW(2-MOD(I,2),NISE+(I-1)/2)                                     VARI-164
      JJ=IABS(J)                                                        VARI-165
      IF (J.NE.JJ) GO TO 3                                              VARI-166
      J1=I                                                              VARI-167
      J2=I                                                              VARI-168
      GO TO 4                                                           VARI-169
    3 J1=JJ+1                                                           VARI-170
      J2=J1+NW(2-MOD(JJ,2),NISE+(JJ-1)/2)-1                             VARI-171
    4 DO 34 K=J1,J2                                                     VARI-172
      J=NW(2-MOD(K,2),NISE+(K-1)/2)                                     VARI-173
      IF ((J.GT.1000).AND.(J.LE.10000)) GO TO 10                        VARI-174
C OPTICAL MODEL AND FOLDING PARAMETERS                                  VARI-175
      IF (J.GT.10000) J=J-9000                                          VARI-176
      IF (LO(7)) GO TO 5                                                VARI-177
      NVAT=42*NPP                                                       VARI-178
      IF (J.GT.NVAT) GO TO 67                                           VARI-179
      IF (MOD(J,42).EQ.25) GO TO 68                                     VARI-180
      GO TO 8                                                           VARI-181
C EXTERNAL OPTICAL PARAMETERS                                           VARI-182
    5 IF (J.GE.NW(2,NPOT+1)) GO TO 69                                   VARI-183
      I1=NW(1,NPOT)-2                                                   VARI-184
      DO 6 L=1,I1                                                       VARI-185
      IF (J.GE.NW(1,NPOT+L).AND.J.LE.NW(2,NPOT+L)) GO TO 7              VARI-186
    6 CONTINUE                                                          VARI-187
      GO TO 70                                                          VARI-188
    7 IF (L.EQ.1) GO TO 8                                               VARI-189
C NOT A FOLDING PARAMETER.                                              VARI-190
      M=NW(1,NPOT+L)                                                    VARI-191
      N=NW(2,NPOT+M-2)                                                  VARI-192
      IV=NW(2,NPOT+N-1)                                                 VARI-193
      IF (IV.GE.1.AND.IV.LE.6) GO TO 8                                  VARI-194
      IF ((IV.LT.0).AND.(J-M.GT.1)) GO TO 71                            VARI-195
      IF ((IV.NE.7).AND.(IV.NE.8)) GO TO 8                              VARI-196
      JY=J                                                              VARI-197
      IF (NW(1,NPOT+NW(2,NPOT+M-2)).NE.1) JY=JY-1                       VARI-198
      IF (M.GT.JY) GO TO 72                                             VARI-199
      IF ((IV.EQ.7).AND.(MOD(JY-M,3).NE.0)) GO TO 73                    VARI-200
      IF ((IV.EQ.8).AND.(MOD(JY+10-M,11).LT.3)) GO TO 74                VARI-201
    8 IF (KF.NE.0) GO TO 9                                              VARI-202
      IF (K.EQ.J1) RAP=DW(NXX+I-1)/DW(NPOT+J-1)                         VARI-203
      DW(NPOT+J-1)=RAP*DW(NPOT+J-1)                                     VARI-204
      GO TO 34                                                          VARI-205
    9 DW(NXX+I-1)=DW(NPOT+J-1)                                          VARI-206
      GO TO 35                                                          VARI-207
   10 J=J-1000                                                          VARI-208
      IF (J.GT.1000) GO TO 15                                           VARI-209
      IF (J.GT.8) GO TO 75                                              VARI-210
C  DEFORMATIONS FOR A GIVEN POTENTIAL                                   VARI-211
      LX(1)=.TRUE.                                                      VARI-212
      DO 11 K1=1,NBET                                                   VARI-213
      IF (LO(1).AND.LO(3).AND.NW(2,NBETA+9*K1-1).NE.0) GO TO 11         VARI-214
      IF (DW(NBETA+9*K1+J-10).NE.0.D0) GO TO 12                         VARI-215
   11 CONTINUE                                                          VARI-216
      GO TO 76                                                          VARI-217
   12 IF (KF.EQ.0) GO TO 13                                             VARI-218
      DW(NXX+I-1)=DW(NBETA+9*K1+J-10)                                   VARI-219
      GO TO 35                                                          VARI-220
   13 IF (K.EQ.J1) RAP=DW(NXX+I-1)/DW(NBETA+9*K1+J-10)                  VARI-221
      DO 14 L=K1,NBET                                                   VARI-222
      IF (LO(1).AND.LO(3).AND.NW(2,NBETA+9*K1-1).NE.0) GO TO 14         VARI-223
      DW(NBETA+9*L+J-10)=DW(NBETA+9*L+J-10)*RAP                         VARI-224
   14 CONTINUE                                                          VARI-225
      GO TO 34                                                          VARI-226
   15 J=J-1000                                                          VARI-227
      IF (J.GT.1000) GO TO 20                                           VARI-228
C DEFORMATIONS FOR A GIVEN MULTIPOLE                                    VARI-229
      IF (J.GT.NBET) GO TO 77                                           VARI-230
      LX(1)=.TRUE.                                                      VARI-231
      DO 16 K1=1,8                                                      VARI-232
      IF (DW(NBETA+9*J+K1-10).NE.0.D0) GO TO 17                         VARI-233
   16 CONTINUE                                                          VARI-234
      GO TO 78                                                          VARI-235
   17 IF (KF.EQ.0) GO TO 18                                             VARI-236
      DW(NXX+I-1)=DW(NBETA+9*J+K1-10)                                   VARI-237
      GO TO 35                                                          VARI-238
   18 IF (K.EQ.J1) RAP=DW(NXX+I-1)/DW(NBETA+9*J+K1-10)                  VARI-239
      DO 19 L=K1,8                                                      VARI-240
   19 DW(NBETA+9*J+L-10)=DW(NBETA+9*J+L-10)*RAP                         VARI-241
      GO TO 34                                                          VARI-242
   20 J=J-1000                                                          VARI-243
      IF (J.GT.1000) GO TO 21                                           VARI-244
C INDIVIDUAL DEFORMATIONS                                               VARI-245
      NBT=10*NBET                                                       VARI-246
      IF (J.GT.NBT) GO TO 79                                            VARI-247
      K1=1+(J-1)/10                                                     VARI-248
      K2=1+MOD(J-1,10)                                                  VARI-249
      IF (K2.GT.8) GO TO 80                                             VARI-250
      LX(1)=.TRUE.                                                      VARI-251
      JK=NBETA+J-1-(J-10)/10                                            VARI-252
      GO TO 32                                                          VARI-253
   21 J=J-1000                                                          VARI-254
      IF (J.GT.1000) GO TO 22                                           VARI-255
C  NUCLEAR PARAMETERS                                                   VARI-256
      IF (J.GT.NVA) GO TO 81                                            VARI-257
      LO(111)=.TRUE.                                                    VARI-258
      JK=NPAA+J-1                                                       VARI-259
      GO TO 32                                                          VARI-260
   22 J=J-1000                                                          VARI-261
      IF (J.GT.1000) GO TO 23                                           VARI-262
C  NUCLEAR MATRIX ELEMENTS                                              VARI-263
      NME=(NIVQ-NT)/3                                                   VARI-264
      IF (J.GT.NME) GO TO 82                                            VARI-265
      LX(2)=.TRUE.                                                      VARI-266
      JK=NT+3*J-1                                                       VARI-267
      GO TO 32                                                          VARI-268
   23 J=J-1000                                                          VARI-269
      IF (J.GT.1000) GO TO 30                                           VARI-270
C  SPIN-ORBIT AND H. F. PARAMETRISATION                                 VARI-271
      IF ((.NOT.LO(81)).AND.(J.GT.6)) GO TO 83                          VARI-272
      IF (J.GT.16) GO TO 27                                             VARI-273
      IF (J.LE.6) GO TO 24                                              VARI-274
      LX(4)=LX(4).OR.(J.LE.11)                                          VARI-275
      LX(5)=LX(5).OR.(J.GT.11)                                          VARI-276
      IF (LO(82).AND.J.GT.9) GO TO 84                                   VARI-277
      IF ((.NOT.LO(86)).AND.(J.GT.11)) GO TO 85                         VARI-278
      IF (LO(82)) GO TO 25                                              VARI-279
      IF ((J.NE.9).AND.(J.LE.11).AND.LO(87)) GO TO 86                   VARI-280
      IF (LO(87)) GO TO 25                                              VARI-281
      IF (J.EQ.7) GO TO 87                                              VARI-282
      IF ((AZ(8).NE.0.D0).AND.(J.GT.8).AND.(J.LE.11)) GO TO 87          VARI-283
      IF ((AZ(8).EQ.0.D0).AND.(J.EQ.8)) GO TO 87                        VARI-284
      GO TO 25                                                          VARI-285
   24 IF (.NOT.LO(4)) GO TO 88                                          VARI-286
      LX(3)=.TRUE.                                                      VARI-287
   25 IF (KF.NE.0) GO TO 26                                             VARI-288
      IF (K.EQ.J1) RAP=DW(NXX+I-1)/AZ(J)                                VARI-289
      AZ(J)=RAP*AZ(J)                                                   VARI-290
      GO TO 34                                                          VARI-291
   26 DW(NXX+I-1)=AZ(J)                                                 VARI-292
      GO TO 35                                                          VARI-293
C LEVEL DENSITY PARAMETERS.                                             VARI-294
   27 J3=J-16                                                           VARI-295
      IF (J3.GT.6*NCONS) GO TO 28                                       VARI-296
      IF (J3.LE.6*(NCONS-NCONT)) LX(6)=.TRUE.                           VARI-297
      IF (J3.GT.6*(NCONS-NCONT)) LX(7)=.TRUE.                           VARI-298
      JK=NSCN+J3-(J3-7)/6                                               VARI-299
      GO TO 32                                                          VARI-300
   28 J3=J3-7*NCONS                                                     VARI-301
      IF (J3.GT.NRD) GO TO 29                                           VARI-302
      LX(8)=.TRUE.                                                      VARI-303
      JK=NGAM+J3-1                                                      VARI-304
      GO TO 32                                                          VARI-305
   29 J3=J3-NRD                                                         VARI-306
      IF (J3.GT.NFISS) GO TO 89                                         VARI-307
      LX(9)=.TRUE.                                                      VARI-308
      JK=NFIS+J-1                                                       VARI-309
      GO TO 32                                                          VARI-310
C DISPERSION PARAMETERS.                                                VARI-311
   30 J=J-1000                                                          VARI-312
      JN=(J-1)/13                                                       VARI-313
      JJ=J-13*JN                                                        VARI-314
      IF (JJ.GT.NPP-1) GO TO 90                                         VARI-315
      IF (NW(2,NIPP+17*JN+2).NE.0) GO TO 91                             VARI-316
      IF ((JJ.LT.6).OR.(JJ.GE.12)) GO TO 31                             VARI-317
      IF ((NW(1,NIPP+17*JN+1).EQ.0).AND.(J.LE.8)) GO TO 92              VARI-318
      IF ((NW(2,NIPP+17*JN+1).EQ.0).AND.(J.GT.8)) GO TO 92              VARI-319
   31 JK=NIPP+2+JJ+3*JN                                                 VARI-320
      LO(113)=.TRUE.                                                    VARI-321
   32 IF (KF.EQ.0) GO TO 33                                             VARI-322
      DW(NXX+I-1)=DW(JK)                                                VARI-323
      GO TO 35                                                          VARI-324
   33 IF (K.EQ.J1) RAP=DW(NXX+I-1)/DW(JK)                               VARI-325
      DW(JK)=RAP*DW(JK)                                                 VARI-326
   34 CONTINUE                                                          VARI-327
   35 CONTINUE                                                          VARI-328
      LO(112)=LX(4).OR.LX(5).OR.LX(6)                                   VARI-329
      LO(114)=LX(7)                                                     VARI-330
   36 IF (KF.EQ.0) GO TO 38                                             VARI-331
      WRITE (MW,1000) NW(1,NIW+1),NW(2,NIW),CHI2,(YY(I),I=1,3)          VARI-332
      WRITE (MW,1001)                                                   VARI-333
      DO 37 I=1,NCOLR                                                   VARI-334
      IF (DW(NTGX+7*I-1).NE.0.D0) WRITE (MW,1002) I,(DW(NTGX+7*I+J-8),J=VARI-335
     13,7)                                                              VARI-336
   37 CONTINUE                                                          VARI-337
      WRITE (MW,1003) (I,DW(NXX+I-1),I=1,NREC)                          VARI-338
      IF (LO(78)) WRITE (MW,1004) (I,DW(NRES+I-1),I=1,NTOT)             VARI-339
      RETURN                                                            VARI-340
   38 IF (KE.EQ.1) RETURN                                               VARI-341
C PRINTING OF RESULTS.                                                  VARI-342
      WRITE (MW,1005) NW(1,NIW+1),NW(2,NIW),CHI2,KE,YY(1)               VARI-343
      IF (KE*(KE-3).NE.0.OR.DW(NRC).LE.0.D0) GO TO 40                   VARI-344
      WRITE (MW,1006) NW(2,NIW+1),YY(3),(DW(NRC+I-1),I=1,NREC)          VARI-345
      WRITE (MW,1007) (DW(NRC+NREC+I-1),I=1,NREC)                       VARI-346
      WRITE (MW,1008)                                                   VARI-347
      L=2*NREC+NRC-1                                                    VARI-348
      DO 39 I=1,NREC                                                    VARI-349
      K=L+1                                                             VARI-350
      L=L+I                                                             VARI-351
   39 WRITE (MW,1009) (DW(J),J=K,L)                                     VARI-352
      GO TO 41                                                          VARI-353
   40 WRITE (MW,1010)                                                   VARI-354
   41 IF (KE.EQ.0) WRITE (MW,1011)                                      VARI-355
      IF (KE.EQ.2) WRITE (MW,1012)                                      VARI-356
      IF (KE.EQ.3) WRITE (MW,1013)                                      VARI-357
      IF (KE.EQ.4) WRITE (MW,1014)                                      VARI-358
      IF (KE.EQ.5) WRITE (MW,1015) NW(2,NIW+1)                          VARI-359
      IF (KE.EQ.6) WRITE (MW,1016) NW(1,NIW+1),NW(2,NIW+1)              VARI-360
      IF (KE.EQ.7) WRITE (MW,1017) NTOT,NREC                            VARI-361
C END OF THE SEARCH                                                     VARI-362
      LO(116)=.FALSE.                                                   VARI-363
      LO(118)=.TRUE.                                                    VARI-364
      IF (LO(76)) GO TO 43                                              VARI-365
      DO 42 I=51,65                                                     VARI-366
   42 LO(I)=LO(I+85)                                                    VARI-367
   43 LO(120)=.NOT.(LO(53).OR.LO(55).OR.LO(56).OR.LO(57).OR.LO(58).OR.LOVARI-368
     1(60).OR.LO(62).OR.LO(64))                                         VARI-369
      IF (LO(41)) LO(120)=LO(120).AND.(.NOT.LO(65))                     VARI-370
      LO(119)=LO(120).AND.(.NOT.(LO(51).OR.LO(64).OR.LO(65)))           VARI-371
      LO(120)=CHI2.EQ.CHI2M.AND.LO(120)                                 VARI-372
C OUTPUT OF FINAL PARAMETERS                                            VARI-373
      WRITE (MW,1018) (LO(J),J=1,100)                                   VARI-374
      IF (LO(7)) GO TO 46                                               VARI-375
      NPO=NPOT-1                                                        VARI-376
C OUTPUT OF USUAL POTENTIALS.                                           VARI-377
      DO 45 J=1,NPP                                                     VARI-378
      IJ=IABS(NW(1,NIPP+15*J-15))                                       VARI-379
      AM3=DW(NWV+18*IJ-17)**.33333333333333D0                           VARI-380
      IF (LO(16)) AM3=AM3+DW(NWV+18*IJ-18)**.33333333333333D0           VARI-381
      WRITE (MW,1019) J,AM3                                             VARI-382
      DO 44 I=1,8                                                       VARI-383
   44 R0(I)=DW(NPO+4*I-2)/AM3                                           VARI-384
      WRITE (MW,1020) (DW(NPO+4*I-3),DW(NPO+4*I-2),R0(I),DW(NPO+4*I-1),DVARI-385
     1W(NPO+4*I),I=1,8),DW(NPO+33)                                      VARI-386
      IF (LO(17)) WRITE (MW,1021) (DW(NPO+I),I=34,42)                   VARI-387
   45 NPO=NPO+1                                                         VARI-388
      GO TO 58                                                          VARI-389
C OUTPUT OF EXTERNAL POTENTIALS.                                        VARI-390
   46 NVMA=NW(1,NPOT)                                                   VARI-391
      WRITE (MW,1022)                                                   VARI-392
   47 IF (NVMA.GE.NW(1,NPOT+1)) GO TO 57                                VARI-393
      I1=NW(1,NVMA+NPOT-1)                                              VARI-394
      IV=NW(2,NVMA+NPOT-1)                                              VARI-395
      IF (IV.EQ.16) GO TO 56                                            VARI-396
      IT1=MOD(I1-1,8)+1                                                 VARI-397
      J1=(I1-1)/8                                                       VARI-398
      I2=NW(1,NPOT+I1-7)                                                VARI-399
      I3=NW(2,NPOT+I1-7)                                                VARI-400
      IF (IV.GE.0) GO TO 48                                             VARI-401
      WRITE (MW,1023) IT1,J1,I2,I3,DW(NPOT+I2-1),DW(NPOT+I2)            VARI-402
      GO TO 56                                                          VARI-403
   48 IF (IV.LT.9) GO TO 49                                             VARI-404
      WRITE (MW,1024) IT1,J1,IV,I2,I3,(DW(NPOT+I-1),I=I2,I3)            VARI-405
      GO TO 56                                                          VARI-406
   49 IF (IV.GT.6) GO TO 55                                             VARI-407
      WRITE (MW,1025) IT1,J1,IV,I2,I3,(DW(NPOT+I-1),I=I2,I3)            VARI-408
      NST=NW(1,NPOT+NVMA+1)                                             VARI-409
      IF (NST.GT.0) GO TO 52                                            VARI-410
      K=IABS(NST)                                                       VARI-411
      EX=DW(NWV+18*K-17)**.33333333333333D0                             VARI-412
      EY=EX                                                             VARI-413
      IF (LO(16)) EX=EX+DW(NWV+18*K-18)**.33333333333333D0              VARI-414
      EY=EY/EX                                                          VARI-415
      EX=DW(NPOT+I2)/EX                                                 VARI-416
      EZ=DW(NPOT+I2-1)                                                  VARI-417
      IF ((.NOT.LO(16)).OR.J1.LE.NPP) GO TO 51                          VARI-418
      ITYZ=IV                                                           VARI-419
      IF (ITYZ.GE.5) ITYZ=ITYZ-4                                        VARI-420
      ITYW=1                                                            VARI-421
      IF (IT1.LE.6) GO TO 50                                            VARI-422
      ITI=7*(J1-NPP)                                                    VARI-423
      ITYW=ITYW*NW(2-MOD(ITI,2),NIVY+(ITI-1)/2)                         VARI-424
   50 IF (LO(6)) ITYW=ITYW-1                                            VARI-425
      IF (ITYZ.GT.1) EZ=EZ/EY**((ITYZ-1)*ITYW)                          VARI-426
   51 WRITE (MW,1026) EZ,EX                                             VARI-427
   52 IF ((IV.NE.5).AND.(IV.NE.6)) GO TO 56                             VARI-428
      WRITE (MW,1027) I2,I3,(DW(NPOT+I-1),I=I2,I3)                      VARI-429
      IF (NST.GT.0) GO TO 56                                            VARI-430
      WRITE (MW,1028)                                                   VARI-431
      NMB=0                                                             VARI-432
   53 NMA=I2                                                            VARI-433
      I4=I2+7                                                           VARI-434
      DO 54 I=I2,I4                                                     VARI-435
      NMB=NMB+1                                                         VARI-436
      J=I-NMA                                                           VARI-437
      IF (IT1.LT.7) J=0                                                 VARI-438
      IF (.NOT.LO(6)) J=J+1                                             VARI-439
   54 R0(NMB)=DW(NPOT+I-1)/EY**J                                        VARI-440
      WRITE (MW,1009) (R0(I),I=1,NMB)                                   VARI-441
      I2=I2+8                                                           VARI-442
      IF (I2.LE.I3) GO TO 53                                            VARI-443
      GO TO 56                                                          VARI-444
   55 WRITE (MW,1025) IT1,J1,IV,I2,I3,(DW(NPOT+I-1),I=I2,I3)            VARI-445
   56 NVMA=I3+1                                                         VARI-446
      GO TO 47                                                          VARI-447
   57 NVMB=NW(2,NPOT+1)                                                 VARI-448
      IF (NVMA.LT.NVMB) WRITE (MW,1029) (DW(NPOT+I-1),I=NVMA,NVMB)      VARI-449
   58 IF (.NOT.LX(1)) GO TO 66                                          VARI-450
C OUTPUT OF DEFORMATIONS.                                               VARI-451
      WRITE (MW,1030) (I,NW(1,NBETA+9*I-1),NW(2,NBETA+9*I-1),(DW(NBETA+9VARI-452
     1*I+J-10),J=1,8),I=1,NBET)                                         VARI-453
      IF (.NOT.LO(16)) GO TO 66                                         VARI-454
      DM=DW(NWV+1)**.33333333333333D0/(DW(NWV)**.33333333333333D0+DW(NWVVARI-455
     1+1)**.33333333333333D0)                                           VARI-456
      WRITE (MW,1031)                                                   VARI-457
      DO 65 I=1,NBET                                                    VARI-458
      K1=0                                                              VARI-459
      K2=0                                                              VARI-460
      IF (LO(3)) GO TO 60                                               VARI-461
   59 K1=1                                                              VARI-462
      K2=NW(1,NBETA+9*I-1)                                              VARI-463
      GO TO 62                                                          VARI-464
   60 IF (LO(1)) GO TO 61                                               VARI-465
      K1=I-1                                                            VARI-466
      K2=K1*NW(1,NBETA+9*I-1)                                           VARI-467
      GO TO 62                                                          VARI-468
   61 IF (NW(2,NBETA+9*I-1).EQ.0) GO TO 59                              VARI-469
   62 IF (.NOT.LO(6)) GO TO 63                                          VARI-470
      K2=K2-K1                                                          VARI-471
      K1=0                                                              VARI-472
   63 DO 64 J=1,6                                                       VARI-473
   64 R0(J)=DW(NBETA+9*I+J-10)/DM**K1                                   VARI-474
      R0(7)=DW(NBETA+9*I-3)/DM**K2                                      VARI-475
      R0(8)=DW(NBETA+9*I-2)/DM**K2                                      VARI-476
   65 WRITE (MW,1032) I,NW(1,NBETA+9*I-1),NW(2,NBETA+9*I-1),R0          VARI-477
C OUTPUT OF OTHER PARAMETERS IN SEARCH                                  VARI-478
   66 IF (LO(111)) WRITE (MW,1033) (DW(NPAA+J-1),I=1,NVA)               VARI-479
      IF (LX(2)) WRITE (MW,1034) (DW(I+2),I=NT,NIVQ,3)                  VARI-480
      IF (LX(3)) WRITE (MW,1035) (AZ(I),I=1,6)                          VARI-481
      IF (LX(4)) WRITE (MW,1036) (AZ(I),I=7,11)                         VARI-482
      IF (LX(5)) WRITE (MW,1037) (AZ(I),I=12,16)                        VARI-483
      IF (LX(6)) WRITE (MW,1038) (I,(DW(NSCN+7*I+J-8),J=1,7),I=1,1)     VARI-484
      IF (LX(7)) WRITE (MW,1038) (I,(DW(NSCN+7*I+J-8),J=1,7),I=1+NCONS-NVARI-485
     1CONT,NCONS)                                                       VARI-486
      IF (LX(8)) WRITE (MW,1039) (I,DW(NGAM+I-1),I=1,NRD)               VARI-487
      IF (LX(9)) WRITE (MW,1040) (I,(DW(NFIS+2*I+J-3),J=1,2),I=1,NFISS) VARI-488
      IF (LO(113)) WRITE (MW,1041) (I,(DW(NIPP+15*I+J-16),J=4,16),I=1,NPVARI-489
     1P)                                                                VARI-490
      RETURN                                                            VARI-491
   67 WRITE (MW,1042) NVAT                                              VARI-492
      GO TO 93                                                          VARI-493
   68 WRITE (MW,1043)                                                   VARI-494
      GO TO 93                                                          VARI-495
   69 WRITE (MW,1044) NW(2,NPOT+1)                                      VARI-496
      GO TO 93                                                          VARI-497
   70 WRITE (MW,1045) I,J,(NW(1,NPOT+L),NW(2,NPOT+L),L=1,I1)            VARI-498
      GO TO 93                                                          VARI-499
   71 WRITE (MW,1046)                                                   VARI-500
      GO TO 93                                                          VARI-501
   72 WRITE (MW,1047)                                                   VARI-502
      GO TO 93                                                          VARI-503
   73 WRITE (MW,1048)                                                   VARI-504
      GO TO 93                                                          VARI-505
   74 WRITE (MW,1049)                                                   VARI-506
      GO TO 93                                                          VARI-507
   75 WRITE (MW,1050)                                                   VARI-508
      GO TO 93                                                          VARI-509
   76 WRITE (MW,1051) I,NW(2-MOD(K,2),NISE+(K-1)/2),J                   VARI-510
      GO TO 93                                                          VARI-511
   77 WRITE (MW,1052) NBET                                              VARI-512
      GO TO 93                                                          VARI-513
   78 WRITE (MW,1053) I,NW(2-MOD(K,2),NISE+(K-1)/2),J                   VARI-514
      GO TO 93                                                          VARI-515
   79 WRITE (MW,1054) NBT                                               VARI-516
      GO TO 93                                                          VARI-517
   80 WRITE (MW,1055)                                                   VARI-518
      GO TO 93                                                          VARI-519
   81 WRITE (MW,1056) NVA                                               VARI-520
      GO TO 93                                                          VARI-521
   82 WRITE (MW,1057) NME                                               VARI-522
      GO TO 93                                                          VARI-523
   83 WRITE (MW,1058)                                                   VARI-524
      GO TO 93                                                          VARI-525
   84 WRITE (MW,1059)                                                   VARI-526
      GO TO 93                                                          VARI-527
   85 WRITE (MW,1060)                                                   VARI-528
      GO TO 93                                                          VARI-529
   86 WRITE (MW,1061)                                                   VARI-530
      GO TO 93                                                          VARI-531
   87 WRITE (MW,1062)                                                   VARI-532
      GO TO 93                                                          VARI-533
   88 WRITE (MW,1063)                                                   VARI-534
      GO TO 93                                                          VARI-535
   89 WRITE (MW,1064)                                                   VARI-536
      GO TO 93                                                          VARI-537
   90 WRITE (MW,1065) NPP                                               VARI-538
      GO TO 93                                                          VARI-539
   91 WRITE (MW,1066)                                                   VARI-540
      GO TO 93                                                          VARI-541
   92 WRITE (MW,1067)                                                   VARI-542
   93 WRITE (MW,1068) I,NW(2-MOD(K,2),NISE+(K-1)/2),J                   VARI-543
      STOP                                                              VARI-544
 1000 FORMAT (/' RUN',I4,'   MAX =',I4,'   ***** CHI2 =',D18.10,' *****'VARI-545
     1,5X,'W(1) =',F10.2,5X,'W(2) =',F5.2,5X,'W(3) =',F10.5)            VARI-546
 1001 FORMAT (/21X,'WEIGHT',12X,'EXP. NORM.',10X,'ERR. NORM.',10X,'CALC.VARI-547
     1 NORM.',12X,'CHI2')                                               VARI-548
 1002 FORMAT (5X,I5,1P,5D20.6)                                          VARI-549
 1003 FORMAT (/' *** VARIABLES'//(6(1X,I3,1P,D16.6)))                   VARI-550
 1004 FORMAT (/' *** FUNCTIONS'//(6(1X,I3,1P,D16.6)))                   VARI-551
 1005 FORMAT (/' RUN',I4,'   MAX =',I4,'   ***** CHI2 =',D18.10,' *****'VARI-552
     1,5X,'KE =',I2,5X,'W(1) =',F12.4)                                  VARI-553
 1006 FORMAT (//' STANDARD ERRORS (VARIANCE AT BEST FIT EQUAL TO DEGREE VARI-554
     1OF FREEDOM:',I6,'.RENORMALISATION FACTOR',D15.6,' )'/(1P,8D15.6)) VARI-555
 1007 FORMAT (/' ERROR ENHANCEMENTS (MULTI/SINGLE VARIABLE ERROR)'/(1P,8VARI-556
     1D15.6))                                                           VARI-557
 1008 FORMAT (/' ERROR CORRELATION MATRIX:')                            VARI-558
 1009 FORMAT (1P,8D15.6)                                                VARI-559
 1010 FORMAT (//' NO INFORMATION ON ERRORS.')                           VARI-560
 1011 FORMAT (//' SEARCH ENDED WITHOUT ERRORS.')                        VARI-561
 1012 FORMAT (//' SEARCH INTERRUPTED BY USER.')                         VARI-562
 1013 FORMAT (//' SEARCH ENDED BY NUMBER OF EVALUATIONS.')              VARI-563
 1014 FORMAT (//' SEARCH ENDED FOR ROUNDING ERRORS.')                   VARI-564
 1015 FORMAT (//' SEARCH ENDED BECAUSE THE FUNCTIONS DO NOT DEPEND ON THVARI-565
     1E VARIABLE',I6)                                                   VARI-566
 1016 FORMAT (//' SEARCH ENDED BECAUSE VARIABLES',2I6,' ARE USELESS IN PVARI-567
     1REPARATORY CALLS.')                                               VARI-568
 1017 FORMAT (//' SEARCH ENDED BECAUSE THE NUMBER OF PARAMETERS',I4,' ISVARI-569
     1 LARGER THAN THE NUMBER OF DATA',I4)                              VARI-570
 1018 FORMAT ('1'/' ******* FINAL RESULTS *******'//' **** FIRST CONTROLVARI-571
     1 CARD ****',2X,'1 ',9(' 1'),' 2 ',9(' 2'),' 3 ',9(' 3'),' 4 ',9(' VARI-572
     24'),' 5'/11X,5('  1 2 3 4 5 6 7 8 9 0')/11X,5(1X,10L2)//' *** SECOVARI-573
     3ND CONTROL CARD ****',2X,'1 ',9(' 1'),' 2 ',9(' 2'),' 3 ',9(' 3'),VARI-574
     4' 4 ',9(' 4'),' 5'/11X,5('  1 2 3 4 5 6 7 8 9 0')/11X,5(1X,10L2)/)VARI-575
 1019 FORMAT (/' OPTICAL POTENTIALS  **',I3,' **     REDUCED RADIUS MULTVARI-576
     1IPLIED BY  ',D15.6/)                                              VARI-577
 1020 FORMAT (' VOLUME/SCALAR REAL POTENTIAL'/' DEPTH',F12.6,' MEV   RADVARI-578
     1IUS',F10.6,' FERMI (REDUCED VALUE',F9.6,')   DIFFUSENESS',F9.6,' FVARI-579
     2ERMI AT THE POWER (1+',F9.6,')'/' VOLUME/SCALAR IMAGINARY POTENTIAVARI-580
     3L'/' DEPTH',F12.6,' MEV   RADIUS',F10.6,' FERMI (REDUCED VALUE',F9VARI-581
     4.6,')   DIFFUSENESS',F9.6,' FERMI AT THE POWER (1+',F9.6,')'/' SURVARI-582
     5FACE/VECTOR REAL POTENTIAL'/' DEPTH',F12.6,' MEV   RADIUS',F10.6,'VARI-583
     6 FERMI (REDUCED VALUE',F9.6,')   DIFFUSENESS',F9.6,' FERMI AT THE VARI-584
     7POWER (1+',F9.6,')'/' SURFACE/VECTOR IMAGINARY POTENTIAL'/' DEPTH'VARI-585
     8,F12.6,' MEV   RADIUS',F10.6,' FERMI (REDUCED VALUE',F9.6,')   DIFVARI-586
     9FUSENESS',F9.6,' FERMI AT THE POWER (1+',F9.6,')'/' REAL SPIN-ORBIVARI-587
     AT/TENSOR POTENTIAL'/' DEPTH',F12.6,' MEV   RADIUS',F10.6,' FERMI (VARI-588
     BREDUCED VALUE',F9.6,')   DIFFUSENESS',F9.6,' FERMI AT THE POWER (1VARI-589
     C+',F9.6,')'/' IMAGINARY SPIN-ORBIT/TENSOR POTENTIAL'/' DEPTH',F12.VARI-590
     D6,' MEV   RADIUS',F10.6,' FERMI (REDUCED VALUE',F9.6,')   DIFFUSENVARI-591
     EESS',F9.6,' FERMI AT THE POWER (1+',F9.6,')'/' COULOMB POTENTIAL'/VARI-592
     F'  CHARGE PRODUCT',F7.0,'  RADIUS',F10.6,' FERMI (REDUCED VALUE',FVARI-593
     G9.6,')   DIFFUSENESS',F9.6,' FERMI AT THE POWER (1+',F9.6,')'/' SPVARI-594
     HIN-ORBIT COULOMB POTENTIAL'/' DEPTH',F12.6,' MEV   RADIUS',F10.6,'VARI-595
     I FERMI (REDUCED VALUE',F9.6,')   DIFFUSENESS',F9.6,' FERMI AT THE VARI-596
     JPOWER (1+',F9.6,')'/14X,'THIRD CHARGE PARAMETER',F9.6)            VARI-597
 1021 FORMAT (/' *** FOLDING MODEL ***'/' REAL PART',10X,'V =',F10.4,6X,VARI-598
     1'R =',F10.4,6X,'A =',F10.4/' IMAGINARY PART',6X,'V =', F10.4,6X,'RVARI-599
     2 =',F10.4,6X,'A =',F10.4/' COULOMB PART',7X,'V =',F10.4,6X,'R =',FVARI-600
     310.4,6X,'A =',F10.4)                                              VARI-601
 1022 FORMAT (/' **** EXTERNAL POTENTIAL PARAMETERS ****'/)             VARI-602
 1023 FORMAT (' (',I1,',',I2,') GIVEN BY POINTS WITH THE PARAMETERS FROMVARI-603
     1',I4,' TO',I4/' STRENGTH =',1P,D15.6,19X,'SCALE =',D15.6)         VARI-604
 1024 FORMAT (' (',I1,',',I2,') TYPE',I3,' FROM',I4,' TO',I4/1P,2D15.6/(VARI-605
     11P,8D15.6))                                                       VARI-606
 1025 FORMAT (' (',I1,',',I2,') TYPE',I3,' FROM',I4,' TO',I4,4X,6D14.6/(VARI-607
     17X,8D14.6))                                                       VARI-608
 1026 FORMAT (' VALUES READ:',F12.6,3X,F9.6)                            VARI-609
 1027 FORMAT (' DEFORMATIONS: PARAMETERS FROM',I4,' TO ',I4,10X,4D18.6/(VARI-610
     152X,4D18.6))                                                      VARI-611
 1028 FORMAT (' DEFORMATIONS READ:')                                    VARI-612
 1029 FORMAT (' FOLDING PARAMETERS:',3F10.5/(20X,3F10.5))               VARI-613
 1030 FORMAT (/' BETA(I,J) FOR  L   K',9X,'V',9X,'W',8X,'VS',8X,'WS',7X,VARI-614
     1'VSO',7X,'WSO',6X,'COUL S.O. COUL'/(5X,I5,5X,I2,2X,I2,2X,8F10.5)) VARI-615
 1031 FORMAT (/10X,'WITHOUT HEAVY ION CORRECTION:'/)                    VARI-616
 1032 FORMAT (5X,I5,5X,I2,2X,I2,2X,8F10.5)                              VARI-617
 1033 FORMAT (/' *** NUCLEAR VARIABLES ***'//(1X,6F20.6))               VARI-618
 1034 FORMAT (/' *** NUCLEAR MATRIX ELEMENTS ***'//(1X,6F20.6))         VARI-619
 1035 FORMAT (/' *** SPIN-ORBIT PARAMETRISATION ***'//1X,6F20.6/)       VARI-620
 1036 FORMAT (/' *** HAUSER-FESHBACH CORRECTIONS ***'//1X,5F20.6/)      VARI-621
 1037 FORMAT (/' *** GIANT DIPOLE RESONANCE PARAM. ***'//1X,5F20.6/)    VARI-622
 1038 FORMAT (/' *** DENSITY OF STATES ***'//(1X,I3,'  SA:',D13.6,6X,'UXVARI-623
     1:',D13.6,5X,'TAU:',D13.6,6X,'SG:',D13.6/28X,'E0:',D13.6,6X,'EX:',DVARI-624
     213.6,7X,'Z:',F5.0))                                               VARI-625
 1039 FORMAT (/' *** GAMMA TRANSMISSION COEFFICIENTS ***'//(1X,4(I5,F20.VARI-626
     16)))                                                              VARI-627
 1040 FORMAT (/' *** FISSION DATA ***'//(1X,2(I5,2F20.6)))              VARI-628
 1041 FORMAT (/' *** DISPERSION RELATION PARAMETERS ***'//(' POTENTIAL',VARI-629
     1I5/1X,6F20.6/1X,6F20.6/1X,F20.6))                                 VARI-630
 1042 FORMAT (' NO OPTICAL PARAMETER INDEX LARGER THAN',I3)             VARI-631
 1043 FORMAT (' NO CHANGE OF COULOMB CHARGES.')                         VARI-632
 1044 FORMAT (' NO EXTERNAL OPTICAL PARAMETER LARGER THAN',I5)          VARI-633
 1045 FORMAT (' VARIABLE',I3,' DEFINED BY',I5,' CANNOT BE USED.'/' THE PVARI-634
     1ARAMETER WHICH CAN BE SEARCHED ARE BETWEEN:'/(4(5X,I6,' TO',I6))) VARI-635
 1046 FORMAT (' FOR A FORM FACTOR GIVEN BY POINTS, THE TWO MULTIPLICATIVVARI-636
     1E FACTOR ONLY CAN BE CHANGED.')                                   VARI-637
 1047 FORMAT (' THE COUPLING OF TWO PARTICLE STATES CANNOT BE CHANGED.')VARI-638
 1048 FORMAT (' ONLY THE OSCILLATOR PARAMETER CAN BE CHANGED FOR A LAGUEVARI-639
     1RRE POLYNOMIAL.')                                                 VARI-640
 1049 FORMAT (' MASSES AND PRODUCT OF CHARGES CANNOT BE CHANGED FOR A SIVARI-641
     1NGLE PARTICLE STATE.')                                            VARI-642
 1050 FORMAT (' NO INDEX FOR DEFORMATION OF A GIVEN POTENTIAL LARGER THAVARI-643
     1N 8.')                                                            VARI-644
 1051 FORMAT (' THE VARIABLE',I3,' CANNOT BE USED BECAUSE ',I5,'  IS A VVARI-645
     1ARIATION OF THE DEFORMATIONS OF POTENTIAL',I3,' WHICH ARE ZERO.') VARI-646
 1052 FORMAT (' NO INDEX FOR MULTIPOLE LARGER THAN',I3)                 VARI-647
 1053 FORMAT (' THE VARIABLE',I3,' CANNOT BE USED BECAUSE ',I5,'  IS A VVARI-648
     1ARIATION FOR A GIVEN MULTIPOLE',I3,' AND THEY ARE ZERO.')         VARI-649
 1054 FORMAT (' NO INDEX FOR INDIVIDUAL MULTIPOLE LARGER THAN',I3)      VARI-650
 1055 FORMAT (' NO INDEX FOR INDIVIDUAL MULTIPOLE EQUAL TO 0 OR 9 MODULOVARI-651
     1 10.')                                                            VARI-652
 1056 FORMAT (' NO INDEX FOR NUCLEAR PARAMETER LARGER THAN',I3)         VARI-653
 1057 FORMAT (' NO INDEX FOR NUCLEAR MATRIX ELEMENT LARGER THAN',I3)    VARI-654
 1058 FORMAT (' NO HAUSER-FESHBACH CORRECTIONS.')                       VARI-655
 1059 FORMAT (' PARAMETER NOT USED FOR SIMPLIFIED COMPOUND NUCLEUS.')   VARI-656
 1060 FORMAT (' NO GAMMA EMISSION IN COMPOUND NUCLEUS.')                VARI-657
 1061 FORMAT (' NO WIDTH FLUCTUATION IN COMPOUND NUCLEUS.')             VARI-658
 1062 FORMAT (' PARAMETER NOT USED IN THIS KIND OF COMPOUND NUCLEUS.')  VARI-659
 1063 FORMAT (' NO PARAMETRISATION OF SOIN-ORBIT DEFORMATION.')         VARI-660
 1064 FORMAT (' THIS LEVEL DENSITY PARAMETER IS NOT USED.')             VARI-661
 1065 FORMAT (' THERE ARE ONLY',I3,' SETS OF DISPERSION PARAMETERS.')   VARI-662
 1066 FORMAT (' THE COEFFICIENTS OF THIS DISPERSION RELATION SHOULD BE RVARI-663
     1EAD EACH TIME AND CANNOT BE IN SEARCH.')                          VARI-664
 1067 FORMAT (' THIS COEFFICIENT OF DISPERSION RELATION IS NOT USED.')  VARI-665
 1068 FORMAT (' VARIABLE',I3,' DEFINED BY',I5,I4,' CANNOT BE USED.'///' VARI-666
     1IN VARI  ...  STOP  ...')                                         VARI-667
      END                                                               VARI-668
