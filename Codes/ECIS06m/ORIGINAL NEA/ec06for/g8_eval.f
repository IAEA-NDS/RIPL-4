C 29/10/06                                                      ECIS06  EVAL-000
      SUBROUTINE EVAL(NW,DW,CM,LO)                                      EVAL-001
C THIS SUBROUTINE CHANGES SOME PARAMETERS TO DO A NEW CALCULATION.      EVAL-002
C THE VALUES READ HERE ARE ABSOLUTE CHANGES ( NEX=0 ),                  EVAL-003
C RELATIVE CHANGES ( NEX>0 ) OR PERCENTAGES ( NEX<0 ).                  EVAL-004
C IF THE LABORATORY ENERGY IS CHANGED, COULOMB FUNCTIONS AND REDUCED    EVAL-005
C NUCLEAR MATRIX ELEMENTS ARE RECALCULATED EVEN IF THEY ARE NOT MODIFIEDEVAL-006
C IF NUCLEAR PARAMETERS ARE CHANGED, REDUCED NUCLEAR MATRIX ELEMENTS AREEVAL-007
C CALCULATED AGAIN. IN ALL THE OTHERS CASES, THE COMPUTATION RESTARTS   EVAL-008
C WITH THE COMPUTATION OF POTENTIALS.                                   EVAL-009
C INDEXES FOR PARAMETERS ARE THE ONES USED IN SEARCH (SEE VARI).        EVAL-010
C HOWEVER 0 MEANS ENERGY IN THE LABORATORY SYSTEM.                      EVAL-011
C INPUT:     NW:      WORKING AREA FOR INTEGERS.                        EVAL-012
C            DW:      WORKING AREA FOR DOUBLE PRECISION IN EQUIVALENCE  EVAL-013
C                     BY CALL WITH NW.                                  EVAL-014
C            CM:      NUCLEAR MASS.                                     EVAL-015
C            LO:      LOGICAL CONTROLS:                                 EVAL-016
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).EVAL-017
C               LO(3)  =.TRUE. ANHARMONIC VIBRATIONAL OR ASYMMETRIC     EVAL-018
C                              ROTATIONAL MODEL.                        EVAL-019
C               LO(4)  =.TRUE. PARAMETRISED SPIN-ORBIT DEFORMATION.     EVAL-020
C               LO(6)  =.TRUE. USE DEFORMATION LENGTHS.                 EVAL-021
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    EVAL-022
C               LO(8)  =.TRUE. RELATIVISTIC KINEMATICS.                 EVAL-023
C               LO(16) =.TRUE. HEAVY-ION DEFINITION OF REDUCED RADII ANDEVAL-024
C                              DEFORMATIONS.                            EVAL-025
C               LO(32) =.TRUE. AUTOMATIC SEARCH ON SOME PARAMETERS.     EVAL-026
C               LO(34) =.TRUE. NEXT CALCULATION CHANGING ENERGY AND/OR  EVAL-027
C                              SOME PARAMETERS.                         EVAL-028
C               LO(75) =.TRUE. NO COMPLETE OUTPUT AT THE FIRST RUN OF A EVAL-029
C                              SEARCH.                                  EVAL-030
C               LO(76) =.TRUE. LO(51) TO LO(65) ARE ALWAYS USED.        EVAL-031
C               LO(81) =.TRUE. HAUSER-FESHBACH CORRECTIONS.             EVAL-032
C               LO(82) =.TRUE. OLD SIMPLIFIED COMPOUND NUCLEUS.         EVAL-033
C               LO(86) =.TRUE. GAMMA EMISSION IN COMPOUND NUCLEUS.      EVAL-034
C               LO(87) =.TRUE. NO WIDTH FLUCTUATIONS.                   EVAL-035
C               LO(111)=.TRUE. NUCLEAR PARAMETERS ARE CHANGED IN SEARCH.EVAL-036
C               LO(112)=.TRUE. SPIN-ORBIT OR COMPOUND NUCLEUS PARAMETERSEVAL-037
C                              ARE CHANGED IN SEARCH.                   EVAL-038
C               LO(113)=.TRUE. DISPERSION RELATION IS CHANGED IN SEARCH.EVAL-039
C               LO(114)=.TRUE. COMPOUND NUCLEUS CONTINUUM IS CHANGED IN EVAL-040
C                              SEARCH.                                  EVAL-041
C               LO(115)=.TRUE. FIRST COMPUTATION FOR THIS ENERGY.       EVAL-042
C               LO(116)=.TRUE. NO OUTPUT.                               EVAL-043
C                                                                       EVAL-044
C FOR THE COMMON  /ADDRE/, /COUPL/, /INTEG/ AND /TITRE/ SEE CALC.       EVAL-045
C FOR THE COMMON  /DCHI2/, AND /NCOMP/ SEE CALX.                        EVAL-046
C                                                                       EVAL-047
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /ADDRE/:                     EVAL-048
C  NWV:       NON INTEGER VALUES FOR THE CHANNELS.                      EVAL-049
C  NIPP:      DISPERSION PARAMETERS.                                    EVAL-050
C  NSCN:      LEVEL DENSITY DESCRIPTION.                                EVAL-051
C  NPOT:      OPTICAL POTENTIAL PARAMETERS.                             EVAL-052
C  NBETA:     DEFORMATION PARAMETERS.                                   EVAL-053
C  NIW:       INTEGER WORKING FIELD FOR THE SEARCH.                     EVAL-054
C  NT:        TABLE OF REDUCED NUCLEAR MATRIX ELEMENTS.                 EVAL-055
C  NIVQ:      TABLE OF MULTIPOLES.                                      EVAL-056
C  NIVY:      TABLE OF FORM FACTOR IDENTIFICATION IVY (FOR COMPUTATION).EVAL-057
C  NCX:       FIRST FREE ADDRESS FOR COMPUTATION OF POTENTIALS.         EVAL-058
C             HERE, FIRST FREE ADDRESS FOR INPUT.                       EVAL-059
C   USED:     NWV,NIPP,NSCN,NPOT,NBETA,NIW,NT,NIVQ,NIVY,NCX.            EVAL-060
C                                                                       EVAL-061
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     EVAL-062
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             EVAL-063
C  NVA:       NUMBER OF NUCLEAR PARAMETERS.                             EVAL-064
C   USED:     NPP,NVA.                                                  EVAL-065
C                                                                       EVAL-066
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCHI2/:                     EVAL-067
C  CHI2M:     MINIMUM CHI-SQUARE IN THE SEARCH.                         EVAL-068
C  YY(1):     STEP SIZE IN THE SEARCH.                                  EVAL-069
C  YY(3):     VARIOUS MEANINGS.  SEE FITE.                              EVAL-070
C   DEFINED:  CHI2M,YY.                                                 EVAL-071
C                                                                       EVAL-072
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /INTEG/:                     EVAL-073
C  IDMT:      TOTAL WORKING FIELD LENGTH AS SINGLE PRECISION.           EVAL-074
C  NBET:      NUMBER OF DIFFERENT DEFORMATIONS (VIBRATIONS+ROTATIONS).  EVAL-075
C   USED:     IDMT,NBET.                                                EVAL-076
C                                                                       EVAL-077
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /NCOMP/:                     EVAL-078
C  NFISS:     NUMBER OF FISSION TRANSMISSION COEFFICIENTS.              EVAL-079
C  NRD:       NUMBER OF GAMMA TRANSMISSION COEFFICIENTS.                EVAL-080
C  NCONT:     NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS.                 EVAL-081
C  NCONS:     NUMBER OF LEVEL DENSITIES NEEDED.                         EVAL-082
C  NCOLX:     TOTAL NUMBER OF LEVELS WITHOUT DISCRETISATION.            EVAL-083
C  AZ(I):     SPIN-ORBIT PARAMETRISATION FOR I-1 TO 6,                  EVAL-084
C             HAUSER FESHBACH PARAMETERS FOR J=7 TO 11,                 EVAL-085
C             GIANT DIPOLE RESONANCE FOR I=12 TO 16.                    EVAL-086
C   DEFINED:  AZ.                                                       EVAL-087
C   USED:     NFISS,NRD,NCONT,NCONS,NCOLX,AZ.                           EVAL-088
C                                                                       EVAL-089
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /TITRE/:                     EVAL-090
C  TITLE(18): TITLE OF THE RUN PRINTED AS HEADING OF RESULTS.           EVAL-091
C   DEFINED:  TITLE.                                                    EVAL-092
C                                                                       EVAL-093
C MEANING OF INDEX AND LOGICAL RETURNED:                                EVAL-094
C      1-1000      OPTICAL MODEL,FOLDING PARAMETERS.                    EVAL-095
C   1001-2000      DEFORMATIONS FOR A GIVEN POTENTIAL.                  EVAL-096
C   2001-3000      DEFORMATIONS FOR A GIVEN MULTIPOLE.                  EVAL-097
C   3001-4000      INDIVIDUAL DEFORMATION.                              EVAL-098
C   4001-5000      NUCLEAR MODEL PARAMETER.             - LO(111)=.TRUE.EVAL-099
C   5001-6000      NUCLEAR MATRIX ELEMENT.                              EVAL-100
C   6001-7000      SPIN-ORBIT AND COMPOUND NUCLEUS PARAMETRISATION.     EVAL-101
C                  SPIN-ORBIT PARAMETRISATION:                          EVAL-102
C                  BZ1, BZ2, BZ3, BZ4, BZ5              - LO(112)=.TRUE.EVAL-103
C                  TGO, BN, FNUG, EGD, GGD              - LO(112)=.TRUE.EVAL-104
C                  SA, UX, TAU, SG, E0, EX FOR GAMMA    - LO(112)=.TRUE.EVAL-105
C                  SA, UX, TAU, SG, E0, EX FOR CONTINUUM- LO(114)=.TRUE.EVAL-106
C                  GAMMA TRANSMISSION FACTORS                           EVAL-107
C                  FISSION TRANSMISSION COEFFICIENT, DEGREE OF FREEDOM  EVAL-108
C   7001-8000      DISPERSION RELATIONS PARAMETRISATION - LO(113)=.TRUE.EVAL-109
C  10001-99999     EXTERNAL OPTICAL MODEL (PARAMETERS ABOVE 1000)       EVAL-110
C IT STOPS THE CALCULATION FOR AN INDEX OF PARAMETER OUT OF RANGE       EVAL-111
C                                                                       EVAL-112
C***********************************************************************EVAL-113
      IMPLICIT REAL*8 (A-H,O-Z)                                         EVAL-114
      LOGICAL LO(150),LX                                                EVAL-115
      DIMENSION NW(2,*),DW(*)                                           EVAL-116
      CHARACTER*4 AA(2),TITLE                                           EVAL-117
      CHARACTER*8 BB(2,2)                                               EVAL-118
      COMMON /ADDRE/ NIPH,NJIT,NWV,NIPP,NPAR,NPAA,NSCN,NFIS,NGAM,NNIV,NPEVAL-119
     1OT,NBETA,NFM,NTGX,NDONN,NNVI,NDE,NISE,NRC,NIW,NNWI,NCC,MCC,NXA,NAMEVAL-120
     21,NFAC,NFG,NXG,NRES,NXX,NT,NIVQ,NIVY,NIVZ,NCOI,MIPI,NXD,MWV,NIXT,NEVAL-121
     3TY,NSM,NTX,NRY,NRCO,NRDO,NVC1,NVC2,NNC,NCX                        EVAL-122
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   EVAL-123
      COMMON /DCHI2/ CHI2,CHI2M,YY(3)                                   EVAL-124
      COMMON /INOUT/ MR,MW,MS                                           EVAL-125
      COMMON /INTEG/ IDMT,NPLACE,NCOLL,NJMAX,ITERM,JDM,JIT,KMIN,KMAX,NCOEVAL-126
     1LS,NCOLT,NBET,LMX,LMAX1,NLT,ISM,NJC,JTX,KCC,MS1,MS2,KBA,KAB,KBC,JTEVAL-127
     2H,NCOLR,NREC,NTOT,LMAX2,KE,ITEMM,KXT,LMAX3,NRZ,NTZ,IPM,IPK,MCM(2),EVAL-128
     3NCT(6)                                                            EVAL-129
      COMMON /NCOMP/ NSP(3),NFISS,NRD,NCONT,NCOJ,NCONS,NIE,NCOLX,NDP,NDQEVAL-130
     1,ACN1,ACN2,AZ(18)                                                 EVAL-131
      COMMON /TITRE/ TITLE(18)                                          EVAL-132
      DATA AA /' NEW','LAST'/                                           EVAL-133
      DATA BB /'INCREASE','S )     ','NEW VALU','ES )    '/             EVAL-134
      NVAT=42*NPP                                                       EVAL-135
      LO(34)=.FALSE.                                                    EVAL-136
      LX=.FALSE.                                                        EVAL-137
      READ (MR,1000) LO(34),LX,NIN,NEX,YY(1)                            EVAL-138
      IF (LX) READ (MR,1001) TITLE                                      EVAL-139
      NX=MIN0(1,MAX0(NEX,-1))+2                                         EVAL-140
      I1=1                                                              EVAL-141
      I2=1                                                              EVAL-142
      IF (.NOT.LO(34)) I1=2                                             EVAL-143
      IF (NX.EQ.2) I2=2                                                 EVAL-144
      WRITE (MW,1002) AA(I1),NIN,BB(1,I2),BB(2,I2)                      EVAL-145
      IF (NX.EQ.1) WRITE (MW,1003)                                      EVAL-146
      DO 1 I=111,120                                                    EVAL-147
    1 LO(I)=.FALSE.                                                     EVAL-148
      IF (NCX+2*NIN.GT.IDMT) CALL MEMO('EVAL',IDMT,NCX+2*NIN)           EVAL-149
      READ (MR,1004) (NW(1,2*I+NCX-2),I=1,NIN)                          EVAL-150
      READ (MR,1005) (DW(2*I+NCX-1),I=1,NIN)                            EVAL-151
      DO 50 I=1,NIN                                                     EVAL-152
      NW(1,2*I+NCX-2)=MAX0(0,NW(1,2*I+NCX-2))                           EVAL-153
      J=MOD(NW(1,2*I+NCX-2),1000)                                       EVAL-154
      IK=(NW(1,2*I+NCX-2)+999)/1000                                     EVAL-155
      EY=DW(2*I+NCX-1)                                                  EVAL-156
      IF (IK.NE.0) GO TO 7                                              EVAL-157
C ENERGY IN THE LABORATORY SYSTEM.                                      EVAL-158
      LO(115)=.TRUE.                                                    EVAL-159
      DW(NWV+12)=EY                                                     EVAL-160
      EX=DW(NWV+2)                                                      EVAL-161
      IF (NX.EQ.1) EY=.01D0*EY*EX                                       EVAL-162
      IF (NX.NE.2) EY=EX+EY                                             EVAL-163
      IF (LO(8)) GO TO 2                                                EVAL-164
      DW(NWV+2)=EY*DW(NWV+1)/(DW(NWV)+DW(NWV+1))                        EVAL-165
      GO TO 3                                                           EVAL-166
C RELATIVISTIC C.-M. ENERGY ECM=SQRT((M1+M2)**2+2*M2*ELAB))-M1-M2.      EVAL-167
    2 DW(NWV+2)=CM*(DSQRT((DW(NWV)+DW(NWV+1))**2+2.D0*DW(NWV+1)*EY/CM)-DEVAL-168
     1W(NWV)-DW(NWV+1))                                                 EVAL-169
    3 WRITE (MW,1006) I,J,EY,DW(2*I+NCX-1),DW(NWV+2),EX                 EVAL-170
      IF (NCOLX.EQ.1) GO TO 50                                          EVAL-171
      DO 5 J=2,NCOLX                                                    EVAL-172
      DW(NWV+22*J-20)=DW(NWV+22*J-20)+DW(NWV+2)-EX                      EVAL-173
      IF (LO(8)) GO TO 4                                                EVAL-174
      DW(NWV+22*J-10)=DW(NWV+22*J-20)*(DW(NWV+22*J-22)+DW(NWV+22*J-21))/EVAL-175
     1DW(NWV+22*J-21)                                                   EVAL-176
      GO TO 5                                                           EVAL-177
    4 DW(NWV+22*J-10)=DW(NWV+22*J-20)*(DW(NWV+22*J-20)/(2.D0*CM)+DW(NWV+EVAL-178
     122*J-22)+DW(NWV+22*J-21))/DW(NWV+22*J-21)                         EVAL-179
    5 CONTINUE                                                          EVAL-180
      GO TO 50                                                          EVAL-181
    6 IF (NEX.LT.0) EY=.01D0*EY*EX                                      EVAL-182
      EY=EX+EY                                                          EVAL-183
      GO TO ( 18 , 23 , 29 , 34 , 36 , 38 , 40 , 49 ),IK                EVAL-184
    7 IF (IK.GT.10) GO TO 8                                             EVAL-185
      IF (IK.GT.8) GO TO 52                                             EVAL-186
      GO TO ( 9 , 19 , 25 , 31 , 35 , 37 , 39 , 48 ),IK                 EVAL-187
    8 IK=1                                                              EVAL-188
      J=NW(1,2*I+NCX-2)-9000                                            EVAL-189
    9 IF (LO(7)) GO TO 11                                               EVAL-190
C OPTICAL MODEL AND FOLDING PARAMETERS.                                 EVAL-191
      N=1+MOD(J-1,42)                                                   EVAL-192
      M=1+(J-1)/42                                                      EVAL-193
      K=IABS(NW(1,NIPP+17*M-17))                                        EVAL-194
      IF (J.GT.NVAT.OR.N.EQ.25) GO TO 53                                EVAL-195
      IF (MOD(N,4).NE.2.OR.N.GT.32) GO TO 10                            EVAL-196
      EX=DW(NWV+22*K-21)**.33333333333333D0                             EVAL-197
      IF (LO(16)) EX=EX+DW(NWV+22*K-22)**.33333333333333D0              EVAL-198
      EY=EY*EX                                                          EVAL-199
   10 EX=DW(NPOT+J-1)                                                   EVAL-200
      IF (NX.NE.2) GO TO 6                                              EVAL-201
      GO TO 18                                                          EVAL-202
C EXTERNAL OPTICAL PARAMETERS.                                          EVAL-203
   11 IF (J.GE.NW(2,NPOT+1)) GO TO 54                                   EVAL-204
      I1=NW(1,NPOT)-2                                                   EVAL-205
      DO 12 L=1,I1                                                      EVAL-206
      IF (J.GE.NW(1,NPOT+L).AND.J.LE.NW(2,NPOT+L)) GO TO 13             EVAL-207
   12 CONTINUE                                                          EVAL-208
      GO TO 55                                                          EVAL-209
   13 EX=DW(NPOT+J-1)                                                   EVAL-210
      IF (L.EQ.1) GO TO 17                                              EVAL-211
C NOT A FOLDING PARAMETER.                                              EVAL-212
      M=NW(1,NPOT+L)                                                    EVAL-213
      N=NW(2,NPOT+M-2)                                                  EVAL-214
      IV=NW(2,NPOT+N-1)                                                 EVAL-215
      IF (IV.LT.1.OR.IV.GT.6) GO TO 16                                  EVAL-216
      NST=NW(1,NPOT+N+1)                                                EVAL-217
      IF (NST.GT.0) GO TO 18                                            EVAL-218
      IT1=MOD(NW(1,NPOT+N-1)-1,8)+1                                     EVAL-219
      J1=(NW(1,NPOT+N-1)-1)/8                                           EVAL-220
      J2=4                                                              EVAL-221
      IF (IT1.GT.6) J2=5                                                EVAL-222
      IF (J-M.GT.1.AND.J-M.LT.J2) GO TO 17                              EVAL-223
      K=IABS(NST)                                                       EVAL-224
      FX=DW(NWV+22*K-21)**.33333333333333D0                             EVAL-225
      FY=FX                                                             EVAL-226
      IF (LO(16)) FX=FX+DW(NWV+22*K-22)**.33333333333333D0              EVAL-227
C TRANSFORMATION OF DEPTH AND RADIUS.                                   EVAL-228
      IF (J-M.EQ.1) EY=EY*FX                                            EVAL-229
      FY=FY/FX                                                          EVAL-230
      IF (J.GT.M+2) GO TO 15                                            EVAL-231
      IF ((.NOT.LO(16)).OR.J1.LE.NPP.OR.J.NE.M) GO TO 17                EVAL-232
      ITYZ=IV                                                           EVAL-233
      IF (ITYZ.GE.5) ITYZ=ITYZ-4                                        EVAL-234
      ITYW=1                                                            EVAL-235
      IF (IT1.LE.6) GO TO 14                                            EVAL-236
      ITI=7*(J1-NPP)                                                    EVAL-237
      ITYW=ITYW*NW(2-MOD(ITI,2),NIVY+(ITI-1)/2)                         EVAL-238
   14 IF (LO(6)) ITYW=ITYW-1                                            EVAL-239
      IF (ITYZ.GT.1) EY=EY*FY**((ITYZ-1)*ITYW)                          EVAL-240
      GO TO 17                                                          EVAL-241
C TRANSFORMATION OF DEFORMATIONS.                                       EVAL-242
   15 J3=J-M-J2                                                         EVAL-243
      IF (IT1.LT.7) J3=0                                                EVAL-244
      IF (.NOT.LO(6)) J3=J3+1                                           EVAL-245
      EY=EY*FY**J3                                                      EVAL-246
      GO TO 17                                                          EVAL-247
   16 IF ((IV.LT.0).AND.(J-M.GT.1)) GO TO 56                            EVAL-248
      IF ((IV.NE.7).AND.(IV.NE.8)) GO TO 17                             EVAL-249
      JY=J                                                              EVAL-250
      IF (NW(1,NPOT+NW(2,NPOT+M-2)).NE.1) JY=JY-1                       EVAL-251
      IF (M.GT.JY) GO TO 57                                             EVAL-252
      IF ((IV.EQ.7).AND.(MOD(JY-M,3).NE.0)) GO TO 58                    EVAL-253
      IF ((IV.EQ.8).AND.(MOD(JY+10-M,11).LT.3)) GO TO 59                EVAL-254
   17 IF (NX.NE.2) GO TO 6                                              EVAL-255
   18 IF (.NOT.LO(7)) WRITE (MW,1007) I,J,N,M,EY,DW(2*I+NCX-1),DW(NPOT+JEVAL-256
     1-1)                                                               EVAL-257
      IF (LO(7)) WRITE (MW,1008) I,J,EY,DW(2*I+NCX-1),DW(NPOT+J-1)      EVAL-258
      DW(NPOT+J-1)=EY                                                   EVAL-259
      GO TO 50                                                          EVAL-260
   19 IF (J.GT.8) GO TO 60                                              EVAL-261
C DEFORMATIONS FOR A GIVEN POTENTIAL.                                   EVAL-262
      K2=J                                                              EVAL-263
      DO 20 K1=1,NBET                                                   EVAL-264
      IF (LO(1).AND.LO(3).AND.NW(2,NBETA+9*K1-1).NE.0) GO TO 20         EVAL-265
      IF (DW(NBETA+K2+9*K1-10).NE.0.D0) GO TO 21                        EVAL-266
   20 CONTINUE                                                          EVAL-267
      GO TO 61                                                          EVAL-268
   21 EY=EY/DW(NBETA+K2+9*K1-10)                                        EVAL-269
      K3=1                                                              EVAL-270
      IF ((NX.NE.1).AND.LO(16)) GO TO 32                                EVAL-271
   22 EX=1                                                              EVAL-272
      IF (NX.NE.2) GO TO 6                                              EVAL-273
   23 WRITE (MW,1009) I,NW(1,2*I+NCX-2),K2,K1,EY,DW(2*I+NCX-1)          EVAL-274
      DO 24 L=K1,NBET                                                   EVAL-275
      IF (LO(1).AND.LO(3).AND.NW(2,NBETA+9*K1-1).NE.0) GO TO 24         EVAL-276
      EX=DW(NBETA+K2+9*L-10)*EY                                         EVAL-277
      WRITE (MW,1010) K2,L,EX,DW(NBETA+K2+9*L-10)                       EVAL-278
      DW(NBETA+K2+9*L-10)=EX                                            EVAL-279
   24 CONTINUE                                                          EVAL-280
      GO TO 50                                                          EVAL-281
C DEFORMATIONS FOR A GIVEN MULTIPOLE.                                   EVAL-282
   25 K1=J                                                              EVAL-283
      IF (J.GT.NBET) GO TO 62                                           EVAL-284
      DO 26 K2=1,8                                                      EVAL-285
      IF (DW(NBETA+K2+9*K1-10).NE.0.D0) GO TO 27                        EVAL-286
   26 CONTINUE                                                          EVAL-287
      GO TO 63                                                          EVAL-288
   27 EY=EY/DW(NBETA+K2+9*K1-10)                                        EVAL-289
      K3=2                                                              EVAL-290
      IF ((NX.NE.1).AND.LO(16)) GO TO 22                                EVAL-291
   28 EX=1                                                              EVAL-292
      IF (NX.NE.2) GO TO 6                                              EVAL-293
   29 WRITE (MW,1009) I,NW(1,2*I+NCX-2),K2,K1,EY,DW(2*I+NCX-1)          EVAL-294
      DO 30 L=K2,8                                                      EVAL-295
      EX=DW(NBETA+L+9*K1-10)*EY                                         EVAL-296
      WRITE (MW,1010) L,K1,EX,DW(NBETA+L+9*K1-10)                       EVAL-297
   30 DW(NBETA+L+9*K1-10)=EX                                            EVAL-298
      GO TO 50                                                          EVAL-299
C INDIVIDUAL DEFORMATIONS.                                              EVAL-300
   31 NBT=10*NBET                                                       EVAL-301
      IF (J.GT.NBT) GO TO 64                                            EVAL-302
      K1=1+(J-1)/10                                                     EVAL-303
      K2=1+MOD(J-1,10)                                                  EVAL-304
      IF (K2.GT.8) GO TO 65                                             EVAL-305
      K3=3                                                              EVAL-306
      IF (.NOT.LO(16)) GO TO 33                                         EVAL-307
C SEARCH OF THE CORRECTIONS FOR HEAVY IONS.                             EVAL-308
   32 EX=DW(NWV+1)**.33333333333333D0/(DW(NWV+1)**.33333333333333D0+DW(NEVAL-309
     1WV)**.33333333333333D0)                                           EVAL-310
      JK=1                                                              EVAL-311
      IF (.NOT.LO(1).AND.LO(3)) JK=K1-1                                 EVAL-312
      K=JK                                                              EVAL-313
      IF (K2.GT.6) K=K*NW(1,NBETA+9*K1-1)                               EVAL-314
      IF (LO(6)) K=K-JK                                                 EVAL-315
      EY=EY*EX**K                                                       EVAL-316
      IF (K3-2) 22 , 28 , 33                                            EVAL-317
   33 EX=DW(NBETA+K2+9*K1-10)                                           EVAL-318
      IF (NX.NE.2) GO TO 6                                              EVAL-319
   34 WRITE (MW,1011) I,NW(1,2*I+NCX-2),K2,K1,EY,DW(2*I+NCX-1),DW(NBETA+EVAL-320
     1K2+9*K1-10)                                                       EVAL-321
      DW(NBETA+K2+9*K1-10)=EY                                           EVAL-322
      GO TO 50                                                          EVAL-323
   35 IF (J.GT.NVA) GO TO 66                                            EVAL-324
      EX=DW(NPAA+J-1)                                                   EVAL-325
      IF (NX.NE.2) GO TO 6                                              EVAL-326
   36 WRITE (MW,1012) I,NW(1,2*I+NCX-2),J,EY,DW(2*I+NCX-1),DW(NPAA+J-1) EVAL-327
      DW(NPAA+J-1)=EY                                                   EVAL-328
      GO TO 50                                                          EVAL-329
C  NUCLEAR MATRIX ELEMENTS.                                             EVAL-330
   37 I3=(NIVQ-NT)/3                                                    EVAL-331
      IF (J.GT.I3) GO TO 67                                             EVAL-332
      EX=DW(NT+3*J-3)                                                   EVAL-333
      LO(111)=.TRUE.                                                    EVAL-334
      IF (NX.NE.2) GO TO 6                                              EVAL-335
   38 WRITE (MW,1013) I,NW(1,2*I+NCX-2),J,EY,DW(2*I+NCX-1),DW(NT+3*J-3) EVAL-336
      DW(NT+3*J-3)=EY                                                   EVAL-337
      GO TO 50                                                          EVAL-338
C SPIN-ORBIT PARAMETRISATION.                                           EVAL-339
   39 IF (J.GT.6) GO TO 41                                              EVAL-340
      IF (.NOT.LO(4)) GO TO 68                                          EVAL-341
      EX=AZ(J)                                                          EVAL-342
      IF (NX.NE.2) GO TO 6                                              EVAL-343
   40 IF (J.GT.6) GO TO 43                                              EVAL-344
      WRITE (MW,1014) I,NW(1,2*I+NCX-2),J,EY,DW(2*I+NCX-1),AZ(J)        EVAL-345
      GO TO 45                                                          EVAL-346
C HAUSER-FESBACH CORRECTION.                                            EVAL-347
   41 IF (.NOT.LO(81)) GO TO 69                                         EVAL-348
      IF (J.GT.16) GO TO 46                                             EVAL-349
      IF (LO(82).AND.J.GT.9) GO TO 70                                   EVAL-350
      IF ((.NOT.LO(86)).AND.J.GT.11) GO TO 70                           EVAL-351
      IF (LO(82)) GO TO 42                                              EVAL-352
      IF (J.NE.9.AND.J.LE.11.AND.LO(87)) GO TO 70                       EVAL-353
      IF (LO(87)) GO TO 42                                              EVAL-354
      IF (J.EQ.7) GO TO 70                                              EVAL-355
      IF (AZ(8).NE.0.D0.AND.J.GT.8.AND.J.LE.11) GO TO 70                EVAL-356
      IF (AZ(8).EQ.0.D0.AND.J.EQ.8) GO TO 70                            EVAL-357
   42 EX=AZ(J)                                                          EVAL-358
      IF (NX.NE.2) GO TO 6                                              EVAL-359
   43 IF (J.GT.16) GO TO 47                                             EVAL-360
      LO(112)=.TRUE.                                                    EVAL-361
      K=J-6                                                             EVAL-362
      IF (K.GT.6) GO TO 44                                              EVAL-363
      WRITE (MW,1015) I,NW(1,2*I+NCX-2),K,EY,DW(2*I+NCX-1),AZ(J)        EVAL-364
      GO TO 45                                                          EVAL-365
   44 K=J-5                                                             EVAL-366
      WRITE (MW,1016) I,NW(1,2*I+NCX-2),K,EY,DW(2*I+NCX-1),AZ(J)        EVAL-367
   45 AZ(J)=EY                                                          EVAL-368
      GO TO 50                                                          EVAL-369
   46 JJ=J-16+(J-16)/6                                                  EVAL-370
      EX=DW(NSCN+JJ-1)                                                  EVAL-371
      IF (NX.NE.2) GO TO 6                                              EVAL-372
   47 K=1+(J-17)/NCONS                                                  EVAL-373
      L=1+MOD(J-17,NCONS)                                               EVAL-374
      IF (K.EQ.NCONS-NCONT) LO(112)=.TRUE.                              EVAL-375
      IF (K.GT.NCONS-NCONT) LO(114)=.TRUE.                              EVAL-376
      WRITE (MW,1017) I,NW(1,2*I+NCX-2),L,K,EY,DW(2*I+NCX-1),DW(NSCN+JJ-EVAL-377
     11)                                                                EVAL-378
      DW(NSCN+JJ-1)=EY                                                  EVAL-379
      GO TO 50                                                          EVAL-380
C DISPERSION PARAMETERS.                                                EVAL-381
   48 IF (.NOT.LO(10)) GO TO 71                                         EVAL-382
      K=(J-1)/13+1                                                      EVAL-383
      IF (K.GT.NPP) GO TO 72                                            EVAL-384
      IF (NW(2,NIPP+17*K-14).NE.0) GO TO 73                             EVAL-385
      L=J-13*K+17                                                       EVAL-386
      LO(113)=.TRUE.                                                    EVAL-387
      EX=DW(NIPP+17*K+L-17)                                             EVAL-388
      IF ((EX.EQ.0.D0).OR.(EY.EQ.0.D0)) WRITE (MW,1018) EX,EY           EVAL-389
      IF (NX.NE.2) GO TO 6                                              EVAL-390
   49 WRITE (MW,1019) I,NW(1,2*I+NCX-2),L,K,EY,DW(2*I+NCX-1),DW(NIPP+17*EVAL-391
     1K+L-17)                                                           EVAL-392
      DW(NIPP+17*K+L-17)=EY                                             EVAL-393
   50 CONTINUE                                                          EVAL-394
      IF (.NOT.LO(32)) RETURN                                           EVAL-395
      CHI2M=1.D35                                                       EVAL-396
      IF (YY(1).EQ.0.D0) YY(1)=20.D0                                    EVAL-397
      YY(3)=0.D0                                                        EVAL-398
      NW(2,NIW)=NW(2,NIW)-NW(1,NIW+1)                                   EVAL-399
      NW(1,NIW+1)=1                                                     EVAL-400
      WRITE (MW,1020) NW(2,NIW),YY(1)                                   EVAL-401
      IF (NW(2,NIW).LE.0) GO TO 75                                      EVAL-402
      II=51                                                             EVAL-403
      IF (LO(76).OR.(.NOT.LO(75))) II=59                                EVAL-404
      DO 51 I=II,65                                                     EVAL-405
      LO(I+85)=LO(I)                                                    EVAL-406
   51 LO(I)=.FALSE.                                                     EVAL-407
      LO(116)=.TRUE.                                                    EVAL-408
      RETURN                                                            EVAL-409
   52 WRITE (MW,1021)                                                   EVAL-410
      GO TO 74                                                          EVAL-411
   53 WRITE (MW,1022) NVAT                                              EVAL-412
      GO TO 74                                                          EVAL-413
   54 WRITE (MW,1023) NW(2,NPOT+1)                                      EVAL-414
      GO TO 74                                                          EVAL-415
   55 WRITE (MW,1024) I,EY,J,(NW(1,NPOT+L),NW(2,NPOT+L),L=1,I1)         EVAL-416
      GO TO 74                                                          EVAL-417
   56 WRITE (MW,1025)                                                   EVAL-418
      GO TO 74                                                          EVAL-419
   57 WRITE (MW,1026)                                                   EVAL-420
      GO TO 74                                                          EVAL-421
   58 WRITE (MW,1027)                                                   EVAL-422
      GO TO 74                                                          EVAL-423
   59 WRITE (MW,1028)                                                   EVAL-424
      GO TO 74                                                          EVAL-425
   60 WRITE (MW,1029)                                                   EVAL-426
      GO TO 74                                                          EVAL-427
   61 WRITE (MW,1030) I,NW(1,2*I+NCX-2),K2                              EVAL-428
      GO TO 74                                                          EVAL-429
   62 WRITE (MW,1031) NBET                                              EVAL-430
      GO TO 74                                                          EVAL-431
   63 WRITE (MW,1032) I,NW(1,2*I+NCX-2),K1                              EVAL-432
      GO TO 74                                                          EVAL-433
   64 WRITE (MW,1033) NBT                                               EVAL-434
      GO TO 74                                                          EVAL-435
   65 WRITE (MW,1034)                                                   EVAL-436
      GO TO 74                                                          EVAL-437
   66 WRITE (MW,1035) NVA                                               EVAL-438
      GO TO 74                                                          EVAL-439
   67 WRITE (MW,1036) I3                                                EVAL-440
      GO TO 74                                                          EVAL-441
   68 WRITE (MW,1037)                                                   EVAL-442
      GO TO 74                                                          EVAL-443
   69 WRITE (MW,1038)                                                   EVAL-444
      GO TO 74                                                          EVAL-445
   70 WRITE (MW,1039)                                                   EVAL-446
      GO TO 74                                                          EVAL-447
   71 WRITE (MW,1040)                                                   EVAL-448
      GO TO 74                                                          EVAL-449
   72 WRITE (MW,1041) NPP                                               EVAL-450
      GO TO 74                                                          EVAL-451
   73 WRITE (MW,1042)                                                   EVAL-452
   74 WRITE (MW,1043) I,NW(1,2*I+NCX-2),IK,J,DW(2*I+NCX-1)              EVAL-453
   75 WRITE (MW,1044)                                                   EVAL-454
      STOP                                                              EVAL-455
 1000 FORMAT (2L1,I3,I5,F10.5)                                          EVAL-456
 1001 FORMAT (18A4)                                                     EVAL-457
 1002 FORMAT ('1',A4,' COMPUTATION WITH',I4,'  NEW PARAMETERS',10X,'( INEVAL-458
     1PUT OF ',2A8)                                                     EVAL-459
 1003 FORMAT (20X,'*** VALUES GIVEN AS PERCENTAGES ***')                EVAL-460
 1004 FORMAT (14I5)                                                     EVAL-461
 1005 FORMAT (7F10.5)                                                   EVAL-462
 1006 FORMAT (2X,I2,' PARAM.',I5,3X,'NEW LAB. ENERGY',1P,D15.8,' (',D15.EVAL-463
     18,')',3X,'CENTRE OF MASS ENERGY',D15.6,5X,'OLD VALUE',D15.6)      EVAL-464
 1007 FORMAT (2X,I2,' PARAM.',I5,3X,'VALUE OF V-OPTICAL(',I2,',',I2,') =EVAL-465
     1',1P,D15.6,' (',D15.6,')',3X,'OLD VALUE',D15.6)                   EVAL-466
 1008 FORMAT (2X,I2,' PARAM.',I5,'TH VALUE OF EXTERNAL POTENTIAL =',1P,DEVAL-467
     115.6,' (',D15.6,')',3X,'OLD VALUE',D15.6)                         EVAL-468
 1009 FORMAT (2X,I2,' PARAM.',I5,3X,'PROPORTIONAL TO BETA(',I2,',',I2,')EVAL-469
     1  WITH RATIO',1P,D15.6,' (',D15.6,')')                            EVAL-470
 1010 FORMAT (15X,'BETA(',I2,',',I2,') =',1P,D15.6,5X,'OLD VALUE',D15.6)EVAL-471
 1011 FORMAT (2X,I2,' PARAM.',I5,3X,'BETA(',I2,',',I2,') =',1P,D15.6,' (EVAL-472
     1',E15.6,')',3X,'OLD VALUE',D15.6)                                 EVAL-473
 1012 FORMAT (2X,I2,' PARAM.',I5,3X,'VAR(',I2,') =',1P,D15.6,' (',D15.6,EVAL-474
     1')',3X,'OLD VALUE',D15.6)                                         EVAL-475
 1013 FORMAT (2X,I2,' PARAM.',I5,3X,'T(4,',I3,') =',1P,D15.6,' (',D15.6,EVAL-476
     1')',3X,'OLD VALUE',D15.6)                                         EVAL-477
 1014 FORMAT (2X,I2,' PARAM.',I5,3X,'AZ(',I1,') =',1P,D15.6,' (',D15.6,'EVAL-478
     1)',3X,'OLD VALUE',D15.6)                                          EVAL-479
 1015 FORMAT (2X,I2,' PARAM.',I5,3X,'BZ(',I1,') =',1P,D15.6,' (',D15.6,'EVAL-480
     1)',3X,'OLD VALUE',D15.6)                                          EVAL-481
 1016 FORMAT (2X,I2,' PARAM.',I5,3X,'GIANT RESONANCE(',I1,') =',1P,D15.6EVAL-482
     1,' (',D15.6,')',3X,'OLD VALUE',D15.6)                             EVAL-483
 1017 FORMAT (2X,I2,' PARAM.',I5,3X,'SCN(',I3,',',I3,') =',1P,D15.6,' ('EVAL-484
     1,D15.6,')',3X,'OLD VALUE',D15.6)                                  EVAL-485
 1018 FORMAT (' NO VERIFICATION OF VALIDITY FOR THE DISPERSION PARAMETEREVAL-486
     1',1P,D15.6,' CHANGED INTO',D15.6)                                 EVAL-487
 1019 FORMAT (2X,I2,' PARAM.',I5,3X,'PIP(',I3,',',I3,') =',1P,D15.6,' ('EVAL-488
     1,D15.6,')',3X,'OLD VALUE',D15.6)                                  EVAL-489
 1020 FORMAT (' NEW MAXIMUM NUMBER OF RUNS:',I6,10X,'STARTING SCALE',F10EVAL-490
     1.2)                                                               EVAL-491
 1021 FORMAT (' NO VARIABLE IS DEFINED BY INDEX BETWEEN 8001 AND 10000.'EVAL-492
     1)                                                                 EVAL-493
 1022 FORMAT (' NO OPTICAL PARAMETER INDEX LARGER THAN',I3,' AND NO CHANEVAL-494
     1GE OF COULOMB CHARGES.')                                          EVAL-495
 1023 FORMAT (' NO EXTERNAL OPTICAL PARAMETER LARGER THAN',I5)          EVAL-496
 1024 FORMAT (' VARIABLE',I3,1P,D15.6,' DEFINED BY',I5,' CANNOT BE USED'EVAL-497
     1/' THE PARAMETER WHICH CAN BE SEARCHED ARE BETWEEN:'/(4(5X,I6,' TOEVAL-498
     2',I6)))                                                           EVAL-499
 1025 FORMAT (' FOR A FORM FACTOR GIVEN BY POINTS, THE TWO MULTIPLICATIVEVAL-500
     1E FACTOR ONLY CAN BE CHANGED.')                                   EVAL-501
 1026 FORMAT (' THE COUPLING OF TWO PARTICLE STATES CANNOT BE CHANGED.')EVAL-502
 1027 FORMAT (' ONLY THE OSCILLATOR PARAMETER CAN BE CHANGED FOR A LAGUEEVAL-503
     1RRE POLYNOMIAL.')                                                 EVAL-504
 1028 FORMAT (' MASSES AND PRODUCT OF CHARGES CANNOT BE CHANGED FOR A SIEVAL-505
     1NGLE PARTICLE STATE.')                                            EVAL-506
 1029 FORMAT (' NOT AN INDEX FOR DEFORMATION OF A GIVEN POTENTIAL.')    EVAL-507
 1030 FORMAT (' THE VARIABLE',I3,' CANNOT BE USED BECAUSE ',I5,'  IS A VEVAL-508
     1ARIATION OF THE DEFORMATIONS OF POTENTIAL',I3,' WHICH ARE ZERO.') EVAL-509
 1031 FORMAT (' NO INDEX FOR MULTIPOLE LARGER THAN',I3)                 EVAL-510
 1032 FORMAT (' THE VARIABLE',I3,' CANNOT BE USED BECAUSE ',I5,'  IS A VEVAL-511
     1ARIATION FOR A GIVEN MULTIPOLE',I3,' AND THEY ARE ZERO.')         EVAL-512
 1033 FORMAT (' NO INDEX FOR INDIVIDUAL MULTIPOLE LARGER THAN',I3)      EVAL-513
 1034 FORMAT (' NO INDEX FOR INDIVIDUAL MULTIPOLE EQUAL TO 0 OR 9 MODULOEVAL-514
     1 10')                                                             EVAL-515
 1035 FORMAT (' NO INDEX FOR NUCLEAR LEVEL PARAMETER LARGER THAN',I3)   EVAL-516
 1036 FORMAT (' NO INDEX FOR NUCLEAR MATRIX ELEMENT LARGER THAN',I3)    EVAL-517
 1037 FORMAT (' THE SPIN-ORBIT IS NOT PARAMETRISED.')                   EVAL-518
 1038 FORMAT (' NO HAUSER-FESBACH CORRECTION TO CROSS-SECTIONS.')       EVAL-519
 1039 FORMAT (' PARAMETER NOT USED FOR THESE HAUSER-FESBACH CORRECTION TEVAL-520
     1O CROSS-SECTIONS.')                                               EVAL-521
 1040 FORMAT (' DISPERSION RELATIONS WERE NOT USED IN PREVIOUS COMPUTATIEVAL-522
     1ON.')                                                             EVAL-523
 1041 FORMAT (' THERE ARE ONLY',I3,' DIFFERENT POTENTIALS.')            EVAL-524
 1042 FORMAT (' THE COEFFICIENTS OF THIS DISPERSION RELATION SHOULD BE REVAL-525
     1EAD EACH TIME.')                                                  EVAL-526
 1043 FORMAT (' VARIABLE',I3,' DEFINED BY',I5,' (',I2,I5,')',' WITH VALUEVAL-527
     1E ',D15.6,' CANNOT BE USED.')                                     EVAL-528
 1044 FORMAT (//' IN EVAL  ...  STOP ...')                              EVAL-529
      END                                                               EVAL-530
