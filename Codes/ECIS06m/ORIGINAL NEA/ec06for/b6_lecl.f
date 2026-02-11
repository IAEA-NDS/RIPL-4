C 24/01/07                                                      ECIS06  LECL-000
      SUBROUTINE LECL(NCOLX,NCOLL,NCONT,IDT,LO,IPI,IPH,WV,IPP,NPA,PA,NA,LECL-001
     1NB,NIMAX,NBET)                                                    LECL-002
C INPUT OF LEVEL DESCRIPTIONS.                                          LECL-003
C INPUT:     NCOLX:   TOTAL NUMBER OF NUCLEAR STATES                    LECL-004
C            NCOLL:   NUMBER OF COUPLED NUCLEAR STATES                  LECL-005
C            NCONT:   NUMBER OF CONTINUUM FOR COMPOUND NUCLEUS          LECL-006
C            IDT:     FREE SPACE IN NPA,PA                              LECL-007
C            LO(I):   LOGICAL CONTROLS:                                 LECL-008
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).LECL-009
C               LO(2)  =.TRUE. SECOND ORDER VIBRATIONAL OR CONSTRAINED  LECL-010
C                              ASYMMETRIC ROTATIONAL MODEL.             LECL-011
C               LO(3)  =.TRUE. ANHARMONIC VIBRATIONAL OR ASYMMETRIC     LECL-012
C                              ROTATIONAL MODEL.                        LECL-013
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    LECL-014
C               LO(8)  =.TRUE. RELATIVISTIC KINEMATICS.                 LECL-015
C               LO(10) =.TRUE. DISPERSION RELATIONS FOR POTENTIALS.     LECL-016
C               LO(15) =.TRUE. REDUCED NUCLEAR MATRIX ELEMENTS READ.    LECL-017
C               LO(18) =.TRUE. PROJECTILE-TARGET ANTISYMMETRISATION.    LECL-018
C               LO(98) =.TRUE. Q ADDED TO THE MASS OF RESIDUAL NUCLEUS  LECL-019
C                              OR OUTGOING PARTICLE.                    LECL-020
C               LO(100)=.TRUE. DIRAC EQUATION.                          LECL-021
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    LECL-022
C               LO(122)=.TRUE. IDENTICAL PARTICLES WITHOUT SPIN.        LECL-023
C               LO(123)=.TRUE. IDENTICAL PARTICLES WITH SPIN.           LECL-024
C OUTPUT:    IPI:     DESCRIPTION OF LEVELS IN IPI(J,*):                LECL-025
C                     PARITY (0 FOR + AND 1 FOR -) FOR J=1.             LECL-026
C                     MULTIPLICITY OF INCIDENT PARTICLE FOR J=2.        LECL-027
C                     MULTIPLICITY OF THE TARGET FOR J=3.               LECL-028
C                     PRODUCT OF CHARGES FOR J=4.                       LECL-029
C                     INDEX OF POTENTIAL FOR J=5.                       LECL-030
C            IPH:     DESCRIPTION OF VIBRATIONAL LEVELS. FOR VIBRATIONALLECL-031
C                     MODEL, IPH(I,J) IS THE NUMBER OF PHONONS FOR      LECL-032
C                     I=1, J=1,NCOLL. IF IPH(1,J) IS 3, THE STATE       LECL-033
C                     IS A MIXTURE OF 1 AND 2-PHONONS STATES.           LECL-034
C                     FOR ROTATIONAL MODEL, IPH(1,J) IS 1 FOR A         LECL-035
C                     VIBRATIONAL BAND AND 2 FOR A MIXTURE WITH THE     LECL-036
C                     GROUND STATE BAND. SINGLE PHONON OR ADDRESS OF    LECL-037
C                     THE DESCRIPTION OF MIXED STATES IN IPH(2,J).      LECL-038
C            WV:      DESCRIPTION OF THE LEVELS IN WV(J,*):             LECL-039
C                     MASS OF INCIDENT PARTICLE FOR J=1.                LECL-040
C                     MASS OF THE TARGET FOR J=2.                       LECL-041
C                     ENERGY IN THE CENTRE OF MASS IN MEV FOR J=3.      LECL-042
C                     (TEMPORARY) PARTICLE EXCITATION FRACTION FOR J=4. LECL-043
C                     ENERGY IN THE LABORATORY SYSTEM IN MEV FOR J=13.  LECL-044
C            IPP(I,J):FIRST LEVEL USING POTENTIAL J FOR I=1 (TEMPORARY).LECL-045
C                     -1 TO READ DISPERSION PARAMETERS FOR I=2.         LECL-046
C            NPA,PA:  STORAGE OF NUCLEAR PARAMETERS EQUIVALENT BY CALL. LECL-047
C            NA:      NUMBER OF INTEGER PARAMETER INFORMATIONS.         LECL-048
C            NB:      STORAGE IN ROAM FOR ASYMMETRIC ROTATION.          LECL-049
C            NIMAX:   TWICE MAXIMUM SUM OF SPINS TARGET+PARTICLE.       LECL-050
C            NBET:    NUMBER OF DIFFERENT PHONONS.                      LECL-051
C            LO:      LOGICAL CONTROLS, DEFINED HERE:                   LECL-052
C                     LO(122) AND LO(123).                              LECL-053
C                                                                       LECL-054
C FOR THE COMMON  /COUPL/ AND /DCONS/ SEE CALC.                         LECL-055
C FOR THE COMMON  /POTE2/ SEE REDM.                                     LECL-056
C                                                                       LECL-057
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     LECL-058
C  CM:        ATOMIC MASS IN MEV.                                       LECL-059
C   USED:     CM.                                                       LECL-060
C                                                                       LECL-061
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     LECL-062
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             LECL-063
C  NVA:       NUMBER OF NUCLEAR PARAMETERS.                             LECL-064
C   DEFINED:  NVA.                                                      LECL-065
C   USED:     NPP.                                                      LECL-066
C                                                                       LECL-067
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     LECL-068
C  NPX:       NUMBER OF POTENTIALS TAKING INTO ACCOUNT DISPERSION.      LECL-069
C   DEFINED:  NTX.                                                      LECL-070
C                                                                       LECL-071
C***********************************************************************LECL-072
      IMPLICIT REAL*8 (A-H,O-Z)                                         LECL-073
      LOGICAL LO(150)                                                   LECL-074
      CHARACTER*1 SIGM,SPI                                              LECL-075
      DIMENSION IPI(11,*),IPH(2,*),WV(22,*),IPP(34,*),NPA(2,*),PA(*)    LECL-076
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   LECL-077
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            LECL-078
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         LECL-079
      COMMON /INOUT/ MR,MW,MS                                           LECL-080
      DATA SIGM /'-'/                                                   LECL-081
C OUTPUT OF NUCLEAR MODEL.                                              LECL-082
      IF (LO(7)) GO TO 3                                                LECL-083
      IF (LO(1)) GO TO 2                                                LECL-084
      IF (LO(3)) GO TO 1                                                LECL-085
      IF (.NOT.LO(2)) WRITE (MW,1000)                                   LECL-086
      IF (LO(2)) WRITE (MW,1001)                                        LECL-087
      GO TO 4                                                           LECL-088
    1 WRITE (MW,1002)                                                   LECL-089
      IF (.NOT.LO(15)) GO TO 29                                         LECL-090
      GO TO 4                                                           LECL-091
    2 IF (.NOT.LO(3)) WRITE (MW,1003)                                   LECL-092
      IF (LO(3)) WRITE (MW,1004)                                        LECL-093
      GO TO 4                                                           LECL-094
    3 WRITE (MW,1005)                                                   LECL-095
C INPUT OF LEVEL DESCRIPTION.                                           LECL-096
    4 NVA=0                                                             LECL-097
      NA=0                                                              LECL-098
      NB=0                                                              LECL-099
      NBET=0                                                            LECL-100
      NIMAX=0                                                           LECL-101
      DO 5 IV=1,NPP                                                     LECL-102
      IPP(1,IV)=-1                                                      LECL-103
    5 IPP(2,IV)=-1                                                      LECL-104
      IF (LO(98).AND.LO(8)) WRITE (MW,1006)                             LECL-105
      DO 27 IV=1,NCOLX                                                  LECL-106
      READ (MR,1007) SP2,N,K,SPI,E,SP1,WV(1,IV),WV(2,IV),SP3,WV(4,IV)   LECL-107
      IPI(2,IV)=1+IDINT(2.00001D0*SP1)                                  LECL-108
      IPI(3,IV)=1+IDINT(2.00001D0*SP2)                                  LECL-109
      IPI(4,IV)=IDINT(SP3)                                              LECL-110
      IF (IV.NE.1) GO TO 8                                              LECL-111
      WV(13,1)=E                                                        LECL-112
      IF (K.EQ.0) K=1                                                   LECL-113
      IF (LO(18).AND.((IPI(3,1).NE.IPI(2,1)).OR.(SPI.EQ.SIGM).OR.LO(100)LECL-114
     1)) GO TO 30                                                       LECL-115
      LO(122)=LO(18).AND.IPI(2,1)*IPI(3,1).EQ.1                         LECL-116
      LO(123)=LO(18).AND.IPI(2,1)*IPI(3,1).NE.1                         LECL-117
      IF (LO(8)) GO TO 6                                                LECL-118
      WV(3,1)=E*WV(2,1)/(WV(1,1)+WV(2,1))                               LECL-119
      WV(4,1)=DSQRT(CK*WV(3,1)*WV(1,1)*WV(2,1)/(WV(1,1)+WV(2,1)))       LECL-120
      GO TO 7                                                           LECL-121
C RELATIVISTIC C.-M. ENERGY ECM=SQRT((M1+M2)**2+2*M2*ELAB))-M1-M2.      LECL-122
    6 WV(3,1)=CM*(DSQRT((WV(1,1)+WV(2,1))**2+2.D0*WV(2,1)*E/CM)-WV(1,1)-LECL-123
     1WV(2,1))                                                          LECL-124
      WV(4,1)=DSQRT(0.125D0*CK*WV(3,1)*(WV(3,1)/CM+2.D0*WV(1,1)+2.D0*WV(LECL-125
     12,1))*(WV(3,1)/CM+2.D0*WV(1,1))*(WV(3,1)/CM+2.D0*WV(2,1)))/(WV(3,1LECL-126
     2)/CM+WV(1,1)+WV(2,1))                                             LECL-127
    7 AM3=WV(2,1)**.33333333333333D0                                    LECL-128
      BM3=WV(1,1)**.33333333333333D0                                    LECL-129
      WRITE (MW,1008) WV(2,1),SP3,AM3,WV(1,1),SP1,BM3                   LECL-130
      WRITE (MW,1009) E,WV(3,1)                                         LECL-131
      E=0.D0                                                            LECL-132
      GO TO 13                                                          LECL-133
    8 IF (WV(1,IV).EQ.0.D0) GO TO 9                                     LECL-134
      AM3=WV(2,IV)**.33333333333333D0                                   LECL-135
      BM3=WV(1,IV)**.33333333333333D0                                   LECL-136
      WRITE (MW,1008) WV(2,IV),SP3,AM3,WV(1,IV),SP1,BM3                 LECL-137
      GO TO 10                                                          LECL-138
    9 WV(1,IV)=WV(1,IV-1)                                               LECL-139
      WV(2,IV)=WV(2,IV-1)                                               LECL-140
      IPI(2,IV)=IPI(2,IV-1)                                             LECL-141
      IPI(4,IV)=IPI(4,IV-1)                                             LECL-142
      SP2=0.5D0*DFLOAT(IPI(3,IV)-1)                                     LECL-143
      SP3=IPI(4,IV)                                                     LECL-144
   10 WV(3,IV)=WV(3,1)-E                                                LECL-145
      IF (.NOT.(LO(98).AND.LO(8))) GO TO 11                             LECL-146
      EP=E*WV(13,IV)                                                    LECL-147
      ET=E-EP                                                           LECL-148
      WV(1,IV)=WV(1,1)+(WV(2,1)-WV(2,IV))+EP/CM                         LECL-149
      WV(2,IV)=WV(2,1)+(WV(1,1)-WV(1,IV))+ET/CM                         LECL-150
   11 IF (K.EQ.0) K=IPI(5,IV-1)                                         LECL-151
      IF (LO(8)) GO TO 12                                               LECL-152
      WV(13,IV)=WV(3,IV)*(WV(1,IV)+WV(2,IV))/WV(2,IV)                   LECL-153
      GO TO 13                                                          LECL-154
   12 WV(13,IV)=WV(3,IV)*(WV(3,IV)/(2.D0*CM)+WV(1,IV)+WV(2,IV))/WV(2,IV)LECL-155
   13 IF (LO(109).AND.IPI(2,IV).NE.2) GO TO 31                          LECL-156
      IF (K.GT.NPP.OR.K.LE.0) GO TO 32                                  LECL-157
      IPI(5,IV)=K                                                       LECL-158
      IF (IPP(1,K).EQ.-1) IPP(1,K)=IV                                   LECL-159
      NIMAX=MAX0(NIMAX,IPI(2,IV)+IPI(3,IV))                             LECL-160
      IPI(1,IV)=0                                                       LECL-161
      IF (SPI.EQ.SIGM) IPI(1,IV)=1                                      LECL-162
      IF (IV.LE.NCOLL) GO TO 15                                         LECL-163
C UNCOUPLED STATES FOR COMPOUND NUCLEUS.                                LECL-164
      IF (NCONT+IV.GT.NCOLX) GO TO 14                                   LECL-165
      WRITE (MW,1010) IV,SP2,SPI,E,K                                    LECL-166
      GO TO 27                                                          LECL-167
C CONTINUUM FOR COMPOUND NUCLEUS.                                       LECL-168
   14 WRITE (MW,1011) IV,SP2,SPI,E,K                                    LECL-169
      GO TO 27                                                          LECL-170
   15 IF (LO(7)) GO TO 17                                               LECL-171
      NNPA=0                                                            LECL-172
      NNVA=0                                                            LECL-173
      IF ((LO(1).EQV.(N.EQ.0)).OR.(LO(3).AND.(.NOT.LO(1)))) GO TO 16    LECL-174
C INPUT OF PHONONS IN HARMONIC VIBRATIONAL OR SYMMETRIC ROTATIONAL MODELLECL-175
      READ (MR,1012) IPH(1,IV),IPH(2,IV),I2,I3                          LECL-176
      NBET=MAX0(NBET,IPH(2,IV),I2,I3)                                   LECL-177
      M=IPH(1,IV)+1                                                     LECL-178
      GO TO ( 17 , 19 , 20 , 20 ) , M                                   LECL-179
C GROUND STATE.                                                         LECL-180
   16 IPH(1,IV)=0                                                       LECL-181
      IF (LO(1)) GO TO 18                                               LECL-182
   17 WRITE (MW,1013) IV,SP2,SPI,E,K                                    LECL-183
      GO TO 27                                                          LECL-184
   18 WRITE (MW,1014) IV,SP2,SPI,E,K                                    LECL-185
      IF (.NOT.LO(3)) GO TO 27                                          LECL-186
      NNVA=(IPI(3,IV)-1)/4                                              LECL-187
      GO TO 21                                                          LECL-188
C 1 PHONON STATE.                                                       LECL-189
   19 IF (.NOT.LO(1)) WRITE (MW,1015) IV,SP2,SPI,E,K,IPH(2,IV)          LECL-190
      IF (LO(1)) WRITE (MW,1016) IV,SP2,SPI,E,K,IPH(2,IV)               LECL-191
      GO TO 27                                                          LECL-192
   20 NNPA=IPH(1,IV)-1                                                  LECL-193
      NNVA=M-2                                                          LECL-194
      IF (.NOT.LO(3)) NNVA=NNVA+1                                       LECL-195
   21 IF (NA+NVA+NNPA+2.GT.IDT) CALL MEMO('LECL',NA+NVA+NNPA+2,IDT)     LECL-196
      IF ((NVA.EQ.0).OR.(NNPA.EQ.0)) GO TO 23                           LECL-197
      DO 22 I=NVA,1,-1                                                  LECL-198
   22 PA(NNPA+NA+I)=PA(NA+I)                                            LECL-199
   23 IF (LO(3)) GO TO 26                                               LECL-200
      WRITE (MW,1013) IV,SP2,SPI,E,K                                    LECL-201
C 2 PHONON STATE.                                                       LECL-202
      NA=NA+1                                                           LECL-203
      NPA(1,NA)=IPH(2,IV)                                               LECL-204
      IF (LO(1)) GO TO 25                                               LECL-205
      IPH(2,IV)=NA                                                      LECL-206
      NPA(2,NA)=I2                                                      LECL-207
      IF (M.EQ.4) GO TO 24                                              LECL-208
      WRITE (MW,1017) NPA(1,NA),NPA(2,NA)                               LECL-209
      GO TO 27                                                          LECL-210
C MIXTURE OF 1 PHONON AND 2 PHONON STATE - INPUT OF MIXING PARAMETER.   LECL-211
   24 NA=NA+1                                                           LECL-212
      NPA(1,NA)=I3                                                      LECL-213
   25 NVA=NVA+1                                                         LECL-214
      NPA(2,NA)=NVA                                                     LECL-215
      READ (MR,1018) PA(NA+NVA)                                         LECL-216
      B1=0.0174532925199433D0*PA(NA+NVA)                                LECL-217
      C1=DCOS(B1)                                                       LECL-218
      C3=DSIN(B1)                                                       LECL-219
      IF (.NOT.LO(1)) WRITE (MW,1019) PA(NA+NVA),C3,NPA(1,NA-1),I2,C1,I3LECL-220
      IF (LO(1)) WRITE (MW,1020) PA(NA+NVA),C3,NPA(1,NA),C1             LECL-221
      GO TO 27                                                          LECL-222
C ASYMMETRIC ROTATIONAL MODEL - INPUT OF MIXING PARAMETERS.             LECL-223
   26 IF (.NOT.LO(1)) GO TO 27                                          LECL-224
      IPH(1,IV)=NNVA                                                    LECL-225
      IPH(2,IV)=NVA                                                     LECL-226
      NB=NB+NNVA+1                                                      LECL-227
      IF (NNVA.EQ.0) GO TO 27                                           LECL-228
      READ (MR,1018) (PA(NA+NVA+J),J=1,NNVA)                            LECL-229
      WRITE (MW,1021) (PA(NA+NVA+J),J=1,NNVA)                           LECL-230
      NVA=NVA+NNVA                                                      LECL-231
   27 IF (LO(1).AND.LO(2).AND.LO(3)) NVA=MAX0(NVA,5)                    LECL-232
      NPX=0                                                             LECL-233
      DO 28 I=1,NCOLL                                                   LECL-234
      IPI(11,I)=IPI(5,I)                                                LECL-235
      IF (LO(10)) IPI(11,I)=I                                           LECL-236
   28 NPX=MAX0(NPX,IPI(11,I))                                           LECL-237
      RETURN                                                            LECL-238
   29 WRITE (MW,1022)                                                   LECL-239
      GO TO 33                                                          LECL-240
   30 WRITE (MW,1023)                                                   LECL-241
      GO TO 33                                                          LECL-242
   31 SP1=0.5D0*DFLOAT(IPI(2,IV)-1)                                     LECL-243
      WRITE (MW,1024) SP1                                               LECL-244
      GO TO 33                                                          LECL-245
   32 WRITE (MW,1025) K,NPP                                             LECL-246
   33 WRITE (MW,1026)                                                   LECL-247
      STOP                                                              LECL-248
 1000 FORMAT (/' FIRST ORDER VIBRATIONAL MODEL.'/)                      LECL-249
 1001 FORMAT (/' SECOND ORDER VIBRATIONAL MODEL.'/)                     LECL-250
 1002 FORMAT (/' ANHARMONIC VIBRATIONAL MODEL.'/)                       LECL-251
 1003 FORMAT (/' SYMMETRIC ROTATIONAL MODEL.'/)                         LECL-252
 1004 FORMAT (/' ASYMMETRIC ROTATIONAL MODEL.'/)                        LECL-253
 1005 FORMAT (/' EXTERNAL FORM FACTOR MODEL.'/)                         LECL-254
 1006 FORMAT (/' MASSES OF INPUT FOR EXCITED STATES MODIFIED BY "Q" VALULECL-255
     1E.'/)                                                             LECL-256
 1007 FORMAT (F5.2,2I2,A1,6F10.5)                                       LECL-257
 1008 FORMAT (/' TARGET',14X,'MASS =',F10.5,11X,'PRODUCT OF CHARGES =',FLECL-258
     16.0,11X,'AT**1/3 =',1P,D15.6,0P/' INCIDENT PARTICLE',3X,'MASS =',FLECL-259
     210.5,3X,3X,'SPIN =',F4.1,3X,'AP**1/3 =',1P,D15.6)                 LECL-260
 1009 FORMAT (10X,'ENERGY(LAB) =',1P,D15.6,' MEV',10X,'ENERGY(C. M.) =',LECL-261
     1D15.6,' MEV.'/)                                                   LECL-262
 1010 FORMAT (' N =',I3,' -   SPIN =',F4.1,A1,' EXCITATION ENERGY =',F8.LECL-263
     14,' MEV',10X,'POTENTIAL',I5,6X,'***** UNCOUPLED STATE *****')     LECL-264
 1011 FORMAT (' N =',I3,' -   SPIN =',F4.1,A1,' EXCITATION ENERGY =',F8.LECL-265
     14,' MEV',10X,'POTENTIAL',I5,6X,'***** START OF A CONTINUUM *****')LECL-266
 1012 FORMAT (14I5)                                                     LECL-267
 1013 FORMAT (' N =',I3,' -   SPIN =',F4.1,A1,' EXCITATION ENERGY =',F8.LECL-268
     14,' MEV',10X,'POTENTIAL',I5)                                      LECL-269
 1014 FORMAT (' N =',I3,' -   SPIN =',F4.1,A1,' EXCITATION ENERGY =',F8.LECL-270
     14,' MEV',10X,'POTENTIAL',I5,6X,'GROUND STATE BAND.')              LECL-271
 1015 FORMAT (' N =',I3,' -   SPIN =',F4.1,A1,' EXCITATION ENERGY =',F8.LECL-272
     14,' MEV',10X,'POTENTIAL',I5,6X,' PHONON STATE WITH PHONON',I3)    LECL-273
 1016 FORMAT (' N =',I3,' -   SPIN =',F4.1,A1,' EXCITATION ENERGY =',F8.LECL-274
     14,' MEV',10X,'POTENTIAL',I5,6X,'VIBRATIONAL BAND OF PHONON',I3)   LECL-275
 1017 FORMAT (23X,'2 PHONONS STATE, WITH PHONONS',I3,' AND',I3)         LECL-276
 1018 FORMAT (7F10.5)                                                   LECL-277
 1019 FORMAT (23X,'MIXING OF 1 AND 2 PHONON STATES WITH',F9.3,' DEGREES'LECL-278
     1/23X,F10.5,' 2 PHONONS STATE, WITH PHONONS',I3,' AND',I3,'  + ',F1LECL-279
     20.5,' 1 PHONON STATE WITH PHONON',I3)                             LECL-280
 1020 FORMAT (23X,'MIXING OF VIBRATIONAL AND GROUND BANDS WITH',F9.3,' DLECL-281
     1EGREES'/23X,F10.5,' VIBRATIONAL BAND OF PHONON',I3,' AND',F10.5,' LECL-282
     2GROUND STATE BAND.')                                              LECL-283
 1021 FORMAT (23X,'BAND MIXING COEFF.',5F10.5)                          LECL-284
 1022 FORMAT (' THE NUCLEAR REDUCED MATRIX ELEMENTS MUST BE READ IN THISLECL-285
     1 MODEL.')                                                         LECL-286
 1023 FORMAT (' PROJECTILE-TARGET ANTISYMMETRISATION VALID ONLY FOR SPINLECL-287
     1 OF PARTICLE EQUAL TO SPIN OF TARGET AND POSITIVE PARITY'/49X,'ANDLECL-288
     2 SCHROEDINGER FORMALISM.')                                        LECL-289
 1024 FORMAT (' PARTICLE SPIN',F5.1,' NOT ALLOWED FOR DIRAC EQUATION.') LECL-290
 1025 FORMAT (' POTENTIAL',I3,' WILL NOT BE READ.TOTAL NUMBER IS:',I3)  LECL-291
 1026 FORMAT (//' IN LECL  ...  STOP  ...')                             LECL-292
      END                                                               LECL-293
