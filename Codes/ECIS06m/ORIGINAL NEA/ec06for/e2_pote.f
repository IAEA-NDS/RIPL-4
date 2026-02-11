C 07/03/07                                                      ECIS06  POTE-000
      SUBROUTINE POTE(BETA,NBETA,IVQ,IVY,IVZ,VAL,NVAL,WV,WVM,ISM,NCOLL,NPOTE-001
     1COLD,NCOLT,IDT,IPI,IPIM,IPP,MCM,LO,V,VCO,VDO,TL,NP,Q)             POTE-002
C COMPUTATION AND OUTPUT OF FORM FACTORS.                               POTE-003
C INPUT:     BETA:    DEFORMATIONS OF THE TARGET.                       POTE-004
C            NBETA:   EQUIVALENT BY CALL WITH BETA.                     POTE-005
C            IVQ:     TABLE OF ANGULAR MOMENTA OF FORM FACTORS.         POTE-006
C            IVY:     TABLE FOR COMPUTATION OF FORM FACTORS.            POTE-007
C            IVZ:     TABLE FOR USE OF FORM FACTORS.                    POTE-008
C            VAL,NVAL:OPTICAL MODEL AND FOLDING PARAMETERS.             POTE-009
C            WV,WVM:  DESCRIPTION OF THE CHANNELS (SEE CALX).           POTE-010
C            ISM:     NUMBER OF POINTS.                                 POTE-011
C            NCOLL:   NUMBER OF COUPLED NUCLEAR STATES.                 POTE-012
C            NCOLD:   TOTAL NUMBER OF STATES WITHOUT CONTINUUM.         POTE-013
C            NCOLT:   TOTAL NUMBER OF STATES WITH CONTINUUM.            POTE-014
C            IDT:     MAXIMUM WORKING FIELD Q.                          POTE-015
C            IPI,IPIM:DESCRIPTION OF THE CHANNELS (SEE CALX).           POTE-016
C            IPP:     CROSS-REFERENCE TO POTENTIALS IN IPI(1,*).        POTE-017
C            MCM:     ANGULAR MOMENTUM LIMITS ON COULOMB CORRECTIONS.   POTE-018
C            LO:      LOGICAL CONTROLS:                                 POTE-019
C               LO(1)  =.TRUE. ROTATIONAL MODEL-(.F.:VIBRATIONAL MODEL).POTE-020
C               LO(3)  =.TRUE. ANHARMONIC VIBRATIONAL OR ASYMMETRIC     POTE-021
C                              ROTATIONAL MODEL.                        POTE-022
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    POTE-023
C               LO(10) =.TRUE. DISPERSION RELATIONS FOR POTENTIALS.     POTE-024
C               LO(11) =.TRUE. DEFORMED COULOMB POTENTIAL.              POTE-025
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            POTE-026
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      POTE-027
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. POTE-028
C               LO(17) =.TRUE. FOLDING MODEL.                           POTE-029
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   POTE-030
C               LO(51) =.TRUE. OUTPUT OF POTENTIALS.                    POTE-031
C               LO(62) =.TRUE. POTENTIALS WRITTEN ON FILE 62.           POTE-032
C               LO(74) =.TRUE. OUTPUT OF TIME IN DIFFERENT STEPS.       POTE-033
C               LO(93) =.TRUE. NO RECOIL CORRECTION FOR REACTIONS.      POTE-034
C               LO(93) =.TRUE. NO RECOIL CORRECTION FOR REACTIONS WITH  POTE-035
C                              SMALL DIFFERENCES (LESS THAN .5) OF      POTE-036
C                              TARGET MASSES.                           POTE-037
C               LO(99) =.TRUE. SCHROEDINGER EQUIVALENT TO DIRAC         POTE-038
C                              EQUATION.                                POTE-039
C               LO(100)=.TRUE. DIRAC EQUATION.                          POTE-040
C               LO(101)=.TRUE. THERE IS A REAL SPIN-ORBIT POTENTIAL.    POTE-041
C               LO(102)=.TRUE. THERE IS AN IMAGINARY SPIN-ORBIT         POTE-042
C                              POTENTIAL.                               POTE-043
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. POTE-044
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    POTE-045
C               LO(117)=.TRUE. FOR ALL CALCULATIONS EXCEPT THE FIRST.   POTE-046
C               LO(120)=.TRUE. OUTPUT AND LAST CALCULATION BEST ONE.    POTE-047
C               LO(129)=.TRUE. REAL SPIN-ORBIT OR DIRAC EQUATION.       POTE-048
C               LO(130)=.TRUE. IMAGINARY SPIN-ORBIT OR DIRAC EQUATION.  POTE-049
C               LO(133)=.TRUE. STORE SCALAR AND COULOMB POTENTIAL       POTE-050
C                              INDEPENDENTLY.                           POTE-051
C OUTPUT:    TL:      TRANSMISSION COEFFICIENTS FOR UNCOUPLED STATES.   POTE-052
C            V(I,J):  POTENTIAL AND TRANSITION FORM FACTORS.            POTE-053
C            VCO(2,I):STRENGTH OF COULOMB POTENTIAL TAILS, SCALAR IN    POTE-054
C                     VCO(1,I), SPIN-ORBIT IN VCO(2,I).                 POTE-055
C            VDO(2,I):STRENGTH OF COULOMB TRANSITION TAILS, SCALAR IN   POTE-056
C                     VDO(1,I), SPIN-ORBIT IN VDO(2,I).                 POTE-057
C            IDT:     RETURNS THE PLACE USED.                           POTE-058
C WORKING AREA:                                                         POTE-059
C            NP,Q:    IN EQUIVALENCE BY CALL.                           POTE-060
C                                                                       POTE-061
C***********************************************************************POTE-062
C     DISPOSITION OF FORM-FACTORS IN THE AREA V(I,*):                   POTE-063
C                                                                       POTE-064
C                     FOR SCHROEDINGER EQUATIONS:                       POTE-065
C     AFTER ITY(1) (=0) REAL CENTRAL POTENTIALS IN ORDER IPI(11,*)      POTE-066
C     AFTER ITY(3) REAL SPIN-ORBIT FORM FACTORS IN THE SAME ORDER.      POTE-067
C     AFTER ITY(5) REAL TRANSITION FORM FACTORS IN ORDER IVZ(1,*).      POTE-068
C     AFTER ITY(7) REAL DERIVED SPIN-ORBIT FORM FACTORS IN ORDER        POTE-069
C           IVZ(1,*) WITH THE R**-2 FORM FACTOR AT THE SAME ADDRESS+INLSPOTE-070
C     AFTER ITY(9) COULOMB CENTRAL POTENTIAL IF NEEDED SEPARATED IN THE POTE-071
C           ORDER OF REAL POTENTIALS.                                   POTE-072
C     AFTER ITY(10) COULOMB SPIN-ORBIT POTENTIAL IF NEEDED SEPARATED IN POTE-073
C           THE ORDER OF SPIN-ORBIT POTENTIALS.                         POTE-074
C     IMAGINARY POTENTIALS AND FORM FACTOR FOLLOW WITH THE SAME ORDER   POTE-075
C           FROM ITY(2), ITY(4), ITY(6) AND ITY(8)                      POTE-076
C                                                                       POTE-077
C                     FOR DIRAC EQUATIONS:                              POTE-078
C     AFTER ITY(1) (=0) FOR EACH LEVEL, 14 FORM FACTORS V(*,I,*) WHICH  POTE-079
C           ARE CENTRAL POTENTIAL FOR I=1,2   SPIN-ORBIT POTENTIAL FOR  POTE-080
C           I=3,4  D=E+M+V-W FOR I=5,6   E-M-V-W FOR I=7,8  D**(1/2) FORPOTE-081
C           I=9,10 D**(-1) FOR I=11,12  AND  - D( TENSOR ) FOR I=13,14  POTE-082
C     AFTER ITX(6) THE INVT SETS OF COMPLEX SCALAR + VECTOR FORM FACTORSPOTE-083
C           FOLLOWED BY THE INSL SETS OF TENSOR FORM FACTORS            POTE-084
C***********************************************************************POTE-085
C                                                                       POTE-086
C FOR THE COMMONS /COUPL/ AND /DCONS/ SEE CALC.                         POTE-087
C FOR THE COMMONS /POTE1/ AND /POTE2/ SEE REDM.                         POTE-088
C                                                                       POTE-089
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /DCONS/:                     POTE-090
C  CHB:       PLANCK CONSTANT /(2*PI) IN MEV*FERMI.                     POTE-091
C  CCZ:       COULOMB ALPHA CONSTANT.                                   POTE-092
C   USED:     CHB,CCZ.                                                  POTE-093
C                                                                       POTE-094
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     POTE-095
C  IQM:       MAXIMUM L-VALUE OF DEFORMATION IN ROTATIONAL MODEL.       POTE-096
C  NBT1:      NUMBER OF PHONONS (VIBRATIONS).                           POTE-097
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             POTE-098
C   USED:     IQM,NBT1,NPP.                                             POTE-099
C                                                                       POTE-100
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE1/:                     POTE-101
C  ITX(16):   STARTING ADDRESS OF DIFFERENT FORM FACTORS.               POTE-102
C             FOR SCHROEDINGER EQUATION, ITX(I)+1 IS THE STARTING       POTE-103
C             ADDRESS OF THE FORM FACTOR READ IN EXTP WITH ITYP=I       POTE-104
C             (POTENTIAL FOR I=1 TO 8, TRANSITION FOR I=9 TO 16).       POTE-105
C             FOR DIRAC EQUATIONS, ITX(1)=0,                            POTE-106
C             ITX(2)+1=ADDRESS OFF FIRST TRANSITION FORM FACTOR,        POTE-107
C             ITX(7)=ADDRESS OF LAST TRANSITION FORM FACTOR,            POTE-108
C             ITX(3)=ADDRESS OF LAST TEMPORARY CENTRAL POTENTIAL,       POTE-109
C             ITX(4)=ITX(7)-24,ITX(5)=ITX(3)-11,ITX(6)=ITX(2)-4.        POTE-110
C             ALL ARE USED FOR SCHROEDINGER, THE FIRST 8 FOR DIRAC.     POTE-111
C  IMAX:      MAXIMUM ANGULAR MOMENTUM.                                 POTE-112
C  INTC:      NUMBER OF FORM FACTORS WITHOUT DEFORMED SPIN-ORBIT        POTE-113
C             INCLUDING CORRECTION TERMS.                               POTE-114
C  INLS:      NUMBER OF SPIN-ORBIT FORM FACTORS NOT TAKING INTO ACCOUNT POTE-115
C             MULTIPLICATION BY 2.                                      POTE-116
C  INVD:      IDEM FOR COULOMB SPIN-ORBIT.                              POTE-117
C  ITXM:      TOTAL NUMBER OF FORM FACTORS.                             POTE-118
C   USED:     ITX,IMAX,INTC,INLS,INVD,ITXM.                             POTE-119
C                                                                       POTE-120
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     POTE-121
C  ITY(1):    STARTING ADDRESS OF REAL CENTRAL POTENTIAL (IT IS 0).     POTE-122
C  ITY(2):    STARTING ADDRESS OF IMAGINARY CENTRAL POTENTIAL.          POTE-123
C  ITY(3):    STARTING ADDRESS OF REAL SPIN-ORBIT POTENTIAL.            POTE-124
C  ITY(4):    STARTING ADDRESS OF IMAGINARY SPIN-ORBIT POTENTIAL.       POTE-125
C  ITY(5):    STARTING ADDRESS OF REAL CENTRAL TRANSITION.              POTE-126
C  ITY(6):    STARTING ADDRESS OF IMAGINARY CENTRAL TRANSITION.         POTE-127
C  ITY(7):    STARTING ADDRESS OF REAL SPIN-ORBIT TRANSITION.           POTE-128
C  ITY(8):    STARTING ADDRESS OF IMAGINARY SPIN-ORBIT TRANSITION.      POTE-129
C  ITY(9):    STARTING ADDRESS OF COULOMB CENTRAL POTENTIAL.            POTE-130
C  ITY(10):   STARTING ADDRESS OF COULOMB SPIN-ORBIT POTENTIAL.         POTE-131
C  ITY(11):   STARTING ADDRESS OF COULOMB CENTRAL TRANSITION.           POTE-132
C  ITY(12):   STARTING ADDRESS OF COULOMB SPIN-ORBIT TRANSITION.        POTE-133
C        ITY(2)=14*NCOLL AND ITY(5)=0 ONLY ARE USED FOR DIRAC EQUATIONS.POTE-134
C  INVT:      NUMBER OF TRANSITION FORM FACTORS WITHOUT SPIN-ORBIT.     POTE-135
C  INTV:      SAME AS INVT, TAKING INTO ACCOUNT DISPERSION.             POTE-136
C  INSL:      NUMBER OF SPIN-ORBIT FORM TRANSITION FACTORS NOT TAKING   POTE-137
C             INTO ACCOUNT MULTIPLICATION BY 2.                         POTE-138
C  NPX:       NUMBER OF POTENTIALS TAKING INTO ACCOUNT DISPERSION.      POTE-139
C   USED:     ITY,INVT,INTV,INSL,NPX.                                   POTE-140
C                                                                       POTE-141
C***********************************************************************POTE-142
      IMPLICIT REAL*8 (A-H,O-Z)                                         POTE-143
      LOGICAL LO(150)                                                   POTE-144
      DIMENSION BETA(9,*),NBETA(18,*),V(ISM,*),IVQ(3,*),IVY(7,*),IVZ(4,*POTE-145
     1),NP(2,*),ILO(16),Q(*),VAL(42,*),NVAL(*),VCO(2,*),VDO(2,*),WV(22,*POTE-146
     2),WVM(22,*),IPI(11,*),IPIM(11,*),IPP(34,*),TL(*),MCM(*),IPIX(11),PPOTE-147
     3GN(10),XGN(10),W(24),WVX(22)                                      POTE-148
      CHARACTER*4 BB(2)                                                 POTE-149
      CHARACTER*8 AA(3,8)                                               POTE-150
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   POTE-151
      COMMON /DCONS/ CM,CHB,CZ,CMB,CCZ,CK,XZ                            POTE-152
      COMMON /POTE1/ ITX(16),IMAX,INTC,INLS,INVC,INVD,ITXM              POTE-153
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         POTE-154
      COMMON /INOUT/ MR,MW,MS                                           POTE-155
      DATA PGN( 1),PGN( 2) / 1.52753387130726D-01, 1.49172986472604D-01/POTE-156
      DATA PGN( 3),PGN( 4) / 1.42096109318382D-01, 1.31688638449177D-01/POTE-157
      DATA PGN( 5),PGN( 6) / 1.18194531961518D-01, 1.01930119817240D-01/POTE-158
      DATA PGN( 7),PGN( 8) / 8.32767415767047D-02, 6.26720483341091D-02/POTE-159
      DATA PGN( 9),PGN(10) / 4.06014298003869D-02, 1.76140071391521D-02/POTE-160
      DATA XGN( 1),XGN( 2) / 7.65265211334973D-02, 2.27785851141645D-01/POTE-161
      DATA XGN( 3),XGN( 4) / 3.73706088715420D-01, 5.10867001950827D-01/POTE-162
      DATA XGN( 5),XGN( 6) / 6.36053680726515D-01, 7.46331906460151D-01/POTE-163
      DATA XGN( 7),XGN( 8) / 8.39116971822219D-01, 9.12234428251326D-01/POTE-164
      DATA XGN( 9),XGN(10) / 9.63971927277914D-01, 9.93128599185095D-01/POTE-165
      DATA BB /'    ','LAST'/                                           POTE-166
      DATA AA/'      RE','AL VOLUM','E/SCALAR',' IMAGINA','RY VOLUM','E/POTE-167
     1SCALAR','     REA','L SURFAC','E/VECTOR','  IMAGIN','. SURFAC','E/POTE-168
     2VECTOR','  REAL S','PIN-ORBI','T/TENSOR',' IMAG. S','PIN-ORBI','T/POTE-169
     3TENSOR','        ','        ',' COULOMB','      SP','IN-ORBIT',' CPOTE-170
     4OULOMB'/                                                          POTE-171
      DATA ILO /62,62,62,62,101,102,62,103,62,12,62,12,13,14,11,19/     POTE-172
      IDY=0                                                             POTE-173
      IF (LO(120)) GO TO 4                                              POTE-174
      ITM=ITXM                                                          POTE-175
      IF (LO(100)) ITM=ITM-ITX(7)                                       POTE-176
      IF (LO(7)) GO TO 3                                                POTE-177
      INVZ=INTC                                                         POTE-178
      INY=1                                                             POTE-179
      IF (LO(109)) INY=3                                                POTE-180
      DO 2 IP=1,NPP                                                     POTE-181
      K=IABS(IPP(1,IP))                                                 POTE-182
      INVY=INVZ+INY                                                     POTE-183
      N1=(INVY+1)/2+1                                                   POTE-184
      N2=N1+MAX0(IMAX,IQM)+1                                            POTE-185
      ID1=INVY+27                                                       POTE-186
      IDX=10                                                            POTE-187
      IF (LO(3)) IDX=36                                                 POTE-188
      IF (.NOT.LO(1)) IDX=1                                             POTE-189
      IDX=MAX0(IDX*ID1+10*INVY,6*ISM)+N2                                POTE-190
      IF (IDX.GT.IDT) CALL MEMO('POTE',IDT,IDX)                         POTE-191
      IK=IP                                                             POTE-192
      IF (LO(100)) IK=ITX(4)+24*IP+1                                    POTE-193
      CALL ROTP(BETA,NBETA,IK,IVY,VAL(1,IP),ID1,ISM,WV(1,K),INVZ,INVY,INPOTE-194
     1Y,PGN,XGN,LO,V,V(1,ITX(3)+1),NP,Q(N1),Q(N2),Q(10*INVY+N2))        POTE-195
C  TO USE THE FOLDING MODEL, THE SPIN-ORBIT POTENTIALS AND THE FIRST    POTE-196
C SPIN-ORBIT NON DIAGONAL FORM FACTORS ARE NOT DERIVED ( V(R) INSTEAD OFPOTE-197
C (1/R)*(D/DR)(V(R)) )  FOR SCHROEDINGER EQUATIONS. FOR DIRAC EQUATIONS POTE-198
C THE DERIVATIVES ARE COMPUTED BUT NOT USED.                            POTE-199
      IDY=MAX0(IDX,IDY)                                                 POTE-200
      IF (LO(74)) CALL HORA                                             POTE-201
      IF (.NOT.LO(17)) GO TO 2                                          POTE-202
C  *******  FOLDING MODEL  **********                                   POTE-203
      IST=5+ISM+IDINT(2.D0*(DABS(VAL(41,IP))+2.D0*DABS(VAL(42,IP)))/WV(8POTE-204
     1,K))                                                              POTE-205
      IF (VAL(40,IP)*VAL(41,IP).EQ.0.D0) IST=ISM+5                      POTE-206
C THE WORKING SPACE IN P IS SHIFTED FROM V WITH THE ORIGINS IN TABLE IT.POTE-207
      NNF=1+ISM*ITM                                                     POTE-208
      IDZ=NNF+4*ITM                                                     POTE-209
      IDX=IDZ+4*IST                                                     POTE-210
      IF (INVZ.NE.0) IDX=IDX+2*IST*IMAX                                 POTE-211
      IDY=MAX0(IDY,IDX)                                                 POTE-212
      IF (IDX.GT.IDT) CALL MEMO('POTE',IDT,IDX)                         POTE-213
      DO 1 J=1,NNF                                                      POTE-214
    1 Q(J)=0.D0                                                         POTE-215
      CALL FOLD(V(1,IK),Q,VAL(34,IP),3,IP,ISM,IST,IVY,INVZ,Q(IDZ),PGN,XGPOTE-216
     1N,WV,NP(1,NNF),LO)                                                POTE-217
      IF (LO(74)) CALL HORA                                             POTE-218
    2 INVZ=0                                                            POTE-219
      GO TO 4                                                           POTE-220
    3 NNF=2*ITM                                                         POTE-221
      IF (NNF.GT.IDT) CALL MEMO('POTE',IDT,NNF)                         POTE-222
      IDZ=IDT-NNF                                                       POTE-223
      IK=1                                                              POTE-224
      IF (LO(100)) IK=ITX(7)+1                                          POTE-225
      CALL STDP(V(1,IK),IVY,ISM,VAL,NVAL,IDZ,IDX,Q(NNF+1),WV,PGN,XGN,NPPPOTE-226
     1,NP,Q(NNF+1),LO)                                                  POTE-227
      IDY=MAX0(IDY,IDX+2*NNF)                                           POTE-228
    4 IF (.NOT.LO(62)) GO TO 14                                         POTE-229
      NPT=8*(NPP+INTC)                                                  POTE-230
      IF (.NOT.LO(101)) NPT=NPT-NPP                                     POTE-231
      IF (.NOT.LO(102)) NPT=NPT-NPP                                     POTE-232
      IF (.NOT.LO(103)) NPT=NPT-NPP                                     POTE-233
      IF (.NOT.LO(11)) NPT=NPT-INTC                                     POTE-234
      IF (.NOT.LO(12)) NPT=NPT-2*INTC                                   POTE-235
      IF (.NOT.LO(13)) NPT=NPT-INTC                                     POTE-236
      IF (.NOT.LO(14)) NPT=NPT-INTC                                     POTE-237
      IF (.NOT.LO(19)) NPT=NPT-INTC                                     POTE-238
      WRITE (62,1000) WV(1,1),WV(13,1),WV(2,1),IPI(4,1),NPT             POTE-239
      IST=2*(ISM/2)                                                     POTE-240
      IDX=IST                                                           POTE-241
      IDY=MAX0(IDY,IDX)                                                 POTE-242
      IF (IDX.GT.IDT) CALL MEMO('POTE',IDT,IDX)                         POTE-243
C PUNCH THE FORM FACTORS.                                               POTE-244
      DO 5 I=1,IST                                                      POTE-245
    5 Q(I)=I*WV(8,1)                                                    POTE-246
      IR=0                                                              POTE-247
      DO 8 IP=1,NPP                                                     POTE-248
      DO 7 IJ=1,8                                                       POTE-249
      IS=1                                                              POTE-250
      IL=ILO(IJ)                                                        POTE-251
      IF (.NOT.LO(IL)) GO TO 7                                          POTE-252
      JI=IP+ITX(IJ)                                                     POTE-253
      IF (LO(100)) JI=ITX(7)+IJ+24*(IP-1)                               POTE-254
      IR=IR+1                                                           POTE-255
      WRITE (MW,1001) (AA(J,IJ),J=1,3),IP,IR                            POTE-256
      WRITE (62,1001) (AA(J,IJ),J=1,3),IP,IR                            POTE-257
      DO 6 I=1,IST,2                                                    POTE-258
      IF (I.EQ.IST-1) IS=2                                              POTE-259
      WRITE (MW,1002) Q(I),V(I,JI),Q(I+1),V(I+1,JI),BB(IS),IR           POTE-260
    6 WRITE (62,1003) Q(I),V(I,JI),Q(I+1),V(I+1,JI),BB(IS),IR           POTE-261
    7 CONTINUE                                                          POTE-262
    8 CONTINUE                                                          POTE-263
      IF (INTC.EQ.0) GO TO 14                                           POTE-264
      DO 13 IQ=1,INTC                                                   POTE-265
      DO 12 IJ=9,16                                                     POTE-266
      IS=1                                                              POTE-267
      IL=ILO(IJ)                                                        POTE-268
      IF (.NOT.LO(IL)) GO TO 12                                         POTE-269
      IP=IQ                                                             POTE-270
      IF (IJ.EQ.13.OR.IJ.EQ.14) IP=IVY(3,IQ)+INLS                       POTE-271
      IF (IJ.GE.15) IP=IVY(IJ-11,IQ)                                    POTE-272
      IF (IP.LE.0) GO TO 12                                             POTE-273
      JI=IP+ITX(IJ)                                                     POTE-274
      IF (LO(100)) JI=ITX(7)+24*NPP+11*IQ+IJ-19                         POTE-275
      IR=IR+1                                                           POTE-276
      WRITE (MW,1004) (AA(J,IJ-8),J=1,3),IQ,IR                          POTE-277
      WRITE (62,1004) (AA(J,IJ-8),J=1,3),IQ,IR                          POTE-278
      IF (.NOT.LO(100).AND.(IJ.GE.13).AND.(IJ.NE.15)) GO TO 10          POTE-279
      DO 9 I=1,IST,2                                                    POTE-280
      IF (I.EQ.IST-1) IS=2                                              POTE-281
      WRITE (MW,1002) Q(I),V(I,JI),Q(I+1),V(I+1,JI),BB(IS),IR           POTE-282
    9 WRITE (62,1003) Q(I),V(I,JI),Q(I+1),V(I+1,JI),BB(IS),IR           POTE-283
      GO TO 12                                                          POTE-284
   10 IF (IJ.EQ.16) JI=JI+INVD                                          POTE-285
      DO 11 I=1,IST,2                                                   POTE-286
      A1=-V(I,JI)*Q(I)**2                                               POTE-287
      A2=-V(I+1,JI)*Q(I+1)**2                                           POTE-288
      IF (I.EQ.IST-1) IS=2                                              POTE-289
      WRITE (MW,1002) Q(I),A1,Q(I+1),A2,BB(IS),IR                       POTE-290
   11 WRITE (62,1003) Q(I),A1,Q(I+1),A2,BB(IS),IR                       POTE-291
   12 CONTINUE                                                          POTE-292
   13 CONTINUE                                                          POTE-293
   14 IF (LO(120)) GO TO 68                                             POTE-294
      IF (INTC.EQ.INVT) GO TO 26                                        POTE-295
C ADDITION OF CORRECTION TERMS OF FORM FACTORS.                         POTE-296
      DO 25 K=1,INVT                                                    POTE-297
      IF (IVY(6,K).LE.0) GO TO 25                                       POTE-298
      K1=IVY(6,K)                                                       POTE-299
      L=IVY(7,K)+2                                                      POTE-300
      DO 24 N=1,8                                                       POTE-301
      IL=ILO(N+8)                                                       POTE-302
      IF ((IL.NE.62).AND.(.NOT.LO(IL))) GO TO 24                        POTE-303
      L1=L                                                              POTE-304
      M=0                                                               POTE-305
      IF (LO(100)) GO TO 17                                             POTE-306
      IF (N.GT.6) GO TO 16                                              POTE-307
      IF (N.GT.4) GO TO 15                                              POTE-308
      K2=ITX(N+8)+K                                                     POTE-309
      K3=ITX(N+8)+K1                                                    POTE-310
      GO TO 19                                                          POTE-311
   15 IF (IVY(3,K).EQ.0) GO TO 24                                       POTE-312
      K2=ITX(N+8)+IVY(3,K)+INLS                                         POTE-313
      K3=ITX(N+8)+IVY(3,K1)+INLS                                        POTE-314
      L1=L1+2                                                           POTE-315
      M=-INLS                                                           POTE-316
      GO TO 19                                                          POTE-317
   16 IF (IVY(N-3,K).EQ.0) GO TO 24                                     POTE-318
      K2=ITX(N+8)+IVY(N-3,K)                                            POTE-319
      K3=ITX(N+8)+IVY(N-3,K1)                                           POTE-320
      IF (N.EQ.8) M=INVD                                                POTE-321
      GO TO 18                                                          POTE-322
   17 IF ((N.GT.4).AND.(N.NE.7)) M=MIN0(11-N,4)                         POTE-323
      K2=ITX(5)+11*K+N                                                  POTE-324
      K3=K2+11*(K1-K)                                                   POTE-325
      IF (N.LE.6) GO TO 19                                              POTE-326
   18 CI=V(ISM,K3)                                                      POTE-327
      CR=V(ISM,K2)                                                      POTE-328
      GO TO 21                                                          POTE-329
   19 DD=0.D0                                                           POTE-330
      CI=0.D0                                                           POTE-331
      CR=0.D0                                                           POTE-332
      DO 20 IS=1,ISM                                                    POTE-333
      DD=DD+WV(8,1)                                                     POTE-334
      CI=CI+V(IS,K3)*DD**L1                                             POTE-335
   20 CR=CR+V(IS,K2)*DD**L1                                             POTE-336
   21 IF ((CI.EQ.0.D0).OR.(CR.EQ.0.D0)) GO TO 24                        POTE-337
      CR=CR/CI                                                          POTE-338
      IF (.NOT.LO(117)) WRITE (MW,1005) K,N,CR                          POTE-339
   22 DO 23 IS=1,ISM                                                    POTE-340
   23 V(IS,K2)=V(IS,K2)-CR*V(IS,K3)                                     POTE-341
      IF (M.EQ.0) GO TO 24                                              POTE-342
      K2=K2+M                                                           POTE-343
      K3=K3+M                                                           POTE-344
      M=0                                                               POTE-345
      GO TO 22                                                          POTE-346
   24 CONTINUE                                                          POTE-347
   25 CONTINUE                                                          POTE-348
   26 RM=ISM*WV(8,1)                                                    POTE-349
      IF (LO(74)) CALL HORA                                             POTE-350
C STORAGE OF OPTICAL POTENTIALS.                                        POTE-351
      DO 44 KJ=1,NCOLT                                                  POTE-352
      JJ=MOD(KJ+NCOLL-1,NCOLT)+1                                        POTE-353
      IF (JJ.GT.NCOLD) GO TO 29                                         POTE-354
      DO 27 I=1,11                                                      POTE-355
   27 IPIX(I)=IPI(I,JJ)                                                 POTE-356
      DO 28 I=1,22                                                      POTE-357
   28 WVX(I)=WV(I,JJ)                                                   POTE-358
      GO TO 32                                                          POTE-359
   29 DO 30 I=1,11                                                      POTE-360
   30 IPIX(I)=IPIM(I,JJ-NCOLD)                                          POTE-361
      DO 31 I=1,22                                                      POTE-362
   31 WVX(I)=WVM(I,JJ-NCOLD)                                            POTE-363
   32 XP=CHB/(2.D0*WVX(6))                                              POTE-364
      XQ=XP**2                                                          POTE-365
      IF (LO(100)) GO TO 38                                             POTE-366
C POTENTIALS FOR SCHROEDINGER EQUATION.                                 POTE-367
      IF (JJ.GT.NCOLL) GO TO 33                                         POTE-368
      I1=IPIX(11)                                                       POTE-369
      I2=I1+ITY(2)                                                      POTE-370
      I3=I1+ITY(3)                                                      POTE-371
      I4=I1+ITY(4)                                                      POTE-372
      I5=I1+ITY(9)                                                      POTE-373
      I6=I1+ITY(10)                                                     POTE-374
      IF (LO(133)) GO TO 34                                             POTE-375
      GO TO 36                                                          POTE-376
   33 I1=1                                                              POTE-377
      I2=2                                                              POTE-378
      I3=3                                                              POTE-379
      I4=4                                                              POTE-380
      I5=5                                                              POTE-381
      I6=6                                                              POTE-382
      IF (LO(133)) GO TO 34                                             POTE-383
      I5=1                                                              POTE-384
      I6=3                                                              POTE-385
      GO TO 36                                                          POTE-386
   34 DO 35 IS=1,ISM                                                    POTE-387
      IF (LO(103)) V(IS,I6)=0.D0                                        POTE-388
   35 V(IS,I5)=0.D0                                                     POTE-389
   36 J=IPIX(5)                                                         POTE-390
      DO 37 IS=1,ISM                                                    POTE-391
      V(IS,I1)=V(IS,J+ITX(1))*(1.D0+WVX(16))+V(IS,J+ITX(3))*(1.D0+WVX(19POTE-392
     1))+WVX(15)*V(IS,J+ITX(2))+WVX(18)*V(IS,J+ITX(4))                  POTE-393
      V(IS,I2)=V(IS,J+ITX(2))*(1.D0+WVX(14))+V(IS,J+ITX(4))*(1.D0+WVX(17POTE-394
     1))                                                                POTE-395
      IF (LO(129)) V(IS,I3)=V(IS,J+ITX(5))*(1.D0+WVX(22))+WVX(21)*V(IS,JPOTE-396
     1+ITX(6))                                                          POTE-397
      IF (LO(130)) V(IS,I4)=V(IS,J+ITX(6))*(1.D0+WVX(20))               POTE-398
      IF (LO(103)) V(IS,I6)=V(IS,I6)+V(IS,J+ITX(8))*XQ                  POTE-399
   37 V(IS,I5)=V(IS,I5)-V(IS,J+ITX(7))                                  POTE-400
      IF (LO(103)) VCO(2,I1)=V(ISM,J+ITX(8))*RM**3/WVX(6)               POTE-401
      IF (.NOT.LO(99)) GO TO 43                                         POTE-402
      VCO(1,I1)=RM*V(ISM,J+ITX(7))*CHB/(2.D0*WVX(7)*WVX(6))             POTE-403
      VCO(2,I1)=VCO(2,I1)+VCO(1,I1)*CHB/(WVX(6)+WVX(7))                 POTE-404
      GO TO 43                                                          POTE-405
C POTENTIALS FOR DIRAC EQUATION:                                        POTE-406
C AT THIS POINT:   AFTER V(1,ITX(7)),                                   POTE-407
C V(*,N,1) SCALAR POTENTIAL V FOR N=1,2 VECTOR POTENTIAL W FOR N=3,4    POTE-408
C          TENSOR POTENTIAL T FOR N=5,6 COULOMB POTENTIALS FOR N=7,8    POTE-409
C          THEIR FIRST DERIVATIVES FOR N=9,16                           POTE-410
C          THEIR SECOND DERIVATIVES FOR N=17,24                         POTE-411
C THE FIRST DERIVATIVES OF NON-COULOMB POTENTIALS HAVE A REVERSED SIGN. POTE-412
   38 M=ITX(4)+24*IPIX(5)                                               POTE-413
      J=JJ                                                              POTE-414
      IF (J.GT.NCOLL) J=1                                               POTE-415
      N=14*J-14                                                         POTE-416
      RR=0.D0                                                           POTE-417
      DO 42 IS=1,ISM                                                    POTE-418
      RR=RR+WVX(8)                                                      POTE-419
      DO 39 K=1,24                                                      POTE-420
   39 W(K)=V(IS,K+M)                                                    POTE-421
      IF (.NOT.LO(10)) GO TO 41                                         POTE-422
      DO 40 K=1,24,8                                                    POTE-423
      W(K)=(1.D0+WVX(16))*W(K)+WVX(15)*W(K+1)                           POTE-424
      W(K+1)=W(K+1)+WVX(14)*W(K+1)                                      POTE-425
      W(K+2)=(1.D0+WVX(19))*W(K+2)+WVX(18)*W(K+3)                       POTE-426
      W(K+3)=W(K+3)+WVX(17)*W(K+3)                                      POTE-427
      W(K+4)=(1.D0+WVX(22))*W(K+4)+WVX(21)*W(K+5)                       POTE-428
   40 W(K+5)=(1.D0+WVX(20))*W(K+5)                                      POTE-429
   41 W(3)=W(3)+W(7)                                                    POTE-430
      W(11)=W(11)+W(15)                                                 POTE-431
      W(19)=W(19)+W(23)                                                 POTE-432
      W(13)=(W(13)+W(16))/WVX(6)                                        POTE-433
      W(21)=(W(21)+W(24))/WVX(6)                                        POTE-434
      W(14)=W(14)/WVX(6)                                                POTE-435
      W(22)=W(22)/WVX(6)                                                POTE-436
      W(15)=WVX(6)+WVX(7)+W(1)-W(3)                                     POTE-437
      W(16)=W(2)-W(4)                                                   POTE-438
      V(IS,N+7)=(WVX(7)-WVX(6)-W(1)-W(3))/CHB                           POTE-439
      V(IS,N+8)=-(W(2)+W(4))/CHB                                        POTE-440
      DD=W(15)**2+W(16)**2                                              POTE-441
      W(5)=W(15)/DD                                                     POTE-442
      W(6)=-W(16)/DD                                                    POTE-443
      DD=DSQRT(DSQRT(DD)/(WVX(6)+WVX(7)))                               POTE-444
      AR=.5D0*DATAN2(W(16),W(15))                                       POTE-445
      V(IS,N+9)=DD*DCOS(AR)                                             POTE-446
      V(IS,N+10)=DD*DSIN(AR)                                            POTE-447
      AR=W(11)-W(9)                                                     POTE-448
      AI=W(12)-W(10)                                                    POTE-449
      BR=W(17)-W(19)                                                    POTE-450
      BI=W(18)-W(20)                                                    POTE-451
      CR=AR*W(5)-AI*W(6)+W(13)                                          POTE-452
      CI=AR*W(6)+AI*W(5)+W(14)                                          POTE-453
      ER=BR*W(5)-BI*W(6)+W(13)*(CR+CR-W(13))-W(14)*(CI+CI-W(14))-W(21)+2POTE-454
     1.D0*CR/RR                                                         POTE-455
      EI=BR*W(6)+BI*W(5)+W(13)*(CI+CI-W(14))+W(14)*(CR+CR-W(13))-W(22)+2POTE-456
     1.D0*CI/RR                                                         POTE-457
      V(IS,N+1)=-W(1)-W(3)*WVX(7)/WVX(6)-(W(1)*W(1)-W(2)*W(2)-W(3)*W(3)+POTE-458
     1W(4)*W(4))/(2.D0*WVX(6))-(.75D0*(CR**2-CI**2)-.5D0*ER)*XP*CHB     POTE-459
      V(IS,N+2)=-W(2)-W(4)*WVX(7)/WVX(6)-(W(1)*W(2)-W(3)*W(4))/WVX(6)-(1POTE-460
     1.5D0*CR*CI-.5D0*EI)*XP*CHB                                        POTE-461
      V(IS,N+3)=.5D0*CR*XP*CHB/RR                                       POTE-462
      V(IS,N+4)=.5D0*CI*XP*CHB/RR                                       POTE-463
      V(IS,N+5)=W(15)/CHB                                               POTE-464
      V(IS,N+6)=W(16)/CHB                                               POTE-465
      V(IS,N+11)=W(5)*CHB                                               POTE-466
      V(IS,N+12)=W(6)*CHB                                               POTE-467
      V(IS,N+13)=W(13)/CHB                                              POTE-468
   42 V(IS,N+14)=W(14)/CHB                                              POTE-469
      VCO(1,J)=W(7)*RR/CHB                                              POTE-470
      VCO(2,J)=RR*(W(7)/(WVX(6)+WVX(7))+W(8)/WVX(6))                    POTE-471
   43 IF (JJ.LE.NCOLL) GO TO 44                                         POTE-472
      NX1=2*ISM+5                                                       POTE-473
      NX2=NX1+10*(IPIX(10)+2)                                           POTE-474
      HH=WV(8,1)                                                        POTE-475
      IF (NX2.GT.IDT) CALL MEMO('POTE',IDT,NX2)                         POTE-476
      IF ((.NOT.LO(93)).OR.(DABS(WV(2,1)-WVX(2)).GT.0.5D0)) HH=HH*WV(2,1POTE-477
     1)/WVX(2)                                                          POTE-478
      CALL TLNC(HH,IPIX,WVX,TL(IPIX(8)+1),ISM,Q,Q(NX1),IPIX(10)+2,NP(1,NPOTE-479
     1X2),Q(NX2),IDT-NX2,V,VCO,LO)                                      POTE-480
   44 CONTINUE                                                          POTE-481
      IF (INTV.EQ.0) GO TO 68                                           POTE-482
C STORAGE OF TRANSITION FORM FACTORS.                                   POTE-483
      DO 67 J=1,INTV                                                    POTE-484
      IJ=IVZ(1,J)                                                       POTE-485
      IF (.NOT.LO(100)) GO TO 45                                        POTE-486
      IK=ITX(6)+4*J                                                     POTE-487
      KI=ITX(5)+11*IJ                                                   POTE-488
   45 L=IVY(4,IJ)                                                       POTE-489
      IF (L.EQ.0) GO TO 50                                              POTE-490
      K=IVY(2,IJ)                                                       POTE-491
      IF (IVQ(2,K).GE.0) GO TO 50                                       POTE-492
C MAGNETIC COULOMB TRANSITION.                                          POTE-493
      DD=2.D0*(WV(8,1)/WV(10,1))**2/CHB                                 POTE-494
      IF (LO(100)) GO TO 47                                             POTE-495
      DO 46 IS=1,ISM                                                    POTE-496
      IF (LO(12)) V(IS,J+ITY(6))=0.D0                                   POTE-497
   46 V(IS,J+ITY(11))=-V(IS,L+ITX(15))*DD                               POTE-498
      IF (IVY(7,IJ).LE.MCM(1)) VDO(1,J)=V(ISM,L+ITX(15))*DD*RM**(IVY(7,IPOTE-499
     1J)+1)                                                             POTE-500
      GO TO 67                                                          POTE-501
   47 DO 49 IS=1,ISM                                                    POTE-502
      DO 48 M=1,4                                                       POTE-503
   48 V(IS,IK+M)=0.D0                                                   POTE-504
   49 V(IS,IK+3)=V(IS,KI+7)*DD                                          POTE-505
      IF (IVY(7,IJ).LE.MCM(1)) VDO(1,J)=V(ISM,KI+7)*DD*RM**(IVY(7,IJ)+1)POTE-506
      GO TO 67                                                          POTE-507
C CENTRAL TRANSITION FORM-FACTOR.                                       POTE-508
   50 I1=IVZ(2,J)/(NCOLL+1)                                             POTE-509
      I2=MOD(IVZ(2,J),NCOLL+1)                                          POTE-510
      DO 51 K=1,9                                                       POTE-511
      WVX(K)=WV(K+13,1)                                                 POTE-512
      IF (I2.NE.0) WVX(K)=.5D0*(WV(K+13,I1)+WV(K+13,I2))                POTE-513
      IF (.NOT.LO(10)) WVX(K)=0.D0                                      POTE-514
   51 CONTINUE                                                          POTE-515
      IF (LO(100)) GO TO 54                                             POTE-516
      DO 52 IS=1,ISM                                                    POTE-517
      IF (LO(12)) V(IS,J+ITY(6))=V(IS,IJ+ITX(10))+V(IS,IJ+ITX(12))+WVX(1POTE-518
     1)*V(IS,IJ+ITX(10))+WVX(4)*V(IS,IJ+ITX(12))                        POTE-519
   52 V(IS,J+ITY(5))=V(IS,IJ+ITX(9))*(1.D0+WVX(3))+V(IS,IJ+ITX(11))*(1.DPOTE-520
     10+WVX(6))+WVX(2)*V(IS,IJ+ITX(10))+WVX(5)*V(IS,IJ+ITX(12))         POTE-521
      IF (L.EQ.0) GO TO 59                                              POTE-522
C COULOMB TRANSITION FORM FACTOR FOR SCHROEDINGER EQUATION.             POTE-523
      DO 53 IS=1,ISM                                                    POTE-524
   53 V(IS,J+ITY(11))=V(IS,J+ITY(11))-V(IS,L+ITX(15))                   POTE-525
      IF (IVY(7,IJ).LE.MCM(1)) VDO(1,J)=V(ISM,L+ITX(15))*RM**(IVY(7,IJ)+POTE-526
     11)                                                                POTE-527
      GO TO 59                                                          POTE-528
   54 DO 55 IS=1,ISM                                                    POTE-529
      V(IS,IK+1)=V(IS,KI+1)+WVX(2)*V(IS,KI+2)                           POTE-530
      V(IS,IK+2)=V(IS,KI+2)+WVX(1)*V(IS,KI+2)                           POTE-531
      V(IS,IK+3)=V(IS,KI+3)+WVX(4)*V(IS,KI+4)                           POTE-532
   55 V(IS,IK+4)=V(IS,KI+4)+WVX(3)*V(IS,KI+4)                           POTE-533
      IF (L.EQ.0) GO TO 57                                              POTE-534
C COULOMB TRANSITION FORM FACTOR FOR DIRAC EQUATION.                    POTE-535
      DO 56 IS=1,ISM                                                    POTE-536
   56 V(IS,IK+3)=V(IS,IK+3)+V(IS,KI+7)                                  POTE-537
      IF (IVY(7,IJ).LE.MCM(1)) VDO(1,J)=V(ISM,KI+7)*RR**(IVY(7,IJ)+1)   POTE-538
   57 DO 58 IS=1,ISM                                                    POTE-539
      CR=V(IS,IK+1)                                                     POTE-540
      CI=V(IS,IK+2)                                                     POTE-541
      V(IS,IK+1)=XP*(V(IS,IK+3)+CR)                                     POTE-542
      V(IS,IK+2)=XP*(V(IS,IK+4)+CI)                                     POTE-543
      V(IS,IK+3)=XP*(V(IS,IK+3)-CR)                                     POTE-544
   58 V(IS,IK+4)=XP*(V(IS,IK+4)-CI)                                     POTE-545
   59 IF (J.GT.INVT) GO TO 67                                           POTE-546
C SPIN-ORBIT TRANSITION FORM FACTOR FOR SCHROEDINGER EQUATION.          POTE-547
      L=IVY(3,IJ)                                                       POTE-548
      IF (L.EQ.0) GO TO 67                                              POTE-549
      IF (LO(100)) GO TO 64                                             POTE-550
      DO 60 IS=1,ISM                                                    POTE-551
      V(IS,L+ITY(7))=V(IS,L+ITX(13))*(1.D0+WVX(9))+V(IS,L+ITX(14))*WVX(8POTE-552
     1)                                                                 POTE-553
   60 V(IS,L+ITY(7)+INSL)=V(IS,L+INLS+ITX(13))*(1.D0+WVX(9))+V(IS,L+ITX(POTE-554
     114))*WVX(8)                                                       POTE-555
      IF (.NOT.LO(14)) GO TO 62                                         POTE-556
      DO 61 IS=1,ISM                                                    POTE-557
      V(IS,L+ITY(8))=V(IS,L+ITX(14))*(1.D0+WVX(7))                      POTE-558
   61 V(IS,L+ITY(8)+INSL)=V(IS,L+INLS+ITX(14))*(1.D0+WVX(7))            POTE-559
   62 K=IVY(5,IJ)                                                       POTE-560
      IF (K.EQ.0) GO TO 67                                              POTE-561
      DO 63 IS=1,ISM                                                    POTE-562
      V(IS,L+ITY(12))=V(IS,L+ITY(12))+V(IS,K+ITX(16))*XQ                POTE-563
   63 V(IS,L+ITY(12)+INSL)=V(IS,L+ITY(12)+INSL)+V(IS,K+INVD+ITX(16))*XQ POTE-564
      IF (IVY(7,IJ).LE.MCM(2)) VDO(2,L)=-V(ISM,K+INVD+ITX(16))*RM**(IVY(POTE-565
     17,IJ)+3)*XQ                                                       POTE-566
      GO TO 67                                                          POTE-567
C SPIN-ORBIT TRANSITION FORM FACTOR FOR DIRAC EQUATION.                 POTE-568
   64 IL=ITX(6)+4*(L+INTV)                                              POTE-569
      RR=0.D0                                                           POTE-570
      DO 65 IS=1,ISM                                                    POTE-571
      RR=RR+WV(8,1)                                                     POTE-572
      V(IS,IL+1)=-XQ*V(IS,KI+9)*(1.D0+WVX(6))                           POTE-573
      V(IS,IL+2)=-XQ*V(IS,KI+10)*(1.D0+WVX(5))                          POTE-574
      V(IS,IL+3)=XQ*V(IS,KI+5)/RR*(1.D0+WVX(6))                         POTE-575
   65 V(IS,IL+4)=XQ*V(IS,KI+6)/RR*(1.D0+WVX(5))                         POTE-576
      IF (IVY(5,IJ).EQ.0) GO TO 67                                      POTE-577
      RR=0.D0                                                           POTE-578
      DO 66 IS=1,ISM                                                    POTE-579
      RR=RR+WV(8,1)                                                     POTE-580
      V(IS,IL+1)=V(IS,IL+1)-XQ*V(IS,KI+11)                              POTE-581
   66 V(IS,IL+3)=V(IS,IL+3)+XQ*V(IS,KI+8)/RR                            POTE-582
      IF (IVY(7,IJ).LE.MCM(2)) VDO(2,L)=XP*V(ISM,KI+8)*RR**(IVY(7,IJ)+1 POTE-583
     1)                                                                 POTE-584
   67 CONTINUE                                                          POTE-585
C OUTPUT OF POTENTIALS.                                                 POTE-586
   68 IF (.NOT.LO(51)) GO TO 77                                         POTE-587
      IF (LO(100)) GO TO 74                                             POTE-588
      DO 69 J=1,NPX                                                     POTE-589
      WRITE (MW,1006) J                                                 POTE-590
      WRITE (MW,1007) (I,V(I,J+ITY(1)),V(I,J+ITY(2)),I=1,ISM)           POTE-591
      IF (LO(133)) WRITE (MW,1008) (I,V(I,J+ITY(9)),I=1,ISM)            POTE-592
      IF (LO(129)) WRITE (MW,1009) (I,V(I,J+ITY(3)),I=1,ISM)            POTE-593
      IF (LO(130)) WRITE (MW,1010) (I,V(I,J+ITY(4)),I=1,ISM)            POTE-594
      IF (LO(133).AND.LO(103)) WRITE (MW,1011) (I,V(I,J+ITY(10)),I=1,ISMPOTE-595
     1)                                                                 POTE-596
   69 CONTINUE                                                          POTE-597
      IF (ITY(2).EQ.ITY(5)) GO TO 72                                    POTE-598
      NVX=ITY(2)-ITY(5)                                                 POTE-599
      WRITE (MW,1012)                                                   POTE-600
      NVB=MIN0(NVX,6)                                                   POTE-601
      DO 70 I=1,ISM                                                     POTE-602
      WRITE (MW,1013) I,(J,V(I,J+ITY(5)),J=1,NVB)                       POTE-603
      IF (NVB.NE.NVX) WRITE (MW,1014) (J,V(I,J+ITY(5)),J=7,NVX)         POTE-604
   70 CONTINUE                                                          POTE-605
      ITMX=ITX(1)                                                       POTE-606
      IF (LO(133)) ITMX=ITY(9)                                          POTE-607
      IF (ITMX.EQ.ITY(6)) GO TO 72                                      POTE-608
      MVX=ITMX-ITY(6)                                                   POTE-609
      WRITE (MW,1015)                                                   POTE-610
      NVB=MIN0(MVX,6)                                                   POTE-611
      DO 71 I=1,ISM                                                     POTE-612
      WRITE (MW,1013) I,(J,V(I,J+ITY(6)),J=1,NVB)                       POTE-613
      IF (NVB.NE.MVX) WRITE (MW,1014) (J,V(I,J+ITY(6)),J=7,MVX)         POTE-614
   71 CONTINUE                                                          POTE-615
   72 IF (.NOT.LO(133)) GO TO 77                                        POTE-616
      IF (ITX(1).EQ.ITY(11)) GO TO 77                                   POTE-617
      MVX=ITX(1)-ITY(11)                                                POTE-618
      WRITE (MW,1016)                                                   POTE-619
      NVB=MIN0(MVX,6)                                                   POTE-620
      DO 73 I=1,ISM                                                     POTE-621
      WRITE (MW,1013) I,(J,V(I,J+ITY(11)),J=1,NVB)                      POTE-622
      IF (NVB.NE.MVX) WRITE (MW,1014) (J,V(I,J+ITY(11)),J=7,MVX)        POTE-623
   73 CONTINUE                                                          POTE-624
      GO TO 77                                                          POTE-625
   74 N=0                                                               POTE-626
      DO 75 L=1,NCOLL                                                   POTE-627
      WRITE (MW,1017) L                                                 POTE-628
      WRITE (MW,1018) (I,(V(I,N+J),J=1,6),I=1,ISM)                      POTE-629
      WRITE (MW,1019) L                                                 POTE-630
      WRITE (MW,1018) (I,(V(I,N+J),J=7,12),I=1,ISM)                     POTE-631
      IF (LO(101).OR.LO(103)) WRITE (MW,1020) L,(I,(V(I,N+J),J=13,14),I=POTE-632
     11,ISM)                                                            POTE-633
   75 N=N+14                                                            POTE-634
      L=INTV+INSL                                                       POTE-635
      IF (L.EQ.0) GO TO 77                                              POTE-636
      N=ITX(2)                                                          POTE-637
      DO 76 K=1,L                                                       POTE-638
      IF (K.LE.INTV) WRITE (MW,1021) K                                  POTE-639
      IF (K.GT.INTV) WRITE (MW,1022) K                                  POTE-640
      WRITE (MW,1023) (I,(V(I,N+J),J=1,4),I=1,ISM)                      POTE-641
   76 N=N+4                                                             POTE-642
   77 IDT=IDY                                                           POTE-643
      RETURN                                                            POTE-644
 1000 FORMAT ('<POTENTI.>',F10.2,1P,D20.8,0P,F10.2,2I5)                 POTE-645
 1001 FORMAT (2X,3A8,' POTENTIAL ***',I2,' ***',5X,I5)                  POTE-646
 1002 FORMAT (2X,F10.5,1P,D20.6,0P,F10.5,1P,D20.6,A4,I6)                POTE-647
 1003 FORMAT (F10.5,1P,D20.6,0P,F10.5,1P,D20.6,A4,I6)                   POTE-648
 1004 FORMAT (2X,3A8,' TRANSITION POTENTIAL ***',I2,' ***',5X,I5)       POTE-649
 1005 FORMAT (' FORM FACTOR',I4,'  POTENTIAL',I3,'  RATIO OF CORRECTING POTE-650
     1DEFORMATION',1P,D15.6)                                            POTE-651
 1006 FORMAT (//' POTENTIAL ****',I2,' ****')                           POTE-652
 1007 FORMAT (//' CENTRAL POTENTIAL:'//(3(5X,I5,1P,2D14.5,' I')))       POTE-653
 1008 FORMAT (//' COULOMB POTENTIAL:'//(6(2X,I4,1P,D14.5)))             POTE-654
 1009 FORMAT (//' REAL SPIN-ORBIT POTENTIAL:'//(6(2X,I4,1P,D14.5)))     POTE-655
 1010 FORMAT (//' IMAGINARY SPIN-ORBIT POTENTIAL:'//(6(2X,I4,1P,D14.5)))POTE-656
 1011 FORMAT (//' COULOMB SPIN-ORBIT POTENTIAL:'//(6(2X,I4,1P,D14.5)))  POTE-657
 1012 FORMAT (//' REAL MULTIPOLES:'/)                                   POTE-658
 1013 FORMAT (4X,I4,6(I4,1P,D14.5))                                     POTE-659
 1014 FORMAT ((8X,6(I4,1P,D14.5)))                                      POTE-660
 1015 FORMAT (//' IMAGINARY MULTIPOLES:'/)                              POTE-661
 1016 FORMAT (//' COULOMB MULTIPOLES:'/)                                POTE-662
 1017 FORMAT (//' POTENTIALS OF SCHROEDINGER''S EQUATION FOR CHANNEL',I3POTE-663
     1/15X,'CENTRAL',27X,'SPIN-ORBIT',28X,'D(R)')                       POTE-664
 1018 FORMAT (2X,I3,1P,2D14.7,' I',5X,2D14.7,' I',5X,2D14.7,' I')       POTE-665
 1019 FORMAT (/17X,'E(R)',28X,'SQRT(D(R))',25X,'D(R)**(-1)',3X,'FOR CHANPOTE-666
     1NEL',I3)                                                          POTE-667
 1020 FORMAT (/45X,'TENSOR POTENTIAL',36X,'FOR CHANNEL',I3/(3(2X,I3,1P,2POTE-668
     1D14.7,' I')))                                                     POTE-669
 1021 FORMAT (//' SCALAR AND VECTOR MULTIPOLES',I6/)                    POTE-670
 1022 FORMAT (//' TENSOR MULTIPOLES',I6/)                               POTE-671
 1023 FORMAT (5X,I5,1P,2D26.7,' I',10X,2D16.7,' I')                     POTE-672
      END                                                               POTE-673
