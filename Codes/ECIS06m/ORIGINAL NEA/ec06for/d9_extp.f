C 28/02/07                                                      ECIS06  EXTP-000
      SUBROUTINE EXTP(NCOLL,NCOLT,WV,NIV,IQ,IVY,IVQ,IPI,FAC,VAL,NVAL,LL,EXTP-001
     1IPP,PIP,IDT,LO)                                                   EXTP-002
C INPUT AND SETUP OF EXTERNAL FORM FACTORS, ELASTIC AND TRANSITIONS.    EXTP-003
C INPUT:     NCOLL:   NUMBER OF COUPLED NUCLEAR STATES.                 EXTP-004
C            NCOLT:   NUMBER OF NUCLEAR STATES INCLUDING UNCOUPLED ONES.EXTP-005
C            WV:      MASSES OF PARTICLE AND TARGET IN WV(1/2,*).       EXTP-006
C            NIV:     ADDRESSES IN TABLE OF NUCLEAR MATRIX ELEMENTS.    EXTP-007
C            IQ:      TABLE OF REDUCED NUCLEAR MATRIX ELEMENTS.         EXTP-008
C            IVY:     TABLE OF FORM FACTORS (SEE REDM).                 EXTP-009
C            IVQ:     TABLE OF ANGULAR MOMENTA (SEE REDM).              EXTP-010
C            IPI(5,*):REFERENCE TO POTENTIALS.                          EXTP-011
C            FAC:     TABLE OF LOGARITHM OF FACTORIALS.                 EXTP-012
C            IPP,PIP: DISPERSION RELATIONS, EQUIVALENT BY CALL.         EXTP-013
C            IDT:     SIZE OF AVAILABLE WORKING SPACE.                  EXTP-014
C            LO(I):   LOGICAL CONTROLS:                                 EXTP-015
C               LO(6)  =.TRUE. USE DEFORMATION LENGTHS.                 EXTP-016
C               LO(10) =.TRUE. DISPERSION RELATIONS FOR POTENTIALS.     EXTP-017
C               LO(12) =.TRUE. DEFORMED IMAGINARY POTENTIAL.            EXTP-018
C               LO(14) =.TRUE. DEFORMED IMAGINARY SPIN-ORBIT OR TENSOR. EXTP-019
C               LO(16) =.TRUE. HEAVY-ION DEFINITION OF REDUCED RADII ANDEXTP-020
C                              DEFORMATIONS.                            EXTP-021
C               LO(73) =.TRUE. NO OUTPUT OF EXTERNAL POTENTIALS WHEN    EXTP-022
C                              THEY ARE READ.                           EXTP-023
C               LO(99) =.TRUE. SCHROEDINGER EQUIVALENT TO DIRAC         EXTP-024
C                              EQUATION.                                EXTP-025
C               LO(100)=.TRUE. DIRAC EQUATION.                          EXTP-026
C               LO(101)=.TRUE. THERE IS A REAL SPIN-ORBIT POTENTIAL.    EXTP-027
C               LO(102)=.TRUE. THERE IS AN IMAGINARY SPIN-ORBIT         EXTP-028
C                              POTENTIAL.                               EXTP-029
C OUTPUT:    VAL,NVAL:IN EQUIVALENCE BY CALL, FOR PARAMETERS READ HERE. EXTP-030
C                     NVA(1) FIRST ADDRESS AFTER LL.                    EXTP-031
C                     NVA(2) NUMBER OF SETS OF FOLDING PARAMETERS.      EXTP-032
C                     NVA(3) FIRST ADDRESS OF FOLDING PARAMETERS.       EXTP-033
C                     NVA(4) LAST ADDRESS OF FOLDING PARAMETERS.        EXTP-034
C            LL:      ADDRESSES IN VA IN EQUIVALENCE WITH NVA(1,3):     EXTP-035
C                     LL(1,ITYP,K) FIRST ADDRESS OF PARAMETERS WHICH CANEXTP-036
C                     BE VARIED FOR ITYP AND FORM FACTOR K, 1 FOR FORM  EXTP-037
C                     FACTORS NOT USED.                                 EXTP-038
C                     LL(2,ITYP,K) LAST ADDRESS FOR ITYP AND FORM FACTOREXTP-039
C                     K, -1 FOR FORM FACTORS NOT USED.                  EXTP-040
C                                                                       EXTP-041
C FOR THE COMMON  /POTE1/ SEE REDM.                                     EXTP-042
C                                                                       EXTP-043
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE1/:                     EXTP-044
C  INTC:      NUMBER OF FORM FACTORS WITHOUT DEFORMED SPIN-ORBIT        EXTP-045
C             INCLUDING CORRECTION TERMS.                               EXTP-046
C  INLS:      NUMBER OF SPIN-ORBIT FORM FACTORS NOT TAKING INTO ACCOUNT EXTP-047
C             MULTIPLICATION BY 2.                                      EXTP-048
C  INVC:      NUMBER OF COULOMB TRANSITION FORM FACTORS.                EXTP-049
C  INVD:      IDEM FOR COULOMB SPIN-ORBIT.                              EXTP-050
C   USED:     INTC,INLS,INVC,INVD.                                      EXTP-051
C                                                                       EXTP-052
C***********************************************************************EXTP-053
C ITYP 1 REAL VOLUME OR DIRAC SCALAR POTENTIAL.                         EXTP-054
C      2 IMAGINARY VOLUME OR DIRAC SCALAR POTENTIAL.                    EXTP-055
C      3 REAL SURFACE OR DIRAC VECTOR POTENTIAL.                        EXTP-056
C      4 IMAGINARY SURFACE OR DIRAC VECTOR POTENTIAL.                   EXTP-057
C      5 REAL SPIN-ORBIT OR DIRAC TENSOR POTENTIAL.                     EXTP-058
C      6 IMAGINARY SPIN-ORBIT OR DIRAC TENSOR POTENTIAL.                EXTP-059
C      7 COULOMB POTENTIAL.                                             EXTP-060
C      8 COULOMB SPIN-ORBIT POTENTIAL.                                  EXTP-061
C                                                                       EXTP-062
C L1,L2 CONTROLS THE LEVEL TO POTENTIAL ASSIGNMENT.                     EXTP-063
C ML    IS 0 FOR THE POTENTIAL OR THE PLACE OF THE TRANSITION FORM.     EXTP-064
C       FACTOR IN THE SEQUENCE OF REDUCED MATRIX ELEMENTS.              EXTP-065
C       THE FORM FACTORS CONTAIN THE DEFORMATION.                       EXTP-066
C L1X,L2X,MLX  FORM FACTOR COPIED TO L1,L2,ML.                          EXTP-067
C                                                                       EXTP-068
C THE SPIN-ORBIT TRANSITION FORM FACTOR TO BE READ IS THE SECOND ONE    EXTP-069
C MULTIPLIED BY R**2. (MULTIPOLE OF AN ORDINARY WOODS-SAXON POTENTIAL). EXTP-070
C                                                                       EXTP-071
C ITYPX  -1  WOODS-SAXON POTENTIAL.                                     EXTP-072
C        -2  FIRST DERIVATIVE   MULTIPLIED BY R/SQRT(4*PI).             EXTP-073
C        -3  SECOND DERIVATIVE  MULTIPLIED BY R**2/(8*PI).              EXTP-074
C        -4  THIRD DERIVATIVE   MULTIPLIED BY R**3/(48*PI**(3/2)).      EXTP-075
C        -5  DEFORMED WOODS-SAXON POTENTIAL.                            EXTP-076
C        -6  DERIVATIVE OF DEFORMED WOODS-SAXON POTENTIAL.              EXTP-077
C        -7  LAGUERRE POLYNOMIAL.                                       EXTP-078
C        -8  SOLUTION IN REAL WOODS-SAXON POTENTIAL.                    EXTP-079
C        -9  BESSEL EXPANSION.                                          EXTP-080
C       -10  LAGUERRE EXPANSION.                                        EXTP-081
C L2X GIVES THE NUMBER OF DEFORMATIONS, OF NODES, OF BESSEL FUNCTIONS,  EXTP-082
C                      OF LAGUERRE POLYNOMIALS.                         EXTP-083
C MLX IS THE L-VALUE OF BESSEL OR LAGUERRE EXPANSION, OR NUMBER OF BOUNDEXTP-084
C FUNCTIONS: =0 OR 1 FOR ONE WITH THE QUANTUM NUMBERS OF THE TRANSITION,EXTP-085
C            =2 FOR TWO FUNCTIONS WITH THE SAME ITYPX,                  EXTP-086
C            =3 WHEN ITYPX=-8 FOR LAGUERRE POLYNOMIAL AS THE SECOND ONE.EXTP-087
C MLX AND -L1X ARE QUANTUM NUMBERS OF VIBRATIONAL BAND:                 EXTP-088
C -L1X GIVES THE MULTIPLICATION OF STEP IN COMPUTING BOUND FUNCTION,    EXTP-089
C -L1X IS THE ORDER OF DERIVATION OF BESSEL OR LAGUERRE EXPANSION.      EXTP-090
C                                                                       EXTP-091
C ALLOWED VALUES OF ITYPX FOR STANDARD POTENTIALS                       EXTP-092
C ***********************************************                       EXTP-093
C **** ML IS 0 ****                                                     EXTP-094
C ITYPX =    -1       -2,-4    -5       -6       -7,-8    -9,-10        EXTP-095
C ITYP = 1   YES      NO       YES      NO       NO       YES           EXTP-096
C ITYP = 2   YES      NO       YES      NO       NO       YES           EXTP-097
C ITYP = 3   YES      NO       YES      NO       NO       YES           EXTP-098
C ITYP = 4   YES      NO       YES      NO       NO       YES           EXTP-099
C ITYP = 5   YES      NO       YES      NO       NO       YES           EXTP-100
C ITYP = 6   YES      NO       YES      NO       NO       YES           EXTP-101
C ITYP = 7   YES      NO       YES      NO       NO       YES           EXTP-102
C ITYP = 8   YES      NO       YES      NO       NO       YES           EXTP-103
C **** ML IS NOT 0 ****                                                 EXTP-104
C ITYP = 1   YES      YES      YES      YES      YES      YES           EXTP-105
C ITYP = 2   YES      YES      YES      YES      YES      YES           EXTP-106
C ITYP = 3   YES      YES      YES      YES      NO       YES           EXTP-107
C ITYP = 4   YES      YES      YES      YES      NO       YES           EXTP-108
C ITYP = 5   YES      YES      YES      YES      NO       YES           EXTP-109
C ITYP = 6   YES      YES      YES      YES      NO       YES           EXTP-110
C ITYP = 7   NO       YES      YES      YES      NO       YES           EXTP-111
C ITYP = 8   NO       YES      YES      YES      NO       YES           EXTP-112
C NUMBER OF PARAMETERS TO STORE:                                        EXTP-113
C ITYPX = -1,-4       -5        -6        -7        -8        -9,-10    EXTP-114
C INTEGER 6           7         8         15        15        9         EXTP-115
C REAL*8  4           4+L2X     4+L2X     1         23        2+L2X     EXTP-116
C THE ODD NUMBER OF INTEGER NEEDED IS ROUNDED TO NEXT EVEN VALUE        EXTP-117
C THERE IS ONE MORE FOR COULOMB POTENTIALS AND ITYPX=-1 TO -6.          EXTP-118
C THERE ARE 5 MORE PARAMETERS FOR ITYPX=-7 AND MLX=2.                   EXTP-119
C THERE ARE 16 OR 6 MORE PARAMETERS FOR ITYPX=-8 AND MLX=2 OR MLX=3.    EXTP-120
C                                                                       EXTP-121
C FOR THE COMMON  /COUPL/ SEE CALC.                                     EXTP-122
C                                                                       EXTP-123
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     EXTP-124
C  NFA:       NUMBER OF LOGARITHMS OF FACTORIALS.                       EXTP-125
C  NPP:       NUMBER OF OPTICAL POTENTIALS.                             EXTP-126
C   USED:     NFA,NPP.                                                  EXTP-127
C                                                                       EXTP-128
C***********************************************************************EXTP-129
      IMPLICIT REAL*8 (A-H,O-Z)                                         EXTP-130
      LOGICAL LO(150)                                                   EXTP-131
      DIMENSION WV(22,*),NIV(NCOLL,NCOLL,*),IQ(6,*),IVY(7,*),IVQ(3,*),IPEXTP-132
     1I(11,*),FAC(*),VAL(*),NVAL(*),LL(2,8,*),IPP(2,18,*),PIP(18,1),ITZ(EXTP-133
     210)                                                               EXTP-134
      CHARACTER*8 AA(3,8)                                               EXTP-135
      CHARACTER*4 IERM,LAST                                             EXTP-136
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   EXTP-137
      COMMON /INOUT/ MR,MW,MS                                           EXTP-138
      COMMON /POTE1/ ITX(16),IMAX,INTC,INLS,INVC,INVD,ITXM              EXTP-139
      DATA AA /'      RE','AL VOLUM','E/SCALAR',' IMAGINA','RY VOLUM','EEXTP-140
     1/SCALAR','     REA','L SURFAC','E/VECTOR','  IMAGIN','. SURFAC','EEXTP-141
     2/VECTOR','  REAL S','PIN-ORBI','T/TENSOR',' IMAG. S','PIN-ORBI','TEXTP-142
     3/TENSOR','        ','        ',' COULOMB','      SP','IN-ORBIT',' EXTP-143
     4COULOMB'/                                                         EXTP-144
      DATA ITZ /2,2,2,2,3,4,7,7,4,4/                                    EXTP-145
      DATA IERM /'LAST'/                                                EXTP-146
      IF (LO(99)) GO TO 44                                              EXTP-147
      NPX=NPP+INTC                                                      EXTP-148
      NMA=3+8*NPX                                                       EXTP-149
      NVAL(1)=NMA                                                       EXTP-150
      IF (NMA.GT.IDT) CALL MEMO('EXTP',IDT,NMA)                         EXTP-151
C SETTING THE CONTROLS OF INPUT                                         EXTP-152
      NTOT=0                                                            EXTP-153
      DO 9 I=1,NPX                                                      EXTP-154
      DO 1 J=1,8                                                        EXTP-155
      LL(1,J,I)=-1                                                      EXTP-156
    1 LL(2,J,I)=-1                                                      EXTP-157
      IF (I.GT.NPP) GO TO 2                                             EXTP-158
      IF (.NOT.LO(101)) LL(1,5,I)=1                                     EXTP-159
      IF (.NOT.LO(102)) LL(1,6,I)=1                                     EXTP-160
      IF (.NOT.LO(101)) LL(1,8,I)=1                                     EXTP-161
      GO TO 7                                                           EXTP-162
    2 K=I-NPP                                                           EXTP-163
      L=IVY(2,K)                                                        EXTP-164
      IF (IVQ(2,L).GE.0) GO TO 4                                        EXTP-165
      DO 3 L=1,6                                                        EXTP-166
    3 LL(1,L,I)=1                                                       EXTP-167
      LL(1,8,I)=1                                                       EXTP-168
      GO TO 7                                                           EXTP-169
    4 IF (LO(12)) GO TO 5                                               EXTP-170
      LL(1,2,I)=1                                                       EXTP-171
      LL(1,4,I)=1                                                       EXTP-172
    5 IF (LO(100)) GO TO 6                                              EXTP-173
      IF (K.GT.INLS) LL(1,5,I)=1                                        EXTP-174
      IF ((.NOT.LO(14)).OR.K.GT.INLS) LL(1,6,I)=1                       EXTP-175
      IF (K.GT.INVC) LL(1,7,I)=1                                        EXTP-176
      IF (K.GT.INVD) LL(1,8,I)=1                                        EXTP-177
      GO TO 7                                                           EXTP-178
    6 IF (IVY(3,K).EQ.0) LL(1,5,I)=1                                    EXTP-179
      IF ((.NOT.LO(14)).OR.IVY(3,K).EQ.0) LL(1,6,I)=1                   EXTP-180
      IF (IVY(4,K).EQ.0) LL(1,7,I)=1                                    EXTP-181
      IF (IVY(5,K).EQ.0) LL(1,8,I)=1                                    EXTP-182
C COUNT OF FORM FACTORS TO BE READ                                      EXTP-183
    7 DO 8 J=1,8                                                        EXTP-184
      IF (LL(1,J,I).EQ.-1) NTOT=NTOT+1                                  EXTP-185
    8 CONTINUE                                                          EXTP-186
    9 CONTINUE                                                          EXTP-187
      NFOLT=0                                                           EXTP-188
   10 IF (NTOT.LE.0) GO TO 42                                           EXTP-189
C INPUT OF A FORM FACTOR                                                EXTP-190
      READ (MR,1000,ERR= 45 ) L1,L2,ML,ITYP,L1X,L2X,MLX,ITYPX,NST,NFOLD,EXTP-191
     1MINT                                                              EXTP-192
      IF (ITYP.GT.8.OR.IABS(NST).GT.NCOLT.OR.ITYP.LE.0) GO TO 46        EXTP-193
      NFOLT=MAX0(NFOLT,NFOLD)                                           EXTP-194
      IF (IABS(ML).NE.0) GO TO 11                                       EXTP-195
C TEST FOR A POTENTIAL                                                  EXTP-196
      IF (L1.NE.L2.OR.L1.GT.NCOLT) GO TO 47                             EXTP-197
      LZ=IPI(5,L1)                                                      EXTP-198
      IF (NST.EQ.0) NST=L1                                              EXTP-199
      IF (.NOT.LO(73)) WRITE (MW,1001) (AA(K,ITYP),K=1,3),LZ            EXTP-200
      GO TO 14                                                          EXTP-201
C TEST FOR A TRANSITION                                                 EXTP-202
   11 IF (L1.GT.NCOLL.OR.L2.GT.NCOLL) GO TO 48                          EXTP-203
      LY=NIV(L1,L2,1)+IABS(ML)-1                                        EXTP-204
      IF (LY.GT.NIV(L1,L2,2)) GO TO 49                                  EXTP-205
      LX=IQ(1,LY)                                                       EXTP-206
      IF (ML.GT.0) GO TO 12                                             EXTP-207
      LX=IVY(6,LX)                                                      EXTP-208
      IF (LX.LE.0) GO TO 50                                             EXTP-209
   12 IF (LO(100)) GO TO 13                                             EXTP-210
      IF ((ITYP.EQ.5).OR.(ITYP.EQ.6)) LX=IVY(3,LX)                      EXTP-211
      IF (ITYP.GE.7) LX=IVY(ITYP-3,LX)                                  EXTP-212
   13 LZ=LX+NPP                                                         EXTP-213
      IF (NST.EQ.0) NST=1                                               EXTP-214
      IF (.NOT.LO(73)) WRITE (MW,1002) (AA(K,ITYP),K=1,3),L1,L2,ML      EXTP-215
C TEST IT IS NOT ALREADY READ AND THAT THERE IS PLACE                   EXTP-216
   14 IF (LL(1,ITYP,LZ).NE.-1) GO TO 51                                 EXTP-217
      IF (NMA+31.GT.IDT) CALL MEMO('EXTP',IDT,NMA+32)                   EXTP-218
      NTOT=NTOT-1                                                       EXTP-219
      NNA=2*NMA-1                                                       EXTP-220
      NVAL(NNA)=8*LZ+ITYP                                               EXTP-221
      NVAL(NNA+2)=NFOLD                                                 EXTP-222
      IF (ITYPX.LT.0.OR.L1X.LE.0) GO TO 21                              EXTP-223
      NVAL(NNA+5)=NMA                                                   EXTP-224
      NMD=NMA+3                                                         EXTP-225
      LL(1,ITYP,LZ)=NMD                                                 EXTP-226
      VAL(NMD)=1.D0                                                     EXTP-227
      IF (MINT.LT.0) READ (MR,1003,ERR= 59 ) VAL(NMD)                   EXTP-228
C USE OF ALREADY STORED FORM FACTOR/TEST EXISTENCE OF COPIED FORM FACTOREXTP-229
      IF (ITYP.NE.ITYPX) GO TO 52                                       EXTP-230
      IF (ML.NE.0) GO TO 15                                             EXTP-231
      IF (L1X.NE.L2X.OR.L1X.GT.NCOLT) GO TO 53                          EXTP-232
      LZX=IPI(5,L1X)                                                    EXTP-233
      IF (.NOT.LO(73)) WRITE (MW,1004) L1X,ITYPX,L1                     EXTP-234
      GO TO 17                                                          EXTP-235
   15 IF (L1X.GT.NCOLL.OR.L2X.GT.NCOLL) GO TO 54                        EXTP-236
      LYY=NIV(L1X,L2X,1)+MLX-1                                          EXTP-237
      IF (LYY.GT.NIV(L1X,L2X,2)) GO TO 55                               EXTP-238
      LXX=IQ(1,LYY)                                                     EXTP-239
      IF (LO(100)) GO TO 16                                             EXTP-240
      IF ((ITYP.EQ.5).OR.(ITYP.EQ.6)) LXX=IVY(3,LXX)                    EXTP-241
      IF (ITYP.GE.7) LXX=IVY(ITYP-3,LXX)                                EXTP-242
   16 LZX=LXX+NPP                                                       EXTP-243
      IF (.NOT.LO(73)) WRITE (MW,1005) L1X,L2X,MLX,ITYPX,L1,L2,ML       EXTP-244
   17 NMC=LL(1,ITYP,LZX)                                                EXTP-245
      IF (NMC.EQ.-1) GO TO 56                                           EXTP-246
C USE OF ALREADY STORED FORM FACTOR/TEST POSSIBILITY TO COPY            EXTP-247
      NMB=NVAL(2*NMC-2)                                                 EXTP-248
      ITYPX=-NVAL(2*NMB)                                                EXTP-249
      IF (NVAL(2*NMB+1).EQ.0.AND.NFOLD.NE.0) NVAL(2*NMB+1)=-1           EXTP-250
      IF (NVAL(2*NMB+1).NE.0.AND.NFOLD.EQ.0) NVAL(2*NMA+1)=-1           EXTP-251
      IF (ITYPX.NE.-16) GO TO 18                                        EXTP-252
      LZX=(NVAL(2*NMB+2)-1)/8                                           EXTP-253
      GO TO 17                                                          EXTP-254
   18 IF (ITYPX.GT.-5.OR.ML.EQ.0) GO TO 20                              EXTP-255
      IF (ITYPX.NE.-7.AND.ITYPX.NE.-8) GO TO 19                         EXTP-256
      LL(1,ITYP+2,LZ)=1                                                 EXTP-257
      NTOT=NTOT-1                                                       EXTP-258
      NMB=NMB+ITZ(-ITYPX)                                               EXTP-259
      IF (NVAL(2*NMB-1).EQ.1) GO TO 19                                  EXTP-260
      IF (NVAL(2*NMB+2).NE.1.OR.NVAL(2*NMB+6).NE.1) GO TO 20            EXTP-261
      K=IQ(2,LY)                                                        EXTP-262
      VAL(NMD)=DCGS(2*IVQ(1,K),NVAL(2*NMB+3),NVAL(2*NMB+7),FAC,NFA)*DSQREXTP-263
     1T(2.D0*IVQ(1,K)+1.D0)/VAL(NMC)                                    EXTP-264
      IF (VAL(NMD).EQ.0.D0) GO TO 57                                    EXTP-265
      GO TO 20                                                          EXTP-266
   19 IF (IQ(2,LY).NE.IQ(2,LYY)) GO TO 58                               EXTP-267
   20 NVAL(NNA+1)=16                                                    EXTP-268
      NVAL(NNA+3)=8*LZX+ITYP                                            EXTP-269
      NVAL(NNA+4)=NMA                                                   EXTP-270
      NMB=NMA+2                                                         EXTP-271
      IF (.NOT.LO(73)) WRITE (MW,1006) NMA,NMB,(NVAL(NNA+I),I=0,4),NMD,VEXTP-272
     1AL(NMD)                                                           EXTP-273
      LL(2,ITYP,LZ)=NMD                                                 EXTP-274
      NMA=NMD+1                                                         EXTP-275
      GO TO 41                                                          EXTP-276
   21 NVAL(NNA+3)=MINT                                                  EXTP-277
      NVAL(NNA+4)=NST                                                   EXTP-278
      IF (ITYPX.LT.0) GO TO 23                                          EXTP-279
      NVAL(NNA+5)=NMA                                                   EXTP-280
      LL(1,ITYP,LZ)=NMA+3                                               EXTP-281
C FORM FACTOR READ FROM CARDS                                           EXTP-282
      READ (MR,1003,ERR= 59 ) VAL(NMA+3),VAL(NMA+4)                     EXTP-283
      IF (VAL(NMA+3).EQ.0.D0) VAL(NMA+3)=1.D0                           EXTP-284
      IF (VAL(NMA+4).EQ.0.D0) VAL(NMA+4)=1.D0                           EXTP-285
      NM=NMA+5                                                          EXTP-286
   22 NMB=NM+3                                                          EXTP-287
      IF (NMB.GT.IDT) CALL MEMO('EXTP',IDT,NMB)                         EXTP-288
      READ (MR,1007,ERR= 59 ) (VAL(I),I=NM,NMB),LAST                    EXTP-289
      NM=NMB+1                                                          EXTP-290
      IF (LAST.NE.IERM) GO TO 22                                        EXTP-291
      NVN=NMA+2                                                         EXTP-292
      NM1=NVN+1                                                         EXTP-293
      NM2=NVN+2                                                         EXTP-294
      NM3=NVN+3                                                         EXTP-295
      NVAL(2*NMA)=1-(NMB-NVN)/2                                         EXTP-296
      IF (.NOT.LO(73)) WRITE (MW,1008) NMA,NVN,(NVAL(NNA+I-1),I=1,6),NM1EXTP-297
     1,NM2,VAL(NM1),VAL(NM2),NM3,NMB,(VAL(I),I=NM3,NMB)                 EXTP-298
      IF (NVAL(2*NMA).GT.-4) GO TO 60                                   EXTP-299
      LL(2,ITYP,LZ)=NMB                                                 EXTP-300
      NMA=NMB+1                                                         EXTP-301
      GO TO 41                                                          EXTP-302
C STANDARD FORM FACTORS                                                 EXTP-303
   23 IF (ITYPX.LT.-10) GO TO 61                                        EXTP-304
      ITYPY=-ITYPX                                                      EXTP-305
      IF (ITYPY.LT.9) GO TO 24                                          EXTP-306
      IF ((ML.NE.0).AND.(MLX.EQ.0)) MLX=IVY(7,LX)                       EXTP-307
      IF (MLX.LT.0) MLX=0                                               EXTP-308
      GO TO 25                                                          EXTP-309
   24 IF (((ITYPY-1)*(ITYPY-5).NE.0).AND.ML.EQ.0) GO TO 62              EXTP-310
      IF ((ITYPY.EQ.1).AND.(ML.NE.0).AND.(ITYP.GE.7)) GO TO 63          EXTP-311
      IF ((ITYPY.GE.7).AND.(ITYP.GT.2)) GO TO 64                        EXTP-312
   25 NVAL(NNA+1)=ITYPY                                                 EXTP-313
      NMB=NMA+ITZ(ITYPY)                                                EXTP-314
      IF (ITYPY.LT.5) GO TO 26                                          EXTP-315
      NVAL(NNA+5)=L2X                                                   EXTP-316
      NVAL(NNA+6)=MLX                                                   EXTP-317
      NVAL(NNA+7)=-L1X                                                  EXTP-318
      IF (ITYPY.EQ.7.OR.ITYPY.EQ.8) GO TO 32                            EXTP-319
      NVAL(2*NMB-1)=NMA                                                 EXTP-320
   26 NVAL(2*NMB)=NMA                                                   EXTP-321
      LL(1,ITYP,LZ)=NMB+1                                               EXTP-322
      IF (LO(73)) GO TO 27                                              EXTP-323
      IF (ITYPY.GE.9) WRITE (MW,1009) NMA,NMB,(NVAL(NNA+I-1),I=1,9)     EXTP-324
      IF (ITYPY.EQ.6) WRITE (MW,1010) NMA,NMB,(NVAL(NNA+I-1),I=1,9)     EXTP-325
      IF (ITYPY.EQ.5) WRITE (MW,1011) NMA,NMB,(NVAL(NNA+I-1),I=1,7)     EXTP-326
      IF (ITYPY.LT.5) WRITE (MW,1012) NMA,NMB,(NVAL(NNA+I-1),I=1,6)     EXTP-327
   27 IF (ITYPY.GE.9) GO TO 29                                          EXTP-328
      NMA=NMB+1                                                         EXTP-329
      NMB=NMB+4                                                         EXTP-330
      IF (ITYP.GT.6) NMB=NMB+1                                          EXTP-331
      READ (MR,1003,ERR= 59 ) (VAL(I),I=NMA,NMB)                        EXTP-332
      IF (ITYPY.GT.6) GO TO 28                                          EXTP-333
      K=IABS(NST)                                                       EXTP-334
      IF (K.EQ.NST) GO TO 28                                            EXTP-335
      IF (.NOT.LO(73)) WRITE (MW,1013) VAL(NMA),VAL(NMA+1)              EXTP-336
      EX=WV(2,K)**.33333333333333D0                                     EXTP-337
      EY=EX                                                             EXTP-338
      IF (LO(16)) EX=EX+WV(1,K)**.33333333333333D0                      EXTP-339
      EY=EY/EX                                                          EXTP-340
      VAL(NMA+1)=VAL(NMA+1)*EX                                          EXTP-341
      IF ((.NOT.LO(16)).OR.ML.EQ.0) GO TO 28                            EXTP-342
      ITYZ=ITYPY                                                        EXTP-343
      IF (ITYZ.GE.5) ITYZ=ITYZ-4                                        EXTP-344
      ITYW=1                                                            EXTP-345
      K=IQ(2,LY)                                                        EXTP-346
      IF (ITYP.GT.6) ITYW=ITYW*IVQ(1,K)                                 EXTP-347
      IF (LO(6)) ITYW=ITYW-1                                            EXTP-348
      IF (ITYZ.GT.1) VAL(NMA)=VAL(NMA)*EY**((ITYZ-1)*ITYW)              EXTP-349
   28 IF (.NOT.LO(73)) WRITE (MW,1014) NMA,NMB,(VAL(I),I=NMA,NMB)       EXTP-350
      IF (ITYPY.LT.5) GO TO 40                                          EXTP-351
   29 NMA=NMB+1                                                         EXTP-352
      NMC=NMA                                                           EXTP-353
      IF (ITYPY.GE.9) NMC=NMC+2                                         EXTP-354
      NMB=NMC+L2X-1                                                     EXTP-355
      IF (NMB.GT.IDT) CALL MEMO('EXTP',IDT,NMB)                         EXTP-356
      READ (MR,1003,ERR= 59 ) (VAL(I),I=NMA,NMB)                        EXTP-357
      IF ((ITYPY.GE.9).AND.(VAL(NMA).EQ.0.D0)) VAL(NMA)=1.D0            EXTP-358
      IF (NST.GT.0.OR.ITYPY.GE.9) GO TO 31                              EXTP-359
      IF (.NOT.LO(73)) WRITE (MW,1015) (VAL(I),I=NMA,NMB)               EXTP-360
      DO 30 I=NMA,NMB                                                   EXTP-361
      J=I-NMA                                                           EXTP-362
      IF (ITYP.LT.7) J=0                                                EXTP-363
      IF (.NOT.LO(6)) J=J+1                                             EXTP-364
   30 VAL(I)=VAL(I)*EY**J                                               EXTP-365
   31 IF (LO(73)) GO TO 40                                              EXTP-366
      IF (ITYPY.GE.9) WRITE (MW,1016) NMA,NMB,VAL(NMA),VAL(NMA+1),L2X,(IEXTP-367
     1,VAL(I),I=NMC,NMB)                                                EXTP-368
      IF (ITYPY.LT.9) WRITE (MW,1017) NMA,NMB,L2X,(I,VAL(I),I=NMA,NMB)  EXTP-369
      GO TO 40                                                          EXTP-370
   32 LL(1,ITYP+2,LZ)=1                                                 EXTP-371
      NTOT=NTOT-1                                                       EXTP-372
      IF (LO(100)) GO TO 65                                             EXTP-373
      NVAL(NNA+7)=-L1X                                                  EXTP-374
      IF (MLX.EQ.0) MLX=1                                               EXTP-375
      IF ((MLX.LT.0.OR.MLX.GT.3).OR.(MLX.EQ.3.AND.ITYPY.EQ.7)) GO TO 66 EXTP-376
      NVAL(NNA+6)=MLX                                                   EXTP-377
      MLY=MLX                                                           EXTP-378
      NMD=NMA+2                                                         EXTP-379
      IF (.NOT.LO(73)) WRITE (MW,1018) NMA,NMD,(NVAL(NNA+I-1),I=1,6)    EXTP-380
      NMD=NMD+1                                                         EXTP-381
      K=IQ(2,LY)                                                        EXTP-382
      IF (MLX.GT.1) GO TO 34                                            EXTP-383
      DO 33 I=1,3                                                       EXTP-384
   33 NVAL(NNA+7+I)=IVQ(I,K)                                            EXTP-385
      NVAL(NNA+11)=NMA                                                  EXTP-386
      NMB=NMD+2                                                         EXTP-387
      IF (.NOT.LO(73)) WRITE (MW,1019) NMD,NMB,(NVAL(NNA+I),I=6,11)     EXTP-388
      NMC=NMB                                                           EXTP-389
      GO TO 36                                                          EXTP-390
   34 READ (MR,1000) (NVAL(NNA+I),I=7,14),NVC                           EXTP-391
      IF (NVAL(NNA+9).EQ.1.AND.NVC.EQ.0) NVAL(NNA+9)=-3                 EXTP-392
      NMB=NMA+7                                                         EXTP-393
      NVAL(NNA+15)=NMA                                                  EXTP-394
      NMC=NMB+1                                                         EXTP-395
      VAL(NMC)=1.D0                                                     EXTP-396
      IF (NVAL(NNA+9).NE.1.OR.NVAL(NNA+13).NE.1.OR.NVC.EQ.0) GO TO 35   EXTP-397
      VAL(NMC)=DCGS(2*IVQ(1,K),NVAL(NNA+10),NVAL(NNA+14),FAC,NFA)*DSQRT(EXTP-398
     1DFLOAT(2*IVQ(1,K)+1))                                             EXTP-399
      IF (VAL(NMC).EQ.0.D0) GO TO 57                                    EXTP-400
   35 IF (.NOT.LO(73)) WRITE (MW,1020) NMD,NMB,(NVAL(NNA+I),I=6,15),NMC,EXTP-401
     1VAL(NMC)                                                          EXTP-402
   36 NMC=NMC+1                                                         EXTP-403
      LL(1,ITYP,LZ)=NMB+1                                               EXTP-404
      IF (ITYPY.EQ.8) GO TO 39                                          EXTP-405
   37 READ (MR,1003,ERR= 59 ) (VAL(NMC+I-1),I=1,3*MLX)                  EXTP-406
      IF (VAL(NMC+1).EQ.0.D0) VAL(NMC+1)=WV(2,L1)                       EXTP-407
      IF (VAL(NMC+2).EQ.0.D0) VAL(NMC+2)=1.D0                           EXTP-408
      IF (MLX.EQ.1) GO TO 38                                            EXTP-409
      IF (VAL(NMC+3).EQ.0.D0) VAL(NMC+3)=VAL(NMC)                       EXTP-410
      IF (VAL(NMC+4).EQ.0.D0) VAL(NMC+4)=WV(2,L1)                       EXTP-411
      IF (VAL(NMC+5).EQ.0.D0) VAL(NMC+5)=1.D0                           EXTP-412
   38 NMB=NMC+3*MLX-1                                                   EXTP-413
      IF (.NOT.LO(73)) WRITE (MW,1021) NMC,NMB,(VAL(I),I=NMC,NMB)       EXTP-414
      GO TO 40                                                          EXTP-415
   39 NMB=NMC+10                                                        EXTP-416
      READ (MR,1003,ERR= 59 ) (VAL(I),I=NMC,NMB)                        EXTP-417
      IF (VAL(NMC+1).EQ.0.D0) VAL(NMC+1)=WV(2,L1)                       EXTP-418
      IF (VAL(NMC+2).EQ.0.D0) VAL(NMC+2)=1.D0                           EXTP-419
      IF (VAL(NMC+4).EQ.0.D0) VAL(NMC+4)=35.D0                          EXTP-420
      IF (.NOT.LO(73)) WRITE (MW,1022) NMC,NMB,(VAL(I),I=NMC,NMB)       EXTP-421
      IF (MLX.EQ.1.AND.MLY.EQ.2.AND.VAL(NMC-10).LT.0.D0.AND.VAL(NMC).LT.EXTP-422
     10.D0) GO TO 67                                                    EXTP-423
      NMC=NMB+1                                                         EXTP-424
      MLX=MLX-2                                                         EXTP-425
      IF (MLX.EQ.1) GO TO 37                                            EXTP-426
      IF (MLX.EQ.0) GO TO 39                                            EXTP-427
   40 LL(2,ITYP,LZ)=NMB                                                 EXTP-428
      NMA=NMB+1                                                         EXTP-429
   41 IF ((.NOT.LO(10)).OR.LZ.GT.NPP) GO TO 10                          EXTP-430
      NMC=LL(1,ITYP,LZ)                                                 EXTP-431
      EX=VAL(NMC)                                                       EXTP-432
      IF (IPP(1,2,LZ).NE.0.AND.EX.EQ.0.D0.AND.ITYP.EQ.2) GO TO 68       EXTP-433
      IF (IPP(2,2,LZ).NE.0.AND.EX.EQ.0.D0.AND.ITYP.EQ.4) GO TO 69       EXTP-434
      IF (PIP(15,LZ).NE.0.D0.AND.EX.EQ.0.D0.AND.ITYP.EQ.1) GO TO 70     EXTP-435
      IF (PIP(7,LZ).NE.0.D0.AND.EX.EQ.0.D0.AND.ITYP.EQ.5) GO TO 71      EXTP-436
      IF (PIP(8,LZ).NE.0.D0.AND.EX.EQ.0.D0.AND.ITYP.EQ.6) GO TO 72      EXTP-437
      GO TO 10                                                          EXTP-438
   42 NVAL(2)=NFOLT                                                     EXTP-439
      NVAL(3)=NMA                                                       EXTP-440
      NVAL(4)=NMA                                                       EXTP-441
      IF (NFOLT.EQ.0) RETURN                                            EXTP-442
      IF (.NOT.LO(73)) WRITE (MW,1023) NFOLD                            EXTP-443
      DO 43 I=1,NFOLT                                                   EXTP-444
      READ (MR,1003,ERR= 77 ) VAL(NMA),VAL(NMA+1),VAL(NMA+2)            EXTP-445
      IF (.NOT.LO(73)) WRITE (MW,1024) I,VAL(NMA),VAL(NMA+1),VAL(NMA+2) EXTP-446
   43 NVAL(4)=NVAL(4)+3                                                 EXTP-447
      RETURN                                                            EXTP-448
   44 WRITE (MW,1025)                                                   EXTP-449
      GO TO 73                                                          EXTP-450
   45 WRITE (MW,1026)                                                   EXTP-451
      GO TO 73                                                          EXTP-452
   46 WRITE (MW,1027) ITYP,NST,NCOLT                                    EXTP-453
      GO TO 73                                                          EXTP-454
   47 WRITE (MW,1028) ML,L1,L2,NCOLT                                    EXTP-455
      GO TO 73                                                          EXTP-456
   48 WRITE (MW,1029) ML,L1,L2,NCOLL                                    EXTP-457
      GO TO 73                                                          EXTP-458
   49 WRITE (MW,1030) ML,L1,L2                                          EXTP-459
      GO TO 73                                                          EXTP-460
   50 WRITE (MW,1031) L1,L2,ML                                          EXTP-461
      GO TO 73                                                          EXTP-462
   51 WRITE (MW,1032) L1,L2,ML,ITYP                                     EXTP-463
      GO TO 73                                                          EXTP-464
   52 WRITE (MW,1033) L1,L2,ML,ITYP,L1X,L2X,MLX,ITYPX                   EXTP-465
      GO TO 73                                                          EXTP-466
   53 WRITE (MW,1034) MLX,L1X,L2X,NCOLT                                 EXTP-467
      GO TO 73                                                          EXTP-468
   54 WRITE (MW,1035) MLX,L1X,L2X,NCOLL                                 EXTP-469
      GO TO 73                                                          EXTP-470
   55 WRITE (MW,1036) MLX,L1X,L2X,NCOLL                                 EXTP-471
      GO TO 73                                                          EXTP-472
   56 WRITE (MW,1037) MLX,L1X,L2X,ITYPX                                 EXTP-473
      GO TO 73                                                          EXTP-474
   57 WRITE (MW,1038) NVAL(NNA+8),NVAL(NNA+12),K                        EXTP-475
      GO TO 73                                                          EXTP-476
   58 WRITE (MW,1039) ITYPX,IQ(2,LYY),IQ(2,LY)                          EXTP-477
      GO TO 73                                                          EXTP-478
   59 WRITE (MW,1040) ITYP,LZ                                           EXTP-479
      GO TO 73                                                          EXTP-480
   60 K=-NVAL(2*NMA)                                                    EXTP-481
      WRITE (MW,1041) K                                                 EXTP-482
      GO TO 73                                                          EXTP-483
   61 WRITE (MW,1042) ITYPX                                             EXTP-484
      GO TO 73                                                          EXTP-485
   62 WRITE (MW,1043) ITYPX,ITYP                                        EXTP-486
      GO TO 73                                                          EXTP-487
   63 WRITE (MW,1044) ITYPX,ITYP                                        EXTP-488
      GO TO 73                                                          EXTP-489
   64 WRITE (MW,1045) ITYPX,ITYP                                        EXTP-490
      GO TO 73                                                          EXTP-491
   65 WRITE (MW,1046) ITYPX                                             EXTP-492
      GO TO 73                                                          EXTP-493
   66 WRITE (MW,1047) MLX,ITYPX                                         EXTP-494
      GO TO 73                                                          EXTP-495
   67 WRITE (MW,1048) VAL(NMC-10),VAL(NMC)                              EXTP-496
      GO TO 73                                                          EXTP-497
   68 WRITE (MW,1049) LZ,IPP(1,2,LZ),EX                                 EXTP-498
      GO TO 73                                                          EXTP-499
   69 WRITE (MW,1050) LZ,IPP(2,2,LZ),EX                                 EXTP-500
      GO TO 73                                                          EXTP-501
   70 WRITE (MW,1051) LZ,LZ,PIP(15,LZ),EX                               EXTP-502
      GO TO 73                                                          EXTP-503
   71 WRITE (MW,1052) LZ,LZ,PIP(7,LZ),EX                                EXTP-504
      GO TO 73                                                          EXTP-505
   72 WRITE (MW,1053) LZ,LZ,PIP(8,LZ),EX                                EXTP-506
   73 IF (NTOT.EQ.0) GO TO 78                                           EXTP-507
      WRITE (MW,1054) NTOT,NPP                                          EXTP-508
      DO 76 J=1,NPX                                                     EXTP-509
      DO 75 I=1,8                                                       EXTP-510
      IF (LL(1,I,J).NE.-1) GO TO 75                                     EXTP-511
      IF (J.GT.NPP) GO TO 74                                            EXTP-512
      WRITE (MW,1055) (AA(K,I),K=1,3),J                                 EXTP-513
      GO TO 75                                                          EXTP-514
   74 WRITE (MW,1056) (AA(K,I),K=1,3),J-NPP                             EXTP-515
   75 CONTINUE                                                          EXTP-516
   76 CONTINUE                                                          EXTP-517
      GO TO 78                                                          EXTP-518
   77 WRITE (MW,1057) I,NFOLD                                           EXTP-519
   78 WRITE (MW,1058)                                                   EXTP-520
      STOP                                                              EXTP-521
 1000 FORMAT (12I5)                                                     EXTP-522
 1001 FORMAT (/3A8,' POTENTIAL NR(',I2,')')                             EXTP-523
 1002 FORMAT (/3A8,' TRANSITION POTENTIAL FROM LEVEL(',I2,') TO LEVEL(',EXTP-524
     1I2,') AND THE ORDER ML =',I2)                                     EXTP-525
 1003 FORMAT (7F10.5)                                                   EXTP-526
 1004 FORMAT (' THE ELASTIC POTENTIAL NR(',I2,') TYP(',I1,') IS COPIED TEXTP-527
     1O NR(',I2,')')                                                    EXTP-528
 1005 FORMAT (' TRANSITION POTENTIAL L1(',I2,') L2(',I2,') ML(',I2,') TYEXTP-529
     1P(',I2,') IS COPIED TO L1(',I2,') L2(',I2,') ML(',I2,')')         EXTP-530
 1006 FORMAT (' USING PARAMETERS',I6,' TO',I6,' FOR COPY TO:',I5,5X,'(COEXTP-531
     1PY):',I3,5X,'FOLD =',I3,5X,'FROM:',I5,5X,'START:',I5/' USING PARAMEXTP-532
     2ETER',I6,' FOR MULTIPLICATIVE FACTOR',D15.8)                      EXTP-533
 1007 FORMAT (F10.5,F20.10,F10.5,F20.10,A4)                             EXTP-534
 1008 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:'/I5,5X,'NB POIEXTP-535
     1NTS:',I5,5X,'FOLD:',I2,5X,'INTG:',I2,5X,'STEP:',I2,5X,'START:',I5/EXTP-536
     25X,'PARAMETERS',I6,' TO',I6,' TO STORE:'/20X,'STRENGTH:',F15.6,10XEXTP-537
     3,'SCALE:',F15.6/' AND PARAMETERS',I6,' TO',I6,' TO STORE:'/4(6X,'REXTP-538
     4ADIUS',7X,'POTENTIAL',2X)/(4(2X,0P,F10.5,3X,1P,D15.7)))           EXTP-539
 1009 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:'/I5,3X,'TYPE:'EXTP-540
     1,I2,3X,'FOLD:',I2,3X,'INTG:',I2,3X,'STEP:',I2,I6,' FUNCTIONS  L:',EXTP-541
     2I3,I4,' DERIVATIONS  START:',I5)                                  EXTP-542
 1010 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:'/I5,3X,'TYPE:'EXTP-543
     1,I2,3X,'FOLD:',I2,3X,'INTG:',I2,3X,'STEP:',I2,I6,' DEFORMATIONS   EXTP-544
     2LBET:',I5,5X,'KBET:',I5,5X,'START:',I5)                           EXTP-545
 1011 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:'/I5,3X,'TYPE:'EXTP-546
     1,I2,3X,'FOLD:',I2,3X,'INTG:',I2,3X,'STEP:',I2,I6,' DEFORMATIONS',5EXTP-547
     2X,'START:',I5)                                                    EXTP-548
 1012 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:',I4,3X,'TYPE:'EXTP-549
     1,I2,3X,'FOLD:',I2,3X,'INTG:',I2,3X,'STEP:',I2,5X,'START:',I5)     EXTP-550
 1013 FORMAT (' VALUES READ:',F12.6,3X,F9.6)                            EXTP-551
 1014 FORMAT (' USING PARAMETERS',I6,' TO',I6,' FOR'/' DEPTH',F12.6,' MEEXTP-552
     1V  RADIUS',F10.6,' F  DIFFUSENESS',F9.6,' F AT THE POWER (1.+',F9.EXTP-553
     26,')',2X,F9.6,' (3RD COUL. PARM.)'/)                              EXTP-554
 1015 FORMAT (' DEFORMATIONS READ:',8F10.5/(19X,8F10.5))                EXTP-555
 1016 FORMAT (' USING PARAMETERS',I6,' TO',I6,' FOR STRENGTH:',F12.6,10XEXTP-556
     1,'SCALE:',F12.6/' AND TO STORE',I5,' BESSEL FUNCTION OR LEGENDRE PEXTP-557
     2OLYNOMIALS COEFFICIENTS:'/(6(3X,I5,F10.5)))                       EXTP-558
 1017 FORMAT (' USING PARAMETERS',I6,' TO',I6,' FOR',I5,' DEFORMATIONS:'EXTP-559
     1/(6(3X,I5,F10.5)))                                                EXTP-560
 1018 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:',I4,3X,'TYPE:'EXTP-561
     1,I2,3X,'FOLD:',I2,3X,'INTG:',I2,3X,'STEP:',I2,5X,'MULT:',I3)      EXTP-562
 1019 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:',I2,' FUNCTIONEXTP-563
     1S   N =',I2,3X,'L =',I3,3X,'2*S =',I2,3X,'2*J =',I3,5X,'START:',I5EXTP-564
     2)                                                                 EXTP-565
 1020 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:'/I5,' FUNCTIONEXTP-566
     1S',2(2X,'N =',I2,3X,'L =',I3,3X,'2*S =',I2,3X,'2*J =',I3),5X,'STAREXTP-567
     2T:',I5/' USING PARAMETER',I6,' FOR MULTIPLICATIVE FACTOR',D18.8)  EXTP-568
 1021 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:'/(' OSCILLATOREXTP-569
     1 PARAMETER',F10.6,2X,'TOTAL MASS',F12.6,2X,'PARTICLE MASS',F12.6))EXTP-570
 1022 FORMAT (' USING PARAMETERS',I6,' TO',I6,' TO STORE:'/' **** BINDINEXTP-571
     1G ENERGY',F12.6,' MEV ****',2X,'TOTAL MASS',F12.6,2X,'PARTICLE MASEXTP-572
     2S',F12.6,2X,'PRODUCT OF CHARGES',F8.2/' SEARCH ON DEPTH OF REAL POEXTP-573
     3TENTIAL FROM',F12.6,' WITH REDUCED RADIUS',F10.6,' FERMI AND DIFFUEXTP-574
     4SENESS',F10.6,' FERMI'/' SPIN-ORBIT POTENTIAL  DEPTH:',F12.6,' MEVEXTP-575
     5  RADIUS:',F10.6,' F  DIFFUSENESS:',F9.6,' F  COULOMB RADIUS:',F10EXTP-576
     6.6,' F'/)                                                         EXTP-577
 1023 FORMAT (/2X,I3,' SETS OF FOLDING PARAMETERS')                     EXTP-578
 1024 FORMAT (2X,I2,5X,' V =',F12.6,5X,' R =',F9.6,5X,' A =',F9.6)      EXTP-579
 1025 FORMAT (' EXTERNAL FORM FACTORS NOT ALLOWED WITH SCHROEDINDER EQUIEXTP-580
     1VALENT OF DIRAC EQUATION')                                        EXTP-581
 1026 FORMAT (' INPUT ERROR IN THE FIRST CARD DEFINING THE POTENTIAL')  EXTP-582
 1027 FORMAT (/' ITYP =',I5,' NOT ALLOWED OR NST =',I5,' LARGER THAN NCOEXTP-583
     1LT =',I3)                                                         EXTP-584
 1028 FORMAT (' WITH ML =',I2,' L1 =',I3,' IS NOT EQUAL TO L2 =',I3,' OREXTP-585
     1 IS LARGER THAN NCOLT =',I3)                                      EXTP-586
 1029 FORMAT (' WITH ML =',I3,' L1 =',I3,' OR L2 =',I3,' IS LARGER THAN EXTP-587
     1NCOLL =',I3)                                                      EXTP-588
 1030 FORMAT (' ML =',I3,' TOO LARGE BETWEEN LEVELS L1 =',I3,' AND L2 ='EXTP-589
     1,I3)                                                              EXTP-590
 1031 FORMAT (' THERE IS NO CORRECTION TERM OF THE FORM FACTOR FOR L1 ='EXTP-591
     1,I3,'  L2 =',I3,'  ML =',I3)                                      EXTP-592
 1032 FORMAT (' THE FORM FACTOR FOR L1 =',I3,'  L2 =',I3,'  ML =',I3,'  EXTP-593
     1AND ITYP =',I2,' IS ALREADY READ OR DOES NOT HAVE TO BE READ')    EXTP-594
 1033 FORMAT (/' DIFFERENT TYPES NOT ALLOWED TO COPY FOR FORM FACTORS L1EXTP-595
     1/L2/ML/ITYP/L1X/L2X/MLX/ITYPX :'/30X,8I5)                         EXTP-596
 1034 FORMAT (' WITH MLX =',I2,' L1X =',I3,' IS NOT EQUAL TO L2X =',I3,'EXTP-597
     1 OR IS LARGER THAN NCOLT =',I3)                                   EXTP-598
 1035 FORMAT (' WITH MLX =',I3,' L1X =',I3,' OR L2X =',I3,' IS LARGER THEXTP-599
     1AN NCOLL =',I3)                                                   EXTP-600
 1036 FORMAT (' MLX =',I3,' TOO LARGE BETWEEN LEVELS L1X =',I3,' AND L2XEXTP-601
     1 =',I5)                                                           EXTP-602
 1037 FORMAT (' FORM FACTOR DEFINED BY L1X =',I3,'  L2X =',I3,'  MLX =',EXTP-603
     1I3,'  ITYPX =',I3,' NOT YET DEFINED')                             EXTP-604
 1038 FORMAT (' NO PARTICLE-HOLE COUPLING WITH 2*JP =',I3,' AND 2*JH =',EXTP-605
     1I3,' TO L =',I3)                                                  EXTP-606
 1039 FORMAT (' COPY WITH ITYPX =',I3,' NOT ALLOWED BETWEEN TRANSITIONS EXTP-607
     1WITH QUANTUM NUMBERS',I2,' AND',I2)                               EXTP-608
 1040 FORMAT (' INPUT ERROR FOR THE POTENTIAL (',I2,',',I2,')')         EXTP-609
 1041 FORMAT (' NUMBER OF DATA',I4,' LESS THAN 4 WHICH ARENECESSARY TO IEXTP-610
     1NTERPOLATE.')                                                     EXTP-611
 1042 FORMAT (/' ITYPX =',I5,' NOT ALLOWED FOR STANDARD FORM FACTORS')  EXTP-612
 1043 FORMAT (' ITYPX,ITYP =',2I5,'  NOT ALLOWED: CENTRAL POTENTIAL CANNEXTP-613
     1OT BE DERIVATIVE OR BOUND FUNCTION')                              EXTP-614
 1044 FORMAT (' ITYPX,ITYP =',2I5,'  NOT ALLOWED: COULOMB TRANSITION FOREXTP-615
     1M FACTOR MUST BE DERIVATIVE')                                     EXTP-616
 1045 FORMAT (' ITYPX,ITYP =',2I5,'  NOT ALLOWED: BOUND STATE FUNCTION CEXTP-617
     1AN BE ONLY REAL OR IMAGINARY TRANSITION POTENTIAL')               EXTP-618
 1046 FORMAT (' ITYPX =',I3,'  NOT ALLOWED IN DIRAC FORMALISM')         EXTP-619
 1047 FORMAT (' L2X =',I3,'  NOT ALLOWED WITH ITYPX =',I4)              EXTP-620
 1048 FORMAT (' FOR ITYP=-8 AND TWO FUNCTIONS, BOTH ARE UNBOUNDED',2F12.EXTP-621
     16)                                                                EXTP-622
 1049 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED FOR DISPERSION RELATIEXTP-623
     1ONS BECAUSE THE VOLUME IMAGINARY STRENGTH IS 0:'/10X,'NV =',I3,'  EXTP-624
     1WV =',D13.6)                                                      EXTP-625
 1050 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED FOR DISPERSION RELATIEXTP-626
     1ONS BECAUSE THE SURFACE IMAGINARY STRENGTH IS 0:'/10X,'NS =',I3,' EXTP-627
     2 WS =',D13.6)                                                     EXTP-628
 1051 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED WITH VARIATION OF THEEXTP-629
     1 HARTREE-FOCK POTENTIAL OF WHICH THE STRENGTH IS 0':/10X,'PIP(15,'EXTP-630
     2,I2,') =',D13.6,'  V =',D13.6)                                    EXTP-631
 1052 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED WITH VARIATION OF THEEXTP-632
     1 REAL SPIN-ORBIT OF WHICH THE STRENGTH IS 0':/10X,'PIP(7,',I2,') =EXTP-633
     2',D13.6,'  VLS =',D13.6)                                          EXTP-634
 1053 FORMAT (' THE POTENTIAL',I4,' CANNOT BE USED WITH VARIATION OF THEEXTP-635
     1 IMAGINARY SPIN-ORBIT OF WHICH THE STRENGTH IS 0':/10X,'PIP(8,',I2EXTP-636
     2,') =',D13.6,'  WLS =',D13.6)                                     EXTP-637
 1054 FORMAT (I6,' FORM FACTORS TO READ, FOR (ITYP, N  )  ( POTENTIALS FEXTP-638
     1OR N SMALLER THAN',I2,' TRANSITIONS AFTER)')                      EXTP-639
 1055 FORMAT (10X,3A8,' MISSING FOR POTENTIAL',I4)                      EXTP-640
 1056 FORMAT (10X,3A8,' MISSING FOR TRANSITION',I4)                     EXTP-641
 1057 FORMAT (' INPUT ERROR FOR THE',I4,'TH SET OF THE',I4,' SETS OF FOLEXTP-642
     1DING PARAMETERS TO BE READ')                                      EXTP-643
 1058 FORMAT (//' IN EXTP  ...  STOP  ...')                             EXTP-644
      END                                                               EXTP-645
