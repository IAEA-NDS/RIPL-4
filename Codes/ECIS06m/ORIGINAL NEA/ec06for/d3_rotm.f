C 04/07/06                                                      ECIS06  ROTM-000
      SUBROUTINE ROTM(NIV,IQ,T,IPI,NCOLL,IPH,NBETA,NVAR,VAR,IT,FAC,IDT,LROTM-001
     1O)                                                                ROTM-002
C NUCLEAR REDUCED MATRIX ELEMENTS FOR THE SYMMETRIC ROTATIONAL MODEL.   ROTM-003
C INPUT:     IPI(J,*):PARITY OF NUCLEAR STATES FOR J=1,                 ROTM-004
C                     MULTIPLICITY OF THE PARTICLE FOR J=2,             ROTM-005
C                     MULTIPLICITY OF THE TARGET FOR J=3.               ROTM-006
C            NCOLL:   NUMBER OF COUPLED STATES                          ROTM-007
C            IPH:     DESCRIPTION OF VIBRATIONAL MODEL (SEE VIBM).      ROTM-008
C            NBETA:   L QUANTUM NUMBERS OF DEFORMATIONS IN NBETA(17,*), ROTM-009
C                     K QUANTUM NUMBERS OF DEFORMATIONS IN NBETA(18,*),.ROTM-010
C            NVAR,VAR:EQUIVALENT BY CALL, MIXTURE COEFFICIENTS OF BANDS.ROTM-011
C                     MODEL.                                            ROTM-012
C            FAC:     TABLE AND NUMBER OF LOGARITHMS OF FACTORIALS.     ROTM-013
C            IDT:     LENGTH AVAILABLE IN THIS SUBROUTINE.              ROTM-014
C            LO(I):   LOGICAL CONTROLS:                                 ROTM-015
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      ROTM-016
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   ROTM-017
C               LO(117)=.TRUE. FOR ALL CALCULATIONS EXCEPT THE FIRST.   ROTM-018
C OUTPUT:    NIV:     NUMBER AND REFERENCES OF INTERACTION FORM FACTORS ROTM-019
C                     BETWEEN EACH CHANNELS. NIV(I1,I2,K): FIRST I OF   ROTM-020
C                     T(3,I) FOR THE PAIR OF NUCLEAR STATES I1,I2 FOR   ROTM-021
C                     K=1 AND LAST ONE FOR K=2.                         ROTM-022
C            IQ,T:    EQUIVALENT BY CALL OF IQ(6,I) AND T(3,I).         ROTM-023
C            IT:      NUMBER OF NUCLEAR MATRIX ELEMENTS IN IQ,T. FOR A  ROTM-024
C                     GIVEN VALUE OF I:                                 ROTM-025
C                     IQ(1,I): REFERENCE TO TABLE OF FORM FACTORS,      ROTM-026
C                              THE CONTROL NUMBER FOR THE FORM FACTOR   ROTM-027
C                              WITH THE VIBRATION N1 AND THE            ROTM-028
C                              MULTIPOLARITY L IS 1000*(L+1)+N1.        ROTM-029
C                     IQ(2,I): REFERENCE TO TABLE OF ANGULAR MOMENTA,   ROTM-030
C                     IQ(3,I): ADDRESS OF THE ASSOCIATED SPIN-ORBIT FORMROTM-031
C                              FACTOR OR 0,                             ROTM-032
C                     IQ(4,I): UNUSED,                                  ROTM-033
C                     T(3,I):  REDUCED NUCLEAR MATRIX ELEMENT:          ROTM-034
C                     (IP||Q||I)=SQRT(2*I+1)*CG(I,IQ,IP,K,0,K),         ROTM-035
C                     (IP+V||Q||I)=SQRT(2*I+1)*CG(I,IQ,IP,K,V,K+V),     ROTM-036
C                     (IP+V||Q||I+V)=SQRT(2*I+1)*CG(I,IQ,IP,K+V,0,K+V). ROTM-037
C                                                                       ROTM-038
C FOR THE COMMON  /COUPL/ SEE CALC.                                     ROTM-039
C                                                                       ROTM-040
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     ROTM-041
C  IQMAX:     MAXIMUM L-EXPANSION IN ROTATIONAL MODEL.                  ROTM-042
C  NFA:       NUMBER OF LOGARITHMS OF FACTORIALS.                       ROTM-043
C  NSPIN:     TWICE THE K-VALUE OF THE ROTATIONAL BAND.                 ROTM-044
C  NVA:       NUMBER OF NUCLEAR PARAMETERS.                             ROTM-045
C   USED:     IQMAX,NFA,NSPIN,NVA.                                      ROTM-046
C                                                                       ROTM-047
C***********************************************************************ROTM-048
      IMPLICIT REAL*8 (A-H,O-Z)                                         ROTM-049
      LOGICAL LO(150)                                                   ROTM-050
      DIMENSION NIV(NCOLL,NCOLL,*),IQ(6,*),T(3,*),IPI(11,*),IPH(2,*),NBEROTM-051
     1TA(18,*),NVAR(2,*),VAR(*),FAC(*),AA(2,2)                          ROTM-052
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   ROTM-053
      COMMON /INOUT/ MR,MW,MS                                           ROTM-054
      IT=1                                                              ROTM-055
      NSP=0                                                             ROTM-056
      IF (LO(13).OR.LO(19)) NSP=10                                      ROTM-057
      DO 24 I1=1,NCOLL                                                  ROTM-058
      IA1=IPI(3,I1)-1                                                   ROTM-059
      AA(1,1)=1.D0                                                      ROTM-060
      AA(1,2)=0.D0                                                      ROTM-061
      IF (IPH(1,I1).LE.1) GO TO 1                                       ROTM-062
      JVAR=IPH(2,I1)                                                    ROTM-063
      IVAR=NVAR(2,JVAR)                                                 ROTM-064
      AY=1.74532925D-02*VAR(IVAR)                                       ROTM-065
      AA(1,1)=DCOS(AY)                                                  ROTM-066
      AA(1,2)=DSIN(AY)                                                  ROTM-067
      IF (.NOT.LO(117)) WRITE (MW,1000) I1,VAR(IVAR),AA(1,1),AA(1,2)    ROTM-068
    1 AA(2,1)=AA(1,1)                                                   ROTM-069
      AA(2,2)=AA(1,2)                                                   ROTM-070
      DO 23 I2=I1,NCOLL                                                 ROTM-071
      IA2=IPI(3,I2)-1                                                   ROTM-072
      IF (I1.EQ.I2) GO TO 2                                             ROTM-073
      AA(2,1)=1.D0                                                      ROTM-074
      AA(2,2)=0.D0                                                      ROTM-075
      IF (IPH(1,I2).LE.1) GO TO 2                                       ROTM-076
      IVAR=IPH(2,I2)                                                    ROTM-077
      JVAR=NVAR(2,IVAR)                                                 ROTM-078
      IF (JVAR.GT.NVA) GO TO 25                                         ROTM-079
      AY=1.74532925D-02*VAR(JVAR)                                       ROTM-080
      AA(2,1)=DCOS(AY)                                                  ROTM-081
      AA(2,2)=DSIN(AY)                                                  ROTM-082
    2 NIV(I1,I2,1)=IT                                                   ROTM-083
      L1=IPH(1,I1)                                                      ROTM-084
      AX=AA(1,1)                                                        ROTM-085
      IF (L1.GT.1) L1=0                                                 ROTM-086
      IF (DABS(AX).LT.1.D-6) GO TO 12                                   ROTM-087
    3 L2=IPH(1,I2)                                                      ROTM-088
      IF (L2.GT.1) L2=0                                                 ROTM-089
      AY=AX*AA(2,1)                                                     ROTM-090
      IF (DABS(AY).LT.1.D-6) GO TO 11                                   ROTM-091
    4 IF (L1.NE.L2) GO TO 7                                             ROTM-092
      NSJ=NSPIN                                                         ROTM-093
      IF (L1.EQ.0) GO TO 5                                              ROTM-094
      N3=IPH(2,I1)                                                      ROTM-095
      IF (N3.NE.IPH(2,I2)) GO TO 11                                     ROTM-096
      NSJ=2*NBETA(18,N3)+NSJ                                            ROTM-097
    5 IF (IPI(1,I1).NE.IPI(1,I2)) GO TO 26                              ROTM-098
      IF ((IABS(NSJ).GT.IA2).OR.(IABS(NSJ).GT.IA1)) GO TO 27            ROTM-099
      IQ1=MIN0(IQMAX,(IA1+IA2)/2)                                       ROTM-100
      IQ2=MAX0(2,IABS(IA1-IA2)/2)                                       ROTM-101
      IF (2*(IQ2/2).NE.IQ2) IQ2=IQ2+1                                   ROTM-102
      IF (IQ2.GT.IQ1) GO TO 11                                          ROTM-103
      DO 6 IQZ=IQ2,IQ1,2                                                ROTM-104
      IF (2*IT.GT.IDT) CALL MEMO('ROTM',IDT,2*IT)                       ROTM-105
      IQ(1,IT)=1000*(IQZ+1)                                             ROTM-106
      IQ(2,IT)=IQZ                                                      ROTM-107
      IQ(3,IT)=NSP                                                      ROTM-108
      T(3,IT)=DFLOAT(1-MOD(IQZ,4))*DSQRT(IA1+1.D0)*DJCG(IA1,2*IQZ,IA2,NSROTM-109
     1J,0,FAC,NFA)*AY                                                   ROTM-110
      IF (DABS(T(3,IT)).GT.1.D-6) IT=IT+1                               ROTM-111
    6 CONTINUE                                                          ROTM-112
      GO TO 11                                                          ROTM-113
    7 IF (L1.GT.L2) GO TO 8                                             ROTM-114
      N3=IPH(2,I2)                                                      ROTM-115
      GO TO 9                                                           ROTM-116
C  TRANSPOSITION.                                                       ROTM-117
    8 AY=AY*DFLOAT(1-MOD(IPI(3,I1)+IPI(3,I2)+2*(IPI(1,I1)+IPI(1,I2)+1),4ROTM-118
     1))                                                                ROTM-119
      N3=IPH(2,I1)                                                      ROTM-120
C 0 PHONONS -1 PHONON.                                                  ROTM-121
    9 NSJ=2*NBETA(18,N3)                                                ROTM-122
      IF (MOD(IPI(1,I1)+IPI(1,I2)+NBETA(17,N3),2).NE.0) GO TO 28        ROTM-123
      IF (IABS(NSJ+NSPIN).GT.IA2) GO TO 27                              ROTM-124
      IQ1=(IA1+IA2)/2                                                   ROTM-125
      IQ2=MAX0(IABS(NSJ),IABS(IA1-IA2))/2                               ROTM-126
      IF (MOD(IQ2+NBETA(17,N3),2).NE.0) IQ2=IQ2+1                       ROTM-127
      IF (IQ2.GT.IQ1) GO TO 11                                          ROTM-128
      DO 10 IQZ=IQ2,IQ1,2                                               ROTM-129
      IF (3*IT.GT.IDT) CALL MEMO('ROTM',IDT,3*IT)                       ROTM-130
      IQ(1,IT)=1000*(IQZ+1)+N3                                          ROTM-131
      IQ(2,IT)=IQZ                                                      ROTM-132
      IQ(3,IT)=NSP                                                      ROTM-133
      T(3,IT)=DFLOAT(1-MOD(IQZ+IPI(1,I1)+IPI(1,I2),4))*DSQRT(IA1+1.D0)*DROTM-134
     1JCG(IA1,2*IQZ,IA2,NSPIN,NSJ,FAC,NFA)*AY                           ROTM-135
      IF (DABS(T(3,IT)).GT.1.D-6) IT=IT+1                               ROTM-136
   10 CONTINUE                                                          ROTM-137
   11 IF ((AA(2,2).EQ.0.D0).OR.(L2.EQ.1)) GO TO 12                      ROTM-138
      L2=1                                                              ROTM-139
      AY=AX*AA(2,2)                                                     ROTM-140
      GO TO 4                                                           ROTM-141
   12 IF ((AA(1,2).EQ.0.D0).OR.(L1.EQ.1)) GO TO 13                      ROTM-142
      L1=1                                                              ROTM-143
      AX=AA(1,2)                                                        ROTM-144
      GO TO 3                                                           ROTM-145
   13 ITI=NIV(I1,I2,1)                                                  ROTM-146
      IF (ITI.EQ.IT) GO TO 22                                           ROTM-147
      ITF=IT-1                                                          ROTM-148
      IF (ITF.EQ.ITI) GO TO 17                                          ROTM-149
      IT1=ITI+1                                                         ROTM-150
      IT=ITI                                                            ROTM-151
      DO 16 I=IT1,ITF                                                   ROTM-152
      DO 14 J=ITI,IT                                                    ROTM-153
      IF ((IQ(1,I).NE.IQ(1,J)).OR.(IQ(2,I).NE.IQ(2,J)).OR.(IQ(3,I).NE.IQROTM-154
     1(3,J))) GO TO 14                                                  ROTM-155
      T(3,J)=T(3,J)+T(3,I)                                              ROTM-156
      GO TO 16                                                          ROTM-157
   14 CONTINUE                                                          ROTM-158
      IT=IT+1                                                           ROTM-159
      DO 15 K=1,3                                                       ROTM-160
   15 IQ(K,IT)=IQ(K,I)                                                  ROTM-161
      T(3,IT)=T(3,I)                                                    ROTM-162
   16 CONTINUE                                                          ROTM-163
      IT=IT+1                                                           ROTM-164
   17 ITF=IT-1                                                          ROTM-165
      DO 21 I=ITI,ITF                                                   ROTM-166
      IF (IT.LE.I) GO TO 22                                             ROTM-167
   18 IF (DABS(T(3,I)).GT.1.D-12) GO TO 21                              ROTM-168
      IT1=I+1                                                           ROTM-169
      IT=IT-1                                                           ROTM-170
      IF (IT1.GT.IT) GO TO 22                                           ROTM-171
      DO 20 J=IT1,IT                                                    ROTM-172
      DO 19 K=1,3                                                       ROTM-173
   19 IQ(K,J-1)=IQ(K,J)                                                 ROTM-174
   20 T(3,J-1)=T(3,J)                                                   ROTM-175
      GO TO 18                                                          ROTM-176
   21 CONTINUE                                                          ROTM-177
   22 NIV(I1,I2,2)=IT-1                                                 ROTM-178
   23 CONTINUE                                                          ROTM-179
   24 CONTINUE                                                          ROTM-180
      IT=IT-1                                                           ROTM-181
      RETURN                                                            ROTM-182
   25 WRITE (MW,1001) JVAR,NVA                                          ROTM-183
      GO TO 29                                                          ROTM-184
   26 WRITE (MW,1002) I1,I2                                             ROTM-185
      GO TO 29                                                          ROTM-186
   27 WRITE (MW,1003) I1,I2                                             ROTM-187
      GO TO 29                                                          ROTM-188
   28 WRITE (MW,1004) I1,I2,N3                                          ROTM-189
      GO TO 29                                                          ROTM-190
   29 WRITE (MW,1005)                                                   ROTM-191
      STOP                                                              ROTM-192
 1000 FORMAT (' STATE',I4,F15.5,' DEGREES      AMPLITUDES =',F15.7,' GROROTM-193
     1UND STATE BAND AND',F15.5,' VIBRATIONAL BAND.')                   ROTM-194
 1001 FORMAT (' NUMBER OF VARIABLES USED:',I5,5X,'EXCEEDS NUMBER OF VARIROTM-195
     2ABLES READ:',I6)                                                  ROTM-196
 1002 FORMAT (/' PARITIES OF STATES',I4,'  AND',I4,'  INCORRECT FOR THE ROTM-197
     1ROTATIONAL MODEL.')                                               ROTM-198
 1003 FORMAT (' TOO LARGE MAGNETIC QUANTUM NUMBER BETWEEN LEVELS',I4,' AROTM-199
     1ND',I4)                                                           ROTM-200
 1004 FORMAT (/' PARITIES OF STATES',I4,'  AND',I4,'  INCORRECT FOR THE ROTM-201
     1ROTATIONAL MODEL WITH THE VIBRATION',I4)                          ROTM-202
 1005 FORMAT (//' IN ROTM  ...  STOP  ...')                             ROTM-203
      END                                                               ROTM-204
