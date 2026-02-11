C 28/02/07                                                      ECIS06  ROAM-000
      SUBROUTINE ROAM(NIV,IQ,T,IPI,NCOLL,IT,BETA,IPH,VAR,VA,FAC,IDT,LO) ROAM-001
C NUCLEAR REDUCED MATRIX ELEMENTS FOR THE ASYMMETRIC ROTATIONAL MODEL.  ROAM-002
C INPUT:     IPI(J,*):PARITY OF NUCLEAR STATES FOR J=1,                 ROAM-003
C                     MULTIPLICITY OF THE PARTICLE FOR J=2,             ROAM-004
C                     MULTIPLICITY OF THE TARGET FOR J=3.               ROAM-005
C            NCOLL:   NUMBER OF COUPLED STATES                          ROAM-006
C            BETA:    DEFORMATIONS FOR L=2 ARE IN BETA(*,1),            ROAM-007
C                     GAMMA VALUES FOR L=2 ARE IN BETA(*,2),            ROAM-008
C                     DEFORMATIONS FOR L=4 ARE IN BETA(*,3),            ROAM-009
C                     ... ETC (SEE DESCRIPTION OF INPUT).               ROAM-010
C            IPH(I,*):NUMBER OF NUCLEAR PARAMETERS FOR I=1,             ROAM-011
C                     THEIR ADDRESS FOR I=2.                            ROAM-012
C            VAR:     EQUIVALENT BY CALL, PARAMETERS FOR SOME MODELS.   ROAM-013
C            FAC:     TABLE AND NUMBER OF LOGARITHMS OF FACTORIALS.     ROAM-014
C            IDT:     LENGTH AVAILABLE IN THIS SUBROUTINE.              ROAM-015
C            LO(I):   LOGICAL CONTROLS:                                 ROAM-016
C               LO(2)  =.TRUE. SECOND ORDER VIBRATIONAL OR CONSTRAINED  ROAM-017
C                              ASYMMETRIC ROTATIONAL MODEL.             ROAM-018
C               LO(13) =.TRUE. DEFORMED REAL SPIN-ORBIT OR TENSOR.      ROAM-019
C               LO(19) =.TRUE. DEFORMED COULOMB SPIN-ORBIT POTENTIAL.   ROAM-020
C               LO(117)=.TRUE. FOR ALL CALCULATIONS EXCEPT THE FIRST.   ROAM-021
C OUTPUT:    NIV:     NUMBER AND REFERENCES OF INTERACTION FORM FACTORS ROAM-022
C                     BETWEEN EACH CHANNELS. NIV(I1,I2,K): FIRST I OF   ROAM-023
C                     T(3,I) FOR THE PAIR OF NUCLEAR STATES I1,I2 FOR   ROAM-024
C                     K=1 AND LAST ONE FOR K=2.                         ROAM-025
C            IQ,T:    EQUIVALENT BY CALL OF IQ(6,I) AND T(3,I).         ROAM-026
C            IT:      NUMBER OF NUCLEAR MATRIX ELEMENTS IN IQ,T. FOR A  ROAM-027
C                     GIVEN VALUE OF I:                                 ROAM-028
C                     IQ(1,I): REFERENCE TO TABLE OF FORM FACTORS,      ROAM-029
C                              THE CONTROL NUMBER FOR FORM FACTOR       ROAM-030
C                              (L=2,K=0) IS 2,FOR (L=2,K=2) IS 3, FOR   ROAM-031
C                              (L=4,K=0) IS 4,FOR (L=4,K=2) IS 5 ...    ROAM-032
C                              AND SO ON...                             ROAM-033
C                     IQ(2,I): REFERENCE TO TABLE OF ANGULAR MOMENTA,   ROAM-034
C                     IQ(3,I): ADDRESS OF THE ASSOCIATED SPIN-ORBIT FORMROAM-035
C                              FACTOR OR 0,                             ROAM-036
C                     IQ(4,I): UNUSED,                                  ROAM-037
C                     T(3,I):  REDUCED NUCLEAR MATRIX ELEMENT:          ROAM-038
C                     (IP||Q(IQ,KQ)||I)=SQRT(2*I+1)*CG(I,IQ,IP,K,KQ,KP) ROAM-039
C            VA:      BAND MIXING COEFFICIENTS FOR L=2.                 ROAM-040
C                                                                       ROAM-041
C FOR THE COMMON  /COUPL/ SEE CALC.                                     ROAM-042
C                                                                       ROAM-043
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /COUPL/:                     ROAM-044
C  IQMAX:     MAXIMUM L-EXPANSION IN ROTATIONAL MODEL.                  ROAM-045
C  NFA:       NUMBER OF LOGARITHMS OF FACTORIALS.                       ROAM-046
C   USED:     IQMAX,NFA.                                                ROAM-047
C                                                                       ROAM-048
C***********************************************************************ROAM-049
      IMPLICIT REAL*8 (A-H,O-Z)                                         ROAM-050
      LOGICAL LO(150)                                                   ROAM-051
      DIMENSION NIV(NCOLL,NCOLL,*),IQ(6,*),T(3,*),IPI(11,*),BETA(9,*),IPROAM-052
     1H(2,*),VAR(*),VA(*),FAC(*)                                        ROAM-053
      COMMON /COUPL/ IQM,IQMAX,NBT1,NFA,NPP,NSPIN,NVA                   ROAM-054
      COMMON /INOUT/ MR,MW,MS                                           ROAM-055
      IVAR=0                                                            ROAM-056
      I2=1                                                              ROAM-057
      IF (.NOT.LO(2)) GO TO 4                                           ROAM-058
C VAR(1) IS THE ANGLE GAMMA OF DAVYDOV AND FILIPPOV MODEL.              ROAM-059
C THE SECOND LEVEL MUST BE A THE FIRST 2+ AND THE THIRD ONE,IF PRESENT, ROAM-060
C THE SECOND 2+ .THEIR MIXING COEFFICIENTS ARE COMPUTED FROM VAR(1).    ROAM-061
C THERE CAN BE ANY NUMBER OF LEVELS AFTER THE THIRD ONE.                ROAM-062
      G=VAR(1)-BETA(1,2)                                                ROAM-063
      IF (.NOT.LO(117)) WRITE (MW,1000) G,(BETA(I,2),I=1,8)             ROAM-064
      DO 1 I=1,8                                                        ROAM-065
    1 BETA(I,2)=BETA(I,2)+G                                             ROAM-066
      IF (.NOT.LO(117)) WRITE (MW,1001) (BETA(I,2),I=1,8)               ROAM-067
      G=0.0174532925199433D0*VAR(1)                                     ROAM-068
      G1=DSIN(G)                                                        ROAM-069
      G2=DCOS(G)                                                        ROAM-070
      G3=DSIN(3.D0*G)                                                   ROAM-071
      G4=DCOS(3.D0*G)                                                   ROAM-072
      G=DSQRT(9.D0-8.D0*G3*G3)                                          ROAM-073
      G5=-(G1*G3+3.D0*G2*G4+G)                                          ROAM-074
      G6=3.D0*G1*G4-G2*G3                                               ROAM-075
      VA(1)=1.D0                                                        ROAM-076
      IF (NCOLL.EQ.1) GO TO 2                                           ROAM-077
      G=DSQRT(G6*G6+G5*G5)                                              ROAM-078
      VA(2)=G5/G                                                        ROAM-079
      VA(3)=G6/G                                                        ROAM-080
      IF (NCOLL.EQ.2) GO TO 2                                           ROAM-081
      VA(4)=-VA(3)                                                      ROAM-082
      VA(5)=VA(2)                                                       ROAM-083
    2 I2=MIN0(NCOLL,3)                                                  ROAM-084
      DO 3 I1=1,I2                                                      ROAM-085
      IF (I1.EQ.1.AND.IPI(3,I1).NE.1) GO TO 14                          ROAM-086
      IF (I1.EQ.1) GO TO 3                                              ROAM-087
      IF (IPI(3,I1).NE.5) GO TO 14                                      ROAM-088
      IF (.NOT.LO(117)) WRITE (MW,1002) I1,VA(2*I1-2),VA(2*I1-1)        ROAM-089
    3 CONTINUE                                                          ROAM-090
      IVAR=I2-1                                                         ROAM-091
      I2=4                                                              ROAM-092
      IF (I2.GT.NCOLL) GO TO 7                                          ROAM-093
C LEVELS NOT RELATED TO THE GAMMA DEFORMATION.                          ROAM-094
    4 DO 6 I1=I2,NCOLL                                                  ROAM-095
      IF (IPI(3,I1).EQ.3) GO TO 16                                      ROAM-096
      IVAR=IPH(2,I1)                                                    ROAM-097
      IF (MOD(IPI(3,I1),4).NE.1) VAR(IVAR+1)=90.D0                      ROAM-098
      K1=IVAR+I1                                                        ROAM-099
      VA(K1)=1.D0                                                       ROAM-100
      IF (IPH(1,I1).EQ.0) GO TO 6                                       ROAM-101
      K2=K1+IPH(1,I1)-1                                                 ROAM-102
      DO 5 K=K1,K2                                                      ROAM-103
      IVAR=IVAR+1                                                       ROAM-104
      G=0.0174532925199433D0*VAR(IVAR)                                  ROAM-105
      VA(K+1)=VA(K)*DSIN(G)                                             ROAM-106
    5 VA(K)=VA(K)*DCOS(G)                                               ROAM-107
      IF (.NOT.LO(117)) WRITE (MW,1002) I1,(VA(K),K=K1,K2),VA(K2+1)     ROAM-108
    6 CONTINUE                                                          ROAM-109
C COMPUTATION OF RED. MAT. ELE. FOR I=/<IP AT THE FIRST CALL.           ROAM-110
    7 IT=1                                                              ROAM-111
      NSP=0                                                             ROAM-112
      IF (LO(13).OR.LO(19)) NSP=10                                      ROAM-113
      DO 13 I1=1,NCOLL                                                  ROAM-114
      NT=IPH(1,I1)+1                                                    ROAM-115
      NV=IPH(2,I1)+I1-1                                                 ROAM-116
      IA1=IPI(3,I1)-1                                                   ROAM-117
      NX=IA1/2-2*(IA1/4)                                                ROAM-118
      DO 12 I2=I1,NCOLL                                                 ROAM-119
      IF (IPI(1,I1).NE.IPI(1,I2)) GO TO 15                              ROAM-120
      NIV(I1,I2,1)=IT                                                   ROAM-121
      MT=IPH(1,I2)+1                                                    ROAM-122
      MV=IPH(2,I2)+I2-1                                                 ROAM-123
      IA2=IPI(3,I2)-1                                                   ROAM-124
      MX=IA2/2-2*(IA2/4)                                                ROAM-125
      IQ1=MIN0(IQMAX,(IA1+IA2)/2,8)                                     ROAM-126
      IQ2=MAX0(2,2*((IABS(IA1-IA2)/2+1)/2))                             ROAM-127
      IF (IQ2.GT.IQ1) GO TO 11                                          ROAM-128
C COMPUTATION OF REDUCED MATRIX ELEMENTS.                               ROAM-129
      DO 10 IQZ=IQ2,IQ1,2                                               ROAM-130
      FS=DFLOAT(IQZ-4*((IQZ+2)/4)+1)                                    ROAM-131
      IQY=IQZ/2+1                                                       ROAM-132
      DO 9 MQZ=1,IQY                                                    ROAM-133
      IF (2*IT.GT.IDT) CALL MEMO('ROAM',IDT,2*IT)                       ROAM-134
      IQ(1,IT)=1000*((IQY*(IQY-1))/2+MQZ)                               ROAM-135
      B=0.D0                                                            ROAM-136
      IQ(2,IT)=IQZ                                                      ROAM-137
      IQ(3,IT)=NSP                                                      ROAM-138
      MQ=2*MQZ-2                                                        ROAM-139
      DO 8 N1=1,NT                                                      ROAM-140
      IF (N1-NX.EQ.0) GO TO 8                                           ROAM-141
      N=2*N1-2                                                          ROAM-142
      LM=N-MQ                                                           ROAM-143
      IF (IABS(LM).GE.2*MT) GO TO 8                                     ROAM-144
      L1=1+IABS(LM)/2                                                   ROAM-145
      FQ=DJCG(IA1,2*IQZ,IA2,2*N,-2*MQ,FAC,NFA)*VA(NV+N1)*VA(MV+L1)      ROAM-146
      IF (LM.LT.0.AND.MX.EQ.1) FQ=-FQ                                   ROAM-147
      IF (MQ*N*LM.NE.0) FQ=.7071068D0*FQ                                ROAM-148
      B=B+FQ                                                            ROAM-149
      IF (N*MQ.EQ.0) GO TO 8                                            ROAM-150
      LM=MQ+N                                                           ROAM-151
      IF (LM.GE.2*MT) GO TO 8                                           ROAM-152
      L1=1+LM/2                                                         ROAM-153
      FQ=DJCG(IA1,2*IQZ,IA2,2*N,2*MQ,FAC,NFA)*VA(NV+N1)*VA(MV+L1)       ROAM-154
      IF (MQ*N*LM.NE.0) FQ=.7071068D0*FQ                                ROAM-155
      B=B+FQ                                                            ROAM-156
    8 CONTINUE                                                          ROAM-157
      IF (MQ.NE.0) B=.7071068D0*B                                       ROAM-158
      T(3,IT)=DSQRT(DFLOAT(IA1+1))*B*FS                                 ROAM-159
      IF (DABS(T(3,IT)).GT.1.D-6) IT=IT+1                               ROAM-160
    9 CONTINUE                                                          ROAM-161
   10 CONTINUE                                                          ROAM-162
   11 NIV(I1,I2,2)=IT-1                                                 ROAM-163
   12 CONTINUE                                                          ROAM-164
   13 CONTINUE                                                          ROAM-165
      IT=IT-1                                                           ROAM-166
      RETURN                                                            ROAM-167
   14 WRITE (MW,1003)                                                   ROAM-168
      GO TO 17                                                          ROAM-169
   15 WRITE (MW,1004) I1,I2                                             ROAM-170
      GO TO 17                                                          ROAM-171
   16 WRITE (MW,1005)                                                   ROAM-172
   17 WRITE (MW,1006)                                                   ROAM-173
      STOP                                                              ROAM-174
 1000 FORMAT (/' FOR CONSTRAINED ASYMMETRIC ROTATIONAL MODEL THE BETA(I,ROAM-175
     12) ARE INCREASED BY',F10.5//23X,'V',9X,'W',8X,'VS',8X,'WS',7X,'VSOROAM-176
     2',7X,'WSO',6X,'COUL S.O. COUL'/' INITIAL VALUES ',8F10.5)         ROAM-177
 1001 FORMAT (' MODIFIED VALUES',8F10.5)                                ROAM-178
 1002 FORMAT (/' BAND MIXING COEFFICIENTS FOR THE LEVEL',I4/(6D20.7))   ROAM-179
 1003 FORMAT (/' LEVEL ORDER INCORRECT FOR LINK BETWEEN DEFORMATION AND ROAM-180
     1BAND MIXING'/' USE 0+-2+-2+ AND THEN THE OTHER LEVELS.')          ROAM-181
 1004 FORMAT (/' PARITIES OF STATES',I4,'  AND',I4,'  INCORRECT FOR THE ROAM-182
     1ROTATIONAL MODEL.')                                               ROAM-183
 1005 FORMAT (/' NO SPIN-1 STATE IN THIS MODEL.')                       ROAM-184
 1006 FORMAT (//' IN ROAM  ...  STOP  ...')                             ROAM-185
      END                                                               ROAM-186
