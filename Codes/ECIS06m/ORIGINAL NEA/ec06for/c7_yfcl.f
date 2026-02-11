C 12/01/07                                                      ECIS06  YFCL-000
      SUBROUTINE YFCL(ETA,RO,U,UP,V,VP,SIGMA,IDIV,NN)                   YFCL-001
C COMPUTATION OF THE COULOMB FUNCTIONS FOR L=0.                         YFCL-002
C CLENSHAW METHOD.  SAME ARGUMENTS AS FCZ0, PLUS NN.                    YFCL-003
C INPUT:     ETA,RO:  COULOMB PARAMETER AND RADIUS.                     YFCL-004
C            SIGMA:   COULOMB PHASE SHIFT FOR L=0.                      YFCL-005
C            NN:      0 FOR ASYMPTOTIC, 1 FOR SERIES AT THE ORIGIN.     YFCL-006
C OUTPUT:    U,UP:    REGULAR FUNCTION AND DERIVATIVE.                  YFCL-007
C            V,VP:    IRREGULAR FUNCTION AND DERIVATIVE.                YFCL-008
C            IDIV:    POWER OF 10.                                      YFCL-009
C***********************************************************************YFCL-010
      IMPLICIT REAL*8 (A-H,O-Z)                                         YFCL-011
      DATA PI /3.1415926535897932D0/                                    YFCL-012
      IDIV=0                                                            YFCL-013
      ETAP=ETA+ETA                                                      YFCL-014
      RO2=RO*RO                                                         YFCL-015
      IF (NN.EQ.1) GO TO 4                                              YFCL-016
C          CLENSHAW ASYMPTOTIC.                                         YFCL-017
      E2=ETA*ETA                                                        YFCL-018
      M=IDINT(40.D0+ETA/4.D0)                                           YFCL-019
      M=2*(M/2)                                                         YFCL-020
      I=1                                                               YFCL-021
      R=DFLOAT(M)                                                       YFCL-022
      D4=4.D0*R                                                         YFCL-023
      TM1=1.D0                                                          YFCL-024
      T=4.D0*ETA/RO-1.D0                                                YFCL-025
      Z=T+T                                                             YFCL-026
      DO 1 J=1,M                                                        YFCL-027
      TP1=Z*T-TM1                                                       YFCL-028
      TM1=T                                                             YFCL-029
    1 T=TP1                                                             YFCL-030
      T=TM1                                                             YFCL-031
      GR=0.D0                                                           YFCL-032
      GPR=0.D0                                                          YFCL-033
      SCR=0.D0                                                          YFCL-034
      EP1R=0.D0                                                         YFCL-035
      ER=0.D0                                                           YFCL-036
      DP1R=0.D0                                                         YFCL-037
      DR=1.D-25                                                         YFCL-038
      C1R=R+1.D0                                                        YFCL-039
      C0R=E2-R*(R+1.D0)                                                 YFCL-040
      GI=0.D0                                                           YFCL-041
      GPI=0.D0                                                          YFCL-042
      SCI=0.D0                                                          YFCL-043
      EP1I=0.D0                                                         YFCL-044
      EI=0.D0                                                           YFCL-045
      DP1I=0.D0                                                         YFCL-046
      DI=0.D0                                                           YFCL-047
      C1I=-3.D0*ETA                                                     YFCL-048
      C0I=-(R+R+1.D0)*ETA                                               YFCL-049
C          BACKWARDS RECURRENCE.                                        YFCL-050
    2 SCZ=C0R*C0R+C0I*C0I                                               YFCL-051
      C2R=C1R*DR-C1I*DI+0.5D0*(DP1R+ER+EP1R)-ETA*DP1I                   YFCL-052
      C2I=C1R*DI+C1I*DR+0.5D0*(DP1I+EI+EP1I)+ETA*DP1R                   YFCL-053
      CR=(C2R*C0R+C2I*C0I)/SCZ                                          YFCL-054
      CI=(C0R*C2I-C0I*C2R)/SCZ                                          YFCL-055
      GR=GR+CR*T                                                        YFCL-056
      GPR=GPR+DR*T                                                      YFCL-057
      SCR=SCR+I*CR                                                      YFCL-058
      GI=GI+CI*T                                                        YFCL-059
      GPI=GPI+DI*T                                                      YFCL-060
      SCI=SCI+I*CI                                                      YFCL-061
      IF (R.EQ.0.D0) GO TO 3                                            YFCL-062
      EM1R=D4*DR+EP1R                                                   YFCL-063
      DM1R=D4*CR+DP1R                                                   YFCL-064
      EP1R=ER                                                           YFCL-065
      ER=EM1R                                                           YFCL-066
      DP1R=DR                                                           YFCL-067
      DR=DM1R                                                           YFCL-068
      EM1I=D4*DI+EP1I                                                   YFCL-069
      DM1I=D4*CI+DP1I                                                   YFCL-070
      EP1I=EI                                                           YFCL-071
      EI=EM1I                                                           YFCL-072
      DP1I=DI                                                           YFCL-073
      DI=DM1I                                                           YFCL-074
      TM1=Z*T-TP1                                                       YFCL-075
      TP1=T                                                             YFCL-076
      T=TM1                                                             YFCL-077
      C0R=C0R+R+R                                                       YFCL-078
      C0I=C0I+ETAP                                                      YFCL-079
      C1R=C1R-1.D0                                                      YFCL-080
      I=-I                                                              YFCL-081
      D4=D4-4.D0                                                        YFCL-082
      R=R-1.D0                                                          YFCL-083
      GO TO 2                                                           YFCL-084
    3 SCR=SCR-0.5D0*CR                                                  YFCL-085
      GR=GR-0.5D0*CR                                                    YFCL-086
      GPR=GPR-0.5D0*DR                                                  YFCL-087
      SCI=SCI-0.5D0*CI                                                  YFCL-088
      GI=GI-0.5D0*CI                                                    YFCL-089
      GPI=GPI-0.5D0*DI                                                  YFCL-090
      Z=SIGMA+RO-ETA*DLOG(RO+RO)                                        YFCL-091
      SCZ=SCR*SCR+SCI*SCI                                               YFCL-092
      CR=DCOS(Z)                                                        YFCL-093
      CI=DSIN(Z)                                                        YFCL-094
      DR=(CR*SCR+CI*SCI)/SCZ                                            YFCL-095
      DI=(CI*SCR-CR*SCI)/SCZ                                            YFCL-096
      SCI=1.D0-ETA/RO                                                   YFCL-097
      SCR=ETAP/RO2                                                      YFCL-098
      CR=-GI*SCI-GPR*SCR                                                YFCL-099
      CI=GR*SCI-GPI*SCR                                                 YFCL-100
      VP=DR*CR-DI*CI                                                    YFCL-101
      UP=DR*CI+DI*CR                                                    YFCL-102
      V=DR*GR-DI*GI                                                     YFCL-103
      U=DR*GI+DI*GR                                                     YFCL-104
      RETURN                                                            YFCL-105
C          SERIES AT THE ORIGIN.                                        YFCL-106
    4 PIETA=PI*ETA                                                      YFCL-107
      IF (DABS(PIETA).GT.36.D0) GO TO 5                                 YFCL-108
      P=DSQRT((DEXP(2.D0*PIETA)-1.D0)/(2.D0*PIETA))                     YFCL-109
      GO TO 7                                                           YFCL-110
    5 IF (PIETA.GT.0.D0) GO TO 6                                        YFCL-111
      P=1.D0/DSQRT(-PIETA-PIETA)                                        YFCL-112
      GO TO 7                                                           YFCL-113
    6 Z=34.588776394910686D0                                            YFCL-114
      IDIV=IDINT(PIETA/Z)                                               YFCL-115
      P=DEXP(PIETA-IDIV*Z)/DSQRT(PIETA+PIETA)                           YFCL-116
      IDIV=15*IDIV                                                      YFCL-117
    7 Z1=ETAP*(PSI(ETA)+.6931471805599453D0)                            YFCL-118
      U0=0.D0                                                           YFCL-119
      U1=RO                                                             YFCL-120
      V0=1.D0                                                           YFCL-121
      V1=Z1*RO                                                          YFCL-122
      U=U0+U1                                                           YFCL-123
      V=V0+V1                                                           YFCL-124
      UP=1.D0                                                           YFCL-125
      VP=Z1                                                             YFCL-126
      XN=2.D0                                                           YFCL-127
      DO 9 N=2,10000                                                    YFCL-128
      XN1=XN*(XN-1.D0)                                                  YFCL-129
      U2=(ETAP*RO*U1-RO2*U0)/XN1                                        YFCL-130
      U=U+U2                                                            YFCL-131
      V2=(ETAP*RO*V1-RO2*V0-ETAP*(XN+XN-1.D0)*U2)/XN1                   YFCL-132
      V=V+V2                                                            YFCL-133
      UP=UP+XN*U2/RO                                                    YFCL-134
      VP=VP+XN*V2/RO                                                    YFCL-135
      IF (DABS(U2).GT.1.D-14*DABS(U)) GO TO 8                           YFCL-136
      IF (DABS(V2).LE.1.D-14*DABS(V)) GO TO 10                          YFCL-137
    8 U0=U1                                                             YFCL-138
      U1=U2                                                             YFCL-139
      V0=V1                                                             YFCL-140
      V1=V2                                                             YFCL-141
    9 XN=XN+1.D0                                                        YFCL-142
   10 PP=V+ETAP*U*DLOG(RO)                                              YFCL-143
      W=U/P                                                             YFCL-144
      WP=UP/P                                                           YFCL-145
      V=P*PP                                                            YFCL-146
      VP=P*(VP+ETAP*(UP*DLOG(RO)+U/RO))                                 YFCL-147
      U=W                                                               YFCL-148
      UP=WP                                                             YFCL-149
      RETURN                                                            YFCL-150
      END                                                               YFCL-151
