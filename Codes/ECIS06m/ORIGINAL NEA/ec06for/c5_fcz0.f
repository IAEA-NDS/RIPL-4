C 12/01/07                                                      ECIS06  FCZ0-000
      SUBROUTINE FCZ0(ETA,RO,F0,FP0,G0,GP0,IEXP,SIGMA)                  FCZ0-001
C COMPUTATION OF THE COULOMB FUNCTIONS FOR L=0.                         FCZ0-002
C INPUT:     ETA,RO:  COULOMB PARAMETER AND RADIUS.                     FCZ0-003
C OUTPUT:    F,FP:    REGULAR FUNCTION AND DERIVATIVE.                  FCZ0-004
C            G,GP:    IRREGULAR FUNCTION AND DERIVATIVE.                FCZ0-005
C            IEXP:    POWER OF 10.                                      FCZ0-006
C            SIGMA:   COULOMB PHASE SHIFT FOR L=0.                      FCZ0-007
C***********************************************************************FCZ0-008
      IMPLICIT REAL*8 (A-H,O-Z)                                         FCZ0-009
      DATA PI /3.1415926535897932D0/                                    FCZ0-010
      SIGMA=SIGM(ETA)                                                   FCZ0-011
      IEXP=0                                                            FCZ0-012
      IF (ETA.LE.28.D0.AND.ETA.GE.-8.D0) GO TO 1                        FCZ0-013
      CALL YFRI(ETA,RO,F0,FP0,G0,GP0,IEXP,SIGMA)                        FCZ0-014
      RETURN                                                            FCZ0-015
    1 IF (ETA.NE.0.D0) GO TO 2                                          FCZ0-016
C BESSEL FUNCTIONS.                                                     FCZ0-017
      F0=DSIN(RO)                                                       FCZ0-018
      G0=DCOS(RO)                                                       FCZ0-019
      FP0=G0                                                            FCZ0-020
      GP0=-F0                                                           FCZ0-021
      RETURN                                                            FCZ0-022
    2 BORNE=1.666666666666667D0*DABS(ETA)+7.5D0                         FCZ0-023
      IF (RO.LT.BORNE) GO TO 3                                          FCZ0-024
      CALL YFAS(ETA,RO,F0,FP0,G0,GP0,SIGMA)                             FCZ0-025
      RETURN                                                            FCZ0-026
    3 IF (ETA.GE.10.D0) GO TO 4                                         FCZ0-027
      IF (ETA.LE.0.D0) GO TO 5                                          FCZ0-028
      IF (RO-2.D0) 15 , 5 , 5                                           FCZ0-029
    4 IF (ETA.GT.(5.D0*RO+6.D1)/7.D0) GO TO 15                          FCZ0-030
C RS=-1 FOR NORMALISATION AT THE ORIGIN,RS=1 AT RO=BORNE(THE BOUNDARY). FCZ0-031
    5 IF (ETA.LT.2.5D0) GO TO 6                                         FCZ0-032
      RS=1.D0                                                           FCZ0-033
      CALL YFAS(ETA,BORNE,F0,FP0,G0,GP0,SIGMA)                          FCZ0-034
      GO TO 8                                                           FCZ0-035
    6 RS=-1.D0                                                          FCZ0-036
C          CLENSHAW ORIGIN.                                             FCZ0-037
      IF (ETA) 7 , 8 , 8                                                FCZ0-038
    7 N=IDINT(-0.5D0*ETA+5.D0)                                          FCZ0-039
      GO TO 9                                                           FCZ0-040
    8 N=IDINT(ETA/5.D0+5.D0)                                            FCZ0-041
    9 N=10*(N/2+1)                                                      FCZ0-042
      TM1=1.D0                                                          FCZ0-043
      T=2.D0*RO/BORNE-1.D0                                              FCZ0-044
      X=T+T                                                             FCZ0-045
      DO 10 I=1,N                                                       FCZ0-046
      TP1=X*T-TM1                                                       FCZ0-047
      TM1=T                                                             FCZ0-048
   10 T=TP1                                                             FCZ0-049
      T=TM1                                                             FCZ0-050
      A1=1.D-30                                                         FCZ0-051
      A2=0.D0                                                           FCZ0-052
      B1=0.D0                                                           FCZ0-053
      B2=A1                                                             FCZ0-054
      S=1.D0                                                            FCZ0-055
      SA=0.D0                                                           FCZ0-056
      SB=0.D0                                                           FCZ0-057
      Z1=0.D0                                                           FCZ0-058
      Z1P=0.D0                                                          FCZ0-059
      AP12=0.D0                                                         FCZ0-060
      AP11=0.D0                                                         FCZ0-061
      BP11=0.D0                                                         FCZ0-062
      Z2=0.D0                                                           FCZ0-063
      Z2P=0.D0                                                          FCZ0-064
      AP22=0.D0                                                         FCZ0-065
      AP21=0.D0                                                         FCZ0-066
      BP21=0.D0                                                         FCZ0-067
      A0=8.D0*ETA/BORNE-1.D0                                            FCZ0-068
      BD=4.D0/(BORNE*BORNE)                                             FCZ0-069
      B0=BD*DFLOAT(N+2)                                                 FCZ0-070
      B4=BD*DFLOAT(N-1)                                                 FCZ0-071
      R4=4.D0*DFLOAT(N)                                                 FCZ0-072
C          DOWNWARDS RECURSION.                                         FCZ0-073
   11 AM11=A0*(A1-AP11)+AP12-B0*B1-B4*BP11                              FCZ0-074
      AM21=A0*(A2-AP21)+AP22-B0*B2-B4*BP21                              FCZ0-075
      SA=SA+S*A1                                                        FCZ0-076
      SB=SB+S*A2                                                        FCZ0-077
      Z1=Z1+A1*T                                                        FCZ0-078
      Z1P=Z1P+B1*T                                                      FCZ0-079
      Z2=Z2+A2*T                                                        FCZ0-080
      Z2P=Z2P+B2*T                                                      FCZ0-081
      IF (R4.EQ.0.D0) GO TO 12                                          FCZ0-082
      BM11=R4*A1+BP11                                                   FCZ0-083
      AP12=AP11                                                         FCZ0-084
      AP11=A1                                                           FCZ0-085
      A1=AM11                                                           FCZ0-086
      BP11=B1                                                           FCZ0-087
      B1=BM11                                                           FCZ0-088
      BM21=R4*A2+BP21                                                   FCZ0-089
      AP22=AP21                                                         FCZ0-090
      AP21=A2                                                           FCZ0-091
      A2=AM21                                                           FCZ0-092
      BP21=B2                                                           FCZ0-093
      B2=BM21                                                           FCZ0-094
      B4=B4-BD                                                          FCZ0-095
      B0=B0-BD                                                          FCZ0-096
      R4=R4-4.D0                                                        FCZ0-097
      S=S*RS                                                            FCZ0-098
      TM1=X*T-TP1                                                       FCZ0-099
      TP1=T                                                             FCZ0-100
      T=TM1                                                             FCZ0-101
      GO TO 11                                                          FCZ0-102
   12 A=AP21-AM21                                                       FCZ0-103
      B=AM11-AP11                                                       FCZ0-104
      SA=A*SA+B*SB                                                      FCZ0-105
      A1=A*A1+B*A2                                                      FCZ0-106
      B1=A*B1+B*B2                                                      FCZ0-107
      Z1=A*Z1+B*Z2                                                      FCZ0-108
      Z1P=A*Z1P+B*Z2P                                                   FCZ0-109
      SA=(SA-0.5D0*A1)/T                                                FCZ0-110
      Z1=Z1-0.5D0*A1                                                    FCZ0-111
      Z1P=Z1P-0.5D0*B1                                                  FCZ0-112
      IF (RS.LT.0.D0) GO TO 13                                          FCZ0-113
      S=F0/(BORNE*SA)                                                   FCZ0-114
      GO TO 14                                                          FCZ0-115
   13 PIETA=PI*ETA                                                      FCZ0-116
      S=DEXP(PIETA)                                                     FCZ0-117
      S=DEXP(-PIETA/2.D0)*DSQRT(2.D0*PIETA/(S-1.D0/S))/SA               FCZ0-118
   14 F0=S*RO*Z1                                                        FCZ0-119
      FP0=S*(Z1+RO*Z1P/BORNE)                                           FCZ0-120
      GO TO 18                                                          FCZ0-121
C   SERIES REGULAR AT THE ORIGIN.                                       FCZ0-122
   15 PI=3.141592653589793D0                                            FCZ0-123
      RO2=RO*RO                                                         FCZ0-124
      ETAP=ETA+ETA                                                      FCZ0-125
      PIETA=PI*ETA                                                      FCZ0-126
      B=DEXP(PIETA)                                                     FCZ0-127
      B=DEXP(0.5D0*PIETA)*DSQRT((B-1.D0/B)/(2.D0*PIETA))                FCZ0-128
      U0=0.D0                                                           FCZ0-129
      U1=RO                                                             FCZ0-130
      U=U0+U1                                                           FCZ0-131
      UP=1.D0                                                           FCZ0-132
      XN=2.D0                                                           FCZ0-133
      DO 16 N=2,10000                                                   FCZ0-134
      XN1=XN*(XN-1.D0)                                                  FCZ0-135
      U2=(ETAP*RO*U1-RO2*U0)/XN1                                        FCZ0-136
      U=U+U2                                                            FCZ0-137
      UP=UP+XN*U2/RO                                                    FCZ0-138
      IF (DABS(U1)+DABS(U2).LT.1.D-10*DABS(U)) GO TO 17                 FCZ0-139
      U0=U1                                                             FCZ0-140
      U1=U2                                                             FCZ0-141
   16 XN=XN+1.D0                                                        FCZ0-142
   17 F0=U/B                                                            FCZ0-143
      FP0=UP/B                                                          FCZ0-144
   18 CALL YFIR(ETA,RO,G0,GP0,SIGMA)                                    FCZ0-145
      RETURN                                                            FCZ0-146
      END                                                               FCZ0-147
