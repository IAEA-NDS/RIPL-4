C 12/01/07                                                      ECIS06  YFIR-000
      SUBROUTINE YFIR(ETA,RO,G0,GP0,SIGMA)                              YFIR-001
C COMPUTATION OF THE COULOMB IRREGULAR FUNCTION FOR L=0.                YFIR-002
C INPUT:     ETA,RO:  COULOMB PARAMETER AND RADIUS.                     YFIR-003
C            SIGMA:   COULOMB PHASE SHIFT FOR L=0.                      YFIR-004
C OUTPUT:    GO,GPO:  REGULAR FUNCTION AND DERIVATIVE.                  YFIR-005
C***********************************************************************YFIR-006
      IMPLICIT REAL*8(A-H,O-Z)                                          YFIR-007
      DATA PI /3.1415926535897932D0/                                    YFIR-008
      IF (ETA.LE.0.D0) GO TO 9                                          YFIR-009
      IF ((RO.LE.(54.D0-ETA)/80.D0).OR.(ETA.LE.22.D0.AND.RO.LE.(30.D0-ETYFIR-010
     1A)/20.D0).OR.(ETA.LE.18.D0.AND.RO.LE.0.075D0*(26.D0-ETA)).OR.(ETA.YFIR-011
     2LE.10.D0.AND.RO.LE..15D0*(18.D0-ETA)).OR.(ETA.LE.3.D0.AND.RO.LE.2.YFIR-012
     325D0+7.35D0*(3.D0-ETA))) GO TO 5                                  YFIR-013
C   TAYLOR SERIES STARTING AT RAU0.                                     YFIR-014
    1 RAU0=1.666666666666667D0*DABS(ETA)+7.5D0                          YFIR-015
      CALL YFAS(ETA,RAU0,F0,FP0,G0,GP0,SIGMA)                           YFIR-016
      X=RAU0-RO                                                         YFIR-017
      X2=X*X                                                            YFIR-018
      X3=X*X2                                                           YFIR-019
      UNR=1.D0/RAU0                                                     YFIR-020
      ETR0=1.D0-2.D0*ETA*UNR                                            YFIR-021
      U0=G0                                                             YFIR-022
      U1=-X*GP0                                                         YFIR-023
      U2=-0.5D0*ETR0*X2*U0                                              YFIR-024
      S=U0+U1+U2                                                        YFIR-025
      V1=U1/X                                                           YFIR-026
      V2=2.D0*U2/X                                                      YFIR-027
      T=V1+V2                                                           YFIR-028
      XN=3.D0                                                           YFIR-029
      DO 3 N=3,10000                                                    YFIR-030
      XN1=XN-1.D0                                                       YFIR-031
      XN1=XN*XN1                                                        YFIR-032
      U3=X*U2*UNR*(1.D0-2.D0/XN)-ETR0*U1*X2/XN1+X3*U0*UNR/XN1           YFIR-033
      S=S+U3                                                            YFIR-034
      V3=XN*U3/X                                                        YFIR-035
      T=T+V3                                                            YFIR-036
      IF (DABS(U3).GT.1.D-11*DABS(S)) GO TO 2                           YFIR-037
      IF (DABS(V3).LE.1.D-11*DABS(T)) GO TO 4                           YFIR-038
    2 U0=U1                                                             YFIR-039
      U1=U2                                                             YFIR-040
      U2=U3                                                             YFIR-041
    3 XN=XN+1.D0                                                        YFIR-042
    4 G0=S                                                              YFIR-043
      GP0=-T                                                            YFIR-044
      RETURN                                                            YFIR-045
C   SERIES ORIGIN.                                                      YFIR-046
    5 RO2=RO*RO                                                         YFIR-047
      ETAP=ETA+ETA                                                      YFIR-048
      PIETA=PI*ETA                                                      YFIR-049
      PIETA2=0.5D0*PIETA                                                YFIR-050
      B=DEXP(PIETA)                                                     YFIR-051
      B=DEXP(PIETA2)*DSQRT((B-1.D0/B)/(2.D0*PIETA))                     YFIR-052
      C1=ETAP*(PSI(ETA)+.6931471805599453D0)                            YFIR-053
      U0=0.D0                                                           YFIR-054
      U1=RO                                                             YFIR-055
      V0=1.D0                                                           YFIR-056
      V1=C1*RO                                                          YFIR-057
      U=U0+U1                                                           YFIR-058
      V=V0+V1                                                           YFIR-059
      UP=1.D0                                                           YFIR-060
      VP=C1                                                             YFIR-061
      XN=2.D0                                                           YFIR-062
      DO 7 N=2,10000                                                    YFIR-063
      XN1=XN*(XN-1.D0)                                                  YFIR-064
      U2=(ETAP*RO*U1-RO2*U0)/XN1                                        YFIR-065
      U=U+U2                                                            YFIR-066
      V2=(ETAP*RO*V1-RO2*V0-ETAP*(XN+XN-1.D0)*U2)/XN1                   YFIR-067
      V=V+V2                                                            YFIR-068
      UP=UP+XN*U2/RO                                                    YFIR-069
      VP=VP+XN*V2/RO                                                    YFIR-070
      IF (DABS(U2).GT.1.D-14*DABS(U)) GO TO 6                           YFIR-071
      IF (DABS(V2).LE.1.D-14*DABS(V)) GO TO 8                           YFIR-072
    6 U0=U1                                                             YFIR-073
      U1=U2                                                             YFIR-074
      V0=V1                                                             YFIR-075
      V1=V2                                                             YFIR-076
    7 XN=XN+1.D0                                                        YFIR-077
    8 GP=V+ETAP*U*DLOG(RO)                                              YFIR-078
      G0=B*GP                                                           YFIR-079
      GP0=B*(VP+ETAP*(UP*DLOG(RO)+U/RO))                                YFIR-080
      RETURN                                                            YFIR-081
    9 IF (RO.LE.0.5D0*ETA+9.D0) GO TO 5                                 YFIR-082
      GO TO 1                                                           YFIR-083
      END                                                               YFIR-084
