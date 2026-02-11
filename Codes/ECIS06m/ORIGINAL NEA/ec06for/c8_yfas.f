C 13/11/05                                                      ECIS06  YFAS-000
      SUBROUTINE YFAS(ETA,RAU,FO,FPO,GO,GPO,SIGO)                       YFAS-001
C COMPUTATION OF THE COULOMB FUNCTIONS FOR L=0.                         YFAS-002
C     ASYMPTOTIC EXPANSION. SAME ARGUMENTS AS FCZ0 EXCEPT FOR IEXP.     YFAS-003
C INPUT:     ETA,RAU: COULOMB PARAMETER AND RADIUS.                     YFAS-004
C            SIGO:    COULOMB PHASE SHIFT FOR L=0.                      YFAS-005
C OUTPUT:    FO,FPO:  REGULAR FUNCTION AND DERIVATIVE.                  YFAS-006
C            GO,GPO:  IRREGULAR FUNCTION AND DERIVATIVE.                YFAS-007
C***********************************************************************YFAS-008
      IMPLICIT REAL*8 (A-H,O-Z)                                         YFAS-009
      TRB=0.D0                                                          YFAS-010
      RAU2=RAU+RAU                                                      YFAS-011
      ETAC=ETA*ETA                                                      YFAS-012
      N=0                                                               YFAS-013
      PS=1.D0                                                           YFAS-014
      GS=0.D0                                                           YFAS-015
      PT=0.D0                                                           YFAS-016
      GT=1.D0-ETA/RAU                                                   YFAS-017
      SF=PS                                                             YFAS-018
      SG=GS                                                             YFAS-019
      SPF=PT                                                            YFAS-020
      SPG=GT                                                            YFAS-021
    1 DENOM=DFLOAT(N+1)*RAU2                                            YFAS-022
      AN=DFLOAT(N+N+1)*ETA/DENOM                                        YFAS-023
      BN=(ETAC-DFLOAT(N*(N+1)))/DENOM                                   YFAS-024
      PS1=AN*PS-BN*PT                                                   YFAS-025
      GS1=AN*GS-BN*GT-PS1/RAU                                           YFAS-026
      PT1=AN*PT+BN*PS                                                   YFAS-027
      GT1=AN*GT+BN*GS-PT1/RAU                                           YFAS-028
      SF=SF+PS1                                                         YFAS-029
      SG=SG+GS1                                                         YFAS-030
      SPF=SPF+PT1                                                       YFAS-031
      SPG=SPG+GT1                                                       YFAS-032
      N=N+1                                                             YFAS-033
      IF (DABS(PS1).GT.TRB) TRB=DABS(PS1)                               YFAS-034
      IF (DABS(PS1).LT.1.D-10*TRB.OR.BN.LT.-1.D0) GO TO 2               YFAS-035
      PS=PS1                                                            YFAS-036
      GS=GS1                                                            YFAS-037
      PT=PT1                                                            YFAS-038
      GT=GT1                                                            YFAS-039
      GO TO 1                                                           YFAS-040
    2 TETAO=RAU-ETA*DLOG(RAU2)+SIGO                                     YFAS-041
      TRA=DSIN(TETAO)                                                   YFAS-042
      TRB=DCOS(TETAO)                                                   YFAS-043
      GO=SF*TRB-SPF*TRA                                                 YFAS-044
      GPO=SG*TRB-SPG*TRA                                                YFAS-045
      FO=SPF*TRB+SF*TRA                                                 YFAS-046
      FPO=SPG*TRB+SG*TRA                                                YFAS-047
      RETURN                                                            YFAS-048
      END                                                               YFAS-049
