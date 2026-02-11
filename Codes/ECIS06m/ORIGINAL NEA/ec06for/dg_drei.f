C 11/01/07                                                      ECIS06  DREI-000
      FUNCTION DREI(AR)                                                 DREI-001
C REAL EXPONENTIAL INTEGRAL FUNCTION MULTIPLIED BY EXPONENTIAL.         DREI-002
C INPUT:     AR:      ARGUMENT.                                         DREI-003
C***********************************************************************DREI-004
      IMPLICIT REAL*8 (A-H,O-Z)                                         DREI-005
      DREI=0.D0                                                         DREI-006
      IF (AR.EQ.0.D0) RETURN                                            DREI-007
      IF (DABS(AR)+18.5D0.GE.32.D0) GO TO 3                             DREI-008
C SERIES EXPANSION.                                                     DREI-009
      DREI=-.57721566490153D0-DLOG(DABS(AR))                            DREI-010
      YR=1.D0                                                           DREI-011
      DO 1 M=1,2000                                                     DREI-012
      AJ=DFLOAT(M)                                                      DREI-013
      YR=-YR*AR/AJ                                                      DREI-014
      IF (DABS(YR).LT.1.D-15*DABS(AR)) GO TO 2                          DREI-015
    1 DREI=DREI-YR/AJ                                                   DREI-016
    2 DREI=DREI*DEXP(AR)                                                DREI-017
      RETURN                                                            DREI-018
C CONTINUED FRACTION.                                                   DREI-019
    3 DO 4 I=1,20                                                       DREI-020
      AJ=DFLOAT(21-I)                                                   DREI-021
      DREI=DREI+AR                                                      DREI-022
      DREI=AJ/DREI                                                      DREI-023
    4 DREI=AJ/(DREI+1.D0)                                               DREI-024
      DREI=1.D0/(DREI+AR)                                               DREI-025
      RETURN                                                            DREI-026
      END                                                               DREI-027
