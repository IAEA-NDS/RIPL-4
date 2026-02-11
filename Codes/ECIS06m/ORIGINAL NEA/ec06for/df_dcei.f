C 14/04/06                                                      ECIS06  DCEI-000
      SUBROUTINE DCEI(AR,AI,ZR,ZI)                                      DCEI-001
C COMPLEX EXPONENTIAL INTEGRAL FUNCTION MULTIPLIED BY EXPONENTIAL FOR   DCEI-002
C COMPLEX ARGUMENT: Z = EXP(X)*EI(X).                                   DCEI-003
C INPUT:     AR,AI:   REAL AND IMAGINARY PART OF THE ARGUMENT X.        DCEI-004
C OUTPUT:    ZR,ZI:   REAL AND IMAGINARY PART OF THE RESULT Z.          DCEI-005
C***********************************************************************DCEI-006
      IMPLICIT REAL*8 (A-H,O-Z)                                         DCEI-007
      ZR=0.D0                                                           DCEI-008
      ZI=0.D0                                                           DCEI-009
      IF ((AR.EQ.0.D0).AND.(AI.EQ.0.D0)) RETURN                         DCEI-010
      IF (DABS(AR)+18.5D0.GE.32.D0) GO TO 3                             DCEI-011
      IF (DSQRT(1024.D0-(AR+18.5D0)**2)/1.665D0.LT.DABS(AI)) GO TO 3    DCEI-012
C SERIES EXPANSION.                                                     DCEI-013
      ZR=-.57721566490153D0-DLOG(AR**2+AI**2)*0.5D0                     DCEI-014
      ZI=-DATAN2(AI,AR)                                                 DCEI-015
      YR=1.D0                                                           DCEI-016
      YI=0.D0                                                           DCEI-017
      DO 1 M=1,2000                                                     DCEI-018
      AJ=DFLOAT(M)                                                      DCEI-019
      YZ=YR                                                             DCEI-020
      YR=-(YZ*AR-YI*AI)/AJ                                              DCEI-021
      YI=-(YZ*AI+YI*AR)/AJ                                              DCEI-022
      IF (YR**2+YI**2.LT.1.D-30*(AR**2+AI**2)) GO TO 2                  DCEI-023
      ZR=ZR-YR/AJ                                                       DCEI-024
    1 ZI=ZI-YI/AJ                                                       DCEI-025
    2 YR=DEXP(AR)                                                       DCEI-026
      YI=YR*ZI                                                          DCEI-027
      YR=YR*ZR                                                          DCEI-028
      ZR=YR*DCOS(AI)-YI*DSIN(AI)                                        DCEI-029
      ZI=YI*DCOS(AI)+YR*DSIN(AI)                                        DCEI-030
      RETURN                                                            DCEI-031
C CONTINUED FRACTION.                                                   DCEI-032
    3 DO 4 I=1,20                                                       DCEI-033
      AJ=DFLOAT(21-I)                                                   DCEI-034
      ZR=ZR+AR                                                          DCEI-035
      ZI=ZI+AI                                                          DCEI-036
      Z=ZR**2+ZI**2                                                     DCEI-037
      ZR=AJ*ZR/Z                                                        DCEI-038
      ZI=-AJ*ZI/Z                                                       DCEI-039
      ZR=ZR+1.D0                                                        DCEI-040
      Z=ZR**2+ZI**2                                                     DCEI-041
      ZR=AJ*ZR/Z                                                        DCEI-042
    4 ZI=-AJ*ZI/Z                                                       DCEI-043
      ZR=ZR+AR                                                          DCEI-044
      ZI=ZI+AI                                                          DCEI-045
      Z=ZR**2+ZI**2                                                     DCEI-046
      ZR=ZR/Z                                                           DCEI-047
      ZI=-ZI/Z                                                          DCEI-048
      RETURN                                                            DCEI-049
      END                                                               DCEI-050
