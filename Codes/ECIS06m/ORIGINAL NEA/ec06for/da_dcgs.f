C 27/06/06                                                      ECIS06  DCGS-000
      FUNCTION DCGS(L,J1,J2,FAC,NFA)                                    DCGS-001
C COMPUTATION OF SINGULAR CLEBSCH-GORDAN COEFFICIENTS.                  DCGS-002
C INPUT:     L,J1,J2: DOUBLE VALUE OF ANGULAR MOMENTA.                  DCGS-003
C            FAC:     TABLE OF LOGARITHM OF FACTORIALS.                 DCGS-004
C            NFA:     LENGTH OF FAC.                                    DCGS-005
C OUTPUT:                                                               DCGS-006
C                                                                       DCGS-007
C                                 /~~~~~~~~~~~~~~~~   ( J2   L   J1  )  DCGS-008
C  CGS(L,J1,J2) = (-)**(J1+1/2)  / (2*J1+1)*(2*J2+1)  (              )  DCGS-009
C                                                     (-1/2  0   1/2 )  DCGS-010
C  WHEN J1 AND J2 ARE HALF-INTEGERS,                                    DCGS-011
C                              ( J2  L  J1 )                            DCGS-012
C  CGS(L,J1,J2) = (-)**(J1-1)  (           )                            DCGS-013
C                              (  0  0  0  )                            DCGS-014
C  WHEN J1 AND J2  ARE INTEGERS.                                        DCGS-015
C                                                                       DCGS-016
C***********************************************************************DCGS-017
      IMPLICIT REAL*8 (A-H,O-Z)                                         DCGS-018
      DIMENSION FAC(*)                                                  DCGS-019
      IF (L.LT.0.OR.J1.LT.0.OR.J2.LT.0) GO TO 1                         DCGS-020
      LL=L+J1+J2+2                                                      DCGS-021
      LT=LL/2                                                           DCGS-022
      IF (2*LT.NE.LL) GO TO 1                                           DCGS-023
      IF (LT.GE.NFA) GO TO 1                                            DCGS-024
      L1=LT-J1                                                          DCGS-025
      L2=LT-J2                                                          DCGS-026
      L3=LT-L                                                           DCGS-027
      IF (L1.LE.0.OR.L2.LE.0.OR.L3.LE.0) GO TO 1                        DCGS-028
      L4=LT+1                                                           DCGS-029
      N1=(L1+1)/2                                                       DCGS-030
      N2=(L2+1)/2                                                       DCGS-031
      N3=(L3+1)/2                                                       DCGS-032
      N4=(L4+1)/2                                                       DCGS-033
      DCGS=DEXP(FAC(N4)-FAC(N1)-FAC(N2)-FAC(N3)-.5D0*(FAC(L4)-FAC(L1)-FADCGS-034
     1C(L2)-FAC(L3)))                                                   DCGS-035
      IF (2*(N1/2).NE.N1) DCGS=-DCGS                                    DCGS-036
      IF (N4-N1-N2-N3+1) 3 , 2 , 1                                      DCGS-037
    1 DCGS=0.D0                                                         DCGS-038
    2 DCGS=2.D0*DCGS                                                    DCGS-039
    3 RETURN                                                            DCGS-040
      END                                                               DCGS-041
