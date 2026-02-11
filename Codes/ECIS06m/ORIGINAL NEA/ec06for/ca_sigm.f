C 02/11/05                                                      ECIS06  SIGM-000
      FUNCTION SIGM(ETA)                                                SIGM-001
C  COULOMB PHASE SHIFT SIGMA0.                                          SIGM-002
C***********************************************************************SIGM-003
      IMPLICIT REAL*8 (A-H,O-Z)                                         SIGM-004
      DIMENSION C(2,13)                                                 SIGM-005
      DATA C1,C2,C3,C4,C5,C6,C7,C /8.333333333333333D-2,-2.7777777777777SIGM-006
     178D-3,7.936507936507937D-4,-5.952380952380952D-4,8.417508417508417SIGM-007
     2D-4,-1.917526917526918D-3,6.41025641025641D-3,1.D-16,1.4D-15,-5.4DSIGM-008
     3-15,-2.07D-14,5.1D-13,-3.6968D-12,7.7823D-12,1.043427D-10,-1.18127SIGM-009
     446D-9,5.0020075D-9,6.116095D-9,-2.056338417D-7,1.133027232D-6,-1.2SIGM-010
     5504934821D-6,-2.01348547807D-5,1.280502823882D-4,-2.152416741149D-SIGM-011
     64,-1.1651675918591D-3,7.218943246663D-3,-9.621971527877D-3,-4.2197SIGM-012
     77345555443D-2,.1665386113822915D0,-4.20026350340952D-2,-.655878071SIGM-013
     85202538D0,.5772156649015329D0,1.D0/                               SIGM-014
      IF (DABS(ETA).GT.1.D-16) GO TO 1                                  SIGM-015
      SIGM=-C(1,13)*ETA                                                 SIGM-016
      GO TO 5                                                           SIGM-017
    1 E=ETA*ETA                                                         SIGM-018
      IF (E.GT.1.D0) GO TO 3                                            SIGM-019
      X=C(1,1)                                                          SIGM-020
      Y=C(2,1)                                                          SIGM-021
      DO 2 I=2,13                                                       SIGM-022
      X=C(1,I)-E*X                                                      SIGM-023
    2 Y=C(2,I)-E*Y                                                      SIGM-024
      SIGM=-DATAN2(ETA*X,Y)                                             SIGM-025
      GO TO 5                                                           SIGM-026
    3 L=1                                                               SIGM-027
      IF (E.LT.64.D0) L=2.D0+DSQRT(64.D0-E)                             SIGM-028
      Z=DFLOAT(L)                                                       SIGM-029
      X=DSQRT(Z*Z+E)                                                    SIGM-030
      Y=DATAN2(ETA,Z)                                                   SIGM-031
      E=1.D0/(Z*Z+E)                                                    SIGM-032
      SIGM=ETA*(DLOG(X)-1.D0)+(Z-.5D0)*Y-(C1*DSIN(Y)+E*(C2*DSIN(3.D0*Y)+SIGM-033
     1E*(C3*DSIN(5.D0*Y)+E*(C4*DSIN(7.D0*Y)+E*(C5*DSIN(9.D0*Y)+E*(C6*DSISIGM-034
     2N(11.D0*Y)+E*C7*DSIN(13.D0*Y)))))))/X                             SIGM-035
      IF (L.EQ.1) GO TO 5                                               SIGM-036
      J=L-1                                                             SIGM-037
      DO 4 I=1,J                                                        SIGM-038
      Z=Z-1.D0                                                          SIGM-039
    4 SIGM=SIGM-DATAN2(ETA,Z)                                           SIGM-040
    5 RETURN                                                            SIGM-041
      END                                                               SIGM-042
