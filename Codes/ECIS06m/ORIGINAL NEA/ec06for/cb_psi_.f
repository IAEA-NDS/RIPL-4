C 12/11/05                                                      ECIS06  PSI_-000
      FUNCTION PSI(ETA)                                                 PSI_-001
C  REAL PART OF PSI(1-I*ETA)-PSI(1)-PSI(I)                              PSI_-002
C  WITH PSI(Z) = LOGARITHMIC DERIVATIVE OF THE GAMMA FUNCTION           PSI_-003
C***********************************************************************PSI_-004
      IMPLICIT REAL*8 (A-H,O-Z)                                         PSI_-005
      DIMENSION B(13)                                                   PSI_-006
      DATA C/.5772156649015329D0/,C1/8.333333333333333D-2/,C2/-8.3333333PSI_-007
     133333333D-3/,C3/3.968253968253968D-3/,C4/-4.166666666666667D-3/,C5PSI_-008
     2/7.575757575757576D-3/,C6/-2.109279609279609D-2/,B/7.4507117898354PSI_-009
     328D-9,2.980350351465228D-8,1.192199259653110D-7,4.769329867878064DPSI_-010
     4-7,1.908212716553938D-6,7.637197637899762D-6,3.058823630702049D-5,PSI_-011
     51.227133475784893D-4,4.941886041194665D-4,2.008392826082234D-3,8.3PSI_-012
     649277381922945D-3,3.692775514337036D-2,2.020569031595948D-1/      PSI_-013
      IF (DABS(ETA).GT.1.D-8) GO TO 1                                   PSI_-014
      PSI=C-1.D0                                                        PSI_-015
      GO TO 5                                                           PSI_-016
    1 E=ETA*ETA                                                         PSI_-017
      IF (E.GT..25D0) GO TO 3                                           PSI_-018
      X=B(1)                                                            PSI_-019
      DO 2 I=2,13                                                       PSI_-020
    2 X=B(I)-E*X                                                        PSI_-021
      PSI=C-1.D0/(1.D0+E)+E*X                                           PSI_-022
      GO TO 5                                                           PSI_-023
    3 L=1                                                               PSI_-024
      IF (E.LT.64.D0) L=2+IDINT(DSQRT(64.D0-E))                         PSI_-025
      X=DSQRT(DFLOAT(L*L)+E)                                            PSI_-026
      Y=DATAN(ETA/DFLOAT(L))                                            PSI_-027
      E=1.D0/(DFLOAT(L*L)+E)                                            PSI_-028
      PSI=DLOG(X)-.5D0*DFLOAT(L)*E-E*(C1*DCOS(2.D0*Y)+E*(C2*DCOS(4.D0*Y)PSI_-029
     1+E*(C3*DCOS(6.D0*Y)+E*(C4*DCOS(8.D0*Y)+E*(C5*DCOS(10.D0*Y)+E*(C6*DPSI_-030
     2COS(12.D0*Y)+E*C1*DCOS(14.D0*Y)))))))+C+C-1.D0                    PSI_-031
      IF (L.EQ.1) GO TO 5                                               PSI_-032
      J=L-1                                                             PSI_-033
      E=ETA*ETA                                                         PSI_-034
      DO 4 I=1,J                                                        PSI_-035
    4 PSI=PSI-1.D0/(DFLOAT(I)+E/DFLOAT(I))                              PSI_-036
    5 RETURN                                                            PSI_-037
      END                                                               PSI_-038
