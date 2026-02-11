C 12/02/06                                                      ECIS06  PADE-000
      SUBROUTINE PADE(R,P,MM,BRE,BIM,EITER,VLR,VLI,LO)                  PADE-001
C  PADE APPROXIMANT OF TYPE I, CONTINUED FRACTION.                      PADE-002
C INPUT:     R:       PARTIAL TAYLOR SUMS.                              PADE-003
C            MM:      NUMBER OF COMPONENTS OF P.                        PADE-004
C            EITER:   CONVERGENCE CRITERION.                            PADE-005
C            VLR,VLI: VARIABLE   HERE TAKEN AS (1.D0,0.D0).             PADE-006
C            LO(I):   LOGICAL CONTROLS:                                 PADE-007
C               LO(23) =.TRUE. NO USE OF PADE AND SHIFT TO USUAL COUPLEDPADE-008
C                              EQUATIONS WHEN THERE IS NO CONVERGENCE.  PADE-009
C               LO(57) =.TRUE. OUTPUT PHASE-SHIFTS AT EACH ITERATION.   PADE-010
C               LO(105)=.TRUE. CONVERGENCE OBTAINED FOR THIS EQUATION.  PADE-011
C               LO(106)=.TRUE. WHEN THE ITERATION IS NOT THE LAST ONE   PADE-012
C                              ALLOWED.                                 PADE-013
C OUTPUT:    BRE,BIM: RESULT IF CONVERGENCE IS OBTAINED.                PADE-014
C WORKING AREA:                                                         PADE-015
C            P:       AT LEAST TWICE AS LONG AS R.                      PADE-016
C***********************************************************************PADE-017
      IMPLICIT REAL*8 (A-H,O-Z)                                         PADE-018
      DIMENSION R(2,MM),P(2,2,MM)                                       PADE-019
      LOGICAL LO(150)                                                   PADE-020
      COMMON /INOUT/ MR,MW,MS                                           PADE-021
      LO(105)=.TRUE.                                                    PADE-022
      IF (MM.GT.3) GO TO 1                                              PADE-023
      LO(105)=.FALSE.                                                   PADE-024
      RETURN                                                            PADE-025
C TAYLOR COEFFICIENTS.                                                  PADE-026
    1 MT=MM+1                                                           PADE-027
      P(1,1,1)=R(1,1)                                                   PADE-028
      P(1,2,1)=0.D0                                                     PADE-029
      P(2,1,1)=R(2,1)                                                   PADE-030
      P(2,2,1)=0.D0                                                     PADE-031
      DO 2 I=2,MM                                                       PADE-032
      P(1,1,I)=R(1,I)-R(1,I-1)                                          PADE-033
      P(1,2,I)=0.D0                                                     PADE-034
      P(2,1,I)=R(2,I)-R(2,I-1)                                          PADE-035
    2 P(2,2,I)=0.D0                                                     PADE-036
C DECOMPOSITION LOOP.                                                   PADE-037
      NM=MM-1                                                           PADE-038
      DO 6 NA=1,NM                                                      PADE-039
      ML=MM-NA                                                          PADE-040
      ZR=P(1,1,NA)                                                      PADE-041
      ZI=P(2,1,NA)                                                      PADE-042
      IF (DABS(ZR)+DABS(ZI).GE.1.D-8) GO TO 4                           PADE-043
      P(1,2,MT-NA)=0.D0                                                 PADE-044
      P(2,2,MT-NA)=0.D0                                                 PADE-045
      DO 3 J=1,ML                                                       PADE-046
      P(1,1,NA+J)=P(1,1,NA+J)-P(1,1,NA)*P(1,2,J)+P(2,1,NA)*P(2,2,J)     PADE-047
    3 P(2,1,NA+J)=P(2,1,NA+J)-P(2,1,NA)*P(1,2,J)-P(1,1,NA)*P(2,2,J)     PADE-048
      GO TO 6                                                           PADE-049
    4 P(1,2,MT-NA)=1.D0                                                 PADE-050
      DO 5 J=1,ML                                                       PADE-051
      ZX=P(1,1,NA)**2+P(2,1,NA)**2                                      PADE-052
      ZR=(P(1,1,NA+J)*P(1,1,NA)+P(2,1,NA+J)*P(2,1,NA))/ZX               PADE-053
      ZI=(P(2,1,NA+J)*P(1,1,NA)-P(1,1,NA+J)*P(2,1,NA))/ZX               PADE-054
      P(1,1,NA+J)=P(1,2,J)-ZR                                           PADE-055
      P(2,1,NA+J)=P(2,2,J)-ZI                                           PADE-056
      P(1,2,J)=ZR                                                       PADE-057
    5 P(2,2,J)=ZI                                                       PADE-058
    6 CONTINUE                                                          PADE-059
      AR1=1.D30                                                         PADE-060
C TRUNCATED CONTINUED FRACTIONS.                                        PADE-061
      DO 9  N=1,NM                                                      PADE-062
      IP=MM+1-N                                                         PADE-063
      P(1,2,N)=P(1,1,IP)                                                PADE-064
      P(2,2,N)=P(2,1,IP)                                                PADE-065
      DO 8  I=2,IP                                                      PADE-066
      ID=IP+1-I                                                         PADE-067
      IF (P(1,2,MT-ID).GT.0.D0) GO TO 7                                 PADE-068
      AR=VLR*P(1,2,N)-VLI*P(2,2,N)+P(1,1,ID)                            PADE-069
      AI=VLR*P(2,2,N)+VLI*P(1,2,N)+P(2,1,ID)                            PADE-070
      P(1,2,N)=AR                                                       PADE-071
      P(2,2,N)=AI                                                       PADE-072
      GO TO 8                                                           PADE-073
    7 AR=1.D0+VLR*P(1,2,N)-VLI*P(2,2,N)                                 PADE-074
      AI=VLR*P(2,2,N)+VLI*P(1,2,N)                                      PADE-075
      ZX=AR*AR+AI*AI                                                    PADE-076
      P(1,2,N)=(P(1,1,ID)*AR+P(2,1,ID)*AI)/ZX                           PADE-077
      P(2,2,N)=(P(2,1,ID)*AR-P(1,1,ID)*AI)/ZX                           PADE-078
    8 CONTINUE                                                          PADE-079
      IF (N.EQ.1) GO TO 9                                               PADE-080
C SEARCH FOR SMALLEST DIFFERENCE.                                       PADE-081
      ZR=P(1,2,N)-P(1,2,N-1)                                            PADE-082
      ZI=P(2,2,N)-P(2,2,N-1)                                            PADE-083
      AI1=DMAX1(DABS(ZR),DABS(ZI))                                      PADE-084
      IF (AI1.GE.AR1) GO TO 9                                           PADE-085
      NN=N                                                              PADE-086
      AR1=AI1                                                           PADE-087
    9 CONTINUE                                                          PADE-088
      LO(105)=AR1.LE.EITER                                              PADE-089
      IF (LO(57)) WRITE (MW,1000) LO(105),MM,NN,P(1,2,NN),P(2,2,NN),P(1,PADE-090
     12,NN-1),P(2,2,NN-1)                                               PADE-091
      IF ((.NOT.LO(105)).AND.(LO(106).OR.LO(23))) RETURN                PADE-092
      LO(105)=.TRUE.                                                    PADE-093
      BRE=P(1,2,NN-1)                                                   PADE-094
      BIM=P(2,2,NN-1)                                                   PADE-095
      IF (NN.EQ.2) RETURN                                               PADE-096
      BRE=0.5D0*(BRE+P(1,2,NN))                                         PADE-097
      BIM=0.5D0*(BIM+P(2,2,NN))                                         PADE-098
      RETURN                                                            PADE-099
 1000 FORMAT (' PADE',5X,L3,5X,'ITER =',I3,5X,'N =',I3,5X,4D15.8)       PADE-100
      END                                                               PADE-101
