C 11/01/07                                                      ECIS06  DLPE-000
      FUNCTION DLPE(EL,EX,EF,CN)                                        DLPE-001
C CONTRIBUTION TO THE REAL VOLUME POTENTIAL OF THE LARGE POSITIVE       DLPE-002
C ENERGIES.                                                             DLPE-003
C INPUT:     EA:      LARGE ENERGY STARTING VALUE ABOVE ZERO ENERGY.    DLPE-004
C            EX:      ENERGY.                                           DLPE-005
C            EF:      FERMI ENERGY.                                     DLPE-006
C            CN:      EXPONENTIAL DECREASE IN |EX-EF-EA| FOR LARGE      DLPE-007
C                     ENERGY TERMS, NEGATIVE FOR DECREASE IN            DLPE-008
C                     SQRT|EX-EF-EA|.                                   DLPE-009
C***********************************************************************DLPE-010
      IMPLICIT REAL*8 (A-H,O-Z)                                         DLPE-011
      DATA PI /3.1415926535897932D0/                                    DLPE-012
      DLPE=0.D0                                                         DLPE-013
      IF (EX.EQ.EF) RETURN                                              DLPE-014
      FF=DSQRT(DABS(EF))                                                DLPE-015
      FL=DSQRT(DABS(EL))                                                DLPE-016
      FX=DSQRT(DABS(EX))                                                DLPE-017
      IF (CN.NE.0.D0) GO TO 5                                           DLPE-018
      DLPE=2.D0*FF*DATAN2(FF,FL)+0.5D0*EL*FL/EF*DLOG(1.D0-EF/EL)-1.5D0*FDLPE-019
     1L*DLOG(DABS(EL-EF))                                               DLPE-020
      IF (EX.GT.0.D0) GO TO 3                                           DLPE-021
      IF (DABS(EX).GT.EL*1.D-5) GO TO 1                                 DLPE-022
      DLPE=DLPE+.5D0*FL*(1.D0+EX/EL/2.D0+(EX/EL)**2/3.D0)               DLPE-023
      GO TO 2                                                           DLPE-024
    1 DLPE=DLPE-0.5D0*EL*FL/EX*DLOG(DABS(1-EX/EL))                      DLPE-025
    2 DLPE=DLPE-2.D0*FX*(DATAN2(FX,FL))+1.5D0*FL*DLOG(DABS(EL-EX))      DLPE-026
      GO TO 16                                                          DLPE-027
    3 IF (DABS(EX).GT.EL*1.D-5) GO TO 4                                 DLPE-028
      DLPE=DLPE+.5D0*FL*(1.D0+EX/EL/2.D0+(EX/EL)**2/3.D0)-FX*(DLOG(DABS(DLPE-029
     1FL-FX)/(FL+FX)))+1.5D0*FL*(DLOG(DABS(EL-EX)))                     DLPE-030
      GO TO 16                                                          DLPE-031
    4 DLPE=DLPE+(FX+1.5D0*FL-0.5D0*EL*FL/EX)*DLOG(FL+FX)+EL*FL/EX*DLOG(FDLPE-032
     1L)                                                                DLPE-033
      IF (DABS(EX-EL).GT.1.D-3) DLPE=DLPE-(FX-1.5D0*FL+0.5D0*EL*FL/EX)*DDLPE-034
     1LOG(DABS(FL-FX))                                                  DLPE-035
      GO TO 16                                                          DLPE-036
    5 IF (CN.GT.0.D0) GO TO 9                                           DLPE-037
      A1=DREI(-CN*FL)                                                   DLPE-038
      A2=DREI(CN*(FX-FL))                                               DLPE-039
      IF (EX.GT.0.D0) GO TO 6                                           DLPE-040
      CALL DCEI(-CN*FL,-CN*FX,A3,A4)                                    DLPE-041
      DLPE=2.D0*FX*A4-3.D0*FL*A3                                        DLPE-042
      IF (EX.GT.-1.D-4) GO TO 7                                         DLPE-043
      DLPE=DLPE+EL*FL/EX*(A3-A1)                                        DLPE-044
      GO TO 8                                                           DLPE-045
    6 A3=DREI(-CN*(FL+FX))                                              DLPE-046
      DLPE=FX*(A2-A3)-1.5D0*FL*(A2+A3)                                  DLPE-047
      IF (EX.LT.1.D-4) GO TO 7                                          DLPE-048
      DLPE=DLPE+.5D0*EL*FL/EX*(A2+A3-2.D0*A1)                           DLPE-049
      GO TO 8                                                           DLPE-050
C USING UP TO THE FOURTH DERIVATIVE OF THE EXPONENTIAL FUNCTION.        DLPE-051
    7 DLPE=DLPE+.5D0*EL*FL*((A1*CN**2+CN/FL+1.D0/EL)*(1.D0+CN**2*EX/12.DDLPE-052
     10)+(CN/FL+3.D0/EL)/EL*EX/6.D0)                                    DLPE-053
    8 CALL DCEI(-CN*FL,-CN*FF,B1,B2)                                    DLPE-054
      DLPE=DLPE-FF*2.D0*B2-EL*FL/EF*(B1-A1)+3.D0*FL*B1                  DLPE-055
      GO TO 16                                                          DLPE-056
    9 IF (EX.LT.0.D0) FX=0.D0                                           DLPE-057
      AV=FX-1.5D0*FL                                                    DLPE-058
      AW=.5D0*EL*FL                                                     DLPE-059
      A2=DREI(CN*(EL-EF))                                               DLPE-060
      A3=DREI(CN*(EL-EX))                                               DLPE-061
      A8=DREI(CN*(EL))                                                  DLPE-062
      A4=(A8-A2)*AW/EF+(A3-A2)*AV                                       DLPE-063
      IF (DABS(EX).GT.1.D-4) GO TO 10                                   DLPE-064
      A4=A4-((A8*CN-1.D0/EL)*(1.D0-.5D0*EX*CN)-0.5D0/EL**2*EX)*AW       DLPE-065
      GO TO 11                                                          DLPE-066
   10 A4=A4-(A8-A3)*AW/EX                                               DLPE-067
   11 A5=A4                                                             DLPE-068
      IMAX=32                                                           DLPE-069
      DO 14 J=1,15                                                      DLPE-070
      A1=0.D0                                                           DLPE-071
      GX=1.D0/EL/DFLOAT(IMAX)                                           DLPE-072
      GY=GX*1.5D0                                                       DLPE-073
      GU=0.D0                                                           DLPE-074
      DO 13 I=1,IMAX                                                    DLPE-075
      IF (I.EQ.IMAX) GY=GX/2.D0                                         DLPE-076
      GU=GU+GX                                                          DLPE-077
      IF (DABS((1.D0-GU*EF)*GU).LT.1.D-6) GO TO 13                      DLPE-078
      IF (EX.GT.0.D0) GO TO 12                                          DLPE-079
      A1=A1+GY*(DSQRT(1.D0/GU))*(1.D0/(1.D0-GU*EX)-1.D0/(1.D0-GU*EF))*DEDLPE-080
     1XP(-CN*(1/GU-EL))/GU                                              DLPE-081
      GO TO 13                                                          DLPE-082
   12 A1=A1+GY/(DSQRT(1.D0/GU)+FX)*DEXP(-CN*(1/GU-EL))/GU**2            DLPE-083
      A1=A1-GY*(DSQRT(1.D0/GU)-FX)/(1.D0-GU*EF)*DEXP(-CN*(1/GU-EL))/GU  DLPE-084
   13 GY=GX                                                             DLPE-085
      DLPE=A1+A4                                                        DLPE-086
      IF (DABS(1.D0-A5/DLPE).LT.1.D-5) GO TO 15                         DLPE-087
      A5=DLPE                                                           DLPE-088
   14 IMAX=2*IMAX                                                       DLPE-089
   15 DLPE=(DLPE-A5)/2+DLPE                                             DLPE-090
   16 DLPE=DLPE/PI                                                      DLPE-091
      RETURN                                                            DLPE-092
      END                                                               DLPE-093
