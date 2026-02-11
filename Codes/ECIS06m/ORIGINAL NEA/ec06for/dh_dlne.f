C 11/01/07                                                      ECIS06  DLNE-000
      FUNCTION DLNE(I2,NV,EA,EF,EP,B,CV,EX)                             DLNE-001
C CONTRIBUTION TO THE REAL VOLUME POTENTIAL OF THE LARGE NEGATIVE       DLNE-002
C ENERGIES.                                                             DLNE-003
C INPUT:     I2:      POWER FOR NEGATIVE CORRECTION, NEGATIVE VALUE     DLNE-004
C                     TO TAKE INTO ACCOUNT THE IMAGINARY POTENTIAL      DLNE-005
C                     ONLY AT THE FIRST POINT.                          DLNE-006
C            NV:      POWER FOR VOLUME OR SCALAR POTENTIAL.             DLNE-007
C            EA:      LARGE ENERGY STARTING VALUE ABOVE FERMI ENERGY.   DLNE-008
C            EF:      FERMI ENERGY.                                     DLNE-009
C            EP:      THRESHOLD ENERGY.                                 DLNE-010
C            B:       CONSTANT FOR VOLUME OR SCALAR POTENTIAL.          DLNE-011
C            CV:      EXPONENTIAL DECREASE IN |EX-EF-EA| FOR LARGE      DLNE-012
C                     ENERGY TERMS, NEGATIVE FOR DECREASE IN            DLNE-013
C                     SQRT|EX-EF-EA|.                                   DLNE-014
C            EX:      ENERGY.                                           DLNE-015
C***********************************************************************DLNE-016
      IMPLICIT REAL*8 (A-H,O-Z)                                         DLNE-017
      DATA PI /3.1415926535897932D0/                                    DLNE-018
      N2=IABS(I2)                                                       DLNE-019
      M2=N2/2                                                           DLNE-020
      MV=NV/2                                                           DLNE-021
      EL=EA-EF                                                          DLNE-022
      AG1X=(EL+EX)**N2/((EL+EX)**N2+EA**N2)                             DLNE-023
      AG2X=(2*EF-EP-EX)**NV/((2*EF-EP-EX)**NV+B**NV)                    DLNE-024
      AG2F=(EF-EP)**NV/((EF-EP)**NV+B**NV)                              DLNE-025
      IF (CV.NE.0.D0) GO TO 1                                           DLNE-026
      IF (I2.LT.0) DLNE=.5D0*DLOG(EA)-AG1X*DLOG(EX+EL)                  DLNE-027
      IF (I2.GT.0) DLNE=.5D0*AG2F*DLOG(EA)-AG1X*AG2X*DLOG(EX+EL)        DLNE-028
      GO TO 3                                                           DLNE-029
    1 IF (CV.GT.0.D0) GO TO 2                                           DLNE-030
      CALL DCEI(0.D0,-CV*DSQRT(EX+EL),BRF,BIF)                          DLNE-031
      CALL DCEI(0.D0,-CV*DSQRT(EA),CRF,CIF)                             DLNE-032
      IF (I2.LT.0) DLNE=2*AG1X*BRF-CRF                                  DLNE-033
      IF (I2.GT.0) DLNE=2*AG1X*AG2X*BRF-AG2F*CRF                        DLNE-034
      GO TO 3                                                           DLNE-035
    2 IF (I2.LT.0) DLNE=AG1X*DREI(CV*(EX+EL))-.5D0*DREI(CV*EA)          DLNE-036
      IF (I2.GT.0) DLNE=(AG1X*AG2X*DREI(CV*(EX+EL))-.5D0*AG2F*DREI(CV*EADLNE-037
     1))                                                                DLNE-038
    3 BN=PI/DFLOAT(2*M2)                                                DLNE-039
      KB=1                                                              DLNE-040
      MP=M2                                                             DLNE-041
      MQ=MV                                                             DLNE-042
      AMP=-EL                                                           DLNE-043
      BMP=EA                                                            DLNE-044
      AMQ=2*EF-EP                                                       DLNE-045
      BMQ=B                                                             DLNE-046
    4 DO 9 I=1,MP                                                       DLNE-047
      CN=(2*I-1)*BN                                                     DLNE-048
      AZR=DCOS(CN)*BMP/MP                                               DLNE-049
      AZI=DSIN(CN)*BMP/MP                                               DLNE-050
      BPR=AMP+BMP*DCOS(CN)                                              DLNE-051
      BPI=BMP*DSIN(CN)                                                  DLNE-052
C RESIDUE AZR+I*AZI, POLE BPR+I*BPI                                     DLNE-053
      CTX=(BPR-EX)**2+BPI**2                                            DLNE-054
      CTF=(BPR-EF)**2+BPI**2                                            DLNE-055
      CR=(BPR-EX)/CTX-(BPR-EF)/CTF                                      DLNE-056
      CI=BPI*(1/CTF-1/CTX)                                              DLNE-057
      BR=CR*AZR-CI*AZI                                                  DLNE-058
      BI=CR*AZI+CI*AZR                                                  DLNE-059
C RESIDUE*(1/(P1(I)-EX)-1/(P1(I)-EF))  BR+I*BI                          DLNE-060
      IF (I2.LT.0) GO TO 5                                              DLNE-061
      AT=(((AMQ-BPR)**2+BPI**2)/BMQ**2)**MQ                             DLNE-062
      CN=DATAN2(-BPI,AMQ-BPR)*MQ*2.D0                                   DLNE-063
      ANR=AT*DCOS(CN)                                                   DLNE-064
      ANI=AT*DSIN(CN)                                                   DLNE-065
      CT=AT**2+2*AT*DCOS(CN)+1.D0                                       DLNE-066
      ER=(AT*DCOS(CN)+AT**2)/CT                                         DLNE-067
      EI=AT*DSIN(CN)/CT                                                 DLNE-068
C VALUE OF THE OTHER FUNCTION AT THE POLE ER+I*EI                       DLNE-069
    5 IF (CV.NE.0.D0) GO TO 6                                           DLNE-070
      CR=-DLOG((EL+BPR)**2+BPI**2)*.5D0                                 DLNE-071
      CI=-DATAN2(BPI,EL+BPR)                                            DLNE-072
      GO TO 8                                                           DLNE-073
    6 IF (CV.GT.0.D0) GO TO 7                                           DLNE-074
      AT=CV*((EL+BPR)**2+BPI**2)**.25D0                                 DLNE-075
      CN=DATAN2(BPI,EL+BPR)/2.D0                                        DLNE-076
      CALL DCEI(-AT*DSIN(CN),AT*DCOS(CN),CR,CI)                         DLNE-077
      CALL DCEI(AT*DSIN(CN),-AT*DCOS(CN),DR,DI)                         DLNE-078
      CR=CR+DR                                                          DLNE-079
      CI=CI+DI                                                          DLNE-080
      GO TO 8                                                           DLNE-081
    7 CALL DCEI(CV*(EL+BPR),CV*BPI,CR,CI)                               DLNE-082
    8 DR=(CR*BR-CI*BI)                                                  DLNE-083
      DI=(CR*BI+CI*BR)                                                  DLNE-084
      IF (I2.GT.0) DR=DR*ER-DI*EI                                       DLNE-085
    9 DLNE=DLNE+DR                                                      DLNE-086
      IF (KB.EQ.2) GO TO 11                                             DLNE-087
      IF (I2.LT.0) GO TO 10                                             DLNE-088
      KB=KB+1                                                           DLNE-089
      BN=PI/DFLOAT(NV)                                                  DLNE-090
      MP=MV                                                             DLNE-091
      MQ=M2                                                             DLNE-092
      AMP=2*EF-EP                                                       DLNE-093
      BMP=B                                                             DLNE-094
      AMQ=-EL                                                           DLNE-095
      BMQ=EA                                                            DLNE-096
      GO TO 4                                                           DLNE-097
   10 DLNE=DLNE*(EF-EP+EA)**NV/((EF-EP+EA)**NV+B**NV)                   DLNE-098
   11 DLNE=DLNE/PI                                                      DLNE-099
      RETURN                                                            DLNE-100
      END                                                               DLNE-101
