C 11/01/07                                                      ECIS06  DISP-000
      SUBROUTINE DISP(IPI,WV,IPP,PIP,NVAL,VAL,NIE,LO)                   DISP-001
C COMPUTATION OF THE DEPTHS OF REAL SURFACE AND VOLUME CONTRIBUTIONS    DISP-002
C OBTAINED BY DISPERSION RELATIONS FROM THE IMAGINARY VOLUME AND        DISP-003
C SURFACE POTENTIALS. COMPUTATION OF THE CORRECTIONS TO THE IMAGINARY   DISP-004
C POTENTIALS DUE TO DIFFERENCES OF ENERGY. VARIATION OF HARTREE-FOCK ANDDISP-005
C SPIN-ORBIT POTENTIALS.                                                DISP-006
C INPUT:     IPI(5,I):POTENTIAL FOR THE LEVEL I.                        DISP-007
C            WV(3,I): CENTRE OF MASS ENERGY OF THE LEVEL I.             DISP-008
C            IPP,PIP: EQUIVALENT BY CALL, PARAMETERS OF THE DISPERSION  DISP-009
C                     RELATIONS (SEE INPUT DESCRIPTION).                DISP-010
C            NVAL,VAL:OPTICAL POTENTIALS.                               DISP-011
C            NIE:     TOTAL NUMBER OF LEVELS, NEGATIVE VALUE FOR        DISP-012
C                     COMPOUND CONTINUA.                                DISP-013
C            LO(I):   LOGICAL CONTROLS:                                 DISP-014
C               LO(7)  =.TRUE. MATRIX ELEMENT AND FORM FACTORS READ.    DISP-015
C               LO(10) =.TRUE. DISPERSION RELATIONS FOR POTENTIALS.     DISP-016
C               LO(109)=.TRUE. FOR DIRAC POTENTIALS.                    DISP-017
C               LO(116)=.TRUE. NO OUTPUT.                               DISP-018
C OUTPUT:    WV(J,I): CORRECTION TO POTENTIALS:                         DISP-019
C                     TO VOLUME/SCALAR IMAGINARY POTENTIAL FOR J=14,    DISP-020
C                     TO VOLUME/SCALAR REAL POTENTIAL FOR J=15,         DISP-021
C                     TO SURFACE/VECTOR IMAGINARY POTENTIAL FOR J=16,   DISP-022
C                     TO SURFACE/VECTOR REAL POTENTIAL FOR J=17,        DISP-023
C                     TO SP-O/TENSOR IMAGINARY POTENTIAL FOR J=18,      DISP-024
C                     TO SP-O/TENSOR REAL POTENTIAL FOR J=19.           DISP-025
C***********************************************************************DISP-026
      IMPLICIT REAL*8 (A-H,O-Z)                                         DISP-027
      LOGICAL LO(150)                                                   DISP-028
      DIMENSION IPI(11,*),WV(22,*),IPP(34,*),PIP(17,*),NVAL(2,*),VAL(42,DISP-029
     1*),WW(9),AV(3,2,3),VV(3),EN(4),WX(9,4)                            DISP-030
      COMMON /INOUT/ MR,MW,MS                                           DISP-031
C INPUT OF DISPERSION COEFFICIENTS.                                     DISP-032
      IF (LO(10).AND.(IPP(2,1).EQ.-1)) CALL LDIS(WV,IPP,PIP,NVAL,VAL,LO)DISP-033
      NCOLT=IABS(NIE)                                                   DISP-034
      IPZ=0                                                             DISP-035
      IF (.NOT.LO(116).AND.LO(10)) WRITE (MW,1000)                      DISP-036
C LOOP ON LEVELS.                                                       DISP-037
      DO 24 IV=1,NCOLT                                                  DISP-038
      DO 1 I=1,9                                                        DISP-039
    1 WW(I)=0.D0                                                        DISP-040
      IPI(11,IV)=MAX0(1,IPI(11,IV))                                     DISP-041
      IP=IPI(5,IV)                                                      DISP-042
      IF (.NOT.LO(10)) GO TO 22                                         DISP-043
      IZ=3                                                              DISP-044
      IF (IPP(1,IP).GE.0) IZ=13                                         DISP-045
      IF (IPP(6,IP).EQ.0) GO TO 7                                       DISP-046
      IF (NCOLT.EQ.NIE) GO TO 6                                         DISP-047
C INPUT OF PARAMETERS FOR A CONTINUUM.                                  DISP-048
      IF (IP.EQ.IPZ) GO TO 2                                            DISP-049
      READ (MR,1001) EN                                                 DISP-050
      READ (MR,1002) WX                                                 DISP-051
      IF (.NOT.LO(116)) WRITE (MW,1003) (EN(I),(WX(J,I),J=1,9),I=1,4)   DISP-052
C INTERPOLATION OF PARAMETERS FOR A CONTINUUM.                          DISP-053
    2 Y=WV(IZ,IV)                                                       DISP-054
      DO 5 L=1,4                                                        DISP-055
      X=0.D0                                                            DISP-056
      DO 3 K=1,4                                                        DISP-057
      IF (K.EQ.L) GO TO 3                                               DISP-058
      X=X*(Y-EN(K))/(EN(L)-EN(K))                                       DISP-059
    3 CONTINUE                                                          DISP-060
      DO 4 I=1,9                                                        DISP-061
    4 WW(I)=WW(I)+X*WX(I,L)                                             DISP-062
    5 CONTINUE                                                          DISP-063
      GO TO 22                                                          DISP-064
C INPUT OF PARAMETERS FOR A STATE.                                      DISP-065
    6 READ (MR,1002) WW                                                 DISP-066
      GO TO 22                                                          DISP-067
C SEARCH ON DEPTHS OF REAL VOLUME, IMAGINARY VOLUME AND IMAGINARY       DISP-068
C SPIN-ORBIT POTENTIALS                                                 DISP-069
    7 IF (LO(7)) GO TO 8                                                DISP-070
      VV(1)=VAL(1,IP)                                                   DISP-071
      VV(2)=VAL(5,IP)                                                   DISP-072
      VV(3)=VAL(21,IP)                                                  DISP-073
      GO TO 9                                                           DISP-074
    8 VV(1)=VAL(NVAL(1,8*IP-5),1)                                       DISP-075
      VV(2)=VAL(NVAL(1,8*IP-4),1)                                       DISP-076
      VV(3)=VAL(NVAL(1,8*IP),1)                                         DISP-077
    9 EF=PIP(5,IP)                                                      DISP-078
      EP=PIP(6,IP)                                                      DISP-079
      EO=EP-EF                                                          DISP-080
      YY=PIP(4,IP)-EF                                                   DISP-081
      Y=WV(IZ,IV)-EF                                                    DISP-082
C SPIN-ORBIT OR TENSOR POTENTIAL.                                       DISP-083
      IF (PIP(9,IP).NE.0.D0) WW(7)=-PIP(9,IP)*(Y-YY)/VV(3)              DISP-084
      IF (PIP(8,IP).NE.0.D0) WW(9)=DEXP(-PIP(8,IP)*(Y-YY))-1.D0         DISP-085
C WORKING ARRAY AV(I,J,K) WITH I=1,2,3 FOR VOLUME/SCALAR/VECTOR         DISP-086
C FORM-FACTOR, J=1,2 FOR TWO TERMS, K=1 FOR REFERENCE POTENTIAL,        DISP-087
C K=2 FOR ACTUAL POTENTIAL, K=3 FOR INTERMEDIATE RESULTS.               DISP-088
      DO 21 I=1,3                                                       DISP-089
      NM=IPP(I+2,IP)                                                    DISP-090
      IF (NM.EQ.0) GO TO 19                                             DISP-091
      N=IABS(NM)                                                        DISP-092
      DO 14 K=1,2                                                       DISP-093
      B=PIP(3*I+K+6,IP)                                                 DISP-094
      DO 10 L=1,2                                                       DISP-095
   10 AV(I,K,L)=0.D0                                                    DISP-096
      IF (DABS(YY).GT.EO) AV(I,K,1)=(DABS(YY)-EO)**N/((DABS(YY)-EO)**N+BDISP-097
     1**N)                                                              DISP-098
      IF (DABS(Y).GT.EO) AV(I,K,2)=(DABS(Y)-EO)**N/((DABS(Y)-EO)**N+B**NDISP-099
     1)                                                                 DISP-100
      IF (LO(109).OR.(I.NE.2)) GO TO 11                                 DISP-101
      IF ((N.GT.0).AND.((PIP(14,IP).NE.0.D0).OR.(PIP(15,IP).NE.0.D0))) GDISP-102
     1O TO 20                                                           DISP-103
   11 AV(I,K,3)=DVWI(Y,EO,N,B,0.D0,0.D0)                                DISP-104
      IF ((I.EQ.3).OR.(N.EQ.NM)) GO TO 15                               DISP-105
      IF (K.EQ.1) GO TO 14                                              DISP-106
C SUM OR DIFFERENCE OF TWO VOLUME TERMS.                                DISP-107
      DO 13 J=1,3                                                       DISP-108
   13 AV(I,1,J)=PIP(3*I+9,IP)*AV(I,1,J)+(1.D0-PIP(3*I+9,IP))*AV(I,2,J)  DISP-109
   14 CONTINUE                                                          DISP-110
   15 IF ((N.LT.0).OR.(I.EQ.2)) GO TO 18                                DISP-111
C VOLUME CORRECTION FOR LARGE ENERGY CONTRIBUTIONS.                     DISP-112
      EA=PIP(7,IP)                                                      DISP-113
      IF (EA.EQ.0.D0) GO TO 18                                          DISP-114
      CN=PIP(12,IP)                                                     DISP-115
      IF (IPP(2,IP).EQ.0) GO TO 16                                      DISP-116
C LARGE NEGATIVE ENERGIES.                                              DISP-117
      N2=IABS(IPP(2,IP))                                                DISP-118
      IF (YY.LT.-EA) AV(1,1,1)=AV(1,1,1)*(1.D0-(YY+EA)**N2/((YY+EA)**N2+DISP-119
     1EA**N2)*DEXP(CN*(DSQRT(EA+EF)-DSQRT(-YY+EF))))                    DISP-120
      IF (Y.LT.-EA) AV(1,1,2)=AV(1,1,2)*(1.D0-(Y+EA)**N2/((Y+EA)**N2+EA*DISP-121
     1*N2)*DEXP(CN*(DSQRT(EA+EF)-DSQRT(-Y+EF))))                        DISP-122
      AV(1,1,3)=AV(1,1,3)+DLNE(IPP(2,IP),N,EA,EF,EP,B,CN,WV(IZ,IV))     DISP-123
   16 IF (PIP(11,IP).EQ.0.D0) GO TO 18                                  DISP-124
C LARGE POSITIVE ENERGY TERM.                                           DISP-125
      EL=EF+EA                                                          DISP-126
      ELL=DSQRT(EL)                                                     DISP-127
      AV(1,2,1)=0.D0                                                    DISP-128
      AV(1,2,2)=0.D0                                                    DISP-129
      IF (YY.GT.EA) AV(1,2,1)=(DSQRT(YY+EF)+ELL*(EL/(YY+EF)-3.D0)/2.D0)*DISP-130
     1DEXP(-CN*(DSQRT(YY+EF)-ELL))                                      DISP-131
      IF (Y.GT.EA) AV(1,2,2)=(DSQRT(Y+EF)+ELL*(EL/(Y+EF)-3.D0)/2.D0)*DEXDISP-132
     1P(-CN*(DSQRT(Y+EF)-ELL))                                          DISP-133
      AV(1,2,3)=DLPE(EL,Y+EF,EF,CN)                                     DISP-134
      XX=(VV(2)-AV(1,2,1)*PIP(11,IP))/AV(1,1,1)                         DISP-135
      DO 17 K=1,3                                                       DISP-136
   17 AV(1,1,K)=XX*AV(1,1,K)+PIP(11,IP)*AV(1,2,K)                       DISP-137
   18 WW(3*I-2)=AV(I,1,2)/AV(I,1,1)-1.D0                                DISP-138
      IF ((I.NE.3).OR.(NM.GT.0)) WW(3*I-1)=AV(I,1,3)/AV(I,1,1)          DISP-139
   19 IF ((I.EQ.1).AND.(PIP(17,IP).NE.0.D0)) WW(3)=DEXP(-PIP(17,IP)*(Y-YDISP-140
     1Y))-1.D0                                                          DISP-141
      GO TO 21                                                          DISP-142
C SURFACE FORM-FACTOR.                                                  DISP-143
   20 CS=PIP(14,IP)                                                     DISP-144
      CN=PIP(15,IP)                                                     DISP-145
      AV(2,1,1)=AV(2,1,1)*DEXP(-CN*YY-CS*DABS(YY))                      DISP-146
      AV(2,1,2)=AV(2,1,2)*DEXP(-CN*Y-CS*DABS(Y))                        DISP-147
      WW(4)=AV(2,1,2)/AV(2,1,1)-1.D0                                    DISP-148
      WW(5)=DVWI(Y,EO,N,B,CS,CN)/AV(2,1,1)                              DISP-149
   21 CONTINUE                                                          DISP-150
   22 DO 23 I=1,9                                                       DISP-151
   23 WV(I+13,IV)=WW(I)                                                 DISP-152
      IPZ=IP                                                            DISP-153
      IF (.NOT.LO(116).AND.LO(10)) WRITE (MW,1004) IP,IV,PIP(4,IP),WV(IZDISP-154
     1,IV),WW                                                           DISP-155
   24 CONTINUE                                                          DISP-156
      RETURN                                                            DISP-157
 1000 FORMAT (' POTENTIAL REFER. ENERGY  ACTUAL ENERGY',10X,'VOLUME/SCALDISP-158
     1AR CORRECTION',21X,'SURFACE/VECTOR CORRECTION'/5X,'LEVEL',37X,'SPIDISP-159
     2N-ORBIT/TENSOR CORRECTION'/38X,2('   IMAGINARY      DISPERSION    DISP-160
     3 REAL',8X))                                                       DISP-161
 1001 FORMAT (4F10.5)                                                   DISP-162
 1002 FORMAT (6F10.5/3F10.5)                                            DISP-163
 1003 FORMAT (' ENERGY FOR INTERPOLATION',20X,'COEFFICIENTS'/1P,7D15.6/6DISP-164
     10X,3D15.6)                                                        DISP-165
 1004 FORMAT (1X,I4,I5,1P,2D14.6,6D15.8/38X,3D15.8)                     DISP-166
      END                                                               DISP-167
