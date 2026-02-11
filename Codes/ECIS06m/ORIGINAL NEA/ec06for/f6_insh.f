C 28/02/07                                                      ECIS06  INSH-000
      SUBROUTINE INSH(P,IC,I1,KT,FAM,X,KAB,ISM,V,NAT,AT,LMD,NVI,LO)     INSH-001
C  E. C. I. S. METHOD: INTEGRATION OF A SINGLE HOMOGENEOUS EQUATION BY  INSH-002
C  THE NUMEROV METHOD       - SCHROEDINGER EQUATION -                   INSH-003
C INPUT:     IC:      CHANNEL NUMBER OF THE EQUATION.                   INSH-004
C            FAM(*,I):MATCHING VALUES FOR I=1 TO 4,                     INSH-005
C                     WAVE NUMBER FOR I=6,                              INSH-006
C                     CONSTANTS OF THIS EQUATION FOR I=7 TO 12.         INSH-007
C            KAB:     MAXIMUM NUMBER OF COUPLED CHANNELS.               INSH-008
C            ISM:     NUMBER OF RADIAL POINTS.                          INSH-009
C            V:       REAL, IMAGINARY POTENTIALS AND COUPLINGS.         INSH-010
C            NAT,AT:  COUPLING COEFFICIENTS.                            INSH-011
C            LMD:     FIRST DIMENSION OF TABLES NAT AND AT.             INSH-012
C            NVI:     TABLE OF ADDRESSES IN NAT,AT.                     INSH-013
C            LO(I):   LOGICAL CONTROLS:                                 INSH-014
C               LO(26) =.TRUE. INTEGRATION STABILISED FOR LONG RANGE    INSH-015
C                              CONSTANT POTENTIAL.                      INSH-016
C               LO(27) =.TRUE. NUMEROV'S METHOD FOR SINGLE EQUATIONS.   INSH-017
C               LO(29) =.TRUE. NO DIAGONAL TERMS IN SECOND MEMBER.      INSH-018
C               LO(30) =.TRUE. PURE DWBA CALCULATION.                   INSH-019
C               LO(101)=.TRUE. THERE IS A REAL SPIN-ORBIT POTENTIAL.    INSH-020
C               LO(102)=.TRUE. THERE IS AN IMAGINARY SPIN-ORBIT         INSH-021
C                              POTENTIAL.                               INSH-022
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. INSH-023
C               LO(121)=.TRUE. OPTICAL MODEL WITHOUT COUPLING.          INSH-024
C               LO(129)=.TRUE. REAL SPIN-ORBIT OR DIRAC EQUATION.       INSH-025
C               LO(130)=.TRUE. IMAGINARY SPIN-ORBIT OR DIRAC EQUATION.  INSH-026
C               LO(133)=.TRUE. STORE SCALAR AND COULOMB POTENTIAL       INSH-027
C                              INDEPENDENTLY.                           INSH-028
C OUTPUT:    P(*,I):  REAL/IMAGINARY REGULAR SOLUTION FOR I=1 AND 2,    INSH-029
C                     REAL/IMAGINARY IRREGULAR SOLUTION FOR I=3 AND 4.  INSH-030
C            FAM(*,I):REAL/IMAGINARY PART OF ZERO'S ORDER SCATTERING    INSH-031
C                     COEFFICIENTS FOR I=9 AND 10.                      INSH-032
C            KT:    THE SOLUTION IS NEGLIGIBLE FOR THE FIRST KT POINTS. INSH-033
C WORKING AREA:                                                         INSH-034
C            X:     FOR THE INTEGRATION.                                INSH-035
C                                                                       INSH-036
C FOR THE COMMON  /CONVE/ SEE CALX.                                     INSH-037
C FOR THE COMMON  /POTE2/ SEE REDM.                                     INSH-038
C                                                                       INSH-039
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     INSH-040
C  BJM:       CONVERGENCE COEFFICIENT OF IMAGINARY POTENTIAL.           INSH-041
C  ACONV:     CONVERGENCE CRITERION FOR POTENTIAL AND FUNCTION.         INSH-042
C   USED:     BJM,ACONV.                                                INSH-043
C                                                                       INSH-044
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     INSH-045
C  ITY(2):    STARTING ADDRESS OF IMAGINARY CENTRAL POTENTIAL.          INSH-046
C  ITY(3):    STARTING ADDRESS OF REAL SPIN-ORBIT POTENTIAL.            INSH-047
C  ITY(4):    STARTING ADDRESS OF IMAGINARY SPIN-ORBIT POTENTIAL.       INSH-048
C  ITY(9):    STARTING ADDRESS OF COULOMB CENTRAL POTENTIAL.            INSH-049
C  ITY(10):   STARTING ADDRESS OF COULOMB SPIN-ORBIT POTENTIAL.         INSH-050
C   USED:     ITY(2),ITY(3),ITY(4),ITY(9),ITY(10).                      INSH-051
C                                                                       INSH-052
C***********************************************************************INSH-053
      IMPLICIT REAL*8 (A-H,O-Z)                                         INSH-054
      LOGICAL LO(150)                                                   INSH-055
      DIMENSION P(ISM,4),FAM(KAB,12),X(2,*),V(ISM,*),NAT(2*LMD,*),AT(LMDINSH-056
     1,*),NVI(KAB,KAB,2)                                                INSH-057
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       INSH-058
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         INSH-059
      BCONV=ACONV                                                       INSH-060
C FOR CLOSED CHANNELS WHEN GREEN'S FUNCTION IS USED.                    INSH-061
      IF (FAM(IC,8).LT.0.D0) BCONV=1.D-15                               INSH-062
C COMPUTATION OF THE REGULAR SOLUTION.                                  INSH-063
      I2=I1+ITY(2)                                                      INSH-064
      DO 1 IS=1,ISM                                                     INSH-065
      X(1,IS+2)=FAM(IC,8)-FAM(IC,10)/DFLOAT(IS*IS)+FAM(IC,7)*V(IS,I1)   INSH-066
    1 X(2,IS+2)=FAM(IC,7)*V(IS,I2)*(1.D0+BJM)                           INSH-067
      IF ((.NOT.LO(101)).OR.(FAM(IC,9).EQ.0.D0))  GO TO 4               INSH-068
      I2=I1+ITY(3)                                                      INSH-069
      DO 2 IS=1,ISM                                                     INSH-070
    2 X(1,IS+2)=X(1,IS+2)+FAM(IC,9)*V(IS,I2)                            INSH-071
      IF (.NOT.LO(102))  GO TO 4                                        INSH-072
      I2=I1+ITY(4)                                                      INSH-073
      DO 3 IS=1,ISM                                                     INSH-074
    3 X(2,IS+2)=X(2,IS+2)+FAM(IC,9)*V(IS,I2)                            INSH-075
    4 IF (.NOT.LO(133)) GO TO 8                                         INSH-076
      IF (FAM(IC,11).EQ.0.D0) GO TO 6                                   INSH-077
      I2=I1+ITY(9)                                                      INSH-078
      DO 5 IS=1,ISM                                                     INSH-079
    5 X(1,IS+2)=X(1,IS+2)+FAM(IC,11)*V(IS,I2)                           INSH-080
    6 IF (FAM(IC,12).EQ.0.D0) GO TO 8                                   INSH-081
      I2=I1+ITY(10)                                                     INSH-082
      DO 7 IS=1,ISM                                                     INSH-083
    7 X(1,IS+2)=X(1,IS+2)+FAM(IC,12)*V(IS,I2)                           INSH-084
    8 IF (.NOT.LO(29)) GO TO 14                                         INSH-085
      K1=NVI(IC,IC,1)                                                   INSH-086
      K2=NVI(IC,IC,2)                                                   INSH-087
      IF (K1.GT.K2) GO TO 14                                            INSH-088
      DO 13 K=K1,K2                                                     INSH-089
      KT=NAT(1,K)                                                       INSH-090
      DO 9 IS=1,ISM                                                     INSH-091
    9 X(1,IS+2)=X(1,IS+2)+AT(2,K)*V(IS,KT)                              INSH-092
      IF (NAT(2,K).EQ.0) GO TO 11                                       INSH-093
      KT=NAT(2,K)+ITY(2)                                                INSH-094
      DO 10 IS=1,ISM                                                    INSH-095
   10 X(2,IS+2)=X(2,IS+2)+AT(2,K)*V(IS,KT)                              INSH-096
   11 IF (.NOT.LO(133)) GO TO 13                                        INSH-097
      IF (AT(3,K).EQ.0.D0) GO TO 13                                     INSH-098
      KT=NAT(1,K)+ITY(9)                                                INSH-099
      DO 12 IS=1,ISM                                                    INSH-100
   12 X(1,IS+2)=X(1,IS+2)+AT(3,K)*V(IS,KT)                              INSH-101
   13 CONTINUE                                                          INSH-102
   14 IF (LO(27)) GO TO 16                                              INSH-103
C MODIFIED NUMEROV METHOD.                                              INSH-104
      DO 15 IS=1,ISM                                                    INSH-105
      A=X(1,IS+2)**2-X(2,IS+2)**2                                       INSH-106
      IF (LO(26)) A=A-X(1,IS+2)**3/30.D0                                INSH-107
      X(2,IS+2)=X(2,IS+2)*(1.D0-X(1,IS+2)/6.D0)                         INSH-108
   15 X(1,IS+2)=X(1,IS+2)-A/12.D0                                       INSH-109
      GO TO 18                                                          INSH-110
C NUMEROV METHOD.                                                       INSH-111
   16 DO 17 IS=1,ISM                                                    INSH-112
      B=(12.D0+X(1,IS+2))**2+X(2,IS+2)**2                               INSH-113
      A=12.D0*(X(1,IS+2)*(12.D0+X(1,IS+2))+X(2,IS+2)**2)/B              INSH-114
      IF (LO(26)) A=A+X(1,IS+2)**3/240.D0                               INSH-115
      X(1,IS+2)=A                                                       INSH-116
   17 X(2,IS+2)=144.D0*X(2,IS+2)/B                                      INSH-117
   18 X(1,1)=0.D0                                                       INSH-118
      X(2,1)=0.D0                                                       INSH-119
      X(1,2)=1.D-15                                                     INSH-120
      X(2,2)=0.D0                                                       INSH-121
      DO 20 IS=1,ISM                                                    INSH-122
      P(IS,3)=X(1,IS+2)                                                 INSH-123
      P(IS,4)=X(2,IS+2)                                                 INSH-124
      HX=X(1,IS+1)*X(1,IS+2)-X(2,IS+1)*X(2,IS+2)                        INSH-125
      HY=X(2,IS+1)*X(1,IS+2)+X(1,IS+1)*X(2,IS+2)                        INSH-126
      X(1,IS+2)=X(1,IS+1)+X(1,IS+1)-X(1,IS)-HX                          INSH-127
      X(2,IS+2)=X(2,IS+1)+X(2,IS+1)-X(2,IS)-HY                          INSH-128
      IF (DABS(X(1,IS+2)).LT.1.D15) GO TO 20                            INSH-129
C  RENORMALISATION OF LARGE FUNCTION VALUES.                            INSH-130
      J=IS+2                                                            INSH-131
      DO 19 I=2,J                                                       INSH-132
      X(1,I)=X(1,I)*1.D-30                                              INSH-133
   19 X(2,I)=X(2,I)*1.D-30                                              INSH-134
   20 CONTINUE                                                          INSH-135
C END OF INTEGRATION.                                                   INSH-136
C MATCHING.                                                             INSH-137
      BRE=X(1,ISM)*FAM(IC,4)-FAM(IC,3)*X(1,ISM+2)                       INSH-138
      BIM=X(2,ISM)*FAM(IC,4)-FAM(IC,3)*X(2,ISM+2)                       INSH-139
      HX=X(1,ISM)*FAM(IC,2)-FAM(IC,1)*X(1,ISM+2)                        INSH-140
      HY=X(2,ISM)*FAM(IC,2)-FAM(IC,1)*X(2,ISM+2)                        INSH-141
      IF (FAM(IC,8).LT.0.D0) GO TO 21                                   INSH-142
      BIM=BIM+HX                                                        INSH-143
      BRE=BRE-HY                                                        INSH-144
   21 A=BRE*BRE+BIM*BIM                                                 INSH-145
      BRE=-BRE/A                                                        INSH-146
      BIM=BIM/A                                                         INSH-147
      FAM(IC,9)=BRE*HX-BIM*HY                                           INSH-148
      FAM(IC,10)=HX*BIM+HY*BRE                                          INSH-149
      IF (LO(121)) RETURN                                               INSH-150
      BRE=BRE/12.D0                                                     INSH-151
      BIM=BIM/12.D0                                                     INSH-152
C NORMALISATION OF THE REGULAR SOLUTION.                                INSH-153
      DO 22 IS=1,ISM                                                    INSH-154
      HX=X(1,IS)+10.D0*X(1,IS+1)+X(1,IS+2)                              INSH-155
      HY=X(2,IS)+10.D0*X(2,IS+1)+X(2,IS+2)                              INSH-156
      P(IS,1)=HX*BRE-HY*BIM                                             INSH-157
   22 P(IS,2)=HX*BIM+HY*BRE                                             INSH-158
C SEARCH OF THE FIRST NON NEGLIGIBLE VALUE.                             INSH-159
      DO 23 KT=1,ISM                                                    INSH-160
      IF ((DABS(P(KT,1))+DABS(P(KT,2))).GT.BCONV) GO TO 24              INSH-161
   23 CONTINUE                                                          INSH-162
   24 IF (LO(30)) RETURN                                                INSH-163
      HX=0.D0                                                           INSH-164
      HY=0.D0                                                           INSH-165
      IF (FAM(IC,8).LT.0.D0) GO TO 25                                   INSH-166
      HX=FAM(IC,1)                                                      INSH-167
      HY=FAM(IC,2)                                                      INSH-168
C COMPUTATION OF THE IRREGULAR SOLUTION - STARTING VALUES.              INSH-169
C  THE TWO LAST POINTS OF IMAGINARY POTENTIAL ARE ZERO.                 INSH-170
   25 X(2,ISM)=HX/FAM(IC,6)                                             INSH-171
      X(2,ISM+2)=HY/FAM(IC,6)                                           INSH-172
      X(1,ISM)=FAM(IC,3)/FAM(IC,6)                                      INSH-173
      X(1,ISM+2)=FAM(IC,4)/FAM(IC,6)                                    INSH-174
      BRE=2.D0-P(ISM,3)                                                 INSH-175
      X(1,ISM+1)=(X(1,ISM)+X(1,ISM+2))/BRE                              INSH-176
      X(2,ISM+1)=(X(2,ISM)+X(2,ISM+2))/BRE                              INSH-177
      J1=ISM-KT                                                         INSH-178
      IF (J1.LE.0) GO TO 27                                             INSH-179
C INTEGRATION.                                                          INSH-180
      DO 26 JS=1,J1                                                     INSH-181
      IS=ISM-JS                                                         INSH-182
      HX=X(1,IS+1)*P(IS,3)-X(2,IS+1)*P(IS,4)                            INSH-183
      HY=X(2,IS+1)*P(IS,3)+X(1,IS+1)*P(IS,4)                            INSH-184
      X(1,IS)=X(1,IS+1)+X(1,IS+1)-X(1,IS+2)-HX                          INSH-185
   26 X(2,IS)=X(2,IS+1)+X(2,IS+1)-X(2,IS+2)-HY                          INSH-186
C COMPUTATION OF IRREGULAR SOLUTION.                                    INSH-187
   27 DO 28 IS=KT,ISM                                                   INSH-188
      P(IS,3)=(X(1,IS)+10.D0*X(1,IS+1)+X(1,IS+2))/12.D0                 INSH-189
   28 P(IS,4)=(X(2,IS)+10.D0*X(2,IS+1)+X(2,IS+2))/12.D0                 INSH-190
      RETURN                                                            INSH-191
      END                                                               INSH-192
