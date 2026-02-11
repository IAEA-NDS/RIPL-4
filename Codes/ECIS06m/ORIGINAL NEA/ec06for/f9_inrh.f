C 27/06/06                                                      ECIS06  INRH-000
      SUBROUTINE INRH(P,IC,KT,FAM,X,NC,ISM,CC,VV,LO)                    INRH-001
C  E. C. I. S. METHOD: INTEGRATION OF A SINGLE HOMOGENEOUS EQUATION BY  INRH-002
C  THE NUMEROV METHOD   - DIRAC EQUATION -                              INRH-003
C INPUT:     IC:      CHANNEL NUMBER OF THE EQUATION.                   INRH-004
C            FAM(*,I):MATCHING VALUES FOR I=1 TO 4,                     INRH-005
C                     CONSTANTS OF THIS EQUATION FOR I=6 TO 10.         INRH-006
C            KAB:     MAXIMUM NUMBER OF COUPLED CHANNELS.               INRH-007
C            NC:      NUMBER OF COUPLED CHANNELS.                       INRH-008
C            ISM:     NUMBER OF RADIAL POINTS.                          INRH-009
C            CC:      EIGENVALUE OF L.S+1                               INRH-010
C            VV:      POTENTIALS, D(R), ....                            INRH-011
C            LO(I):   LOGICAL CONTROLS:                                 INRH-012
C               LO(26) =.TRUE. INTEGRATION STABILISED FOR LONG RANGE    INRH-013
C                              CONSTANT POTENTIAL.                      INRH-014
C               LO(27) =.TRUE. NUMEROV'S METHOD FOR SINGLE EQUATIONS.   INRH-015
C               LO(121)=.TRUE. OPTICAL MODEL WITHOUT COUPLING.          INRH-016
C OUTPUT:    P(*,I):  REGULAR SOLUTION FOR I=1, 2, 5 AND 6,             INRH-017
C                     IRREGULAR SOLUTION FOR I=3, 4, 7 AND 8,           INRH-018
C                     LARGE COMPONENTS FOR I=1 TO 4,                    INRH-019
C                     SMALL COMPONENT FOR I=5 TO 8,                     INRH-020
C                     REAL/IMAGINARY IRREGULAR SOLUTION FOR I=3 AND 4.  INRH-021
C            FAM(*,I):REAL/IMAGINARY PART OF ZERO'S ORDER SCATTERING    INRH-022
C                     COEFFICIENTS FOR I=9 AND 10.                      INRH-023
C            KT:    THE SOLUTION IS NEGLIGIBLE FOR THE FIRST KT POINTS. INRH-024
C WORKING AREA:                                                         INRH-025
C            X:     FOR THE INTEGRATION.                                INRH-026
C                                                                       INRH-027
C FOR THE COMMON  /CONVE/ SEE CALX.                                     INRH-028
C                                                                       INRH-029
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /CONVE/:                     INRH-030
C  H:         STEP SIZE FOR INTEGRATION.                                INRH-031
C  ACONV:     CONVERGENCE CRITERION FOR POTENTIAL AND FUNCTION.         INRH-032
C   USED:     H,ACONV.                                                  INRH-033
C                                                                       INRH-034
C***********************************************************************INRH-035
      IMPLICIT REAL*8 (A-H,O-Z)                                         INRH-036
      LOGICAL LO(150)                                                   INRH-037
      DIMENSION P(ISM,8),FAM(NC,10),X(2,*),VV(ISM,14)                   INRH-038
      COMMON /CONVE/ H,BJM,EITER,ACONV,CONJ,HCONV                       INRH-039
C COMPUTATION OF THE REGULAR SOLUTION.                                  INRH-040
      DO 1 IS=1,ISM                                                     INRH-041
      X(1,IS+2)=FAM(IC,10)/DFLOAT(IS)**2-FAM(IC,8)-FAM(IC,7)*VV(IS,1)-FAINRH-042
     1M(IC,9)*VV(IS,3)                                                  INRH-043
    1 X(2,IS+2)=-FAM(IC,7)*VV(IS,2)-FAM(IC,9)*VV(IS,4)                  INRH-044
      IF (LO(27)) GO TO 3                                               INRH-045
C MODIFIED NUMEROV METHOD.                                              INRH-046
      DO 2 IS=1,ISM                                                     INRH-047
      A=X(1,IS+2)**2-X(2,IS+2)**2                                       INRH-048
      IF (LO(26)) A=A+X(1,IS+2)**3/30.D0                                INRH-049
      X(2,IS+2)=X(2,IS+2)*(1.D0+X(1,IS+2)/6.D0)                         INRH-050
    2 X(1,IS+2)=X(1,IS+2)+A/12.D0                                       INRH-051
      GO TO 5                                                           INRH-052
C NUMEROV METHOD.                                                       INRH-053
    3 DO 4 IS=1,ISM                                                     INRH-054
      B=(12.D0-X(1,IS+2))**2+X(2,IS+2)**2                               INRH-055
      A=12.D0*(X(1,IS+2)*(12.D0-X(1,IS+2))-X(2,IS+2)**2)/B              INRH-056
      IF (LO(26)) A=A-X(1,IS+2)**3/240.D0                               INRH-057
      X(1,IS+2)=A                                                       INRH-058
    4 X(2,IS+2)=144.D0*X(2,IS+2)/B                                      INRH-059
    5 X(1,1)=0.D0                                                       INRH-060
      X(2,1)=0.D0                                                       INRH-061
      X(1,2)=1.D-15                                                     INRH-062
      X(2,2)=0.D0                                                       INRH-063
      DO 7 IS=1,ISM                                                     INRH-064
      P(IS,3)=X(1,IS+2)                                                 INRH-065
      P(IS,4)=X(2,IS+2)                                                 INRH-066
      HX=X(1,IS+1)*X(1,IS+2)-X(2,IS+1)*X(2,IS+2)                        INRH-067
      HY=X(2,IS+1)*X(1,IS+2)+X(1,IS+1)*X(2,IS+2)                        INRH-068
      X(1,IS+2)=X(1,IS+1)+X(1,IS+1)-X(1,IS)+HX                          INRH-069
      X(2,IS+2)=X(2,IS+1)+X(2,IS+1)-X(2,IS)+HY                          INRH-070
      IF (DABS(X(1,IS+2)).LT.1.D15) GO TO 7                             INRH-071
C RENORMALISATION OF LARGE FUNCTION VALUES.                             INRH-072
      J=IS+2                                                            INRH-073
      DO 6 I=2,J                                                        INRH-074
      X(1,I)=X(1,I)*1.D-30                                              INRH-075
    6 X(2,I)=X(2,I)*1.D-30                                              INRH-076
    7 CONTINUE                                                          INRH-077
C END OF INTEGRATION.                                                   INRH-078
C MATCHING WITH TWO VALUES.                                             INRH-079
      BRE=X(1,ISM)*FAM(IC,4)-FAM(IC,3)*X(1,ISM+2)                       INRH-080
      BIM=X(2,ISM)*FAM(IC,4)-FAM(IC,3)*X(2,ISM+2)                       INRH-081
      HX=X(1,ISM)*FAM(IC,2)-FAM(IC,1)*X(1,ISM+2)                        INRH-082
      HY=X(2,ISM)*FAM(IC,2)-FAM(IC,1)*X(2,ISM+2)                        INRH-083
      BIM=BIM+HX                                                        INRH-084
      BRE=BRE-HY                                                        INRH-085
      BRR=BRE*BRE+BIM*BIM                                               INRH-086
      BRE=-BRE/BRR                                                      INRH-087
      BIM=BIM/BRR                                                       INRH-088
      FAM(IC,9)=BRE*HX-BIM*HY                                           INRH-089
      FAM(IC,10)=HX*BIM+HY*BRE                                          INRH-090
      IF (LO(121)) RETURN                                               INRH-091
      BRE=BRE/12.D0                                                     INRH-092
      BIM=BIM/12.D0                                                     INRH-093
C NORMALISATION OF THE REGULAR SOLUTION.                                INRH-094
C AND SEARCH OF THE FIRST NON NEGLIGIBLE VALUE.                         INRH-095
      DO 8 IS=1,ISM                                                     INRH-096
      BRR=X(1,IS)+10.D0*X(1,IS+1)+X(1,IS+2)                             INRH-097
      HY=X(2,IS)+10.D0*X(2,IS+1)+X(2,IS+2)                              INRH-098
      HX=BRR*BRE-HY*BIM                                                 INRH-099
      HY=HY*BRE+BRR*BIM                                                 INRH-100
      P(IS,1)=HX*VV(IS,9)-HY*VV(IS,10)                                  INRH-101
    8 P(IS,2)=HX*VV(IS,10)+HY*VV(IS,9)                                  INRH-102
      DO 9 KT=1,ISM                                                     INRH-103
      IF ((DABS(P(KT,1))+DABS(P(KT,2))).GT.ACONV) GO TO 10              INRH-104
    9 CONTINUE                                                          INRH-105
   10 HX=FAM(IC,1)                                                      INRH-106
      HY=FAM(IC,2)                                                      INRH-107
C COMPUTATION OF THE IRREGULAR SOLUTION - STARTING VALUES.              INRH-108
C THE LAST POINT OF IMAGINARY POTENTIAL IS NEGLECTED.                   INRH-109
      X(2,ISM)=HX/FAM(IC,6)                                             INRH-110
      X(2,ISM+2)=HY/FAM(IC,6)                                           INRH-111
      X(1,ISM)=FAM(IC,3)/FAM(IC,6)                                      INRH-112
      X(1,ISM+2)=FAM(IC,4)/FAM(IC,6)                                    INRH-113
      BRE=2.D0+P(ISM,3)                                                 INRH-114
      X(1,ISM+1)=(X(1,ISM)+X(1,ISM+2))/BRE                              INRH-115
      X(2,ISM+1)=(X(2,ISM)+X(2,ISM+2))/BRE                              INRH-116
      KR=MIN0(KT,ISM-6)                                                 INRH-117
      I1=ISM-KR                                                         INRH-118
C INTEGRATION.                                                          INRH-119
      DO 11 JS=1,I1                                                     INRH-120
      IS=ISM-JS                                                         INRH-121
      HX=X(1,IS+1)*P(IS,3)-X(2,IS+1)*P(IS,4)                            INRH-122
      HY=X(2,IS+1)*P(IS,3)+X(1,IS+1)*P(IS,4)                            INRH-123
      X(1,IS)=X(1,IS+1)+X(1,IS+1)-X(1,IS+2)+HX                          INRH-124
   11 X(2,IS)=X(2,IS+1)+X(2,IS+1)-X(2,IS+2)+HY                          INRH-125
C COMPUTATION OF IRREGULAR SOLUTION.                                    INRH-126
      DO 12 IS=KR,ISM                                                   INRH-127
      BRE=(X(1,IS)+10.D0*X(1,IS+1)+X(1,IS+2))/12.D0                     INRH-128
      BIM=(X(2,IS)+10.D0*X(2,IS+1)+X(2,IS+2))/12.D0                     INRH-129
      P(IS,3)=BRE*VV(IS,9)-BIM*VV(IS,10)                                INRH-130
   12 P(IS,4)=BRE*VV(IS,10)+BIM*VV(IS,9)                                INRH-131
      ISM3=ISM-3                                                        INRH-132
      HH=H*60.D0                                                        INRH-133
      DO 22 L=1,4                                                       INRH-134
      DO 21 IS=KT,ISM                                                   INRH-135
      IF (IS-KR.GT.2) GO TO 16                                          INRH-136
      IF (IS-KR-1)  13 , 14 , 15                                        INRH-137
   13 P(IS,L+4)=(-147.D0*P(IS,L)+360.D0*P(IS+1,L)-450.D0*P(IS+2,L)+400.DINRH-138
     10*P(IS+3,L)-225.D0*P(IS+4,L)+72.D0*P(IS+5,L)-10.D0*P(IS+6,L))/HH  INRH-139
      GO TO 21                                                          INRH-140
   14 P(IS,L+4)=(-10.D0*P(KR,L)-77.D0*P(IS,L)+150.D0*P(IS+1,L)-100.D0*P(INRH-141
     1IS+2,L)+50.D0*P(IS+3,L)-15.D0*P(IS+4,L)+2.D0*P(IS+5,L))/HH        INRH-142
      GO TO 21                                                          INRH-143
   15 P(IS,L+4)=(2.D0*P(KR,L)-24.D0*P(IS-1,L)-35.D0*P(IS,L)+80.D0*P(IS+1INRH-144
     1,L)-30.D0*P(IS+2,L)+8.D0*P(IS+3,L)-P(IS+4,L))/HH                  INRH-145
      GO TO 21                                                          INRH-146
   16 IF (IS.GT.ISM3) GO TO 17                                          INRH-147
      P(IS,L+4)=(45.D0*(P(IS+1,L)-P(IS-1,L))-9.D0*(P(IS+2,L)-P(IS-2,L))+INRH-148
     1P(IS+3,L)-P(IS-3,L))/HH                                           INRH-149
      GO TO 21                                                          INRH-150
   17 IF (IS-ISM+1) 18 , 19 , 20                                        INRH-151
   18 P(ISM-2,L+4)=(P(ISM-6,L)-8.D0*P(ISM-5,L)+30.D0*P(ISM-4,L)-80.D0*P(INRH-152
     1ISM3,L)+35.D0*P(ISM-2,L)+24.D0*P(ISM-1,L)-2.D0*P(ISM,L))/HH       INRH-153
      GO TO 21                                                          INRH-154
   19 P(ISM-1,L+4)=(-2.D0*P(ISM-6,L)+15.D0*P(ISM-5,L)-50.D0*P(ISM-4,L)+1INRH-155
     100.D0*P(ISM3,L)-150.D0*P(ISM-2,L)+77.D0*P(ISM-1,L)+10.D0*P(ISM,L))INRH-156
     2/HH                                                               INRH-157
      GO TO 21                                                          INRH-158
   20 P(ISM,L+4)=(10.D0*P(ISM-6,L)-72.D0*P(ISM-5,L)+225.D0*P(ISM-4,L)-40INRH-159
     10.D0*P(ISM3,L)+450.D0*P(ISM-2,L)-360.D0*P(ISM-1,L)+147.D0*P(ISM,L)INRH-160
     2)/HH                                                              INRH-161
   21 CONTINUE                                                          INRH-162
   22 CONTINUE                                                          INRH-163
      R=H*DFLOAT(KT-1)                                                  INRH-164
      DO 23 IS=KT,ISM                                                   INRH-165
      R=R+H                                                             INRH-166
      BR=P(IS,1)*(CC/R+VV(IS,13))-P(IS,5)-P(IS,2)*VV(IS,14)             INRH-167
      BI=P(IS,2)*(CC/R+VV(IS,13))-P(IS,6)+P(IS,1)*VV(IS,14)             INRH-168
      CR=P(IS,3)*(CC/R+VV(IS,13))-P(IS,7)-P(IS,4)*VV(IS,14)             INRH-169
      CI=P(IS,4)*(CC/R+VV(IS,13))-P(IS,8)+P(IS,3)*VV(IS,14)             INRH-170
      P(IS,5)=BR*VV(IS,11)-BI*VV(IS,12)                                 INRH-171
      P(IS,6)=BI*VV(IS,11)+BR*VV(IS,12)                                 INRH-172
      P(IS,7)=CR*VV(IS,11)-CI*VV(IS,12)                                 INRH-173
   23 P(IS,8)=CI*VV(IS,11)+CR*VV(IS,12)                                 INRH-174
      IF (KT.EQ.1) RETURN                                               INRH-175
      KR=KT-1                                                           INRH-176
      DO 25 IS=1,KR                                                     INRH-177
      DO 24 L=1,8                                                       INRH-178
   24 P(IS,L)=0.D0                                                      INRH-179
   25 CONTINUE                                                          INRH-180
      RETURN                                                            INRH-181
      END                                                               INRH-182
