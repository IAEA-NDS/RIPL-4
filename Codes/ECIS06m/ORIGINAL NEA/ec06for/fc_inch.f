C 08/03/07                                                      ECIS06  INCH-000
      SUBROUTINE INCH(V,MC,NAT,AT,LMD,NVI,FAM,Y,X,W,ISM,KAB,NC,NCIN,NML,INCH-001
     1JSX,KR,LO)                                                        INCH-002
C  STANDARD INTEGRATION OF THE COUPLED  EQUATIONS BY THE NUMEROV METHOD.INCH-003
C INPUT:     V:       POTENTIALS AND COUPLINGS.                         INCH-004
C            MC(*,J): ADDRESS OF CENTRAL POTENTIALS FOR J=4.            INCH-005
C            NAT,AT:  GEOMETRICAL COEFFICIENTS.                         INCH-006
C            LMD:     FIRST DIMENSION OF TABLES NAT AND AT.             INCH-007
C            NVI:     ADDRESSES IN THE TABLE NAT,AT.                    INCH-008
C            FAM:     MATCHING VALUES AND CONSTANTS OF EQUATIONS.       INCH-009
C            Y:       COULOMB CORRECTIONS.                              INCH-010
C            ISM:     NUMBER OF RADIAL POINTS.                          INCH-011
C            KAB:     MAXIMUM NUMBER OF COUPLED CHANNELS.               INCH-012
C            NC:      NUMBER OF COUPLED CHANNELS.                       INCH-013
C            NCIN:    NUMBER OF SOLUTIONS.                              INCH-014
C            NML:     MAXIMUM NUMBER OF POINTS WHERE THE COUPLING       INCH-015
C                     POTENTIALS CAN BE CALCULATED AT A TIME.           INCH-016
C            JSX:     PERIODICITY OF THE SCHMIDT'S ORTHOGONALISATION    INCH-017
C                     IF LO(42) IS .TRUE..                              INCH-018
C            LO(I):   LOGICAL CONTROLS:                                 INCH-019
C               LO(26) =.TRUE. INTEGRATION STABILISED FOR LONG RANGE    INCH-020
C                              CONSTANT POTENTIAL.                      INCH-021
C               LO(42) =.TRUE. SCHMIDT'S ORTHOGONALISATION OF SOLUTIONS INCH-022
C                              IN USUAL COUPLED EQUATIONS.              INCH-023
C               LO(44) =.TRUE. COULOMB CORRECTIONS.                     INCH-024
C               LO(57) =.TRUE. OUTPUT PHASE-SHIFTS AT EACH ITERATION.   INCH-025
C               LO(74) =.TRUE. OUTPUT OF TIME IN DIFFERENT STEPS.       INCH-026
C               LO(103)=.TRUE. THERE IS A COULOMB SPIN-ORBIT POTENTIAL. INCH-027
C               LO(129)=.TRUE. REAL SPIN-ORBIT OR DIRAC EQUATION.       INCH-028
C               LO(130)=.TRUE. IMAGINARY SPIN-ORBIT OR DIRAC EQUATION.  INCH-029
C               LO(133)=.TRUE. STORE SCALAR AND COULOMB POTENTIAL       INCH-030
C                              INDEPENDENTLY.                           INCH-031
C OUTPUT:    X:       SCATTERING COEFFICIENTS MULTIPLIED KF/KI          INCH-032
C                     REAL PART IN X(IC,IC',5), IMAGINARY PART IN       INCH-033
C                     X(IC,IC',2) FOR INCOMING CHANNEL IC'.             INCH-034
C WORKING AREAS:                                                        INCH-035
C            W:       REAL/IMAGINARY COUPLING POTENTIALS.               INCH-036
C            X:       NUMEROV RECURRENCE.                               INCH-037
C            KR:      WORKING FIELD FOR LINS.                           INCH-038
C                                                                       INCH-039
C FOR THE COMMON  /POTE2/ SEE REDM.                                     INCH-040
C                                                                       INCH-041
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /POTE2/:                     INCH-042
C  ITY(2):    STARTING ADDRESS OF IMAGINARY CENTRAL POTENTIAL.          INCH-043
C  ITY(3):    STARTING ADDRESS OF REAL SPIN-ORBIT POTENTIAL.            INCH-044
C  ITY(4):    STARTING ADDRESS OF IMAGINARY SPIN-ORBIT POTENTIAL.       INCH-045
C  ITY(9):    STARTING ADDRESS OF COULOMB CENTRAL POTENTIAL.            INCH-046
C  ITY(10):   STARTING ADDRESS OF COULOMB SPIN-ORBIT POTENTIAL.         INCH-047
C   USED:     ITY(2),ITY(3),ITY(4),ITY(9),ITY(10).                      INCH-048
C                                                                       INCH-049
C***********************************************************************INCH-050
      IMPLICIT REAL*8 (A-H,O-Z)                                         INCH-051
      LOGICAL LO(150)                                                   INCH-052
      DIMENSION V(ISM,*),MC(KAB,6),NAT(2*LMD,*),AT(LMD,*),NVI(KAB,KAB,2)INCH-053
     1,FAM(KAB,12),Y(KAB,KAB,4),X(KAB,KAB,6),W(2,NC,NC,*),KR(*)         INCH-054
      COMMON /INOUT/ MR,MW,MS                                           INCH-055
      COMMON /POTE2/ ITY(12),INVT,INTV,INSL,NPX                         INCH-056
      IF (LO(74)) CALL HORA                                             INCH-057
      JR=0                                                              INCH-058
      DO 3 I=1,NC                                                       INCH-059
      DO 2 J=1,NC                                                       INCH-060
      DO 1 K=1,6                                                        INCH-061
    1 X(J,I,K)=0.D0                                                     INCH-062
    2 CONTINUE                                                          INCH-063
    3 X(I,I,3)=1.D-15                                                   INCH-064
C  RADIAL INTEGRATION LOOP.                                             INCH-065
      DO 40 JS=1,ISM                                                    INCH-066
      IS=MOD(JS-1,NML)+1                                                INCH-067
      IF (JS.LE.JR) GO TO 26                                            INCH-068
      JI=JR+1                                                           INCH-069
      M=JS-JI                                                           INCH-070
      JR=MIN0(ISM,JR+NML)                                               INCH-071
C THE POTENTIALS ARE FIRST CALCULATED IN W(1,IC,IC',I) AND              INCH-072
C W(2,IC,IC',I) STARTING WITH I=2   FOR IC LARGER OR EQUAL TO IC'.      INCH-073
      DO 19 L=1,NC                                                      INCH-074
      DO 18 J=L,NC                                                      INCH-075
      DO 4 I=JI,JR                                                      INCH-076
      W(1,J,L,I+1)=0.D0                                                 INCH-077
    4 W(2,J,L,I+1)=0.D0                                                 INCH-078
      IF (L.NE.J) GO TO 12                                              INCH-079
C OPTICAL MODEL CONTRIBUTION.                                           INCH-080
      I1=MC(L,4)                                                        INCH-081
      I2=I1+ITY(2)                                                      INCH-082
      DO 5 I=JI,JR                                                      INCH-083
      W(1,J,L,I+1)=FAM(L,8)-FAM(L,10)/DFLOAT(I+M)**2+FAM(L,7)*V(I+M,I1) INCH-084
    5 W(2,J,L,I+1)=FAM(L,7)*V(I+M,I2)                                   INCH-085
      IF (.NOT.LO(129)) GO TO 8                                         INCH-086
      I2=I1+ITY(3)                                                      INCH-087
      DO 6 I=JI,JR                                                      INCH-088
    6 W(1,J,L,I+1)=W(1,J,L,I+1)+FAM(L,9)*V(I+M,I2)                      INCH-089
      IF (.NOT.LO(130)) GO TO 8                                         INCH-090
      I2=I1+ITY(4)                                                      INCH-091
      DO 7 I=JI,JR                                                      INCH-092
    7 W(2,J,L,I+1)=W(2,J,L,I+1)+FAM(L,9)*V(I+M,I2)                      INCH-093
    8 IF (.NOT.LO(133)) GO TO 12                                        INCH-094
      IF (FAM(L,11).EQ.0.D0) GO TO 10                                   INCH-095
      I2=I1+ITY(9)                                                      INCH-096
      DO 9 I=JI,JR                                                      INCH-097
    9 W(1,J,L,I+1)=W(1,J,L,I+1)+FAM(L,11)*V(I+M,I2)                     INCH-098
   10 IF (.NOT.LO(103)) GO TO 12                                        INCH-099
      I2=I1+ITY(10)                                                     INCH-100
      DO 11 I=JI,JR                                                     INCH-101
   11 W(1,J,L,I+1)=W(1,J,L,I+1)+FAM(L,12)*V(I+M,I2)                     INCH-102
C COUPLED CHANNEL CONTRIBUTION.                                         INCH-103
   12 K1=NVI(J,L,1)                                                     INCH-104
      K2=NVI(J,L,2)                                                     INCH-105
      IF (K1.GT.K2) GO TO 18                                            INCH-106
      DO 17 K=K1,K2                                                     INCH-107
      KT=NAT(1,K)                                                       INCH-108
      KU=NAT(2,K)                                                       INCH-109
      DO 13 I=JI,JR                                                     INCH-110
   13 W(1,J,L,I+1)=W(1,J,L,I+1)+AT(2,K)*V(I+M,KT)                       INCH-111
      IF (NAT(2,K).EQ.0) GO TO 15                                       INCH-112
      KU=NAT(2,K)+ITY(2)                                                INCH-113
      DO 14 I=JI,JR                                                     INCH-114
   14 W(2,J,L,I+1)=W(2,J,L,I+1)+AT(2,K)*V(I+M,KU)                       INCH-115
   15 IF (.NOT.LO(133)) GO TO 17                                        INCH-116
      IF (AT(3,K).EQ.0.D0) GO TO 17                                     INCH-117
      KU=KT+ITY(11)                                                     INCH-118
      DO 16 I=JI,JR                                                     INCH-119
   16 W(1,J,L,I+1)=W(1,J,L,I+1)+AT(3,K)*V(I+M,KU)                       INCH-120
   17 CONTINUE                                                          INCH-121
   18 CONTINUE                                                          INCH-122
   19 CONTINUE                                                          INCH-123
      JJ=JR-JI+1                                                        INCH-124
      DO 25 I=1,JJ                                                      INCH-125
      J=I+1                                                             INCH-126
C SYMMETRISATION OF OLD VALUES.                                         INCH-127
      DO 21 K=1,NC                                                      INCH-128
      DO 20 L=K,NC                                                      INCH-129
      W(1,K,L,J)=W(1,L,K,J)                                             INCH-130
   20 W(2,K,L,J)=W(2,L,K,J)                                             INCH-131
   21 CONTINUE                                                          INCH-132
C COMPUTATION OF V +V*V/12. FOR L LARGER OR EQUAL TO K.                 INCH-133
      DO 24 K=1,NC                                                      INCH-134
      DO 23 L=K,NC                                                      INCH-135
      BRE=0.D0                                                          INCH-136
      BIM=0.D0                                                          INCH-137
      DO 22 N=1,NC                                                      INCH-138
      BRE=BRE-W(1,N,K,J)*W(1,N,L,J)+W(2,N,K,J)*W(2,N,L,J)               INCH-139
   22 BIM=BIM-W(2,N,K,J)*W(1,N,L,J)-W(1,N,K,J)*W(2,N,L,J)               INCH-140
C SYMMETRISATION.                                                       INCH-141
      W(1,L,K,I)=W(1,L,K,J)+BRE/12.D0                                   INCH-142
      W(2,L,K,I)=W(2,L,K,J)+BIM/12.D0                                   INCH-143
      W(1,K,L,I)=W(1,L,K,I)                                             INCH-144
   23 W(2,K,L,I)=W(2,L,K,I)                                             INCH-145
      IF (LO(26)) W(1,K,K,I)=W(1,K,K,I)+W(1,K,K,J)**3/360.D0            INCH-146
   24 CONTINUE                                                          INCH-147
   25 CONTINUE                                                          INCH-148
   26 DO 33 I=1,NC                                                      INCH-149
      DO 27 J=1,NC                                                      INCH-150
      X(I,J,1)=X(I,J,2)                                                 INCH-151
      X(I,J,2)=X(I,J,3)                                                 INCH-152
      X(I,J,4)=X(I,J,5)                                                 INCH-153
   27 X(I,J,5)=X(I,J,6)                                                 INCH-154
   28 DO 30 J=1,NC                                                      INCH-155
      HX=0.D0                                                           INCH-156
      HY=0.D0                                                           INCH-157
      DO 29 K=1,NC                                                      INCH-158
      HX=HX+W(1,K,J,IS)*X(I,K,2)-W(2,K,J,IS)*X(I,K,5)                   INCH-159
   29 HY=HY+W(2,K,J,IS)*X(I,K,2)+W(1,K,J,IS)*X(I,K,5)                   INCH-160
      X(I,J,3)=X(I,J,2)+X(I,J,2)-X(I,J,1)-HX                            INCH-161
   30 X(I,J,6)=X(I,J,5)+X(I,J,5)-X(I,J,4)-HY                            INCH-162
      IF (DABS(X(I,I,3)).LT.1.D15) GO TO 33                             INCH-163
      DO 32 JI=1,6                                                      INCH-164
      DO 31 IJ=1,NC                                                     INCH-165
   31 X(I,IJ,JI)=X(I,IJ,JI)*1.D-30                                      INCH-166
   32 CONTINUE                                                          INCH-167
      GO TO 28                                                          INCH-168
   33 CONTINUE                                                          INCH-169
      IF (.NOT.(LO(42).AND.(MOD(JS,JSX).EQ.0).AND.(JS.NE.ISM))) GO TO 40INCH-170
C  SCHMIDT ORTHOGONALISATION PROCEDURE EVERY JSX STEPS.                 INCH-171
      DO 39 I=1,NC                                                      INCH-172
      IF (I.EQ.1) GO TO 37                                              INCH-173
      IM=I-1                                                            INCH-174
      DO 36 K=1,IM                                                      INCH-175
      X(I,K,1)=0.D0                                                     INCH-176
      X(I,K,4)=0.D0                                                     INCH-177
      DO 34 J=1,NC                                                      INCH-178
      X(I,K,1)=X(I,K,1)+X(K,J,2)*X(I,J,2)-X(K,J,5)*X(I,J,5)             INCH-179
   34 X(I,K,4)=X(I,K,4)+X(K,J,5)*X(I,J,2)+X(K,J,2)*X(I,J,5)             INCH-180
      BRE=X(I,K,1)*X(K,K,1)-X(I,K,4)*X(K,K,4)                           INCH-181
      BIM=X(I,K,1)*X(K,K,4)+X(I,K,4)*X(K,K,1)                           INCH-182
      DO 35 J=1,NC                                                      INCH-183
      X(I,J,2)=X(I,J,2)-BRE*X(K,J,2)+BIM*X(K,J,5)                       INCH-184
      X(I,J,5)=X(I,J,5)-BRE*X(K,J,5)-BIM*X(K,J,2)                       INCH-185
      X(I,J,3)=X(I,J,3)-BRE*X(K,J,3)+BIM*X(K,J,6)                       INCH-186
   35 X(I,J,6)=X(I,J,6)-BRE*X(K,J,6)-BIM*X(K,J,3)                       INCH-187
   36 CONTINUE                                                          INCH-188
   37 BRE=0.D0                                                          INCH-189
      BIM=0.D0                                                          INCH-190
      DO 38 J=1,NC                                                      INCH-191
      BRE=BRE+X(I,J,2)**2-X(I,J,5)**2                                   INCH-192
   38 BIM=BIM+2.D0*X(I,J,2)*X(I,J,5)                                    INCH-193
      BIM=BIM/BRE                                                       INCH-194
      X(I,I,1)=1.D0/(BRE*(1.D0+BIM**2))                                 INCH-195
   39 X(I,I,4)=-X(I,I,1)*BIM                                            INCH-196
   40 CONTINUE                                                          INCH-197
      IF (LO(74)) CALL HORA                                             INCH-198
C  MATRICES OF PSEUDO-WRONSKIANS FOR THE MATCHING CONDITION.            INCH-199
      DO 42 I=1,NC                                                      INCH-200
      DO 41 J=1,NC                                                      INCH-201
      X(I,J,2)=(X(I,J,1)*FAM(J,2)-FAM(J,1)*X(I,J,3))                    INCH-202
      X(I,J,5)=(X(I,J,6)*FAM(J,1)-FAM(J,2)*X(I,J,4))                    INCH-203
      X(I,J,1)=(X(I,J,3)*FAM(J,3)-FAM(J,4)*X(I,J,1))                    INCH-204
      X(I,J,4)=(X(I,J,4)*FAM(J,4)-FAM(J,3)*X(I,J,6))                    INCH-205
      IF (FAM(J,8).LT.0.D0) GO TO 41                                    INCH-206
      X(I,J,1)=X(I,J,1)-X(I,J,5)                                        INCH-207
      X(I,J,4)=X(I,J,4)+X(I,J,2)                                        INCH-208
   41 CONTINUE                                                          INCH-209
   42 CONTINUE                                                          INCH-210
C  COMPLEX LINEAR SYSTEM OF EQUATIONS.                                  INCH-211
      NCX=NCIN                                                          INCH-212
      IF (LO(44)) NCX=NC                                                INCH-213
      CALL LINS(X(1,1,4),KAB,X,KAB,X(1,1,5),KAB,X(1,1,2),KAB,NC,NCX,KR,IINCH-214
     1ER)                                                               INCH-215
      IF (IER.NE.0) GO TO 48                                            INCH-216
      IF (.NOT.LO(44)) RETURN                                           INCH-217
C BUILDING THE LINEAR SYSTEM FOR S-MATRIX WITH COULOMB CORRECTIONS.     INCH-218
      DO 45 I=1,NC                                                      INCH-219
      DO 44 J=1,NC                                                      INCH-220
      X(I,J,1)=Y(J,I,1)                                                 INCH-221
      X(I,J,3)=X(I,J,2)                                                 INCH-222
      X(I,J,4)=Y(J,I,2)                                                 INCH-223
      X(I,J,6)=-Y(J,I,1)+X(I,J,5)                                       INCH-224
      DO 43 K=1,NC                                                      INCH-225
      X(I,J,1)=X(I,J,1)-(Y(J,K,1)-Y(J,K,4))*X(I,K,2)+(Y(J,K,3)+Y(J,K,2))INCH-226
     1*X(I,K,5)                                                         INCH-227
      X(I,J,3)=X(I,J,3)-Y(J,K,3)*X(I,K,2)-Y(J,K,1)*X(I,K,5)             INCH-228
      X(I,J,4)=X(I,J,4)-(Y(J,K,1)-Y(J,K,4))*X(I,K,5)-(Y(J,K,2)+Y(J,K,3))INCH-229
     1*X(I,K,2)                                                         INCH-230
   43 X(I,J,6)=X(I,J,6)-Y(J,K,3)*X(I,K,5)+Y(J,K,1)*X(I,K,2)             INCH-231
   44 CONTINUE                                                          INCH-232
   45 X(I,I,4)=X(I,I,4)+1.D0                                            INCH-233
C TRANSFER OF THE SECOND MEMBERS IN X(1,1,5) AND X(1,1,2).              INCH-234
      DO 47 I=1,NC                                                      INCH-235
      DO 46 J=1,NC                                                      INCH-236
      Y(I,J,2)=X(I,J,2)                                                 INCH-237
      X(I,J,2)=X(I,J,3)                                                 INCH-238
      Y(I,J,1)=X(I,J,5)                                                 INCH-239
   46 X(I,J,5)=X(I,J,6)                                                 INCH-240
   47 CONTINUE                                                          INCH-241
      CALL LINS(X(1,1,4),KAB,X,KAB,X(1,1,5),KAB,X(1,1,2),KAB,NC,NCIN,KR INCH-242
     1,IER)                                                             INCH-243
      IF (IER.NE.0) GO TO 48                                            INCH-244
      IF (LO(57)) WRITE (MW,1000) ((I,J,Y(I,J,1),Y(I,J,2),X(I,J,5),X(I,JINCH-245
     1,2),J=1,NCIN),I=1,NC)                                             INCH-246
      RETURN                                                            INCH-247
   48 WRITE (MW,1001) IER                                               INCH-248
      STOP                                                              INCH-249
 1000 FORMAT (/25X,'UNCORRECTED VALUES',30X,'CORRECTED VALUES'/(2X,2I3,4INCH-250
     1D25.15))                                                          INCH-251
 1001 FORMAT (' RETURN FROM LINS WITH IER =',I2///' IN INCH  ...  STOP  INCH-252
     1...')                                                             INCH-253
      END                                                               INCH-254
