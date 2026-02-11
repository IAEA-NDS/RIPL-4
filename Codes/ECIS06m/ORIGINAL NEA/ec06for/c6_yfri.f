C 13/11/05                                                      ECIS06  YFRI-000
      SUBROUTINE YFRI(ETA,RAU,FO,FPO,GO,GPO,IDIV,SIGMA)                 YFRI-001
C COMPUTATION OF THE COULOMB FUNCTIONS FOR L=0.                         YFRI-002
C RICCATI METHOD.  SAME ARGUMENTS AS FCZ0.                              YFRI-003
C INPUT:     ETA,RAU: COULOMB PARAMETER AND RADIUS.                     YFRI-004
C            SIGMA:   COULOMB PHASE SHIFT FOR L=0.                      YFRI-005
C OUTPUT:    FO,FPO:  REGULAR FUNCTION AND DERIVATIVE.                  YFRI-006
C            GO,GPO:  IRREGULAR FUNCTION AND DERIVATIVE.                YFRI-007
C            IDIV:    POWER OF 10.                                      YFRI-008
C***********************************************************************YFRI-009
      IMPLICIT REAL*8 (A-H,O-Z)                                         YFRI-010
      DIMENSION Q(5),QP(5)                                              YFRI-011
C        RICCATI COEFFICIENTS.                                          YFRI-012
      DATA G61,G62,G63,G64,G65,G66,G67,G68,G69,G610,G611/ 1.159057617187YFRI-013
     15D-2,3.863525390625D-2,4.6600341796875D-2,4.8583984375D-2,1.156514YFRI-014
     248567708D0,5.6874755859375D0,1.323888288225445D1,1.713083224826384YFRI-015
     3D1,1.269003295898436D1,5.05523681640625D0,8.42539464010415D-1/    YFRI-016
      DATA G81,G82,G83,G84,G85,G86,G87,G88,G89,G810,G811,G812,G813,G814,YFRI-017
     1G815/ 1.851092066083633D-2,8.63842964172363D-2,1.564757823944092D-YFRI-018
     21,1.430139541625977D-1,1.924622058868408D-1,8.500803152720129D0,7.YFRI-019
     3265429720878595D1,3.057942376817972D2,7.699689544836672D2,1.254157YFRI-020
     4054424285D3,1.361719536066055D3,9.831831171035763D2,4.547869927883YFRI-021
     5148D2,1.222640538215636D2,1.455524450256709D1/                    YFRI-022
      DATA GP61,GP62,GP63,GP64,GP65,GP66/ 0.289764404296875D-1,0.2318115YFRI-023
     1234375D0,0.8056640625D0,1.6015625D0,0.3046875D0,5.625D0/          YFRI-024
      DATA GP81,GP82,GP83,GP84,GP85,GP86,GP87,GP88/ 0.647882223129272D-1YFRI-025
     1,0.6910743713378906D0,0.3322952270507811D1,0.94830322265625D1,17.6YFRI-026
     296533203125D0,34.787109375D0,50.203125D0,78.75D0/                 YFRI-027
      DATA Q/ 0.4959570165D-1,0.8888888889D-2,0.2455199181D-2,0.91089580YFRI-028
     161D-3,0.2534684115D-3/                                            YFRI-029
      DATA QP,HO,HPO/ 0.1728260369D0,0.3174603174D-3,0.3581214850D-2,0.3YFRI-030
     1117824680D-3,0.9073966427D-3,2*0.D0/                              YFRI-031
      ETAC=ETA*ETA                                                      YFRI-032
      ETA2=ETA+ETA                                                      YFRI-033
      ETARO=ETA*RAU                                                     YFRI-034
      IND=0                                                             YFRI-035
      JND=0                                                             YFRI-036
      IG=0                                                              YFRI-037
      IDIV=0                                                            YFRI-038
      IF (ETA.GT.0.D0) GO TO 1                                          YFRI-039
      IF (-ETARO-14.0625D0) 3 , 15 , 15                                 YFRI-040
    1 IF (DABS(RAU-ETA2).LE.1.D-9) GO TO 14                             YFRI-041
      IF (RAU-ETA2) 6 , 14 , 2                                          YFRI-042
    2 IF (RAU-ETA2-2.D1*(ETA**0.25D0)) 4 , 15 , 15                      YFRI-043
    3 NN=1                                                              YFRI-044
      GO TO 5                                                           YFRI-045
    4 NN=0                                                              YFRI-046
    5 CALL YFCL(ETA,RAU,FO,FPO,GO,GPO,SIGMA,IDIV,NN)                    YFRI-047
      RETURN                                                            YFRI-048
    6 IF (ETARO.LE.12.D0) GO TO 3                                       YFRI-049
      TRA=ETA2-6.75D0*(ETA**0.4D0)                                      YFRI-050
      IF (RAU.LE.TRA) GO TO 7                                           YFRI-051
      IND=1                                                             YFRI-052
      JND=1                                                             YFRI-053
      RO=RAU                                                            YFRI-054
      RAU=TRA                                                           YFRI-055
      RAU0=TRA                                                          YFRI-056
C             RICCATI   1.                                              YFRI-057
    7 X=RAU/ETA2                                                        YFRI-058
      U=(1.D0-X)/X                                                      YFRI-059
      X2=X*X                                                            YFRI-060
      RU=DSQRT(U)                                                       YFRI-061
      RX=DSQRT(X)                                                       YFRI-062
      TRE=1.D0/(U*RU*ETA2)                                              YFRI-063
      TRB=TRE*TRE                                                       YFRI-064
      FI=(DSQRT((1.D0-X)*X)+DATAN2(RX,DSQRT(1.D0-RX*RX))-1.5707963267948YFRI-065
     197D0)*ETA2                                                        YFRI-066
      TR1=-0.25D0*DLOG(U)                                               YFRI-067
      TR2=-((9.D0*U+6.D0)*U+5.D0)/48.D0                                 YFRI-068
      TR3=((((-3.D0*U-4.D0)*U+6.D0)*U+12.D0)*U+5.D0)/64.D0              YFRI-069
      TR4=-((((((U+2.D0)*945.D0*U+1395.D0)*U+12300.D0)*U+25191.D0)*U+198YFRI-070
     190.D0)*U+5525.D0)/46080.D0                                        YFRI-071
      TR5=((((((((-27.D0*U-72.D0)*U-68.D0)*U+360.D0)*U+2190.D0)*U+4808.DYFRI-072
     10)*U+5148.D0)*U+2712.D0)*U+565.D0)/2048.D0                        YFRI-073
      TR6=-((((((((((G61*U+G62)*U+G63)*U+G64)*U+G65)*U+G66)*U+G67)*U+G68YFRI-074
     1)*U+G69)*U+G610)*U+G611)                                          YFRI-075
      TR7=((((((((((((-81.D0*U-324.D0)*U-486.D0)*U-404.D0)*U+4509.D0)*U+YFRI-076
     152344.D0)*U+233436.D0)*U+567864.D0)*U+838521.D0)*U+775884.D0)*U+44YFRI-077
     21450.D0)*U+141660.D0)*U+19675.D0)/6144.D0                         YFRI-078
      TR8=(((((((((((((G81*U+G82)*U+G83)*U+G84)*U+G85)*U+G86)*U+G87)*U+GYFRI-079
     188)*U+G89)*U+G810)*U+G811)*U+G812)*U+G813)*U+G814)*U+G815         YFRI-080
      FI=FI+TRE*(TR2+TRB*(TR4+TRB*(TR6+TRB*TR8)))                       YFRI-081
      PSI=-FI                                                           YFRI-082
      TRA=TR1+TRB*(TR3+TRB*(TR5+TRB*TR7))                               YFRI-083
      FI=FI+TRA                                                         YFRI-084
      PSI=PSI+TRA                                                       YFRI-085
      FIP=RU*ETA2                                                       YFRI-086
      TRA=1.D0/(X2*U)                                                   YFRI-087
      TR1=0.25D0                                                        YFRI-088
      TRE=TRE/(X2*X2*U)                                                 YFRI-089
      TRB=TRB/(X2*X2)                                                   YFRI-090
      TR2=-(8.D0*X-3.D0)/32.D0                                          YFRI-091
      TR3=((24.D0*X-12.D0)*X+3.D0)/64.D0                                YFRI-092
      TR4=(((-1536.D0*X+704.D0)*X-336.D0)*X+63.D0)/2048.D0              YFRI-093
      TR5=((((1920.D0*X-576.D0)*X+504.D0)*X-180.D0)*X+27.D0)/1024.D0    YFRI-094
      TR6=((((-GP66*X+GP65)*X-GP64)*X+GP63)*X-GP62)*X+GP61              YFRI-095
      TR7=-((((((-40320.D0*X-10560.D0)*X-13248.D0)*X+7560.D0)*X-3132.D0)YFRI-096
     1*X+756.D0)*X-81.D0)/2048.D0                                       YFRI-097
      TR8=-(((((((GP88*X+GP87)*X+GP86)*X-GP85)*X+GP84)*X-GP83)*X+GP82)*XYFRI-098
     1-GP81)                                                            YFRI-099
      FIP=FIP+TRE*(TR2+TRB*(TR4+TRB*(TR6+TRB*TR8)))                     YFRI-100
      PSIP=-FIP                                                         YFRI-101
      TRA=TRA*(TR1+TRB*(TR3+TRB*(TR5+TRB*TR7)))                         YFRI-102
      FIP=FIP+TRA                                                       YFRI-103
      PSIP=PSIP+TRA                                                     YFRI-104
      XXX=34.588776394910686D0                                          YFRI-105
      INDG=IDINT(PSI/XXX)                                               YFRI-106
      IDIV=15*INDG                                                      YFRI-107
      IF (INDG.EQ.0) GO TO 8                                            YFRI-108
      PSI=PSI-XXX*INDG                                                  YFRI-109
      FI=FI+XXX*INDG                                                    YFRI-110
    8 FO=0.5D0*DEXP(FI)                                                 YFRI-111
      GO=DEXP(PSI)                                                      YFRI-112
      FPO=FO*FIP/ETA2                                                   YFRI-113
      GPO=GO*PSIP/ETA2                                                  YFRI-114
      IF (JND.EQ.0) RETURN                                              YFRI-115
      RAU=RO                                                            YFRI-116
      GO=FO                                                             YFRI-117
      GPO=FPO                                                           YFRI-118
    9 X=RAU0-RO                                                         YFRI-119
      X2=X*X                                                            YFRI-120
      X3=X*X2                                                           YFRI-121
      UNR=1.D0/RAU0                                                     YFRI-122
      ETR0=1.D0-2.D0*ETA*UNR                                            YFRI-123
      U0=GO                                                             YFRI-124
      U1=-X*GPO                                                         YFRI-125
      U2=-0.5D0*ETR0*X2*U0                                              YFRI-126
      S=U0+U1+U2                                                        YFRI-127
      V1=U1/X                                                           YFRI-128
      V2=2.D0*U2/X                                                      YFRI-129
      T=V1+V2                                                           YFRI-130
      XN=3.D0                                                           YFRI-131
      DO 11 N=3,10000                                                   YFRI-132
      XN1=XN-1.D0                                                       YFRI-133
      XN1=XN*XN1                                                        YFRI-134
      U3=X*U2*UNR*(1.D0-2.D0/XN)-ETR0*U1*X2/XN1+X3*U0*UNR/XN1           YFRI-135
      S=S+U3                                                            YFRI-136
      V3=XN*U3/X                                                        YFRI-137
      T=T+V3                                                            YFRI-138
      IF (DABS(U3).GT.1.D-10*DABS(S)) GO TO 10                          YFRI-139
      IF (DABS(V3).LE.1.D-10*DABS(T)) GO TO 12                          YFRI-140
   10 U0=U1                                                             YFRI-141
      U1=U2                                                             YFRI-142
      U2=U3                                                             YFRI-143
   11 XN=XN+1.D0                                                        YFRI-144
   12 IF (IG.EQ.0) GO TO 13                                             YFRI-145
      GO=S                                                              YFRI-146
      GPO=-T                                                            YFRI-147
      FO=HO                                                             YFRI-148
      FPO=HPO                                                           YFRI-149
      RETURN                                                            YFRI-150
   13 HO=S                                                              YFRI-151
      HPO=-T                                                            YFRI-152
   14 ET0=ETA**(0.1666666666666667D0)                                   YFRI-153
      ETAD=ETAC*ETAC                                                    YFRI-154
      ET=ETA**(0.6666666666666667D0)                                    YFRI-155
      ET1=ET*ET                                                         YFRI-156
      ET2=ET1*ET1                                                       YFRI-157
      ET3=ET2*ET                                                        YFRI-158
      ET4=ETAD*ET                                                       YFRI-159
      ET5=ET4*ET                                                        YFRI-160
      FO=1.D0-Q(1)/ET1-Q(2)/ETAC-Q(3)/ET3-Q(4)/ETAD-Q(5)/ET5            YFRI-161
      GO=1.D0+Q(1)/ET1-Q(2)/ETAC+Q(3)/ET3-Q(4)/ETAD+Q(5)/ET5            YFRI-162
      FPO=1.D0+QP(1)/ET+QP(2)/ETAC+QP(3)/ET2+QP(4)/ETAD+QP(5)/ET4       YFRI-163
      GPO=1.D0-QP(1)/ET+QP(2)/ETAC-QP(3)/ET2+QP(4)/ETAD-QP(5)/ET4       YFRI-164
      FO=0.7063326373D0*ET0*FO                                          YFRI-165
      GO=1.223404016D0*ET0*GO                                           YFRI-166
      FPO=0.4086957323D0*FPO/ET0                                        YFRI-167
      GPO=-0.7078817734D0*GPO/ET0                                       YFRI-168
      IDIV=0                                                            YFRI-169
      IF (IND.EQ.0) RETURN                                              YFRI-170
      IG=1                                                              YFRI-171
      RAU0=ETA2                                                         YFRI-172
      GO TO 9                                                           YFRI-173
C        RICCATI 2 AND 3.                                               YFRI-174
   15 X=ETA2/RAU                                                        YFRI-175
      X2=X*X                                                            YFRI-176
      U=1.D0-X                                                          YFRI-177
      RU=DSQRT(U)                                                       YFRI-178
      U3=U*U*U                                                          YFRI-179
      TRD=1.D0/(U3*ETA2*ETA2)                                           YFRI-180
      TRC=X2*TRD                                                        YFRI-181
      TRE=1.D0/(U*RU*ETA2)                                              YFRI-182
      FI=-0.25D0*DLOG(U)                                                YFRI-183
      TRB=TRD/64.D0                                                     YFRI-184
      TR3=(((3.D0*U-4.D0)*U-6.D0)*U+12.D0)*U-5.D0                       YFRI-185
      TR5=((((((((-27.D0*U+72.D0)*U-68.D0)*U-360.D0)*U+2190.D0)*U-4808.DYFRI-186
     10)*U+5148.D0)*U-2712.D0)*U+565.D0)/32.D0                          YFRI-187
      TR7=((((((((((((81.D0*U-324.D0)*U+486.D0)*U-404.D0)*U-4509.D0)*U+5YFRI-188
     12344.D0)*U-233436.D0)*U+567864.D0)*U-838521.D0)*U+775884.D0)*U-441YFRI-189
     2450.D0)*U+141660.D0)*U-19675.D0)/96.D0                            YFRI-190
      FI=FI+TRB*(TR3+TRD*(TR5+TRD*TR7))                                 YFRI-191
      FIP=0.25D0/U                                                      YFRI-192
      TRB=3.D0*TRC/(64.D0*U)                                            YFRI-193
      TR3=(4.D0-X)*X-8.D0                                               YFRI-194
      TR5=((((9.D0*X-60.D0)*X+168.D0)*X-192.D0)*X+640.D0)/16.D0         YFRI-195
      TR7=((((((-27.D0*X+252.D0)*X-1044.D0)*X+2520.D0)*X-4416.D0)*X-3520YFRI-196
     1.D0)*X-13440.D0)/32.D0                                            YFRI-197
      FIP=FIP+TRB*(TR3+TRC*(TR5+TRC*TR7))                               YFRI-198
      TRA=DABS((RU-1.D0)/(RU+1.D0))                                     YFRI-199
      PSI=(0.5D0*DLOG(TRA)+RU/X)*ETA2+0.785398163397448D0               YFRI-200
      TR2=-((9.D0*U-6.D0)*U+5.D0)/48.D0                                 YFRI-201
      TR4=((((((U-2.D0)*945.D0*U+1395.D0)*U-12300.D0)*U+25191.D0)*U-1989YFRI-202
     10.D0)*U+5525.D0)/46080.D0                                         YFRI-203
      TR6=(((((((((-G61*U+G62)*U-G63)*U+G64)*U-G65)*U+G66)*U-G67)*U+G68)YFRI-204
     1*U-G69)*U+G610)*U-G611                                            YFRI-205
      TR8=(((((((((((((G81*U-G82)*U+G83)*U-G84)*U+G85)*U-G86)*U+G87)*U-GYFRI-206
     188)*U+G89)*U-G810)*U+G811)*U-G812)*U+G813)*U-G814)*U+G815         YFRI-207
      PSI=PSI+TRE*(TR2+TRD*(TR4+TRD*(TR6+TRD*TR8)))                     YFRI-208
      PSIP=-RU*ETA2/X2                                                  YFRI-209
      TRB=TRE*X/U                                                       YFRI-210
      TR2=(3.D0*X-8.D0)/32.D0                                           YFRI-211
      TR4=-(((63.D0*X-336.D0)*X+704.D0)*X-1536.D0)/2048.D0              YFRI-212
      TR6=((((GP61*X-GP62)*X+GP63)*X-GP64)*X+GP65)*X-GP66               YFRI-213
      TR8=((((((-GP81*X+GP82)*X-GP83)*X+GP84)*X-GP85)*X+GP86)*X+GP87)*X+YFRI-214
     1GP88                                                              YFRI-215
      PSIP=PSIP+TRB*(TR2+TRC*(TR4+TRC*(TR6+TRC*TR8)))                   YFRI-216
      TRA=DEXP(FI)                                                      YFRI-217
      FO=TRA*DSIN(PSI)                                                  YFRI-218
      GO=TRA*DCOS(PSI)                                                  YFRI-219
      IF (ETA.GT.0.D0) GO TO 16                                         YFRI-220
      TRA=FO                                                            YFRI-221
      FO=-GO                                                            YFRI-222
      GO=TRA                                                            YFRI-223
   16 TRA=-ETA2/(RAU*RAU)                                               YFRI-224
      FPO=(FIP*FO+PSIP*GO)*TRA                                          YFRI-225
      GPO=(FIP*GO-PSIP*FO)*TRA                                          YFRI-226
      RETURN                                                            YFRI-227
      END                                                               YFRI-228
