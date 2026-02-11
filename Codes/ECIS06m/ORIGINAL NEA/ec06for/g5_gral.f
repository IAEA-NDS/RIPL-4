C 28/06/06                                                      ECIS06  GRAL-000
      SUBROUTINE GRAL(TGR,GRR,FRR,MTHETA,MF,CM,ID1,LT1,LT2)             GRAL-001
C PLOTS OF CROSS-SECTION AND POLARISATIONS WITH AND WITHOUT EXPERIMENTALGRAL-002
C DATA.                                                                 GRAL-003
C INPUT:     TGR:     ANGLES.                                           GRAL-004
C            GRR:     CALCULATED VALUES.                                GRAL-005
C            FRR:     EXPERIMENTAL VALUES.                              GRAL-006
C            MTHETA:  NUMBER OF ANGLES.                                 GRAL-007
C            MF,CM:   EQUIVALENT BY CALL, DESCRIPTION OF DATA (SEE OBSE)GRAL-008
C                     MF(5,J) IS THE NUMBER OF POWER OF 10 IN A LINE    GRAL-009
C                     FOR CROSS SECTIONS AND MUST BE 1 TO PLOT A        GRAL-010
C                     POLARISATION (ALWAYS BETWEEN -1 AND 1).           GRAL-011
C                     CM(I,J) IS A TITLE FOR I=6,10.                    GRAL-012
C            ID1:     NUMBER OF DIFFERENT POLARISATIONS.                GRAL-013
C            LT1:     .TRUE. FOR POLARISATION,                          GRAL-014
C                     .FALSE. FOR CROSS SECTION.                        GRAL-015
C            LT2:     .TRUE. IF THERE ARE NO EXPERIMENTAL DATA.         GRAL-016
C***********************************************************************GRAL-017
      IMPLICIT REAL*8 (A-H,O-Z)                                         GRAL-018
      LOGICAL LT1,LT2                                                   GRAL-019
      DIMENSION TGR(MTHETA),GRR(ID1,MTHETA),FRR(MTHETA),MF(10,*)        GRAL-020
      CHARACTER*1 VGR(112),GRA(17)                                      GRAL-021
      CHARACTER*4 CM(10,*)                                              GRAL-022
      COMMON /INOUT/ MR,MW,MS                                           GRAL-023
      DATA GRA /'+','*',',',':','-','1','2','3','4','5','6','7','8','9',GRAL-024
     1'0',' ','.'/                                                      GRAL-025
      IF (LT1) GO TO 14                                                 GRAL-026
C CROSS SECTIONS.                                                       GRAL-027
      IF ((MF(5,1).LE.0).OR.(MF(5,1).GT.50)) RETURN                     GRAL-028
      NA2=100/MF(5,1)                                                   GRAL-029
      AA2=DFLOAT(NA2)                                                   GRAL-030
      WRITE (MW,1000)                                                   GRAL-031
      WRITE (MW,1001) GRA(1),(CM(J,1),J=6,10)                           GRAL-032
      IF (.NOT.LT2) WRITE (MW,1002)                                     GRAL-033
      AA3=AA2*DFLOAT((MF(5,1)+1)/2)                                     GRAL-034
      AA1=AA2/DLOG(10.D0)                                               GRAL-035
      NA1=101-MF(5,1)*NA2                                               GRAL-036
      A5=101.5D0                                                        GRAL-037
      DO 1 I=1,112                                                      GRAL-038
    1 VGR(I)=GRA(16)                                                    GRAL-039
      DO 2 I=NA1,112,NA2                                                GRAL-040
    2 VGR(I)=GRA(6)                                                     GRAL-041
      WRITE (MW,1003) VGR                                               GRAL-042
      DO 3 I=1,112                                                      GRAL-043
    3 VGR(I)=GRA(17)                                                    GRAL-044
      DO 4 I=NA1,112,NA2                                                GRAL-045
    4 VGR(I)=GRA(2)                                                     GRAL-046
      WRITE (MW,1003) VGR                                               GRAL-047
      DO 11 I=1,MTHETA                                                  GRAL-048
      DO 5 J=2,111                                                      GRAL-049
    5 VGR(J)=GRA(16)                                                    GRAL-050
      DO 6 J=NA1,111,NA2                                                GRAL-051
    6 VGR(J)=GRA(17)                                                    GRAL-052
      IF (GRR(1,I).LT.0.D0) GO TO 27                                    GRAL-053
      A4=AA1*DLOG(GRR(1,I))+A5                                          GRAL-054
      NG=0                                                              GRAL-055
    7 IF (NG.GT.20) GO TO 28                                            GRAL-056
      NG=NG+1                                                           GRAL-057
      IF (A4.GT.2.D0) GO TO 8                                           GRAL-058
      A4=A4+AA3                                                         GRAL-059
      A5=A5+AA3                                                         GRAL-060
      WRITE (MW,1004)                                                   GRAL-061
      GO TO 7                                                           GRAL-062
    8 IF (A4.LT.112.D0) GO TO 9                                         GRAL-063
      A4=A4-AA3                                                         GRAL-064
      A5=A5-AA3                                                         GRAL-065
      WRITE (MW,1005)                                                   GRAL-066
      GO TO 8                                                           GRAL-067
    9 L=IDINT(A4)                                                       GRAL-068
      VGR(L)=GRA(1)                                                     GRAL-069
      IF (LT2) GO TO 10                                                 GRAL-070
C EXPERIMENTAL DATA.                                                    GRAL-071
      IF (FRR(I).LT.0.D0) GO TO 29                                      GRAL-072
      L=IDINT(AA1*DLOG(FRR(I))+A5)                                      GRAL-073
      IF (L.GT.1.AND.L.LT.112) VGR(L)=GRA(2)                            GRAL-074
   10 WRITE (MW,1006) TGR(I),VGR                                        GRAL-075
   11 CONTINUE                                                          GRAL-076
      DO 12 I=1,112                                                     GRAL-077
   12 VGR(I)=GRA(17)                                                    GRAL-078
      DO 13 I=NA1,112,NA2                                               GRAL-079
   13 VGR(I)=GRA(2)                                                     GRAL-080
      WRITE (MW,1003) VGR                                               GRAL-081
      RETURN                                                            GRAL-082
C POLARISATION.                                                         GRAL-083
   14 IF (LT2) GO TO 15                                                 GRAL-084
      WRITE (MW,1000)                                                   GRAL-085
      WRITE (MW,1001) GRA(1),(CM(J,1),J=6,10)                           GRAL-086
      WRITE (MW,1002)                                                   GRAL-087
      GO TO 17                                                          GRAL-088
   15 NK=0                                                              GRAL-089
      DO 16 K=1,ID1                                                     GRAL-090
      IF (MF(5,K).NE.1) GO TO 16                                        GRAL-091
      KN=MOD(NK,15)+1                                                   GRAL-092
      NK=NK+1                                                           GRAL-093
      IF (NK.EQ.1) WRITE (MW,1000)                                      GRAL-094
      WRITE (MW,1001) GRA(KN),(CM(J,K),J=6,10)                          GRAL-095
   16 CONTINUE                                                          GRAL-096
      IF (NK.EQ.0) RETURN                                               GRAL-097
   17 WRITE (MW,1007)                                                   GRAL-098
      DO 18 I=1,112                                                     GRAL-099
   18 VGR(I)=GRA(16)                                                    GRAL-100
      DO 19 I=1,103                                                     GRAL-101
   19 VGR(I)=GRA(17)                                                    GRAL-102
      DO 20 I=2,102,10                                                  GRAL-103
   20 VGR(I)=GRA(2)                                                     GRAL-104
      WRITE (MW,1003) VGR                                               GRAL-105
      DO 24 I=1,MTHETA                                                  GRAL-106
      DO 21 M=2,102                                                     GRAL-107
   21 VGR(M)=GRA(16)                                                    GRAL-108
      VGR(52)=GRA(17)                                                   GRAL-109
      NK=0                                                              GRAL-110
      DO 22 K=1,ID1                                                     GRAL-111
      IF (MF(5,K).NE.1) GO TO 22                                        GRAL-112
      NK=MOD(NK,15)+1                                                   GRAL-113
      M=IDINT(50.D0*GRR(K,I)+52.5D0)                                    GRAL-114
      IF (M.GT.1.AND.M.LE.102) VGR(M)=GRA(NK)                           GRAL-115
   22 CONTINUE                                                          GRAL-116
      IF (LT2) GO TO 23                                                 GRAL-117
      M=IDINT(50.D0*FRR(I)+52.5D0)                                      GRAL-118
      IF (M.GT.1.AND.M.LE.102) VGR(M)=GRA(2)                            GRAL-119
   23 WRITE (MW,1006) TGR(I),VGR                                        GRAL-120
   24 CONTINUE                                                          GRAL-121
      DO 25 I=1,103                                                     GRAL-122
   25 VGR(I)=GRA(17)                                                    GRAL-123
      DO 26 I=2,102,10                                                  GRAL-124
   26 VGR(I)=GRA(2)                                                     GRAL-125
      WRITE (MW,1003) VGR                                               GRAL-126
      RETURN                                                            GRAL-127
   27 WRITE (MW,1008) I,GRR(1,I)                                        GRAL-128
      GO TO 30                                                          GRAL-129
   28 WRITE (MW,1009) I,GRR(1,I)                                        GRAL-130
      GO TO 30                                                          GRAL-131
   29 WRITE (MW,1010) I,FRR(I)                                          GRAL-132
   30 WRITE (MW,1011)                                                   GRAL-133
      RETURN                                                            GRAL-134
 1000 FORMAT ('1')                                                      GRAL-135
 1001 FORMAT (45X,A1,4X,5A4/)                                           GRAL-136
 1002 FORMAT (45X,'+    CALCULATED VALUE'/45X,'*    EXPERIMENTAL VALUE'/GRAL-137
     1)                                                                 GRAL-138
 1003 FORMAT (7X,112A1)                                                 GRAL-139
 1004 FORMAT (40(' -*'))                                                GRAL-140
 1005 FORMAT (40(' *-'))                                                GRAL-141
 1006 FORMAT (F7.2,112A1)                                               GRAL-142
 1007 FORMAT (//6X,' -1',7X,' -.8',6X,' -.6',6X,' -.4',6X,' -.2',7X,' 0'GRAL-143
     1,8X,' .2',7X,' .4',7X,' .6',7X,' .8',7X,' 1')                     GRAL-144
 1008 FORMAT (' FOR THE',I4,'TH LINE, CALCULATED VALUE =',D15.8)        GRAL-145
 1009 FORMAT (' FOR THE',I4,'TH LINE, MORE THAN 20 TRANSLATIONS. CALCULAGRAL-146
     1TED VALUE =',D15.8)                                               GRAL-147
 1010 FORMAT (' FOR THE',I4,'TH LINE, EXPERIMENTAL VALUE =',D15.8)      GRAL-148
 1011 FORMAT (//' **** GRAPH CANCELED ****')                            GRAL-149
      END                                                               GRAL-150
