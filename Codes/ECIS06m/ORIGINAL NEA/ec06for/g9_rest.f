C 21/08/06                                                      ECIS06  REST-000
      SUBROUTINE REST(KF,NW,DW,IDMX,LO)                                 REST-001
C IF KF=0                                                               REST-002
C IT SAVES ON TAPE MS ALL THE DATA NEEDED TO RE-START A SEARCH STOPPED  REST-003
C BY THE NUMBER OF EVALUATIONS IF LO(35)=.TRUE.                         REST-004
C IF KF.NE.0                                                            REST-005
C IT READS ON TAPE MS TO RE-START THE SEARCH                            REST-006
C IT IS CALLED ONLY IF LO(36)=.TRUE.                                    REST-007
C ******** UNLESS THE MAIN SUBROUTINE IS CHANGED, MS = 8 ***************REST-008
C INPUT:     KE:      0 TO WRITE ON MS, ANY VALUE TO READ ON MS.        REST-009
C            NW,DW:   WORKING ARRAY.                                    REST-010
C            IDMX:    SIZE OF THE WORKING ARRAY.                        REST-011
C            LO(I):   LOGICAL CONTROLS:                                 REST-012
C               LO(35) =.TRUE. SEARCH SAVED ON TAPE 8 IF CORRECTLY ENDEDREST-013
C                              OR STOPPED BY THE NUMBER OF EVALUATIONS. REST-014
C                                                                       REST-015
C THE COMMON /ADDRE/ IS USED IN CALC, CALX, CAL1, VARI, EVAL AND REST.  REST-016
C THE COMMON /ANGUL/ IS USED IN CALX, LECT, SCHE, LCSP, RESU AND REST.  REST-017
C THE COMMON /CONVE/ IS USED IN CALC, CALX, CAL1, INTI, INSI, INTR      REST-018
C                               INRI AND REST.                          REST-019
C THE COMMON /COUPL/ IS USED IN CALC, CALX, LECL, LECT, LECD, COLF,     REST-020
C                               REDM, VIBM, ROTM, ROAM, EXTP, CAL1,     REST-021
C                               POTE, ROTP, ROTD, STDP, VARI, EVAL      REST-022
C                               AND REST.                               REST-023
C THE COMMON /DCHI2/ IS USED IN CALC, CALX, RESU, VARI, EVAL, REST      REST-024
C                               AND FITE.                               REST-025
C THE COMMON /DCONS/ IS USED IN ECIS, CALC, LECL, LECT, COLF, KHCO,     REST-026
C                               CONU, CAL1, POTE, ROTP, ROTZ, STDP,     REST-027
C                               STBF, MTCH, SCAM, SCHE, LSCP, RESU,     REST-028
C                               SCAT AND REST.                          REST-029
C THE COMMON /INTEG/ IS USED IN CALC, CALX, CAL1, VARI, EVAL AND REST.  REST-030
C THE COMMON /NCOMP/ IS USED IN CALC, CALX, LECT, CONU, GGDR, CAL1,     REST-031
C                               QUAN, SCAM, SCHE, RESU, VARI, EVAL,     REST-032
C                               AND REST.                               REST-033
C THE COMMON /POTE1/ IS USED IN CALC, REDM, EXTP, POTE, ROTP, STDP,     REST-034
C                               FOLD AND REST.                          REST-035
C THE COMMON /POTE2/ IS USED IN CALC, LECL, REDM, CAL1, POTE, QUAN,     REST-036
C                               MTCH AND REST.                          REST-037
C THE COMMON /TITRE/ IS USED IN CALC, CALX, RESU, EVAL AND REST.        REST-038
C                                                                       REST-039
C SIGNIFICANCE OF THE QUANTITIES IN COMMON /TITRE/:                     REST-040
C  TITLE(18): TITLE OF THE RUN PRINTED AS HEADING OF RESULTS.           REST-041
C   DEFINED:  TITLE.                                                    REST-042
C   USED:     TITLE.                                                    REST-043
C                                                                       REST-044
C***********************************************************************REST-045
      IMPLICIT REAL*8 (A-H,O-Z)                                         REST-046
      LOGICAL LO(150)                                                   REST-047
      DIMENSION NW(2,*),DW(*),T(18)                                     REST-048
      CHARACTER*4 T,TITLE                                               REST-049
      COMMON /ADDRE/ N1(51)                                             REST-050
      COMMON /ANGUL/ N2(14)                                             REST-051
      COMMON /CONVE/ A3(6)                                              REST-052
      COMMON /COUPL/ N3(7)                                              REST-053
      COMMON /DCHI2/ A4(5)                                              REST-054
      COMMON /DCONS/ A5(7)                                              REST-055
      COMMON /INOUT/ MR,MW,MS                                           REST-056
      COMMON /INTEG/ N7(46)                                             REST-057
      COMMON /NCOMP/ N8(52)                                             REST-058
      COMMON /POTE1/ N9(22)                                             REST-059
      COMMON /POTE2/ NA(12)                                             REST-060
      COMMON /TITRE/ TITLE(18)                                          REST-061
      IDMY=IDMX                                                         REST-062
      REWIND MS                                                         REST-063
      IF (KF.NE.0) GO TO 1                                              REST-064
C WRITES COMMONS FIRST AND W AFTER                                      REST-065
      WRITE (MS) IDMX,LO,N1,N2,A3,N3,A4,A5,N7,N8,N9,NA,TITLE            REST-066
C N1(46) IS NRCO,N7(2) IS NPLACE,N7(1) IS IDMT                          REST-067
      I1=N1(46)                                                         REST-068
      WRITE (MW,1000) MS,I1                                             REST-069
      WRITE (MS) (DW(I),I=1,I1)                                         REST-070
      REWIND MS                                                         REST-071
      RETURN                                                            REST-072
C READ COMMONS,COMPUTES LIMITS OF W AND READ W                          REST-073
    1 READ (MS) I,LO,N1,N2,A3,N3,A4,A5,N7,N8,N9,NA,T                    REST-074
      I1=N1(46)                                                         REST-075
      N7(1)=IDMY                                                        REST-076
      WRITE (MW,1001) TITLE,T,N7(1)                                     REST-077
      IF (N7(2).GT.N7(1)) CALL MEMO('REST',N7(1),N7(2))                 REST-078
      READ (MS) (DW(I),I=1,I1)                                          REST-079
      REWIND MS                                                         REST-080
C READS ON 5 IF SAVE HAS TO BE DONE AGAIN AND THE NEW NUMBER OF EVAL.   REST-081
C ANY MODIFICATION OF THE SEARCH CAN BE READ AT THIS PLACE IF THEY      REST-082
C DO NOT SPOIL THE SEARCH                                               REST-083
      LO(35)=.FALSE.                                                    REST-084
      READ (MR,1002) LO(35),N,ECH,RAP                                   REST-085
C N1(20) IS NIW, N7(27) IS NREC                                         REST-086
      IF (N.EQ.0) N=N7(34)+2                                            REST-087
      M=N1(20)                                                          REST-088
      NW(2,M)=NW(2,M)+N                                                 REST-089
      WRITE (MW,1003) MS,I1                                             REST-090
      WRITE (MW,1004) N,NW(2,M)                                         REST-091
      IF (LO(35)) WRITE (MW,1005) MS                                    REST-092
      IF (.NOT.LO(35)) WRITE (MW,1006) MS                               REST-093
      IF (ECH.LT.1.D0) GO TO 2                                          REST-094
      WRITE (MW,1008) A4(3),ECH                                         REST-095
      A4(3)=ECH                                                         REST-096
    2 IF (RAP.LT.1.D0) GO TO 3                                          REST-097
      WRITE (MW,1009) A4(4),RAP                                         REST-098
      A4(4)=RAP                                                         REST-099
    3 RETURN                                                            REST-100
 1000 FORMAT (///15H OUTPUT ON TAPE,I3,24H OF DW FROM DW(1) TO DW(,I6,1HREST-101
     1)///)                                                             REST-102
 1001 FORMAT ('1'//1X,18A4//' RESTART OF COMPUTATION'//' LAST TITLE : ''REST-103
     1''',18A4,''''''//' AVAILABLE WORKING SPACE',I8/)                  REST-104
 1002 FORMAT (L1,I4,5X,2F10.5)                                          REST-105
 1003 FORMAT (' INPUT FROM TAPE',I3,' OF DW FROM DW(1) TO DW(',I6,')'//)REST-106
 1004 FORMAT (/' MAXIMUM NUMBER OF EVALUATIONS INCREASED BY',I4,5X,'NEW REST-107
     1VALUE',I6)                                                        REST-108
 1005 FORMAT (' SAVE ON TAPE',I3,' IF NECESSARY.')                      REST-109
 1006 FORMAT (' NO SAVE ON TAPE',I3)                                    REST-110
 1008 FORMAT (' CHANGE OF ECH FROM',F10.5,' TO',F10.5)                  REST-111
 1009 FORMAT (' CHANGE OF RAP FROM',F10.5,' TO',F10.5)                  REST-112
      END                                                               REST-113
