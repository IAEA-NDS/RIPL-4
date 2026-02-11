C 07/02/06                                                      ECIS06  HORA-000
      SUBROUTINE HORA                                                   HORA-001
C GIVES THE TIME ELAPSED SINCE THE FIRST CALL AND SINCE THE LAST CALL.  HORA-002
C***********************************************************************HORA-003
      DIMENSION A(2)                                                    HORA-004
      COMMON /INOUT/ MR,MW,MS                                           HORA-005
      DATA M0,M2 /2*0/                                                  HORA-006
      SAVE M0,M2                                                        HORA-007
      M1=INT(1000.*ETIME(A))                                            HORA-008
      N5=M1-M2                                                          HORA-009
      NT=M1-M0                                                          HORA-010
      M2=M1                                                             HORA-011
      NS=NT/1000                                                        HORA-012
      NT=NT-1000*NS                                                     HORA-013
      NM=NS/60                                                          HORA-014
      NS=NS-60*NM                                                       HORA-015
      NH=NM/60                                                          HORA-016
      NM=NM-60*NH                                                       HORA-017
      WRITE (MW,1000) NH,NM,NS,NT,N5                                    HORA-018
      RETURN                                                            HORA-019
 1000 FORMAT (' *** TOTAL TIME ***',I3,'H',I3,'MN',I3,'S',I4,' MS',10X,'HORA-020
     1DIFFERENCE SINCE LAST CALL',I9,' MILLISECONDS.')                  HORA-021
      END                                                               HORA-022
