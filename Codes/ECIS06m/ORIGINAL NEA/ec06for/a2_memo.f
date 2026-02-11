C 19/11/05                                                      ECIS06  MEMO-000
      SUBROUTINE MEMO(NAME,IDMT,NPLACE)                                 MEMO-001
C CONTROL OF MEMORY AREA.                                               MEMO-002
C INPUT:     NAME:    NAME OF CALLING SUBROUTINE.                       MEMO-003
C            IDMT:    PREVIOUS SPACE.                                   MEMO-004
C            NPLACE:  REQUESTED SPACE.                                  MEMO-005
C                                                                       MEMO-006
C THIS SUBROUTINE CAN BE REPLACED BY A SUBROUTINE INCREASING THE COMMON MEMO-007
C IF IT IS POSSIBLE.                                                    MEMO-008
C***********************************************************************MEMO-009
      CHARACTER*4 NAME                                                  MEMO-010
      COMMON /INOUT/ MR,MW,MS                                           MEMO-011
      IF (IDMT.GE.NPLACE) RETURN                                        MEMO-012
      WRITE (MW,1000) IDMT,NPLACE,NAME                                  MEMO-013
      STOP                                                              MEMO-014
 1000 FORMAT (' NOT ENOUGH PLACE .....   ',I10,' MEMORIES ALLOWED',I10,'MEMO-015
     1 MEMORIES REQUESTED.'///' IN ',A4,'  ...  STOP  ...')             MEMO-016
      END                                                               MEMO-017
