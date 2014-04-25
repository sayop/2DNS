!> \file: SimulationSetup.F90
!> \author: Sayop Kim

MODULE RestartDataOut_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE
   INTEGER, PARAMETER :: IOunit = 10
   INTEGER :: ios
CONTAINS
!-----------------------------------------------------------------------------!
  SUBROUTINE IORestartID(id)
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: NADV, INADV
    USE TimeIntegration_m, ONLY: T
    !Input/Output Restart ID
    !id = 1: read RESTART.DATA
    !id = 2: write RESTART.DATA
    INTEGER :: id
    CHARACTER(LEN=12) :: fileName
    CHARACTER(LEN=10) :: inputVar

    fileName = 'RESTART.DATA'

    IF(id .EQ. 1) THEN
      OPEN(IOunit, FILE = fileName, FORM = 'FORMATTED', ACTION = 'READ', &
           STATUS = 'OLD', IOSTAT = ios)
      READ(IOunit,'(A6,I8)') inputVar,inadv
      READ(IOunit,'(A6,g15.8)') inputVar,T
      CLOSE(IOunit)
    ELSE
      OPEN(IOunit, File = fileName, FORM = 'FORMATTED', ACTION = 'WRITE')
      WRITE(IOunit,'(A6,I8)') 'NADV=',nadv
      WRITE(IOunit,'(A6,g15.8)') 'TIME=',T
      CLOSE(IOunit)
    END IF

  END SUBROUTINE IORestartID

!-----------------------------------------------------------------------------!
  SUBROUTINE ReadRestartData()
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: imin, jmin, imax, jmax, V
    USE SimulationSetup_m, ONLY: SetTransformedVariables, SetPrimativeVariables
    REAL(KIND=wp) :: dummy
    INTEGER :: i, j

    OPEN(IOunit, FILE = 'output.tec', FORM = 'FORMATTED', ACTION = 'READ', &
         STATUS = 'OLD', IOSTAT = ios)
    READ(IOunit,*)
    READ(IOunit,*)
    READ(IOunit,*)
    READ(IOunit,*)

    DO j = jmin, jmax
      DO i = imin, imax
        READ(IOunit,'(9g17.8)') dummy, dummy, dummy, V(1,i,j), V(2,i,j), &
                                V(3,i,j), V(4,i,j), V(5,i,j), dummy
      END DO
    END DO
    CLOSE(IOunit)
    CALL SetTransformedVariables()
    CALL SetPrimativeVariables()
  END SUBROUTINE ReadRestartData
END MODULE RestartDataOut_m
