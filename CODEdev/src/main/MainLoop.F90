!> \file: MainLoop.F90
!> \author: Sayop Kim

MODULE MainLoop_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE

CONTAINS

!-----------------------------------------------------------------------------!
  SUBROUTINE MainLoop()
!-----------------------------------------------------------------------------!
    USE io_m, ONLY: WriteErrorLog, WriteDataOut
    USE TimeIntegration_m
    USE SimulationVars_m, ONLY: nadv, nmax, imin, imax, jmin, jmax, V, U, RMSerr, &
                                ifinish
    USE SimulationSetup_m, ONLY: SetBoundaryConditions, SetTransformedVariables, &
                                 SetPrimativeVariables

    INTEGER :: i,j
    TimeLoop: DO nadv = 1, nmax
      CALL SetBoundaryConditions()
      CALL SetTimeStep()
      CALL TimeIntegration()
      CALL SetPrimativeVariables()
      ! Update Contravariant velocities that need to be updated
      ! in SetBoundaryConditions
      !CALL SetTransformedVariables()
      CALL CheckConvergence()
      IF(IFINISH .EQ. 1 .OR. nadv .EQ. nmax) THEN
        WRITE(*,'(A27,I6,A10,g15.6)') 'NORMAL TERMINATION AT NADV=',nadv, &
                                      'and RMSerr=',RMSerr
        return
      ELSE
        WRITE(*,'(A5,I6,A3,g15.6,A4,g15.6,A8,g15.6)') 'NADV=',nadv, &
                                                     'T=',t, 'DT=', dt, &
                                                     'RMSerr=',RMSerr
        CALL WriteErrorLog(nadv)
      END IF
    END DO TimeLoop
    CALL WriteErrorLog(nadv)
    CALL WriteDataOut()

  END SUBROUTINE MainLoop
END MODULE MainLoop_m

