!> \file: main.F90
!> \author: Sayop Kim

PROGRAM main
   USE Parameters_m, ONLY: wp
   USE ReadGrid_m, ONLY: INGRID
   USE MainLoop_m
   USE SimulationSetup_m, ONLY: InitializeGridArrays, SetInitialConditions, &
                                SetTransformedVariables
   USE GridJacobian_m, ONLY: InverseGridMetricsArrays, GridMetricsArrays
   USE io_m

   IMPLICIT NONE
   CHARACTER(LEN=filenameLength) :: outputfile = 'output.tec'
   REAL(KIND=wp) :: START, FINISH

   CALL CPU_TIME(START)
   CALL ReadInput()
   CALL InitializeGridArrays()
   CALL INGRID()
   CALL InverseGridMetricsArrays()
   CALL GridMetricsArrays()
   CALL SetInitialConditions()
   CALL MainLoop()
   CALL WriteTecPlot(outputfile,'"I","J","JACOBIAN","rho","u","v","pres", &
                                 "temp","MACH"')
   CALL CPU_TIME(FINISH)
   WRITE(*,*) "TOTAL CPU TIME:",FINISH-START," SECONDS"

END PROGRAM main
