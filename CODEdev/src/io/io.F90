!> \file: io.F90
!> \author: Sayop Kim
!> \brief: Provides routines to read input and write output

MODULE io_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE

   ! cgammainit: initial heat capacity ratio: Cp/Cv
   REAL(KIND=wp) :: rhoinit, uinit, vinit, cgammainit, pinit, tinit
   INTEGER, PARAMETER :: IOunit = 10, filenameLength = 64
   CHARACTER(LEN=50) :: prjTitle
   CHARACTER(LEN=filenameLength) :: gridFile
CONTAINS

!-----------------------------------------------------------------------------!
   SUBROUTINE ReadInput()
!-----------------------------------------------------------------------------!
!  Read input files for transformation 1:
!-----------------------------------------------------------------------------!
     USE SimulationVars_m, ONLY: imax, jmax, ngl, nmax, errLimit, &
                                 dens_ref, temp_ref, vel_ref, length_ref, &
                                 cp, ivisc
     USE TimeIntegration_m, ONLY: CFL
     USE AUSMPWplus_m, ONLY: alpha, epsil, kappa, limiter
     USE BCvisc_m, ONLY: wallTemp
     IMPLICIT NONE
     INTEGER :: ios
     CHARACTER(LEN=10) :: inputVar

     OPEN(IOunit, FILE = 'input.dat', FORM = 'FORMATTED', ACTION = 'READ', &
           STATUS = 'OLD', IOSTAT = ios)
     IF(ios /= 0) THEN
        WRITE(*,'(a)') ""
        WRITE(*,'(a)') "Fatal error: Could not open the input data file."
        RETURN
     ELSE
        WRITE(*,'(a)') ""
        WRITE(*,'(a)') "Reading input file for transformation 1"
     ENDIF

     READ(IOunit,*)
     READ(IOunit,'(a)') prjTitle
     WRITE(*,'(4a)') 'Project Title:', '"',TRIM(prjTitle),'"'
     READ(IOunit,*) inputVar, imax
     WRITE(*,'(a,i6)') inputVar,imax
     READ(IOunit,*) inputVar, jmax
     WRITE(*,'(a,i6)') inputVar,jmax
     READ(IOunit,*) inputVar, ngl
     WRITE(*,'(a,i6)') inputVar, ngl
     READ(IOunit,*) inputVar, gridFile
     WRITE(*,'(a,2x,a)') inputVar, gridFile

     !Read initial conditions
     READ(IOunit,*)
     READ(IOunit,*) inputVar, rhoinit
     WRITE(*,'(a,g15.6)') inputVar, rhoinit
     READ(IOunit,*) inputVar, uinit
     WRITE(*,'(a,g15.6)') inputVar, uinit
     READ(IOunit,*) inputVar, vinit
     WRITE(*,'(a,g15.6)') inputVar, vinit
     READ(IOunit,*) inputVar, pinit
     WRITE(*,'(a,g15.6)') inputVar, pinit
     READ(IOunit,*) inputVar, tinit
     WRITE(*,'(a,g15.6)') inputVar, tinit
     READ(IOunit,*) inputVar, cgammainit
     WRITE(*,'(a,g15.6)') inputVar, cgammainit

     !Read free stream conditions
     READ(IOunit,*)
     READ(IOunit,*) inputVar, dens_ref
     WRITE(*,'(a,g15.6)') inputVar, dens_ref
     READ(IOunit,*) inputVar, temp_ref
     WRITE(*,'(a,g15.6)') inputVar, temp_ref
     READ(IOunit,*) inputVar, vel_ref
     WRITE(*,'(a,g15.6)') inputVar, vel_ref
     READ(IOunit,*) inputVar, length_ref
     WRITE(*,'(a,g15.6)') inputVar, length_ref
     READ(IOunit,*) inputVar, CP
     WRITE(*,'(a,g15.6)') inputVar, CP

     !Read Simulation parameters
     READ(IOunit,*)
     READ(IOunit,*) inputVar, nmax
     WRITE(*,'(a,i6)') inputVar, nmax
     READ(IOunit,*) inputVar, CFL
     WRITE(*,'(a,g15.6)') inputVar, CFL
     READ(IOunit,*) inputVar, errLimit
     WRITE(*,'(a,g15.6)') inputVar, errLimit

     !Read AUSMPW+ parameters
     READ(IOunit,*)
     READ(IOunit,*) inputVar, alpha
     WRITE(*,'(a,g15.6)') inputVar, alpha
     READ(IOunit,*) inputVar, epsil
     WRITE(*,'(a,i6)') inputVar, epsil
     READ(IOunit,*) inputVar, limiter
     WRITE(*,'(a,i6)') inputVar, limiter
     READ(IOunit,*) inputVar, kappa
     WRITE(*,'(a,g15.6)') inputVar, kappa

     !Viscous solver parameters
     READ(IOunit,*)
     READ(IOunit,*) inputVar, ivisc
     WRITE(*,'(a,i6)') inputVar, ivisc
     READ(IOunit,*) inputVar, wallTemp
     WRITE(*,'(a,g15.6)') inputVar, wallTemp


  END SUBROUTINE ReadInput

!-----------------------------------------------------------------------------!
   SUBROUTINE WriteTecPlot(fileName,varList)
!-----------------------------------------------------------------------------!
!  Write Tecplot file
!-----------------------------------------------------------------------------!
     USE SimulationVars_m, ONLY: imin, jmin, imax, jmax, ires, jres, xp, &
                                 V, RHO, UVEL, VVEL, PRES, TEMP,  MACH
     USE GridJacobian_m!, ONLY: JACOBIAN
     IMPLICIT NONE
     CHARACTER(LEN=filenameLength), INTENT(IN) :: fileName
     CHARACTER(LEN=*), INTENT(IN) :: varList
     INTEGER :: i, j

     OPEN(IOunit, File = fileName, FORM = 'FORMATTED', ACTION = 'WRITE')
     ! writes the two line TECPLOT header
     WRITE(IOunit,'(a)') 'Title="' // TRIM(prjTitle) // '"'
     WRITE(IOunit,'(a)') 'Variables=' // TRIM(varList)

     WRITE(IOunit,'(a)') ""
     WRITE(IOunit,'(a,i6,a,i6,a)') 'Zone I=', ires, ', J=', jres, ', F=POINT'

     DO j = jmin, jmax
        DO i = imin, imax
!           WRITE(IOunit,'(9g15.6)') xp(1,i,j), xp(2,i,j), JACOBIAN(I,J), &
!                                    V(1,i,j), V(2,i,j), V(3,i,j), V(4,i,j), &
!                                    V(5,i,j), MACH(i,j)
           WRITE(IOunit,'(9g15.6)') xp(1,i,j), xp(2,i,j), JACOBIAN(I,J), &
                                    RHO(i,j), UVEL(i,j), VVEL(i,j), PRES(i,j), &
                                    TEMP(i,j), MACH(i,j)
        ENDDO
     ENDDO
     CLOSE(IOunit)

   END SUBROUTINE WriteTecPlot   

!-----------------------------------------------------------------------------!
  SUBROUTINE WriteErrorLog(iter)
!-----------------------------------------------------------------------------!
!  Write Error Log file
!-----------------------------------------------------------------------------!
  USE SimulationVars_m, ONLY: RMSerr
  IMPLICIT NONE
  INTEGER :: iter
  CHARACTER(LEN=filenameLength) :: fileName = 'ErrorLog.dat'

  IF(iter .EQ. 1) THEN
    OPEN(IOunit, File = fileName, FORM = 'FORMATTED', ACTION = 'WRITE')
    WRITE(IOunit,'(A,2A10)') '#','Iteration','RMS error'
  ELSE
    OPEN(IOunit, File = fileName, FORM = 'FORMATTED', ACTION = 'WRITE', &
         POSITION = 'APPEND')
  END IF

    WRITE(IOunit,'(i6,g15.6)') iter, RMSerr
    CLOSE(IOunit)
  END SUBROUTINE WriteErrorLog

!-----------------------------------------------------------------------------!
  SUBROUTINE WriteDataOut()
!-----------------------------------------------------------------------------!
!  Write primative data out
!-----------------------------------------------------------------------------!
  USE SimulationVars_m, ONLY: V, imin, imax, jmin, jmax, xp, MACH

  CHARACTER(LEN=filenameLength) :: fileName = 'WallData.dat'
  INTEGER :: i, j

  OPEN(IOunit, File = fileName, FORM = 'FORMATTED', ACTION = 'WRITE')
  WRITE(IOunit,'(A,6A15)') '#','x','Density','u','v','p','Mach'
  j = jmin
  DO i = imin, imax
    WRITE(IOunit,'(6g15.6)') XP(1,i,j), V(1,i,j), V(2,i,j), V(3,i,j), V(4,i,j), &
                             MACH(i,j)
  END DO
  END SUBROUTINE WriteDataOut
END MODULE io_m
