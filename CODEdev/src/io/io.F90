!> \file: io.F90
!> \author: Sayop Kim
!> \brief: Provides routines to read input and write output

MODULE io_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE

   ! cgammainit: initial heat capacity ratio: Cp/Cv
   REAL(KIND=wp) :: rhoinit, uinit, vinit, cgammainit, pinit, tinit
   REAL(KIND=wp) :: xLoc
   INTEGER, PARAMETER :: IOunit = 10, filenameLength = 64
   CHARACTER(LEN=50) :: prjTitle
   CHARACTER(LEN=filenameLength) :: gridFile
CONTAINS

!-----------------------------------------------------------------------------!
   SUBROUTINE ReadInput()
!-----------------------------------------------------------------------------!
!  Read input files for transformation 1:
!-----------------------------------------------------------------------------!
     USE SimulationVars_m, ONLY: imax, jmax, ngl, nmax, restart, errLimit, &
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
     READ(IOunit,*) inputVar, restart
     WRITE(*,'(a,i6)') inputVar, restart
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

     !Read BL profile data collecting point
     READ(IOunit,*)
     READ(IOunit,*) inputVar, xLoc
     WRITE(*,'(a,g15.6)') inputVar, xLoc

     CLOSE(IOunit)
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
           WRITE(IOunit,'(9g17.8)') xp(1,i,j), xp(2,i,j), JACOBIAN(I,J), &
                                    V(1,i,j), V(2,i,j), V(3,i,j), V(4,i,j), &
                                    V(5,i,j), MACH(i,j)
!           WRITE(IOunit,'(9e17.10)') xp(1,i,j), xp(2,i,j), JACOBIAN(I,J), &
!                                    RHO(i,j), UVEL(i,j), VVEL(i,j), PRES(i,j), &
!                                    TEMP(i,j), MACH(i,j)
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

    WRITE(IOunit,'(i6,e15.8)') iter, RMSerr
    CLOSE(IOunit)
  END SUBROUTINE WriteErrorLog

!-----------------------------------------------------------------------------!
  SUBROUTINE WriteDataOut()
!-----------------------------------------------------------------------------!
!  Write primative data out
!-----------------------------------------------------------------------------!
  USE SimulationVars_m, ONLY: V, imin, imax, jmin, jmax, xp, MACH, Temp

  CHARACTER(LEN=filenameLength) :: fileName1 = 'WallData.dat'
  CHARACTER(LEN=filenameLength) :: fileName2 = 'BoundaryLayer.dat'
  INTEGER :: i, j, iLeft, iRight
  REAL(KIND=wp) :: distL, distR, NDUvel, DTemp

  OPEN(IOunit, File = fileName1, FORM = 'FORMATTED', ACTION = 'WRITE')
  WRITE(IOunit,'(A,6A15)') '#','x','Density','u','v','p','Mach'
  j = jmin
  DO i = imin, imax
    WRITE(IOunit,'(6e17.8)') XP(1,i,j), V(1,i,j), V(2,i,j), V(3,i,j), V(4,i,j), &
                             MACH(i,j)
  END DO
  CLOSE(IOunit)

  OPEN(IOunit, File = fileName2, FORM = 'FORMATTED', ACTION = 'WRITE')
  WRITE(IOunit,'(A,A6,g15.6)') '#','xLoc=', xLoc
  WRITE(IOunit,'(A,3A15)') '#','y','NonDim_Uvel','Dim_Temp'
  !Find nearest i=const line to xLoc
  DO j = jmin, jmax
    DO i = imin, imax - 1
      IF(XP(1,i,j) .LE. xLoc .AND. XP(1,i+1,j) .GT. xLoc) THEN
        iLeft = i
        iRight = i + 1
        distL = xLoc - XP(1,iLeft,j)
        distR = XP(1,iRight,j) - xLoc
        !Linear interpolation
        NDUvel = V(2,iLeft,j) + distL / (distL + distR) * ( &
                                V(2,iRight,j) - V(2,iLeft,j) )
        DTemp = Temp(iLeft,j) + distL / (distL + distR) * ( &
                                Temp(iRight,j) - Temp(iLeft,j) )
        EXIT
      END IF
    END DO
    WRITE(IOunit,'(3e17.8)') XP(2,i,j), NDUvel, DTemp
  END DO
  CLOSE(IOunit)
  END SUBROUTINE WriteDataOut
END MODULE io_m
