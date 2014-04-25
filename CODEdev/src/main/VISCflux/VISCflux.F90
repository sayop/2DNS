!> \file: SimulationSetup.F90
!> \author: Sayop Kim

MODULE VISCflux_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE

CONTAINS
!-----------------------------------------------------------------------------!
  SUBROUTINE VISCflux()
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: imin, jmin, imax, jmax, incell, jncell, &
                                DF, V
    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY, JACOBIAN

    INTEGER i, j, n
    REAL(KIND=wp), DIMENSION(4,incell,jncell) :: FPhalf, GPhalf
    REAL(KIND=wp) :: mx, my, uhalf, vhalf, tau_xx, tau_yy, tau_xy, &
                     JACOBIANavg

    FPhalf = 0.0_wp
    GPhalf = 0.0_wp
    !
    ! Set transformed flux along j-const lines at every (i+1/2) points
    !
    DO j = jmin + 1, jmax - 1
      DO i = imin, imax - 1
        JACOBIANavg = 0.5_wp * (JACOBIAN(i+1,j) + JACOBIAN(i,j))
        mx = 0.5_wp * ( PIPX(i+1,j) + PIPX(i,j) )
        my = 0.5_wp * ( PIPY(i+1,j) + PIPY(i,j) )
        uhalf = 0.5_wp * ( V(2,i+1,j) + V(2,i,j) )
        vhalf = 0.5_wp * ( V(3,i+1,j) + V(3,i,j) )
        tau_xx = ShearStress(i,j,1,1,'i')
        tau_yy = ShearStress(i,j,2,2,'i')
        tau_xy = ShearStress(i,j,1,2,'i')
        FPhalf(1,i,j) = 0.0_wp
        FPhalf(2,i,j) = mx * tau_xx + my * tau_xy
        FPhalf(3,i,j) = mx * tau_xy + my * tau_yy
        FPhalf(4,i,j) = mx * ( uhalf * tau_xx + vhalf * tau_xy - &
                               HeatFlux(i,j,1,'i') ) + &
                        my * ( uhalf * tau_xy + vhalf * tau_yy - &
                               HeatFlux(i,j,2,'i') )                       
        DO n = 1, 4
          FPhalf(n,i,j) = FPhalf(n,i,j) / JACOBIANavg
        END DO
      END DO
    END DO
    
    !
    ! Set transformed flux in along i-const lines at every (j+1/2) points
    !
    DO i = imin + 1, imax - 1
      DO j = jmin, jmax - 1
        JACOBIANavg = 0.5_wp * (JACOBIAN(i,j+1) + JACOBIAN(i,j))
        mx = 0.5_wp * ( PJPX(i,j+1) + PJPX(i,j) )
        my = 0.5_wp * ( PJPY(i,j+1) + PJPY(i,j) )
        uhalf = 0.5_wp * ( V(2,i,j+1) + V(2,i,j) )
        vhalf = 0.5_wp * ( V(3,i,j+1) + V(3,i,j) )
        tau_xx = ShearStress(i,j,1,1,'j')
        tau_yy = ShearStress(i,j,2,2,'j')
        tau_xy = ShearStress(i,j,1,2,'j')
        GPhalf(1,i,j) = 0.0_wp
        GPhalf(2,i,j) = mx * tau_xx + my * tau_xy
        GPhalf(3,i,j) = mx * tau_xy + my * tau_yy
        GPhalf(4,i,j) = mx * ( uhalf * tau_xx + vhalf * tau_xy - &
                               HeatFlux(i,j,1,'j') ) + &
                        my * ( uhalf * tau_xy + vhalf * tau_yy - &
                               HeatFlux(i,j,2,'j') )
        DO n = 1, 4
          GPhalf(n,i,j) = GPhalf(n,i,j) / JACOBIANavg
        END DO
      END DO
    END DO

    ! Update DF vector
    DO i = imin + 1, imax - 1
      DO j = jmin + 1, jmax - 1
        DO n = 1, 4
          DF(n,i,j) = DF(n,i,j) - ( &
                      FPhalf(n,i,j) - FPhalf(n,i-1,j) + &
                      GPhalf(n,i,j) - GPhalf(n,i,j-1) )
        END DO
      END DO
    END DO

  END SUBROUTINE VISCflux

  FUNCTION SutherlandLaw(temp,isu) RESULT(SU)
    !Return dimensional viscosity and thermal conductivity
    INTEGER :: i, j, isu
    REAL(KIND=wp) :: SU, C1, C2, C3, C4, temp

    C1 = 1.458E-06_wp
    C2 = 110.4_wp
    C3 = 2.495E-03_wp
    C4 = 194.0_wp 

    IF(isu .EQ. 1) THEN
      ! Set dimensional viscosity
      SU = C1 * temp ** 1.5 / (temp + C2)
    ELSE IF(isu .EQ. 2) THEN
      ! Set dimensional thermal conductivity
      SU = C3 * temp ** 1.5 / (temp + C4)
    END IF
  END FUNCTION SutherlandLaw

  FUNCTION ShearStress(i,j,i1,i2,axis) RESULT(TAU)
    !Calculate ShearStress at the half point (i+1/2,j) with 'i'
    !or (i,j+1/2) with 'j'
    USE SimulationVars_m, ONLY: V, TEMP, RE_REF, MU_REF
    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY

    INTEGER :: i, j, i1, i2, ip, jp
    REAL :: DUDX, DUDY, DVDX, DVDY, &
            DUDI, DUDJ, DVDI, DVDJ, MU, &
            TAU
    CHARACTER(LEN=1) :: axis

    IF(axis .EQ. 'i') THEN
      ip = 1
      jp = 0
      DUDI = V(2,i+1,j) - V(2,i,j)
      DUDJ = 0.25_wp * ( V(2,i+1,j+1) + V(2,i,j+1) - &
                         V(2,i+1,j-1) - V(2,i,j-1) )
      DVDI = V(3,i+1,j) - V(3,i,j)
      DVDJ = 0.25_wp * ( V(3,i+1,j+1) + V(3,i,j+1) - &
                         V(3,i+1,j-1) - V(3,i,j-1) )
    ELSEIF(axis .EQ. 'j') THEN
      ip = 0
      jp = 1
      DUDI = 0.25_wp * ( V(2,i+1,j+1) + V(2,i+1,j) - &
                         V(2,i-1,j+1) - V(2,i-1,j) )
      DUDJ = V(2,i,j+1) - V(2,i,j)
      DVDI = 0.25_wp * ( V(3,i+1,j+1) + V(3,i+1,j) - &
                         V(3,i-1,j+1) - V(3,i-1,j) )
      DVDJ = V(3,i,j+1) - V(3,i,j)
    END IF


    IF(i1 .EQ. i2) THEN
      !Set DUDX
      DUDX = 0.5_wp * ( (PIPX(i+ip,j+jp) + PIPX(i,j)) * DUDI + &
                        (PJPX(i+ip,j+jp) + PJPX(i,j)) * DUDJ )
      !Set DVDY
      DVDY = 0.5_wp * ( (PIPY(i+ip,j+jp) + PIPY(i,j)) * DVDI + &
                        (PJPY(i+ip,j+jp) + PJPY(i,j)) * DVDJ )
      IF(i1 .EQ. 1) THEN
        !Set TAU_xx
        TAU = (2.0_wp * DUDX - DVDY) / 1.5_wp
      ELSE
        !Set TAU_yy
        TAU = (2.0_wp * DVDY - DUDX) / 1.5_wp
      END IF
    ELSE
      !Set DUDY
      DUDY = 0.5_wp * ( (PIPY(i+ip,j+jp) + PIPY(i,j)) * DUDI + &
                        (PJPY(i+ip,j+jp) + PJPY(i,j)) * DUDJ )
      !Set DVDX
      DVDX = 0.5_wp * ( (PIPX(i+ip,j+jp) + PIPX(i,j)) * DVDI + &
                        (PJPX(i+ip,j+jp) + PJPX(i,j)) * DVDJ )
      !Set TAU_xy
      TAU = DUDY + DVDX
    END IF

    !Dimensionless viscosity coefficient
    MU = 0.5_wp * ( SutherlandLaw(TEMP(i+ip,j+jp),1) + &
                    SutherlandLaw(TEMP(i,j),1) ) / MU_REF

    !Final update for TAU
    TAU = MU * TAU / RE_REF 
  END FUNCTION ShearStress

  FUNCTION HeatFlux(i,j,id,axis) RESULT(HF)
    !Calculate HeatFlux at the half point (i+1/2,j) with 'i'
    !or (i,j+1/2) with 'j'
    USE SimulationVars_m, ONLY: V, TEMP, RE_REF, MU_REF, MACH_REF, CP, &
                                cgamma
    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY
    REAL(KIND=wp) :: HF, PR, DTDI, DTDJ, DTDX, DTDY, MU, K
    INTEGER :: i, j, id, ip, jp
    CHARACTER(LEN=1) :: axis

    IF(axis .EQ. 'i') THEN
      ip = 1
      jp = 0
      !DTDI = TEMP(i+1,j) - TEMP(i,j)
      DTDI = V(5,i+1,j) - V(5,i,j)
      !DTDJ = 0.25_wp * ( TEMP(i+1,j+1) + TEMP(i,j+1) - &
      !                   TEMP(i+1,j-1) - TEMP(i,j-1) )
      DTDJ = 0.25_wp * ( V(5,i+1,j+1) + V(5,i,j+1) - &
                         V(5,i+1,j-1) - V(5,i,j-1) )
    ELSE
      ip = 0
      jp = 1
      !DTDJ = TEMP(i,j+1) - TEMP(i,j)
      DTDJ = V(5,i,j+1) - V(5,i,j)
      !DTDI = 0.25_wp * ( TEMP(i+1,j+1) + TEMP(i+1,j) - &
      !                   TEMP(i-1,j-1) - TEMP(i-1,j) )
      DTDI = 0.25_wp * ( V(5,i+1,j+1) + V(5,i+1,j) - &
                         V(5,i-1,j-1) - V(5,i-1,j) )
    END IF

    !Dimensional viscosity
    MU = 0.5_wp * ( SutherlandLaw(TEMP(i+ip,j+jp),1) + &
                    SutherlandLaw(TEMP(i,j),1) )
    !Dimensional thermal conductivity
    K = 0.5_wp * ( SutherlandLaw(TEMP(i+ip,j+jp),2) + &
                   SutherlandLaw(TEMP(i,j),2) ) 


    !Set Prandtl number
    PR = CP * MU / K

    !NonDimensionalize MU
    MU = MU / MU_REF

    IF(id .EQ. 1) THEN
      !Set HeatFlux in x-direction: Q_x
      !Set DTDX
      DTDX = 0.5_wp * ( (PIPX(i+ip,j+jp) + PIPX(i,j)) * DTDI + &
                        (PJPX(i+ip,j+jp) + PJPX(i,j)) * DTDJ )
      HF = -MU * DTDX / ( (cgamma - 1.0_wp) * MACH_REF ** 2 * RE_REF * PR )
    ELSE
      !Set HeatFlux in y-direction: Q_y
      !Set DTDY
      DTDY = 0.5_wp * ( (PIPY(i+ip,j+jp) + PIPY(i,j)) * DTDI + &
                        (PJPY(i+ip,j+jp) + PJPY(i,j)) * DTDJ )
      HF = -MU * DTDY / ( (cgamma - 1.0_wp) * MACH_REF ** 2 * RE_REF * PR )
    END IF
  END FUNCTION HeatFlux
END MODULE VISCflux_m
