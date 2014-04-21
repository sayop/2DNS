!> \file: AUSMPWplus.F90
!> \author: Sayop Kim

MODULE AUSMPWplus_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE
   REAL(KIND=wp) :: alpha, epsil, kappa
   INTEGER :: limiter

CONTAINS
!-----------------------------------------------------------------------------!
  SUBROUTINE SetAUSMPWplus()
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: imin, jmin, imax, jmax, incell, jncell, &
                                U, V, DF, cgamma, nadv
    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY, JACOBIAN

    INTEGER i, j, n
    ! FPpHalf: Transformed flux in i-direction at (i+1/2,j)
    ! FPmHalf: Transformed flux in i-direction at (i-1/2,j)
    ! GPpHalf: Transformed flux in j-direction at (i,j+1/2)
    ! GPmHalf: Transformed flux in j-direction at (i,j-1/2)
    !REAL(KIND=wp), DIMENSION(4) :: FPpHalf, FPmHalf, GPpHalf, GPmHalf
    ! UL: Left-extrapolated state vector at (i+1/2,j) or (i,j+1/2)
    ! UR: Right-extrapolated state vector at (i+1/2,j) or (i,j+1/2)
    !REAL(KIND=wp), DIMENSION(4) :: UL, UR

    ! FPhalf: Transformed flux in i-direction at i+1/2
    ! GPhalf: Transformed flux in j-direction at j+1/2
    REAL(KIND=wp), DIMENSION(4,incell,jncell) :: FPhalf, GPhalf
    REAL(KIND=wp) :: XIX, XIY, ETAX, ETAY, JACOBIANavg, A1, A2, &
                     rr, uu, vv, Pl, Pr, h0L, h0R, h0norm, &
                     rrl, rrr, uul, uur, vvl, vvr, &
                     utildL, vtildL, utildR, vtildR, Cavg, &
                     MtildL, MtildR, MtildLplus, MtildRminus, &
                     Pplus, Pminus, Pmin, fl, fr, Mavg
                     
    REAL(KIND=wp), DIMENSION(4) :: UL, UR, VL, VR, XEPl, XEPr
    REAL(KIND=wp), DIMENSION(2) :: UCL, UCR

    !-----------------------------------
    !  
    !  o-------o-------o-------o-------o
    !  |       |       |       |       |
    !  |       G       G       G       |
    !  |       |       |       |       |
    !  o---F---o---F---o---F---o---F---o
    !  |       |       |       |       |
    !  |       G       G       G       |
    !  |       |       |       |       |
    !  o---F---o---F---o---F---o---F---o
    !  |       |       |       |       |
    !  |       G       G       G       |
    !  |       |       |       |       |
    !  o-------o-------o-------o-------o
    !
    ! o: Grid points
    ! F: node for FP flux variables
    ! G: node for GP flux variables

    DF = 0.0_wp
    FPhalf = 0.0_wp
    GPhalf = 0.0_wp
    XEPl = 0.0_wp
    XEPr = 0.0_wp
    !
    ! Set transformed flux in along j-const lines at every (i+1/2) points
    !
    DO j = jmin + 1, jmax - 1
      DO i = imin, imax - 1
        ! Arithmatic average of grid metrics at i + 1/2
        XIX = 0.5_wp * (PIPX(i,j) + PIPX(i+1,j))
        XIY = 0.5_wp * (PIPY(i,j) + PIPY(i+1,j))
        ETAX = 0.5_wp * (PJPX(i,j) + PJPX(i+1,j))
        ETAY = 0.5_wp * (PJPY(i,j) + PJPY(i+1,j))
        JACOBIANavg = 0.5_wp * (JACOBIAN(i,j) + JACOBIAN(i+1,j))
        A1 = sqrt(XIX ** 2 + XIY ** 2)
        A2 = sqrt(ETAX ** 2 + ETAY ** 2)
        ! Set u- and v-velcotiy at i + 1/2
        ! First calculate left-extrapolated state vector and convert
        UL = LeftExtrapolateU('i',i,j)
        rr = UL(1)
        uu = UL(2) / rr
        vv = UL(3) / rr
        Pl = (cgamma - 1.0_wp) * ( UL(4) - &
             0.5_wp * rr * (uu ** 2 + vv ** 2) )
!        VL = LeftExtrapolateV('i',i,j)
!        rr = VL(1)
!        uu = VL(2)
!        vv = VL(3)
!        Pl = VL(4)
        h0L = Pl * cgamma / (rr * (cgamma - 1.0_wp)) + &
              0.5_wp * (uu ** 2 + vv ** 2)
        !UCL = LeftExtrapolateUC('i',i,j)
        !utildL = UCL(1) / A1
        !vtildL = UCL(2) / A1
        utildL = (XIX * uu + XIY * vv) / A1
        !vtildL = (ETAX * uu + ETAY * vv) / A1
        vtildL = (ETAX * uu + ETAY * vv) / A2
        !utildL = XIX * uu / A1
        ! Then calculate right-extrapolated state vector from i+1 point
        UR = RightExtrapolateU('i',i,j)
        rr = UR(1)
        uu = UR(2) / rr
        vv = UR(3) / rr
        Pr = (cgamma - 1.0_wp) * ( UR(4) - &
             0.5_wp * rr * (uu ** 2 + vv ** 2) )
!        VR = RightExtrapolateV('j',i,j)
!        rr = VR(1)
!        uu = VR(2)
!        vv = VR(3)
!        Pr = VR(4)
        h0R = Pr * cgamma / (rr * (cgamma - 1.0_wp)) + &
              0.5_wp * (uu ** 2 + vv ** 2)
        !UCR = RightExtrapolateUC('i',i,j)
        !utildR = UCR(1) / A1
        !vtildR = UCR(2) / A1
        utildR = (XIX * uu + XIY * vv) / A1
        !vtildR = (ETAX * uu + ETAY * vv) / A1
        vtildR = (ETAX * uu + ETAY * vv) / A2
        !utildR = XIX * uu / A1
        ! Calculate stagnation enthalpy normal to the interface
        h0norm = 0.5_wp * (h0L + h0R - 0.5_wp * (vtildL ** 2 + vtildR ** 2))
        ! Calculate cell averaged speed of sound
        ! Use Cavg for calculating Cs and then update to Cavg
        ! Use rr temporarily to store half of sum of utildL and utildR
        Cavg = sqrt(2.0_wp * h0norm * (cgamma - 1.0_wp) / (cgamma + 1.0_wp))
        rr = 0.5_wp * (utildL + utildR)
        IF(rr .GE. 0.0) THEN
          Cavg = Cavg ** 2 / max(abs(utildL),Cavg)
        ELSE
          Cavg = Cavg ** 2 / max(abs(utildR),Cavg)
        END IF
        ! Calculate cell face Mach numbers
        MtildL = utildL / Cavg
        MtildR = utildR / Cavg
        ! Calculate split Mach numbers and pressures
        IF(abs(MtildL) .LE. 1.0) THEN
          MtildLplus = 0.25_wp * (MtildL + 1.0_wp) ** 2
          Pplus = MtildLplus * (2.0_wp - MtildL) + alpha * MtildL * &
                  (MtildL ** 2 - 1.0_wp) ** 2
        ELSE
          MtildLplus = 0.5_wp * (MtildL + abs(MtildL))
          Pplus = 0.5_wp * (1.0_wp + sign(1.0_wp,MtildL))
        END IF
        IF(abs(MtildR) .LE. 1.0) THEN
          MtildRminus = -0.25_wp * (MtildR - 1.0_wp) ** 2
          Pminus = -MtildRminus * (2.0_wp + MtildR) - alpha * MtildR * &
                   (MtildR ** 2 - 1.0_wp) ** 2
        ELSE
          MtildRminus = 0.5_wp * (MtildR - abs(MtildR))
          Pminus = 0.5_wp * (1.0_wp - sign(1.0_wp,MtildR))
        END IF
        ! Set pressure weighting terms
        ! Use rr temporarily to store Ps and w
        rr = Pplus * Pl + Pminus * Pr
        Pmin = min(V(4,i,j-1), V(4,i,j+1), V(4,i+1,j-1), V(4,i+1,j+1))
        IF(rr .LE. 0.0) THEN
          fl = 0.0_wp
          fr = 0.0_wp
        ELSE
          fl = (Pl / rr - 1.0_wp) * min(1.0_wp,Pmin/min(Pl,Pr)) ** 2
          fr = (Pr / rr - 1.0_wp) * min(1.0_wp,Pmin/min(Pl,Pr)) ** 2
        END IF
        rr = 1.0_wp - min(Pl/Pr, Pr/Pl) ** 3
        ! Apply the weightings to MtildLplus, MtildRminus
        ! Now MtildLplus and MtildRminus become averaged Mach number
        ! on the cell face
        Mavg = MtildLplus + MtildRminus
        IF(Mavg .GE. 0.0) THEN
          MtildLplus = MtildLplus + MtildRminus * ((1.0_wp - rr) * &
                       (1.0_wp + fr) - fl)
          MtildRminus = MtildRminus * rr * (1.0_wp + fr)
        ELSE
          MtildRminus = MtildRminus + MtildLplus * ((1.0_wp - rr) * &
                        (1.0_wp + fl) - fr)
          MtildLplus = MtildLplus * rr * (1.0_wp + fl)
        END IF
        ! Update FPhalf vector
        ! Use rr as a temporary variable to store coefficients
        !XEPl = Pl * LeftExtrapolateXE('i',i,j)
        !XEPr = Pr * RightExtrapolateXE('i',i,j)
        XEPl = 0.0_wp
        XEPl(2) = Pl * XIX
        XEPl(3) = Pl * XIY
        XEPr = 0.0_wp
        XEPr(2) = Pr * XIX
        XEPr(3) = Pr * XIY

        DO n = 1, 4
          FPhalf(n,i,j) = ( MtildLplus * Cavg * A1 * UL(n) + &
                            MtildRminus * Cavg * A1 * UR(n) + &
                            Pplus * XEPl(n) + Pminus * XEPr(n) ) / &
                            JACOBIANavg
        END DO
      END DO
    END DO

    XEPl = 0.0_wp
    XEPr = 0.0_wp
    !
    ! Set transformed flux in along i-const lines at every (j+1/2) points
    !
    DO i = imin + 1, imax - 1
      DO j = jmin, jmax - 1
        ! Arithmatic average of grid metrics at i + 1/2
        XIX = 0.5_wp * (PIPX(i,j) + PIPX(i,j+1))
        XIY = 0.5_wp * (PIPY(i,j) + PIPY(i,j+1))
        ETAX = 0.5_wp * (PJPX(i,j) + PJPX(i,j+1))
        ETAY = 0.5_wp * (PJPY(i,j) + PJPY(i,j+1))
        JACOBIANavg = 0.5_wp * (JACOBIAN(i,j) + JACOBIAN(i,j+1))
        A1 = sqrt(XIX ** 2 + XIY ** 2)
        A2 = sqrt(ETAX ** 2 + ETAY ** 2)
        ! Set u- and v-velcotiy at j + 1/2
        ! First calculate left-extrapolated state vector and convert
        UL = LeftExtrapolateU('j',i,j)
        rr = UL(1)
        uu = UL(2) / rr
        vv = UL(3) / rr
        Pl = (cgamma - 1.0_wp) * ( UL(4) - &
             0.5_wp * rr * (uu ** 2 + vv ** 2) )
!        VL = LeftExtrapolateV('j',i,j)
!        rr = VL(1)
!        uu = VL(2)
!        vv = VL(3)
!        Pl = VL(4)
        h0L = Pl * cgamma / (rr * (cgamma - 1.0_wp)) + &
              0.5_wp * (uu ** 2 + vv ** 2)
        !UCL = LeftExtrapolateUC('j',i,j)
        !utildL = UCL(1) / A1
        !vtildL = UCL(2) / A1
        !utildL = (XIX * uu + XIY * vv) / A1
        utildL = (XIX * uu + XIY * vv) / A1
        vtildL = (ETAX * uu + ETAY * vv) / A2
        !utildL = XIX * uu / A1
        !vtildL = ETAY * vv / A1
        ! Then calculate right-extrapolated state vector from i+1 point
        UR = RightExtrapolateU('j',i,j)
        rr = UR(1)
        uu = UR(2) / rr
        vv = UR(3) / rr
        Pr = (cgamma - 1.0_wp) * ( UR(4) - &
             0.5_wp * rr * (uu ** 2 + vv ** 2) )
!        VR = RightExtrapolateV('j',i,j)
!        rr = VR(1)
!        uu = VR(2)
!        vv = VR(3)
!        Pr = VR(4)
        h0R = Pr * cgamma / (rr * (cgamma - 1.0_wp)) + &
              0.5_wp * (uu ** 2 + vv ** 2)
        !UCR = RightExtrapolateUC('j',i,j)
        !utildR = UCR(1) / A1
        !vtildR = UCR(2) / A1
        !utildR = (XIX * uu + XIY * vv) / A1
        utildR = (XIX * uu + XIY * vv) / A1
        vtildR = (ETAX * uu + ETAY * vv) / A2
        !utildR = XIX * uu / A1
        !vtildR = ETAY * vv / A1
        ! Calculate stagnation enthalpy normal to the interface
        h0norm = 0.5_wp * (h0L + h0R - 0.5_wp * (utildL ** 2 + utildR ** 2))
        ! Calculate cell averaged speed of sound
        ! Use Cavg for calculating Cs and then update to Cavg
        ! Use rr temporarily to store half of sum of vtildL and vtildR
        Cavg = sqrt(2.0_wp * h0norm * (cgamma - 1.0_wp) / (cgamma + 1.0_wp))
        rr = 0.5_wp * (vtildL + vtildR)
        IF(rr .GE. 0.0) THEN
          Cavg = Cavg ** 2 / max(abs(vtildL),Cavg)
        ELSE
          Cavg = Cavg ** 2 / max(abs(vtildR),Cavg)
        END IF
        ! Calculate cell face Mach numbers
        MtildL = vtildL / Cavg
        MtildR = vtildR / Cavg
        ! Calculate split Mach numbers and pressures
        IF(abs(MtildL) .LE. 1.0) THEN
          MtildLplus = 0.25_wp * (MtildL + 1.0_wp) ** 2
          Pplus = MtildLplus * (2.0_wp - MtildL) + alpha * MtildL * &
                  (MtildL ** 2 - 1.0_wp) ** 2
        ELSE
          MtildLplus = 0.5_wp * (MtildL + abs(MtildL))
          Pplus = 0.5_wp * (1.0_wp + sign(1.0_wp,MtildL))
        END IF
        IF(abs(MtildR) .LE. 1.0) THEN
          MtildRminus = -0.25_wp * (MtildR - 1.0_wp) ** 2
          Pminus = -MtildRminus * (2.0_wp + MtildR) - alpha * MtildR * &
                   (MtildR ** 2 - 1.0_wp) ** 2
        ELSE
          MtildRminus = 0.5_wp * (MtildR - abs(MtildR))
          Pminus = 0.5_wp * (1.0_wp - sign(1.0_wp,MtildR))
        END IF
        ! Set pressure weighting terms
        ! Use rr temporarily to store Ps and w
        rr = Pplus * Pl + Pminus * Pr
        Pmin = min(V(4,i-1,j), V(4,i+1,j), V(4,i-1,j+1), V(4,i+1,j+1))
        IF(rr .LE. 0.0) THEN
          fl = 0.0_wp
          fr = 0.0_wp
        ELSE
          fl = (Pl / rr - 1.0_wp) * min(1.0_wp,Pmin/min(Pl,Pr)) ** 2
          fr = (Pr / rr - 1.0_wp) * min(1.0_wp,Pmin/min(Pl,Pr)) ** 2
        END IF
        rr = 1.0_wp - min(Pl/Pr, Pr/Pl) ** 3
        ! Apply the weightings to MtildLplus, MtildRminus
        ! Now MtildLplus and MtildRminus become averaged Mach number
        ! on the cell face
        Mavg = MtildLplus + MtildRminus
        IF(Mavg .GE. 0.0) THEN
          MtildLplus = MtildLplus + MtildRminus * ((1.0_wp - rr) * &
                       (1.0_wp + fr) - fl)
          MtildRminus = MtildRminus * rr * (1.0_wp + fr)
        ELSE
          MtildRminus = MtildRminus + MtildLplus * ((1.0_wp - rr) * &
                        (1.0_wp + fl) - fr)
          MtildLplus = MtildLplus * rr * (1.0_wp + fl)
        END IF
        ! Update GPhalf vector
        ! Use rr as a temporary variable to store coefficients
        !XEPl = Pl * LeftExtrapolateXE('j',i,j)
        !XEPr = Pr * RightExtrapolateXE('j',i,j)
        XEPl = 0.0_wp
        XEPl(2) = Pl * ETAX
        XEPl(3) = Pl * ETAY
        XEPr = 0.0_wp
        XEPr(2) = Pr * ETAX
        XEPr(3) = Pr * ETAY
        DO n = 1, 4
          GPhalf(n,i,j) = ( MtildLplus * Cavg * A2 * UL(n) + &
                            MtildRminus * Cavg * A2 * UR(n) + &
                            Pplus * XEPl(n) + Pminus * XEPr(n) ) / &
                            JACOBIANavg
        END DO
      END DO
    END DO

    ! Update DF vector
    DO i = imin + 1, imax - 1
      DO j = jmin + 1, jmax - 1
        DO n = 1, 4
          DF(n,i,j) = FPhalf(n,i,j) - FPhalf(n,i-1,j) + &
                      GPhalf(n,i,j) - GPhalf(n,i,j-1)
        END DO
      END DO
    END DO

  END SUBROUTINE SetAUSMPWplus

  FUNCTION LeftExtrapolateU(axis,i,j) RESULT(UL)
    USE SimulationVars_m, ONLY: U, imin, jmin
    IMPLICIT NONE
    INTEGER :: i, j, n, im, jm, ip, jp
    CHARACTER(LEN=1) :: axis
    REAL(KIND=wp), DIMENSION(4) :: UL
    ! DelUmHalf: Delta_U_i-1/2 = U_i - U_i-1
    ! DelUpHalf: Delta_U_i+1/2 = U_i+1 - U_i
    ! phi1 = phi(rL)
    ! phi2 = phi(1/rL)
    REAL(KIND=wp) :: DelUmHalf, DelUpHalf, rL, phi1, phi2, beta

    IF(axis .EQ. 'i') THEN
      im = -1
      jm = 0
      ip = 1
      jp = 0
    ELSEIF(axis .EQ. 'j') THEN
      im = 0
      jm = -1
      ip = 0
      jp = 1
    END IF
    beta = (3.0_wp - kappa) / (1.0_wp - kappa)
    DO n = 1, 4
      IF(i .EQ. imin .OR. j .EQ. jmin) THEN
        UL(n) = U(n,i,j)
        CYCLE
      END IF
      DelUmHalf = U(n,i,j) - U(n,i+im,j+im)
      DelUpHalf = U(n,i+ip,j+jp) - U(n,i,j)
      SELECT CASE(limiter)
      CASE(0)
        phi1 = DelUmHalf
        phi2 = DelUpHalf
      CASE(1)
        !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
        phi1 = minmod(DelUmHalf, beta*DelUpHalf)
        phi2 = minmod(DelUpHalf, beta*DelUmHalf)
      END SELECT
      UL(n) = U(n,i,j) + 0.25_wp * epsil * (&
              (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
    END DO
  END FUNCTION LeftExtrapolateU

!  FUNCTION LeftExtrapolateV(axis,i,j) RESULT(VL)
!    USE SimulationVars_m, ONLY: V, imin, jmin
!    IMPLICIT NONE
!    INTEGER :: i, j, n, im, jm, ip, jp
!    CHARACTER(LEN=1) :: axis
!    REAL(KIND=wp), DIMENSION(4) :: VL
!    ! DelUmHalf: Delta_U_i-1/2 = U_i - U_i-1
!    ! DelUpHalf: Delta_U_i+1/2 = U_i+1 - U_i
!    ! phi1 = phi(rL)
!    ! phi2 = phi(1/rL)
!    REAL(KIND=wp) :: DelVmHalf, DelVpHalf, rL, phi1, phi2, beta
!
!    IF(axis .EQ. 'i') THEN
!      im = -1
!      jm = 0
!      ip = 1
!      jp = 0
!    ELSEIF(axis .EQ. 'j') THEN
!      im = 0
!      jm = -1
!      ip = 0
!      jp = 1
!    END IF
!    beta = (3.0_wp - kappa) / (1.0_wp - kappa)
!    DO n = 1, 4
!      IF(i .EQ. imin .OR. j .EQ. jmin) THEN
!        VL(n) = V(n,i,j)
!        CYCLE
!      END IF
!      DelVmHalf = V(n,i,j) - V(n,i+im,j+im)
!      DelVpHalf = V(n,i+ip,j+jp) - V(n,i,j)
!      SELECT CASE(limiter)
!      CASE(0)
!        phi1 = DelVmHalf
!        phi2 = DelvpHalf
!      CASE(1)
!        !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
!        phi1 = minmod(DelVmHalf, beta*DelVpHalf)
!        phi2 = minmod(DelvpHalf, beta*DelVmHalf)
!      END SELECT
!      VL(n) = V(n,i,j) + 0.25_wp * epsil * (&
!              (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
!    END DO
!  END FUNCTION LeftExtrapolateV

!  FUNCTION LeftExtrapolateUC(axis,i,j) RESULT(UCL)
!    USE SimulationVars_m, ONLY: UCON, imin, jmin
!    IMPLICIT NONE
!    INTEGER :: i, j, n, im, jm, ip, jp
!    CHARACTER(LEN=1) :: axis
!    REAL(KIND=wp), DIMENSION(2) :: UCL
!    REAL(KIND=wp) :: DelUmHalf, DelUpHalf, rL, phi1, phi2, beta
!
!    IF(axis .EQ. 'i') THEN
!      im = -1
!      jm = 0
!      ip = 1
!      jp = 0
!    ELSEIF(axis .EQ. 'j') THEN
!      im = 0
!      jm = -1
!      ip = 0
!      jp = 1
!    END IF
!    beta = (3.0_wp - kappa) / (1.0_wp - kappa)
!    DO n = 1, 2
!      IF(i .EQ. imin .OR. j .EQ. jmin) THEN
!        UCL(n) = UCON(n,i,j)
!        CYCLE
!      END IF
!      DelUmHalf = UCON(n,i,j) - UCON(n,i+im,j+im)
!      DelUpHalf = UCON(n,i+ip,j+jp) - UCON(n,i,j)
!      SELECT CASE(limiter)
!      CASE(0)
!        phi1 = DelUmHalf
!        phi2 = DelUpHalf
!      CASE(1)
!        !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
!        phi1 = minmod(DelUmHalf, beta*DelUpHalf)
!        phi2 = minmod(DelUpHalf, beta*DelUmHalf)
!      END SELECT
!      UCL(n) = UCON(n,i,j) + 0.25_wp * epsil * (&
!               (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
!    END DO
!  END FUNCTION LeftExtrapolateUC

!  FUNCTION LeftExtrapolateXE(axis,i,j) RESULT(XE)
!    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY
!    USE SimulationVars_m, ONLY: imin, jmin
!    IMPLICIT NONE
!    INTEGER :: i, j, n, im, jm, ip, jp
!    CHARACTER(LEN=1) :: axis
!    ! XE(1) = X(4) = 0.0
!    ! XE(2): XI_x in i-direction or ETA_x in j-direction
!    ! XE(3): XI_y in i-direction or ETA_y in j-direction
!    REAL(KIND=wp), DIMENSION(4) :: XE
!    REAL(KIND=wp) :: DelUmHalf, DelUpHalf, rL, phi1, phi2, beta
!
!    XE = 0.0_wp
!    beta = (3.0_wp - kappa) / (1.0_wp - kappa)
!    IF(axis .EQ. 'i') THEN
!      IF(i .EQ. imin) THEN
!        ! left boundary
!        XE(2) = PIPX(i,j)
!        XE(3) = PIPY(i,j)
!      ELSE
!        im = -1
!        jm = 0
!        ip = 1
!        jp = 0
!        ! Set XI_X
!        DelUmHalf = PIPX(i,j) - PIPX(i+im,j+im)
!        DelUpHalf = PIPX(i+ip,j+jp) - PIPX(i,j)
!        SELECT CASE(limiter)
!        CASE(0)
!          phi1 = DelUmHalf
!          phi2 = DelUpHalf
!        CASE(1)
!          !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
!          phi1 = minmod(DelUmHalf, beta*DelUpHalf)
!          phi2 = minmod(DelUpHalf, beta*DelUmHalf)
!        END SELECT
!        XE(2) = PIPX(i,j) + 0.25_wp * epsil * (&
!                (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
!        ! Set XI_Y
!        DelUmHalf = PIPY(i,j) - PIPY(i+im,j+im)
!        DelUpHalf = PIPY(i+ip,j+jp) - PIPY(i,j)
!        SELECT CASE(limiter)
!        CASE(0) 
!          phi1 = DelUmHalf
!          phi2 = DelUpHalf
!        CASE(1)
!          !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
!          phi1 = minmod(DelUmHalf, beta*DelUpHalf)
!          phi2 = minmod(DelUpHalf, beta*DelUmHalf)
!        END SELECT
!        XE(3) = PIPY(i,j) + 0.25_wp * epsil * (&
!                (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
!      END IF
!    ELSEIF(axis .EQ. 'j') THEN
!      IF(j .EQ. jmin) THEN
!        ! left boundary
!        XE(2) = PJPX(i,j)
!        XE(3) = PJPY(i,j)
!      ELSE
!        im = 0
!        jm = -1
!        ip = 0
!        jp = 1
!        ! Set ETA_X
!        DelUmHalf = PJPX(i,j) - PJPX(i+im,j+im)
!        DelUpHalf = PJPX(i+ip,j+jp) - PJPX(i,j)
!        SELECT CASE(limiter)
!        CASE(0) 
!          phi1 = DelUmHalf
!          phi2 = DelUpHalf
!        CASE(1)
!          !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
!          phi1 = minmod(DelUmHalf, beta*DelUpHalf)
!          phi2 = minmod(DelUpHalf, beta*DelUmHalf)
!        END SELECT
!        XE(2) = PJPX(i,j) + 0.25_wp * epsil * (&
!                (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
!        ! Set ETA_Y
!        DelUmHalf = PJPY(i,j) - PJPY(i+im,j+im)
!        DelUpHalf = PJPY(i+ip,j+jp) - PJPY(i,j)
!        SELECT CASE(limiter)
!        CASE(0) 
!          phi1 = DelUmHalf
!          phi2 = DelUpHalf
!        CASE(1)
!          !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
!          phi1 = minmod(DelUmHalf, beta*DelUpHalf)
!          phi2 = minmod(DelUpHalf, beta*DelUmHalf)
!        END SELECT
!        XE(3) = PJPY(i,j) + 0.25_wp * epsil * (&
!                (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
!      END IF
!    END IF
!  END FUNCTION LeftExtrapolateXE

  FUNCTION RightExtrapolateU(axis,i,j) RESULT(UR)
    USE SimulationVars_m, ONLY: U, imax, jmax
    IMPLICIT NONE
    INTEGER :: i, j, n, ip, jp, ipp, jpp
    CHARACTER(LEN=1) :: axis
    REAL(KIND=wp), DIMENSION(4) :: UR
    REAL(KIND=wp) :: DelUpHalf, DelUppHalf, rR, phi1, phi2, beta

    IF(axis .EQ. 'i') THEN
      ip = 1
      ipp = 2
      jp = 0
      jpp = 0
    ELSEIF(axis .EQ. 'j') THEN
      ip = 0
      ipp = 0
      jp = 1
      jpp = 2
    END IF
    beta = (3.0_wp - kappa) / (1.0_wp - kappa)
    DO n = 1, 4
      IF(i+ip .EQ. imax .OR. j+jp .EQ. jmax) THEN
        ! Right boundary
        UR(n) = U(n,i+ip,j+jp)
        CYCLE
      END IF
      DelUpHalf = U(n,i+ip,j+jp) - U(n,i,j)
      DelUppHalf = U(n,i+ipp,j+jpp) - U(n,i+ip,j+jp)
      SELECT CASE(limiter)
      CASE(0)
        phi1 = DelUpHalf
        phi2 = DelUppHalf
      CASE(1)
        !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
        phi1 = minmod(DelUpHalf, beta*DelUppHalf)
        phi2 = minmod(DelUppHalf, beta*DelUpHalf)
      END SELECT
      UR(n) = U(n,i+ip,j+jp) - 0.25_wp * epsil * (&
              (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
    END DO
  END FUNCTION RightExtrapolateU

  FUNCTION RightExtrapolateV(axis,i,j) RESULT(VR)
    USE SimulationVars_m, ONLY: V, imax, jmax
    IMPLICIT NONE
    INTEGER :: i, j, n, ip, jp, ipp, jpp
    CHARACTER(LEN=1) :: axis
    REAL(KIND=wp), DIMENSION(4) :: VR
    REAL(KIND=wp) :: DelVpHalf, DelVppHalf, rR, phi1, phi2, beta

    IF(axis .EQ. 'i') THEN
      ip = 1
      ipp = 2
      jp = 0
      jpp = 0
    ELSEIF(axis .EQ. 'j') THEN
      ip = 0
      ipp = 0
      jp = 1
      jpp = 2
    END IF
    beta = (3.0_wp - kappa) / (1.0_wp - kappa)
    DO n = 1, 4
      IF(i+ip .EQ. imax .OR. j+jp .EQ. jmax) THEN
        ! Right boundary
        VR(n) = V(n,i+ip,j+jp)
        CYCLE
      END IF
      DelVpHalf = V(n,i+ip,j+jp) - V(n,i,j)
      DelVppHalf = V(n,i+ipp,j+jpp) - V(n,i+ip,j+jp)
      SELECT CASE(limiter)
      CASE(0)
        phi1 = DelVpHalf
        phi2 = DelVppHalf
      CASE(1)
        !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
        phi1 = minmod(DelVpHalf, beta*DelVppHalf)
        phi2 = minmod(DelVppHalf, beta*DelVpHalf)
      END SELECT
      VR(n) = V(n,i+ip,j+jp) - 0.25_wp * epsil * (&
              (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
    END DO
  END FUNCTION RightExtrapolateV

  FUNCTION RightExtrapolateUC(axis,i,j) RESULT(UCR)
    USE SimulationVars_m, ONLY: UCON, imax, jmax
    IMPLICIT NONE
    INTEGER :: i, j, n, ip, jp, ipp, jpp
    CHARACTER(LEN=1) :: axis
    REAL(KIND=wp), DIMENSION(2) :: UCR
    REAL(KIND=wp) :: DelUpHalf, DelUppHalf, rR, phi1, phi2, beta

    IF(axis .EQ. 'i') THEN
      ip = 1
      ipp = 2
      jp = 0
      jpp = 0
    ELSEIF(axis .EQ. 'j') THEN
      ip = 0
      ipp = 0
      jp = 1
      jpp = 2
    END IF
    beta = (3.0_wp - kappa) / (1.0_wp - kappa)
    DO n = 1, 2
      IF(i+ip .EQ. imax .OR. j+jp .EQ. jmax) THEN
        ! Right boundary
        UCR(n) = UCON(n,i+ip,j+jp)
        CYCLE
      END IF
      DelUpHalf = UCON(n,i+ip,j+jp) - UCON(n,i,j)
      DelUppHalf = UCON(n,i+ipp,j+jpp) - UCON(n,i+ip,j+jp)
      SELECT CASE(limiter)
      CASE(0)
        phi1 = DelUpHalf
        phi2 = DelUppHalf
      CASE(1)
        !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
        phi1 = minmod(DelUpHalf, beta*DelUppHalf)
        phi2 = minmod(DelUppHalf, beta*DelUpHalf)
      END SELECT
      UCR(n) = UCON(n,i+ip,j+jp) - 0.25_wp * epsil * (&
              (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
    END DO
  END FUNCTION RightExtrapolateUC

  FUNCTION RightExtrapolateXE(axis,i,j) RESULT(XE)
    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY
    USE SimulationVars_m, ONLY: imax, jmax
    IMPLICIT NONE
    INTEGER :: i, j, n, ip, jp, ipp, jpp
    CHARACTER(LEN=1) :: axis
    ! XE(1) = XE(4) = 0.0
    ! XE(2): XI_x in i-direction or ETA_x in j-direction
    ! XE(3): XI_y in i-direction or ETA_y in j-direction
    REAL(KIND=wp), DIMENSION(4) :: XE
    REAL(KIND=wp) :: DelUpHalf, DelUppHalf, rR, phi1, phi2, beta

    XE = 0.0_wp
    beta = (3.0_wp - kappa) / (1.0_wp - kappa)
    IF(axis .EQ. 'i') THEN
      ip = 1
      ipp = 2
      jp = 0
      jpp = 0
      IF(i+ip .EQ. imax) THEN
        XE(2) = PIPX(i+ip,j+jp)
        XE(3) = PIPY(i+ip,j+jp)
      ELSE
        ! Set XI_X
        DelUpHalf = PIPX(i+ip,j+jp) - PIPX(i,j)
        DelUppHalf = PIPX(i+ipp,j+jpp) - PIPX(i+ip,j+jp)
        SELECT CASE(limiter)
        CASE(0)
          phi1 = DelUpHalf
          phi2 = DelUppHalf
        CASE(1)
          !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
          phi1 = minmod(DelUpHalf, beta*DelUppHalf)
          phi2 = minmod(DelUppHalf, beta*DelUpHalf)
        END SELECT
        XE(2) = PIPX(i+ip,j+jp) - 0.25_wp * epsil * (&
                (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
        ! Set XI_Y
        DelUpHalf = PIPY(i+ip,j+jp) - PIPY(i,j)
        DelUppHalf = PIPY(i+ipp,j+jpp) - PIPY(i+ip,j+jp)
        SELECT CASE(limiter)
        CASE(0)
          phi1 = DelUpHalf
          phi2 = DelUppHalf
        CASE(1)
          !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
          phi1 = minmod(DelUpHalf, beta*DelUppHalf)
          phi2 = minmod(DelUppHalf, beta*DelUpHalf)
        END SELECT
        XE(3) = PIPY(i+ip,j+jp) - 0.25_wp * epsil * (&
                (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
      END IF
    ELSEIF(axis .EQ. 'j') THEN
      ip = 0
      ipp = 0
      jp = 1
      jpp = 2
      IF(j+jp .EQ. jmax) THEN
        XE(2) = PJPX(i+ip,j+jp)
        XE(3) = PJPY(i+ip,j+jp)
      ELSE
        ! Set ETA_X
        DelUpHalf = PJPX(i+ip,j+jp) - PJPX(i,j)
        DelUppHalf = PJPX(i+ipp,j+jpp) - PJPX(i+ip,j+jp)
        SELECT CASE(limiter)
        CASE(0)
          phi1 = DelUpHalf
          phi2 = DelUppHalf
        CASE(1)
          !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
          phi1 = minmod(DelUpHalf, beta*DelUppHalf)
          phi2 = minmod(DelUppHalf, beta*DelUpHalf)
        END SELECT
        XE(2) = PJPX(i+ip,j+jp) - 0.25_wp * epsil * (&
                (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
        ! Set ETA_Y
        DelUpHalf = PJPY(i+ip,j+jp) - PJPY(i,j)
        DelUppHalf = PJPY(i+ipp,j+jpp) - PJPY(i+ip,j+jp)
        SELECT CASE(limiter)
        CASE(0)
          phi1 = DelUpHalf
          phi2 = DelUppHalf
        CASE(1)
          !rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
          phi1 = minmod(DelUpHalf, beta*DelUppHalf)
          phi2 = minmod(DelUppHalf, beta*DelUpHalf)
        END SELECT
        XE(3) = PJPY(i+ip,j+jp) - 0.25_wp * epsil * (&
                (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * phi2)
      END IF
    END IF
  END FUNCTION RightExtrapolateXE

  FUNCTION minmod(xx,yy) RESULT(phi)
    IMPLICIT NONE
    REAL(KIND=wp) :: xx, yy, phi

    phi = sign(1.0_wp,xx) * max(0.0_wp, min(abs(xx),sign(1.0_wp,xx)*yy))
  END FUNCTION minmod
END MODULE AUSMPWplus_m
