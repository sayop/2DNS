!> \file: AUSMPWplus.F90
!> \author: Sayop Kim

MODULE AUSMPWplus_m
   USE Parameters_m, ONLY: wp

   IMPLICIT NONE
   REAL(KIND=wp) :: alpha, kappa
   INTEGER :: limiter, epsil

CONTAINS
!-----------------------------------------------------------------------------!
  SUBROUTINE SetAUSMPWplus()
!-----------------------------------------------------------------------------!
    USE SimulationVars_m, ONLY: imin, jmin, imax, jmax, incell, jncell, &
                                U, V, DF, cgamma, nadv
    USE GridJacobian_m, ONLY: PIPX, PIPY, PJPX, PJPY, JACOBIAN

    INTEGER i, j, n
    ! UL: Left-extrapolated state vector at (i+1/2,j) or (i,j+1/2)
    ! UR: Right-extrapolated state vector at (i+1/2,j) or (i,j+1/2)
    !REAL(KIND=wp), DIMENSION(4) :: UL, UR

    ! FPhalf: Transformed flux in i-direction at i+1/2
    ! GPhalf: Transformed flux in j-direction at j+1/2
    REAL(KIND=wp), DIMENSION(4,incell,jncell) :: FPhalf, GPhalf
    REAL(KIND=wp) :: XIX, XIY, ETAX, ETAY, JACOBIANavg, A1, A2, &
                     rr, uu, vv, Pl, Pr, h0L, h0R, h0norm, &
                     utildL, vtildL, utildR, vtildR, Cavg, &
                     MtildL, MtildR, MtildLplus, MtildRminus, &
                     MbarLplus, MbarRminus, &
                     Pplus, Pminus, Pmin, fl, fr, Mavg, Ps
                     
    REAL(KIND=wp), DIMENSION(4) :: UL, UR, XEPl, XEPr

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
    !
    ! Set transformed flux along j-const lines at every (i+1/2) points
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
        h0L = Pl * cgamma / (rr * (cgamma - 1.0_wp)) + &
              0.5_wp * (uu ** 2 + vv ** 2)
        utildL = (XIX * uu + XIY * vv) / A1
        vtildL = (ETAX * uu + ETAY * vv) / A2
        ! Then calculate right-extrapolated state vector from i+1 point
        UR = RightExtrapolateU('i',i,j)
        rr = UR(1)
        uu = UR(2) / rr
        vv = UR(3) / rr
        Pr = (cgamma - 1.0_wp) * ( UR(4) - &
             0.5_wp * rr * (uu ** 2 + vv ** 2) )
        h0R = Pr * cgamma / (rr * (cgamma - 1.0_wp)) + &
              0.5_wp * (uu ** 2 + vv ** 2)
        utildR = (XIX * uu + XIY * vv) / A1
        vtildR = (ETAX * uu + ETAY * vv) / A2
        ! Calculate stagnation enthalpy normal to the interface
        h0norm = 0.5_wp * (h0L + h0R - 0.5_wp * (vtildL ** 2 + vtildR ** 2))
        ! Calculate cell averaged speed of sound
        ! Use Cavg for calculating Cs and then update to Cavg
        ! Use rr temporarily to store half of sum of utildL and utildR
        Cavg = sqrt(2.0_wp * h0norm * (cgamma - 1.0_wp) / (cgamma + 1.0_wp))
        rr = 0.5_wp * (utildL + utildR)
        IF(rr .GE. 0.0_wp) THEN
          Cavg = Cavg ** 2 / max(abs(utildL),Cavg)
        ELSE
          Cavg = Cavg ** 2 / max(abs(utildR),Cavg)
        END IF
        ! Calculate cell face Mach numbers
        MtildL = utildL / Cavg
        MtildR = utildR / Cavg
        ! Calculate split Mach numbers and pressures
        IF(abs(MtildL) .LE. 1.0_wp) THEN
          MtildLplus = 0.25_wp * (MtildL + 1.0_wp) ** 2
          Pplus = MtildLplus * (2.0_wp - MtildL) + alpha * MtildL * &
                  (MtildL ** 2 - 1.0_wp) ** 2
        ELSE
          MtildLplus = 0.5_wp * (MtildL + abs(MtildL))
          Pplus = 0.5_wp * (1.0_wp + sign(1.0_wp,MtildL))
        END IF
        IF(abs(MtildR) .LE. 1.0_wp) THEN
          MtildRminus = -0.25_wp * (MtildR - 1.0_wp) ** 2
          Pminus = -MtildRminus * (2.0_wp + MtildR) - alpha * MtildR * &
                   (MtildR ** 2 - 1.0_wp) ** 2
        ELSE
          MtildRminus = 0.5_wp * (MtildR - abs(MtildR))
          Pminus = 0.5_wp * (1.0_wp - sign(1.0_wp,MtildR))
        END IF
        ! Set pressure weighting terms
        ! Use rr temporarily to store Ps and w
        Ps = Pplus * Pl + Pminus * Pr
        Pmin = min(V(4,i,j-1), V(4,i,j+1), V(4,i+1,j-1), V(4,i+1,j+1))
        IF(Ps .LE. 0.0_wp) THEN
          fl = 0.0_wp
          fr = 0.0_wp
        ELSE
          fl = (Pl / Ps - 1.0_wp) * min(1.0_wp,Pmin/min(Pl,Pr)) ** 2
          fr = (Pr / Ps - 1.0_wp) * min(1.0_wp,Pmin/min(Pl,Pr)) ** 2
        END IF
        rr = 1.0_wp - min(Pl/Pr, Pr/Pl) ** 3
        ! Apply the weightings to MtildLplus, MtildRminus
        ! Now MtildLplus and MtildRminus become averaged Mach number
        ! on the cell face
        Mavg = MtildLplus + MtildRminus
        IF(Mavg .GE. 0.0_wp) THEN
          MbarLplus = MtildLplus + MtildRminus * ((1.0_wp - rr) * &
                       (1.0_wp + fr) - fl)
          MbarRminus = MtildRminus * rr * (1.0_wp + fr)
        ELSE
          MbarRminus = MtildRminus + MtildLplus * ((1.0_wp - rr) * &
                        (1.0_wp + fl) - fr)
          MbarLplus = MtildLplus * rr * (1.0_wp + fl)
        END IF
        ! Update FPhalf vector
!        XEPl = 0.0_wp
!        XEPl(2) = Pl * XIX
!        XEPl(3) = Pl * XIY
!        XEPr = 0.0_wp
!        XEPr(2) = Pr * XIX
!        XEPr(3) = Pr * XIY
        FPhalf(1,i,j) = MbarLplus * Cavg * A1 * UL(1) + &
                        MbarRminus * Cavg * A1 * UR(1)
        FPhalf(2,i,j) = MbarLplus * Cavg * A1 * UL(2) + &
                        MbarRminus * Cavg * A1 * UR(2) + &
                        Ps * XIX
        FPhalf(3,i,j) = MbarLplus * Cavg * A1 * UL(3) + &
                        MbarRminus * Cavg * A1 * UR(3) + &
                        Ps * XIY
        FPhalf(4,i,j) = MbarLplus * Cavg * A1 * (UL(4) + Pl) + &
                        MbarRminus * Cavg * A1 * (UR(4) + Pr)
        DO n = 1, 4
!          FPhalf(n,i,j) = ( MbarLplus * Cavg * A1 * UL(n) + &
!                            MbarRminus * Cavg * A1 * UR(n) + &
!                            Pplus * XEPl(n) + Pminus * XEPr(n) ) / &
!                            JACOBIANavg
           FPhalf(n,i,j) = FPhalf(n,i,j) / JACOBIANavg
        END DO
      END DO
    END DO

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
        h0L = Pl * cgamma / (rr * (cgamma - 1.0_wp)) + &
              0.5_wp * (uu ** 2 + vv ** 2)
        utildL = (XIX * uu + XIY * vv) / A1
        vtildL = (ETAX * uu + ETAY * vv) / A2
        ! Then calculate right-extrapolated state vector from i+1 point
        UR = RightExtrapolateU('j',i,j)
        rr = UR(1)
        uu = UR(2) / rr
        vv = UR(3) / rr
        Pr = (cgamma - 1.0_wp) * ( UR(4) - &
             0.5_wp * rr * (uu ** 2 + vv ** 2) )
        h0R = Pr * cgamma / (rr * (cgamma - 1.0_wp)) + &
              0.5_wp * (uu ** 2 + vv ** 2)
        utildR = (XIX * uu + XIY * vv) / A1
        vtildR = (ETAX * uu + ETAY * vv) / A2
        ! Calculate stagnation enthalpy normal to the interface
        h0norm = 0.5_wp * (h0L + h0R - 0.5_wp * (utildL ** 2 + utildR ** 2))
        ! Calculate cell averaged speed of sound
        ! Use Cavg for calculating Cs and then update to Cavg
        ! Use rr temporarily to store half of sum of vtildL and vtildR
        Cavg = sqrt(2.0_wp * h0norm * (cgamma - 1.0_wp) / (cgamma + 1.0_wp))
        rr = 0.5_wp * (vtildL + vtildR)
        IF(rr .GE. 0.0_wp) THEN
          Cavg = Cavg ** 2 / max(abs(vtildL),Cavg)
        ELSE
          Cavg = Cavg ** 2 / max(abs(vtildR),Cavg)
        END IF
        ! Calculate cell face Mach numbers
        MtildL = vtildL / Cavg
        MtildR = vtildR / Cavg
        ! Calculate split Mach numbers and pressures
        IF(abs(MtildL) .LE. 1.0_wp) THEN
          MtildLplus = 0.25_wp * (MtildL + 1.0_wp) ** 2
          Pplus = MtildLplus * (2.0_wp - MtildL) + alpha * MtildL * &
                  (MtildL ** 2 - 1.0_wp) ** 2
        ELSE
          MtildLplus = 0.5_wp * (MtildL + abs(MtildL))
          Pplus = 0.5_wp * (1.0_wp + sign(1.0_wp,MtildL))
        END IF
        IF(abs(MtildR) .LE. 1.0_wp) THEN
          MtildRminus = -0.25_wp * (MtildR - 1.0_wp) ** 2
          Pminus = -MtildRminus * (2.0_wp + MtildR) - alpha * MtildR * &
                   (MtildR ** 2 - 1.0_wp) ** 2
        ELSE
          MtildRminus = 0.5_wp * (MtildR - abs(MtildR))
          Pminus = 0.5_wp * (1.0_wp - sign(1.0_wp,MtildR))
        END IF
        ! Set pressure weighting terms
        Ps = Pplus * Pl + Pminus * Pr
        Pmin = min(V(4,i-1,j), V(4,i+1,j), V(4,i-1,j+1), V(4,i+1,j+1))
        IF(Ps .LE. 0.0_wp) THEN
          fl = 0.0_wp
          fr = 0.0_wp
        ELSE
          fl = (Pl / Ps - 1.0_wp) * min(1.0_wp,Pmin/min(Pl,Pr)) ** 2
          fr = (Pr / Ps - 1.0_wp) * min(1.0_wp,Pmin/min(Pl,Pr)) ** 2
        END IF
        rr = 1.0_wp - min(Pl/Pr, Pr/Pl) ** 3
        ! Apply the weightings to MtildLplus, MtildRminus
        ! Now MtildLplus and MtildRminus become averaged Mach number
        ! on the cell face
        Mavg = MtildLplus + MtildRminus
        IF(Mavg .GE. 0.0_wp) THEN
          MbarLplus = MtildLplus + MtildRminus * ((1.0_wp - rr) * &
                       (1.0_wp + fr) - fl)
          MbarRminus = MtildRminus * rr * (1.0_wp + fr)
        ELSE
          MbarRminus = MtildRminus + MtildLplus * ((1.0_wp - rr) * &
                        (1.0_wp + fl) - fr)
          MbarLplus = MtildLplus * rr * (1.0_wp + fl)
        END IF
        ! Update GPhalf vector
!        XEPl = 0.0_wp
!        XEPl(2) = Pl * ETAX
!        XEPl(3) = Pl * ETAY
!        XEPr = 0.0_wp
!        XEPr(2) = Pr * ETAX
!        XEPr(3) = Pr * ETAY
        GPhalf(1,i,j) = MbarLplus * Cavg * A2 * UL(1) + &
                        MbarRminus * Cavg * A2 * UR(1)
        GPhalf(2,i,j) = MbarLplus * Cavg * A2 * UL(2) + &
                        MbarRminus * Cavg * A2 * UR(2) + &
                        Ps * ETAX
        GPhalf(3,i,j) = MbarLplus * Cavg * A2 * UL(3) + &
                        MbarRminus * Cavg * A2 * UR(3) + &
                        Ps * ETAY
        GPhalf(4,i,j) = MbarLplus * Cavg * A2 * (UL(4) + Pl) + &
                        MbarRminus * Cavg * A2 * (UR(4) + Pr)
        DO n = 1, 4
!          GPhalf(n,i,j) = ( MbarLplus * Cavg * A2 * UL(n) + &
!                            MbarRminus * Cavg * A2 * UR(n) + &
!                            Pplus * XEPl(n) + Pminus * XEPr(n) ) / &
!                            JACOBIANavg
           GPhalf(n,i,j) = GPhalf(n,i,j) / JACOBIANavg
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
    USE SimulationVars_m, ONLY: U, imin, jmin, imax, jmax
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
    !beta = (3.0_wp - kappa) / (1.0_wp - kappa)
    DO n = 1, 4
      IF(i .EQ. imin .OR. j .EQ. jmin) THEN
        ! Left boundary
        UL(n) = U(n,i,j)
        CYCLE
      END IF
      DelUmHalf = U(n,i,j) - U(n,i+im,j+jm)
      DelUpHalf = U(n,i+ip,j+jp) - U(n,i,j)
      rL = (DelUpHalf + 1.0E-6) / (DelUmHalf + 1.0E-6)
      SELECT CASE(limiter)
      CASE(0)
        phi1 = 1.0_wp
        phi2 = 1.0_wp
      CASE(1)
        IF(rL .GE. 0.0_wp) THEN
          phi1 = minmod(rL, 1.0_wp)
          phi2 = minmod(1.0_wp / rL, 1.0_wp)
        ELSE
          phi1 = 0.0_wp
          phi2 = 0.0_wp
        END IF
      END SELECT
      UL(n) = U(n,i,j) + 0.25_wp * epsil * DelUmHalf * ( &
              (1.0_wp - kappa) * phi1 + (1.0_wp + kappa) * rL * phi2 )
    END DO
  END FUNCTION LeftExtrapolateU


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
    !beta = (3.0_wp - kappa) / (1.0_wp - kappa)
    DO n = 1, 4
      IF(i+ip .EQ. imax .OR. j+jp .EQ. jmax) THEN
        ! Right boundary
        UR(n) = U(n,i+ip,j+jp)
        CYCLE
      END IF
      DelUpHalf = U(n,i+ip,j+jp) - U(n,i,j)
      DelUppHalf = U(n,i+ipp,j+jpp) - U(n,i+ip,j+jp)
      rR = (DelUpHalf + 1.0E-6) / (DelUppHalf + 1.0E-6)
      SELECT CASE(limiter)
      CASE(0)
        phi1 = 1.0_wp
        phi2 = 1.0_wp
      CASE(1)
        IF(rR .GE. 0.0_wp) THEN
          phi1 = minmod(rR, 1.0_wp)
          phi2 = minmod(1.0_wp / rR, 1.0_wp)
        ELSE
          phi1 = 0.0_wp
          phi2 = 0.0_wp
        END IF
      END SELECT
      UR(n) = U(n,i+ip,j+jp) - 0.25_wp * epsil * DelUppHalf * ( &
              (1.0_wp + kappa) * rR * phi2 + (1.0_wp - kappa) * phi1 )
    END DO
  END FUNCTION RightExtrapolateU


  FUNCTION minmod(xx,yy) RESULT(phi)
    IMPLICIT NONE
    REAL(KIND=wp) :: xx, yy, phi

    phi = sign(1.0_wp,xx) * max(0.0_wp, min(abs(xx),sign(1.0_wp,xx)*yy))
  END FUNCTION minmod
END MODULE AUSMPWplus_m
