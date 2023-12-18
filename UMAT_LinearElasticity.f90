SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, &
                DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, &
                DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, &
                NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, &
                DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC) 

  !DEC$ ATTRIBUTES DLLExport,StdCall,reference :: umat
  use, intrinsic :: iso_c_binding
  implicit double precision(a-h,o-z)
  parameter (nprecd=2)

  CHARACTER(C_CHAR) :: CMNAME
  REAL(C_DOUBLE), DIMENSION(NTENS) :: STRESS, STATEV, DDSDDT, DRPLDE, STRAN, DSTRAN, PREDEF, DPRED, Sig, dSig
  REAL(C_DOUBLE), DIMENSION(NTENS, NTENS) :: DDSDDE, DFGRD0, DFGRD1
  REAL(C_DOUBLE), DIMENSION(3) :: COORDS, DROT, TIME
  REAL(C_DOUBLE), DIMENSION(NPROPS) :: PROPS
  REAL(C_DOUBLE) :: SSE, SPD, SCD, RPL, DTIME, TEMP, DTEMP, PNEWDT
  INTEGER(C_INT) :: NDI, NSHR, NTENS, NSTATEV, NPROPS, CELENT, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
  REAL(C_DOUBLE) :: G, ENU, one, two, FAC, D1, D2, DSTRANVOL
  INTEGER(C_INT) :: i

  ! Arguments:
  !          I/O  Type
  !  PROPS    I   R()  : List with model parameters
  !  DSTRAN   I   R()  : Strain increment
  !  DDSDDE   O   R(,) : Material stiffness matrix
  !  STRESS  I/O  R()  : stresses
  !  STATEV  I/O  R()  : state variables

  ! Contents of PROPS(2)
  !  1 : G       shear modulus
  !  2 : ENU     Poisson's ratio

  G = PROPS(1)
  ENU = PROPS(2)
  one = 1.0d0
  two = 2.0d0

  ! calculate elastic stress increment (DSigE = elastic stiffness D * strain increment DEps)
  FAC = two * G / (one - two * ENU)
  D1 = FAC * (one - ENU)
  D2 = FAC * ENU
  DSTRANVOL = DSTRAN(1) + DSTRAN(2) + DSTRAN(3)
  dSig(1) = (D1 - D2) * DSTRAN(1) + D2 * DSTRANVOL
  dSig(2) = (D1 - D2) * DSTRAN(2) + D2 * DSTRANVOL
  dSig(3) = (D1 - D2) * DSTRAN(3) + D2 * DSTRANVOL
  dSig(4) = G * DSTRAN(4)

  IF (NTENS == 6) then
    dSig(5) = G * DSTRAN(5)
    dSig(6) = G * DSTRAN(6)
  END IF

  ! elastic stress
  Sig = STRESS + dSig

  ! stress state parameters update
  DO i = 1, NTENS
    STRESS(i) = Sig(i)
  END DO

  DDSDDE = 0.0
  DDSDDE(1:3, 1:3) = D2
  DDSDDE(1, 1) = D1
  DDSDDE(2, 2) = D1
  DDSDDE(3, 3) = D1
  DDSDDE(4, 4) = G
  IF (NTENS == 6) then
    DDSDDE(5,5) = G
    DDSDDE(6,6) = G
  END IF

  RETURN

END SUBROUTINE UMAT
