*USER SUBROUTINES
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATEV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      implicit real(8) (a-h,o-z)
!
      CHARACTER*(*) CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATEV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
!
! Arguments:
!          I/O  Type
!  PROPS    I   R()  : List with model parameters
!  DSTRAN   I   R()  : Strain increment
!  DDSDDE   O   R(,) : Material stiffness matrix
!  STRESS  I/O  R()  : stresses
!  STATEV  I/O  R()  : state variables
!
!
!---  Local variables
!
      Dimension DE(6,6), dSig(6), Prs_E(3), Prs(3),
     *          xN1(3), xN2(3), xN3(3),Sig(6)

    ! this statement needs to be activated if a dll is going to be created on Windows
    !_!DEC$ ATTRIBUTES DLLExport,StdCall,reference :: umat

!
!     Mohr-Coulomb modelwith tension cut-off
!     date 6 february 2015
!     Hans Teunissen
!
! Contents of PROPS(6) MC
!  1 : G       shear modulus
!  2 : ENU     Poisson's ratio
!  3 : C       Cohesion
!  4 : Phi     Friction angle (degrees)
!  5 : Psi     Dilation angle (degrees)
!  6 : Tens    Allowable tensile stress
!
        Rad  = 45d0 / datan(1d0)
*
* ... start correction routine
*
        call AbaqusToPlaxisStressOrdering(DSTRAN)
        call AbaqusToPlaxisStressOrdering(STRESS)
        ipl     =   0
        G       =   PROPS(1)       ! G
        ENU     =   PROPS(2)       ! nu
        C       =   PROPS(3)       ! C
        Phi     =   PROPS(4) / Rad ! Phi in radians
        Psi     =   PROPS(5) / Rad ! Psi in radians
        sTens   =  -PROPS(6)       ! tensile strength (change sign)
        sPhi    =   Sin(Phi)
        sPsi    =   Sin(Psi)
        cCosPhi =   C*Cos(Phi)
        If (sPhi.Gt.0) Then
          If (sTens.Gt.cCosPhi/sPhi) sTens = cCosPhi/sPhi
        End If
*
        ! Fill elastic material matrix
        F1  = 2*G*(1-ENU)/(1-2*ENU)
        F2  = 2*G*( ENU )/(1-2*ENU)
        Call MZeroR(DE,36)
        Do i=1,3
          Do j=1,3
            DE(i,j) = F2
          End Do
          DE(i,i) = F1
          DE(i+3,i+3) = G
        End Do
*
        ! elastic stress increment
        Call MatVec( DE, 6, DSTRAN, 6, dSig)
        ! elastic stress
        Call AddVec( STRESS, dSig, 1d0, 1d0, 6, Sig )
        ! calculate principal stresses and directions
        iOpt = 1
        Call PrnSig(iOpt, Sig, xN1, xN2, xN3, S1, S2, S3, P, Q)
        ! Sig    : tension     positive
        ! Prs(E) : compression positive
        Prs_E(1) = - S1 ! minus sign
        Prs_E(2) = - S2
        Prs_E(3) = - S3
        iArea = 2
        Call MC_Tens( iArea, G, ENU, sPhi, sPsi, cCosPhi, sTens,
     *                Prs_E, Prs, ipl,
     *                fme,fte,fhe,fm,ft,fh ,xLamM,xLamT)
          If (ipl.Ne.0) Then ! some plasticity
          ! Check Sig1 > Sig2 > Sig3
          If (                 Prs(2).Lt.Prs(3)) iarea = 1   ! Tr. compression
          If (IArea.Eq.2 .And. Prs(1).Lt.Prs(2)) iarea = 3   ! Tr. extension
          If (iArea.Ne.2) Then
            Call MC_Tens( iArea, G, ENU, sPhi, sPsi, cCosPhi, sTens,
     *                    Prs_E, Prs, ipl,
     *                    fme,fte,fhe,fm,ft,fh ,xLamM,xLamT)
          End If
          ! Prs : compression positive
          S1 = - Prs(1) ! minus sign
          S2 = - Prs(2)
          S3  =- Prs(3)
          ! back to Cartesian stresses
          Call CarSig(S1,S2,S3,xN1,xN2,xN3,Sig)
          ! Sig    : tension positive
        End If
*
* ... stress state parameters update
*
        call AbaqusToPlaxisStressOrdering(Sig)
        Do i=1,NTENS
          STRESS(i) = sig(i)
        End Do
        statev(1)=dble(ipl)
        sse=0d0
        spd=0d0
        scd=0d0
*
* ... Tangent stiffness matrix to be returned (done by elastic stiffness) no correction needed!
*
        G       =   PROPS(1)       ! G
        ENU     =   PROPS(2)       ! nu
        F1  = 2*G*(1-ENU)/(1-2*ENU)
        F2  = 2*G*( ENU )/(1-2*ENU)
        Call MZeroR(DDSDDE,36)
        Do i=1,3
          Do j=1,3
            DDSDDE(i,j) = F2
          End Do
          DDSDDE(i,i) = F1
          DDSDDE(i+3,i+3) = G
        End Do
*
* ... end UMAT routine
*
      Return
      End
      subroutine AbaqusToPlaxisStressOrdering(S)
      implicit real(8) (a-h,o-z)
      dimension S(*)
      temp=S(5)
      S(5)=S(6)
      S(6)=temp
      end subroutine AbaqusToPlaxisStressOrdering

C***********************************************************************
      Subroutine MC_Tens( iArea, G, xNu, sPhi, sPsi, cCosPhi, sTens,
     *                    Prs_E, Prs, ipl
     *                    ,fme,fte,fhe,fm,ft,fh ,xLamM,xLamT
     *                   )
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
!
! Routine to solve stresses according to Coulomb / tension criterion
! Caution : compression positive
!
!     Sig = Sig^E - xLamM * DdGmdS - xLamT * DdGtdS
!
!        I/O Type
! iArea   I    I    : 1 : Triax compression corner     Sig1 > Sig2 = Sig3
!                     2 : Regular yield surface        Sig1 > Sig2 > Sig3
!                     3 : Triax extension              Sig1 = Sig2 > Sig3
! G       I    R    : Shear modulus
! xNu     I    R    : Poisson's ratio
! sPhi    I    R    : Sine of friction angle
! sPsi    I    R    : Sine of dilation angle
! cCosPhi I    R    : Cohesion * Cos(phi)
! sTens   I    R    : Allowable tensile stress (negative value)
! Prs_E   I    R(3) : Elastic principal stresses   Sig1E >= Sig2E >= Sig3E
! Prs     O    R(3) : Resulting principal stresses
! ipl     O    I    : Plasticity indicator
!                     0 : elastic
!                     1 : Coulomb surface
!                     2 : Tension surface (also corner with Coulomb)
! Local:
!   fme,fte,fm,ft : Value of yield functions m=Coulomb, t=tension, e=elastic
!   xLamM,xLamT   : Plastic multipliers
!   DdG#dS(3)     : Elastic D matrix times derivative of plastic potential
!   PrsE(3)       : Copy of Prs_E with possible correction for corners
!
      Dimension Prs_E(3),Prs(3)
      Dimension DdGmdS(3),DdGtdS(3),dPrs(3),PrsE(3)
      Do i=1,3
        PrsE(i) = Prs_E(i)
        Prs (i) = Prs_E(i)
      End Do
      ipl = 0
      xLamM = 0
      xLamT = 0
      ft=-1
      fm=-1
      Select Case (IArea)
        Case (1) ! Triax compression ! Sig1 > Sig2 = Sig3
          !!!!!!!!!!!!!!!!!!! COMPRESSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! first correct sig2=sig3
          PrsE(2) = (Prs_E(2) + Prs_E(3) ) / 2
          PrsE(3) = PrsE(2)
          Prs (2) = PrsE(2)
          Prs (3) = PrsE(3)

          fme = ( 2*PrsE(1) - PrsE(2) - PrsE(3) ) /4
     *         -( 2*PrsE(1) + PrsE(2) + PrsE(3) ) /4 * sPhi - cCosPhi
          fte = sTens - ( PrsE(2) + PrsE(3) ) /2
          If (fme.Lt.0 .And. fte.Lt.0 ) Goto 999 ! Return Both elastic
          fmTol = 1d-10*abs(fme)+1d-10
          ftTol = 1d-10*abs(fte)+1d-10

          ! Try Coulomb only
          DdGmdS(1) = G/(1-2*xNu)/2 * ( 2*(1-2*xNu)  -(2      )*sPsi )
          DdGmdS(2) = G/(1-2*xNu)/2 * ( - (1-2*xNu)  -(1+2*xNu)*sPsi )
          DdGmdS(3) = G/(1-2*xNu)/2 * ( - (1-2*xNu)  -(1+2*xNu)*sPsi )
          amm =G/4 * (3  - sPsi - sPhi + sPhi*sPsi*(3+2*xNu)/(1-2*xNu) )
          If (fme.Gt.0) Then
            xLamM = fme / amm  ! assume only coulomb surface
            xLamT = 0
            ! Calculate stresses Prs = Prs - xLamM*DdGmdS
            Call AddVec( PrsE, DdGmdS, 1d0, -xLamM, 3, Prs )
            ! check yield functions
            fm = ( 2*Prs(1) - Prs(2) - Prs(3) ) /4
     *          -( 2*Prs(1) + Prs(2) + Prs(3) ) /4 * sPhi - cCosPhi
            ft = sTens - ( Prs(2) + Prs(3) ) /2
            ! when correct return
            ipl = 1 ! Coulomb
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Goto 999 ! Return
          End If
!         coming here : tension surface active
          ipl = 2 ! tens
          DdGtdS(1) =  -G/(1-2*xNu) * ( 2* xNu  )
          DdGtdS(2) =  -G/(1-2*xNu) * (  1      )
          DdGtdS(3) =  -G/(1-2*xNu) * (  1      )
          atm =G/2 * ( 1 + sPsi*(1+2*xNu)/(1-2*xNu) )
          amt =G/2 * ( 1 + sPhi*(1+2*xNu)/(1-2*xNu) )
          att = G /(1-2*xNu)
          If (fte.Gt.0) Then
            xLamT = fte/att
            xLamM = 0
            ! Prs = PrsE - xLamT * DdGtdS
            Call AddVec( PrsE, DdGtdS, 1d0, -xLamT, 3, Prs )
            fm = ( 2*Prs(1) - Prs(2) - Prs(3) ) /4
     *          -( 2*Prs(1) + Prs(2) + Prs(3) ) /4 * sPhi - cCosPhi
            ft = sTens - ( Prs(2) + Prs(3) ) /2
            ft1= sTens - Prs(1)
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Then ! Goto 999 ! Return
              If (ft1.Gt.0) Prs(1) = sTens
              Goto 999 ! Return
            End If
          End If
          ! combined: coulomb + tension
          ipl = 2 ! tens (+compr)
          Det = amm*att-amt*atm
          xLamM = ( att*fme - amt*fte ) / det
          xLamT = (-atm*fme + amm*fte ) / det
          ! Prs = PrsE - xLamM * DdGmdS - xLamT * DdGtdS
          Call AddVec( DdGmdS, DdGtdS, xLamM, xLamT, 3, dPrs )
          Call AddVec( PrsE, dPrs, 1d0, -1d0, 3, Prs )
          fm = ( 2*Prs(1) - Prs(2) - Prs(3) ) /4
     *        -( 2*Prs(1) + Prs(2) + Prs(3) ) /4 * sPhi - cCosPhi
          ft = sTens - ( Prs(2) + Prs(3) ) /2
          ft1= sTens - Prs(1)
          if (ft1.gt.1d-10) write(10,*)'FT1 positive !!!'
          If (ft1.Gt.0) Then
            Prs(1) = sTens
          End If
          If ( fm .Lt. fmTol .And.
     *         ft .Lt. ftTol       ) Goto 999 ! Return
        Case (2) ! Regular
          !!!!!!!!!!!!!!!!!!! REGULAR !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          fme = ( PrsE(1) - PrsE(3) ) /2
     *         -( PrsE(1) + PrsE(3) ) /2 * sPhi - cCosPhi
          fte = sTens - PrsE(3)
          If (fme.Lt.0 .And. fte.Lt.0 ) Goto 999 ! Return ! Both elastic
          fmTol = 1d-10*abs(fme)+1d-10
          ftTol = 1d-10*abs(fte)+1d-10
          ! Try Coulomb only
          DdGmdS(1) = G/(1-2*xNu) * ( 1-2*xNu -     sPsi )
          DdGmdS(2) = G/(1-2*xNu) * (        -2*xNu*sPsi )
          DdGmdS(3) = G/(1-2*xNu) * (-1+2*xNu -     sPsi )
          amm = G * ( 1 + sPhi*sPsi/(1-2*xNu) )
          If (fme.Gt.0) Then
            xLamM = fme / amm
            xLamT = 0
            ! Prs = PrsE - xLamM * DdGmdS
            Call AddVec( PrsE, DdGmdS, 1d0, -xLamM, 3, Prs )
            fm = ( Prs(1) - Prs(3) ) /2
     *          -( Prs(1) + Prs(3) ) /2 * sPhi - cCosPhi
            ft = sTens - Prs(3)
            ipl = 1 ! Coulomb
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Goto 999 ! Return
          End If
          ! tension surface
          ipl = 2 ! tens
          DdGtdS(1) = 2*G/(1-2*xNu) * ( -  xNu   )
          DdGtdS(2) = 2*G/(1-2*xNu) * ( -  xNu   )
          DdGtdS(3) = 2*G/(1-2*xNu) * ( -(1-xNu) )
          atm =   G * ( 1 + sPsi/(1-2*xNu) )
          amt =   G * ( 1 + sPhi/(1-2*xNu) )
          att = 2*G * ( 1-xNu)/(1-2*xNu)
          If (fte.Gt.0) Then
            xLamT = fte/att
            xLamM = 0
            ! Prs = PrsE - xLamT * DdGtdS
            Call AddVec( PrsE, DdGtdS, 1d0, -xLamT, 3, Prs )
            fm = ( Prs(1) - Prs(3) ) /2
     *          -( Prs(1) + Prs(3) ) /2 * sPhi - cCosPhi
            ft = sTens - Prs(3)
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Goto 999 ! Return
          End If
          ! coulomb + tension
          Det = amm*att-amt*atm
          xLamM = ( att*fme - amt*fte ) / det
          xLamT = (-atm*fme + amm*fte ) / det
          ! Prs = PrsE - xLamM * DdGmdS - xLamT * DdGtdS
          Call AddVec( DdGmdS, DdGtdS, xLamM, xLamT, 3, dPrs )
          Call AddVec( PrsE, dPrs, 1d0, -1d0, 3, Prs )
          fm = ( Prs(1) - Prs(3) ) /2
     *        -( Prs(1) + Prs(3) ) /2 * sPhi - cCosPhi
          ft = sTens - Prs(3)
          ft1= sTens - Prs(1)
          if (ft1.gt.1d-10) write(10,*)'FT1 positive !!!'
          If ( fm .Lt. fmTol .And.
     *         ft .Lt. ftTol       ) Goto 999 ! Return
        Case (3) ! Triax extension
          !!!!!!!!!!!!!!!!!!! EXTENSION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          ! first correct sig1=sig2
          PrsE(2) = (Prs_E(1) + Prs_E(2) ) / 2
          PrsE(1) = PrsE(2)
          Prs (1) = PrsE(1)
          Prs (2) = PrsE(2)

          fme = ( PrsE(1) + PrsE(2) - 2*PrsE(3) ) /4
     *         -( PrsE(1) + PrsE(2) + 2*PrsE(3) ) /4 * sPhi - cCosPhi
          fte = sTens - PrsE(3)
          If (fme.Lt.0 .And. fte.Lt.0 ) Goto 999 ! Return ! Both elastic
          fmTol = 1d-10*abs(fme)+1d-10
          ftTol = 1d-10*abs(fte)+1d-10

          ! Try Coulomb only
          DdGmds(1) = G/(1-2*xNu)/2 * (   (1-2*xNu)  -(1+2*xNu)*sPsi )
          DdGmds(2) = G/(1-2*xNu)/2 * (   (1-2*xNu)  -(1+2*xNu)*sPsi )
          DdGmds(3) = G/(1-2*xNu)/2 * (-2*(1-2*xNu)  -(2      )*sPsi )
          amm =G/4 * (3 + sPsi + sPhi + sPhi*sPsi*(3+2*xNu)/(1-2*xNu) )
          If (fme.Gt.0) Then
            xLamM = fme / amm
            xLamT = 0
            ! Prs = PrsE - xLamM * DdGmdS
            Call AddVec( PrsE, DdGmdS, 1d0, -xLamM, 3, Prs )
            fm = ( Prs(1) + Prs(2) - 2*Prs(3) ) /4
     *          -( Prs(1) + Prs(2) + 2*Prs(3) ) /4 * sPhi - cCosPhi
            ft = sTens -  Prs(3)
            ipl = 1 ! Coulomb
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Goto 999 ! Return
          End If
          ! tension surface
          ipl = 2 ! tens
          DdGtds(1) =  -2*G/(1-2*xNu) * (  xNu  )
          DdGtds(2) =  -2*G/(1-2*xNu) * (  xNu  )
          DdGtds(3) =  -2*G/(1-2*xNu) * ( 1-xNu )

          amt =G * ( 1 + sPhi /(1-2*xNu) )
          atm =G * ( 1 + sPsi /(1-2*xNu) )
          att = 2*G*(1-xNu) /(1-2*xNu)

          If (fte.Gt.0) Then
            xLamT = fte/att
            xLamM = 0
            ! Prs = PrsE - xLamT * DdGtdS
            Call AddVec( PrsE, DdGtdS, 1d0, -xLamT, 3, Prs )
            fm = ( Prs(1) + Prs(2) - 2*Prs(3) ) /4
     *          -( Prs(1) + Prs(2) + 2*Prs(3) ) /4 * sPhi - cCosPhi
            ft = sTens -  Prs(3)

            ft1= sTens - Prs(1)
            If ( fm .Lt. fmTol .And.
     *           ft .Lt. ftTol       ) Then ! Goto 999 ! Return
              If (ft1.Gt.0) Then
                Prs(1) = sTens
                Prs(2) = sTens
              End If
              Goto 999 ! Return
            End If
          End If
          ! coulomb + tension
          Det = amm*att-amt*atm
          xLamM = ( att*fme - amt*fte ) / det
          xLamT = (-atm*fme + amm*fte ) / det
          ! Prs = PrsE - xLamM * DdGmdS - xLamT * DdGtdS
          Call AddVec( DdGmdS, DdGtdS, xLamM, xLamT, 3, dPrs )
          Call AddVec( PrsE, dPrs, 1d0, -1d0, 3, Prs )
          fm = ( 2*Prs(1) - Prs(2) - Prs(3) ) /4
     *        -( 2*Prs(1) + Prs(2) + Prs(3) ) /4 * sPhi - cCosPhi
          ft = sTens - ( Prs(2) + Prs(3) ) /2
          ft1= sTens - Prs(1)
          if (ft1.gt.1d-10) write(10,*)'FT1 positive !!!'
          If ( fm .Lt. fmTol .And.
     *         ft .Lt. ftTol       ) Goto 999 ! Return
        Case Default
          Stop ' incorrect iArea in MC_Tens '
      End Select
      ! Normally, no situation should come here
      write(10,*)'Area: ',iArea
      write(10,*)'fme : ',fme
      write(10,*)'fte : ',fte
      write(10,*)'fm  : ',fm
      write(10,*)'ft  : ',ft
!      Call Error ( ' how did we get here ?' )
!      Stop ' how did we get here ?'
  999 Continue
      If (fm.Gt.1d-10 .Or. ft.Gt.1d-10) Then
        Write(10,902) iarea,fm,ft
  902   Format('fm>0 or ft>0, iarea:',i5,2f10.5)
      End If
      Return

      End ! MC_Tens
      Subroutine MZEROR(R,K)
C
C***********************************************************************
C
C     Function: To make a real array R with dimension K to zero
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J) = 0.0D0
      End Do

      Return
      End


      Subroutine MZEROI(I,K)
C
C***********************************************************************
C
C     Function: To make an integre array I with Dimension K to zero
C
C***********************************************************************
C
      Dimension I(*)

      Do J=1,K
        I(J)=0
      End Do

      Return
      End

      Subroutine SETRVAL(R,K,V)
C
C***********************************************************************
C
C     Function: To fill a real array R with Dimension K with value V
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J)=V
      End Do

      Return
      End

      Subroutine SETIVAL(I,K,IV)
C
C***********************************************************************
C
C     Function: To fill an integer array I with Dimension K with value IV
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension I(*)

      Do J=1,K
        I(J)=IV
      End Do

      Return
      End

      Subroutine COPYIVEC(I1,I2,K)
C
C***********************************************************************
C
C     Function: To copy an integer array I1 with Dimension K to I2
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension I1(*),I2(*)

      Do  J=1,K
        I2(J)=I1(J)
      End Do

      Return
      End

      Subroutine COPYRVEC(R1,R2,K)
C
C***********************************************************************
C
C     Function: To copy a Double array R1 with Dimension K to R2
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension R1(*),R2(*)

      Do J=1,K
        R2(J)=R1(J)
      End Do

      Return
      End


      Logical Function IS0ARR(A,N)
C
C***********************************************************************
C    Function :  To check whether a real array contains only zero values.
C                When an array contains only zero's is might not need to be
C                written to the XXX file.
C                exit Function when first non-zero value occured or when
C                all elements are checked and are zero.
C
C    Input:  A : array to be checked
C            N : number of elements in array that should be checked
C
C    Output : .TRUE.  when all elements are 0
C             .FALSE. when at least one element is not zero
C
C    Called by :  Subroutine TOBXX
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*)
      Is0Arr=.False.
      Do I=1,N
        If ( A(I) .Ne. 0 ) Return
      End Do
      Is0Arr=.True.
      Return
      End

      Logical Function IS0IARR(IARR,N)
C
C***********************************************************************
C    Function :  To check whether a integer array contains only zero values.
C                Similar to IS0ARR
C
C    Input:  IARR : array to be checked
C            N    : number of elements in array that should be checked
C
C    Output : .TRUE.  when all elements are 0
C             .FALSE. when at least one element is not zero
C
C    Called by :  Subroutine TOBXX
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Dimension IARR(*)

      Is0IArr=.False.
      Do I=1,N
        If ( IARR(I) .Ne. 0 ) Return
      End Do
      Is0IArr=.True.
      Return
      End
C***********************************************************************
      Subroutine MulVec(V,F,K)
C***********************************************************************
C
C     Function: To multiply a real vector V with dimension K by F
C
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*)

      Do J=1,K
        V(J)=F*V(J)
      End Do

      Return
      End     ! Subroutine Mulvec
C***********************************************************************
      Subroutine MatVec(xMat,IM,Vec,N,VecR)
C***********************************************************************
C
C     Calculate VecR = xMat*Vec
C
C I   xMat  : (Square) Matrix (IM,*)
C I   Vec   : Vector
C I   N     : Number of rows/colums
C O   VecR  : Resulting vector
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
C***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine MatVec

C***********************************************************************
      Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
C***********************************************************************
C
C     Calculate VecR() = R1*Vec1()+R2*Vec2()
C
C I   Vec1,
C I   Vec2  : Vectors
C I   R1,R2 : Multipliers
C I   N     : Number of rows
C O   VecR  : Resulting vector
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
C***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine AddVec
C
C***********************************************************************
      Double Precision Function DInProd(A,B,N)
C***********************************************************************
C
C     Returns the Inproduct of two vectors
C
C I   A,B  : Two vectors
C I   N    : Used length of vectors
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*),B(*)
C***********************************************************************

      X = 0
      Do I=1,N
        X = X + A(I)*B(I)
      End Do
      DInProd = X
      Return
      End     ! Function DInProd
C
C***********************************************************************
      Subroutine MatMat(xMat1,Id1,xMat2,Id2,nR1,nC2,nC1,xMatR,IdR)
C***********************************************************************
C
C     Calculate xMatR = xMat1*xMat2
C
C I   xMat1 : Matrix (Id1,*)
C I   xMat2 : Matrix (Id2,*)
C I   nR1   : Number of rows in resulting matrix    (= No rows in xMat1)
C I   nC2   : Number of columns in resulting matrix (= No cols in xMat2)
C I   nC1   : Number of columns in matrix xMat1
C             = Number  rows    in matrix xMat2
C O   xMatR : Resulting matrix (IdR,*)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(Id1,*),xMat2(Id2,*),xMatR(IdR,*)
C**********************************************************************

      Do I=1,nR1
        Do J=1,nC2
          X=0
          Do K=1,nC1
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMat

C***********************************************************************
      Subroutine MatMatSq(n, xMat1, xMat2, xMatR)
C***********************************************************************
C
C     Calculate xMatR = xMat1*xMat2 for square matrices, size n
C
C I   n     : Dimension of matrices
C I   xMat1 : Matrix (n,*)
C I   xMat2 : Matrix (n,*)
C O   xMatR : Resulting matrix (n,*)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(n,*),xMat2(n,*),xMatR(n,*)
C**********************************************************************

      Do I=1,n
        Do J=1,n
          X=0
          Do K=1,n
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMatSq

C***********************************************************************
      Subroutine WriVal ( io, C , V )
C***********************************************************************
C
C Write (Double) value to file unit io (when io>0)
C
C***********************************************************************
C
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,V
    1 Format( A,3x, 1x,1p,e12.5)
      Return
      End
C***********************************************************************
      Subroutine WriIVl ( io, C , I )
C***********************************************************************
C
C Write (integer) value to file unit io (when io>0)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,I
    1 Format( A,3x, 1x,I6)
      Return
      End
C***********************************************************************
      Subroutine WriIVc ( io, C , iV , n )
C***********************************************************************
C
C Write (integer) vector to file unit io (when io>0)
C
C***********************************************************************
      Character C*(*)
      Dimension iV(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) (iv(i),i=1,n)
    1 Format( ( 2(3x,5i4) ) )
      Return
      End
C***********************************************************************
      Subroutine WriVec ( io, C , V , n )
C***********************************************************************
C
C Write (Double) vector to file unit io (when io>0)
C 6 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      If (Len_Trim(C).Le.6) Then
        Write(io,2) C,( V(i),i=1,n)
      Else
        Write(io,*) C
        Write(io,1) ( V(i),i=1,n)
      End If
    1 Format( ( 2(1x, 3(1x,1p,e10.3) ) ) )
    2 Format( A, ( T7, 2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
C***********************************************************************
      Subroutine WriVec5( io, C , V , n )
C***********************************************************************
C
C Write (Double) vector to file unit io (when io>0)
C 5 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) ( V(i),i=1,n)
    1 Format( 5(1x,1p,e12.5) )
      Return
      End
C***********************************************************************
      Subroutine WriMat ( io, C , V , nd, nr, nc )
C***********************************************************************
C
C Write (Double) matrix to file unit io (when io>0)
C 6 values per line
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(nd,*)

      If (io.Le.0) Return

      Write(io,*) C
      Do j=1,nr
        Write(io,1) j,( V(j,i),i=1,nc)
      End Do
    1 Format(i4, (  T7,2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
C***********************************************************************

      Subroutine MatInvPiv(Aorig,B,N)
      Implicit Double Precision (A-H,O-Z)
      Dimension Aorig(n,*), B(n,*),A(:,:)
      Allocatable :: A
      Allocate ( A(n,n) ) ! No error checking !!
      Call CopyRVec(AOrig, A, n*n )
      Call MZeroR(B,n*n)
      Do i=1,n
        B(i,i) = 1d0
      End Do
      Do I=1,n
        T=A(I,I)
        iPiv=i
        Do j=i+1,n
          If ( Abs(A(j,i)) .Gt. Abs(A(iPiv,i))  ) iPiv=j
        End Do
        If (iPiv.Ne.i) Then
          Do j=1,n
            x         = A( i  ,j)
            A( i  ,j) = A(iPiv,j)
            A(iPiv,j) = x
            x         = B( i  ,j)
            B( i  ,j) = B(iPiv,j)
            B(iPiv,j) = x
          End Do
          T=A(I,I)
        End If
        Do J=1,n
          A(I,J)=A(I,J)/T
          B(I,J)=B(I,J)/T
        End Do
        Do K=1,n
          If (K.Ne.I) Then
            T=A(K,I)
            Do J=1,n
              A(K,J)=A(K,J)-T*A(I,J)
              B(K,J)=B(K,J)-T*B(I,J)
            End Do
          End If
        End Do
      End Do
      DeAllocate ( A  )
      Return
      End ! MatinvPiv

C***********************************************************************
      Subroutine PrnSig(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension S(*),xN1(*),xN2(*),xN3(*)

      If (iOpt.Eq.1) Then
        Call Eig_3(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
      Else
        Call Eig_3a(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
      End If
      Return
      End
C***********************************************************************
      Subroutine Eig_3(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3),V(3,3),
     *          xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! Set V to unity matrix
      V(1,1) = 1
      V(2,1) = 0
      V(3,1) = 0

      V(1,2) = 0
      V(2,2) = 1
      V(3,2) = 0

      V(1,3) = 0
      V(2,3) = 0
      V(3,3) = 1


      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      it = 0
      itmax = 50
      Do While ( it.Lt.itMax .And.
     *           abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
          If (a(ip,iq) .Ne. 0.0) Then
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            v1p=c*v(1,ip)-s*v(1,iq)
            v2p=c*v(2,ip)-s*v(2,iq)
            v3p=c*v(3,ip)-s*v(3,iq)
            v(1,iq)=s*v(1,ip)+c*v(1,iq)
            v(2,iq)=s*v(2,ip)+c*v(2,iq)
            v(3,iq)=s*v(3,ip)+c*v(3,iq)
            v(1,ip)=v1p
            v(2,ip)=v2p
            v(3,ip)=v3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      ! Sort eigenvalues S1 <= S2 <= S3
      is1 = 1
      is2 = 2
      is3 = 3
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
        it  = is3
        is3 = is2
        is2 = it
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      Do i=1,3
        xN1(i) = v(i,is1) ! first  column
        xN2(i) = v(i,is2) ! second column
        xN3(i) = v(i,is3) ! third  column
      End Do
      Return
      End ! Eig_3

      Subroutine Eig_3a(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3)   !  V(3,3),xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      If (iOpt.Eq.1) Tol = 1d-50*abs_max_s
      it=0
      itmax = 50
      Do While ( it.lt.itmax .And.
     *           abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )

        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
          If (a(ip,iq) .Ne. 0.0) Then         ! ongelijk nul ?
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      Return
      End ! Eig_3a

C
C***********************************************************************
      Logical Function LEqual(A,B,Eps)
C***********************************************************************
C
C     Returns .TRUE.  when two real values are (almost) equal,
C             .FALSE. otherwise
C
C I   A,B  : Two real values to be compared
C I   Eps  : Toleration (Magnitude ~= 1E-5)
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
C***********************************************************************
      LEqual =.True.
      If (A .Eq. B) Return
      If (DAbs(A-B) .LT. 0.5D0*Eps*( DAbs(A) + DAbs(B) + Eps ) )Return
      LEqual =.False.
      Return
      End     ! function LEqual
C
C***********************************************************************
      Subroutine CrossProd(xN1,xN2,xN3)
C***********************************************************************
C
C     Returns cross product of xN1 and xN2
C
C I   xN1,xN2 : Two basic vectors
C O   xN3     : Resulting vector
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*)
C***********************************************************************

      xN3(1) = xN1(2)*xN2(3) - xN1(3)*xN2(2)
      xN3(2) = xN1(3)*xN2(1) - xN1(1)*xN2(3)
      xN3(3) = xN1(1)*xN2(2) - xN1(2)*xN2(1)

      Return
      End     ! Subroutine CrossProd
C
C***********************************************************************
      Double Precision Function ArcSin(X,ie)
C***********************************************************************
C
C     Returns the Arc Sine of X
C
C I   X : Input value
C
C     Note : In stead of using default routine DASIN we use this one
C            because Â³XÂ³ can be slightly beyond 1 and this will give
C            a RTE using DASIN(X)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
C***********************************************************************
      Ie=0
      S = (1-X*X)
!      If (S .Lt. -1E-10) Ie=1
!      If (S .Lt. -1E-10) Write(*,1) X,S
!      If (S .Lt. -1E-10) Write(2,1) X,S
    1 Format(' ArcSin(',1x,1p,e13.5e3,') , S =',1x,1p,e13.5e3)
      If (S.LT.0) S = 0
      S = DSQRT(S)
      ArcSin = DATan2(X,S)
      Return
      End     ! function ArcSin
C
C***********************************************************************
      Subroutine CarSig(S1,S2,S3,xN1,xN2,xN3,SNew)
C***********************************************************************
C
C     Returns the Cartesian stresses using the principal stresses S1..S3
C     and the principal directions
C
C I   S1..S3   : Principal stresses
C I   xN1..xN3 : Principal directions (xNi for Si)
C
C***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*),SNew(*)
      Dimension SM(3,3),T(3,3),TT(3,3),STT(3,3)
C***********************************************************************
C
C**** Fill transformation (rotation) matrix
C
      Do I=1,3
        T(I,1) = xN1(I)
        T(I,2) = xN2(I)
        T(I,3) = xN3(I)
        TT(1,I) = T(I,1)
        TT(2,I) = T(I,2)
        TT(3,I) = T(I,3)
      End Do
!      Call MatTranspose(T,3,TT,3,3,3)

      Call MZeroR(SM,9)
      SM(1,1) = S1
      SM(2,2) = S2
      SM(3,3) = S3
C
C**** SMnew = T*SM*TT
C
      Call MatMat(SM ,3,  TT,3 , 3,3,3 ,STT,3)
      Call MatMat( T ,3, STT,3 , 3,3,3 ,SM ,3)
!     Call MatMatSq(3, SM,  TT, STT )   ! STT = SM*TT
!     Call MatMatSq(3,  T, STT, SM  )   ! SM  =  T*STT
C
C**** Extract cartesian stress vector from stress matrix
C
      Do I=1,3
        SNew(I) = SM(I,I)
      End Do
      SNew(4) = SM(2,1)
      SNew(5) = SM(3,2)
      SNew(6) = SM(3,1)

      Return
      End     ! Subroutine CarSig
C**********************************************************************
      subroutine setveclen(xn,n,xl)
C**********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN(*)
      x=0
      do i=1,n
        x=x+xn(i)**2
      end do
      if (x.Ne.0) Then
        f=xl/sqrt(x)
        do i=1,3
          xn(i)=xn(i)*f
        end do
      end if
      return
      end ! setveclen

C**********************************************************************
C End Of file
C**********************************************************************