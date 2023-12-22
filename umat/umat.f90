SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL, DDSDDT, &
                DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, &
                DTEMP, PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, &
                NSTATEV, PROPS, NPROPS, COORDS, DROT, PNEWDT, CELENT, &
                DFGRD0, DFGRD1, NOEL, NPT, LAYER, KSPT, KSTEP, KINC) 
    
  !DEC$ ATTRIBUTES DLLExport,StdCall,reference :: umat
  use, intrinsic :: iso_c_binding
  use HookesLaw, only: ElasticData
  use MathUtils
  use MohrCoulomb
  use MohrCoulombFlowRule
  
  
  implicit double precision(a-h,o-z)
  parameter (nprecd=2)
        
  CHARACTER(C_CHAR) :: CMNAME
  REAL(C_DOUBLE), DIMENSION(NTENS) :: STRESS, STATEV, DDSDDT, DRPLDE, STRAN, DSTRAN, PREDEF, DPRED, Sig, dSig, princ_stress_vector, dfdsig, dgdsig, sig_elastic
  REAL(C_DOUBLE), DIMENSION(NTENS, NTENS) :: DDSDDE, DFGRD0, DFGRD1, Dep, dgdfT
  REAL(C_DOUBLE), DIMENSION(3) :: COORDS, DROT, TIME, invariants, principal_stresses
  REAL(C_DOUBLE), DIMENSION(NPROPS) :: PROPS
  REAL(C_DOUBLE) :: SSE, SPD, SCD, RPL, DTIME, TEMP, DTEMP, PNEWDT
  INTEGER(C_INT) :: NDI, NSHR, NTENS, NSTATEV, NPROPS, CELENT, NOEL, NPT, LAYER, KSPT, KSTEP, KINC
  REAL(C_DOUBLE) :: G, ENU, one, two, FAC, D1, D2, DSTRANVOL, Dep_denominator, yield_criteria, lambda
  INTEGER(C_INT) :: i
  
  REAL(C_DOUBLE), DIMENSION(NTENS, NTENS) :: elastic_matrix
  real(C_DOUBLE), DIMENSION(3, 3) :: stress_matrix, rotation_matrix
  real(C_DOUBLE), DIMENSION(6) :: full_stress_vector
        
  !REAL(C_DOUBLE), DIMENSION(3,3) :: A
  !real(C_DOUBLE), DIMENSION(3) :: wr, wi
  !real(C_DOUBLE), DIMENSION(3,3) :: vl
  !real(C_DOUBLE), DIMENSION(3,3) :: vr
  !integer(C_INT) :: info, lwork
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
  !  3 : phi
  !  4 : c
  !  5 : psi
        
                    
  type(ElasticData) :: elastic_data
  type(Utils) :: math_utils
  type(YieldSurface) :: yield_surface
  type(FlowRule) :: flow_rule
  !        
  elastic_matrix =  elastic_data%fill_elastic_matrix(PROPS(1), PROPS(2), NTENS)
   
  dSig = matmul(elastic_matrix, DSTRAN)
  
  ! elastic stress
  sig_elastic = STRESS + dSig 
  
  full_stress_vector = 0
  
  do i = 1,NTENS
      full_stress_vector(i) = sig_elastic(i)
  end do
  
  !stress_matrix = math_utils%set_3x3_stress_matrix(Sig, NTENS)
  
  !print *, "STRESS" , STRESS
  
  !call math_utils%eigen3x3(stress_matrix, principal_stresses, rotation_matrix)
  
  
  
  !invariants = math_utils%calculate_stress_invariants(STRESS, NTENS)
  !principal_stresses = math_utils%calculate_principal_stresses(invariants)
  !print *, "STRESS" , STRESS
  !print *, "principal_stresses"  , principal_stresses
  !
  !  
  !princ_stress_vector = 0.0
  !do i = 1,3
  !    princ_stress_vector(i) = principal_stresses(i)
  !end do
  
  invariants = math_utils%calculate_stress_invariants(sig_elastic, NTENS)
  
  !print *, "principal_stresses" , principal_stresses
  yield_criteria = yield_surface%get_invariant_yield_criteria(invariants, PROPS(3),PROPS(4))
  
  if (yield_criteria > 0.0) then
      print *, "yield_criteria" , yield_criteria
      
  else if (yield_criteria < 0.0) then
      print *, "yield_criteria" , yield_criteria
  end if

  if (yield_criteria > 0.0) then
      
  
      dfdsig = yield_surface%get_derivatives(invariants, full_stress_vector, PROPS(3))
      dgdsig = flow_rule%get_derivatives(invariants, full_stress_vector, PROPS(5))
  
  
      !dgdfT = math_utils%vector_outer_product(dgdsig, dfdsig)
      !print *, "dgdfT" , dgdfT
  
      ! todo add hardening parameter
      !Dep_denominator = math_utils%vector_vector_dot_product(dfdsig,math_utils%matrix_vector_product(elastic_matrix, dgdsig)
      Dep_denominator = dot_product(dfdsig,matmul(elastic_matrix, dgdsig))
      !Dep_denominator = ddot(NTENS,dfdsig,1,matmul(elastic_matrix, dgdsig),1)
  
      !Dep = elastic_matrix - math_utils%matrix_matrix_product(math_utils%matrix_matrix_product(elastic_matrix, dgdfT), elastic_matrix) / Dep_denominator
      !Dep = elastic_matrix - matmul(matmul(elastic_matrix, dgdfT), elastic_matrix) / Dep_denominator
      ! Dep = De * dgdsig * dfdsigT *De / ((dfdsigT * De * dgdsig) + H)
  
      ! lambda = dfdsigT dot sigma_el / ( dfdsigT dot De dot dgdsig), elastic potential
      lambda = dot_product(dfdsig, Sig) / Dep_denominator
      
      ! total stress
      Sig = full_stress_vector - lambda * matmul(elastic_matrix, dgdsig)
      
      ! calculate dSigmaEl = D * dEpsilon
      !dSig =  math_utils%matrix_vector_product(Dep, DSTRAN)
      !dSig = matmul(Dep, DSTRAN)
      
      
      DDSDDE = Dep 
  else
      
      DDSDDE = elastic_matrix 
  end if
        
  
  !math_utils%rotate_matrix(
  

        
  ! stress state parameters update
  DO i = 1, NTENS
      STRESS(i) = Sig(i)
  END DO
  
  
  
  
  
        
  RETURN
      
END SUBROUTINE UMAT
