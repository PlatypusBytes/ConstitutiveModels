
    module MatsuokaNakai
    
    use, intrinsic :: iso_c_binding
    use MathUtils
    
    implicit None
    
    type :: YieldSurface
    contains
            procedure :: get_yield_criteria
            procedure :: get_yield_criteria_princ
            procedure :: get_derivatives_princ
            procedure :: get_invariant_yield_criteria
            procedure :: get_derivatives
    end type YieldSurface
        
        
    contains
    
    function get_invariant_yield_criteria(this, invariants, c, phi) result(yield_criteria)
        implicit double precision(a-h,o-z)
        class(YieldSurface), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: invariants(3)
        real(C_DOUBLE) :: yield_criteria
        
        real(C_DOUBLE) :: sigma_m, sigma_dev, theta, c, phi, c_cos_phi, phi_rad
        
        
    
    
    
    
    module MohrCoulomb
    
    use, intrinsic :: iso_c_binding
    use MathUtils
    
    implicit None
    
    type :: YieldSurface
    contains
            procedure :: get_yield_criteria
            procedure :: get_yield_criteria_princ
            procedure :: get_derivatives_princ
            procedure :: get_invariant_yield_criteria
            procedure :: get_derivatives
    end type YieldSurface
        
        
    contains
    
    
    function get_yield_criteria(this, sigma1, sigma2, sigma3, phi, c) result(yield_criteria)
        implicit double precision(a-h,o-z)
        class(YieldSurface), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: sigma1, sigma2, sigma3, phi, c
        real(C_DOUBLE) :: yield_criteria(6)
        
        real(C_DOUBLE) :: c_cos_phi, phi_rad
        
        phi_rad = phi * 4 * atan(1.0) / 180.0
        c_cos_phi = c * cos(phi_rad)
        
        yield_criteria(1) = 1/2 *(sigma2-sigma3) + 1/2 *(sigma2-sigma3)*sin(phi_rad) - c_cos_phi
        yield_criteria(2) = 1/2 *(sigma3-sigma2) + 1/2 *(sigma3-sigma2)*sin(phi_rad) - c_cos_phi
        yield_criteria(3) = 1/2 *(sigma3-sigma1) + 1/2 *(sigma3-sigma1)*sin(phi_rad) - c_cos_phi
        yield_criteria(4) = 1/2 *(sigma1-sigma3) + 1/2 *(sigma1-sigma3)*sin(phi_rad) - c_cos_phi
        yield_criteria(5) = 1/2 *(sigma1-sigma2) + 1/2 *(sigma1-sigma2)*sin(phi_rad) - c_cos_phi
        yield_criteria(6) = 1/2 *(sigma2-sigma1) + 1/2 *(sigma2-sigma1)*sin(phi_rad) - c_cos_phi
        
        
    end function get_yield_criteria
    
    
    function get_yield_criteria_princ(this, principal_stresses, phi, c) result(yield_criteria)
        implicit double precision(a-h,o-z)
        class(YieldSurface), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: principal_stresses(3)
        real(C_DOUBLE), intent(in) :: phi, c
        real(C_DOUBLE) :: yield_criteria
        
        real(C_DOUBLE) :: c_cos_phi, pi, phi_rad
        !real(C_DOUBLE) :: tollerance
        
        phi_rad = phi * 4 * atan(1.0) / 180.0
        
        !print *,"phi_rad", phi_rad
        !tollerance = 1e-10
        c_cos_phi = c * cos(phi_rad )
        
        !print *,"princ_stress", principal_stresses(1), principal_stresses(2), principal_stresses(3)
        yield_criteria = 1/2 *(principal_stresses(1)-principal_stresses(3)) + 1/2 *(principal_stresses(1)+principal_stresses(3))*sin(phi_rad) - c_cos_phi
        
        !print *,"1 -3 ", principal_stresses(1)-principal_stresses(3)
        !print *,"1 +3 ", principal_stresses(1)+principal_stresses(3)
        
        !print *,"yield_criteria", yield_criteria
        !
        !if (abs(sigma2-sigma3) < tollerance) then
        !    yield_criteria = 1/2 *(sigma1-sigma2) + 1/2 *(sigma1-sigma2)*sin(phi) - c_cos_phi
        !els
        !    yield_criteria = 1/2 *(sigma2-sigma3) + 1/2 *(sigma2-sigma3)*sin(phi) - c_cos_phi
        !end if
        
        
        !yield_criteria(1) = 1/2 *(sigma2-sigma3) + 1/2 *(sigma2+sigma3)*sin(phi) - c_cos_phi
        !yield_criteria(2) = 1/2 *(sigma3-sigma2) + 1/2 *(sigma3+sigma2)*sin(phi) - c_cos_phi
        !yield_criteria(3) = 1/2 *(sigma3-sigma1) + 1/2 *(sigma3+sigma1)*sin(phi) - c_cos_phi
        !yield_criteria(4) = 1/2 *(sigma1-sigma3) + 1/2 *(sigma1+sigma3)*sin(phi) - c_cos_phi
        !yield_criteria(5) = 1/2 *(sigma1-sigma2) + 1/2 *(sigma1+sigma2)*sin(phi) - c_cos_phi
        !yield_criteria(6) = 1/2 *(sigma2-sigma1) + 1/2 *(sigma2+sigma1)*sin(phi) - c_cos_phi
        
        
    end function get_yield_criteria_princ
    
    function get_derivatives_princ(this, phi, ntens) result(dyield)
        implicit double precision(a-h,o-z)
        class(YieldSurface), intent(inout) :: this
        real(C_DOUBLE), intent(in) :: phi
        integer(C_INT), intent(in) :: ntens
        
        real(C_DOUBLE) :: dyield(ntens), phi_rad
        
        phi_rad = phi * 4 * atan(1.0) / 180.0
        ! initialize dyield as zero
        dyield = 0
        
        ! derivatives with respect to sigma1
        dyield(1) = 1/2 + 1/2 * sin(phi_rad)
        
        ! derivatives with respect to sigma2
        dyield(2) = -1/2 + 1/2 * sin(phi_rad)
        
        ! derivatives with respect to sigma3
        dyield(3) = -1/2 + 1/2 * sin(phi_rad)
  

    end function get_derivatives_princ
    
    
    function get_derivatives(this, invariants, stress_vector, phi) result(derivatives)
        implicit double precision(a-h,o-z)
        class(YieldSurface), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: invariants(3), stress_vector(6)
        real(C_DOUBLE), intent(in) :: phi
        real(C_DOUBLE) :: derivatives(6)
        real(C_DOUBLE) :: deviatoric_stresses(6), dsigma_m_dsigma(6), dsigma_dev_dsigma(6), dJ3_dsigma(6)
        real(C_DOUBLE) :: sigma_m, sigma_dev, theta, phi_rad
        integer(C_INT) :: i
        
        type(Utils) :: math_utils

        ! chainrule: df/dsigma = df/dsigma_m * dsigma_m/dsigma + df/dsigma_dev * dsigma_dev/dsigma + df/dtheta * dtheta/dsigma
        ! mohr coulomb criteria with invariants = sigma_m* sin(phi) + sigma_dev * (cos(theta)/(3**(1/2)) - sin(theta)*sin(phi)/3) - c*cos(phi)
        ! df/dsigma = C1 * dsigma_m/dsigma + C2 * dsigma_dev/dsigma + C3 * dJ3/dsigma
        ! where: C1 = dF/ dsigma_m, C2 = dF/dsigma_dev - tan(3*theta)/sigma_dev * dF/dTheta, C3 = - sqrt(3)/(2*sigma_dev**3cos(3*theta)*dF/dTheta
        
        phi_rad = phi * 4 * atan(1.0) / 180.0
        
        sigma_dev = invariants(2)
        theta = invariants(3)
        
        ! calculate dF/dInvariants
        dF_dsigma_m = sin(phi_rad)
        dF_dsigma_dev = cos(theta)/(3**(1/2)) - sin(theta)*sin(phi_rad)/3
        dF_dtheta = sigma_dev * (-sin(theta)/(3**(1/2)) - cos(theta)*sin(phi_rad)/3)
        
        ! calculate dInvariants/dsigma
        call math_utils%calculate_derivative_invariants_sigma(invariants, stress_vector, dsigma_m_dsigma, dsigma_dev_dsigma, dJ3_dsigma)

        ! calculate Df/dsigma
        derivatives = dF_dsigma_m * dsigma_m_dsigma + (dF_dsigma_dev - tan(3*theta)/sigma_dev * dF_dtheta) * dsigma_dev_dsigma + &
                        - 3**(1/2)/(2*sigma_dev**3*cos(3*theta)) * dF_dtheta * dJ3_dsigma
        
        end function get_derivatives
    
    function get_invariant_yield_criteria(this, invariants, c, phi) result(yield_criteria)
        implicit double precision(a-h,o-z)
        class(YieldSurface), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: invariants(3)
        real(C_DOUBLE) :: yield_criteria
        
        real(C_DOUBLE) :: sigma_m, sigma_dev, theta, c, phi, c_cos_phi, phi_rad
        
        
        sigma_m = invariants(1)
        sigma_dev = invariants(2)
        theta = invariants(3)
        
        phi_rad = phi * 4 * atan(1.0) / 180.0
        
        c_cos_phi = c * cos(phi_rad)
                       
        yield_criteria = sigma_m * sin(phi_rad) + sigma_dev * (cos(theta)/(3**(1/2)) - sin(theta)*sin(phi_rad)/3) - c_cos_phi
        
        
    end function get_invariant_yield_criteria
    


end Module MohrCoulomb

        
    