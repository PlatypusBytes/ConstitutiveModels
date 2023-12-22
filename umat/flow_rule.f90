module MohrCoulombFlowRule
    
    use, intrinsic :: iso_c_binding
    use MathUtils
    implicit None
    
    type :: FlowRule
    contains
            procedure :: get_flow_rule_princ
            procedure :: get_derivatives_princ
            procedure :: get_derivatives
    end type FlowRule
        
        
    contains
    
    
    !function get_yield_criteria(this, sigma1, sigma2, sigma3, phi, c) result(yield_criteria)
    !    implicit double precision(a-h,o-z)
    !    class(YieldSurface), intent(inout) :: this
    !    
    !    real(C_DOUBLE), intent(in) :: sigma1, sigma2, sigma3, phi, c
    !    real(C_DOUBLE) :: yield_criteria(6)
    !    
    !    real(C_DOUBLE) :: c_cos_phi
    !    
    !    c_cos_phi = c * cos(phi)
    !    
    !    yield_criteria(1) = 1/2 *(sigma2-sigma3) + 1/2 *(sigma2-sigma3)*sin(phi) - c_cos_phi
    !    yield_criteria(2) = 1/2 *(sigma3-sigma2) + 1/2 *(sigma3-sigma2)*sin(phi) - c_cos_phi
    !    yield_criteria(3) = 1/2 *(sigma3-sigma1) + 1/2 *(sigma3-sigma1)*sin(phi) - c_cos_phi
    !    yield_criteria(4) = 1/2 *(sigma1-sigma3) + 1/2 *(sigma1-sigma3)*sin(phi) - c_cos_phi
    !    yield_criteria(5) = 1/2 *(sigma1-sigma2) + 1/2 *(sigma1-sigma2)*sin(phi) - c_cos_phi
    !    yield_criteria(6) = 1/2 *(sigma2-sigma1) + 1/2 *(sigma2-sigma1)*sin(phi) - c_cos_phi
    !    
    !    
    !end function get_yield_criteria
    
    
    function get_flow_rule_princ(this, sigma1, sigma2, sigma3, psi, c) result(flow_rule)
        implicit double precision(a-h,o-z)
        class(FlowRule), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: sigma1, sigma2, sigma3, psi, c
        real(C_DOUBLE) :: flow_rule
        
        
        flow_rule = 1/2 *(sigma1-sigma3) + 1/2 *(sigma1+sigma3)*sin(psi)
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
        
        
    end function get_flow_rule_princ
    
    function get_derivatives_princ(this, psi, ntens) result(dflow)
        implicit double precision(a-h,o-z)
        class(FlowRule), intent(inout) :: this
        real(C_DOUBLE), intent(in) :: psi
        integer(C_INT), intent(in) :: ntens
        real(C_DOUBLE) :: dflow(ntens)
        
        ! initialize dflow as zero
        dflow = 0
        
        ! derivatives with respect to sigma1
        dflow(1) = 1/2 + 1/2 * sin(psi)
        
        ! derivatives with respect to sigma2
        dflow(2) = -1/2 + 1/2 * sin(psi)
        
        ! derivatives with respect to sigma3
        dflow(3) = -1/2 + 1/2 * sin(psi)
  

    end function get_derivatives_princ
    
        

    function get_derivatives(this, invariants, stress_vector, psi) result(derivatives)
        implicit double precision(a-h,o-z)
        class(FlowRule), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: invariants(3), stress_vector(6)
        real(C_DOUBLE), intent(in) :: psi
        real(C_DOUBLE) :: derivatives(6)
        real(C_DOUBLE) :: deviatoric_stresses(6), dsigma_m_dsigma(6), dsigma_dev_dsigma(6), dJ3_dsigma(6)
        real(C_DOUBLE) :: sigma_dev, theta, psi_rad
        integer(C_INT) :: i
        
        type(Utils) :: math_utils

        ! chainrule: df/dsigma = df/dsigma_m * dsigma_m/dsigma + df/dsigma_dev * dsigma_dev/dsigma + df/dtheta * dtheta/dsigma
        ! mohr coulomb criteria with invariants = sigma_m* sin(phi) + sigma_dev * (cos(theta)/(3**(1/2)) - sin(theta)*sin(phi)/3) - c*cos(phi)
        ! df/dsigma = C1 * dsigma_m/dsigma + C2 * dsigma_dev/dsigma + C3 * dJ3/dsigma
        ! where: C1 = dF/ dsigma_m, C2 = dF/dsigma_dev - tan(3*theta)/sigma_dev * dF/dTheta, C3 = - sqrt(3)/(2*sigma_dev**3cos(3*theta)*dF/dTheta
        
        psi_rad = psi * 4*atan(1.0)/180
        
        sigma_dev = invariants(2)
        theta = invariants(3)
        
        ! calculate dF/dInvariants
        dG_dsigma_m = sin(psi_rad)
        dG_dsigma_dev = cos(theta)/(3**(1/2)) - sin(theta)*sin(psi_rad)/3
        dG_dtheta = sigma_dev * (-sin(theta)/(3**(1/2)) - cos(theta)*sin(psi_rad)/3)
        
        ! calculate dInvariants/dsigma
        call math_utils%calculate_derivative_invariants_sigma(invariants, stress_vector, dsigma_m_dsigma, dsigma_dev_dsigma, dJ3_dsigma)

        ! calculate Df/dsigma
        derivatives = dF_dsigma_m * dsigma_m_dsigma + (dF_dsigma_dev - tan(3*theta)/sigma_dev * dF_dtheta) * dsigma_dev_dsigma + &
                        - 3**(1/2)/(2*sigma_dev**3*cos(3*theta)) * dF_dtheta * dJ3_dsigma
        
        end function get_derivatives
    

end Module MohrCoulombFlowRule

        
    