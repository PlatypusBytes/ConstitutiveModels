module MathUtils
    
    use, intrinsic :: iso_c_binding
    implicit None
    
    type :: Utils
    contains
            procedure :: vector_vector_dot_product
            procedure :: matrix_vector_product
            procedure :: matrix_matrix_product
            procedure :: vector_outer_product
            procedure :: calculate_stress_invariants
            procedure :: calculate_principal_stresses
            procedure :: set_3x3_stress_matrix
            procedure :: identity
            procedure :: determinant3x3
            procedure :: determinant2x2
            procedure :: eigen3x3
            procedure :: rotate_matrix
            procedure :: calculate_derivative_invariants_sigma
    end type Utils
        
        
    contains
    
    
    
    function vector_vector_dot_product(this, x, y) result(v)
        implicit double precision(a-h,o-z)
        class(Utils), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: x(:), y(:)
        real (C_DOUBLE) :: v
        integer :: i
        
        v = 0.0
        do i = 1, size(x)
            v = v + x(i) * y(i)
        end do
    
    end function vector_vector_dot_product
        
        
    
    ! matrix vector product
    function matrix_vector_product(this, A, x) result(y)
        implicit double precision(a-h,o-z)
        class(Utils), intent(inout) :: this
    
        real(C_DOUBLE), intent(in) :: A(:,:)
        real(C_DOUBLE), intent(in) :: x(:)
        real(C_DOUBLE) :: y(size(x))
        integer :: i, j

        do i = 1, size(x)
            y(i) = 0.0
            do j = 1, size(A, 2) ! Iterate over the columns of A
                y(i) = y(i) + A(i, j) * x(j)
            end do
        end do
    end function matrix_vector_product
    
    
    function matrix_matrix_product(this, A, B) result(C)
        implicit double precision(a-h,o-z)
        class(Utils), intent(inout) :: this
    
        real(C_DOUBLE), intent(in) :: A(:,:)
        real(C_DOUBLE), intent(in) :: B(:,:)
        real(C_DOUBLE) :: C(size(A, 1), size(B, 2))
        integer :: i, j, k

        do i = 1, size(A, 1)
            do j = 1, size(B, 2)
                C(i, j) = 0.0
                do k = 1, size(A, 2)
                    C(i, j) = C(i, j) + A(i, k) * B(k, j)
                end do
            end do
        end do
    end function matrix_matrix_product
    
    
    ! calculate the stress invariants s, t, theta
    function calculate_stress_invariants(this, stress_vector, NTENS) result(I)
        implicit double precision(a-h,o-z)
        class(Utils), intent(inout) :: this
        
        integer(C_INT), intent(in) :: NTENS
        real(C_DOUBLE), intent(in) :: stress_vector(NTENS)
        real(C_DOUBLE) :: I(3)
        real(C_DOUBLE) :: s, sx, sy, sz,t, J3
        
        ! s invariant, distance origin to pi plane
        s = 1/(3**(1/2)) *(stress_vector(1) + stress_vector(2) + stress_vector(3))
        
        sx = (2*stress_vector(1) - stress_vector(2) - stress_vector(3))/3
        sy = (2*stress_vector(2) - stress_vector(1) - stress_vector(3))/3
        sz = (2*stress_vector(3) - stress_vector(1) - stress_vector(2))/3
        
        ! mean stress        
        I(1) = 1/(3**(1/2)) * s
        
        ! t invariant, perpendicular distance stress point from space diagonal
        if (NTENS == 4) then
            ! 2D case
            t = 1/(3**(1/2)) *((stress_vector(1) - stress_vector(2))**2 + (stress_vector(2) - stress_vector(3))**2 + (stress_vector(3) - stress_vector(1))**2 &
                +6*stress_vector(4)**2 )**(1/2)       
            
            J3 = sx*sy*sz - sz*stress_vector(4)**2 
        else if (NTENS == 6) then
            ! 3D case
            t = 1/(3**(1/2)) *((stress_vector(1) - stress_vector(2))**2 + (stress_vector(2) - stress_vector(3))**2 + (stress_vector(3) - stress_vector(1))**2 &
                +6*stress_vector(4)**2 +6*stress_vector(5)**2 + +6*stress_vector(6)**2)**(1/2)     
            
        
            J3 = sx*sy*sz - sx*stress_vector(4)**2 - sy*stress_vector(5)**2 - sz*stress_vector(6)**2 + 2*stress_vector(4)*stress_vector(5)*stress_vector(6)
        else 
            print *, "ERROR: NTENS must be 4 or 6"
            
        end if
        
        ! deviatoric stress
        I(2) = (3/2)**(1/2) * t
        
        
        ! theta invariant, lode angle, measure of angular position of the stress point within the pi plane
        I(3) = (1/3) * asin(-3* 6**(1/2) * J3/(t**3))
        
    end function calculate_stress_invariants
    
    function calculate_principal_stresses(this, invariants) result(principal_stresses)
        implicit double precision(a-h,o-z)
        class(Utils), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: invariants(3)
        real(C_DOUBLE) :: principal_stresses(3)
        real(C_DOUBLE) :: pi
        
        ! define pi like this for maximum precision
        pi = 4*atan(1.0)
        
        principal_stresses(1) = invariants(1) + 2/3*invariants(2) * sin(invariants(3) - 2*pi/3)
        principal_stresses(2) = invariants(1) + 2/3*invariants(2) * sin(invariants(3))
        principal_stresses(3) = invariants(1) + 2/3*invariants(2) * sin(invariants(3) + 2*pi/3)
                
        
    end function calculate_principal_stresses
    
    
    function vector_outer_product(this, x, y) result(A)
        implicit double precision(a-h,o-z)
        class(Utils), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: x(:)
        real(C_DOUBLE), intent(in) :: y(:)
        real(C_DOUBLE) :: A(size(x), size(y))
        integer :: i, j
        
        A = 0.0
        do i = 1, size(x)
            do j = 1, size(y)
                A(i, j) = x(i) * y(j)
            end do
        end do
        
    end function vector_outer_product
    
    function identity(this, n) result(A)
        implicit double precision(a-h,o-z)
        class(Utils), intent(inout) :: this
        
        integer(C_INT), intent(in) :: n
        real(C_DOUBLE) :: A(n, n)
        integer :: i, j
        
        A = 0.0
        do i = 1, n
            A(i, i) = 1.0
        end do
        
    end function identity
    
    function set_3x3_stress_matrix(this, stress_vector, NTENS) result(stress_matrix)
        implicit double precision(a-h,o-z)
        class(Utils), intent(inout) :: this
        
        integer(C_INT), intent(in) :: NTENS
        real(C_DOUBLE), intent(in) :: stress_vector(NTENS)
        real(C_DOUBLE) :: stress_matrix(3, 3)
        
        stress_matrix(1, 1) = stress_vector(1)
        stress_matrix(2, 2) = stress_vector(2)
        stress_matrix(3, 3) = stress_vector(3)
        
        stress_matrix(1, 2) = stress_vector(4)
        stress_matrix(2, 1) = stress_vector(4)
        
        ! 3D case todo check ordering
        if (NTENS .eq. 6) then
            stress_matrix(1, 3) = stress_vector(5)
            stress_matrix(3, 1) = stress_vector(5)
            stress_matrix(2, 3) = stress_vector(6)
            stress_matrix(3, 2) = stress_vector(6)
        end if
        
    end function set_3x3_stress_matrix
    
    
    function determinant2x2(this, A) result(det)
            implicit double precision(a-h,o-z)
            class (Utils), intent(inout) :: this
            real(C_DOUBLE), intent(in) :: A(2,2)
            real(C_DOUBLE) :: det
           
            det = A(1,1)*A(2,2) - A(1,2)*A(2,1)
            
        end function determinant2x2
    
    
    function determinant3x3(this, A) result(det)
            implicit double precision(a-h,o-z)
            class (Utils), intent(inout) :: this
            
            real(C_DOUBLE), intent(in) :: A(3,3)
            real(C_DOUBLE) :: det
           
            det = A(1,1) * this%determinant2x2(A(2:3, 2:3)) - A(1,2) * this%determinant2x2(A(2:3, [1,3])) + A(1,3) * this%determinant2x2(A(2:3, [1,2]))
            
            
        end function determinant3x3
        
    
    subroutine eigen3x3(this, A, eigenvalues, rotation_matrix)
        implicit double precision(a-h,o-z)
        class (Utils), intent(inout) :: this
    
        real(C_DOUBLE), intent(in) :: A(3,3)
        real(C_DOUBLE), intent(out) :: eigenvalues(3)
        real(C_DOUBLE), intent(out) :: rotation_matrix(3,3)
    
        real(C_DOUBLE) :: trace, q, p1, p2, precision, pi
        real(C_DOUBLE), dimension(3,3) :: B, identity, tmp1, tmp2, tmp3
        real(C_DOUBLE), dimension(3) :: eigen_vector
        integer(C_INT) :: i
        precision = 1e-12
        pi = 4*atan(1.0)
    
    
        !calculate trace
        trace = 0.0
    
        ! todo generalize
        p1 = A(1,2)**2 + A(1,3)**2 + A(2,3)**2
    
        if (abs(p1) < precision) then
            do i = 1, size(A, 1)
                eigenvalues(i) = A(i,i)
            end do
        else
        
            do i = 1, size(A, 1)
               trace = trace + A(i, i)
            end do
            q = trace/size(A, 1)
        
            p2 = (A(1,1) - q)**2 + (A(2,2) - q)**2 + (A(3,3) - q)**2 + 2*p1
            p = sqrt(p2/6)
        
            identity = this%identity(3)
            B = (1/p)*(A - q*identity)
            r = this%determinant3x3(B)/2
        
            if (r <= -1 + precision) then
                phi = pi/3
            else if (r >= 1 - precision) then
                phi = 0
            else
                phi = acos(r)/3
            end if
        
            ! calculate eigenvalues
            eigenvalues(1) = q + 2*p*cos(phi)
            eigenvalues(3) = q + 2*p*cos(phi + (2*pi/3))
            eigenvalues(2) = 3*q - eigenvalues(1) - eigenvalues(3)
        
        
            ! calculate eigenvectors
            tmp1 = A - eigenvalues(1) * identity
            x1 = - (tmp1(1,2) + tmp1(1,3) ) / tmp1(1,1)
        
            eigen_vector(1) = x1
            eigen_vector(2) = 1
            eigen_vector(3) = 1
        
!            norm_vector_1 = 1/sqrt(x1**2 + 2) * eigen_vector
            ! normalize
            rotation_matrix(:,1) = 1/sqrt(x1**2 + 2) * eigen_vector
            
            
            tmp2 = A - eigenvalues(2) * identity
        
            x2 = - (tmp2(2,1) + tmp2(2,3) ) / tmp2(2,2)
        
            eigen_vector(1) = 1
            eigen_vector(2) = x2
            eigen_vector(3) = 1
            
            ! normalize
            rotation_matrix(:,2)  = 1/sqrt(x2**2 + 2) * eigen_vector
        
            tmp3 = A - eigenvalues(3) * identity
        
            x3 = - (tmp3(3,1) + tmp3(3,2) ) / tmp3(3,3)
        
            eigen_vector(1) = 1
            eigen_vector(2) = 1
            eigen_vector(3) = x3
            
            ! normalize
            rotation_matrix(:,3) = 1/sqrt(x3**2 + 2) * eigen_vector
    
        end if
    end subroutine eigen3x3
    
    function rotate_matrix(this, A, rotation_matrix) result(B)
        implicit double precision(a-h,o-z)
        class (Utils), intent(inout) :: this
        
        real(C_DOUBLE), intent(in) :: A(:,:)
        real(C_DOUBLE), intent(in) :: rotation_matrix(3,3)
        real(C_DOUBLE) :: B(size(A, 1), size(A, 2)), AT(size(A, 2), size(A, 1))
        integer :: i, j, k
        
        
        AT = transpose(A)
        
        B = matmul(matmul(A, rotation_matrix), AT)

        
    end function rotate_matrix
    
    
    subroutine calculate_derivative_invariants_sigma(this, invariants, stress_vector, dsigma_m_dsigma, dsigma_dev_dsigma, dJ3_dsigma)
    
        implicit double precision(a-h,o-z)
        class (Utils), intent(inout) :: this
        real(C_DOUBLE), intent(in) :: invariants(3), stress_vector(6)
        real(C_DOUBLE), intent(out) :: dsigma_m_dsigma(6), dsigma_dev_dsigma(6), dJ3_dsigma(6)
        real(C_DOUBLE) :: sx, sy, sz,  tau_xy, tau_yz, tau_xz, sigma_m, sigma_dev
        
        sigma_m = invariants(1)
        sigma_dev = invariants(2)
        
        sx = stress_vector(1) - sigma_m
        sy = stress_vector(2) - sigma_m
        sz = stress_vector(3) - sigma_m
        
        tau_xy = stress_vector(4)
        tau_yz = stress_vector(5)
        tau_xz = stress_vector(6)
        
        dsigma_m_dsigma = 1/3 * [1,1,1,0,0,0]
        dsigma_dev_dsigma = 1/(2*sigma_dev) * [sx, sy, sz, 2*tau_xy, 2*tau_yz, 2*tau_xz]
        dJ3_dsigma = 1/3 * [sy*sz - tau_yz**2, sx*sz - tau_xz**2, sx*sy - tau_xy**2, 2*(tau_yz*tau_xz - sz*tau_xy), 2*(tau_xz*tau_xy - sx*tau_yz), 2*(tau_xy*tau_yz - sy*tau_xz)] + &
                            sigma_dev**2/3 * [1,1,1,0,0,0]
        
    end subroutine calculate_derivative_invariants_sigma

end Module MathUtils

        
    