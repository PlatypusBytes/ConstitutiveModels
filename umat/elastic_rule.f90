module HookesLaw
    use, intrinsic :: iso_c_binding
    implicit double precision(a-h,o-z)
    
    type :: ElasticData
    contains 
        procedure :: fill_elastic_matrix
    end type ElasticData
      
    contains
      
        function fill_elastic_matrix(this, E, nu, ntens) result(elastic_matrix)
            class(ElasticData), intent(inout) :: this
            real(C_DOUBLE), intent(in) :: E, nu
            integer(C_INT), intent(in) :: ntens
            real(C_DOUBLE) :: elastic_matrix(ntens,ntens)
        
            real(C_DOUBLE) :: v1,v2,v3
            
            ! initialize elastic matrix to zero
            elastic_matrix = 0.0
            
            v1 = E*(1-nu)/((1+nu)*(1-2*nu))
            v2 = E*nu/((1+nu)*(1-2*nu))
            v3 = E/(2*(1+nu))
            
            ! 3d matrix
            if (ntens .eq. 6) then
            
                elastic_matrix(1,1) = v1
                elastic_matrix(1,2) = v2
                elastic_matrix(1,3) = v2
                elastic_matrix(2,1) = v2
                elastic_matrix(2,2) = v1
                elastic_matrix(2,3) = v2
                elastic_matrix(3,1) = v2
                elastic_matrix(3,2) = v2
                elastic_matrix(3,3) = v1
                elastic_matrix(4,4) = v3
                elastic_matrix(5,5) = v3
                elastic_matrix(6,6) = v3
            !    
            ! 2d matrix
            else if (ntens .eq. 4) then
                elastic_matrix(1,1) = v1
                elastic_matrix(1,2) = v2
                elastic_matrix(1,3) = v2
                elastic_matrix(2,1) = v2
                elastic_matrix(2,2) = v1
                elastic_matrix(2,3) = v2
                elastic_matrix(3,1) = v2
                elastic_matrix(3,2) = v2
                elastic_matrix(3,3) = v1
                elastic_matrix(4,4) = v3
            else
                print *, "Error: ntens must be 4 or 6"
            end if
            
        end function fill_elastic_matrix       
        
end module HookesLaw

      
      
      
      

