import numpy as np
import cffi

class Utils:

    @staticmethod
    def run_c_umat(umat_file_name, stress, state_vars, strain, dstrain, props, ndi=3, nshr=3, nstatv=1):
        """
        Call the UMAT function from Python using CFFI.

        Temperature and predefined field variables are not used in this function.

        Parameters:
        - umat_file_name: Path to the compiled UMAT shared library (DLL)
        - stress: Initial stress vector [σ11, σ22, σ33, σ12, σ13, σ23]
        - state_vars: State variables array
        - strain: Total strain vector
        - dstrain: Increment in strain vector
        - props: Material properties array
        - ndi: Number of direct stress components
        - nshr: Number of shear stress components
        - nstatv: Number of state variables

        Returns:
        - stress: Updated stress vector
        - ddsdde: Updated Jacobian matrix
        - state_vars: Updated state variables
        """

        # Convert inputs to lists if they are not already
        stress = list(stress)
        state_vars = list(state_vars)
        strain = list(strain)
        dstrain = list(dstrain)
        props = list(props)

        # Initialize CFFI
        ffi = cffi.FFI()

        # Define the C function signature
        ffi.cdef("""
        void umat(double* STRESS, double* STATEV, double* DDSDDE,
                  double* SSE, double* SPD, double* SCD, double* RPL,
                  double* DDSDDT, double* DRPLDE, double* DRPLDT,
                  double* STRAN, double* DSTRAN, double* TIME, double* DTIME,
                  double* TEMP, double* DTEMP, double* PREDEF, double* DPRED,
                  char* CMNAME, int* NDI, int* NSHR, int* NTENS, int* NSTATV,
                  double* PROPS, int* NPROPS, double* COORDS, double* DROT,
                  double* PNEWDT, double* CELENT, double* DFGRD0, double* DFGRD1,
                  int* NOEL, int* NPT, int* LAYER, int* KSPT, int* KSTEP, int* KINC);
        """)

        # Load the library
        # Adjust the path as needed
        lib = ffi.dlopen(umat_file_name)  # Windows

        # number of stress components
        ntens = ndi + nshr

        # Convert numpy arrays to C arrays via CFFI
        c_stress = ffi.new("double[]", stress) # stress vector
        c_statev = ffi.new("double[]", state_vars) # state variables
        c_ddsdde = ffi.new("double[]", ntens * ntens)  # Flat array of the Jacobian matrix

        c_sse = ffi.new("double*", 0.0) # elastic strain energy density
        c_spd = ffi.new("double*", 0.0) # plastic dissipation
        c_scd = ffi.new("double*", 0.0) # creep dissipation
        c_rpl = ffi.new("double*", 0.0) # volumetric heat generation

        c_strain = ffi.new("double[]", strain)
        c_dstrain = ffi.new("double[]", dstrain)

        c_time = ffi.new("double[]", [0.0, 0.0])  # [step time, total time]
        c_dtime = ffi.new("double*", 1.0)

        # Material name (80 characters)
        c_cmname = ffi.new("char[]", b"MY_MATERIAL" + b" " * (80 - len(b"MY_MATERIAL")))

        # Integer parameters
        c_ndi = ffi.new("int*", ndi)
        c_nshr = ffi.new("int*", nshr)
        c_ntens = ffi.new("int*", ntens)
        c_nstatv = ffi.new("int*", nstatv)

        c_props = ffi.new("double[]", props)
        c_nprops = ffi.new("int*", len(props))

        # Rotation matrix (identity)
        c_drot = ffi.new("double[]", [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])

        # Element info
        c_noel = ffi.new("int*", 1)
        c_npt = ffi.new("int*", 1)
        c_layer = ffi.new("int*", 1)
        c_kspt = ffi.new("int*", 1)
        c_kstep = ffi.new("int*", 1)
        c_kinc = ffi.new("int*", 1)

        # Call UMAT function without temperature and predefined field variables
        lib.umat(
            c_stress, c_statev, c_ddsdde,
            c_sse, c_spd, c_scd, c_rpl,
            ffi.NULL, ffi.NULL, ffi.NULL,  # DDSDDT, DRPLDE, DRPLDT (null)
            c_strain, c_dstrain, c_time, c_dtime,
            ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL,  # TEMP, DTEMP, PREDEF, DPRED (null)
            c_cmname, c_ndi, c_nshr, c_ntens, c_nstatv,
            c_props, c_nprops,
            ffi.NULL,  # COORDS (null)
            c_drot,
            ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL,  # PNEWDT, CELENT, DFGRD0, DFGRD1 (null)
            c_noel, c_npt, c_layer, c_kspt, c_kstep, c_kinc
        )

        # Convert C arrays back to numpy arrays
        stress_updated = np.array([c_stress[i] for i in range(ntens)])
        statev_updated = np.array([c_statev[i] for i in range(nstatv)])

        # Reshape the flat DDSDDE array to a 2D matrix
        ddsdde_updated = np.array([[c_ddsdde[i * ntens + j] for j in range(ntens)]
                                   for i in range(ntens)])

        return stress_updated, ddsdde_updated, statev_updated

# Example usage
if __name__ == "__main__":
    # Example material parameters
    E = 200000.0  # Young's modulus
    nu = 0.3  # Poisson's ratio
    props = [E, nu, 0.0]

    # Initial conditions
    stress = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Initial stress
    state_vars = [0.0]  # Initial state variable(s)

    # Apply strain increment
    strain = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]  # Initial strain
    dstrain = [-0.001, 0.0, 0.0, 0.0, 0.0, 0.0]  # Strain increment

    # Call UMAT
    updated_stress, jacobian, updated_state_vars = Utils.run_c_umat(r"C:\software_development\ConstitutiveModels\c_models\linear_elastic_with_vertical_tens_cutoff\linear_elastic_with_vertical_tens_cutoff.dll",
        stress, state_vars, strain, dstrain, props
    )

    print("Updated stress:", updated_stress)
    print("Updated state variables:", updated_state_vars)
    print("Jacobian matrix:")
    print(jacobian)

