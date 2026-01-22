import numpy as np
import cffi

class Utils:

    @staticmethod
    def run_fortran_umat(umat_file_name, stress, state_vars, strain, dstrain, props, time_step, ndi=3, nshr=3):
        """
        Call the UMAT function from Python using CFFI.

        Temperature and predefined field variables are not used in this function.

        Args:
            umat_file_name (str): Path to the compiled UMAT shared library (DLL)
            stress (Sequence[float]): Initial stress vector [σ11, σ22, σ33, σ12, σ13, σ23]
            state_vars (Sequence[float]): State variables array
            strain (Sequence[float]): Total strain vector
            dstrain (Sequence[float]): Increment in strain vector
            props (Sequence[float]): Material properties array
            ndi (int): Number of direct stress components (default is 3 for 3D)
            nshr (int): Number of shear stress components (default is 3 for 3D)

        Returns:
            stress (np.ndarray): Updated stress vector
            ddsdde (np.ndarray): Updated Jacobian matrix
            state_vars (np.ndarray): Updated state variables
        """

        import ctypes
        import numpy as np
        import os

        # 2. The EXACT function name inside the DLL (check with dumpbin/nm)
        # Common examples: 'umat_', 'UMAT', '_umat_'
        possible_function_names = ['umat_', 'UMAT', 'umat','_umat_']  # <<< MUST BE CORRECT >>>

        # 3. Define dimensions/sizes based on your UMAT usage (e.g., 3D solid)
        NDI = ndi
        NSHR = nshr
        NTENS = NDI + NSHR  # 6 for 3D solid
        NSTATV = len(state_vars)  # <<< Set to the actual number of state variables your UMAT uses
        NPROPS = len(props)  # <<< Set to the actual number of material props your UMAT expects
        # --- End Configuration ---

        # Check if DLL exists
        if not os.path.exists(umat_file_name):
            raise FileNotFoundError(f"DLL not found at: {umat_file_name}")

        # --- Load the DLL ---
        try:
            umat_lib = ctypes.CDLL(umat_file_name)
            # if platform.system() == "Windows":
            #     # WinDLL uses stdcall convention by default, which might be needed depending on compiler flags
            #     # CDLL uses cdecl. Often Fortran DLLs might align better with cdecl. Try both if one fails.
            #     umat_lib = ctypes.CDLL(umat_file_name)
            #     # umat_lib = ctypes.WinDLL(dll_path) # Alternative for Windows
            # else:

        except OSError as e:
            print(f"Error loading DLL: {e}")
            print("Ensure the DLL architecture (32/64-bit) matches your Python interpreter.")
            exit()

        # --- Get the function pointer ---
        fortran_function_found = False
        for fortran_func_name in possible_function_names:
            try:
                umat_func = getattr(umat_lib, fortran_func_name)
                fortran_function_found = True
                break
            except AttributeError:
                # Function not found, try the next one
                pass


        if not fortran_function_found:
            raise AttributeError (f"None of the function names  is found in the dll.")

        # --- Define Argument Data Types using ctypes ---
        # This MUST precisely match the Fortran subroutine's interface AND how Fortran passes them
        # Fortran usually passes EVERYTHING by reference (pointer).
        # We assume standard ABAQUS types: DOUBLE PRECISION for reals, default INTEGER for integers.

        c_double_p = ctypes.POINTER(ctypes.c_double)
        c_int_p = ctypes.POINTER(ctypes.c_int)  # Adjust if using non-default integers
        c_char_p = ctypes.c_char_p  # For CMNAME (tricky, see notes)

        # Define the expected argument types for the Fortran function
        # Order matters! Must match the UMAT subroutine definition exactly.
        umat_func.argtypes = [
            # STRESS(NTENS)        - Output/Input Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # STATEV(NSTATV)       - Output/Input Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # DDSDDE(NTENS,NTENS)  - Output Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
            # SSE, SPD, SCD        - Output Scalars (Double Precision) - Passed by reference
            c_double_p, c_double_p, c_double_p,
            # RPL                  - Output Scalar (Double Precision) - Passed by reference
            c_double_p,
            # DDSDDT(NTENS)        - Output Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # DRPLDE(NTENS)        - Output Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # DRPLDT               - Output Scalar (Double Precision) - Passed by reference
            c_double_p,
            # STRAN(NTENS)         - Input Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # DSTRAN(NTENS)        - Input Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # TIME(2)              - Input Array (Double Precision)
            c_double_p, # np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # DTIME                - Input Scalar (Double Precision) - Passed by reference
            c_double_p,
            # TEMP                 - Input Scalar (Double Precision) - Passed by reference
            c_double_p,
            # DTEMP                - Input Scalar (Double Precision) - Passed by reference
            c_double_p,
            # PREDEF(1)            - Input Array (Double Precision) - Size may vary
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # DPRED(1)             - Input Array (Double Precision) - Size may vary
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # CMNAME               - Input Character String - VERY TRICKY - see notes
            c_char_p,  # This might require a hidden length argument depending on compiler
            # NDI, NSHR, NTENS     - Input Scalars (Integer) - Passed by reference
            c_int_p, c_int_p, c_int_p,
            # NSTATV               - Input Scalar (Integer) - Passed by reference
            c_int_p,
            # PROPS(NPROPS)        - Input Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # NPROPS               - Input Scalar (Integer) - Passed by reference
            c_int_p,
            # COORDS(3)            - Input Array (Double Precision)
            c_double_p, #np.ctypeslib.ndpointer(dtype=np.float64, ndim=1, flags='F_CONTIGUOUS'),
            # DROT(3,3)            - Input Array (Double Precision)
            c_double_p, #np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
            # PNEWDT               - Output/Input Scalar (Double Precision) - Passed by reference
            c_double_p,
            # CELENT               - Input Scalar (Double Precision) - Passed by reference
            c_double_p,
            # DFGRD0(3,3)          - Input Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
            # DFGRD1(3,3)          - Input Array (Double Precision)
            np.ctypeslib.ndpointer(dtype=np.float64, ndim=2, flags='F_CONTIGUOUS'),
            # NOEL, NPT, LAYER, KSPT - Input Scalars (Integer) - Passed by reference
            c_int_p, c_int_p, c_int_p, c_int_p,
            # KSTEP, KINC          - Input Scalars (Integer) - Passed by reference
            c_int_p, c_int_p
            # --- Possible Hidden Arguments (Compiler Dependent) ---
            # String length for CMNAME might be added here implicitly by the compiler
            # ctypes.c_size_t # Example if CMNAME length is passed
        ]

        # Fortran subroutines don't typically return a value in the C sense
        umat_func.restype = None

        # --- Prepare Input Arguments ---
        # Allocate NumPy arrays - IMPORTANT: use order='F' for Fortran compatibility
        stress = np.array(stress, dtype=np.float64, order='F')
        statev = np.array(state_vars, dtype=np.float64, order='F')
        ddsdde = np.zeros((NTENS, NTENS), dtype=np.float64, order='F')
        ddsddt = np.zeros(NTENS, dtype=np.float64, order='F')
        drplde = np.zeros(NTENS, dtype=np.float64, order='F')
        stran = np.array(strain, dtype=np.float64, order='F')
        dstran = np.array(dstrain, dtype=np.float64, order='F')
        time_val = None #np.array([0.0, 1.0], dtype=np.float64, order='F')  # Example time step
        predef = np.zeros(1, dtype=np.float64, order='F')
        dpred = np.zeros(1, dtype=np.float64, order='F')
        props = np.array(props, dtype=np.float64, order='F')  # Material properties

        coords = None # np.array([0.0, 0.0, 0.0], dtype=np.float64, order='F')
        drot = None # np.eye(NDI, dtype=np.float64, order='F')
        dfgrd0 = np.eye(NDI, dtype=np.float64, order='F')
        dfgrd1 = np.eye(NDI, dtype=np.float64, order='F')

        # Create ctypes variables for scalars passed by reference
        sse = ctypes.c_double(0.0)
        spd = ctypes.c_double(0.0)
        scd = ctypes.c_double(0.0)
        rpl = ctypes.c_double(0.0)
        drpldt = ctypes.c_double(0.0)
        dtime = ctypes.c_double(0.0)
        temp = ctypes.c_double(0.0)
        dtemp = ctypes.c_double(0.0)
        pnewdt = ctypes.c_double(1.0)
        celent = ctypes.c_double(1.0)

        # Integers passed by reference
        c_ndi = ctypes.c_int(NDI)
        c_nshr = ctypes.c_int(NSHR)
        c_ntens = ctypes.c_int(NTENS)
        c_nstatv = ctypes.c_int(NSTATV)
        c_nprops = ctypes.c_int(NPROPS)
        c_noel = ctypes.c_int(1)
        c_npt = ctypes.c_int(1)
        c_layer = ctypes.c_int(0)
        c_kspt = ctypes.c_int(0)
        c_kstep = ctypes.c_int(time_step)
        c_kinc = ctypes.c_int(1)

        # --- CMNAME Handling (VERY DELICATE) ---
        # Fortran passes fixed-length strings often with a hidden length argument.
        cmname_len = 80  # Standard Abaqus length
        cmname_buf = ctypes.create_string_buffer(cmname_len)
        cmname_buf.value = b"MYMATERIAL".ljust(cmname_len)  # Pad with spaces if needed

        # --- Call the Fortran function ---
        umat_func(
            stress, statev, ddsdde,
            ctypes.byref(sse), ctypes.byref(spd), ctypes.byref(scd),  # Pass scalars by reference
            ctypes.byref(rpl),
            ddsddt, drplde,
            ctypes.byref(drpldt),
            stran, dstran, time_val,
            ctypes.byref(dtime), ctypes.byref(temp), ctypes.byref(dtemp),
            predef, dpred,
            cmname_buf,  # Pass the buffer or simple bytes
            ctypes.byref(c_ndi), ctypes.byref(c_nshr), ctypes.byref(c_ntens),
            ctypes.byref(c_nstatv),
            props,
            ctypes.byref(c_nprops),
            coords, drot,
            ctypes.byref(pnewdt),
            ctypes.byref(celent),
            dfgrd0, dfgrd1,
            ctypes.byref(c_noel), ctypes.byref(c_npt), ctypes.byref(c_layer), ctypes.byref(c_kspt),
            ctypes.byref(c_kstep), ctypes.byref(c_kinc)
            # Add hidden string length argument here if needed , ctypes.c_size_t(cmname_len)
        )

        # transfer fortran output to numpy arrays
        stress = np.ctypeslib.as_array(stress)
        statev = np.ctypeslib.as_array(statev)
        ddsdde = np.ctypeslib.as_array(ddsdde)


        return stress, ddsdde, statev

    @staticmethod
    def run_c_umat(umat_file_name, stress, state_vars, strain, dstrain, props, time_step, ndi=3, nshr=3):
        """
        Call the UMAT function from Python using CFFI.

        Temperature and predefined field variables are not used in this function.

        Args:
            umat_file_name (str): Path to the compiled UMAT shared library (DLL)
            stress (Sequence[float]): Initial stress vector [σ11, σ22, σ33, σ12, σ13, σ23]
            state_vars (Sequence[float]): State variables array
            strain (Sequence[float]): Total strain vector
            dstrain (Sequence[float]): Increment in strain vector
            props (Sequence[float]): Material properties array
            ndi (int): Number of direct stress components (default is 3 for 3D)
            nshr (int): Number of shear stress components (default is 3 for 3D)

        Returns:
            stress (np.ndarray): Updated stress vector
            ddsdde (np.ndarray): Updated Jacobian matrix
            state_vars (np.ndarray): Updated state variables
        """

        # Convert inputs to lists if they are not already
        stress = list(stress)
        state_vars = list(state_vars)
        strain = list(strain)
        dstrain = list(dstrain)
        props = list(props)

        nstatv = len(state_vars)

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
        c_dtime = ffi.new("double*", time_step)

        # Material name (80 characters)
        c_cmname = ffi.new("char[]", b"MY_MATERIAL" + b" " * (80 - len(b"MY_MATERIAL")))

        # Integer parameters
        c_ndi = ffi.new("int*", ndi)
        c_nshr = ffi.new("int*", nshr)
        c_ntens = ffi.new("int*", ntens)
        c_nstatv = ffi.new("int*", nstatv)

        c_props = ffi.new("double[]", props)
        c_nprops = ffi.new("int*", len(props))

        # # Rotation matrix (identity)
        # c_drot = ffi.new("double[]", [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0])

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
            c_sse, c_spd, c_scd, ffi.NULL,
            ffi.NULL, ffi.NULL, ffi.NULL,  # DDSDDT, DRPLDE, DRPLDT (null)
            c_strain, c_dstrain, c_time, c_dtime,
            ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL,  # TEMP, DTEMP, PREDEF, DPRED (null)
            c_cmname, c_ndi, c_nshr, c_ntens, c_nstatv,
            c_props, c_nprops,
            ffi.NULL,  # COORDS (null)
            ffi.NULL, # DROT (null)
            ffi.NULL, ffi.NULL, ffi.NULL, ffi.NULL,  # PNEWDT, CELENT, DFGRD0, DFGRD1 (null)
            c_noel, c_npt, ffi.NULL, ffi.NULL, c_kstep, c_kinc
        )

        # Convert C arrays back to numpy arrays
        stress_updated = np.array([c_stress[i] for i in range(ntens)])
        statev_updated = np.array([c_statev[i] for i in range(nstatv)])

        # Reshape the flat DDSDDE array to a 2D matrix
        ddsdde_updated = np.array([[c_ddsdde[i * ntens + j] for j in range(ntens)]
                                   for i in range(ntens)])

        return stress_updated, ddsdde_updated, statev_updated
