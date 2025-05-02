

int calculate_cap_hardening()
{
    double Ddgdsigma[VOIGTSIZE_3D];

    f = cap_yield_function(q_special, M, p, pc);
    if (f< ZERO_TOL)
    {
    // elastic, converged
        return 1;
    }
    else
    {

        res =  calculate_residual();

        while (iter < max_iter && abs(res) > ZERO_TOL)
        {
            // calculate K cap
            double K_cap = calculate_K_cap();

            if (iter == 0)
            {
            dgdsigma = calculate_dgdsigma();
            }

            dfdh = calculate_dfdh();
            dhdlamda = calculate_dhdlambda();

            dfdsigma = calculate_dfdsigma();


            matrix_vector_multiply(D, dgdsigma, VOIGTSIZE_3D, Ddgdsigma);

            double dfdsigmaDdgdsigma =  vector_dot_product(dfdsigma, Ddgdsigma, VOIGTSIZE_3D);

            dfdlambda = dfdsigmaDdgdsigma - dfdh * dhdlamda;

            // implicit
            d_lambda = d_lambda + res /(dfdlambda); // + viscosity / (d_t + SMALL_VALUE));

            d_eps_p = d_lambda * dgdsigma;
            d_eps_p_vol = calculate_volumetric_strain_3d(d_eps_p);
            eps_p_vol = eps_p_vol_state + d_eps_p_vol;

            princ_stress = princ_stress_elast - d_lambda * Ddgdsigma;
            // resort princ stress
            sort(princ_stress, VOIGTSIZE_3D);

            pc = calculate_pc(eps_p_vol, pc_state);

            pc = max(pc, pc_state);

            f = cap_yield_function(q_special, M, p, pc);

            res = calculate_residual(f);

        }

        if (res < ZERO_TOL)
        {
            // is converged
            state_void = e_current;
            state_pc = pc;
            state_loading_type = "cap";
            is_converged = 1;
        }

    }

}

double calculate_residual(double f)
{
    return f;
}
