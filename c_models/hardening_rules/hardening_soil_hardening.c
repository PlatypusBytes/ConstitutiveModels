

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
        int i = 0;
        while (i < max_iter && abs(res) > ZERO_TOL)
        {
            // calculate K cap
            double K_cap = calculate_K_cap();

            dgdsigma = calculate_dgdsigma();
            // normalize dgdsigma


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

            i++;

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

int calculate_cone_hardening()
{
    f = cone_yield_function();
    res = calculate_residual(f);

    int i = 0;
    while (i < max_iter && abs(res) > ZERO_TOL)
    {
        // calculate cone parameters
        parameters = calculate_cone_parameters();

        double dgdsigma[VOIGTSIZE_3D];

        derivative_shear_yield_function(const double q_a, const double q, const double E_50, const double E_ur, const double gamma_ps, dgdsigma)

        dgdsigma = calculate_cone_dgdsigma();
        // normalize


        dfdh = calculate_dfdh_cone();
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

        f = cone_yield_function();

        res = calculate_residual(f);

        i++;

    }

}

double calculate_residual(double f)
{
    return f;
}


double calculate_dfdh_cone()
{
    return -1.0;
}

double calculate_dhdlambda_cap()
{
 // dpc = -(1+e)/(lambda_star - kappa_star) * p_c * d_eps_p_vol;
 // d_eps_p_vol = (d_eps_1 + d_eps_2 + d_eps_3)
 // d_eps = d_lambda * dgdsigma
 // d_eps_p = d_lambda * sum(dgdsigma) , in principle stress space
 // dpc = -(1+e)/(lambda_star - kappa_star) * p_c * d_lambda * sum(dgdsigma);

    double dpcdlambda = -(1+e)/(lambda_star - kappa_star) * p_c * sum(dgdsigma);

    return dpcdlambda;
}

double calculate_dhdlambda_cone()
{
    double dpcdlambda = calculate_dhdlambda_cap();
    double dhdlamda = 1.0 + dpcdlambda;
    return dhdlamda;
}
