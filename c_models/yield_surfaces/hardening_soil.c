
double shear_yield_function(const double q_a,const double q, const double E_50, const double E_ur, const double gamma_ps)
{


//double f_1 = q_a/E_50 * q /(q_a - q) - 2* q/ E_ur - gamma_ps;

// plaxis manual:

double f_1 = 2 * q_a/E_50 * q /(q_a - q) - 2* q/ E_ur - gamma_ps;

return f_1;

}


void calculate_dqa_dsigma(double dqadsigma[VOIGTSIZE_3D], const double phi_rad, const double R_f)
{
    for (int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        dqadsigma[i] = 0.0;
    }

    dqadsigma[2] = -2.0*sin(phi_rad) / (1.0 - sin(phi_rad)) / R_f;

}

void derivative_shear_yield_function(const double q_a, const double q, const double E_50, const double E_ur, const double gamma_ps, double dfds[VOIGTSIZE_3D])
{
    // derivative of shear yield function

    double dqdsigma[VOIGTSIZE_3D];
    calculate_dq_dsigma_triaxial_state_3d(dqdsigma)


    double dqadsigma[VOIGTSIZE_3D];
    calculate_dqa_dsigma(dqadsigma, phi_rad_peak, R_f);

    double dfdq = (2.0 * q_a*q_a) / (E_50 * pow((q_a - q),2.0)) - 2.0 / E_ur

    double dfdqa = - (2 * q*q_a) / (E_50 * pow((q_a - q),2))

    double dfds[VOIGTSIZE_3D];

    for( int i = 0; i < VOIGTSIZE_3D; ++i)
    {
        dfds[i] = dfdq * dqdsigma[i] + dfdqa * dqadsigma[i];
    }

}

//double cap_yield_function(const double q, const double M, const double r, const double p, const double pc)
//{
//    // Cap yield function
//    double f_2 = q * q / (M*M * r*r) + p*p - pc*pc;
//
//    return f_2;
//}

double cap_yield_function(const double q_special, const double M, const double p, const double pc)
{
    // Cap yield function from plaxis manual
    double f_2 = q_special * q_special / (M*M ) + p*p - pc*pc;

    return f_2;
}

double derivative_cap_yield_function(const double q_special, const double M, const double p,)
{
    // Cap yield function
    dqdsigma = calcukate_dqds();
    dpdsigma = calculate_dpdsigma();


    double dfds = 2.0 * q_special / (M*M) * dqdsigma + 2.0*p * dpdsigma;

    return dfds;
}

double calculate_special_deviatoric_stress(const double sigma_1, const double sigma_2, const double sigma_3, const double phi_rad)
{
    const double alpha = 3 + sin(phi_rad) / (1 - sin(phi_rad));

    // triaxial compression
    if (abs(-sigma_2 - sigma_3) < EPSILON)
    {
        return -(sigma_1-sigma_3);
    }
    if (-sigma_1 > -sigma)

    return sigma_1 + (alpha -1) * sigma_2 - alpha* sigma_3;

}


// plastic shear strain
double calculate_shear_hardening_parameter(const double plastic_strain_1, const double plastic_strain_vol)
{
double gamma_ps = -(2* plastic_strain_1 - plastic_strain_vol) ;
return gamma_ps;
}

double calculate_ultimate_deviatoric_stress(const double c, const double phi_rad,  const double sigma_3)
{
    // follows from mohr coulomb yield criterion
    double q_f = 2 * sin(phi_rad) / (1- sin(phi_rad)) * (sigma_3 +  c /tan(phi_rad));

    return q_f;
    // Calculate the deviatoric ultimate stress based on the given parameters
    *q_a = 2.0 * c * cos(phi_rad) * (1.0 - (sigma_3 / p_ref) * sin(phi_rad)) / (1.0 - m);
}

double calculate_plastic_volumetric_strain_rate(const double dilatancy_angle_mob, const double gamma_ps_rate)
{
    return sin(dilatancy_angle_mob) * gamma_ps_rate;
}

double calculate_mobilised_dilatancy_angle(const double sin_phi_mob, const double friction_angle, const double sin_phi_cv, const double dilatancy_angle)
{

    const double phi = friction_angle;
    const double psi = dilatancy_angle;




    if (sin_phi_mob < 3.0/4.0 * sin(phi))
    {
        return 0.0;
    }
    else if (sin_phi_mob >= sin(phi))
    {
        if (psi >0)
        {
            double sin_psi_mob = max(0.0, (sin_phi_mob-sin_phi_cv)/(1-sin_phi_mob*sin_phi_cv));
            return asin(sin_psi_mob);
        }
        else{
            return psi;
        }
    }
    else if (phi < EPSILON)
    {
        return 0.0;
    }

}

double calculate_mobilised_sin_friction_angle(const double sigma_1, const double sigma_3, const double c, const double phi_rad)
{
    // Calculate the mobilised sin of the friction angle based on the given parameters
    return (sigma_1 - sigma_3) / (sigma_1 + sigma_3 -2.0 * c /tan(phi_rad));

}

double calculate_critical_sin_friction_angle(const double phi_rad, const double psi_rad)
{
    return (sin(phi_rad)- sin(psi_rad))/(1- sin(phi_rad)*sin(psi_rad));
}
