
double shear_yield_function(const double q_a,const double q, const double E_50, const double E_ur, const double gamma_ps)
{


//double f_1 = q_a/E_50 * q /(q_a - q) - 2* q/ E_ur - gamma_ps;

// plaxis manual:
double f_1 = q_a/E_50 * q /(q_a - q) - 2* q/ E_ur - gamma_ps;

return f_1;

}

double cap_yield_function(const double q, const double M, const double r, const double p, const double pc)
{
    // Cap yield function
    double f_2 = q * q / (M*M * r*r) + p*p - pc*pc;

    return f_2;
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

double calculate_mobilized_dilatancy_angle(const double friction_angle_mob, const double friction_angle, const double friction_angle_cv, const double dilatancy_angle)
{
    if (sin(friction_angle_mob) < 3.0/4.0 * sin(friction_angle))
    {
        return 0.0;
    }
    else if (sin(friction_angle_mob) >= sin(friction_angle))
    {
    if dilatancy_angle
        return friction_angle_cv;
    }
    else
    {
        return friction_angle_mob;
    }
}