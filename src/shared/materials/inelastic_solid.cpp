#include "inelastic_solid.h"
#include "base_particles.hpp"

namespace SPH
{
//=================================================================================================//
void HardeningPlasticSolid::initializeLocalParameters(BaseParticles *base_particles)
{
    PlasticSolid::initializeLocalParameters(base_particles);
    base_particles->registerVariable(inverse_plastic_strain_, "InversePlasticRightCauchyStrain",
                                     [&](size_t i) -> Matd
                                     { return Matd::Identity(); });
    base_particles->registerVariable(hardening_parameter_, "HardeningParameter");
    base_particles->addVariableToRestart<Matd>("InversePlasticRightCauchyStrain");
    base_particles->addVariableToRestart<Real>("HardeningParameter");
}
//=================================================================================================//
Matd HardeningPlasticSolid::PlasticConstitutiveRelation(const Matd &F, size_t index_i, Real dt)
{
    Matd be = F * inverse_plastic_strain_[index_i] * F.transpose();
    Matd normalized_be = be * pow(be.determinant(), -OneOverDimensions);
    Real normalized_be_isentropic = normalized_be.trace() * OneOverDimensions;
    Matd deviatoric_PK = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd::Identity());
    Real deviatoric_PK_norm = deviatoric_PK.norm();
    Real trial_function = deviatoric_PK_norm -
                          sqrt_2_over_3_ * (hardening_modulus_ * hardening_parameter_[index_i] + yield_stress_);
    if (trial_function > 0.0)
    {
        Real renormalized_shear_modulus = normalized_be_isentropic * G0_;
        Real relax_increment = 0.5 * trial_function / (renormalized_shear_modulus + hardening_modulus_ / 3.0);
        hardening_parameter_[index_i] += sqrt_2_over_3_ * relax_increment;
        deviatoric_PK -= 2.0 * renormalized_shear_modulus * relax_increment * deviatoric_PK / deviatoric_PK_norm;
        Matd relaxed_be = deviatoric_PK / G0_ + normalized_be_isentropic * Matd::Identity();
        normalized_be = relaxed_be * pow(relaxed_be.determinant(), -OneOverDimensions);
    }
    Matd inverse_F = F.inverse();
    Matd inverse_F_T = inverse_F.transpose();
    inverse_plastic_strain_[index_i] = inverse_F * normalized_be * inverse_F_T;

    return (deviatoric_PK + VolumetricKirchhoff(F.determinant()) * Matd::Identity()) * inverse_F_T;
}
//=================================================================================================//
void ViscousPlasticSolid::initializeLocalParameters(BaseParticles *base_particles)
{
    PlasticSolid::initializeLocalParameters(base_particles);
    base_particles->registerVariable(inverse_plastic_strain_, "InversePlasticRightCauchyStrain",
                                     [&](size_t i) -> Matd
                                     { return Matd::Identity(); });
    base_particles->addVariableToRestart<Matd>("InversePlasticRightCauchyStrain");
}
//=================================================================================================//
Matd ViscousPlasticSolid::PlasticConstitutiveRelation(const Matd &F, size_t index_i, Real dt)
{
    Matd be = F * inverse_plastic_strain_[index_i] * F.transpose();
    Matd normalized_be = be * pow(be.determinant(), -OneOverDimensions);
    Real normalized_be_isentropic = normalized_be.trace() * OneOverDimensions;
    Matd deviatoric_K = DeviatoricKirchhoff(normalized_be - normalized_be_isentropic * Matd::Identity());
    Real deviatoric_K_norm = deviatoric_K.norm();
    Real trial_function = deviatoric_K_norm - sqrt_2_over_3_ * yield_stress_; 
	if (trial_function > 0.0)
    {
        Real renormalized_shear_modulus = normalized_be_isentropic * G0_;
        Real s_Mid = 0.0;
        Real s_Max = deviatoric_K_norm;
        Real s_Min = sqrt_2_over_3_ * yield_stress_;
        Real gfunc = 0.0;
        Real Precision = 1.0e-6;
        Real Relative_Error;
        do
        {
            s_Mid = (s_Max + s_Min) / 2.0;
            gfunc = pow(viscous_modulus_, 1.0 / Herschel_Bulkley_power_) * (s_Mid - deviatoric_K_norm) + 
            2.0 * renormalized_shear_modulus * dt * pow((s_Mid - sqrt_2_over_3_ * yield_stress_), 1.0 / Herschel_Bulkley_power_);
            if (gfunc < 0.0)
            {
                s_Min = s_Mid;
            }
            else
            {
                s_Max = s_Mid;
            }
            Relative_Error = gfunc / deviatoric_K_norm;
        } while (fabs(Relative_Error) >= Precision); 
        deviatoric_K = s_Mid * deviatoric_K / deviatoric_K_norm;
        Matd relaxed_be = deviatoric_K / G0_ + normalized_be_isentropic * Matd::Identity();
        normalized_be = relaxed_be * pow(relaxed_be.determinant(), -OneOverDimensions);
    }
    Matd inverse_F = F.inverse();
    Matd inverse_F_T = inverse_F.transpose();
    be = pow(F.determinant(), (2.0 / 3.0)) * normalized_be;
    inverse_plastic_strain_[index_i] = inverse_F * be * inverse_F_T;

     return (deviatoric_K + VolumetricKirchhoff(F.determinant()) * Matd::Identity()) * inverse_F_T;
}
} // namespace SPH
