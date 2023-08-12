#include "inelastic_dynamics.h"

namespace SPH
{
//=====================================================================================================//
namespace solid_dynamics
{
//=================================================================================================//
PlasticIntegration1stHalf::
    PlasticIntegration1stHalf(BaseInnerRelation &inner_relation) 
    : Integration1stHalf(inner_relation),
    plastic_solid_(DynamicCast<PlasticSolid>(this, elastic_solid_))
{
    numerical_dissipation_factor_ = 0.5;
}
//=================================================================================================//
void PlasticIntegration1stHalf::initialization(size_t index_i, Real dt)
{
    pos_[index_i] += vel_[index_i] * dt * 0.5;
    F_[index_i] += dF_dt_[index_i] * dt * 0.5;
    // �� = ��0 / J ; J = det[F]
    rho_[index_i] = rho0_ / F_[index_i].determinant();
    // obtain Kirchhoff stress tensor
    // Matd inverse_F_T = F_[index_i].inverse().transpose();
    // stress_PK1_B_[index_i] = plastic_solid_.PlasticConstitutiveRelation(F_[index_i], index_i, dt) * inverse_F_T * B_[index_i];
    
    // obtain the first Piola-Kirchhoff stress from Kirchhoff stress tensor
     stress_PK1_B_[index_i] = plastic_solid_.PlasticConstitutiveRelation(F_[index_i], index_i, dt) * B_[index_i]; 
}
//=================================================================================================//
} // namespace solid_dynamics
  //=====================================================================================================//
} // namespace SPH
