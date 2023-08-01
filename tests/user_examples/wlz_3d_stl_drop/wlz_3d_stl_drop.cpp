/**
 * this file is following the 3d_self_contact.cpp
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/Xiaohuangren.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-10.0, -10.0, -10.0);
Vec3d domain_upper_bound(10.0, 10.0, 5.0);
Real resolution_ref = (domain_upper_bound[0] - domain_lower_bound[0]) / 200.0;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
StdVec<Vecd> imported_model_observation_location = {Vec3d(0.0, 0.0, 0.0)};
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
// viscoplastic material
Real physical_viscosity = 1000.0;
Real gravity_g = 2.0;
Real rho0_s = 1.0e3;                                                                                    /* ρ density. kg/m^3 */
Real Bulk_modulus = 1.09e5;                                                                             /* κ bulk modulus. Pa */
Real Shear_modulus = 1.12e4;                                                                            /* μ/G shear modulus. Pa */
Real yield_stress = 0.1;                                                                                /* σ_Y yield stress. Pa */
Real viscous_modulus = 10.0;                                                                            /* η viscosity.  */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa. 此处32487 */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* ν Poisson's ratio. 取值(-1,0.5). 此处0.45 */
Real Herschel_Bulkley_power = 1.0;
//----------------------------------------------------------------------
//	Body shapes used in the case.
//----------------------------------------------------------------------
class Xiaohuangren : public ComplexShape
{
  public:
    explicit Xiaohuangren(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_file, Vecd::Zero(), 1.0);
    }
};
class StationaryPlate : public ComplexShape
{
  public:
    explicit StationaryPlate(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_plate(9.5, 9.5, 0.25);
        Vecd translation_plate(0.0, 0.0, -7.0);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_plate), halfsize_plate);
    }
};
//----------------------------------------------------------------------
//	The main program
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, resolution_ref);
    // Tag for run particle relaxation for the initial body fitted distribution.
    system.setRunParticleRelaxation(false); // false
    // Tag for reload initially relaxed particles.
    system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody xiaohuangren(system, makeShared<Xiaohuangren>("Xiaohuangren"));
    xiaohuangren.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    xiaohuangren.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                        yield_stress, viscous_modulus, Herschel_Bulkley_power);
    (!system.RunParticleRelaxation() && system.ReloadParticles())
        ? xiaohuangren.generateParticles<ParticleGeneratorReload>(io_environment, xiaohuangren.getName())
        : xiaohuangren.generateParticles<ParticleGeneratorLattice>();

    SolidBody stationary_plate(system, makeShared<StationaryPlate>("StationaryPlate"));
    stationary_plate.defineParticlesAndMaterial<SolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    stationary_plate.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation xiaohuangren_inner(xiaohuangren);
    SelfSurfaceContactRelation xiaohuangren_self_contact(xiaohuangren);
    SurfaceContactRelation xiaohuangren_contact(xiaohuangren_self_contact, {&stationary_plate});
    //----------------------------------------------------------------------
    //	check whether run particle relaxation for body fitted particle distribution.
    //----------------------------------------------------------------------
    if (system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        // Random reset the insert body particle position.
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(xiaohuangren);
        // Write the particle reload files.
        ReloadParticleIO write_particle_reload_files(io_environment, xiaohuangren);
        // A  Physics relaxation step.
        relax_dynamics::RelaxationStepInner relaxation_step_inner(xiaohuangren_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_states.writeToFile(0);
        //----------------------------------------------------------------------
        //	Particle relaxation loop.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        // Output particles position for reload.
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	This section define all numerical methods will be used in this case.
    //----------------------------------------------------------------------
    // initialize a time step
    SimpleDynamics<TimeStepInitialization> initialization_with_gravity(xiaohuangren, makeShared<Gravity>(Vecd(0.0, 0.0, -gravity_g)));
    // Corrected configuration for reproducing rigid rotation.
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration(xiaohuangren_inner);
    // Time step size
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(xiaohuangren, 0.1);
    // stress relaxation.
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> stress_relaxation_first_half(xiaohuangren_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(xiaohuangren_inner);
    // Algorithms for solid-solid contacts.
    InteractionDynamics<solid_dynamics::ContactDensitySummation> xiaohuangren_update_contact_density(xiaohuangren_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> xiaohuangren_compute_solid_contact_forces(xiaohuangren_contact);
    InteractionDynamics<solid_dynamics::SelfContactDensitySummation> xiaohuangren_self_contact_density(xiaohuangren_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> xiaohuangren_self_contact_forces(xiaohuangren_self_contact);
    // Damping the velocity field for quasi-static solution
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        coil_damping(0.2, xiaohuangren_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	From here the time stepping begins.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    system.initializeSystemConfigurations();
    // apply initial condition
    corrected_configuration.exec();
    write_states.writeToFile(0);
    // Setup time stepping control parameters.
    int ite = 0;
    Real end_time = 10.0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    // Statistics for computing time.
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: "
                          << dt << "\n";
            }
            initialization_with_gravity.exec();
            // contact dynamics.
            xiaohuangren_self_contact_density.exec();
            xiaohuangren_self_contact_forces.exec();
            xiaohuangren_update_contact_density.exec();
            xiaohuangren_compute_solid_contact_forces.exec();
            // Stress relaxation and damping.
            stress_relaxation_first_half.exec(dt);
            coil_damping.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;

            // update particle neighbor relations for contact dynamics
            xiaohuangren.updateCellLinkedList();
            xiaohuangren_self_contact.updateConfiguration();
            xiaohuangren_contact.updateConfiguration();
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
