#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.02; 
Vec2d ball_center_1(-0.52, 0.0);
Vec2d ball_center_2(0.52, 0.0);
BoundingBox system_domain_bounds(Vec2d(-2.0, -2.0),
                                 Vec2d(2.0, 2.0));
Real ball_radius_1 = 0.5;
Real ball_radius_2 = 0.5;

StdVec<Vecd> ball_observation_location_1 = {ball_center_1};
StdVec<Vecd> ball_observation_location_2 = {ball_center_2};
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real physical_viscosity = 1000.0;
Real rho0_s = 1.0e3;
Real Bulk_modulus = 1.09e5;
Real Shear_modulus = 1.12e4;
Real yield_stress = 0.1;
Real viscous_modulus = 10.0;
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus);
Real Herschel_Bulkley_power = 1.0;

class Ball1InitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit Ball1InitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][0] = 0.05;
    }
};
class Ball2InitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit Ball2InitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

    void update(size_t index_i, Real dt)
    {
        vel_[index_i][0] = -0.05;
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false); // true
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true); // false
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody ball_1(sph_system, makeShared<GeometricShapeBall>(ball_center_1, ball_radius_1, "Ball_1"));
    ball_1.defineBodyLevelSetShape();
    ball_1.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                  yield_stress, viscous_modulus, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? ball_1.generateParticles<ParticleGeneratorReload>(io_environment, ball_1.getName())
        : ball_1.generateParticles<ParticleGeneratorLattice>();

    SolidBody ball_2(sph_system, makeShared<GeometricShapeBall>(ball_center_2, ball_radius_2, "Ball_2"));
    ball_2.defineBodyLevelSetShape();
    ball_2.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                  yield_stress, viscous_modulus, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? ball_2.generateParticles<ParticleGeneratorReload>(io_environment, ball_2.getName())
        : ball_2.generateParticles<ParticleGeneratorLattice>();

    ObserverBody ball_observer_1(sph_system, "BallObserver1");
    ball_observer_1.generateParticles<ObserverParticleGenerator>(ball_observation_location_1);
    ObserverBody ball_observer_2(sph_system, "BallObserver2");
    ball_observer_2.generateParticles<ObserverParticleGenerator>(ball_observation_location_2);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation ball_inner_1(ball_1);
        InnerRelation ball_inner_2(ball_2);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for ball.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> ball_random_particles_1(ball_1);
        relax_dynamics::RelaxationStepInner ball_relaxation_step_inner_1(ball_inner_1);
        SimpleDynamics<RandomizeParticlePosition> ball_random_particles_2(ball_2);
        relax_dynamics::RelaxationStepInner ball_relaxation_step_inner_2(ball_inner_2);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&ball_1, &ball_2});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        ball_random_particles_1.exec(0.25);
        ball_random_particles_2.exec(0.25);
        write_relaxed_particles.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            ball_relaxation_step_inner_1.exec();
            ball_relaxation_step_inner_2.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_relaxed_particles.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation ball_inner_1(ball_1);
    InnerRelation ball_inner_2(ball_2);
    SelfSurfaceContactRelation ball_self_contact_1(ball_1);
    SelfSurfaceContactRelation ball_self_contact_2(ball_2);
    SurfaceContactRelation ball_contact_1(ball_self_contact_1, {&ball_2});
    SurfaceContactRelation ball_contact_2(ball_self_contact_2, {&ball_1});
    ContactRelation ball_observer_contact_1(ball_observer_1, {&ball_1});
    ContactRelation ball_observer_contact_2(ball_observer_2, {&ball_2});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<Ball1InitialCondition> initial_condition_1(ball_1);
    SimpleDynamics<Ball2InitialCondition> initial_condition_2(ball_2);
    InteractionWithUpdate<CorrectedConfigurationInner> ball_corrected_configuration_1(ball_inner_1);
    InteractionWithUpdate<CorrectedConfigurationInner> ball_corrected_configuration_2(ball_inner_2);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_get_time_step_size_1(ball_1, 0.1);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_get_time_step_size_2(ball_2, 0.1);
    /** stress relaxation for the balls. */
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> ball_stress_relaxation_first_half_1(ball_inner_1);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_stress_relaxation_second_half_1(ball_inner_1);
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> ball_stress_relaxation_first_half_2(ball_inner_2);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_stress_relaxation_second_half_2(ball_inner_2);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> ball_update_contact_density_1(ball_contact_1);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_compute_solid_contact_forces_1(ball_contact_1);
    InteractionDynamics<solid_dynamics::ContactDensitySummation> ball_update_contact_density_2(ball_contact_2);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_compute_solid_contact_forces_2(ball_contact_2);
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center_1(ball_1, Vecd(1.0, 1.0));
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center_2(ball_2, Vecd(1.0, 1.0));
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        ball_1_position_damping(0.5, ball_inner_1, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        ball_2_position_damping(0.5, ball_inner_2, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_ball_center_displacement_1("Position", io_environment, ball_observer_contact_1);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_ball_center_displacement_2("Position", io_environment, ball_observer_contact_2);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    ball_corrected_configuration_1.exec();
    ball_corrected_configuration_2.exec();
    initial_condition_1.exec();
    initial_condition_2.exec();
    /** Initial states output. */
    body_states_recording.writeToFile(0);
    write_ball_center_displacement_1.writeToFile(0);
    write_ball_center_displacement_2.writeToFile(0);
    /** Main loop. */
    int ite = 0;
    Real T0 = 3.0;
    Real end_time = T0;
    Real output_interval = 0.02 * T0;
    Real Dt = 0.1 * output_interval;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                initial_condition_1.exec();
                initial_condition_2.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                ball_update_contact_density_1.exec();
                ball_compute_solid_contact_forces_1.exec();
                ball_stress_relaxation_first_half_1.exec(dt);
                ball_1_position_damping.exec(dt);
                constrain_mass_center_1.exec(dt);
                ball_stress_relaxation_second_half_1.exec(dt);

                ball_update_contact_density_2.exec();
                ball_compute_solid_contact_forces_2.exec();
                ball_stress_relaxation_first_half_2.exec(dt);
                ball_2_position_damping.exec(dt);
                constrain_mass_center_2.exec(dt);
                ball_stress_relaxation_second_half_2.exec(dt);

                ball_1.updateCellLinkedList();
                ball_2.updateCellLinkedList();
                ball_contact_1.updateConfiguration();
                ball_contact_2.updateConfiguration();

                ite++;
                Real dt_1 = ball_get_time_step_size_1.exec();
                Real dt_2 = ball_get_time_step_size_2.exec();
                dt = SMIN(dt_1, dt_2);
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                write_ball_center_displacement_1.writeToFile(ite);
            }
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
