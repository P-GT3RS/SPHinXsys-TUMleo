#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.025; /* reference resolution. */
Vec3d ball_center(0.0, 0.0, 1.0);
BoundingBox system_domain_bounds(Vec3d(-1.5, -1.5, -0.5),
                                 Vec3d(1.5, 1.5, 2.0));
StdVec<Vecd> ball_observation_location = {ball_center};
Real ball_radius = 0.5; /* R radius. m */
Real column_radius = 0.15;
Real column_length = 2.0;
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
// viscoplastic
Real gravity_g = 2.0;                                                                                   /* 过大会导致球体直接穿过地面,过小会导致球体刚接触地面就停止计算 */
Real rho0_s = 1.0e3;                                                                                    /* ρ density. kg/m^3 */
Real Bulk_modulus = 1.09e5;                                                                             /* κ bulk modulus. Pa */
Real Shear_modulus = 1.12e4;                                                                            /* μ/G shear modulus. Pa */
Real yield_stress = 0.1;                                                                                /* σ_Y yield stress. Pa */
Real viscous_modulus = 10.0;                                                                            /* η viscosity. 过大会导致理想粘塑性产生明显回弹,过小会导致球体无法产生回弹 */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* ν Poisson's ratio. 取值(-1,0.5). 此处0.45 */
Real Herschel_Bulkley_power = 1.0;                                                                      /* h Herschel_Bulkley_power. */

class Column_1 : public ComplexShape
{
  public:
    explicit Column_1(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column_1(0.0, -0.25, -0.1);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0.0, 0.0), column_radius,
                                       column_length / 2.0, 50, translation_column_1);
    }
};
class Column_2 : public ComplexShape
{
  public:
    explicit Column_2(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_column_2(0.0, 0.25, -0.1);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0.0, 0.0), column_radius,
                                       column_length / 2.0, 50, translation_column_2);
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
    SolidBody column_1(sph_system, makeShared<Column_1>("Column_1"));
    column_1.defineBodyLevelSetShape();
    column_1.defineParticlesAndMaterial<SolidParticles, Solid>(rho0_s, Youngs_modulus, poisson);
    column_1.generateParticles<ParticleGeneratorLattice>();

    SolidBody column_2(sph_system, makeShared<Column_2>("Column_2"));
    column_2.defineBodyLevelSetShape();
    column_2.defineParticlesAndMaterial<SolidParticles, Solid>(rho0_s, Youngs_modulus, poisson);
    column_2.generateParticles<ParticleGeneratorLattice>();

    // 问题:用此方法生成八面体而非球体
    // SolidBody ball(sph_system, makeShared<Ball>("Ball"));
    // ball.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    // ball.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
    //	yield_stress, viscous_modulus, Herschel_Bulkley_power);
    //(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    //	? ball.generateParticles<ParticleGeneratorReload>(io_environment, ball.getName())
    //	: ball.generateParticles<ParticleGeneratorLattice>();

    SolidBody ball(sph_system, makeShared<GeometricShapeBall>(ball_center, ball_radius, "BallBody"));
    ball.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    ball.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                yield_stress, viscous_modulus, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? ball.generateParticles<ParticleGeneratorReload>(io_environment, ball.getName())
        : ball.generateParticles<ParticleGeneratorLattice>();

    ObserverBody ball_observer(sph_system, "BallObserver");
    ball_observer.generateParticles<ObserverParticleGenerator>(ball_observation_location);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation ball_inner(ball);
        InnerRelation column_1_inner(column_1);
        InnerRelation column_2_inner(column_2);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for ball.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> ball_random_particles(ball);
        SimpleDynamics<RandomizeParticlePosition> column_1_random_particles(column_1);
        SimpleDynamics<RandomizeParticlePosition> column_2_random_particles(column_2);
        relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
        relax_dynamics::RelaxationStepInner column_1_relaxation_step_inner(column_1_inner);
        relax_dynamics::RelaxationStepInner column_2_relaxation_step_inner(column_2_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&ball, &column_1, &column_2});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        ball_random_particles.exec(0.25);
        column_1_random_particles.exec(0.25);
        column_2_random_particles.exec(0.25);
        write_relaxed_particles.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            ball_relaxation_step_inner.exec();
            column_1_relaxation_step_inner.exec();
            column_2_relaxation_step_inner.exec();
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
    InnerRelation ball_inner(ball);
    SurfaceContactRelation ball_contact(ball, {&column_1, &column_2});
    ContactRelation ball_observer_contact(ball_observer, {&ball});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<TimeStepInitialization> ball_initialize_timestep(ball, makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g)));
    InteractionWithUpdate<CorrectedConfigurationInner> ball_corrected_configuration(ball_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_get_time_step_size(ball, 0.1);
    /** stress relaxation for the balls. */
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> ball_stress_relaxation_first_half(ball_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_stress_relaxation_second_half(ball_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> ball_update_contact_density(ball_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_compute_solid_contact_forces(ball_contact);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_ball_displacement("Position", io_environment, ball_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    ball_corrected_configuration.exec();
    /** Initial states output. */
    body_states_recording.writeToFile(0);
    write_ball_displacement.writeToFile(0);
    /** Main loop. */
    int ite = 0;
    Real T0 = 2.0;
    Real end_time = T0;
    Real output_interval = 0.01 * T0;
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
                ball_initialize_timestep.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                ball_update_contact_density.exec();
                ball_compute_solid_contact_forces.exec();
                ball_stress_relaxation_first_half.exec(dt);
                ball_stress_relaxation_second_half.exec(dt);

                ball.updateCellLinkedList();
                ball_contact.updateConfiguration();

                ite++;
                dt = ball_get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            write_ball_displacement.writeToFile(ite);
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
}
