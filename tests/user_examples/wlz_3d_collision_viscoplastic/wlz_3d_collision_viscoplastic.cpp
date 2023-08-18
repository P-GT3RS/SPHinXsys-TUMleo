/**
 * @file 	ball_cylinder_collision_3D.cpp
 * @details one ball with viscoplastic material collides with two cylinders.
 * @author 	Liezhao Wu, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"  /* SPHinXsys Library. */
using namespace SPH;    /* Namespace cite here. */
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.0025;/* reference resolution. */
Real DL = 0.1;/* box length. */
Real DH = 0.4;/* box height. */
Real ball_radius = resolution_ref * 20.0; /* ball radius. */
Real cylinder_radius = resolution_ref * 4.0; /* cylinder radius. */
Real cylinder_length = DL;
BoundingBox system_domain_bounds(Vec3d(-DL, -DL, -DH), Vec3d(DL, DL, DH));
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real gravity_g = 2.0;
Real rho0_s = 1.0e3;                                                                                    /* ¦Ñ density. kg/m^3 */
Real Bulk_modulus = 1.09e5;                                                                             /* ¦Ê bulk modulus. Pa */
Real Shear_modulus = 1.12e4;                                                                            /* ¦Ì/G shear modulus. Pa */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* ¦Í Poisson's ratio. */
Real yield_stress = 0.1;                                                                                /* ¦Ò_Y yield stress. Pa */
Real viscosity = 10.0;                                                                                  /* ¦Ç viscosity. */
Real Herschel_Bulkley_power = 1.0;                                                                      /* viscoplastic material. */ 
//----------------------------------------------------------------------
//	Geometric shapes
//----------------------------------------------------------------------
class Ball : public ComplexShape
{
  public:
    explicit Ball(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_ball(0.0, 0.0, 0.5 * DH);
        add<TriangleMeshShapeSphere>(ball_radius, 6, translation_ball);
    }
};
class CylinderOne : public ComplexShape
{
  public:
    explicit CylinderOne(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_cylinder_1(0.0, -0.2 * DL, 0.0);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0.0, 0.0), cylinder_radius,
                                       cylinder_length, 25, translation_cylinder_1);
    }
};
class CylinderTwo : public ComplexShape
{
  public:
    explicit CylinderTwo(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_cylinder_2(0.0, 0.2 * DL, 0.0);
        add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0.0, 0.0), cylinder_radius,
                                       cylinder_length, 25, translation_cylinder_2);
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
    sph_system.setRunParticleRelaxation(false); /* true/false */
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true); 
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody cylinder_1(sph_system, makeShared<CylinderOne>("CylinderOne"));
    cylinder_1.defineBodyLevelSetShape();
    cylinder_1.defineParticlesAndMaterial<SolidParticles, Solid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder_1.generateParticles<ParticleGeneratorReload>(io_environment, cylinder_1.getName())
        : cylinder_1.generateParticles<ParticleGeneratorLattice>();

    SolidBody cylinder_2(sph_system, makeShared<CylinderTwo>("CylinderTwo"));
    cylinder_2.defineBodyLevelSetShape();
    cylinder_2.defineParticlesAndMaterial<SolidParticles, Solid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder_2.generateParticles<ParticleGeneratorReload>(io_environment, cylinder_2.getName())
        : cylinder_2.generateParticles<ParticleGeneratorLattice>();

    SolidBody ball(sph_system, makeShared<Ball>("Ball"));
    ball.defineBodyLevelSetShape();
    ball.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
    yield_stress, viscosity, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? ball.generateParticles<ParticleGeneratorReload>(io_environment, ball.getName())
        : ball.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation ball_inner(ball);
        InnerRelation cylinder_1_inner(cylinder_1);
        InnerRelation cylinder_2_inner(cylinder_2);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for ball.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> ball_random_particles(ball);
        SimpleDynamics<RandomizeParticlePosition> cylinder_1_random_particles(cylinder_1);
        SimpleDynamics<RandomizeParticlePosition> cylinder_2_random_particles(cylinder_2);
        relax_dynamics::RelaxationStepInner ball_relaxation_step_inner(ball_inner);
        relax_dynamics::RelaxationStepInner cylinder_1_relaxation_step_inner(cylinder_1_inner);
        relax_dynamics::RelaxationStepInner cylinder_2_relaxation_step_inner(cylinder_2_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_body_state(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&ball, &cylinder_1, &cylinder_2});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        ball_random_particles.exec(0.25);
        cylinder_1_random_particles.exec(0.25);
        cylinder_2_random_particles.exec(0.25);
        write_body_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            ball_relaxation_step_inner.exec();
            cylinder_1_relaxation_step_inner.exec();
            cylinder_2_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_body_state.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation ball_inner(ball);
    SurfaceContactRelation ball_contact(ball, {&cylinder_1, &cylinder_2});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> ball_initialize_timestep(ball, gravity_ptr);
    InteractionWithUpdate<CorrectedConfigurationInner> ball_corrected_configuration(ball_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_get_time_step_size(ball, 0.05);
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
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    ball_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int ite = 0;
    Real T0 = 1.0;
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