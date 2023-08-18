/**
 * @file 	drop_gravity_3D.cpp
 * @details two balls with viscoplastic and oobleck(shear thickening) material bouncing with a boundary wall
 * @author 	Liezhao Wu, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"                 /* SPHinXsys Library. */
using namespace SPH;                   /* Namespace cite here. */
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real resolution_ref = 0.005;           /* reference resolution. */
Real DL = 0.8;                         /* boundary wall length. */
Real DH = 0.6;                         /* height. */
Real BW = resolution_ref * 4;          /* boundary wall width. */
Real ball_radius = resolution_ref * 12;/* ball radius. */
BoundingBox system_domain_bounds(Vec3d(0.0, 0.0, -BW), Vec3d(0.5 * DL, DL, DH));
//----------------------------------------------------------------------
//	Global parameters on material properties.
//----------------------------------------------------------------------
Real gravity_g = 1.0;
Real rho0_s = 1.0e3;                                                                                    /* �� density. kg/m^3 */
Real Bulk_modulus = 1.09e5;                                                                             /* �� bulk modulus. Pa */
Real Shear_modulus = 1.12e4;                                                                            /* ��/G shear modulus. Pa */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa. */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* �� Poisson's ratio. */
Real yield_stress = 0.1;                                                                                /* ��_Y yield stress. Pa */
Real viscosity = 10.0;                                                                                  /* �� viscosity. */
Real Herschel_Bulkley_power_1 = 1.0;                                                                    /* viscoplastic material. */
Real Herschel_Bulkley_power_2 = 2.8;                                                                    /* oobleck(shear thickening) material. */
//----------------------------------------------------------------------
//	Geometric shapes
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d halfsize_wall(0.5 * DL, DL, BW / 2.0);     
        Vec3d translation_wall(0.0, 0.0, 0.0); 
        add<TriangleMeshShapeBrick>(halfsize_wall, 1, translation_wall);
    }
};
class BallOne : public ComplexShape
{
  public:
    explicit BallOne(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_ball_1(0.25 * DL, 0.3 * DL, 0.75 * DH);
        add<TriangleMeshShapeSphere>(ball_radius, 5, translation_ball_1);
    }
};
class BallTwo : public ComplexShape
{
  public:
    explicit BallTwo(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_ball_2(0.25 * DL, 0.7 * DL, 0.75 * DH);
        add<TriangleMeshShapeSphere>(ball_radius, 5, translation_ball_2);
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
    SolidBody ball_1(sph_system, makeShared<BallOne>("BallOne"));
    ball_1.defineBodyLevelSetShape();
    ball_1.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
    yield_stress, viscosity, Herschel_Bulkley_power_1);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? ball_1.generateParticles<ParticleGeneratorReload>(io_environment, ball_1.getName())
        : ball_1.generateParticles<ParticleGeneratorLattice>();

    SolidBody ball_2(sph_system, makeShared<BallTwo>("BallTwo"));
    ball_2.defineBodyLevelSetShape();
    ball_2.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
    yield_stress, viscosity, Herschel_Bulkley_power_2);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? ball_2.generateParticles<ParticleGeneratorReload>(io_environment, ball_2.getName())
        : ball_2.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation ball_1_inner(ball_1);
        InnerRelation ball_2_inner(ball_2);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for ball.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> ball_1_random_particles(ball_1);
        SimpleDynamics<RandomizeParticlePosition> ball_2_random_particles(ball_2);
        relax_dynamics::RelaxationStepInner ball_1_relaxation_step_inner(ball_1_inner);
        relax_dynamics::RelaxationStepInner ball_2_relaxation_step_inner(ball_2_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_ball_state(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&ball_1, &ball_2});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        ball_1_random_particles.exec(0.25);
        ball_2_random_particles.exec(0.25);
        write_ball_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            ball_1_relaxation_step_inner.exec();
            ball_2_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_ball_state.writeToFile(ite);
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
    InnerRelation ball_1_inner(ball_1);
    SurfaceContactRelation ball_1_contact(ball_1, {&wall_boundary});
    InnerRelation ball_2_inner(ball_2);
    SurfaceContactRelation ball_2_contact(ball_2, {&wall_boundary});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> ball_1_initialize_timestep(ball_1, gravity_ptr);
    SimpleDynamics<TimeStepInitialization> ball_2_initialize_timestep(ball_2, gravity_ptr);
    InteractionWithUpdate<CorrectedConfigurationInner> ball_1_corrected_configuration(ball_1_inner);
    InteractionWithUpdate<CorrectedConfigurationInner> ball_2_corrected_configuration(ball_2_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_1_get_time_step_size(ball_1, 0.05);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> ball_2_get_time_step_size(ball_2, 0.05);
    /** stress relaxation for the ball. */
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> ball_1_stress_relaxation_first_half(ball_1_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_1_stress_relaxation_second_half(ball_1_inner);
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> ball_2_stress_relaxation_first_half(ball_2_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> ball_2_stress_relaxation_second_half(ball_2_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> ball_1_update_contact_density(ball_1_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_1_compute_solid_contact_forces(ball_1_contact);
    InteractionDynamics<solid_dynamics::ContactDensitySummation> ball_2_update_contact_density(ball_2_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> ball_2_compute_solid_contact_forces(ball_2_contact);
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
    ball_1_corrected_configuration.exec();
    ball_2_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int ite = 0;
    Real T0 = 3.0;
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
                ball_1_initialize_timestep.exec();
                ball_2_initialize_timestep.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << " dt: " << dt << "\n";
                }
                ball_1_update_contact_density.exec();
                ball_1_compute_solid_contact_forces.exec();
                ball_1_stress_relaxation_first_half.exec(dt);
                ball_1_stress_relaxation_second_half.exec(dt);

                ball_1.updateCellLinkedList();
                ball_1_contact.updateConfiguration();

                ball_2_update_contact_density.exec();
                ball_2_compute_solid_contact_forces.exec();
                ball_2_stress_relaxation_first_half.exec(dt);
                ball_2_stress_relaxation_second_half.exec(dt);

                ball_2.updateCellLinkedList();
                ball_2_contact.updateConfiguration();

                ite++;
                Real dt_1 = ball_1_get_time_step_size.exec();
                Real dt_2 = ball_2_get_time_step_size.exec();
                dt = SMIN(dt_1, dt_2);
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
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