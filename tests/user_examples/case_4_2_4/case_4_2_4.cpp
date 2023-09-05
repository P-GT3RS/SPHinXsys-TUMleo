/**
 * @file 	4.2.4 Cartoon Collision with Ground.cpp
 * @brief 	two cartoon characters from .stl file with viscoplastic and oobleck material bouncing with a boundary
 * @author 	Liezhao Wu, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/Xiaohuangren.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.2;
Real DH = 0.4;
Real resolution_ref = 0.0025;
Real BW = resolution_ref * 4;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-DL, -DL, -BW), Vec3d(DL, DL, DH));
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real gravity_g = 1.0;
Real rho0_s = 1.0e3;                                                                                    /* ρ density. kg/m^3 */
Real Bulk_modulus = 1.09e5;                                                                             /* κ bulk modulus. Pa */
Real Shear_modulus = 1.12e4;                                                                            /* μ/G shear modulus. Pa */
Real yield_stress = 0.1;                                                                                /* σ_Y yield stress. Pa */
Real viscous_modulus = 10.0;                                                                            /* η viscosity.  */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa. 此处32487 */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* ν Poisson's ratio. 取值(-1,0.5). 此处0.45 */
Real Herschel_Bulkley_power_1 = 1.0;
Real Herschel_Bulkley_power_2 = 2.8;
//----------------------------------------------------------------------
//	Body shapes used in the case.
//----------------------------------------------------------------------
class XiaohuangrenOne : public ComplexShape
{
  public:
    explicit XiaohuangrenOne(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_1(0.0, -0.25 * DL, DH / 2.0);
        add<TriangleMeshShapeSTL>(full_path_to_file, translation_1, 0.01);
    }
};
class XiaohuangrenTwo : public ComplexShape
{
  public:
    explicit XiaohuangrenTwo(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_2(0.0, 0.5 * DL, DH / 2.0);
        add<TriangleMeshShapeSTL>(full_path_to_file, translation_2, 0.01);
    }
};
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wall(0.5 * DL, DL, BW / 2.0);
        Vecd translation_wall(0.0, 0.0, 0.0);
        add<TriangleMeshShapeBrick>(halfsize_wall, 5, translation_wall);
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
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody xiaohuangren_1(sph_system, makeShared<XiaohuangrenOne>("XiaohuangrenOne"));
    xiaohuangren_1.defineBodyLevelSetShape();
    xiaohuangren_1.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                          yield_stress, viscous_modulus, Herschel_Bulkley_power_1);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? xiaohuangren_1.generateParticles<ParticleGeneratorReload>(io_environment, xiaohuangren_1.getName())
        : xiaohuangren_1.generateParticles<ParticleGeneratorLattice>();

    SolidBody xiaohuangren_2(sph_system, makeShared<XiaohuangrenTwo>("XiaohuangrenTwo"));
    xiaohuangren_2.defineBodyLevelSetShape();
    xiaohuangren_2.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                          yield_stress, viscous_modulus, Herschel_Bulkley_power_2);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? xiaohuangren_2.generateParticles<ParticleGeneratorReload>(io_environment, xiaohuangren_2.getName())
        : xiaohuangren_2.generateParticles<ParticleGeneratorLattice>();

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
        InnerRelation xiaohuangren_1_inner(xiaohuangren_1);
        InnerRelation xiaohuangren_2_inner(xiaohuangren_2);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for the insert body.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> xiaohuangren_1_random_particles(xiaohuangren_1);
        SimpleDynamics<RandomizeParticlePosition> xiaohuangren_2_random_particles(xiaohuangren_2);
        relax_dynamics::RelaxationStepInner relaxation_1_step_inner(xiaohuangren_1_inner);
        relax_dynamics::RelaxationStepInner relaxation_2_step_inner(xiaohuangren_2_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_body_state(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&xiaohuangren_1, &xiaohuangren_2});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        xiaohuangren_1_random_particles.exec(0.25);
        xiaohuangren_2_random_particles.exec(0.25);
        relaxation_1_step_inner.SurfaceBounding().exec();
        relaxation_2_step_inner.SurfaceBounding().exec();
        write_body_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            relaxation_1_step_inner.exec();
            relaxation_2_step_inner.exec();
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
    InnerRelation xiaohuangren_1_inner(xiaohuangren_1);
    SelfSurfaceContactRelation xiaohuangren_1_self_contact(xiaohuangren_1);
    SurfaceContactRelation xiaohuangren_1_contact(xiaohuangren_1_self_contact, {&wall_boundary});
    InnerRelation xiaohuangren_2_inner(xiaohuangren_2);
    SelfSurfaceContactRelation xiaohuangren_2_self_contact(xiaohuangren_2);
    SurfaceContactRelation xiaohuangren_2_contact(xiaohuangren_2_self_contact, {&wall_boundary});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> initialization_1_with_gravity(xiaohuangren_1, gravity_ptr);
    SimpleDynamics<TimeStepInitialization> initialization_2_with_gravity(xiaohuangren_2, gravity_ptr);
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration_1(xiaohuangren_1_inner);
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration_2(xiaohuangren_2_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size_1(xiaohuangren_1, 0.1);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size_2(xiaohuangren_2, 0.1);
    /** stress relaxation for the balls. */
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> stress_1_relaxation_first_half(xiaohuangren_1_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_1_relaxation_second_half(xiaohuangren_1_inner);
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> stress_2_relaxation_first_half(xiaohuangren_2_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_2_relaxation_second_half(xiaohuangren_2_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> xiaohuangren_1_update_contact_density(xiaohuangren_1_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> xiaohuangren_1_compute_solid_contact_forces(xiaohuangren_1_contact);
    InteractionDynamics<solid_dynamics::SelfContactDensitySummation> xiaohuangren_1_self_contact_density(xiaohuangren_1_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> xiaohuangren_1_self_contact_forces(xiaohuangren_1_self_contact);

    InteractionDynamics<solid_dynamics::ContactDensitySummation> xiaohuangren_2_update_contact_density(xiaohuangren_2_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> xiaohuangren_2_compute_solid_contact_forces(xiaohuangren_2_contact);
    InteractionDynamics<solid_dynamics::SelfContactDensitySummation> xiaohuangren_2_self_contact_density(xiaohuangren_2_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> xiaohuangren_2_self_contact_forces(xiaohuangren_2_self_contact);

    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center_2(xiaohuangren_2, Vec3d(1.0, 0.0, 0.0));
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
    corrected_configuration_1.exec();
    corrected_configuration_2.exec();
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
                initialization_1_with_gravity.exec();
                initialization_2_with_gravity.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                xiaohuangren_1_self_contact_density.exec();
                xiaohuangren_1_self_contact_forces.exec();
                xiaohuangren_1_update_contact_density.exec();
                xiaohuangren_1_compute_solid_contact_forces.exec();

                stress_1_relaxation_first_half.exec(dt);
                stress_1_relaxation_second_half.exec(dt);

                xiaohuangren_1.updateCellLinkedList();
                xiaohuangren_1_self_contact.updateConfiguration();
                xiaohuangren_1_contact.updateConfiguration();


                xiaohuangren_2_self_contact_density.exec();
                xiaohuangren_2_self_contact_forces.exec();
                xiaohuangren_2_update_contact_density.exec();
                xiaohuangren_2_compute_solid_contact_forces.exec();

                stress_2_relaxation_first_half.exec(dt);
                constrain_mass_center_2.exec(dt);
                stress_2_relaxation_second_half.exec(dt);

                xiaohuangren_2.updateCellLinkedList();
                xiaohuangren_2_self_contact.updateConfiguration();
                xiaohuangren_2_contact.updateConfiguration();

                ite++;
                Real dt_1 = computing_time_step_size_1.exec();
                Real dt_2 = computing_time_step_size_2.exec();
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
