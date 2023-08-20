/**
 * @file 	stl_cartoon_drop_3D.cpp
 * @details a cartoon bodies from .stl file with oobleck material bouncing with a boundary wall
 * @author 	Liezhao Wu, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" /* SPHinXsys Library. */
using namespace SPH;   /* Namespace cite here. */
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/Cartoon.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.5;
Real DH = 0.5;
Real resolution_ref = 0.003;
Real BW = resolution_ref * 4;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-DL, -DL, -BW), Vec3d(DL, DL, DH));
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
Real gravity_g = 2.0;
Real rho0_s = 1.0e3;                                                                                    /* ¦Ñ density. kg/m^3 */
Real Bulk_modulus = 1.09e5;                                                                             /* ¦Ê bulk modulus. Pa */
Real Shear_modulus = 1.12e4;                                                                            /* ¦Ì/G shear modulus. Pa */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* ¦Í Poisson's ratio. */
Real yield_stress = 0.1;                                                                                /* ¦Ò_Y yield stress. Pa */
Real viscosity = 10.0;                                                                                  /* ¦Ç viscosity.  */
Real Herschel_Bulkley_power = 2.8;                                                                    /* viscoplastic material. */
//----------------------------------------------------------------------
//	Body shapes used in the case.
//----------------------------------------------------------------------
class Cartoon : public ComplexShape
{
  public:
    explicit Cartoon(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation_1(0.0, 0.0, 0.75 * DH);
        add<TriangleMeshShapeSTL>(full_path_to_file, translation_1, 0.01);
    }
};
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wall(0.5 * DL, 0.5 * DL, BW / 2.0);
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
    sph_system.setRunParticleRelaxation(false); /* true/false */
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody cartoon(sph_system, makeShared<Cartoon>("Cartoon"));
    cartoon.defineBodyLevelSetShape();
    cartoon.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                     yield_stress, viscosity, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cartoon.generateParticles<ParticleGeneratorReload>(io_environment, cartoon.getName())
        : cartoon.generateParticles<ParticleGeneratorLattice>();

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
        InnerRelation cartoon_inner(cartoon);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for ball.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> random_cartoon_particles(cartoon);
        relax_dynamics::RelaxationStepInner relaxation_step_inner(cartoon_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_states(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, cartoon);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_cartoon_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_states.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
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
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation cartoon_inner(cartoon);
    SelfSurfaceContactRelation cartoon_self_contact(cartoon);
    SurfaceContactRelation cartoon_contact(cartoon_self_contact, {&wall_boundary});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    // initialize a time step
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> initialization_with_gravity(cartoon, gravity_ptr);
    // Corrected configuration for reproducing rigid rotation.
    InteractionWithUpdate<CorrectedConfigurationInner> corrected_configuration(cartoon_inner);
    // Time step size
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_size(cartoon, 0.1);
    // stress relaxation for the ball.
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> stress_relaxation_first_half(cartoon_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(cartoon_inner);
    // Algorithms for solid-solid contacts.
    InteractionDynamics<solid_dynamics::ContactDensitySummation> cartoon_update_contact_density(cartoon_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> cartoon_compute_solid_contact_forces(cartoon_contact);
    InteractionDynamics<solid_dynamics::SelfContactDensitySummation> cartoon_self_contact_density(cartoon_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> cartoon_self_contact_forces(cartoon_self_contact);

    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center(cartoon, Vec3d(1.0, 1.0, 0.0));
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
    // apply initial condition
    corrected_configuration.exec();
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
    Real output_period = 0.01 * T0;
    Real Dt = 0.1 * output_period;
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
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                initialization_with_gravity.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }
                // contact dynamics.
                cartoon_self_contact_density.exec();
                cartoon_self_contact_forces.exec();
                cartoon_update_contact_density.exec();
                cartoon_compute_solid_contact_forces.exec();
                // Stress relaxation and damping.
                stress_relaxation_first_half.exec(dt);
                stress_relaxation_second_half.exec(dt);
                // update particle neighbor relations for contact dynamics
                cartoon.updateCellLinkedList();
                cartoon_self_contact.updateConfiguration();
                cartoon_contact.updateConfiguration();

                ite++;
                dt = computing_time_step_size.exec();
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