/**
 * @file 	shaving_cream_drop_3D.cpp
 * @brief
 * @details
 * @author 	Liezhao Wu, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/shaving_foam.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 0.2; 
Real DH = 0.5;
Real resolution_ref = 0.0025;  /* reference resolution. */
Real BW = resolution_ref * 4; /* width for BCs. */
BoundingBox system_domain_bounds(Vec3d(-DL, -DL, 0.0), Vec3d(DL, DL, DH));

Real start_time = 0.0;
Real end_time = 5.0;
Vec3d displacement(0.2, 0.0, 0.0);
int frequency = 5;
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
// shaving cream
Real rho0_s = 77.7;                                                                                     /* ¦Ñ density. kg/m^3 */
Real Bulk_modulus = 1.09e5;                                                                             /* ¦Ê bulk modulus. Pa */
Real Shear_modulus = 290.0;                                                                             /* ¦Ì/G shear modulus. Pa */
Real yield_stress = 31.9;                                                                               /* ¦Ò_Y yield stress. Pa */
Real viscous_modulus = 27.2;                                                                            /* ¦Ç viscosity. */
Real Herschel_Bulkley_power = 0.22;                                                                     /* h Herschel_Bulkley_power. */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* Poisson's ratio. */
Real gravity_g = 9.8;
//----------------------------------------------------------------------
//	Geometries used in the case.
//----------------------------------------------------------------------
class Cream : public ComplexShape
{
  public:
    explicit Cream(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d translation(0.0, 0.0, DH / 2.0);
        add<TriangleMeshShapeSTL>(full_path_to_file, translation, 0.02);
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
    /* running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false); // true
    /* starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true); // false
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody cream(sph_system, makeShared<Cream>("Cream"));
    cream.defineBodyLevelSetShape();
    cream.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                 yield_stress, viscous_modulus, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cream.generateParticles<ParticleGeneratorReload>(io_environment, cream.getName())
        : cream.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation cream_inner(cream);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> cream_random_particles(cream);
        relax_dynamics::RelaxationStepInner cream_relaxation_step_inner(cream_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_cream_state(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, cream);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        cream_random_particles.exec(0.25);
        write_cream_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            cream_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_cream_state.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of cream particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation cream_inner(cream);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> cream_initialize_gravity(cream, gravity_ptr);
    InteractionWithUpdate<CorrectedConfigurationInner> cream_corrected_configuration(cream_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> cream_get_time_step_size(cream, 0.5);
    /** stress relaxation for the foam. */
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> cream_stress_relaxation_first_half(cream_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> cream_stress_relaxation_second_half(cream_inner);

    Vec3d halfsize(DL, DL, 0.025);
    Vec3d translation(0.0, 0.0, DH / 2.0);
    BodyRegionByParticle wall_boundary(cream, makeShared<TriangleMeshShapeBrick>(halfsize, 10, translation));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> wall_constraint(wall_boundary);
    SimpleDynamics<solid_dynamics::ForthAndBackMotionBodyPart> wall_initialize_motion(wall_boundary, start_time, end_time, displacement, frequency);
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
    cream_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int ite = 0;
    Real T0 = 5.0;
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
                cream_initialize_gravity.exec();
                //wall_initialize_motion.exec(dt);
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                cream_stress_relaxation_first_half.exec(dt);
                wall_constraint.exec(dt);
                cream_stress_relaxation_second_half.exec(dt);

                cream.updateCellLinkedList();

                ite++;
                dt = cream_get_time_step_size.exec();
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
