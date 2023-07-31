/**
 * @file 	particle_relaxation.cpp
 * @brief 	This is the test of using levelset to generate body fitted particles (3D).
 * @details We use this case to test the particle generation and relaxation for a complex geometry.
 *			Before particle generation, we clean the sharp corners of the model.
 * @author 	Yongchuan Yu and Xiangyu Hu
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
Vec3d domain_lower_bound(-4.0, -6.0, -5.0);
Vec3d domain_upper_bound(4.0, 6.0, 5.0);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 12.5;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
StdVec<Vecd> imported_model_observation_location = { Vec3d(0.0, 0.0, 0.0) };
    //----------------------------------------------------------------------
//	Global parameters on material properties
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
//	Geometric shapes
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d halfsize_wall(1.5, 1.5, 0.1);     /* 地面长的一半,宽的一半,厚度的一半. */
        Vec3d translation_wall(0.0, 0.0, -4.5); /* 地面中间层向Z轴负方向移动4.5 */
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall), halfsize_wall);
    }
};
//----------------------------------------------------------------------
//	define a body from the imported model.
//----------------------------------------------------------------------
class SolidBodyFromMesh : public ComplexShape
{
  public:
    explicit SolidBodyFromMesh(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation(0.0, 0.0, 0.0);
        add<TriangleMeshShapeSTL>(full_path_to_file, translation, 1.0);
    }
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, dp_0);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(true); // true
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true); // false
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody imported_model(sph_system, makeShared<SolidBodyFromMesh>("SolidBodyFromMesh"));
    imported_model.defineAdaptation<ParticleRefinementNearSurface>(1.15, 1.0, 3);
    //imported_model.defineBodyLevelSetShape()->correctLevelSetSign()->writeLevelSet(io_environment);
    imported_model.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    imported_model.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
    yield_stress, viscous_modulus, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? imported_model.generateParticles<ParticleGeneratorReload>(io_environment, imported_model.getName())
        : imported_model.generateParticles<ParticleGeneratorMultiResolution>();
    //imported_model.addBodyStateForRecording<Real>("SmoothingLengthRatio");

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    ObserverBody imported_model_observer(sph_system, "ImportedModelObserver");
    imported_model_observer.generateParticles<ObserverParticleGenerator>(imported_model_observation_location);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        AdaptiveInnerRelation imported_model_inner(imported_model);
        //InnerRelation imported_model_inner(imported_model);// Error: reference DynamicCasting class SPH::MultilevelCellLinkedList leads to nullptr!
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for ball.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> imported_model_random_particles(imported_model);
        relax_dynamics::RelaxationStepInner imported_model_relaxation_step_inner(imported_model_inner);
        //SimpleDynamics<relax_dynamics::UpdateSmoothingLengthRatioByShape> update_smoothing_length_ratio(imported_model);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
        //MeshRecordingToPlt cell_linked_list_recording(io_environment, imported_model.getCellLinkedList());
        ReloadParticleIO write_particle_reload_files(io_environment, {&imported_model, &wall_boundary});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        imported_model_random_particles.exec(0.25);
        imported_model_relaxation_step_inner.SurfaceBounding().exec();
        //update_smoothing_length_ratio.exec();
        //imported_model.updateCellLinkedList();
        //cell_linked_list_recording.writeToFile(0);
        write_relaxed_particles.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite_p = 0;
        int relax_step = 1000;
        while (ite_p < relax_step)
        {
            //update_smoothing_length_ratio.exec();
            imported_model_relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                write_relaxed_particles.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    //AdaptiveInnerRelation imported_model_inner(imported_model);
    InnerRelation imported_model_inner(imported_model);//Error: reference DynamicCasting class SPH::MultilevelCellLinkedList leads to nullptr!
    SelfSurfaceContactRelation imported_model_self_contact(imported_model);
    SurfaceContactRelation imported_model_contact(imported_model_self_contact, {&wall_boundary});
    ContactRelation imported_model_observer_contact(imported_model_observer, {&imported_model});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SimpleDynamics<TimeStepInitialization> imported_model_initialize_timestep(imported_model, makeShared<Gravity>(Vec3d(0.0, 0.0, -gravity_g)));
    InteractionWithUpdate<CorrectedConfigurationInner> imported_model_corrected_configuration(imported_model_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> imported_model_get_time_step_size(imported_model, 0.1);
    /** stress relaxation for the ball. */
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> imported_model_stress_relaxation_first_half(imported_model_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> imported_model_stress_relaxation_second_half(imported_model_inner);
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> imported_model_update_contact_density(imported_model_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> imported_model_compute_solid_contact_forces(imported_model_contact);
    InteractionDynamics<solid_dynamics::SelfContactDensitySummation> imported_model_self_contact_density(imported_model_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> imported_model_self_contact_forces(imported_model_self_contact);
    /** Damping for one ball */
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        damping(0.5, imported_model_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_imported_model_displacement("Position", io_environment, imported_model_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    imported_model_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    write_imported_model_displacement.writeToFile(0);
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
                imported_model_initialize_timestep.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << " dt: " << dt << "\n";
                }
                imported_model_self_contact_density.exec();
                imported_model_self_contact_forces.exec();
                imported_model_update_contact_density.exec();
                imported_model_compute_solid_contact_forces.exec();
                imported_model_stress_relaxation_first_half.exec(dt);
                damping.exec(dt);
                imported_model_stress_relaxation_second_half.exec(dt);

                imported_model.updateCellLinkedList();
                imported_model_self_contact.updateConfiguration();
                imported_model_contact.updateConfiguration();

                ite++;
                dt = imported_model_get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                write_imported_model_displacement.writeToFile(ite);
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
    write_imported_model_displacement.testResult();
    return 0;
}
