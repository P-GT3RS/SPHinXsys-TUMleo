#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.0;                /* wall length. */
Real DH = 0.25;               /* wall thickness. */
Real resolution_ref = 0.02;   /* reference resolution. */
Real BW = resolution_ref * 4; /* width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-DL / 2.0 - BW, -DL), Vec2d(DL / 2.0 + BW, DH + BW));
Real foam_radius = 0.5;
Vec2d foam_center(0.0, DH - foam_radius);
StdVec<Vecd> observation_location = {foam_center}; /* observer location. */
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
// shaving cream
Real rho0_s = 77.7;                 /* ¦Ñ density. kg/m^3 */
Real Bulk_modulus = 1.09e5;         /* ¦Ê bulk modulus. Pa */
Real Shear_modulus = 290;           /* ¦Ì/G shear modulus. Pa */
Real yield_stress = 31.9;           /* ¦Ò_Y yield stress. Pa */
Real viscous_modulus = 27.2;        /* ¦Ç viscosity. */
Real Herschel_Bulkley_power = 0.22; /* h Herschel_Bulkley_power. */
//  2D formula
// Real Youngs_modulus = (4.0 * Shear_modulus * Bulk_modulus) / (Bulk_modulus + Shear_modulus);
// Real poisson = (Bulk_modulus - Shear_modulus) / (Bulk_modulus + Shear_modulus);
//  3D formula
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* Poisson's ratio. */
Real gravity_g = 9.8;
std::vector<Vecd> wall_shape{Vecd(-DL / 2.0, 0.0), Vecd(-DL / 2.0, DH),
                             Vecd(DL / 2.0, DH), Vecd(DL / 2.0, 0.0), Vecd(-DL / 2.0, 0.0)};
std::vector<Vecd> foam_base_shape{Vecd(-foam_radius, DH - foam_radius), Vecd(-foam_radius, 0.0),
                                  Vecd(foam_radius, 0.0), Vecd(foam_radius, DH - foam_radius), Vecd(-foam_radius, DH - foam_radius)};
//----------------------------------------------------------------------
//	Geometric shapes
//----------------------------------------------------------------------
class Foam : public MultiPolygonShape
{
  public:
    explicit Foam(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(foam_base_shape, ShapeBooleanOps::add);
        multi_polygon_.addACircle(foam_center, foam_radius, 100, ShapeBooleanOps::add);
    }
};
MultiPolygon createFoamConstrainShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(wall_shape, ShapeBooleanOps::add);
    return multi_polygon;
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
    SolidBody foam(sph_system, makeShared<Foam>("Foam"));
    foam.defineBodyLevelSetShape();
    foam.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
    yield_stress, viscous_modulus, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? foam.generateParticles<ParticleGeneratorReload>(io_environment, foam.getName())
        : foam.generateParticles<ParticleGeneratorLattice>();

    ObserverBody foam_observer(sph_system, "FoamObserver");
    foam_observer.generateParticles<ObserverParticleGenerator>(observation_location);

    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation foam_inner(foam);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> foam_random_particles(foam);
        relax_dynamics::RelaxationStepInner foam_relaxation_step_inner(foam_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_foam_state(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&foam});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        foam_random_particles.exec(0.25);
        write_foam_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            foam_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_foam_state.writeToFile(ite);
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
    InnerRelation foam_inner(foam);
    ContactRelation foam_observer_contact(foam_observer, {&foam});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> foam_initialize_timestep(foam, gravity_ptr);
    InteractionWithUpdate<CorrectedConfigurationInner> foam_corrected_configuration(foam_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> foam_get_time_step_size(foam, 0.2);
    /** stress relaxation for the foam. */
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> foam_stress_relaxation_first_half(foam_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> foam_stress_relaxation_second_half(foam_inner);
    /** Algorithms for solid-solid contact. */
    // clamping a solid body part. This is softer than a direct constraint
    BodyRegionByParticle foam_base(foam, makeShared<MultiPolygonShape>(createFoamConstrainShape()));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_foam_base(foam_base);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        foam_displacement_recording("Position", io_environment, foam_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    foam_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	Initial states output.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile(0);
    foam_displacement_recording.writeToFile(0);
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
                foam_initialize_timestep.exec();
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                foam_stress_relaxation_first_half.exec(dt);
                constraint_foam_base.exec();
                foam_stress_relaxation_second_half.exec(dt);

                foam.updateCellLinkedList();

                ite++;
                dt = foam_get_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                foam_displacement_recording.writeToFile(ite);
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

    foam_displacement_recording.testResult();

    return 0;
}
