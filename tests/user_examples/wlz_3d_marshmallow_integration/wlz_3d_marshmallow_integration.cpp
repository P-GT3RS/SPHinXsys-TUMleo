/**
 * @file 	marshmallow.cpp
 * @brief 	marshmallow with two crackers
 * @details  暂时不稳定
 * @author 	Liezhao Wu, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real resolution_ref = 0.002;
Real DL = 0.05; // cracker length
Real DH = 0.05; // height between crackers
Real BW = resolution_ref * 4.0;
Real marshmallow_radius = resolution_ref * 10.0;
BoundingBox system_domain_bounds(Vec3d(-DL, -DL, -DL),
                                 Vec3d(DL, DL, DL));
Real time_to_apply_boundary_velocity = 1.0;
Real velocity = 0.05;
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
// s'more interior
Real rho0_s = 50.0;                                                                                     /* ρ density. kg/m^3 */
Real Bulk_modulus = 1.09e5;                                                                             /* κ bulk modulus. Pa */
Real Shear_modulus = 80.0;                                                                              /* μ/G shear modulus. Pa */
Real yield_stress = 10.0;                                                                               /* σ_Y yield stress. Pa */
Real viscous_modulus = 16.0;                                                                            /* η viscosity.  */
Real Herschel_Bulkley_power = 0.43;                                                                     /* h Herschel_Bulkley_power. */
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);      /* E Young's modulus. Pa */
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus); /* ν Poisson's ratio. 取值(-1,0.5). */
//----------------------------------------------------------------------
//	Bodies used in the case.
//----------------------------------------------------------------------
class Marshmallow : public ComplexShape
{
  public:
    explicit Marshmallow(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation_marshmallow(0.0, 0.0, 0.0);
        //add<TriangleMeshShapeCylinder>(SimTK::UnitVec3(1.0, 0.0, 0.0),
        //                               marshmallow_radius, DL / 2.0, 20, translation_marshmallow);
        add<TriangleMeshShapeSphere>(marshmallow_radius, 5, translation_marshmallow);
    }
};
class LowerCracker : public ComplexShape
{
  public:
    explicit LowerCracker(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_lower_cracker(DL, DL, 0.5 * BW);
        Vecd translation_lower_cracker(0.0, 0.0, -DH / 2.0);
        add<TriangleMeshShapeBrick>(halfsize_lower_cracker, 20, translation_lower_cracker);
    }
};
class UpperCracker : public ComplexShape
{
  public:
    explicit UpperCracker(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_upper_cracker(DL, DL, 0.5 * BW);
        Vecd translation_upper_cracker(0.0, 0.0, DH / 2.0);
        add<TriangleMeshShapeBrick>(halfsize_upper_cracker, 20, translation_upper_cracker);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class LowerCrackerInitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit LowerCrackerInitialCondition(SPHBody &solid_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(solid_body){};

  protected:
    void update(size_t index_i, Real dt)
    {
        Real current_time = GlobalStaticVariables::physical_time_;
        if (current_time < time_to_apply_boundary_velocity)
        {
            vel_[index_i][2] = velocity;
            pos_[index_i] += vel_[index_i] * dt;
        }
    };
};
class UpperCrackerInitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit UpperCrackerInitialCondition(SPHBody &solid_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(solid_body){};

  protected:
    void update(size_t index_i, Real dt)
    {
        Real current_time = GlobalStaticVariables::physical_time_;
        if (current_time < time_to_apply_boundary_velocity)
        {
            vel_[index_i][2] = -velocity;
            pos_[index_i] += vel_[index_i] * dt;
        }
    };
};
//----------------------------------------------------------------------
//	Body parts usually for impose constraints.
//----------------------------------------------------------------------
class BoundaryGeometry0 : public BodyPartByParticle
{
  public:
    BoundaryGeometry0(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry0::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry0(){};

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][0] < 0.5 * DL + 0.3 * resolution_ref && base_particles_.pos_[index_i][0] > 0.5 * DL - 0.3 * resolution_ref)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class BoundaryGeometry1 : public BodyPartByParticle
{
  public:
    BoundaryGeometry1(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry1::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry1(){};

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][1] < 0.5 * DL + 0.3 * resolution_ref && base_particles_.pos_[index_i][1] > 0.5 * DL - 0.3 * resolution_ref)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class BoundaryGeometry2 : public BodyPartByParticle
{
  public:
    BoundaryGeometry2(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry2::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry2(){};

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][2] < 0.5 * DH + 0.3 * resolution_ref && base_particles_.pos_[index_i][2] > 0.5 * DH - 0.3 * resolution_ref)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false); // true
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true); // false
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody marshmallow(sph_system, makeShared<Marshmallow>("Marshmallow"));
    marshmallow.defineBodyLevelSetShape();
    marshmallow.defineParticlesAndMaterial<ElasticSolidParticles, ViscousPlasticSolid>(rho0_s, Youngs_modulus, poisson,
                                                                                       yield_stress, viscous_modulus, Herschel_Bulkley_power);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? marshmallow.generateParticles<ParticleGeneratorReload>(io_environment, marshmallow.getName())
        : marshmallow.generateParticles<ParticleGeneratorLattice>();

    SolidBody upper_cracker(sph_system, makeShared<UpperCracker>("UpperCracker"));
    upper_cracker.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    upper_cracker.generateParticles<ParticleGeneratorLattice>();

    SolidBody lower_cracker(sph_system, makeShared<LowerCracker>("LowerCracker"));
    lower_cracker.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    lower_cracker.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation marshmallow_inner(marshmallow);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> marshmallow_random_particles(marshmallow);
        relax_dynamics::RelaxationStepInner marshmallow_relaxation_step_inner(marshmallow_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_marshmallow_state(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, marshmallow);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        marshmallow_random_particles.exec(0.25);
        write_marshmallow_state.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            marshmallow_relaxation_step_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_marshmallow_state.writeToFile(ite);
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
    InnerRelation marshmallow_inner(marshmallow);
    SelfSurfaceContactRelation marshmallow_self_contact(marshmallow);
    SurfaceContactRelation marshmallow_cracker_contact(marshmallow_self_contact, {&upper_cracker, &lower_cracker});
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    SimpleDynamics<UpperCrackerInitialCondition> initialization_upper_cracker(upper_cracker);
    SimpleDynamics<LowerCrackerInitialCondition> initialization_lower_cracker(lower_cracker);
    // corrected strong configuration
    InteractionWithUpdate<CorrectedConfigurationInner> marshmallow_corrected_configuration(marshmallow_inner);
    // time step size calculation
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> computing_time_step_marshmallow(marshmallow, 0.2);
    // stress relaxation for the strip
    Dynamics1Level<solid_dynamics::PlasticIntegration1stHalf> stress_relaxation_first_half(marshmallow_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(marshmallow_inner);

    BoundaryGeometry0 boundary_geometry0(marshmallow, "BoundaryGeometry0");
    SimpleDynamics<solid_dynamics::FixedInAxisDirection>
        constrain_holder0(boundary_geometry0, Vec3d(1.0, 0.0, 0.0));
    BoundaryGeometry1 boundary_geometry1(marshmallow, "BoundaryGeometry1");
    SimpleDynamics<solid_dynamics::FixedInAxisDirection>
        constrain_holder1(boundary_geometry1, Vec3d(0.0, 1.0, 0.0));
    BoundaryGeometry2 boundary_geometry2(marshmallow, "BoundaryGeometry2");
    SimpleDynamics<solid_dynamics::FixedInAxisDirection>
        constrain_holder2(boundary_geometry2, Vec3d(0.0, 0.0, 1.0));
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center(marshmallow, Vec3d(1.0, 1.0, 1.0));
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> marshmallow_update_contact_density(marshmallow_cracker_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> marshmallow_compute_solid_contact_forces(marshmallow_cracker_contact);
    //-----------------------------------------------------------------------------
    // Output for particle relaxation.
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp write_marshmallow_states(io_environment, sph_system.real_bodies_);
    //-----------------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //-----------------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    marshmallow_corrected_configuration.exec();
    //-----------------------------------------------------------------------------
    // from here the time stepping begines
    //-----------------------------------------------------------------------------
    write_marshmallow_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int ite = 0;
    Real T0 = 1.0;
    Real end_Time = T0;
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
    while (GlobalStaticVariables::physical_time_ < end_Time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "
                              << dt << "\n";
                }
                initialization_upper_cracker.exec(dt);
                initialization_lower_cracker.exec(dt);

                marshmallow_update_contact_density.exec();
                marshmallow_compute_solid_contact_forces.exec();

                stress_relaxation_first_half.exec(dt);
                constrain_holder0.exec(dt);
                constrain_holder1.exec(dt);
                constrain_holder2.exec(dt);
                constrain_mass_center.exec(dt);
                constrain_holder0.exec(dt);
                constrain_holder1.exec(dt);
                constrain_holder2.exec(dt);
                stress_relaxation_second_half.exec(dt);

                marshmallow.updateCellLinkedList();
                upper_cracker.updateCellLinkedList();
                lower_cracker.updateCellLinkedList();
                marshmallow_cracker_contact.updateConfiguration();

                ite++;
                dt = computing_time_step_marshmallow.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
        }
        TickCount t2 = TickCount::now();
        write_marshmallow_states.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
