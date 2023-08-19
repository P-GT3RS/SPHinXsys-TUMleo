/**
 * @file 	marshmallow.cpp
 * @brief 	marshmallow with two crackers
 * @details
 * @author 	Liezhao Wu, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 9.05; // rubber strip length
Real PH = 3.05; // rubber strip height
Real real_PH = 3.0;
Real initial_distance = 0.0;
// particle spacing, at least three particles
Real resolution_ref = real_PH / 60.0;
Real BW = resolution_ref * 4.0; // boundary width
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(0.0, -initial_distance - BW),
                                 Vec2d(PL, PH + initial_distance + BW));
Real time_to_apply_velocity = 0.85;
Real cracker_velocity = 2.0;
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
Real rho0_s = 1.0;                // reference density
Real Youngs_modulus = 1.0e6;      // reference Youngs modulus
Real poisson = 0.45;              // Poisson ratio
Real physical_viscosity = 1000.0; // physical damping
//----------------------------------------------------------------------
//	Geometries used in the case.
//----------------------------------------------------------------------
std::vector<Vecd> createLowerCrackerShape()
{
    std::vector<Vecd> lower_cracker_shape;
    lower_cracker_shape.push_back(Vecd(0.0, -initial_distance - BW));
    lower_cracker_shape.push_back(Vecd(0.0, -initial_distance));
    lower_cracker_shape.push_back(Vecd(PL, -initial_distance));
    lower_cracker_shape.push_back(Vecd(PL, -initial_distance - BW));
    lower_cracker_shape.push_back(Vecd(0.0, -initial_distance - BW));

    return lower_cracker_shape;
}
std::vector<Vecd> createUpperCrackerShape()
{
    std::vector<Vecd> upper_cracker_shape;
    upper_cracker_shape.push_back(Vecd(0.0, PH + initial_distance));
    upper_cracker_shape.push_back(Vecd(0.0, PH + initial_distance + BW));
    upper_cracker_shape.push_back(Vecd(PL, PH + initial_distance + BW));
    upper_cracker_shape.push_back(Vecd(PL, PH + initial_distance));
    upper_cracker_shape.push_back(Vecd(0.0, PH + initial_distance));

    return upper_cracker_shape;
}
std::vector<Vecd> createMarshmallowShape()
{
    std::vector<Vecd> marshmallow_shape;
    marshmallow_shape.push_back(Vecd(0.0, 0.0));
    marshmallow_shape.push_back(Vecd(0.0, PH));
    marshmallow_shape.push_back(Vecd(PL, PH));
    marshmallow_shape.push_back(Vecd(PL, 0.0));
    marshmallow_shape.push_back(Vecd(0.0, 0.0));

    return marshmallow_shape;
}
//----------------------------------------------------------------------
//	Bodies used in the case.
//----------------------------------------------------------------------
class LowerCracker : public MultiPolygonShape
{
  public:
    explicit LowerCracker(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createLowerCrackerShape(), ShapeBooleanOps::add);
    }
};
class UpperCracker : public MultiPolygonShape
{
  public:
    explicit UpperCracker(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createUpperCrackerShape(), ShapeBooleanOps::add);
    }
};
class Marshmallow : public MultiPolygonShape
{
  public:
    explicit Marshmallow(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createMarshmallowShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class LowerCrackerInitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit LowerCrackerInitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

  protected:
    void update(size_t index_i, Real dt)
    {
        Real current_time = GlobalStaticVariables::physical_time_;
        if (current_time < time_to_apply_velocity)
        {
            vel_[index_i][0] = 0.0;
            vel_[index_i][1] = cracker_velocity;
            pos_[index_i] += vel_[index_i] * dt;
        }
    };
};
class UpperCrackerInitialCondition
    : public solid_dynamics::ElasticDynamicsInitialCondition
{
  public:
    explicit UpperCrackerInitialCondition(SPHBody &sph_body)
        : solid_dynamics::ElasticDynamicsInitialCondition(sph_body){};

  protected:
    void update(size_t index_i, Real dt)
    {
        Real current_time = GlobalStaticVariables::physical_time_;
        if (current_time < time_to_apply_velocity)
        {
            vel_[index_i][0] = 0.0;
            vel_[index_i][1] = -cracker_velocity;
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
    BoundaryGeometry0(SPHBody &sph_body, const std::string &body_part_name)
        : BodyPartByParticle(sph_body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry0::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry0(){};

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][0] < 0.5 * PL + 0.3 * resolution_ref &&
            base_particles_.pos_[index_i][0] > 0.5 * PL - 0.3 * resolution_ref)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class BoundaryGeometry1 : public BodyPartByParticle
{
  public:
    BoundaryGeometry1(SPHBody &sph_body, const std::string &body_part_name)
        : BodyPartByParticle(sph_body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry1::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry1(){};

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][1] < 0.5 * PH + 0.3 * resolution_ref &&
            base_particles_.pos_[index_i][1] > 0.5 * PH - 0.3 * resolution_ref)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//------------------------------------------------------------------------------
// Main program starts here.
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody marshmallow(sph_system, makeShared<Marshmallow>("Marshmallow"));
    marshmallow.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    marshmallow.generateParticles<ParticleGeneratorLattice>();

    SolidBody lower_cracker(sph_system, makeShared<LowerCracker>("LowerCracker"));
    lower_cracker.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    lower_cracker.generateParticles<ParticleGeneratorLattice>();

    SolidBody upper_cracker(sph_system, makeShared<UpperCracker>("UpperCracker"));
    upper_cracker.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    upper_cracker.generateParticles<ParticleGeneratorLattice>();
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
    SimpleDynamics<LowerCrackerInitialCondition> initialization_lower_cracker(lower_cracker);
    SimpleDynamics<UpperCrackerInitialCondition> initialization_upper_cracker(upper_cracker);
    // corrected strong configuration
    InteractionWithUpdate<CorrectedConfigurationInner> marshmallow_corrected_configuration(marshmallow_inner);
    // time step size calculation
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> marshmallow_computing_time_step_size(marshmallow,0.5);
    // stress relaxation for the strip
    Dynamics1Level<solid_dynamics::Integration1stHalfKirchhoff> marshmallow_stress_relaxation_first_half(marshmallow_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> marshmallow_stress_relaxation_second_half(marshmallow_inner);

    BoundaryGeometry0 boundary_geometry0(marshmallow, "BoundaryGeometry0");
    SimpleDynamics<solid_dynamics::FixedInAxisDirection> constrain_holder0(boundary_geometry0, Vecd(0.0, 1.0));
    BoundaryGeometry1 boundary_geometry1(marshmallow, "BoundaryGeometry1");
    SimpleDynamics<solid_dynamics::FixedInAxisDirection> constrain_holder1(boundary_geometry1, Vecd(1.0, 0.0));
    SimpleDynamics<solid_dynamics::ConstrainSolidBodyMassCenter> constrain_mass_center(marshmallow, Vecd(1.0, 1.0));
    /** Algorithms for solid-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> marshmallow_update_contact_density(marshmallow_cracker_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> marshmallow_compute_solid_contact_forces(marshmallow_cracker_contact);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        marshmallow_position_damping(0.99, marshmallow_inner, "Velocity", physical_viscosity);
    //-----------------------------------------------------------------------------
    // Output for particle relaxation.
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
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
    body_states_recording.writeToFile(0);
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
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
                }
                marshmallow_update_contact_density.exec();
                marshmallow_compute_solid_contact_forces.exec();

                marshmallow_stress_relaxation_first_half.exec(dt);
                constrain_holder0.exec(dt);
                constrain_holder1.exec(dt);
                // constrain_mass_center.exec(dt);
                marshmallow_position_damping.exec(dt);
                constrain_holder0.exec(dt);
                constrain_holder1.exec(dt);
                // constrain_mass_center.exec(dt);
                marshmallow_stress_relaxation_second_half.exec(dt);

                initialization_upper_cracker.exec(dt);
                initialization_lower_cracker.exec(dt);
                marshmallow.updateCellLinkedList();
                lower_cracker.updateCellLinkedList();
                upper_cracker.updateCellLinkedList();
                marshmallow_cracker_contact.updateConfiguration();

                ite++;
                dt = marshmallow_computing_time_step_size.exec();
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