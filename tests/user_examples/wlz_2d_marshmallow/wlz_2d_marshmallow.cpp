#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 3.05;                 // cracker length
Real PH = 2.05;                 // height between crackers
Real resolution_ref = 0.02;
Real BW = resolution_ref * 4.0; // cracker width
Real marshmallow_radius = 1.1;
Vec2d marshmallow_center(0.0, 0.0);
BoundingBox system_domain_bounds(Vec2d(-PL - BW, -PH - BW), Vec2d(PL + BW, PH + BW));
Real time_to_apply_velocity = 0.45;
Real cracker_velocity = 0.2;
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
Real physical_viscosity = 1500.0;   
// s'more interior (s'more exterior)
Real rho0_s = 50.0;
Real Bulk_modulus = 1.09e5;
Real Shear_modulus = 80.0;          // 5.0e4
Real yield_stress = 10.0;           // 1000.0
Real viscous_modulus = 16.0;        // 0.1
Real Herschel_Bulkley_power = 0.43; // 1.0
// 2D formula
//Real poisson = (Bulk_modulus - Shear_modulus) / (Bulk_modulus + Shear_modulus);
//Real Youngs_modulus = (4.0 * Shear_modulus * Bulk_modulus) / (Bulk_modulus + Shear_modulus);
// 3D formula
Real poisson = (3.0 * Bulk_modulus - 2.0 * Shear_modulus) / (6.0 * Bulk_modulus + 2.0 * Shear_modulus);
Real Youngs_modulus = (9.0 * Shear_modulus * Bulk_modulus) / (3.0 * Bulk_modulus + Shear_modulus);
//----------------------------------------------------------------------
//	Geometries used in the case.
//----------------------------------------------------------------------
std::vector<Vecd> createLowerCrackerShape()
{
    std::vector<Vecd> lower_cracker_shape;
    lower_cracker_shape.push_back(Vecd(-PL/ 2.0, -PH / 2.0));
    lower_cracker_shape.push_back(Vecd(PL / 2.0, -PH / 2.0));
    lower_cracker_shape.push_back(Vecd(PL / 2.0, -PH / 2.0 - BW));
    lower_cracker_shape.push_back(Vecd(-PL / 2.0, -PH / 2.0 - BW));
    lower_cracker_shape.push_back(Vecd(-PL / 2.0, -PH / 2.0));

    return lower_cracker_shape;
}
std::vector<Vecd> createUpperCrackerShape()
{
    std::vector<Vecd> upper_cracker_shape;
    upper_cracker_shape.push_back(Vecd(-PL / 2.0, PH / 2.0));
    upper_cracker_shape.push_back(Vecd(-PL / 2.0, PH / 2.0 + BW));
    upper_cracker_shape.push_back(Vecd(PL / 2.0, PH / 2.0 + BW));
    upper_cracker_shape.push_back(Vecd(PL / 2.0, PH / 2.0));
    upper_cracker_shape.push_back(Vecd(-PL / 2.0, PH / 2.0));

    return upper_cracker_shape;
}
std::vector<Vecd> createUpperMarshmallowShape()
{
    std::vector<Vecd> upper_marshmallow_shape;
    upper_marshmallow_shape.push_back(Vecd(-PL / 2.0, PH / 2.0));
    upper_marshmallow_shape.push_back(Vecd(-PL / 2.0, marshmallow_radius));
    upper_marshmallow_shape.push_back(Vecd(PL / 2.0, marshmallow_radius));
    upper_marshmallow_shape.push_back(Vecd(PL / 2.0, PH / 2.0));
    upper_marshmallow_shape.push_back(Vecd(-PL / 2.0, PH / 2.0));

    return upper_marshmallow_shape;
}
std::vector<Vecd> createLowerMarshmallowShape()
{
    std::vector<Vecd> lower_marshmallow_shape;
    lower_marshmallow_shape.push_back(Vecd(-PL / 2.0, -PH / 2.0));
    lower_marshmallow_shape.push_back(Vecd(PL / 2.0, -PH / 2.0));
    lower_marshmallow_shape.push_back(Vecd(PL / 2.0, -marshmallow_radius));
    lower_marshmallow_shape.push_back(Vecd(-PL / 2.0, -marshmallow_radius));
    lower_marshmallow_shape.push_back(Vecd(-PL / 2.0, -PH / 2.0));

    return lower_marshmallow_shape;
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
        // 第三项resolution决定圆周的圆滑程度,数字代表由正n边形拟合,数字越大越圆滑.
        // 修改圆形使其两侧接触为线接触.
        multi_polygon_.addACircle(marshmallow_center, marshmallow_radius, 200, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createUpperMarshmallowShape(), ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(createLowerMarshmallowShape(), ShapeBooleanOps::sub);
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
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.setRunParticleRelaxation(false); // false
    sph_system.setReloadParticles(true);        
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

    SolidBody lower_cracker(sph_system, makeShared<LowerCracker>("LowerCracker"));
    lower_cracker.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    lower_cracker.generateParticles<ParticleGeneratorLattice>();

    SolidBody upper_cracker(sph_system, makeShared<UpperCracker>("UpperCracker"));
    upper_cracker.defineParticlesAndMaterial<ElasticSolidParticles, NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    upper_cracker.generateParticles<ParticleGeneratorLattice>();
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
    // SurfaceContactRelation upper_contact(strip_self_contact, { &strip_base_part_two_body });
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    SimpleDynamics<LowerCrackerInitialCondition> initialization_lower_cracker(lower_cracker);
    SimpleDynamics<UpperCrackerInitialCondition> initialization_upper_cracker(upper_cracker);
    // corrected strong configuration
    InteractionWithUpdate<CorrectedConfigurationInner> marshmallow_corrected_configuration(marshmallow_inner);
    // time step size calculation
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> marshmallow_computing_time_step_size(marshmallow);
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
        marshmallow_position_damping(0.5, marshmallow_inner, "Velocity", physical_viscosity);
    //-----------------------------------------------------------------------------
    // outputs
    //-----------------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
    //-----------------------------------------------------------------------------
    //	Setup particle configuration and initial conditions
    //-----------------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    marshmallow_corrected_configuration.exec();
    //-----------------------------------------------------------------------------
    // from here the time stepping begines
    //-----------------------------------------------------------------------------
    write_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int ite = 0;
    Real T0 = 0.5;
    Real End_Time = T0;
    Real D_Time = 0.01 * T0;
    Real Dt = 0.1 * D_Time; /**< Time period for data observing */
    Real dt = 0.0;          // default acoustic time step sizes
    // statistics for computing time
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    // computation loop starts
    while (GlobalStaticVariables::physical_time_ < End_Time)
    {
        Real integration_time = 0.0;
        // integrate time (loop) until the next output time
        while (integration_time < D_Time)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: "<< dt << "\n";
                }

                marshmallow_update_contact_density.exec();
                marshmallow_compute_solid_contact_forces.exec();

                marshmallow_stress_relaxation_first_half.exec(dt);
                constrain_holder0.exec(dt);
                constrain_holder1.exec(dt);
                constrain_mass_center.exec(dt);
                marshmallow_position_damping.exec(dt);
                constrain_holder0.exec(dt);
                constrain_holder1.exec(dt);
                constrain_mass_center.exec(dt);
                marshmallow_stress_relaxation_second_half.exec(dt);

                initialization_upper_cracker.exec(dt);
                initialization_lower_cracker.exec(dt);
                marshmallow.updateCellLinkedList();
                lower_cracker.updateCellLinkedList();
                upper_cracker.updateCellLinkedList();
                // strip_self_contact.updateConfiguration();
                marshmallow_cracker_contact.updateConfiguration();
                // upper_contact.updateConfiguration();

                ite++;
                dt = 0.2 * marshmallow_computing_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
        }
        TickCount t2 = TickCount::now();
        write_body_states.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
