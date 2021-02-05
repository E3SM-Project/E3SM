#ifndef SCREAM_P3_MICROPHYSICS_HPP
#define SCREAM_P3_MICROPHYSICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "physics/p3/p3_main_impl.hpp"
#include "physics/p3/p3_functions.hpp"
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs 

#include <string>

namespace scream
{

/*
 * The class responsible to handle the atmosphere microphysics
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
 *  Note: for now, scream is only going to accommodate P3 as microphysics
*/

  using namespace p3;
  using P3F          = Functions<Real, DefaultDevice>;
  using Spack        = typename P3F::Spack;
  using Smask        = typename P3F::Smask;
  using Pack         = ekat::Pack<Real,Spack::n>;
  using physics_fun  = scream::physics::Functions<Real, DefaultDevice>;

  using view_1d  = typename P3F::view_1d<Real>;
  using view_2d  = typename P3F::view_2d<Spack>;
  using view_2d_const  = typename P3F::view_2d<const Spack>;
  using sview_2d = typename KokkosTypes<DefaultDevice>::template view_2d<Real>;

class P3Microphysics : public AtmosphereProcess
{
public:
  using field_type       = Field<      Real>;
  using const_field_type = Field<const Real>;

  // Constructors
  P3Microphysics (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Microphysics"; }

  // The communicator used by subcomponent
  const ekat::Comm& get_comm () const { return m_p3_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_p3_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real>& field_repo) const;

  // Get the set of required/computed fields
  const std::set<FieldIdentifier>& get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>& get_computed_fields () const { return m_computed_fields; }

  /*--------------------------------------------------------------------------------------------*/
  // Most individual processes have a pre-processing step that constructs needed variables from
  // the set of fields stored in the field manager.  A structure like this defines those operations,
  // which can then be called during run_imple in the main .cpp code.
  // Structure to handle the local generation of data needed by p3_main in run_impl
  struct run_local_vars {
    run_local_vars(const int ncol, const int npack, view_2d_const pmid_, view_2d T_atm_, view_2d_const ast_, view_2d_const zi_) : 
      m_ncol(ncol),
      m_npack(npack),
      // IN
      pmid(pmid_),
      T_atm(T_atm_),
      ast(ast_),
      zi(zi_),
      // OUT
      exner("exner",ncol,npack),
      th_atm("th_atm",ncol,npack),
      cld_frac_l("cld_frac_l",ncol,npack),
      cld_frac_i("cld_frac_i",ncol,npack),
      cld_frac_r("cld_frac_r",ncol,npack),
      dz("dz",ncol,npack)
    {
      // Nothing else to initialize at the moment.
    };
    KOKKOS_INLINE_FUNCTION
    void operator()(const int icol) const {
      int ipack, ivec;
      for (ipack=0;ipack<m_npack;ipack++) {
        // Exner
        const Spack opmid(pmid(icol,ipack));
        const Smask opmid_mask(!isnan(opmid) and opmid>0.0);
        auto oexner = physics_fun::get_exner(opmid,opmid_mask);
        exner(icol,ipack)  = oexner;
        // Potential temperature
        const Spack oT_atm(T_atm(icol,ipack));
        const Smask oT_atm_mask(!isnan(oT_atm) and oT_atm>0.0);
        auto oth = physics_fun::T_to_th(oT_atm,oexner,oT_atm_mask);
        th_atm(icol,ipack) = oth; 
        // Cloud fraction and dz
        const Spack oast(ast(icol,ipack));
        const Smask oasti_mask(!isnan(oast) and oast>mincld);
        cld_frac_l(icol,ipack).set(oasti_mask,oast);
        cld_frac_i(icol,ipack).set(oasti_mask,oast);
        cld_frac_r(icol,ipack).set(oasti_mask,oast);
        Int kstr = ipack==0 ? 1 : 0;  // If ipack == 0 then we need to skip the first index for rain fraction (i.e. TOM)
        for (int kk=kstr;kk<Spack::n;kk++)
        {
          // Hard-coded max-overlap cloud fraction calculation.  Cycle through the layers from top to bottom and determine if the rain fraction needs to
          // be updated to match the cloud fraction in the layer above.  It is necessary to calculate the location of the layer directly above this one,
          // labeled ipack_m1 and ivec_m1 respectively.  Note, the top layer has no layer above it, which is why we have the kstr index in the loop.
          ivec  = kk % Spack::n;
          Int ipack_m1 = (ipack*Spack::n + kk) / Spack::n;
          Int ivec_m1  = (ipack*Spack::n + kk) % Spack::n;
          cld_frac_r(icol,ipack)[kk] = ast(icol,ipack_m1)[ivec_m1]>cld_frac_r(icol,ipack)[kk] ? 
                                              ast(icol,ipack_m1)[ivec_m1] :
                                              cld_frac_r(icol,ipack)[kk];
          // dz is calculated as the difference between the two layer interfaces.  Note that the lower the index the higher the altitude.
          // We also want to make sure we use the top level index for assignment since dz[0] = zi[0]-zi[1], for example.
          dz(icol,ipack_m1)[ivec_m1] = zi(icol,ipack_m1)[ivec_m1]-zi(icol,ipack)[ivec]; 
        }
        //
      }
    }

    int m_ncol, m_npack;
    Real mincld = 0.0001;  // TODO: These should be stored somewhere as more universal constants.  Or maybe in the P3 class hpp
    view_2d_const pmid;
    view_2d       T_atm;
    view_2d_const ast;
    view_2d_const zi;
    view_2d       exner;
    view_2d       th_atm;
    view_2d       cld_frac_l;
    view_2d       cld_frac_i;
    view_2d       cld_frac_r;
    view_2d       dz;
  }; // run_local_vars
  /* --------------------------------------------------------------------------------------------*/

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const util::TimeStamp& t0);
  void run_impl        (const Real dt);
  void finalize_impl   ();

  // Setting the fields in the atmospheric process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

  std::set<FieldIdentifier> m_required_fields;
  std::set<FieldIdentifier> m_computed_fields;

  std::map<std::string,const_field_type>  m_p3_fields_in;
  std::map<std::string,field_type>        m_p3_fields_out;

  template<typename T>
  using view_type = field_type::view_type<T*>;

  template<typename T>
  using host_view_type = field_type::get_view_type<view_type<T>,Host>;

  using host_view_in_type   = host_view_type<const_field_type::RT>;
  using host_view_out_type  = host_view_type<      field_type::RT>;

  std::map<std::string,host_view_in_type>   m_p3_host_views_in;
  std::map<std::string,host_view_out_type>  m_p3_host_views_out;

  std::map<std::string,const Real*>  m_raw_ptrs_in;
  std::map<std::string,Real*>        m_raw_ptrs_out;

  // Used to init some fields. For now, only needed for stand-alone p3 runs
  std::shared_ptr<FieldInitializer>  m_initializer;

  util::TimeStamp   m_current_ts;
  ekat::Comm              m_p3_comm;

  ekat::ParameterList     m_p3_params;

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_levs;

  // Store the structures for each arguement to p3_main;
  P3F::P3PrognosticState prog_state;
  P3F::P3DiagnosticInputs diag_inputs;
  P3F::P3DiagnosticOutputs diag_outputs;
  P3F::P3HistoryOnly history_only;
  P3F::P3Infrastructure infrastructure;
  // Iteration count is internal to P3 and keeps track of the number of times p3_main has been called.
  // infrastructure.it is passed as an arguement to p3_main and is used for identifying which iteration an error occurs. 

}; // class P3Microphysics

} // namespace scream

#endif // SCREAM_P3_MICROPHYSICS_HPP
