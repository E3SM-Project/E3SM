#include <mam4xx/mam4.hpp>
#include <mam_coupling.hpp>

namespace scream::impl {

// performs gas phase chemistry calculations on a single level of a single
// atmospheric column
KOKKOS_INLINE_FUNCTION
void gas_phase_chemistry(Real q[mam_coupling::gas_pcnst()]) {
  constexpr int gas_pcnst = mam_coupling::gas_pcnst();

  // fetch the zenith angle (not its cosine!) in degrees for this column.
  // FIXME: RRTMGP does this in Fortran--shall we also for now?
  Real zenith = 0.0; // [deg]

  // Xform geopotential height from m to km and pressure from Pa to mb
  zsurf = rga * phis;
  zintr = m2km * zi;
  zmid = m2km * (zm + zsurf);
  zint = m2km * (zi + zsurf);

  zint(:ncol,pver+1) = m2km * (zi(:ncol,pver+1) + zsurf(:ncol))
  zintr(:ncol,pver+1)= m2km *  zi(:ncol,pver+1)

  // ... map incoming concentrations to working array
  for (int mm = 0; mm < pcnst; ++mm) {
    int nn = map2chm[mm];
    if (nn > 0) {
      mmr(:ncol,:,nn) = state_q(:ncol,:,mm);
    }
  }

  // ... set atmosphere mean mass
  set_mean_mass(ncol, mbar);

  // ... xform from mmr to vmr
  mmr2vmr(mmr, vmr, mbar, ncol);


  // ... xform water vapor from mmr to vmr and set upper bndy values
  qh2o(:ncol,:) = state_q(:ncol,:,1)
  h2o_to_vmr(state_q(:,:,1), h2ovmr, mbar, ncol);

  // ... set the "invariants"
  setinv(invariants, tfld, h2ovmr, vmr, pmid, ncol, lchnk, pbuf);

  // ... set the column densities at the upper boundary
  set_ub_col(col_delta, vmr, invariants, pdel, ncol, lchnk);

  // ... set rates for "tabular" and user specified reactions
  setrxt(reaction_rates, tfld, ncol);

  // compute the relative humidity
  qsat(tfld, pmid, satv, satq);
  Real relhum = 0.622 * h2ovmr / satq;
  relhum = min_max_bound(0, 1, relhum);

  cwat(:ncol,:pver) = cldw(:ncol,:pver)

  usrrxt(reaction_rates, tfld, invariants, invariants(:,:,indexm), ncol);
  adjrxt(reaction_rates, invariants, invariants(1,1,indexm), ncol);

  //===================================
  // Photolysis rates at time = t(n+1)
  //===================================

  // ... set the column densities
  setcol(col_delta, col_dens);

  // ... calculate the photodissociation rates
  Real esfact = 1.0;
  shr_orb_decl(calday, eccen, mvelpp, lambm0, obliqr, delta, esfact);

    // ... look up the photolysis rates from table
    table_photo(reaction_rates, pmid, pdel, tfld, col_dens, zen_angle,
      asdir, cwat, cldfr, esfact, ncol);

    // ... compute the extraneous frcing at time = t(n+1)
    setext(extfrc, lchnk, ncol, zintr);

    for (int mm = 0; mm < extcnt; ++mm) {
      if (mm != synoz_ndx) {
        for (int kk = 0; kk < pver; ++kk) {
          extfrc(:ncol,kk,mm) = extfrc(:ncol,kk,mm) / invariants(:ncol,kk,indexm)
        }
      }
    }

    // ... Form the washout rates
    sethet(het_rates, pmid, zmid, phis, tfld,
      cmfdqr, prain, nevapr, delt, invariants(:,:,indexm),
      vmr, ncol, lchnk);

    ltrop_sol(:ncol) = 0 // apply solver to all levels

    // save h2so4 before gas phase chem (for later new particle nucleation)
    del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4)

    vmr0(:ncol,:,:) = vmr(:ncol,:,:) // mixing ratios before chemistry changes

    //========================================
    // ... Call the class solution algorithms
    //========================================

    //	... Solve for "Implicit" species
    imp_sol(vmr, reaction_rates, het_rates, extfrc, delt,
      invariants(1,1,indexm), ncol, lchnk, ltrop_sol(:ncol));

    if (convproc_do_aer) {
      vmr2mmr(vmr, mmr_new, mbar, ncol);
      mmr_new(:ncol,:,:) = 0.5_r8*( mmr(:ncol,:,:)+mmr_new(:ncol,:,:) );
      //RCE - mmr_new = average of mmr values before and after imp_sol
      het_diags(het_rates(:ncol,:,:), mmr_new(:ncol,:,:), pdel(:ncol,:), lchnk, ncol );
    }

  // save h2so4 change by gas phase chem (for later new particle nucleation)
  if (ndx_h2so4 > 0) {
     del_h2so4_gasprod(1:ncol,:) = vmr(1:ncol,:,ndx_h2so4) - del_h2so4_gasprod(1:ncol,:)
  }
}

} // namespace scream::impl
