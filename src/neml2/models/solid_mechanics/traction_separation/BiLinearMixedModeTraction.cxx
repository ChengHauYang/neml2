// Copyright 2024, UChicago Argonne, LLC
// All Rights Reserved
// Software Name: NEML2 -- the New Engineering material Model Library, version 2
// By: Argonne National Laboratory
// OPEN SOURCE LICENSE (MIT)
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#include "neml2/models/solid_mechanics/traction_separation/BiLinearMixedModeTraction.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/tanh.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/imap.h"

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() +=
      " using a bilinear mixed-mode damage law with Benzeggagh-Kenane or power-law criterion.";

  options.set_output("damage") = VariableName(STATE, "damage");
  options.set("damage").doc() = "Scalar damage variable";

  options.set_input("displacement_jump_old") = VariableName(OLD_FORCES, "displacement_jump");
  options.set("displacement_jump_old").doc() = "Old displacement jump (used when lagging)";

  options.set_input("time") = VariableName(FORCES, "t");
  options.set("time").doc() = "Current time (for viscous regularization; time step = t - t_old)";

  options.set_parameter<TensorName<Scalar>>("stiffness");
  options.set("stiffness").doc() = "Penalty elastic stiffness K";

  options.set_parameter<TensorName<Scalar>>("mode_I_fracture_energy");
  options.set("mode_I_fracture_energy").doc() = "Mode I critical energy release rate G_Ic";

  options.set_parameter<TensorName<Scalar>>("mode_II_fracture_energy");
  options.set("mode_II_fracture_energy").doc() = "Mode II critical energy release rate G_IIc";

  options.set_parameter<TensorName<Scalar>>("normal_strength");
  options.set("normal_strength").doc() = "Tensile (normal) strength N";

  options.set_parameter<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Shear strength S";

  options.set_parameter<TensorName<Scalar>>("criterion_exponent");
  options.set("criterion_exponent").doc() = "BK or power-law exponent eta";

  options.set<double>("viscosity") = 0.0;
  options.set("viscosity").doc() = "Viscous regularization coefficient (0 = no regularization)";

  options.set<double>("regularization") = 1e-10;
  options.set("regularization").doc() = "Regularization parameter alpha for the smooth Heaviside";

  options.set<bool>("use_power_law_criterion") = false;
  options.set("use_power_law_criterion").doc() =
      "Use power-law mixed-mode criterion instead of Benzeggagh-Kenane";

  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() =
      "Use the old displacement jump when computing mode mixity ratio beta";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() =
      "Use the old displacement jump when computing the effective displacement jump delta_m";

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _K(declare_parameter<Scalar>("K", "stiffness")),
    _GI_c(declare_parameter<Scalar>("GI_c", "mode_I_fracture_energy")),
    _GII_c(declare_parameter<Scalar>("GII_c", "mode_II_fracture_energy")),
    _N(declare_parameter<Scalar>("N", "normal_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "criterion_exponent")),
    _viscosity(options.get<double>("viscosity")),
    _alpha(options.get<double>("regularization")),
    _use_power_law(options.get<bool>("use_power_law_criterion")),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_disp_jump(options.get<bool>("lag_displacement_jump")),
    _damage(declare_output_variable<Scalar>("damage")),
    _damage_old(declare_input_variable<Scalar>(_damage.name().old())),
    _displacement_jump_old(declare_input_variable<Vec>("displacement_jump_old")),
    _t(declare_input_variable<Scalar>("time")),
    _t_old(declare_input_variable<Scalar>(_t.name().old()))
{
}

void
BiLinearMixedModeTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto & jump = _displacement_jump();
  const auto & jump_old = _displacement_jump_old();

  const auto delta_n = jump(0);
  const auto delta_s1 = jump(1);
  const auto delta_s2 = jump(2);

  // -------------------------------------------------------------------------
  // Regularized Heaviside for normal/tangential traction split (used in traction)
  // H(x) = 0.5*(1 + tanh(x/alpha)),  dH/dx = (1 - tanh^2(x/alpha)) / (2*alpha)
  // -------------------------------------------------------------------------
  const auto th_n = tanh(delta_n / _alpha);
  const auto H = (1.0 + th_n) * 0.5;
  const auto dH = (1.0 - th_n * th_n) / (2.0 * _alpha);

  const auto delta_n_pos = H * delta_n;
  const auto delta_n_neg = delta_n - delta_n_pos;
  const auto zero_s = Scalar::zeros_like(delta_s1);
  const auto delta_active = Vec::fill(delta_n_pos, delta_s1, delta_s2);
  const auto delta_inactive = Vec::fill(delta_n_neg, zero_s, zero_s);

  // -------------------------------------------------------------------------
  // Step 1: Mode mixity (beta = delta_s / delta_n when delta_n > 0)
  // -------------------------------------------------------------------------
  const auto & jump_mix = _lag_mode_mixity ? jump_old : jump;
  const auto dn_mix = jump_mix(0);
  const auto ds1_mix = jump_mix(1);
  const auto ds2_mix = jump_mix(2);
  const auto ds_sq_mix = ds1_mix * ds1_mix + ds2_mix * ds2_mix;
  // regularize tangential norm to avoid sqrt singularity
  const auto ds_mix = sqrt(ds_sq_mix + Scalar(_alpha * _alpha, dn_mix.options()));
  const auto has_opening = dn_mix > Scalar::zeros_like(dn_mix);
  const auto safe_dn = where(has_opening, dn_mix, Scalar::ones_like(dn_mix));
  const auto beta = where(has_opening, ds_mix / safe_dn, Scalar::zeros_like(dn_mix));

  // dbeta/ddelta (Vec) — non-zero only when not lagging and delta_n > 0
  const auto zero_v = Vec::zeros(_displacement_jump.options());
  Vec dbeta_ddelta = zero_v;
  if (dout_din && !_lag_mode_mixity)
  {
    const auto dbd0 = where(has_opening, -ds_mix / (safe_dn * safe_dn), zero_s);
    const auto dbd1 = where(has_opening, ds1_mix / (ds_mix * safe_dn), zero_s);
    const auto dbd2 = where(has_opening, ds2_mix / (ds_mix * safe_dn), zero_s);
    dbeta_ddelta = Vec::fill(dbd0, dbd1, dbd2);
  }

  // -------------------------------------------------------------------------
  // Step 2: Damage initiation displacement jump (delta_init)
  // -------------------------------------------------------------------------
  const auto delta_n0 = _N / _K;   // pure-normal initiation gap
  const auto delta_s0 = _S / _K;   // pure-shear  initiation gap
  const auto delta_mixed_sq = delta_s0 * delta_s0 + beta * beta * delta_n0 * delta_n0;
  const auto delta_mixed = sqrt(delta_mixed_sq);
  const auto delta_init_mixed =
      delta_n0 * delta_s0 * sqrt(1.0 + beta * beta) / delta_mixed;
  const auto delta_init = where(has_opening, delta_init_mixed, delta_s0);

  // ddelta_init/dbeta (Scalar) — only when delta_n > 0
  const auto ddelta_init_dbeta_mixed =
      delta_init_mixed * beta *
      (1.0 / (1.0 + beta * beta) - delta_n0 * delta_n0 / delta_mixed_sq);
  // ddelta_init/ddelta (Vec)
  Vec ddelta_init_ddelta = zero_v;
  if (dout_din && !_lag_mode_mixity)
    ddelta_init_ddelta = where(has_opening,
                               ddelta_init_dbeta_mixed,
                               zero_s) *
                         dbeta_ddelta;

  // -------------------------------------------------------------------------
  // Step 3: Full-degradation displacement jump (delta_final)
  // -------------------------------------------------------------------------
  Scalar delta_final;
  Vec ddelta_final_ddelta = zero_v;

  if (!_use_power_law)
  {
    // Benzeggagh-Kenane criterion
    const auto beta_sq_ratio = beta * beta / (1.0 + beta * beta);
    const auto delta_final_mixed =
        2.0 / (_K * delta_init) *
        (_GI_c + (_GII_c - _GI_c) * pow(beta_sq_ratio, _eta));
    const auto delta_final_shear = 2.0 * std::sqrt(2.0) * _GII_c / _S;
    delta_final = where(has_opening, delta_final_mixed, delta_final_shear);

    if (dout_din && !_lag_mode_mixity)
    {
      // d(delta_final)/d(delta_init) * ddelta_init/ddelta
      const auto ddelta_final_ddelta_init = where(has_opening,
                                                  -delta_final_mixed / delta_init,
                                                  zero_s);
      // d(delta_final)/d(beta) * dbeta/ddelta
      const auto dbeta_sq_ratio_dbeta = 2.0 * beta / ((1.0 + beta * beta) * (1.0 + beta * beta));
      const auto ddelta_final_dbeta =
          where(has_opening,
                2.0 / (_K * delta_init) * (_GII_c - _GI_c) * _eta *
                    pow(beta_sq_ratio, _eta - 1.0) * dbeta_sq_ratio_dbeta,
                zero_s);
      ddelta_final_ddelta =
          ddelta_final_ddelta_init * ddelta_init_ddelta + ddelta_final_dbeta * dbeta_ddelta;
    }
  }
  else
  {
    // Power-law criterion
    const auto Gc_mixed =
        pow(1.0 / _GI_c, _eta) + pow(beta * beta / _GII_c, _eta);
    const auto Gc_mixed_term = pow(Gc_mixed, -1.0 / _eta);
    const auto delta_final_mixed =
        (2.0 + 2.0 * beta * beta) / (_K * delta_init) * Gc_mixed_term;
    const auto delta_final_shear = 2.0 * std::sqrt(2.0) * _GII_c / _S;
    delta_final = where(has_opening, delta_final_mixed, delta_final_shear);

    if (dout_din && !_lag_mode_mixity)
    {
      const auto ddelta_final_ddelta_init = where(has_opening,
                                                  -delta_final_mixed / delta_init,
                                                  zero_s);
      const auto dGc_mixed_dbeta =
          _eta * pow(beta * beta / _GII_c, _eta - 1.0) * (2.0 * beta / _GII_c);
      const auto prefactor = (2.0 + 2.0 * beta * beta) / (_K * delta_init);
      const auto dprefactor_dbeta = 4.0 * beta / (_K * delta_init);
      const auto dGc_term_dbeta =
          (-1.0 / _eta) * pow(Gc_mixed, -1.0 / _eta - 1.0) * dGc_mixed_dbeta;
      const auto ddelta_final_dbeta =
          where(has_opening,
                dprefactor_dbeta * Gc_mixed_term + prefactor * dGc_term_dbeta,
                zero_s);
      ddelta_final_ddelta =
          ddelta_final_ddelta_init * ddelta_init_ddelta + ddelta_final_dbeta * dbeta_ddelta;
    }
  }

  // -------------------------------------------------------------------------
  // Step 4: Effective mixed-mode displacement jump (delta_m)
  // -------------------------------------------------------------------------
  const auto & jump_dm = _lag_disp_jump ? jump_old : jump;
  const auto dn_dm = jump_dm(0);
  const auto ds1_dm = jump_dm(1);
  const auto ds2_dm = jump_dm(2);
  // Regularized Heaviside applied to current jump's normal component for Macaulay bracket
  const auto th_dm = tanh(dn_dm / _alpha);
  const auto H_dm = (1.0 + th_dm) * 0.5;
  const auto dH_dm = (1.0 - th_dm * th_dm) / (2.0 * _alpha);
  const auto dn_pos_dm = H_dm * dn_dm;
  const auto delta_m = sqrt(dn_pos_dm * dn_pos_dm + ds1_dm * ds1_dm + ds2_dm * ds2_dm +
                             Scalar(_alpha * _alpha, dn_dm.options()));

  // ddelta_m/ddelta (Vec) — non-zero only when not lagging
  Vec ddelta_m_ddelta = zero_v;
  if (dout_din && !_lag_disp_jump)
  {
    const auto ddn_pos_ddn = H_dm + dn_dm * dH_dm;
    const auto ddm0 = dn_pos_dm * ddn_pos_ddn / delta_m;
    const auto ddm1 = ds1_dm / delta_m;
    const auto ddm2 = ds2_dm / delta_m;
    ddelta_m_ddelta = Vec::fill(ddm0, ddm1, ddm2);
  }

  // -------------------------------------------------------------------------
  // Step 5: Damage (bilinear law + irreversibility + viscous regularization)
  // -------------------------------------------------------------------------
  const auto & d_old = _damage_old();
  const auto dt = _t() - _t_old();

  // Bilinear damage
  const auto d_linear =
      delta_final * (delta_m - delta_init) / (delta_m * (delta_final - delta_init));
  const auto d_trial = where(delta_m > delta_final,
                             Scalar::ones_like(delta_m),
                             where(delta_m > delta_init, d_linear, Scalar::zeros_like(delta_m)));

  // d(d_linear)/d(delta) via quotient rule
  Vec dd_trial_ddelta = zero_v;
  if (dout_din)
  {
    const auto dm_fi = delta_m * (delta_final - delta_init);
    const auto numerator = delta_final * (delta_m - delta_init);
    const auto d_numerator_ddelta = ddelta_final_ddelta * (delta_m - delta_init) +
                                    delta_final * (ddelta_m_ddelta - ddelta_init_ddelta);
    const auto d_denominator_ddelta =
        ddelta_m_ddelta * (delta_final - delta_init) +
        delta_m * (ddelta_final_ddelta - ddelta_init_ddelta);
    const auto dd_linear_ddelta =
        (d_numerator_ddelta * dm_fi - numerator * d_denominator_ddelta) / (dm_fi * dm_fi);

    // Apply branch mask: non-zero only in softening regime (delta_init < delta_m < delta_final)
    dd_trial_ddelta = where(delta_m > delta_final, zero_v,
                            where(delta_m > delta_init, dd_linear_ddelta, zero_v));
  }

  // Viscous regularization
  const auto visc_dt = Scalar(_viscosity, dt.options()) / dt;
  const auto denom = visc_dt + 1.0;
  const auto d_viscous = (d_trial + visc_dt * d_old) / denom;
  const auto dd_viscous_ddelta = dd_trial_ddelta / denom;

  // Irreversibility: enforce d >= d_old
  const auto advancing = d_viscous >= d_old;
  const auto d_final = where(advancing, d_viscous, d_old);
  const auto dd_final_ddelta = where(advancing, dd_viscous_ddelta, zero_v);

  if (out)
    _damage = d_final;

  // -------------------------------------------------------------------------
  // Step 6: Traction
  // T = (1-d)*K*delta_active + K*delta_inactive
  // -------------------------------------------------------------------------
  if (out)
    _traction = (1.0 - d_final) * _K * delta_active + _K * delta_inactive;

  // -------------------------------------------------------------------------
  // Step 7: Traction Jacobian
  // dT/ddelta = (1-d)*K*d(delta_active)/ddelta + K*d(delta_inactive)/ddelta
  //           - K * outer(delta_active, dd/ddelta)
  // -------------------------------------------------------------------------
  if (dout_din)
    if (_displacement_jump.is_dependent())
    {
      // d(delta_n_pos)/d(delta_n) = H + delta_n * dH
      const auto ddn_pos_ddn = H + delta_n * dH;
      const auto ddn_neg_ddn = 1.0 - ddn_pos_ddn;
      // d(delta_active)/d(delta): diagonal R2
      const auto dact = R2::fill(ddn_pos_ddn, Scalar::ones_like(delta_s1), Scalar::ones_like(delta_s2));
      // d(delta_inactive)/d(delta): diagonal R2 with (ddn_neg, 0, 0)
      const auto dinact = R2::fill(ddn_neg_ddn, zero_s, zero_s);

      _traction.d(_displacement_jump) =
          (1.0 - d_final) * _K * dact + _K * dinact +
          (-_K) * outer(delta_active, dd_final_ddelta);

      _damage.d(_displacement_jump) = dd_final_ddelta;
    }
}
} // namespace neml2
