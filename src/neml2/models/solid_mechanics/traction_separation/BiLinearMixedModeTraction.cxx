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
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/cos.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sin.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/where.h"

#include <cmath>

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() +=
      " Mixed-mode bilinear traction-separation law with damage (BK or POWER_LAW criterion).";

  options.set_parameter<TensorName<Scalar>>("penalty_stiffness");
  options.set("penalty_stiffness").doc() = "Penalty stiffness K";

  options.set_parameter<TensorName<Scalar>>("GIc");
  options.set("GIc").doc() = "Critical energy release rate in mode I (normal)";

  options.set_parameter<TensorName<Scalar>>("GIIc");
  options.set("GIIc").doc() = "Critical energy release rate in mode II (shear)";

  options.set_parameter<TensorName<Scalar>>("normal_strength");
  options.set("normal_strength").doc() = "Normal (tensile) strength N";

  options.set_parameter<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Shear strength S";

  options.set_parameter<TensorName<Scalar>>("eta");
  options.set("eta").doc() = "Power-law exponent for mixed-mode criterion";

  options.set_parameter<TensorName<Scalar>>("viscosity");
  options.set("viscosity").doc() = "Viscous regularization coefficient (0 disables)";

  options.set_parameter<TensorName<Scalar>>("alpha");
  options.set("alpha").doc() = "Smoothing length for regularized Heaviside";

  options.set<std::string>("mixed_mode_criterion") = "BK";
  options.set("mixed_mode_criterion").doc() = "Mixed-mode criterion: 'BK' or 'POWER_LAW'";

  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() =
      "Use old displacement jump when computing mode-mixity ratio beta";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() =
      "Use old displacement jump when computing effective displacement jump delta_m";

  options.set_input("dt") = VariableName(STATE, "dt");
  options.set("dt").doc() = "Time step size (used for viscous regularization)";

  options.set_input("damage_old") = VariableName(OLD_STATE, "damage");
  options.set("damage_old").doc() = "Old damage value for irreversibility enforcement";

  options.set_input("jump_old") = VariableName(OLD_STATE, "displacement_jump");
  options.set("jump_old").doc() =
      "Old displacement jump (required when lag_mode_mixity or lag_displacement_jump is true)";

  options.set_output("damage") = VariableName(STATE, "damage");
  options.set("damage").doc() = "Scalar damage variable d in [0, 1]";

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _K(declare_parameter<Scalar>("K", "penalty_stiffness")),
    _GIc(declare_parameter<Scalar>("GIc", "GIc")),
    _GIIc(declare_parameter<Scalar>("GIIc", "GIIc")),
    _N(declare_parameter<Scalar>("N", "normal_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "eta")),
    _visc(declare_parameter<Scalar>("visc", "viscosity")),
    _alpha(declare_parameter<Scalar>("alpha", "alpha")),
    _criterion(options.get<std::string>("mixed_mode_criterion")),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_disp_jump(options.get<bool>("lag_displacement_jump")),
    _dt(declare_input_variable<Scalar>("dt")),
    _jump_old((_lag_mode_mixity || _lag_disp_jump) ? &declare_input_variable<Vec>("jump_old")
                                                   : nullptr),
    _damage_old(declare_input_variable<Scalar>("damage_old")),
    _damage(declare_output_variable<Scalar>("damage"))
{
}

Scalar
BiLinearMixedModeTraction::reg_heaviside(const Scalar & x) const
{
  const Scalar zero = Scalar::zeros_like(x);
  const Scalar one = Scalar::ones_like(x);
  const Scalar mid = 0.5 * (1.0 + sin(M_PI * x / (2.0 * _alpha)));
  return where(x <= -_alpha, zero, where(x < _alpha, mid, one));
}

Scalar
BiLinearMixedModeTraction::reg_heaviside_deriv(const Scalar & x) const
{
  const Scalar zero = Scalar::zeros_like(x);
  const Scalar mid = 0.25 * (M_PI / _alpha) * cos(M_PI * x / (2.0 * _alpha));
  return where((x < _alpha) && (x > -_alpha), mid, zero);
}

void
BiLinearMixedModeTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto g = _jump();
  const auto g_old = _jump_old ? (*_jump_old)() : g;

  // Which jump is used for mixity/effective computations
  const auto delta_mix = _lag_mode_mixity ? g_old : g;
  const auto g_eff = _lag_disp_jump ? g_old : g;

  // ----- Step 1: Mode mixity ratio beta -----
  const Scalar dn_mix = delta_mix(0);
  const Scalar ds_mix = sqrt(delta_mix(1) * delta_mix(1) + delta_mix(2) * delta_mix(2));

  const Scalar beta = where(dn_mix > 0, ds_mix / dn_mix, Scalar::zeros_like(dn_mix));

  Vec dbeta_ddelta = Vec::zeros_like(g);
  if (!_lag_mode_mixity)
  {
    // Zero shear: dbeta/dn non-zero but shear components zero
    const Vec dbeta_active = Vec::fill(-ds_mix / (dn_mix * dn_mix),
                                       delta_mix(1) / (ds_mix * dn_mix),
                                       delta_mix(2) / (ds_mix * dn_mix));
    const Vec dbeta_zero_shear = Vec::fill(-ds_mix / (dn_mix * dn_mix),
                                           Scalar::zeros_like(dn_mix),
                                           Scalar::zeros_like(dn_mix));
    const Vec dbeta_nonzero_shear =
        where(ds_mix > 0, dbeta_active, dbeta_zero_shear);
    dbeta_ddelta = where(dn_mix > 0, dbeta_nonzero_shear, Vec::zeros_like(g));
  }

  // ----- Step 2: Damage initiation displacement jump delta_init -----
  const Scalar delta_n0 = _N / _K;
  const Scalar delta_s0 = _S / _K;

  const Scalar delta_mixed_denom =
      sqrt(delta_s0 * delta_s0 + (beta * delta_n0) * (beta * delta_n0));
  const Scalar delta_init_mixed =
      delta_n0 * delta_s0 * sqrt(1.0 + beta * beta) / delta_mixed_denom;

  Scalar delta_init = where(dn_mix > 0, delta_init_mixed, delta_s0);

  Vec ddelta_init_ddelta = Vec::zeros_like(g);
  if (!_lag_mode_mixity)
  {
    const Scalar ddelta_init_dbeta =
        delta_init_mixed * beta *
        (1.0 / (1.0 + beta * beta) -
         (delta_n0 * delta_n0) / (delta_mixed_denom * delta_mixed_denom));
    ddelta_init_ddelta = where(dn_mix > 0, dbeta_ddelta * ddelta_init_dbeta, Vec::zeros_like(g));
  }

  // ----- Step 3: Full-degradation displacement jump delta_final -----
  const Scalar delta_final_shear = std::sqrt(2.0) * 2.0 * _GIIc / _S;
  const Scalar beta2 = beta * beta;
  const Scalar r = beta2 / (1.0 + beta2);

  Scalar delta_final_mixed = delta_final_shear;
  Vec ddelta_final_ddelta = Vec::zeros_like(g);

  if (_criterion == "BK")
  {
    delta_final_mixed = 2.0 / (_K * delta_init) * (_GIc + (_GIIc - _GIc) * pow(r, _eta));

    if (!_lag_mode_mixity)
    {
      const Scalar ddelta_final_ddelta_init = -delta_final_mixed / delta_init;
      const Scalar dr_dbeta = 2.0 * beta / ((1.0 + beta2) * (1.0 + beta2));
      const Scalar ddelta_final_dbeta = 2.0 / (_K * delta_init) * (_GIIc - _GIc) * _eta *
                                        pow(r, _eta - 1.0) * dr_dbeta;
      const Vec dfinal_ddelta = ddelta_final_ddelta_init * ddelta_init_ddelta +
                                ddelta_final_dbeta * dbeta_ddelta;
      ddelta_final_ddelta = where(dn_mix > 0, dfinal_ddelta, Vec::zeros_like(g));
    }
  }
  else // POWER_LAW
  {
    const Scalar Gc_mixed = pow(1.0 / _GIc, _eta) + pow(beta2 / _GIIc, _eta);
    const Scalar factor = (2.0 + 2.0 * beta2) / (_K * delta_init);
    delta_final_mixed = factor * pow(Gc_mixed, -1.0 / _eta);

    if (!_lag_mode_mixity)
    {
      const Scalar ddelta_final_ddelta_init = -delta_final_mixed / delta_init;
      const Scalar dGc_dbeta = _eta * pow(beta2 / _GIIc, _eta - 1.0) * (2.0 * beta / _GIIc);
      const Scalar dfactor_dbeta = 4.0 * beta / (_K * delta_init);
      const Scalar ddelta_final_dbeta =
          dfactor_dbeta * pow(Gc_mixed, -1.0 / _eta) +
          factor * (-1.0 / _eta) * pow(Gc_mixed, -1.0 / _eta - 1.0) * dGc_dbeta;
      const Vec dfinal_ddelta = ddelta_final_ddelta_init * ddelta_init_ddelta +
                                ddelta_final_dbeta * dbeta_ddelta;
      ddelta_final_ddelta = where(dn_mix > 0, dfinal_ddelta, Vec::zeros_like(g));
    }
  }

  const Scalar delta_final = where(dn_mix > 0, delta_final_mixed, delta_final_shear);

  // ----- Step 4: Effective displacement jump delta_m -----
  const Scalar dn_eff = g_eff(0);
  const Scalar H = reg_heaviside(dn_eff);
  const Scalar dn_pos = H * dn_eff;
  const Scalar ds_eff_sq = g_eff(1) * g_eff(1) + g_eff(2) * g_eff(2);
  const Scalar delta_m = sqrt(dn_pos * dn_pos + ds_eff_sq);

  Vec ddelta_m_ddelta = Vec::zeros_like(g);
  if (!_lag_disp_jump)
  {
    const Scalar dH = reg_heaviside_deriv(dn_eff);
    const Scalar ddn_pos_ddn = H + dH * dn_eff;
    const Vec ddelta_m_num = Vec::fill(dn_pos * ddn_pos_ddn, g_eff(1), g_eff(2));
    ddelta_m_ddelta = where(delta_m > 1e-14, ddelta_m_num / delta_m, Vec::zeros_like(g));
  }

  // ----- Step 5: Damage -----
  const Scalar denom = delta_final - delta_init;
  const Scalar d_mid = delta_final * (delta_m - delta_init) / (delta_m * denom);
  const Scalar d_trial = where(
      delta_m < delta_init, Scalar::zeros_like(delta_m),
      where(delta_m > delta_final, Scalar::ones_like(delta_m), d_mid));

  const Scalar d_old = _damage_old();
  // Boolean masks for branch selection (avoid Double tensors as where conditions)
  const Scalar irrev_mask = d_trial < d_old;
  const Scalar d_irrev = where(irrev_mask, d_old, d_trial);

  // Derivative of d in softening branch
  const Scalar denom_m = delta_m * denom;
  const Vec numer_1 =
      ddelta_final_ddelta * (delta_m - delta_init) +
      delta_final * (ddelta_m_ddelta - ddelta_init_ddelta);
  const Vec numer_2 =
      ddelta_m_ddelta * denom + delta_m * (ddelta_final_ddelta - ddelta_init_ddelta);
  const Vec dd_mid_ddelta =
      numer_1 / denom_m - delta_final * (delta_m - delta_init) * numer_2 / (denom_m * denom_m);

  // Zero derivative outside the softening interval or when irreversible
  Vec dd_ddelta = where(
      delta_m < delta_init || delta_m > delta_final || irrev_mask,
      Vec::zeros_like(g),
      dd_mid_ddelta);

  // Viscous regularization
  const Scalar visc_over_dt = _visc / _dt();
  const Scalar visc_denom = visc_over_dt + 1.0;
  const Scalar damage = (d_irrev + visc_over_dt * d_old) / visc_denom;
  dd_ddelta = dd_ddelta / visc_denom;

  // ----- Step 6: Traction -----
  // Normal split: active for opening (dn > 0), inactive for compression
  const Scalar dn = g(0);
  const Scalar zero = Scalar::zeros_like(dn);
  const Scalar dn_pos_t = where(dn > 0, dn, zero);
  const Scalar dn_neg_t = where(dn < 0, dn, zero);
  const Vec g_active = Vec::fill(dn_pos_t, g(1), g(2));
  const Vec g_inactive = Vec::fill(dn_neg_t, zero, zero);

  if (out)
  {
    _damage = damage;
    _traction = (1.0 - damage) * _K * g_active + _K * g_inactive;
  }

  // ----- Step 7: Derivatives -----
  if (dout_din)
  {
    if (_jump.is_dependent())
    {
      _damage.d(_jump) = dd_ddelta;

      const Scalar one = Scalar::ones_like(dn);
      const Scalar dact_n_ddn = where(dn > 0, one, zero);
      const Scalar dinact_n_ddn = where(dn < 0, one, zero);
      const R2 dg_active_dg = R2::fill(dact_n_ddn, one, one);
      const R2 dg_inactive_dg = R2::fill(dinact_n_ddn, zero, zero);

      R2 dt_dg = (1.0 - damage) * _K * dg_active_dg + _K * dg_inactive_dg;
      dt_dg = dt_dg + (-_K) * outer(g_active, dd_ddelta);

      _traction.d(_jump) = dt_dg;
    }

    if (_damage_old.is_dependent())
    {
      // d(damage)/d(d_old):
      //   d_irrev = where(irrev_mask, d_old, d_trial)
      //   d(d_irrev)/d(d_old) = 1 when frozen, 0 when advancing
      //   damage = (d_irrev + visc_over_dt * d_old) / visc_denom
      //   d(damage)/d(d_old) = (d(d_irrev)/d(d_old) + visc_over_dt) / visc_denom
      //     frozen:   (1 + visc_over_dt) / visc_denom = 1
      //     advancing: visc_over_dt / visc_denom
      const Scalar dd_dd_old =
          where(irrev_mask, Scalar::ones_like(damage), visc_over_dt / visc_denom);
      _damage.d(_damage_old) = dd_dd_old;
      // d(traction)/d(d_old) = -K * g_active * d(damage)/d(d_old)
      _traction.d(_damage_old) = (-_K) * dd_dd_old * g_active;
    }
  }
}
} // namespace neml2
