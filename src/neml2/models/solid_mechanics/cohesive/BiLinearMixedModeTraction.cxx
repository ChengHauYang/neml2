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

#include "neml2/models/solid_mechanics/cohesive/BiLinearMixedModeTraction.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/heaviside.h"
#include "neml2/tensors/functions/macaulay.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/imap.h"
#include "neml2/base/EnumSelection.h"

#include <cmath>

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() =
      "Bilinear mixed-mode traction-separation law with irreversible damage and optional "
      "viscous regularization. Supports BK (Benzeggagh-Kenane) and POWER_LAW mixed-mode "
      "propagation criteria.";

  options.set_parameter<TensorName<Scalar>>("penalty_stiffness");
  options.set("penalty_stiffness").doc() = "Penalty elastic stiffness K";

  options.set_parameter<TensorName<Scalar>>("normal_fracture_energy");
  options.set("normal_fracture_energy").doc() = "Mode I critical energy release rate G_Ic";

  options.set_parameter<TensorName<Scalar>>("shear_fracture_energy");
  options.set("shear_fracture_energy").doc() = "Mode II critical energy release rate G_IIc";

  options.set_parameter<TensorName<Scalar>>("tensile_strength");
  options.set("tensile_strength").doc() = "Tensile (normal) strength N";

  options.set_parameter<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Shear strength S";

  options.set_parameter<TensorName<Scalar>>("power_law_exponent");
  options.set("power_law_exponent").doc() = "Power law exponent eta for the BK or POWER_LAW criterion";

  options.set_parameter<TensorName<Scalar>>("viscosity");
  options.set("viscosity").doc() = "Viscous regularization coefficient (0 = no regularization)";

  EnumSelection crit_sel({"BK", "POWER_LAW"},
                         {static_cast<int>(Criterion::BK),
                          static_cast<int>(Criterion::POWER_LAW)},
                         "BK");
  options.set<EnumSelection>("criterion") = crit_sel;
  options.set("criterion").doc() =
      "Mixed-mode propagation criterion. Options: " + crit_sel.join();

  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() =
      "Use displacement jump from the previous step when computing mode-mixity ratio beta "
      "and the critical displacement jumps delta_init and delta_final; their derivatives "
      "w.r.t. the current jump are then zero (consistent tangent for lagged quantities)";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() =
      "Use displacement jump from the previous step when computing the effective mixed-mode "
      "displacement jump delta_m; derivative w.r.t. current jump is then zero";

  options.set_input("time_step") = VariableName(FORCES, "time_step");
  options.set("time_step").doc() = "Time step size dt (required when viscosity > 0)";

  options.set_output("damage") = VariableName(STATE, "damage");
  options.set("damage").doc() = "Scalar damage variable d in [0, 1]";

  options.set_output("stored_displacement_jump") =
      VariableName(STATE, "internal", "displacement_jump_store");
  options.set("stored_displacement_jump").doc() =
      "Copy of the current displacement jump stored as state for use in the next step";

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _K(declare_parameter<Scalar>("K", "penalty_stiffness")),
    _GI_c(declare_parameter<Scalar>("GIc", "normal_fracture_energy")),
    _GII_c(declare_parameter<Scalar>("GIIc", "shear_fracture_energy")),
    _N(declare_parameter<Scalar>("N", "tensile_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "power_law_exponent")),
    _viscosity(declare_parameter<Scalar>("viscosity", "viscosity")),
    _criterion(options.get<EnumSelection>("criterion").as<Criterion>()),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_displacement_jump(options.get<bool>("lag_displacement_jump")),
    _dt(declare_input_variable<Scalar>("time_step")),
    _damage(declare_output_variable<Scalar>("damage")),
    _damage_old(declare_input_variable<Scalar>(_damage.name().old())),
    _jump_stored(declare_output_variable<Vec>("stored_displacement_jump")),
    _jump_old(declare_input_variable<Vec>(_jump_stored.name().old()))
{
}

void
BiLinearMixedModeTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto & delta = _jump();
  const auto delta_n = delta(0);
  const auto one = Scalar::ones_like(delta_n);
  const auto zero = Scalar::zeros_like(delta_n);

  // ---- Step 1: Mode mixity ----
  // Select jump for mixity computation (lagged or current)
  const auto & delta_sel = _lag_mode_mixity ? _jump_old() : delta;
  const auto sel_n = delta_sel(0);
  const auto sel_s1 = delta_sel(1);
  const auto sel_s2 = delta_sel(2);

  const auto sel_s = neml2::sqrt(sel_s1 * sel_s1 + sel_s2 * sel_s2);
  const auto normal_open = sel_n > zero;
  const auto safe_sel_n = where(normal_open, sel_n, one);
  const auto beta = where(normal_open, sel_s / safe_sel_n, zero);
  const auto beta2 = beta * beta;

  // ---- Step 2: Damage-initiation displacement jump delta_init ----
  const auto delta_n0 = _N / _K;
  const auto delta_s0 = _S / _K;
  const auto delta_mixed_sq = delta_s0 * delta_s0 + beta2 * delta_n0 * delta_n0;
  const auto safe_mixed = neml2::sqrt(
      where(delta_mixed_sq > zero, delta_mixed_sq, one));
  const auto delta_init_mixed = delta_n0 * delta_s0 * neml2::sqrt(one + beta2) / safe_mixed;
  const auto delta_init = where(normal_open, delta_init_mixed, delta_s0);

  // ---- Step 3: Full-degradation displacement jump delta_final ----
  const auto delta_final_shear = std::sqrt(2.0) * 2.0 * _GII_c / _S;

  Scalar delta_final;
  if (_criterion == Criterion::BK)
  {
    // BK criterion: delta_final = 2/(K*delta_init) * (GIc + (GIIc-GIc) * bsr^eta)
    const auto bsr = beta2 / (one + beta2);
    const auto bsr_eta = neml2::pow(bsr, _eta);
    const auto delta_final_bk =
        2.0 * (_GI_c + (_GII_c - _GI_c) * bsr_eta) / (_K * delta_init);
    delta_final = where(normal_open, delta_final_bk, delta_final_shear);
  }
  else
  {
    // POWER_LAW criterion
    const auto GI_c_inv_eta = neml2::pow(one / _GI_c, _eta);
    const auto beta2_GII_eta = neml2::pow(beta2 / _GII_c, _eta);
    const auto Gc_mixed = GI_c_inv_eta + beta2_GII_eta;
    const auto Gc_mixed_inv_eta = neml2::pow(one / Gc_mixed, one / _eta);
    const auto delta_final_pl = 2.0 * (one + beta2) / (_K * delta_init) * Gc_mixed_inv_eta;
    delta_final = where(normal_open, delta_final_pl, delta_final_shear);
  }

  // ---- Step 4: Effective mixed-mode displacement jump delta_m ----
  const auto & delta_eff = _lag_displacement_jump ? _jump_old() : delta;
  const auto delta_n_eff = delta_eff(0);
  const auto delta_s1_eff = delta_eff(1);
  const auto delta_s2_eff = delta_eff(2);

  const auto H_eff = neml2::heaviside(delta_n_eff);
  const auto delta_n_pos_eff = neml2::macaulay(delta_n_eff);
  const auto delta_m = neml2::sqrt(delta_s1_eff * delta_s1_eff + delta_s2_eff * delta_s2_eff +
                                   delta_n_pos_eff * delta_n_pos_eff);

  // ---- Step 5: Damage (bilinear law + irreversibility + viscous regularization) ----
  const auto safe_dm = where(delta_m > zero, delta_m, one);
  const auto safe_df_di =
      where((delta_final - delta_init) > zero, delta_final - delta_init, one);
  const auto d_soft = delta_final * (delta_m - delta_init) / (safe_dm * safe_df_di);
  const auto d_trial = where(delta_m < delta_init, zero,
                              where(delta_m > delta_final, one, d_soft));
  const auto advancing = where(d_trial > _damage_old(), one, zero);
  const auto in_soft =
      one - where(delta_m < delta_init, one, zero) - where(delta_m > delta_final, one, zero);
  const auto d_irrev = where(d_trial > _damage_old(), d_trial, _damage_old());
  const auto visc_denom = _viscosity / _dt() + one;
  const auto d = (d_irrev + _viscosity * _damage_old() / _dt()) / visc_denom;

  // ---- Step 6: Traction (from current jump, always unlagged for traction itself) ----
  const auto H_n = neml2::heaviside(delta_n);
  const auto delta_active_t =
      R2::fill(H_n, Scalar::ones_like(H_n), Scalar::ones_like(H_n)) * delta;
  const auto delta_inactive_t =
      R2::fill(Scalar::ones_like(H_n) - H_n, Scalar::zeros_like(H_n), Scalar::zeros_like(H_n)) *
      delta;

  if (out)
  {
    _damage = d;
    _jump_stored = delta;
    _traction = (one - d) * _K * delta_active_t + _K * delta_inactive_t;
  }

  // ---- Step 7: Jacobians ----
  if (dout_din)
  {
    // d(jump_stored)/d(jump) = I (stored jump equals current jump)
    if (_jump.is_dependent())
      _jump_stored.d(_jump) = imap_v<Vec>(_jump.options());

    // d(damage)/d(old_damage): frozen → (1-advancing), advancing → 0; plus viscous term
    const auto dd_d_dold = (one - advancing + _viscosity / _dt()) / visc_denom;
    _damage.d(_damage_old) = dd_d_dold;

    // d(traction)/d(old_damage)
    _traction.d(_damage_old) = -_K * delta_active_t * dd_d_dold;

    if (_jump.is_dependent())
    {
      const auto dd_dm =
          in_soft * delta_final * delta_init / (safe_df_di * safe_dm * safe_dm);

      // ddelta_m/d(delta): gradient of delta_m w.r.t. the jump used for delta_m
      // = [delta_n_pos_eff, delta_s1_eff, delta_s2_eff] / delta_m
      //   (uses d(macaulay(x))/dx = heaviside(x) for the normal component)
      const auto delta_active_eff =
          R2::fill(H_eff, Scalar::ones_like(H_eff), Scalar::ones_like(H_eff)) * delta_eff;
      const auto ddelta_m_ddelta = delta_active_eff / safe_dm;

      const auto lag_m = _lag_displacement_jump ? zero : one;
      const auto dd_ddelta = advancing * dd_dm / visc_denom * lag_m * ddelta_m_ddelta;

      _traction.d(_jump) =
          (one - d) * _K * R2::fill(H_n, one, one) +
          _K * R2::fill(one - H_n, zero, zero) -
          _K * outer(delta_active_t, dd_ddelta);
    }
  }
}
} // namespace neml2
