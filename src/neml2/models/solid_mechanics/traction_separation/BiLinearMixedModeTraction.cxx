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
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/heaviside.h"
#include "neml2/tensors/functions/macaulay.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/outer.h"

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Bilinear mixed-mode traction-separation law";

  options.set_parameter<TensorName<Scalar>>("penalty_stiffness") = "1e6";
  options.set("penalty_stiffness").doc() = "Penalty elastic stiffness K";

  options.set_parameter<TensorName<Scalar>>("GI_c") = "1.0";
  options.set("GI_c").doc() = "Mode I critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("GII_c") = "1.0";
  options.set("GII_c").doc() = "Mode II critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("normal_strength") = "1.0";
  options.set("normal_strength").doc() = "Tensile (normal) strength N";

  options.set_parameter<TensorName<Scalar>>("shear_strength") = "1.0";
  options.set("shear_strength").doc() = "Shear strength S";

  options.set_parameter<TensorName<Scalar>>("exponent") = "1.0";
  options.set("exponent").doc() = "Power law exponent eta";

  options.set_parameter<TensorName<Scalar>>("viscosity") = "0.0";
  options.set("viscosity").doc() = "Viscous regularization coefficient";

  options.set_parameter<TensorName<Scalar>>("alpha") = "1e-10";
  options.set("alpha").doc() = "Regularization parameter for the Macaulay bracket";

  options.set<std::string>("criterion") = "BK";
  options.set("criterion").doc() = "Mixed mode propagation criterion: BK or POWER_LAW";

  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() = "Use delta_old when computing mixity";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() = "Use delta_old when computing delta_m";

  options.set_input("old_damage") = VariableName(OLD_STATE, "internal", "d");
  options.set("old_damage").doc() = "Damage variable from previous step";

  options.set_output("damage") = VariableName(STATE, "internal", "d");
  options.set("damage").doc() = "Damage variable";

  options.set_input("time") = VariableName(FORCES, "t");
  options.set("time").doc() = "Time";

  options.set_input("old_time") = VariableName(OLD_FORCES, "t");
  options.set("old_time").doc() = "Time from previous step";

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparation(options),
    _K(declare_parameter<Scalar>("penalty_stiffness", "penalty_stiffness", true)),
    _GI_c(declare_parameter<Scalar>("GI_c", "GI_c", true)),
    _GII_c(declare_parameter<Scalar>("GII_c", "GII_c", true)),
    _N(declare_parameter<Scalar>("normal_strength", "normal_strength", true)),
    _S(declare_parameter<Scalar>("shear_strength", "shear_strength", true)),
    _eta(declare_parameter<Scalar>("exponent", "exponent", true)),
    _viscosity(declare_parameter<Scalar>("viscosity", "viscosity", true)),
    _alpha(declare_parameter<Scalar>("alpha", "alpha", true)),
    _criterion(options.get<std::string>("criterion") == "BK" ? 0 : 1),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_disp_jump(options.get<bool>("lag_displacement_jump")),
    _d_old(declare_input_variable<Scalar>("old_damage")),
    _t(declare_input_variable<Scalar>("time")),
    _t_old(declare_input_variable<Scalar>("old_time")),
    _jump_old(declare_input_variable<Vec>(_jump.name().old())),
    _d(declare_output_variable<Scalar>("damage"))
{
}

void
BiLinearMixedModeTraction::set_value(bool out, bool dout_din, bool d2out_din2)
{
  const auto t = _t();
  const auto t_old = _t_old();
  const auto diff_t = t - t_old;
  const auto dt = neml2::where(diff_t > 1e-16, diff_t, Scalar::full(1e-16, t.options()));

  const auto delta = _jump();
  const auto delta_old = _jump_old();

  // 1. Compute mode mixity ratio
  const auto delta_mix = _lag_mode_mixity ? delta_old : delta;
  const auto dn_mix = delta_mix(0);
  const auto ds_mix = neml2::sqrt(delta_mix(1) * delta_mix(1) + delta_mix(2) * delta_mix(2) + _alpha * _alpha);
  const auto beta = neml2::where(dn_mix > 0.0, ds_mix / dn_mix, Scalar::zeros(delta.options()));

  // 2. Compute damage-initiation displacement jump
  const auto delta_normal0 = _N / _K;
  const auto delta_shear0 = _S / _K;
  const auto delta_mixed =
      neml2::sqrt(delta_shear0 * delta_shear0 + (beta * delta_normal0) * (beta * delta_normal0));
  const auto delta_init_opening =
      delta_normal0 * delta_shear0 * neml2::sqrt(1.0 + beta * beta) / delta_mixed;
  const auto delta_init = neml2::where(dn_mix > 0.0, delta_init_opening, delta_shear0);

  // 3. Compute full-degradation displacement jump
  Scalar delta_final;
  if (_criterion == 0) // BK
  {
    const auto beta_sq_ratio = beta * beta / (1.0 + beta * beta);
    const auto delta_final_opening =
        2.0 / (_K * delta_init) * (_GI_c + (_GII_c - _GI_c) * neml2::pow(beta_sq_ratio, _eta));
    const auto delta_final_pure_shear_exact = 1.41421356 * 2.0 * _GII_c / _S;
    delta_final = neml2::where(dn_mix > 0.0, delta_final_opening, delta_final_pure_shear_exact);
  }
  else // POWER_LAW
  {
    const auto Gc_mixed = neml2::pow(1.0 / _GI_c, _eta) + neml2::pow(beta * beta / _GII_c, _eta);
    const auto delta_final_opening =
        (2.0 + 2.0 * beta * beta) / (_K * delta_init) * neml2::pow(Gc_mixed, -1.0 / _eta);
    const auto delta_final_pure_shear_exact = 1.41421356 * 2.0 * _GII_c / _S;
    delta_final = neml2::where(dn_mix > 0.0, delta_final_opening, delta_final_pure_shear_exact);
  }

  // 4. Compute effective displacement jump
  const auto delta_eff_jump = _lag_disp_jump ? delta_old : delta;
  const auto dn_eff = delta_eff_jump(0);
  const auto dn_pos = neml2::macaulay(dn_eff);
  const auto delta_m =
      neml2::sqrt(delta_eff_jump(1) * delta_eff_jump(1) + delta_eff_jump(2) * delta_eff_jump(2) +
                  dn_pos * dn_pos + _alpha * _alpha);

  // 5. Compute damage
  const auto d_trial = neml2::where(
      delta_m < delta_init,
      Scalar::zeros(delta.options()),
      neml2::where(delta_m > delta_final,
                   Scalar::ones(delta.options()),
                   delta_final * (delta_m - delta_init) / (delta_m * (delta_final - delta_init))));

  const auto d_irreversible = neml2::where(d_trial < _d_old(), _d_old(), d_trial);

  // Viscous regularization
  const auto d_viscous = (d_irreversible + _viscosity * _d_old() / dt) / (_viscosity / dt + 1.0);

  if (out)
  {
    _d = d_viscous;

    // 6. Compute traction
    const auto H = neml2::heaviside(delta(0));
    const auto dn_pos_cur = H * delta(0);
    const auto dn_neg_cur = delta(0) - dn_pos_cur;

    const auto delta_active = Vec::fill(dn_pos_cur, delta(1), delta(2));
    const auto delta_inactive = Vec::fill(dn_neg_cur, Scalar::zeros(delta.options()), Scalar::zeros(delta.options()));

    _traction = (1.0 - d_viscous) * _K * delta_active + _K * delta_inactive;
  }

  if (dout_din || d2out_din2)
  {
    // AD
  }
}

void
BiLinearMixedModeTraction::request_AD()
{
  Model::request_AD(_traction, _jump);
  Model::request_AD(_d, _jump);
}

} // namespace neml2
