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
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/tanh.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/inner.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/misc/types.h"

namespace neml2
{
register_NEML2_object(BiLinearMixedModeTraction);

OptionSet
BiLinearMixedModeTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Bilinear mixed-mode traction-separation law.";

  EnumSelection criterion_selection({"BK", "POWER_LAW"}, "BK");
  options.set<EnumSelection>("criterion") = criterion_selection;
  options.set("criterion").doc() = "Mixed mode propagation criterion";

  options.set_parameter<TensorName<Scalar>>("normal_stiffness");
  options.set("normal_stiffness").doc() = "Penalty elastic stiffness";

  options.set_parameter<TensorName<Scalar>>("GI_c");
  options.set("GI_c").doc() = "Mode I critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("GII_c");
  options.set("GII_c").doc() = "Mode II critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("tensile_strength");
  options.set("tensile_strength").doc() = "Tensile (normal) strength";

  options.set_parameter<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Shear strength";

  options.set_parameter<TensorName<Scalar>>("exponent") = TensorName<Scalar>("1.0");
  options.set("exponent").doc() = "Power law exponent (BK or POWER_LAW criterion)";

  options.set_parameter<TensorName<Scalar>>("viscosity") = TensorName<Scalar>("0.0");
  options.set("viscosity").doc() = "Viscous regularization coefficient (0 = no regularization)";

  options.set_parameter<TensorName<Scalar>>("alpha") = TensorName<Scalar>("1e-10");
  options.set("alpha").doc() = "Regularization parameter for the Macaulay bracket";

  options.set_input("time") = VariableName(FORCES, "t");
  options.set("time").doc() = "Time";

  options.set_output("damage") = VariableName(STATE, "internal", "d");
  options.set("damage").doc() = "Damage variable";

  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() = "Use old displacement jump when computing mode mixity";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() =
      "Use old displacement jump when computing effective displacement jump";

  return options;
}

BiLinearMixedModeTraction::BiLinearMixedModeTraction(const OptionSet & options)
  : TractionSeparation(options),
    _criterion(options.get<EnumSelection>("criterion").as<Criterion>()),
    _K(declare_parameter<Scalar>("K", "normal_stiffness")),
    _GI_c(declare_parameter<Scalar>("GI_c", "GI_c")),
    _GII_c(declare_parameter<Scalar>("GII_c", "GII_c")),
    _N(declare_parameter<Scalar>("N", "tensile_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "exponent")),
    _viscosity(declare_parameter<Scalar>("viscosity", "viscosity")),
    _alpha(declare_parameter<Scalar>("alpha", "alpha")),
    _d(declare_output_variable<Scalar>("damage")),
    _d_old(declare_input_variable<Scalar>(_d.name().old())),
    _displacement_jump_old(declare_input_variable<Vec>(_displacement_jump.name().old())),
    _t(declare_input_variable<Scalar>(VariableName(FORCES, "t"))),
    _t_old(declare_input_variable<Scalar>(VariableName(FORCES, "t").old())),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_disp_jump(options.get<bool>("lag_displacement_jump"))
{
}

void
BiLinearMixedModeTraction::request_AD()
{
  _traction.request_AD(_displacement_jump);
  _traction.request_AD(_t);
  _traction.request_AD(_t_old);
  _d.request_AD(_displacement_jump);
  _d.request_AD(_t);
  _d.request_AD(_t_old);
}

void
BiLinearMixedModeTraction::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  // 1. Compute mode mixity ratio
  const auto & delta_mix = _lag_mode_mixity ? _displacement_jump_old : _displacement_jump;
  auto delta_n_mix = delta_mix()(0);
  auto delta_s1_mix = delta_mix()(1);
  auto delta_s2_mix = delta_mix()(2);
  auto delta_s_mix = neml2::sqrt(delta_n_mix * delta_n_mix + delta_s1_mix * delta_s1_mix + delta_s2_mix * delta_s2_mix);

  auto opening = delta_n_mix > 0;
  auto beta = neml2::where(opening, (delta_s_mix / delta_n_mix), Scalar::zeros_like(delta_n_mix));

  // 2. Compute damage-initiation displacement jump
  auto delta_normal0 = _N / _K;
  auto delta_shear0 = _S / _K;

  auto options = _displacement_jump.options();

  auto delta_init = neml2::where(opening,
                                 (delta_normal0 * delta_shear0 * neml2::sqrt(Scalar::fill(1.0, options) + beta * beta) /
                                  neml2::sqrt(delta_shear0 * delta_shear0 +
                                              beta * beta * delta_normal0 * delta_normal0)),
                                 delta_shear0);

  // 3. Compute full-degradation displacement jump
  auto delta_final = Scalar::zeros_like(delta_init);
  if (_criterion == Criterion::BK)
  {
    auto beta_sq_ratio = beta * beta / (Scalar::fill(1.0, options) + beta * beta);
    delta_final = neml2::where(opening,
                               (Scalar::fill(2.0, options) / (_K * delta_init) *
                                (_GI_c + (_GII_c - _GI_c) * neml2::pow(beta_sq_ratio, _eta))),
                               (neml2::sqrt(Scalar::fill(2.0, options)) * 2.0 * _GII_c / _S));
  }
  else // POWER_LAW
  {
    auto Gc_mixed = neml2::pow(Scalar::fill(1.0, options) / _GI_c, _eta) + neml2::pow(beta * beta / _GII_c, _eta);
    delta_final = neml2::where(opening,
                               ((Scalar::fill(2.0, options) + 2.0 * beta * beta) / (_K * delta_init) *
                                neml2::pow(Gc_mixed, Scalar::fill(-1.0, options) / _eta)),
                               (neml2::sqrt(Scalar::fill(2.0, options)) * 2.0 * _GII_c / _S));
  }

  // 4. Compute effective displacement jump
  const auto & delta = _lag_disp_jump ? _displacement_jump_old : _displacement_jump;
  auto delta_n = delta()(0);
  auto delta_s1 = delta()(1);
  auto delta_s2 = delta()(2);

  // Regularized Macaulay bracket
  auto H = Scalar::fill(0.5, options) * (Scalar::fill(1.0, options) + neml2::tanh(delta_n / _alpha));
  auto delta_n_pos = H * delta_n;
  auto delta_m = neml2::sqrt(delta_s1 * delta_s1 + delta_s2 * delta_s2 + delta_n_pos * delta_n_pos);

  // 5. Compute damage
  auto d_trial = Scalar::zeros_like(delta_m);
  auto mask1 = delta_m > delta_init;
  auto mask2 = delta_m > delta_final;
  d_trial = neml2::where(
      mask1, delta_final * (delta_m - delta_init) / (delta_m * (delta_final - delta_init)), Scalar::zeros_like(delta_m));
  d_trial = neml2::where(mask2, Scalar::fill(1.0, options), d_trial);

  // Irreversibility
  auto d_irr = neml2::where(d_trial < _d_old(), _d_old(), d_trial);

  // Viscous regularization
  auto dt = _t() - _t_old();
  auto d_new = neml2::where(dt > 0, (d_irr + _viscosity * _d_old / dt) / (_viscosity / dt + Scalar::fill(1.0, options)), d_irr);
  _d = d_new;

  // 6. Compute traction
  auto delta_n_neg = delta_n - delta_n_pos;
  Vec delta_active = Vec::fill(delta_n_pos, delta_s1, delta_s2);
  Vec delta_inactive = Vec::fill(delta_n_neg, Scalar::zeros_like(delta_n), Scalar::zeros_like(delta_n));

  _traction = (Scalar::fill(1.0, options) - d_new) * _K * delta_active + _K * delta_inactive;
}
} // namespace neml2
