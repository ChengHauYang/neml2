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

#include "neml2/models/solid_mechanics/traction_separation/ExpTractionSeparation.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/where.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Exponential cohesive traction-separation law.";

  options.set_output("effective_displacement_jump") =
      VariableName(STATE, "interface", "effective_displacement_jump");
  options.set("effective_displacement_jump").doc() = "Weighted effective displacement jump vector";

  options.set_output("effective_displacement_jump_max") =
      VariableName(STATE, "internal", "effective_displacement_jump_max");
  options.set("effective_displacement_jump_max").doc() =
      "Historical maximum scalar effective displacement jump";

  options.set_output("damage") = VariableName(STATE, "internal", "interface_damage");
  options.set("damage").doc() = "Scalar exponential interface damage";

  options.set_parameter<TensorName<Scalar>>("critical_energy_release_rate");
  options.set("critical_energy_release_rate").doc() = "Critical energy release rate";

  options.set_parameter<TensorName<Scalar>>("softening_length");
  options.set("softening_length").doc() = "Exponential softening length scale";

  options.set_parameter<TensorName<Scalar>>("tangential_weight") = "1.0";
  options.set("tangential_weight").doc() = "Tangential weighting factor in the effective jump";

  options.set<double>("regularization") = 1e-16;
  options.set("regularization").doc() = "Small positive regularization inside the effective norm";

  options.set<bool>("irreversible_damage") = true;
  options.set("irreversible_damage").doc() = "Whether damage is driven by the historical maximum";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _effective_jump_max_old(declare_input_variable<Scalar>(
        options.get<VariableName>("effective_displacement_jump_max").old())),
    _effective_jump_max(declare_output_variable<Scalar>("effective_displacement_jump_max")),
    _damage(declare_output_variable<Scalar>("damage")),
    _effective_jump(declare_output_variable<Vec>("effective_displacement_jump")),
    _G_c(declare_parameter<Scalar>("G_c", "critical_energy_release_rate")),
    _delta_0(declare_parameter<Scalar>("delta_0", "softening_length")),
    _beta(declare_parameter<Scalar>("beta", "tangential_weight")),
    _eps(options.get<double>("regularization")),
    _irreversible_damage(options.get<bool>("irreversible_damage"))
{
}

void
ExpTractionSeparation::request_AD()
{
  _traction.request_AD(_jump);
  _damage.request_AD(_jump);
  _effective_jump.request_AD(_jump);
  _effective_jump_max.request_AD(_jump);
}

void
ExpTractionSeparation::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto delta_n = _jump()(0);
  const auto delta_s1 = _jump()(1);
  const auto delta_s2 = _jump()(2);
  const auto beta_sqrt = sqrt(_beta);
  auto delta_eff =
      sqrt(delta_n * delta_n + _beta * (delta_s1 * delta_s1 + delta_s2 * delta_s2) + _eps);

  _effective_jump = Vec::fill(delta_n, beta_sqrt * delta_s1, beta_sqrt * delta_s2);

  if (_irreversible_damage)
    delta_eff = where(delta_eff > _effective_jump_max_old, delta_eff, _effective_jump_max_old());

  _effective_jump_max = delta_eff;
  _damage = 1.0 - exp(-delta_eff / _delta_0);
  _traction = (1.0 - _damage) * _G_c / (_delta_0 * _delta_0) * _jump();
}
} // namespace neml2
