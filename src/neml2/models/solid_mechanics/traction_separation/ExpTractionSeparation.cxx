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
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/imap.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() +=
      " following an isotropic exponential softening law d = 1 - exp(-delta_eff/delta0). "
      "Irreversibility is enforced by tracking the historical maximum effective jump.";

  options.set_input("effective_displacement_jump_max") =
      VariableName(STATE, "effective_displacement_jump_max");
  options.set("effective_displacement_jump_max").doc() =
      "Current maximum effective displacement jump (state variable for irreversibility)";

  options.set_output("damage") = VariableName(STATE, "damage");
  options.set("damage").doc() = "Scalar damage variable d in [0,1)";

  options.set_parameter<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Fracture energy G_c";

  options.set_parameter<TensorName<Scalar>>("softening_length");
  options.set("softening_length").doc() = "Softening length scale delta_0";

  options.set_parameter<TensorName<Scalar>>("tangential_weight");
  options.set("tangential_weight").doc() = "Tangential weighting factor beta in the effective jump";

  options.set<double>("regularization") = 1e-16;
  options.set("regularization").doc() = "Small epsilon added inside sqrt to avoid singularity";

  options.set<bool>("irreversible_damage") = true;
  options.set("irreversible_damage").doc() =
      "Enforce irreversibility by using the historical maximum effective jump";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparationLaw(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length")),
    _beta(declare_parameter<Scalar>("beta", "tangential_weight")),
    _eps(options.get<double>("regularization")),
    _irreversible(options.get<bool>("irreversible_damage")),
    _delta_eff_max(declare_output_variable<Scalar>("effective_displacement_jump_max")),
    _delta_eff_max_old(
        declare_input_variable<Scalar>(_delta_eff_max.name().old())),
    _damage(declare_output_variable<Scalar>("damage"))
{
}

void
ExpTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto & jump = _displacement_jump();
  const auto delta_n = jump(0);
  const auto delta_s1 = jump(1);
  const auto delta_s2 = jump(2);

  // Effective scalar displacement jump (with regularizer to avoid sqrt singularity)
  const auto delta_t_sq = delta_s1 * delta_s1 + delta_s2 * delta_s2;
  const auto delta_n_sq = delta_n * delta_n;
  auto delta_eff = sqrt(delta_n_sq + _beta * delta_t_sq + Scalar(_eps, delta_n.options()));

  // Derivative of delta_eff w.r.t. delta (before irreversibility) — a Vec
  // d(delta_eff)/d(delta) = (delta_n, beta*delta_s1, beta*delta_s2) / delta_eff
  const auto dd_eff_ddelta_before =
      Vec::fill(delta_n / delta_eff, _beta * delta_s1 / delta_eff, _beta * delta_s2 / delta_eff);

  // Irreversibility: clamp delta_eff to its historical maximum
  const auto & old_max = _delta_eff_max_old();
  Scalar advancing;
  if (_irreversible)
    advancing = delta_eff > old_max;
  else
    advancing = Scalar::ones_like(delta_eff);

  if (_irreversible)
    delta_eff = where(advancing, delta_eff, old_max);

  // dd_eff_ddelta is zero when frozen at old_max
  const auto zero_vec = Vec::zeros(_displacement_jump.options());
  const auto dd_eff_ddelta = where(advancing, dd_eff_ddelta_before, zero_vec);

  // Update state: max effective jump
  if (out)
    _delta_eff_max = where(advancing, delta_eff, old_max);

  // Damage: d = 1 - exp(-delta_eff / delta0)
  const auto exp_term = exp(-delta_eff / _delta0);
  const auto d = 1.0 - exp_term;

  // Traction prefactor: c = Gc / delta0^2
  const auto c = _Gc / (_delta0 * _delta0);

  if (out)
  {
    _damage = d;
    _traction = (1.0 - d) * c * jump;
  }

  if (dout_din)
    if (_displacement_jump.is_dependent())
    {
      // d(d)/d(delta) = exp_term / delta0 * dd_eff_ddelta
      const auto dd_ddelta = (exp_term / _delta0) * dd_eff_ddelta;

      // d(T)/d(delta) = (1-d)*c*I - c * outer(jump, dd_ddelta)
      const auto I = imap_v<Vec>(_displacement_jump.options());
      _traction.d(_displacement_jump) =
          (1.0 - d) * c * I + (-c) * outer(jump, dd_ddelta);

      // d(damage)/d(delta) = exp_term / delta0 * dd_eff_ddelta
      _damage.d(_displacement_jump) = dd_ddelta;

      // d(delta_eff_max)/d(delta) = dd_eff_ddelta (when advancing)
      _delta_eff_max.d(_displacement_jump) = where(advancing, dd_eff_ddelta_before, zero_vec);
    }
}
} // namespace neml2
