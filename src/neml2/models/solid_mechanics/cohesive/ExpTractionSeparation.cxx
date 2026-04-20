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

#include "neml2/models/solid_mechanics/cohesive/ExpTractionSeparation.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/imap.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/where.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() = "Scalar-damage exponential traction-separation law: "
                  "d = 1 - exp(-delta_eff / delta0), T = (1-d) * (Gc/delta0^2) * delta";

  options.set_parameter<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Critical fracture energy G_c";

  options.set_parameter<TensorName<Scalar>>("softening_length");
  options.set("softening_length").doc() = "Softening length scale delta_0";

  options.set_parameter<TensorName<Scalar>>("tangential_weight");
  options.set("tangential_weight").doc() = "Tangential weighting factor beta in the effective jump";

  options.set<bool>("irreversible_damage") = true;
  options.set("irreversible_damage").doc() =
      "If true, enforce damage irreversibility via historical max effective jump";

  options.set_output("effective_displacement_jump_max") =
      VariableName(STATE, "internal", "delta_eff_max");
  options.set("effective_displacement_jump_max").doc() =
      "Historical maximum effective scalar displacement jump (state variable)";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparationLaw(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length")),
    _beta(declare_parameter<Scalar>("beta", "tangential_weight")),
    _irreversible(options.get<bool>("irreversible_damage")),
    _delta_eff_max(declare_output_variable<Scalar>("effective_displacement_jump_max")),
    _delta_eff_max_old(declare_input_variable<Scalar>(_delta_eff_max.name().old()))
{
}

void
ExpTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto & delta = _jump();
  const auto delta_n = delta(0);
  const auto delta_s1 = delta(1);
  const auto delta_s2 = delta(2);

  // delta_eff = sqrt(delta_n^2 + beta*(delta_s1^2 + delta_s2^2))
  const auto delta_eff_sq =
      delta_n * delta_n + _beta * (delta_s1 * delta_s1 + delta_s2 * delta_s2);
  const auto delta_eff = neml2::sqrt(delta_eff_sq);

  // Apply irreversibility by selecting between advancing and frozen state
  auto delta_eff_active = delta_eff;
  if (_irreversible)
  {
    const auto & old_max = _delta_eff_max_old();
    delta_eff_active = where(delta_eff > old_max, delta_eff, old_max);
  }
  _delta_eff_max = delta_eff_active;

  const auto exp_neg = neml2::exp(-delta_eff_active / _delta0);
  const auto c = _Gc / (_delta0 * _delta0);

  if (out)
    _traction = exp_neg * c * delta;

  if (dout_din)
    if (_jump.is_dependent())
    {
      // J = c * exp_neg * I - (c/delta0) * exp_neg * outer(delta, d_delta_eff_active/d_delta)
      // d(delta_eff)/d(delta) = [delta_n, beta*delta_s1, beta*delta_s2] / delta_eff
      //
      // Guard against delta_eff = 0 (no jump: gradient is zero)
      const auto safe_eff = where(delta_eff > Scalar::zeros_like(delta_eff), delta_eff,
                                  Scalar::ones_like(delta_eff));
      const auto ddelta_eff_ddelta =
          R2::fill(Scalar::ones_like(delta_eff), _beta, _beta) * delta / safe_eff;

      // When irreversible and frozen, the derivative of delta_eff_active w.r.t. delta is zero
      const auto advancing_mask =
          _irreversible
              ? where(delta_eff > _delta_eff_max_old(), Scalar::ones_like(delta_eff),
                      Scalar::zeros_like(delta_eff))
              : Scalar::ones_like(delta_eff);

      const auto I = imap_v<Vec>(_jump.options());
      _traction.d(_jump) =
          c * exp_neg * I -
          (c / _delta0) * exp_neg * outer(delta, advancing_mask * ddelta_eff_ddelta);
      _delta_eff_max.d(_jump) = advancing_mask * ddelta_eff_ddelta;
    }
}
} // namespace neml2
