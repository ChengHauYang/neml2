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

#include "neml2/models/solid_mechanics/ExpTractionSeparation.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/norm_sq.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/stack.h"
#include "neml2/tensors/functions/where.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() = "Exponential traction-separation law with optional irreversible effective jump "
                  "history.";

  options.set<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Fracture energy";

  options.set<TensorName<Scalar>>("softening_length");
  options.set("softening_length").doc() = "Softening length scale";

  options.set<TensorName<Scalar>>("beta");
  options.set("beta").doc() = "Tangential weighting factor";

  options.set_input("old_effective_jump_max") =
      VariableName(OLD_STATE, "internal", "effective_jump_max");
  options.set("old_effective_jump_max").doc() = "Previous maximum effective displacement jump";

  options.set_output("damage") = VariableName(STATE, "internal", "damage");
  options.set("damage").doc() = "Scalar cohesive damage";

  options.set_output("effective_jump_max") = VariableName(STATE, "internal", "effective_jump_max");
  options.set("effective_jump_max").doc() = "Updated maximum effective displacement jump";

  options.set_output("effective_displacement_jump");
  options.set("effective_displacement_jump").doc() = "Optional effective displacement jump vector";

  options.set<bool>("irreversible_damage") = true;
  options.set("irreversible_damage").doc() = "Freeze damage at the historical maximum effective "
                                             "displacement jump";

  options.set<double>("eps") = 1e-16;
  options.set("eps").doc() = "Small regularizer used in the effective jump norm";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparationLaw(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length")),
    _beta(declare_parameter<Scalar>("beta_weight", "beta")),
    _delta_eff_max_old(declare_input_variable<Scalar>("old_effective_jump_max")),
    _damage(declare_output_variable<Scalar>("damage")),
    _delta_eff_max(declare_output_variable<Scalar>("effective_jump_max")),
    _effective_jump(options.get("effective_displacement_jump").user_specified()
                        ? &declare_output_variable<Vec>("effective_displacement_jump")
                        : nullptr),
    _irreversible_damage(options.get<bool>("irreversible_damage")),
    _eps(options.get<double>("eps"))
{
}

void
ExpTractionSeparation::request_AD()
{
  TractionSeparationLaw::request_AD();
  _damage.request_AD(_jump);
  _damage.request_AD(_delta_eff_max_old);
  _delta_eff_max.request_AD(_jump);
  _delta_eff_max.request_AD(_delta_eff_max_old);
  _traction.request_AD(_delta_eff_max_old);
  if (_effective_jump)
    _effective_jump->request_AD(_jump);
}

void
ExpTractionSeparation::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  const auto jump = _jump();
  const Scalar delta_n = jump(0);
  const Scalar delta_t_sq = jump(1) * jump(1) + jump(2) * jump(2);
  const Scalar delta_eff = sqrt(delta_n * delta_n + _beta * delta_t_sq + _eps);
  const Scalar old_max = _delta_eff_max_old();
  const Scalar delta_eff_hist =
      _irreversible_damage ? where(delta_eff > old_max, delta_eff, old_max) : delta_eff;
  const Scalar damage = 1.0 - exp(-delta_eff_hist / _delta0);
  const Scalar prefactor = (1.0 - damage) * _Gc / (_delta0 * _delta0);

  _damage = damage;
  _delta_eff_max = delta_eff_hist;
  _traction = prefactor * jump;

  if (_effective_jump)
    *_effective_jump = Vec(base_stack({delta_n, sqrt(_beta) * jump(1), sqrt(_beta) * jump(2)}));
}
} // namespace neml2
