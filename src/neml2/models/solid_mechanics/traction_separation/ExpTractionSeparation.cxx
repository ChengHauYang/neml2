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

#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/imap.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() = "Exponential cohesive traction-separation law with optional irreversible damage "
                  "based on an effective displacement jump.";

  options.set_parameter<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Fracture energy \\f$ G_c \\f$";

  options.set_parameter<TensorName<Scalar>>("softening_length");
  options.set("softening_length").doc() = "Characteristic displacement scale \\f$ \\delta_0 \\f$";

  options.set_parameter<TensorName<Scalar>>("shear_weight") = "1.0";
  options.set("shear_weight").doc() =
      "Tangential weighting coefficient \\f$ \\beta \\f$ in the effective jump norm";

  options.set<double>("eps") = 1e-16;
  options.set("eps").doc() = "Small regularization added to the effective jump norm";

  options.set<bool>("irreversible_damage") = true;
  options.set("irreversible_damage").doc() =
      "Use the maximum historical effective jump to prevent damage healing";

  options.set_output("damage") = VariableName(STATE, "internal", "damage");
  options.set("damage").doc() = "Cohesive damage variable";

  options.set_output("effective_jump_max") =
      VariableName(STATE, "internal", "effective_displacement_jump_scalar_max");
  options.set("effective_jump_max").doc() = "Maximum effective scalar displacement jump";

  options.set<VariableName>("effective_jump") = VariableName(STATE, "internal", "effective_jump");
  options.set("effective_jump").doc() =
      "Optional weighted displacement jump vector for post-processing";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparationLaw(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length")),
    _beta(declare_parameter<Scalar>("beta", "shear_weight")),
    _eps(options.get<double>("eps")),
    _irreversible_damage(options.get<bool>("irreversible_damage")),
    _damage(declare_output_variable<Scalar>("damage")),
    _effective_jump(options.get<VariableName>("effective_jump").empty()
                        ? nullptr
                        : &declare_output_variable<Vec>("effective_jump")),
    _effective_jump_max(declare_output_variable<Scalar>("effective_jump_max")),
    _effective_jump_max_old(_irreversible_damage
                                ? &declare_input_variable<Scalar>(_effective_jump_max.name().old())
                                : nullptr)
{
}

void
ExpTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto jump = _jump();
  const auto delta_n = jump(0);
  const auto delta_t1 = jump(1);
  const auto delta_t2 = jump(2);
  const auto delta_eff_raw =
      sqrt(delta_n * delta_n + _beta * (delta_t1 * delta_t1 + delta_t2 * delta_t2) + _eps);
  const auto ddelta_eff_raw =
      Vec::fill(delta_n, _beta * delta_t1, _beta * delta_t2) / delta_eff_raw;

  Scalar delta_eff = delta_eff_raw;
  Vec ddelta_eff = ddelta_eff_raw;
  if (_irreversible_damage)
  {
    const auto advancing = delta_eff_raw > (*_effective_jump_max_old)();
    delta_eff = where(advancing, delta_eff_raw, (*_effective_jump_max_old)());
    ddelta_eff = where(advancing, ddelta_eff_raw, Vec::zeros_like(ddelta_eff_raw));
  }

  const auto damage = 1.0 - exp(-delta_eff / _delta0);
  const auto stiffness = _Gc / (_delta0 * _delta0);
  const auto scale = exp(-delta_eff / _delta0) * stiffness;

  if (out)
  {
    _damage = damage;
    _effective_jump_max = delta_eff;
    if (_effective_jump)
      *_effective_jump = Vec::fill(delta_n, sqrt(_beta) * delta_t1, sqrt(_beta) * delta_t2);
    _traction = scale * jump;
  }

  if (dout_din)
  {
    if (_jump.is_dependent())
    {
      _damage.d(_jump) = exp(-delta_eff / _delta0) / _delta0 * ddelta_eff;
      _effective_jump_max.d(_jump) = ddelta_eff;
      if (_effective_jump)
        _effective_jump->d(_jump) = R2::fill(Scalar::ones_like(delta_n), sqrt(_beta), sqrt(_beta));
      _traction.d(_jump) =
          scale * (imap_v<Vec>(_jump.options()) - neml2::outer(jump, ddelta_eff / _delta0));
    }
  }
}
} // namespace neml2
