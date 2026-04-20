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

#include "neml2/models/solid_mechanics/cohesive_zone/ExpTractionSeparation.h"
#include "neml2/tensors/tensors.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/max.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/outer.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Exponential traction separation law.";

  options.set_parameter<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Fracture energy Gc";

  options.set_parameter<TensorName<Scalar>>("softening_length");
  options.set("softening_length").doc() = "Softening length scale delta0";

  options.set_parameter<TensorName<Scalar>>("weighting_factor");
  options.set("weighting_factor").doc() = "Tangential weighting factor beta";
  options.set<TensorName<Scalar>>("weighting_factor") = "1.0";

  options.set<bool>("irreversible_damage") = true;
  options.set("irreversible_damage").doc() = "Whether to use irreversible damage";

  options.set_output("max_effective_displacement_jump") =
      VariableName(STATE, "internal", "delta_max");
  options.set("max_effective_displacement_jump").doc() = "Maximum effective displacement jump";

  options.set_output("damage") = VariableName(STATE, "internal", "d");
  options.set("damage").doc() = "Damage variable";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy", false)),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length", false)),
    _beta(declare_parameter<Scalar>("beta", "weighting_factor", false)),
    _irreversible_damage(options.get<bool>("irreversible_damage")),
    _delta_max(declare_output_variable<Scalar>("max_effective_displacement_jump")),
    _delta_max_old(_irreversible_damage
                       ? &declare_input_variable<Scalar>(_delta_max.name().old())
                       : nullptr),
    _damage(declare_output_variable<Scalar>("damage"))
{
}

void
ExpTractionSeparation::set_value(bool out, bool dout_din, [[maybe_unused]] bool d2out_din2)
{
  auto delta = _delta();
  auto delta_n = delta(0);
  auto delta_t1 = delta(1);
  auto delta_t2 = delta(2);

  auto delta_eff = sqrt(delta_n * delta_n + _beta * (delta_t1 * delta_t1 + delta_t2 * delta_t2));

  auto zero = Scalar(0.0, delta.options());
  auto one = Scalar(1.0, delta.options());

  if (_irreversible_damage)
  {
    auto delta_max_old = (*_delta_max_old)();
    delta_eff = neml2::where(delta_eff > delta_max_old, delta_eff, delta_max_old);
    if (out)
      _delta_max = Tensor(delta_eff);
  }

  auto d = one - exp(-delta_eff / _delta0);
  auto c = _Gc / (_delta0 * _delta0);

  if (out)
  {
    _damage = Tensor(d);
    _t = Tensor((one - d) * c * delta);
  }

  if (dout_din)
  {
    auto ddelta_eff_ddelta = Vec::fill(delta_n, _beta * delta_t1, _beta * delta_t2) / delta_eff;
    if (_irreversible_damage)
    {
      auto delta_max_old = (*_delta_max_old)();
      ddelta_eff_ddelta = neml2::where(delta_eff > delta_max_old, one, zero) * ddelta_eff_ddelta;
    }
    auto dd_ddelta = (one / _delta0) * exp(-delta_eff / _delta0) * ddelta_eff_ddelta;

    if (_delta.is_dependent())
    {
      _t.d(_delta) = Tensor((one - d) * c * R2::identity(_delta.options()));
      _t.d(_delta) += Tensor(-c * outer(delta, dd_ddelta));
      _damage.d(_delta) = Tensor(dd_ddelta);
      _delta_max.d(_delta) = Tensor(ddelta_eff_ddelta);
    }
  }
}
} // namespace neml2
