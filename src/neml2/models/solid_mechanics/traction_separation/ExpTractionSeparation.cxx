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
#include "neml2/tensors/functions/where.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/misc/types.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Exponential softening law for cohesive zone modeling.";

  options.set_parameter<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Fracture energy G_c";

  options.set_parameter<TensorName<Scalar>>("softening_length");
  options.set("softening_length").doc() = "Softening length scale delta_0";

  options.set_parameter<TensorName<Scalar>>("tangential_weight") = TensorName<Scalar>("1.0");
  options.set("tangential_weight").doc() = "Tangential weighting factor beta";

  options.set_parameter<TensorName<Scalar>>("epsilon") = TensorName<Scalar>("1e-16");
  options.set("epsilon").doc() = "Small regularizer for norm";

  options.set<bool>("irreversible_damage") = true;
  options.set("irreversible_damage").doc() = "If true, use historical maximum scalar jump";

  options.set_output("damage") = VariableName(STATE, "internal", "d");
  options.set("damage").doc() = "Damage variable";

  options.set_output("effective_jump_max") = VariableName(STATE, "internal", "delta_eff_max");
  options.set("effective_jump_max").doc() = "Maximum effective scalar jump";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length")),
    _beta(declare_parameter<Scalar>("beta", "tangential_weight")),
    _eps(declare_parameter<Scalar>("eps", "epsilon")),
    _irreversible_damage(options.get<bool>("irreversible_damage")),
    _delta_eff_max(declare_output_variable<Scalar>("effective_jump_max")),
    _delta_eff_max_old(declare_input_variable<Scalar>(_delta_eff_max.name().old())),
    _damage(declare_output_variable<Scalar>("damage"))
{
}

void
ExpTractionSeparation::request_AD()
{
  _traction.request_AD(_displacement_jump);
  _delta_eff_max.request_AD(_displacement_jump);
  _damage.request_AD(_displacement_jump);
}

void
ExpTractionSeparation::set_value(bool out, bool /*dout_din*/, bool /*d2out_din2*/)
{
  if (!out)
    return;

  auto delta_n = _displacement_jump()(0);
  auto delta_t1 = _displacement_jump()(1);
  auto delta_t2 = _displacement_jump()(2);
  auto delta_t_sq = delta_t1 * delta_t1 + delta_t2 * delta_t2;

  auto delta_eff = neml2::sqrt(delta_n * delta_n + _beta * delta_t_sq + _eps);

  if (_irreversible_damage)
  {
    _delta_eff_max = neml2::where(delta_eff > _delta_eff_max_old(), delta_eff, _delta_eff_max_old());
    delta_eff = _delta_eff_max();
  }

  _damage = 1.0 - neml2::exp(-delta_eff / _delta0);
  auto c = _Gc / (_delta0 * _delta0);
  _traction = (1.0 - _damage) * c * _displacement_jump;
}
} // namespace neml2
