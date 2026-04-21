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
#include "neml2/tensors/functions/outer.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Exponential traction-separation law";

  options.set_parameter<TensorName<Scalar>>("fracture_energy") = "1.0";
  options.set("fracture_energy").doc() = "Fracture energy G_c";

  options.set_parameter<TensorName<Scalar>>("softening_length") = "1.0";
  options.set("softening_length").doc() = "Softening length scale delta_0";

  options.set_parameter<TensorName<Scalar>>("weighting_factor") = "1.0";
  options.set("weighting_factor").doc() = "Tangential weighting factor beta";

  options.set_parameter<TensorName<Scalar>>("eps") = "1e-16";
  options.set("eps").doc() = "Small regularizer for norm";

  options.set<bool>("irreversible_damage") = true;
  options.set("irreversible_damage").doc() = "Whether to use historical max delta_eff";

  options.set_input("old_maximum_effective_displacement_jump") =
      VariableName(OLD_STATE, "internal", "delta_eff_max");
  options.set("old_maximum_effective_displacement_jump").doc() =
      "Maximum effective scalar jump from previous step";

  options.set_output("maximum_effective_displacement_jump") =
      VariableName(STATE, "internal", "delta_eff_max");
  options.set("maximum_effective_displacement_jump").doc() = "Maximum effective scalar jump";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _Gc(declare_parameter<Scalar>("fracture_energy", "fracture_energy", true)),
    _delta0(declare_parameter<Scalar>("softening_length", "softening_length", true)),
    _beta(declare_parameter<Scalar>("weighting_factor", "weighting_factor", true)),
    _eps(declare_parameter<Scalar>("eps", "eps", true)),
    _irreversible(options.get<bool>("irreversible_damage")),
    _delta_eff_max_old(declare_input_variable<Scalar>("old_maximum_effective_displacement_jump")),
    _delta_eff_max(declare_output_variable<Scalar>("maximum_effective_displacement_jump"))
{
}

void
ExpTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto dn = _jump()(0);
  const auto dt1 = _jump()(1);
  const auto dt2 = _jump()(2);
  const auto dt_sq = dt1 * dt1 + dt2 * dt2;
  const auto delta_eff_cur = neml2::sqrt(dn * dn + _beta * dt_sq + _eps);

  const auto delta_eff = _irreversible ? neml2::where(delta_eff_cur > _delta_eff_max_old(), delta_eff_cur, _delta_eff_max_old()) : delta_eff_cur;

  if (out)
  {
    _delta_eff_max = delta_eff;
    const auto d = 1.0 - neml2::exp(-delta_eff / _delta0);
    const auto c = _Gc / (_delta0 * _delta0);
    _traction = (1.0 - d) * c * _jump();
  }

  if (dout_din)
  {
    const auto c = _Gc / (_delta0 * _delta0);
    const auto exp_term = neml2::exp(-delta_eff / _delta0);
    const auto d = 1.0 - exp_term;

    if (_jump.is_dependent())
    {
      const auto ddelta_eff_cur_ddelta =
          Vec::fill(dn, _beta * dt1, _beta * dt2) / delta_eff_cur;

      Vec ddelta_eff_ddelta = ddelta_eff_cur_ddelta;
      if (_irreversible)
        ddelta_eff_ddelta = neml2::where(delta_eff_cur > _delta_eff_max_old(),
                                         ddelta_eff_cur_ddelta,
                                         Vec::zeros(_jump.options()));

      const auto dd_ddelta = exp_term * (1.0 / _delta0) * ddelta_eff_ddelta;

      _traction.d(_jump) =
          (1.0 - d) * c * R2::identity(_jump.options()) - c * neml2::outer(_jump(), dd_ddelta);

      _delta_eff_max.d(_jump) = ddelta_eff_ddelta;
    }
  }
}
} // namespace neml2
