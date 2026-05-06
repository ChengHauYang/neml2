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
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/R2.h"
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
  options.doc() +=
      " Exponential cohesive form with irreversible scalar damage driven by the historic maximum "
      "effective displacement jump.";

  options.set_parameter<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Fracture energy G_c";

  options.set_parameter<TensorName<Scalar>>("softening_length_scale");
  options.set("softening_length_scale").doc() = "Softening length scale delta_0";

  options.set_parameter<TensorName<Scalar>>("tangential_weighting_factor");
  options.set("tangential_weighting_factor").doc() =
      "Weighting factor beta on the squared tangential jump";

  options.set<double>("epsilon") = 1e-16;
  options.set("epsilon").doc() =
      "Small regularizer added inside the effective-jump sqrt to keep its "
      "derivative bounded at zero jump";

  options.set_output("effective_displacement_jump_max") =
      VariableName(STATE, "effective_displacement_jump_max");
  options.set("effective_displacement_jump_max").doc() =
      "Historic maximum effective scalar displacement jump";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparation(options),
    _delta_eff_max(declare_output_variable<Scalar>("effective_displacement_jump_max")),
    _delta_eff_max_old(declare_input_variable<Scalar>(_delta_eff_max.name().old())),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length_scale")),
    _beta(declare_parameter<Scalar>("beta", "tangential_weighting_factor")),
    _eps(options.get<double>("epsilon"))
{
}

void
ExpTractionSeparation::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto dn = _delta()(0);
  const auto ds1 = _delta()(1);
  const auto ds2 = _delta()(2);

  const auto delta_eff = neml2::sqrt(dn * dn + _beta * (ds1 * ds1 + ds2 * ds2) + _eps);
  const auto adv_mask = delta_eff > _delta_eff_max_old();
  const auto delta_eff_max = neml2::where(adv_mask, delta_eff, _delta_eff_max_old());

  const auto exp_term = neml2::exp(-delta_eff_max / _delta0);
  const auto d = 1.0 - exp_term;
  const auto c = _Gc / (_delta0 * _delta0);

  if (out)
  {
    _delta_eff_max = delta_eff_max;
    _T = (1.0 - d) * c * _delta();
  }

  if (dout_din)
  {
    const auto zero_v = Vec::zeros_like(_delta());
    const auto zero_s = Scalar::zeros_like(delta_eff);
    const auto one_s = Scalar::ones_like(delta_eff);

    // d(delta_eff)/d(jump) = (1/delta_eff) * (dn, beta*ds1, beta*ds2)
    const auto d_delta_eff_d_jump =
        Vec::fill(dn, _beta * ds1, _beta * ds2) / delta_eff;
    const auto d_delta_eff_max_d_jump =
        neml2::where(adv_mask, d_delta_eff_d_jump, zero_v);
    const auto d_delta_eff_max_d_old = neml2::where(adv_mask, zero_s, one_s);

    // d(d)/d(delta_eff_max) = (1/delta0) * exp(-delta_eff_max / delta0)
    const auto d_d_d_delta_eff_max = exp_term / _delta0;
    const auto d_d_d_jump = d_d_d_delta_eff_max * d_delta_eff_max_d_jump;
    const auto d_d_d_old = d_d_d_delta_eff_max * d_delta_eff_max_d_old;

    // dT/d(jump) = (1-d)*c*I - c * outer(jump, d_d/d_jump)
    const auto I = R2::identity(_delta.options());
    _T.d(_delta) = (1.0 - d) * c * I - c * neml2::outer(_delta(), d_d_d_jump);

    // dT/d(delta_eff_max_old) = -c * d_d/d_old * jump
    _T.d(_delta_eff_max_old) = -c * d_d_d_old * _delta();

    _delta_eff_max.d(_delta) = d_delta_eff_max_d_jump;
    _delta_eff_max.d(_delta_eff_max_old) = d_delta_eff_max_d_old;
  }
}
} // namespace neml2
