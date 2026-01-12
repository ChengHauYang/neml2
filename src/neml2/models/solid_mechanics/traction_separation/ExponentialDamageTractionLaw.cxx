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

#include "neml2/models/solid_mechanics/traction_separation/ExponentialDamageTractionLaw.h"
#include "neml2/tensors/functions/vdot.h"
#include "neml2/tensors/functions/norm.h"

namespace neml2
{

register_NEML2_object(ExponentialDamageTractionLaw);

OptionSet
ExponentialDamageTractionLaw::expected_options()
{
  OptionSet options = CohesiveTractionLaw::expected_options();
  options.doc() += " Exponential traction–separation law with effective separation and damage.";

  options.set_parameter<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Interface fracture energy";

  options.set_parameter<TensorName<Scalar>>("tangential_opening_weight");
  options.set("tangential_opening_weight").doc() =
      "Weights tangential opening relative to normal opening.";

  options.set_parameter<TensorName<Scalar>>("softening_length_scale");
  options.set("softening_length_scale").doc() =
      "Softening length scale controlling damage evolution.";

  return options;
}

ExponentialDamageTractionLaw::ExponentialDamageTractionLaw(const OptionSet & options)
  : CohesiveTractionLaw(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _beta(declare_parameter<Scalar>("beta", "tangential_opening_weight")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length_scale"))
{
}

void
ExponentialDamageTractionLaw::set_value(bool out, bool /*dout_din*/, bool /**/)
{
  const auto g = _jump();
  const auto n = _normal();

  auto delta_n = neml2::vdot(g, n);

  auto g_t = g - delta_n * n;

  auto delta_t = neml2::norm(g_t);

  auto delta_eff = sqrt(delta_n * delta_n + _beta * delta_t * delta_t);

  auto damage = 1.0 - exp(-delta_eff / _delta0);

  auto prefactor = -_Gc / (_delta0 * _delta0) * (1.0 - damage);

  if (out)
  {
    // Treat prefactor and g as ATen tensors and multiply them
    const at::Tensor tprod =
        static_cast<const at::Tensor &>(prefactor) * static_cast<const at::Tensor &>(g);

    // Reconstruct neml2::Vec with the same dynamic/intermediate dimensions as g
    _traction = neml2::Vec(tprod, g.dynamic_dim(), g.intmd_dim());
  }
}

} // namespace neml2
