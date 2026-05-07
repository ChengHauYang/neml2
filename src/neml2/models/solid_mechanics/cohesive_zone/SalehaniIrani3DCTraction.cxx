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

#include "neml2/models/solid_mechanics/cohesive_zone/SalehaniIrani3DCTraction.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() +=
      " Implements the Salehani-Irani 3D coupled exponential law (Salehani & Irani 2018), "
      "in which normal opening enters linearly and the two tangential components enter "
      "quadratically in a single shared exponential envelope.";

  options.add_parameter<Scalar>("normal_gap_at_maximum_normal_traction",
                                "Characteristic normal gap at the peak normal traction.");
  options.add_parameter<Scalar>("tangential_gap_at_maximum_shear_traction",
                                "Characteristic tangential gap at the peak shear traction "
                                "(raw, unscaled).");
  options.add_parameter<Scalar>("maximum_normal_traction", "Peak normal traction.");
  options.add_parameter<Scalar>("maximum_shear_traction", "Peak shear traction.");

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparation(options),
    _delta0_n(declare_parameter<Scalar>("delta0_n", "normal_gap_at_maximum_normal_traction")),
    _delta0_s(declare_parameter<Scalar>("delta0_s", "tangential_gap_at_maximum_shear_traction")),
    _Tmax_n(declare_parameter<Scalar>("Tmax_n", "maximum_normal_traction")),
    _Tmax_s(declare_parameter<Scalar>("Tmax_s", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // TODO: implement the coupled exponential traction. From the spec:
  //   x       = delta_n / delta0_n
  //           + (delta_s1 / (sqrt(2) * delta0_s))^2
  //           + (delta_s2 / (sqrt(2) * delta0_s))^2
  //   T_n     = e        * Tmax_n * (delta_n  / delta0_n)              * exp(-x)
  //   T_si    = sqrt(2e) * Tmax_s * (delta_si / (sqrt(2) * delta0_s))  * exp(-x)
  // Note the sqrt(2) scaling is baked into the shear characteristic gap in MOOSE;
  // decide whether to keep that here or expose the scaled value as the parameter.
  if (out)
    _T = _delta();

  // TODO: implement the full 3x3 Jacobian (off-diagonal terms come from the shared exp(-x)).
  if (dout_din)
  {
    // _T.d(_delta) = ...;
  }
}
} // namespace neml2
