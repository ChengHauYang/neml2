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

#include "neml2/models/solid_mechanics/traction_separation/SalehaniIrani3DCTraction.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/outer.h"

#include <cmath>

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() +=
      " Salehani-Irani (2018) 3D exponential traction-separation law with fully coupled Jacobian.";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction");
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Normal gap at which the normal traction is maximum";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction");
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Shear gap at which the shear traction is maximum (stored internally as sqrt(2) times this "
      "value)";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Peak normal traction";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Peak shear traction";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparationLaw(options),
    _delta_n0(declare_parameter<Scalar>("delta_n0", "normal_gap_at_maximum_normal_traction")),
    _delta_t0(declare_parameter<Scalar>("delta_t0", "tangential_gap_at_maximum_shear_traction")),
    _T_n_max(declare_parameter<Scalar>("T_n_max", "maximum_normal_traction")),
    _T_s_max(declare_parameter<Scalar>("T_s_max", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, bool d2out_din2)
{
  // Internal characteristic gaps; shear is scaled by sqrt(2) per the Salehani-Irani formulation
  const Scalar dn0 = _delta_n0;
  const Scalar ds0 = std::sqrt(2.0) * _delta_t0;

  // Traction prefactors: a_n = e*T_n_max, a_s = sqrt(2*e)*T_s_max
  const double e = std::exp(1.0);
  const Scalar a_n = e * _T_n_max;
  const Scalar a_s = std::sqrt(2.0 * e) * _T_s_max;

  const auto g = _jump();
  const Scalar dn = g(0);
  const Scalar ds1 = g(1);
  const Scalar ds2 = g(2);

  // Normalized components
  const Scalar bn = dn / dn0;
  const Scalar bs1 = ds1 / ds0;
  const Scalar bs2 = ds2 / ds0;

  // Exponent: x = bn + bs1^2 + bs2^2  (normal linear, shear quadratic)
  const Scalar x = bn + bs1 * bs1 + bs2 * bs2;
  const Scalar ex = exp(-x);

  // traction = a_i * b_i * exp(-x)
  const Vec t = Vec::fill(a_n * bn, a_s * bs1, a_s * bs2) * ex;

  if (out)
    _traction = t;

  if (dout_din)
  {
    if (_jump.is_dependent())
    {
      // J_ij = a_i * ex * (dbi/dg_j - bi * dx/dg_j)
      //      = ex * diag(a_n/dn0, a_s/ds0, a_s/ds0) - outer(t, dx/dg)
      const Vec dx_dg = Vec::fill(Scalar::full_like(dn, 1.0) / dn0,
                                  2.0 * ds1 / (ds0 * ds0),
                                  2.0 * ds2 / (ds0 * ds0));
      const R2 J_diag = R2::fill(a_n / dn0, a_s / ds0, a_s / ds0) * ex;
      _traction.d(_jump) = J_diag + (-outer(t, dx_dg));
    }
  }

  if (d2out_din2)
  {
    // not implemented
  }
}
} // namespace neml2
