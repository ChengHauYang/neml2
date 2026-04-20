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
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/outer.h"

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() = "Salehani & Irani (2018) 3D coupled traction-separation law.";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction");
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Normal gap at which maximum normal traction is reached";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction");
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Tangential gap at which maximum shear traction is reached";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Maximum normal traction";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Maximum shear traction";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparation(options),
    _delta_un_max(declare_parameter<Scalar>("delta_un_max", "normal_gap_at_maximum_normal_traction")),
    _delta_us_max(
        declare_parameter<Scalar>("delta_us_max", "tangential_gap_at_maximum_shear_traction")),
    _tn_max(declare_parameter<Scalar>("tn_max", "maximum_normal_traction")),
    _ts_max(declare_parameter<Scalar>("ts_max", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  auto options = _displacement_jump.options();
  auto dn0 = _delta_un_max;
  auto ds0 = neml2::sqrt(Scalar::fill(2.0, options)) * _delta_us_max;

  auto a0 = neml2::exp(Scalar::fill(1.0, options)) * _tn_max;
  auto a12 = neml2::sqrt(Scalar::fill(2.0, options) * neml2::exp(Scalar::fill(1.0, options))) * _ts_max;

  auto dn = _displacement_jump()(0);
  auto ds1 = _displacement_jump()(1);
  auto ds2 = _displacement_jump()(2);

  auto x = (dn / dn0) + (ds1 / ds0) * (ds1 / ds0) + (ds2 / ds0) * (ds2 / ds0);
  auto exp_x = neml2::exp(-x);

  if (out)
  {
    auto tn = a0 * (dn / dn0) * exp_x;
    auto ts1 = a12 * (ds1 / ds0) * exp_x;
    auto ts2 = a12 * (ds2 / ds0) * exp_x;
    _traction = Vec::fill(tn, ts1, ts2);
  }

  if (dout_din)
  {
    if (_displacement_jump.is_dependent())
    {
      auto dx_du0 = Scalar::fill(1.0, options) / dn0;
      auto dx_du1 = Scalar::fill(2.0, options) * ds1 / (ds0 * ds0);
      auto dx_du2 = Scalar::fill(2.0, options) * ds2 / (ds0 * ds0);
      auto dx_du = Vec::fill(dx_du0, dx_du1, dx_du2);

      auto b0 = dn / dn0;
      auto b1 = ds1 / ds0;
      auto b2 = ds2 / ds0;
      auto b = Vec::fill(b0, b1, b2);

      auto a = Vec::fill(a0, a12, a12);

      auto db_du = R2::fill(Scalar::fill(1.0, options) / dn0, Scalar::fill(1.0, options) / ds0, Scalar::fill(1.0, options) / ds0);

      // dTi_duj = a_i * exp_x * (dbi_duj - b_i * dx_duj)
      // dT_du = diag(a) * exp_x * (db_du - outer(b, dx_du))
      _traction.d(_displacement_jump) =
          R2::fill(a(0), a(1), a(2)) * exp_x * (db_du - neml2::outer(b, dx_du));
    }
  }
}
} // namespace neml2
