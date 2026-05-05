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

#include <cmath>

namespace neml2
{
register_NEML2_object(SalehaniIrani3DCTraction);

OptionSet
SalehaniIrani3DCTraction::expected_options()
{
  OptionSet options = TractionSeparation::expected_options();
  options.doc() +=
      " Salehani-Irani exponential coupling: normal component is linear-times-exponential and "
      "shear components are linear-times-exponential with a quadratic exponent contribution.";

  options.set_parameter<TensorName<Scalar>>("normal_gap_at_maximum_normal_traction");
  options.set("normal_gap_at_maximum_normal_traction").doc() =
      "Normal displacement jump at which the normal traction is maximum";

  options.set_parameter<TensorName<Scalar>>("tangential_gap_at_maximum_shear_traction");
  options.set("tangential_gap_at_maximum_shear_traction").doc() =
      "Tangential displacement jump at which the shear traction is maximum (raw input; an "
      "internal sqrt(2) scaling is applied)";

  options.set_parameter<TensorName<Scalar>>("maximum_normal_traction");
  options.set("maximum_normal_traction").doc() = "Peak normal traction";

  options.set_parameter<TensorName<Scalar>>("maximum_shear_traction");
  options.set("maximum_shear_traction").doc() = "Peak shear traction";

  return options;
}

SalehaniIrani3DCTraction::SalehaniIrani3DCTraction(const OptionSet & options)
  : TractionSeparation(options),
    _du0_n_raw(declare_parameter<Scalar>("du0_n", "normal_gap_at_maximum_normal_traction")),
    _du0_t_raw(declare_parameter<Scalar>("du0_t", "tangential_gap_at_maximum_shear_traction")),
    _Tn_max(declare_parameter<Scalar>("Tn_max", "maximum_normal_traction")),
    _Tt_max(declare_parameter<Scalar>("Tt_max", "maximum_shear_traction"))
{
}

void
SalehaniIrani3DCTraction::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  const auto sqrt2 = std::sqrt(2.0);
  const auto e1 = std::exp(1.0);
  const auto sqrt_2e = std::sqrt(2.0 * e1);

  const auto du0_n = _du0_n_raw;
  const auto du0_t = sqrt2 * _du0_t_raw;

  const auto a_n = e1 * _Tn_max;
  const auto a_t = sqrt_2e * _Tt_max;

  const auto dn = _delta()(0);
  const auto ds1 = _delta()(1);
  const auto ds2 = _delta()(2);

  const auto bn = dn / du0_n;
  const auto bs1 = ds1 / du0_t;
  const auto bs2 = ds2 / du0_t;

  const auto x = bn + bs1 * bs1 + bs2 * bs2;
  const auto exp_x = neml2::exp(-x);

  if (out)
    _T = Vec::fill(a_n * bn * exp_x, a_t * bs1 * exp_x, a_t * bs2 * exp_x);

  if (dout_din)
  {
    // dx/d(delta_j)
    const auto dx_dn = 1.0 / du0_n;
    const auto dx_ds1 = 2.0 * ds1 / (du0_t * du0_t);
    const auto dx_ds2 = 2.0 * ds2 / (du0_t * du0_t);

    // db_i/d(delta_j) is 1/du0_i on the diagonal, 0 off-diagonal
    const auto inv_du0_n = 1.0 / du0_n;
    const auto inv_du0_t = 1.0 / du0_t;

    // Row 0 (normal traction)
    const auto j00 = a_n * exp_x * (inv_du0_n - bn * dx_dn);
    const auto j01 = a_n * exp_x * (-bn * dx_ds1);
    const auto j02 = a_n * exp_x * (-bn * dx_ds2);

    // Row 1 (shear-1 traction)
    const auto j10 = a_t * exp_x * (-bs1 * dx_dn);
    const auto j11 = a_t * exp_x * (inv_du0_t - bs1 * dx_ds1);
    const auto j12 = a_t * exp_x * (-bs1 * dx_ds2);

    // Row 2 (shear-2 traction)
    const auto j20 = a_t * exp_x * (-bs2 * dx_dn);
    const auto j21 = a_t * exp_x * (-bs2 * dx_ds1);
    const auto j22 = a_t * exp_x * (inv_du0_t - bs2 * dx_ds2);

    _T.d(_delta) = R2::fill(j00, j01, j02, j10, j11, j12, j20, j21, j22);
  }
}
} // namespace neml2
