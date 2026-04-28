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

#pragma once

#include "neml2/models/solid_mechanics/traction_separation/TractionSeparation.h"

namespace neml2
{
class Scalar;

/**
 * @brief Salehani-Irani 3D coupled exponential traction-separation law.
 *
 * Couples normal (\f$\alpha=1\f$) and shear (\f$\alpha=2\f$) responses through a single
 * exponential factor \f$ \exp(-x) \f$ with
 * \f$ x = \delta_n/\delta_{u,n} + (\delta_{s1}/\delta_{u,t})^2 + (\delta_{s2}/\delta_{u,t})^2 \f$.
 * The shear characteristic length stored internally is \f$ \sqrt{2} \f$ times the user-supplied
 * tangential gap.
 */
class SalehaniIrani3DCTraction : public TractionSeparation
{
public:
  static OptionSet expected_options();

  SalehaniIrani3DCTraction(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Normal characteristic gap (at maximum normal traction)
  const Scalar & _du0_n_raw;

  /// Tangential characteristic gap (at maximum shear traction); raw user input
  const Scalar & _du0_t_raw;

  /// Maximum normal traction
  const Scalar & _Tn_max;

  /// Maximum shear traction
  const Scalar & _Tt_max;
};
} // namespace neml2
