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

#include "neml2/models/solid_mechanics/cohesive_zone/TractionSeparation.h"

namespace neml2
{
class Scalar;

/**
 * @brief Bilinear mixed-mode traction separation law.
 *
 * See design/traction_separation_law/BiLinearMixedModeTraction_spec_simple.md
 */
class BiLinearMixedModeTraction : public TractionSeparation
{
public:
  static OptionSet expected_options();

  BiLinearMixedModeTraction(const OptionSet & options);

protected:
  void set_value(bool out, bool dout_din, bool d2out_din2) override;

  /// Penalty elastic stiffness
  const Scalar & _K;
  /// Mode I critical energy release rate
  const Scalar & _GIc;
  /// Mode II critical energy release rate
  const Scalar & _GIIc;
  /// Tensile strength
  const Scalar & _N;
  /// Shear strength
  const Scalar & _S;
  /// Power law exponent
  const Scalar & _eta;

  /// Mixed mode criterion: BK or POWER_LAW
  const std::string _criterion;

  /// Regularization parameter for Macaulay bracket
  const double _alpha;

  /// Damage (state)
  Variable<Scalar> & _d;

  /// Damage (old state)
  const Variable<Scalar> & _d_old;
};
} // namespace neml2
