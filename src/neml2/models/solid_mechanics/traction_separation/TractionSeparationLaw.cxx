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

#include "neml2/models/solid_mechanics/traction_separation/TractionSeparationLaw.h"

namespace neml2
{
OptionSet
TractionSeparationLaw::expected_options()
{
  OptionSet options = Model::expected_options();
  options.doc() = "Base interface for cohesive traction-separation laws mapping a displacement "
                  "jump to a traction.";

  options.set_input("displacement_jump") = VariableName(FORCES, "interface_displacement_jump");
  options.set("displacement_jump").doc() = "Displacement jump in the local interface frame \\f$ "
                                           "[\\delta_n, \\delta_{t1}, \\delta_{t2}] \\f$";

  options.set_output("traction") = VariableName(STATE, "interface_traction");
  options.set("traction").doc() =
      "Interface traction in the local interface frame \\f$ [T_n, T_{t1}, T_{t2}] \\f$";

  return options;
}

TractionSeparationLaw::TractionSeparationLaw(const OptionSet & options)
  : Model(options),
    _jump(declare_input_variable<Vec>("displacement_jump")),
    _traction(declare_output_variable<Vec>("traction"))
{
}

void
TractionSeparationLaw::diagnose() const
{
  Model::diagnose();
  diagnostic_assert_force(_jump);
  diagnostic_assert_state(_traction);
}
} // namespace neml2
