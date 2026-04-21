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

#include "neml2/models/solid_mechanics/traction_separation/TractionSeparation.h"

namespace neml2
{
OptionSet
TractionSeparation::expected_options()
{
  OptionSet options = Model::expected_options();

  options.set_input("interface_displacement_jump") = VariableName(FORCES, "interface_displacement_jump");
  options.set("interface_displacement_jump").doc() = "Interface displacement jump";

  options.set_output("interface_traction") = VariableName(FORCES, "interface_traction");
  options.set("interface_traction").doc() = "Interface traction";

  return options;
}

TractionSeparation::TractionSeparation(const OptionSet & options)
  : Model(options),
    _jump(declare_input_variable<Vec>("interface_displacement_jump")),
    _traction(declare_output_variable<Vec>("interface_traction"))
{
}
} // namespace neml2
