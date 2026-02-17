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

#include <ATen/ATen.h>

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

  options.set<bool>("irreversible") = false;
  options.set("irreversible").doc() =
      "If true, enforce damage irreversibility using d = max(d, d_old).";

  // Old damage is only required when irreversibility is enabled, but we still define the option
  // here so the name/mount is standardized.
  options.set_input("damage_old") = VariableName(OLD_STATE, "damage");
  options.set("damage_old").doc() = "Old damage value used when irreversibility is enabled.";

  // Damage is ALWAYS produced by this model, so we always declare it as a state output so the
  // driver can store it and provide old_state/damage next step.
  options.set_output("damage") = VariableName(STATE, "damage");
  options.set("damage").doc() = "Current damage value.";

  return options;
}

ExponentialDamageTractionLaw::ExponentialDamageTractionLaw(const OptionSet & options)
  : CohesiveTractionLaw(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _beta(declare_parameter<Scalar>("beta", "tangential_opening_weight")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length_scale")),
    _irreversible(options.get<bool>("irreversible")),
    _damage_old(_irreversible ? &declare_input_variable<Scalar>("damage_old") : nullptr),
    _damage(declare_output_variable<Scalar>("damage"))
{
}

void
ExponentialDamageTractionLaw::set_value(bool out, bool /*dout_din*/, bool /**/)
{
  // _jump() is in the local CZM frame:
  //   g(0) = normal component, g(1), g(2) = tangential components
  const auto g = _jump();

  const auto delta_n = g(0);
  const auto delta_t = sqrt(g(1) * g(1) + g(2) * g(2));

  auto delta_eff = sqrt(delta_n * delta_n + _beta * delta_t * delta_t);

  at::Tensor damage_t = 1.0 - at::exp(-static_cast<const at::Tensor &>(delta_eff) /
                                      static_cast<const at::Tensor &>(_delta0));

  // Optional irreversibility: never allow damage to decrease.
  if (_irreversible)
  {
    const at::Tensor dold = static_cast<const at::Tensor &>((*_damage_old)());
    damage_t = at::maximum(damage_t, dold);
  }

  const Scalar damage(damage_t, g.dynamic_dim(), g.intmd_dim());

  const Scalar prefactor = -_Gc / (_delta0 * _delta0) * (1.0 - damage);

  if (out)
  {
    _damage = damage;          // always store
    _traction = prefactor * g; // NEML2 broadcasting

#ifdef DEBUG
    // taction output
    std::cout << "---traction output---------\n";
    std::cout << "base_dim        = " << _traction.base_dim() << "\n";
    std::cout << "dynamic_dim     = " << _traction.dynamic_dim() << "\n";
    std::cout << "intmd_dim       = " << _traction.intmd_dim() << "\n";

    for (Size i = 0; i < _traction.base_dim(); ++i)
      std::cout << "base_size[" << i << "] = " << _traction.base_size(i) << "\n";

    for (Size i = 0; i < _traction.dynamic_dim(); ++i)
      std::cout << "dynamic_size[" << i << "] = " << _traction.dynamic_size(i) << "\n";

    for (Size i = 0; i < _traction.intmd_dim(); ++i)
      std::cout << "intmd_size[" << i << "] = " << _traction.intmd_size(i) << "\n";

    std::cout << "base_sizes      = " << _traction.base_sizes() << "\n";
    std::cout << "dynamic_sizes   = " << _traction.dynamic_sizes() << "\n";
    std::cout << "intmd_sizes     = " << _traction.intmd_sizes() << "\n";

    // damage output
    std::cout << "---damage output---------\n";
    std::cout << "base_dim        = " << _damage.base_dim() << "\n";
    std::cout << "dynamic_dim     = " << _damage.dynamic_dim() << "\n";
    std::cout << "intmd_dim       = " << _damage.intmd_dim() << "\n";
    for (Size i = 0; i < _damage.base_dim(); ++i)
      std::cout << "base_size[" << i << "] = " << _damage.base_size(i) << "\n";

    for (Size i = 0; i < _damage.dynamic_dim(); ++i)
      std::cout << "dynamic_size[" << i << "] = " << _damage.dynamic_size(i) << "\n";

    for (Size i = 0; i < _damage.intmd_dim(); ++i)
      std::cout << "intmd_size[" << i << "] = " << _damage.intmd_size(i) << "\n";
    std::cout << "base_sizes      = " << _damage.base_sizes() << "\n";
    std::cout << "dynamic_sizes   = " << _damage.dynamic_sizes() << "\n";
    std::cout << "intmd_sizes     = " << _damage.intmd_sizes() << "\n";
#endif
  }
}

} // namespace neml2
