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

#include "neml2/models/solid_mechanics/traction_separation/ExpTractionSeparation.h"
#include "neml2/tensors/R2.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/exp.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/where.h"

namespace neml2
{
register_NEML2_object(ExpTractionSeparation);

OptionSet
ExpTractionSeparation::expected_options()
{
  OptionSet options = TractionSeparationLaw::expected_options();
  options.doc() +=
      " Exponential damage traction-separation law with optional damage irreversibility.";

  options.set_parameter<TensorName<Scalar>>("fracture_energy");
  options.set("fracture_energy").doc() = "Interface fracture energy Gc";

  options.set_parameter<TensorName<Scalar>>("softening_length_scale");
  options.set("softening_length_scale").doc() = "Softening length scale delta0";

  options.set_parameter<TensorName<Scalar>>("tangential_opening_weight");
  options.set("tangential_opening_weight").doc() =
      "Weight beta >= 0 scaling tangential vs. normal effective opening";

  options.set<bool>("irreversible") = false;
  options.set("irreversible").doc() =
      "If true, enforce damage irreversibility: d = max(d, d_old)";

  options.set_input("damage_old") = VariableName(OLD_STATE, "damage");
  options.set("damage_old").doc() = "Old damage value (required when irreversible = true)";

  options.set_output("damage") = VariableName(STATE, "damage");
  options.set("damage").doc() = "Scalar damage d in [0, 1)";

  return options;
}

ExpTractionSeparation::ExpTractionSeparation(const OptionSet & options)
  : TractionSeparationLaw(options),
    _Gc(declare_parameter<Scalar>("Gc", "fracture_energy")),
    _delta0(declare_parameter<Scalar>("delta0", "softening_length_scale")),
    _beta(declare_parameter<Scalar>("beta", "tangential_opening_weight")),
    _irreversible(options.get<bool>("irreversible")),
    _damage_old(_irreversible ? &declare_input_variable<Scalar>("damage_old") : nullptr),
    _damage(declare_output_variable<Scalar>("damage"))
{
}

void
ExpTractionSeparation::set_value(bool out, bool dout_din, bool d2out_din2)
{
  const auto g = _jump();
  const Scalar dn = g(0);
  const Scalar dt1 = g(1);
  const Scalar dt2 = g(2);

  // Effective displacement jump: sqrt(dn^2 + beta*(dt1^2 + dt2^2))
  const Scalar delta_eff_sq = dn * dn + _beta * (dt1 * dt1 + dt2 * dt2);
  const Scalar delta_eff = sqrt(delta_eff_sq);

  // Damage from exponential law: d = 1 - exp(-delta_eff / delta0)
  const Scalar d_trial = 1.0 - exp(-delta_eff / _delta0);

  // Irreversibility: d = max(d_trial, d_old)
  // Precompute d_old and frozen_mask here so they are available for both out and dout_din
  Scalar d_old = Scalar::zeros_like(d_trial);
  Scalar frozen_mask = d_trial < Scalar::zeros_like(d_trial); // boolean false-like placeholder
  Scalar damage = d_trial;
  if (_irreversible)
  {
    d_old = (*_damage_old)();
    frozen_mask = d_trial < d_old;
    damage = where(frozen_mask, d_old, d_trial);
  }

  // Traction prefactor: c = Gc / delta0^2
  const Scalar c = _Gc / (_delta0 * _delta0);

  if (out)
  {
    _damage = damage;
    _traction = (1.0 - damage) * c * g;
  }

  if (dout_din)
  {
    if (_jump.is_dependent())
    {
      // dd/dg = exp(-delta_eff/delta0) / (delta0 * delta_eff) * (dn, beta*dt1, beta*dt2)
      const Scalar exp_term = exp(-delta_eff / _delta0) / (_delta0 * delta_eff);
      Vec dd_dg = exp_term * Vec::fill(dn, _beta * dt1, _beta * dt2);

      if (_irreversible)
        dd_dg = where(frozen_mask, Vec::zeros_like(dd_dg), dd_dg);

      // dt/dg = (1-d)*c*I - c*outer(g, dd/dg)
      _damage.d(_jump) = dd_dg;
      _traction.d(_jump) = (1.0 - damage) * c * R2::identity(_jump.options()) +
                           (-c) * outer(g, dd_dg);
    }

    if (_irreversible && _damage_old->is_dependent())
    {
      // d(damage)/d(d_old) = 1 where frozen, 0 where advancing
      const Scalar dd_dd_old =
          where(frozen_mask, Scalar::ones_like(damage), Scalar::zeros_like(damage));
      _damage.d(*_damage_old) = dd_dd_old;
      // d(traction)/d(d_old) = -c * g * d(damage)/d(d_old)
      _traction.d(*_damage_old) = (-c) * dd_dd_old * g;
    }
  }

  if (d2out_din2)
  {
    // not implemented
  }
}
} // namespace neml2
