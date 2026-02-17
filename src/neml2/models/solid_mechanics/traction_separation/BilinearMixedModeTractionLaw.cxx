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

#include "neml2/models/solid_mechanics/traction_separation/BilinearMixedModeTractionLaw.h"

#include "neml2/misc/assertions.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/norm.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sqrt.h"

#include <ATen/ATen.h>
#include <cmath>

namespace neml2
{

register_NEML2_object(BilinearMixedModeTractionLaw);

OptionSet
BilinearMixedModeTractionLaw::expected_options()
{
  OptionSet options = CohesiveTractionLaw::expected_options();
  options.doc() += " Mixed-mode bilinear traction-separation law with damage (BK or POWER_LAW).";

  // --- parameters ---
  options.set_parameter<TensorName<Scalar>>("penalty_stiffness");
  options.set("penalty_stiffness").doc() = "Penalty stiffness K.";

  options.set_parameter<TensorName<Scalar>>("GIc");
  options.set("GIc").doc() = "Critical energy release rate in normal direction.";

  options.set_parameter<TensorName<Scalar>>("GIIc");
  options.set("GIIc").doc() = "Critical energy release rate in shear direction.";

  options.set_parameter<TensorName<Scalar>>("normal_strength");
  options.set("normal_strength").doc() = "Normal (tensile) strength N.";

  options.set_parameter<TensorName<Scalar>>("shear_strength");
  options.set("shear_strength").doc() = "Shear strength S.";

  options.set_parameter<TensorName<Scalar>>("eta");
  options.set("eta").doc() = "Power law parameter for mixed-mode criterion.";

  options.set_parameter<TensorName<Scalar>>("viscosity") = TensorName<Scalar>("viscosity");
  options.set("viscosity").doc() = "Viscosity for viscous regularization (0 disables).";

  // MOOSE-aligned smoothing length for regularized heaviside
  options.set_parameter<TensorName<Scalar>>("alpha") = TensorName<Scalar>("alpha");
  options.set("alpha").doc() = "Smoothing length for regularized Heaviside (MOOSE style).";

  // criterion stored as a string option (BK or POWER_LAW)
  options.set<std::string>("mixed_mode_criterion") = "BK";
  options.set("mixed_mode_criterion").doc() = "Mixed-mode criterion: 'BK' or 'POWER_LAW'.";

  // lagging
  options.set<bool>("lag_mode_mixity") = true;
  options.set("lag_mode_mixity").doc() = "If true, use old jump to compute mode mixity ratio beta.";

  options.set<bool>("lag_displacement_jump") = false;
  options.set("lag_displacement_jump").doc() =
      "If true, use old jump to compute effective displacement jump delta_m.";

  // dt input for viscous regularization
  options.set_input("dt") = VariableName(STATE, "dt");
  options.set("dt").doc() = "Time step size (required if viscosity > 0).";

  // outputs
  options.set_output("damage") = VariableName(STATE, "damage");
  options.set("damage").doc() = "Scalar damage variable d in [0,1].";

  // Old state inputs
  options.set_input("damage_old") = VariableName(OLD_STATE, "damage");
  options.set("damage_old").doc() = "Old damage value for irreversibility.";

  options.set_input("jump_old") = VariableName(OLD_STATE, "jump");
  options.set("jump_old").doc() = "Old displacement jump for lagging (needed only if "
                                  "lag_mode_mixity or lag_displacement_jump is true).";

  return options;
}

BilinearMixedModeTractionLaw::BilinearMixedModeTractionLaw(const OptionSet & options)
  : CohesiveTractionLaw(options),

    _K(declare_parameter<Scalar>("K", "penalty_stiffness")),
    _GIc(declare_parameter<Scalar>("GIc", "GIc")),
    _GIIc(declare_parameter<Scalar>("GIIc", "GIIc")),
    _N(declare_parameter<Scalar>("N", "normal_strength")),
    _S(declare_parameter<Scalar>("S", "shear_strength")),
    _eta(declare_parameter<Scalar>("eta", "eta")),
    _visc(declare_parameter<Scalar>("visc", "viscosity")),
    _alpha(declare_parameter<Scalar>("alpha", "alpha")),

    _criterion(options.get<std::string>("mixed_mode_criterion")),
    _lag_mode_mixity(options.get<bool>("lag_mode_mixity")),
    _lag_disp_jump(options.get<bool>("lag_displacement_jump")),

    _dt(declare_input_variable<Scalar>("dt")),
    _jump_old((_lag_mode_mixity || _lag_disp_jump) ? &declare_input_variable<Vec>("jump_old")
                                                   : nullptr),
    _damage_old(declare_input_variable<Scalar>("damage_old")),
    _damage(declare_output_variable<Scalar>("damage"))
{
}

/**
 * if x <= -l: 0
 * else if x < l: 0.5*(1 + sin(pi*x/(2*l)))
 * else: 1
 */
Scalar
BilinearMixedModeTractionLaw::regularizedHeaviside(const Scalar & x,
                                                   const Scalar & smoothing_length) const
{
  const at::Tensor xt = static_cast<const at::Tensor &>(x);
  const at::Tensor lt = static_cast<const at::Tensor &>(smoothing_length);

  const at::Tensor zero = at::zeros_like(xt);
  const at::Tensor one = at::ones_like(xt);

  const at::Tensor mid = 0.5 * (1.0 + at::sin(M_PI * xt / (2.0 * lt)));
  const at::Tensor H = at::where(xt <= -lt, zero, at::where(xt < lt, mid, one));

  return Scalar(H, x.dynamic_dim(), x.intmd_dim());
}

/**
 * if |x| < l: 0.25*pi/l*cos(pi*x/(2*l))
 * else: 0
 */
Scalar
BilinearMixedModeTractionLaw::regularizedHeavisideDerivative(const Scalar & x,
                                                             const Scalar & smoothing_length) const
{
  const at::Tensor xt = static_cast<const at::Tensor &>(x);
  const at::Tensor lt = static_cast<const at::Tensor &>(smoothing_length);

  const at::Tensor zero = at::zeros_like(xt);
  const at::Tensor mid = 0.25 * (M_PI / lt) * at::cos(M_PI * xt / (2.0 * lt));
  const at::Tensor dH = at::where((xt < lt) & (xt > -lt), mid, zero);

  return Scalar(dH, x.dynamic_dim(), x.intmd_dim());
}

Scalar
BilinearMixedModeTractionLaw::macauley_pos(const Scalar & x, const Scalar & smoothing_length) const
{
  // <x>_+ = H(x) * x
  return regularizedHeaviside(x, smoothing_length) * x;
}

Scalar
BilinearMixedModeTractionLaw::maximum(const Scalar & a, const Scalar & b) const
{
  const at::Tensor ta = static_cast<const at::Tensor &>(a);
  const at::Tensor tb = static_cast<const at::Tensor &>(b);
  return Scalar(at::maximum(ta, tb), a.dynamic_dim(), a.intmd_dim());
}

Scalar
BilinearMixedModeTractionLaw::minimum(const Scalar & a, const Scalar & b) const
{
  const at::Tensor ta = static_cast<const at::Tensor &>(a);
  const at::Tensor tb = static_cast<const at::Tensor &>(b);
  return Scalar(at::minimum(ta, tb), a.dynamic_dim(), a.intmd_dim());
}

void
BilinearMixedModeTractionLaw::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // _jump() and _jump_old() are in the local CZM frame:
  //   g(0) = normal component, g(1), g(2) = tangential components
  // This matches MOOSE:
  //   interface_jump_local = R^T * interface_jump_global
  // So we must NOT re-project using _normal().

  // Step 1: compute mode mixity ratio, beta
  // Current jump, old jump (if lagging is disabled, old jump is not needed)
  const auto g = _jump();
  const auto g_old = _jump_old ? (*_jump_old)() : g;

  // Which jump used for beta and delta_m
  const auto delta = _lag_mode_mixity ? g_old : g;
  const auto g_eff = _lag_disp_jump ? g_old : g;

  // Mode mixity beta (using delta)
  const Scalar dn_beta = delta(0);
  const Scalar ds_beta = sqrt(delta(1) * delta(1) + delta(2) * delta(2));

  // beta = ds/dn if opening (dn>0), else 0
  const at::Tensor dn_bt = static_cast<const at::Tensor &>(dn_beta);
  const at::Tensor ds_bt = static_cast<const at::Tensor &>(ds_beta);
  const at::Tensor zero_t = at::zeros_like(dn_bt);
  const at::Tensor beta_t = at::where(dn_bt > 0, ds_bt / dn_bt, zero_t);
  const Scalar beta(beta_t, dn_beta.dynamic_dim(), dn_beta.intmd_dim());

  // Step 2: compute critical displacement jump
  const Scalar delta_n0 = _N / _K;
  const Scalar delta_s0 = _S / _K;

  // delta_init
  Scalar delta_init = delta_s0;
  {
    const Scalar delta_mixed = sqrt(delta_s0 * delta_s0 + (beta * delta_n0) * (beta * delta_n0));
    const Scalar di_mixed = delta_n0 * delta_s0 * sqrt(1.0 + beta * beta) / delta_mixed;

    const at::Tensor di_t = at::where(dn_bt > 0,
                                      static_cast<const at::Tensor &>(di_mixed),
                                      static_cast<const at::Tensor &>(delta_s0));
    delta_init = Scalar(di_t, dn_beta.dynamic_dim(), dn_beta.intmd_dim());
  }

  // Step 3: compute final displacement jump
  // delta_final
  Scalar delta_final = std::sqrt(2.0) * 2.0 * _GIIc / _S; // pure shear baseline
  {
    const at::Tensor beta2_t = static_cast<const at::Tensor &>(beta * beta);
    const at::Tensor one_t = at::ones_like(beta2_t);
    const at::Tensor r_t = beta2_t / (one_t + beta2_t); // r = beta^2/(1+beta^2)

    at::Tensor df_mixed_t = static_cast<const at::Tensor &>(delta_final);

    if (_criterion == "BK")
    {
      const Scalar r(r_t, beta.dynamic_dim(), beta.intmd_dim());
      const Scalar term = _GIc + (_GIIc - _GIc) * pow(r, _eta);
      const Scalar df = 2.0 / (_K * delta_init) * term;
      df_mixed_t = static_cast<const at::Tensor &>(df);
    }
    else if (_criterion == "POWER_LAW")
    {
      const Scalar beta2(beta2_t, beta.dynamic_dim(), beta.intmd_dim());
      const Scalar invGIc = 1.0 / _GIc;
      const Scalar invGIIc = 1.0 / _GIIc;

      const Scalar Gc_mixed = pow(invGIc, _eta) + pow(beta2 * invGIIc, _eta);

      const Scalar factor = (2.0 + 2.0 * beta2) / (_K * delta_init);
      const Scalar df = factor * pow(Gc_mixed, -1.0 / _eta);
      df_mixed_t = static_cast<const at::Tensor &>(df);
    }

    const at::Tensor df_t =
        at::where(dn_bt > 0, df_mixed_t, static_cast<const at::Tensor &>(delta_final));
    delta_final = Scalar(df_t, dn_beta.dynamic_dim(), dn_beta.intmd_dim());
  }

  // Step 4: compute effective displacement jump delta_m (using g_eff)
  // Effective displacement jump delta_m (using g_eff)
  const Scalar dn_eff = g_eff(0);
  const Scalar ds_eff = sqrt(g_eff(1) * g_eff(1) + g_eff(2) * g_eff(2));

  // dn_pos = H(dn)*dn (Macaulay regularization)
  const Scalar dn_pos = macauley_pos(dn_eff, _alpha);

  // delta_m = sqrt(ds^2 + dn_pos^2)
  const Scalar delta_m = sqrt(ds_eff * ds_eff + dn_pos * dn_pos);

  // Step 5: compute damage
  // Damage evolution (piecewise + irreversibility + viscosity)
  const at::Tensor dm_t = static_cast<const at::Tensor &>(delta_m);
  const at::Tensor di_t = static_cast<const at::Tensor &>(delta_init);
  const at::Tensor df_t = static_cast<const at::Tensor &>(delta_final);

  const at::Tensor d0 = at::zeros_like(dm_t);
  const at::Tensor d1 = at::ones_like(dm_t);

  // Guard only the physically required denominator and regularize dm_t near 0.
  const auto assert_nonzero = [](const at::Tensor & denom, const char * what)
  {
    const at::Tensor tol = at::full_like(denom, 1e-14);
    const bool bad_denom = at::any(at::abs(denom) <= tol).item<bool>();
    neml_assert(!bad_denom,
                "BilinearMixedModeTractionLaw: ",
                what,
                " is zero/too small. Check inputs and parameters.");
  };
  assert_nonzero(df_t - di_t, "delta_final - delta_init");
  assert_nonzero(dm_t, "effective displacement jump");

  const at::Tensor dmid = df_t * (dm_t - di_t) / dm_t / (df_t - di_t);

  at::Tensor dnew = at::where(dm_t < di_t, d0, at::where(dm_t > df_t, d1, dmid));

  // irreversibility
  const at::Tensor dold = static_cast<const at::Tensor &>(_damage_old());
  dnew = at::maximum(dnew, dold);

  // viscous regularization: (d + visc * d_old / dt) / (visc/dt + 1)
  const at::Tensor visc_t = static_cast<const at::Tensor &>(_visc);
  const at::Tensor dt_t = static_cast<const at::Tensor &>(_dt());
  const at::Tensor dt_safe = dt_t + at::full_like(dt_t, 1e-14);
  const at::Tensor denom = visc_t / dt_safe + 1.0;
  const at::Tensor dvisc = (dnew + visc_t * dold / dt_safe) / denom;

  const Scalar damage(dvisc, delta_m.dynamic_dim(), delta_m.intmd_dim());

  // Step 6: compute traction
  if (out)
  {
    _damage = damage;

    const Scalar dn = g(0);
    const Scalar zero = Scalar::zeros_like(dn);

    const Scalar dn_pos = maximum(dn, zero);
    const Scalar dn_neg = minimum(dn, zero);

    const Vec g_active = Vec::fill(dn_pos, g(1), g(2));
    const Vec g_inactive = Vec::fill(dn_neg, zero, zero);

    // Traction (active/inactive split; compression not degraded)
    const Scalar one_minus_d = 1.0 - damage;
    const Vec t_active = (one_minus_d * _K) * g_active;
    const Vec t_inactive = _K * g_inactive;
    _traction = t_active + t_inactive;

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

  // Step 7: compute derivatives
  if (dout_din)
  {
    // _traction.d(_jump);
  }
}

} // namespace neml2
