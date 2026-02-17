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
#include "neml2/tensors/R2.h"
#include "neml2/tensors/Scalar.h"
#include "neml2/tensors/Vec.h"
#include "neml2/tensors/functions/norm.h"
#include "neml2/tensors/functions/cos.h"
#include "neml2/tensors/functions/outer.h"
#include "neml2/tensors/functions/pow.h"
#include "neml2/tensors/functions/sin.h"
#include "neml2/tensors/functions/sqrt.h"
#include "neml2/tensors/functions/where.h"

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
  const auto zero = Scalar::zeros_like(x);
  const auto one = Scalar::ones_like(x);
  const auto mid = 0.5 * (1.0 + sin(M_PI * x / (2.0 * smoothing_length)));
  return where(x <= -smoothing_length, zero, where(x < smoothing_length, mid, one));
}

/**
 * if |x| < l: 0.25*pi/l*cos(pi*x/(2*l))
 * else: 0
 */
Scalar
BilinearMixedModeTractionLaw::regularizedHeavisideDerivative(const Scalar & x,
                                                             const Scalar & smoothing_length) const
{
  const auto zero = Scalar::zeros_like(x);
  const auto mid = 0.25 * (M_PI / smoothing_length) * cos(M_PI * x / (2.0 * smoothing_length));
  return where((x < smoothing_length) && (x > -smoothing_length), mid, zero);
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
  return where(a > b, a, b);
}

Scalar
BilinearMixedModeTractionLaw::minimum(const Scalar & a, const Scalar & b) const
{
  return where(a < b, a, b);
}

void
BilinearMixedModeTractionLaw::set_value(bool out, bool dout_din, bool /*d2out_din2*/)
{
  // _jump() and _jump_old() are in the local CZM frame:
  //   g(0) = normal component, g(1), g(2) = tangential components
  // This matches MOOSE:
  //   interface_jump_local = R^T * interface_jump_global
  // So we must NOT re-project using _normal().

  const auto assert_nonzero = [](const Scalar & denom, const char * what)
  {
    const at::Tensor dt = static_cast<const at::Tensor &>(denom);
    const at::Tensor tol = at::full_like(dt, 1e-14);
    const bool bad_denom = at::any(at::abs(dt) <= tol).item<bool>();
    neml_assert(!bad_denom,
                "BilinearMixedModeTractionLaw: ",
                what,
                " is zero/too small. Check inputs and parameters.");
  };

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
  const Scalar beta = where(dn_beta > 0, ds_beta / dn_beta, Scalar::zeros_like(dn_beta));

  // dbeta/ddelta is only defined in opening with non-zero shear jump.
  assert_nonzero(ds_beta, "shear jump for mode mixity");
  assert_nonzero(dn_beta, "normal jump for mode mixity");
  const Scalar one = Scalar::ones_like(dn_beta);
  const Vec dbeta_active = Vec::fill(-ds_beta / (dn_beta * dn_beta),
                                     delta(1) / (ds_beta * dn_beta),
                                     delta(2) / (ds_beta * dn_beta));
  const Vec dbeta_ddelta = where(dn_beta > 0, dbeta_active, Vec::zeros_like(delta));

  // Step 2: compute critical displacement jump
  const Scalar delta_n0 = _N / _K;
  const Scalar delta_s0 = _S / _K;

  // delta_init
  Scalar delta_init = delta_s0;
  Vec ddelta_init_ddelta = Vec::zeros_like(delta);
  {
    const Scalar delta_mixed = sqrt(delta_s0 * delta_s0 + (beta * delta_n0) * (beta * delta_n0));
    const Scalar di_mixed = delta_n0 * delta_s0 * sqrt(1.0 + beta * beta) / delta_mixed;

    delta_init = where(dn_beta > 0, di_mixed, delta_s0);

    // ddelta_init/ddelta = ddelta_init/dbeta * dbeta/ddelta
    if (!_lag_mode_mixity)
    {
      const Scalar ddelta_init_dbeta =
          delta_init * beta * (1 / (one + beta * beta) - pow(delta_init / delta_mixed, 2));
      ddelta_init_ddelta =
          where(dn_beta > 0, dbeta_ddelta * ddelta_init_dbeta, Vec::zeros_like(delta));
    }
  }

  // Step 3: compute final displacement jump
  // delta_final
  Scalar delta_final = std::sqrt(2.0) * 2.0 * _GIIc / _S; // pure shear baseline
  Vec ddelta_final_ddelta = Vec::zeros_like(delta);
  {
    const Scalar beta2 = beta * beta;
    const Scalar r = beta2 / (1.0 + beta2); // r = beta^2/(1+beta^2)
    Scalar delta_final_mixed = delta_final;

    if (_criterion == "BK")
    {
      const Scalar term = _GIc + (_GIIc - _GIc) * pow(r, _eta);
      delta_final_mixed = 2.0 / (_K * delta_init) * term;

      if (!_lag_mode_mixity)
      {
        // ddelta_final/ddelta = ddelta_final/ddelta_init * ddelta_init/ddelta +
        //                       ddelta_final/dbeta * dbeta/ddelta
        const Scalar ddelta_final_ddelta_init = -delta_final_mixed / delta_init;
        const Scalar ddelta_final_dbeta = 2.0 / (_K * delta_init) * (_GIIc - _GIc) * _eta *
                                          pow(r, _eta - 1.0) * 2.0 * beta *
                                          (1.0 - pow(beta / (1.0 + beta2), 2.0));
        const Vec ddelta_final_mixed_ddelta =
            ddelta_final_ddelta_init * ddelta_init_ddelta + ddelta_final_dbeta * dbeta_ddelta;
        ddelta_final_ddelta = where(dn_beta > 0, ddelta_final_mixed_ddelta, Vec::zeros_like(delta));
      }
    }
    else if (_criterion == "POWER_LAW")
    {
      const Scalar invGIc = 1.0 / _GIc;
      const Scalar invGIIc = 1.0 / _GIIc;
      const Scalar Gc_mixed = pow(invGIc, _eta) + pow(beta2 * invGIIc, _eta);

      const Scalar factor = (2.0 + 2.0 * beta2) / (_K * delta_init);
      delta_final_mixed = factor * pow(Gc_mixed, -1.0 / _eta);

      if (!_lag_mode_mixity)
      {
        // ddelta_final/ddelta = ddelta_final/ddelta_init * ddelta_init/ddelta +
        //                       ddelta_final/dbeta * dbeta/ddelta
        const Scalar ddelta_final_ddelta_init = -delta_final_mixed / delta_init;
        const Scalar ddelta_final_dbeta =
            delta_final_mixed * 2.0 * beta / (1.0 + beta2) -
            factor * pow(Gc_mixed, -1.0 / _eta - 1.0) *
                (pow(invGIc, _eta - 1.0) + pow(beta2 * invGIIc, _eta - 1.0) * 2.0 * beta * invGIIc);
        const Vec ddelta_final_mixed_ddelta =
            ddelta_final_ddelta_init * ddelta_init_ddelta + ddelta_final_dbeta * dbeta_ddelta;
        ddelta_final_ddelta = where(dn_beta > 0, ddelta_final_mixed_ddelta, Vec::zeros_like(delta));
      }
    }

    delta_final = where(dn_beta > 0, delta_final_mixed, delta_final);
  }

  // Step 4: compute effective displacement jump delta_m (using g_eff)
  // Effective displacement jump delta_m (using g_eff)
  const Scalar dn_eff = g_eff(0);
  const Scalar ds_eff = sqrt(g_eff(1) * g_eff(1) + g_eff(2) * g_eff(2));

  // dn_pos = H(dn)*dn (Macaulay regularization)
  const Scalar dn_eff_pos = macauley_pos(dn_eff, _alpha);

  // delta_m = sqrt(ds^2 + dn_pos^2)
  const Scalar delta_m = sqrt(ds_eff * ds_eff + dn_eff_pos * dn_eff_pos);

  // d(delta_m)/d(delta)
  Vec ddelta_m_ddelta = Vec::zeros_like(g);
  if (!_lag_disp_jump)
  {
    const auto hdiff_smoothing = Scalar::full_like(dn_eff, 1e-6);
    const Scalar ddn_pos_ddn = regularizedHeavisideDerivative(dn_eff, hdiff_smoothing) * dn_eff +
                               regularizedHeaviside(dn_eff, _alpha);
    const Vec ddelta_m_num = Vec::fill(dn_eff_pos * ddn_pos_ddn, g_eff(1), g_eff(2));
    ddelta_m_ddelta = where(delta_m > 1e-14, ddelta_m_num / delta_m, Vec::zeros_like(g));
  }

  // Step 5: compute damage
  // Damage evolution (piecewise + irreversibility + viscosity), aligned with MOOSE.
  const Scalar d0 = Scalar::zeros_like(delta_m);
  const Scalar d1 = Scalar::ones_like(delta_m);

  // Guard only the physically required denominator.
  assert_nonzero(delta_final - delta_init, "delta_final - delta_init");

  const Scalar dmid = delta_final * (delta_m - delta_init) / delta_m / (delta_final - delta_init);
  // if, else if, else in one expression
  const Scalar d_trial = where(delta_m < delta_init, d0, where(delta_m > delta_final, d1, dmid));

  // dd/ddelta in the softening branch, 0 outside and under irreversibility.
  const Scalar denom_mid = delta_m * (delta_final - delta_init);
  const Vec numer_1 = ddelta_final_ddelta * (delta_m - delta_init) +
                      delta_final * (ddelta_m_ddelta - ddelta_init_ddelta);
  const Vec numer_2 = ddelta_m_ddelta * (delta_final - delta_init) +
                      delta_m * (ddelta_final_ddelta - ddelta_init_ddelta);
  const Vec dd_mid =
      numer_1 / denom_mid - delta_final * (delta_m - delta_init) * numer_2 / pow(denom_mid, 2.0);

  // zero derivative outside the softening branch
  const Scalar d_old = _damage_old();
  const Scalar irrev = d_trial < d_old;
  const Scalar d_irrev = where(irrev, d_old, d_trial);
  Vec dd_ddelta =
      where(delta_m < delta_init || delta_m > delta_final || irrev, Vec::zeros_like(g), dd_mid);

  // Viscous regularization.
  assert_nonzero(_dt(), "time step size for viscous regularization");
  const Scalar visc_denom = _visc / _dt() + 1.0;
  const Scalar damage = (d_irrev + _visc * d_old / _dt()) / visc_denom;
  dd_ddelta = dd_ddelta / visc_denom;

  // Shared split terms used by output and derivatives.
  const Scalar dn = g(0);
  const Scalar zero = Scalar::zeros_like(dn);
  const Scalar dn_pos = maximum(dn, zero);
  const Scalar dn_neg = minimum(dn, zero);
  const Vec g_active = Vec::fill(dn_pos, g(1), g(2));
  const Vec g_inactive = Vec::fill(dn_neg, zero, zero);

  // Step 6: compute traction
  if (out)
  {
    _damage = damage;

    // Traction (active/inactive split)
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
    if (_jump.is_dependent())
    {
      _damage.d(_jump) = dd_ddelta;

      // d(traction)/d(delta): split contribution + damage contribution.
      const Scalar one = Scalar::ones_like(dn);

      const Scalar dactive_n_ddn = where(dn > 0, one, zero);
      const Scalar dinactive_n_ddn = where(dn < 0, one, zero);
      const R2 ddelta_active_ddelta = R2::fill(dactive_n_ddn, one, one);
      const R2 ddelta_inactive_ddelta = R2::fill(dinactive_n_ddn, zero, zero);

      R2 dtraction_ddelta =
          (1.0 - damage) * _K * ddelta_active_ddelta + _K * ddelta_inactive_ddelta;

      dtraction_ddelta = dtraction_ddelta - _K * outer(g_active, dd_ddelta);

      _traction.d(_jump) = dtraction_ddelta;
    }
  }
}

} // namespace neml2
