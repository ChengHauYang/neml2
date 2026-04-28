# NEML2 — A Code-to-Physics Reverse Engineering

*A pedagogical reconstruction of the continuum theories, constitutive equations, and numerical formulations that NEML2 actually implements — extracted by reading the C++ source, not the comments.*

---

## Executive Summary

NEML2 is a **multi-physics constitutive modeling library** for the engineering response of materials under combined mechanical, thermal, chemical, and fracture-driving loads. Read at the level of equations rather than classes, the codebase implements the following physical systems:

| Domain | What is solved | Where it lives |
|---|---|---|
| **Continuum solid mechanics** | $\boldsymbol{\sigma}(\boldsymbol{\varepsilon}, \boldsymbol{q}, T, t)$ — constitutive response of the Cauchy / Mandel stress under elasto-(visco-)plasticity, anisotropic elasticity, hardening, and ductile damage | `models/solid_mechanics/` |
| **Crystal plasticity** | Slip-system kinematics on FCC / BCC / HCP crystals; orientation evolution under multiplicative $\boldsymbol{F}=\boldsymbol{F}_e \boldsymbol{F}_p$ | `models/crystallography/`, `models/solid_mechanics/crystal_plasticity/` |
| **Brittle phase-field fracture** | AT1 / AT2 regularized variational fracture with energy decomposition (no-split / vol-dev / spectral) and KKT-irreversibility | `models/phase_field_fracture/` |
| **Unsaturated porous flow & poromechanics** | Porosity-dependent permeability (Kozeny–Carman, power, exponential), van Genuchten / Brooks–Corey retention curves, Biot-type coupling stress | `models/porous_flow/` |
| **Chemo-mechanical reaction kinetics** | Solid-state nucleation/growth (Avrami–Erofeev), shrinking-core (contracting geometry), diffusion-limited reactions, geometry feedback | `models/chemical_reactions/` |
| **1-D finite volume transport** | First-order upwind advective fluxes with Dirichlet BCs and gradient reconstruction | `models/finite_volume/` |
| **Generic time integration & rate algebra** | Backward / forward Euler, finite-difference rates, SR2 invariants, linear combinations | `models/common/`, `models/ImplicitUpdate.cxx` |

The unifying idea is that every physical relation — whether constitutive, kinematic, kinetic, balance, or boundary — is expressed as a **node** of a directed acyclic graph that consumes named tensor quantities and produces named tensor quantities. Both the residual *and the analytic Jacobian* of every equation are written explicitly in the code (no opaque AD black boxes for the leaves), which makes the framework suitable for fully implicit, batched, GPU-residential integration.

> **Key Physical Insight #0 — Physics fidelity through differentiability.** The code preserves analytic Jacobians for every constitutive law it ships. Where regularization is unavoidable (e.g. derivative of a tensor norm at the origin), the regularization parameter is principled (machine $\varepsilon$, Macaulay smoothing) rather than ad-hoc, so consistency between the residual and its tangent is exact to floating point.

---

## 1. The Physics Tree (Mind Map)

```
                              ┌──────────────────────────┐
                              │     User input file (.i) │
                              │   Declares ComposedModel │
                              └────────────┬─────────────┘
                                           │
                       ┌───────────────────┴──────────────────┐
                       │       ComposedModel = DAG of         │
                       │       physical/operational nodes     │
                       └───┬────────────────────────────────┬─┘
                           │                                │
        ┌──────────────────┴──────────┐         ┌───────────┴──────────────────┐
        │  PHYSICAL NODES             │         │  OPERATIONAL NODES           │
        │  (close a constitutive eqn) │         │  (algebra, time, solvers)    │
        └───┬─────────────────────────┘         └───┬──────────────────────────┘
            │                                       │
   ┌────────┼──────────────────────────────┐    ┌───┼──────────────────────────┐
   │        │                              │    │   │                          │
   ▼        ▼                              ▼    ▼   ▼                          ▼
SOLID  CRYSTAL                      PHASE-FIELD  COMMON algebra            DRIVERS
MECH   PLASTICITY                   FRACTURE     • SR2Invariant            • TransientDriver
       │                            │            • LinearCombination       • prescribed forces
       │                            │            • VariableRate (rates)    • predictor/corrector
       │                            │
       │                            ├── DegradationFunction g(d)
       │                            ├── CrackGeometricFunction α(d)  (AT1/AT2)
       │                            ├── StrainEnergyDensity ψ⁺/ψ⁻
       │                            └── KKT complementarity (Fischer–Burmeister)
       │
       ├── Elasticity:    σ = ℂ(K,G,θ_orient) : εᵉ        (iso/cubic/Green-Lagrange)
       ├── Yield:         f(σ̄, σ_y, h, X, φ)            (J2, GTN)
       ├── Flow rule:     ε̇ᵖ = γ̇ Nₘ                     (associative, Perzyna VP)
       ├── Hardening:     ḣ, Ẋ                            (Voce, FA, Chaboche)
       ├── Damage:        φ̇ = (1−φ) tr(ε̇ᵖ)              (Gurson cavitation)
       └── Crystal plasticity:
            • slip-system geometry s^α ⊗ m^α
            • τ^α = σ : Mˢʸᵐ
            • γ̇^α = γ̇₀ |τ^α/g^α|^n sign(τ^α)
            • ġ = θ₀(1 − g/τ_f) Σ|γ̇^β|
            • Lᵖ = Σ γ̇^α (s^α ⊗ m^α), rotated by Q
            • orientation rate: Ω = w − wᵖ + [dᵖ, ε]

   ┌─────────────┬──────────────────────────────┐
   ▼             ▼                              ▼
POROUS FLOW   CHEMICAL REACTIONS         FINITE VOLUME (1-D)
• EffectiveSat• Avrami f = k(1−a)(−ln(1−a))ⁿ   • upwind flux F_{i+½} = v⁺uᵢ + v⁻uᵢ₊₁
• k(φ) laws    • contracting f = k(1−a)ⁿ        • central gradient
• Pc(Sₑ)       • diffusion-limited shrinking    • Dirichlet BC append
• AdvectiveStr.  core (cylindrical channels)
```

Two things to notice from the mind map:

1. **No PDE solver lives inside `models/`.** The library produces *constitutive closures and discretized fluxes*; momentum / energy / mass balances at the continuum level are the host code's responsibility (NEML2 is meant to plug into MOOSE, libMesh, or a Python driver). The `finite_volume/` module is an exception — a deliberately minimal 1-D demonstrator.
2. **The "operational" nodes** (BackwardEuler, ForwardEuler, ImplicitUpdate, VariableRate, SR2Invariant) are *mathematically generic*. They are what turns the bag of pointwise constitutive equations above into a *time-integrable*, *implicitly solvable* system.

---

## 2. Mathematical Foundation

### 2.1 Kinematics and stress measures

NEML2 supports both **small-strain additive** and **large-deformation multiplicative** decompositions.

**Multiplicative split** (large deformation):
$$\boldsymbol{F} = \boldsymbol{F}_e \boldsymbol{F}_p, \qquad \boldsymbol{L} = \dot{\boldsymbol{F}}\boldsymbol{F}^{-1} = \boldsymbol{L}_e + \boldsymbol{L}_p, \qquad \boldsymbol{D} = \mathrm{sym}(\boldsymbol{L}),\;\; \boldsymbol{W} = \mathrm{skew}(\boldsymbol{L})$$

**Green–Lagrange strain** (`elasticity/GreenLagrangeStrain.cxx:39–61`):
$$\boldsymbol{E} = \tfrac{1}{2}(\boldsymbol{C} - \boldsymbol{I}) = \tfrac{1}{2}(\boldsymbol{F}^\top \boldsymbol{F} - \boldsymbol{I})$$

**Stress measures used in the code**

| Symbol | Meaning | Where used |
|---|---|---|
| $\boldsymbol{\sigma}$ | Cauchy stress (true stress) | Outputs from elasticity; inputs to `MandelStress` |
| $\boldsymbol{M}$ | Mandel stress, work-conjugate to $\boldsymbol{L}_p$ in the intermediate config | All flow rules, yield, hardening |
| $\boldsymbol{P}$ | First Piola–Kirchhoff | `AdvectiveStress` (poromechanics) |
| $\boldsymbol{S}$ | Second Piola–Kirchhoff | Hyperelastic models with Green–Lagrange strain |

`IsotropicMandelStress.cxx:49` collapses $\boldsymbol{M} \equiv \boldsymbol{\sigma}$ for isotropic small-strain problems — a simplification that is *exact* in the small-strain regime and a *very good approximation* under finite rotations when the elastic stretch is small (typical of metals).

> **Key Physical Insight #1 — Mandel stress is the pivot.** Every plasticity submodel reads $\boldsymbol{M}$, never $\boldsymbol{\sigma}$. This is the single most important convention in the codebase: it is what makes flow rules frame-objective without any explicit Jaumann / Green–Naghdi rate machinery.

### 2.2 Elasticity

**Linear isotropic** (`elasticity/LinearIsotropicElasticity.cxx:38–62`):
$$\boldsymbol{\sigma} = K\,\mathrm{tr}(\boldsymbol{\varepsilon}_e)\,\boldsymbol{I} + 2G\,\mathrm{dev}(\boldsymbol{\varepsilon}_e), \qquad \mathbb{C} = 3K\,\mathbb{I}_v + 2G\,\mathbb{I}_d$$

Volumetric/deviatoric projection is encoded as the *fourth-order isotropic projectors* $\mathbb{I}_v = \tfrac{1}{3}\boldsymbol{I}\otimes\boldsymbol{I}$ and $\mathbb{I}_d = \mathbb{I}^{\mathrm{sym}} - \mathbb{I}_v$ (`elasticity/IsotropicElasticityTensor.cxx:60`). This is *not* an aesthetic choice — it lets the inverse map (compliance form) be written without inverting any matrix.

**Cubic anisotropic** (`elasticity/CubicElasticityTensor.cxx:62`) parameterizes the stiffness through the three independent cubic constants $C_{11}, C_{12}, C_{44}$ via three symmetry basis tensors:
$$\mathbb{C} = C_1\,\mathbb{I}_1 + C_2\,\mathbb{I}_2 + C_3\,\mathbb{I}_3$$

For crystal plasticity the same $\mathbb{C}$ is rotated by the orientation $\boldsymbol{Q}$:
$$\mathbb{C}_{\text{sample}} = (\boldsymbol{Q}\otimes\boldsymbol{Q}\otimes\boldsymbol{Q}\otimes\boldsymbol{Q}):\mathbb{C}_{\text{crystal}}$$

### 2.3 Yield criteria

**Von Mises (J2)** (`solid_mechanics/YieldFunction.cxx:68–71`):
$$f = \sqrt{\tfrac{2}{3}}\,(\bar{\sigma} - \sigma_y - h), \qquad \bar{\sigma} = \sqrt{\tfrac{3}{2}}\,\|\mathrm{dev}(\boldsymbol{M}-\boldsymbol{X})\|$$

The pre-factor $\sqrt{2/3}$ is what makes the *plastic multiplier* $\dot\gamma$ in the flow rule equal the equivalent plastic strain rate $\dot{\bar\varepsilon}^p$ — a normalization convention that propagates into every hardening law downstream.

**Gurson–Tvergaard–Needleman (GTN) ductile damage** (`solid_mechanics/GTNYieldFunction.cxx:89–91`):
$$f = \left(\frac{\bar\sigma}{\sigma_y + h}\right)^2 + 2 q_1 \phi \cosh\!\left(\frac{q_2}{2}\,\frac{3\sigma_h - \sigma_s}{\sigma_y + h}\right) - (q_3 \phi^2 + 1)$$

with $\phi$ the void volume fraction, $\sigma_h$ the hydrostatic stress, $\sigma_s$ a sintering back-pressure, and $(q_1, q_2, q_3)$ the Tvergaard fit parameters that capture void interaction. The second term *couples hydrostatic stress to plastic flow* — the physics that makes spallation and ductile rupture possible.

### 2.4 Flow rules

**Associative flow** (`solid_mechanics/AssociativePlasticFlow.cxx:62–67`):
$$\dot{\boldsymbol{\varepsilon}}^p = \dot\gamma\,\boldsymbol{N}_M, \qquad \boldsymbol{N}_M = \frac{\partial f}{\partial \boldsymbol{M}}$$

For J2 specifically (`AssociativeJ2FlowDirection.cxx:60–62`):
$$\boldsymbol{N}_M = \tfrac{3}{2}\frac{\mathrm{dev}(\boldsymbol{M})}{\bar\sigma}$$

The denominator $\bar\sigma$ is regularized by `machine_precision()` so the derivative is well-defined at the elastic origin without changing the flow direction anywhere it matters.

**Perzyna viscoplasticity** (`solid_mechanics/PerzynaPlasticFlowRate.cxx:62–66`):
$$\dot\gamma = \left\langle\frac{f}{\eta}\right\rangle^n$$

where $\langle\cdot\rangle = \max(\cdot, 0)$ is the Macaulay bracket (encoded as a Heaviside multiplication in the code). This is the rate-dependent generalization of the von Mises criterion: "yield" becomes a smooth onset rather than a hard switch, and rate sensitivity $n$ controls the steepness of the onset.

> **Key Physical Insight #2 — The Macaulay bracket is the entire physics of rate-dependent plasticity.** It encodes both unilateral activation ("no plastic flow when below yield") and the smoothness needed for an implicit Newton solve. `Heaviside(f)` in the code is doing genuine constitutive work, not numerical filtering.

### 2.5 Hardening

**Isotropic, linear** (`LinearIsotropicHardening.cxx:37–57`): $h = K\bar\varepsilon^p$.

**Isotropic, Voce (saturation)** (`VoceIsotropicHardening.cxx:38–61`):
$$h(\bar\varepsilon^p) = R\bigl[1 - \exp(-d\,\bar\varepsilon^p)\bigr], \quad \frac{\partial h}{\partial \bar\varepsilon^p} = R d\,\exp(-d\,\bar\varepsilon^p)$$

A rate-form variant (`SlopeSaturationVoceIsotropicHardening.cxx:67`):
$$\dot h = \theta_0\!\left(1 - \frac{h}{R}\right)\dot\gamma$$

**Kinematic, linear** (`LinearKinematicHardening.cxx:39–59`): $\boldsymbol{X} = H\,\boldsymbol{K}_p$.

**Frederick–Armstrong (dynamic recovery)** (`FredrickArmstrongPlasticHardening.cxx:40–74`):
$$\dot{\boldsymbol{X}} = \left(\tfrac{2}{3}C\,\boldsymbol{N}_M - g\boldsymbol{X}\right)\dot\gamma$$

The first term is the *direct hardening* (back-stress accumulation along the flow direction); the second is *dynamic recovery* — back-stress that bleeds away in proportion to its own magnitude, which is what produces realistic Bauschinger response and ratcheting saturation.

**Chaboche (with static recovery)** (`ChabochePlasticHardening.cxx:69–74`):
$$\dot{\boldsymbol{X}} = \left(\tfrac{2}{3}C\,\boldsymbol{N}_M - g\boldsymbol{X}\right)\dot\gamma \;-\; A\,\|\boldsymbol{X}\|^{a-1}\boldsymbol{X}$$

The third (subtractive) term is *static recovery* — time-driven, not strain-driven — that captures creep-induced relaxation of the back-stress. The combination $C, g, A, a$ together specify the cyclic stress–strain hysteresis loop shape, mean-stress relaxation, and ratcheting rate — all from one rate ODE per back-stress.

> **Key Physical Insight #3 — A multi-back-stress Chaboche model is just a `ComposedModel` of single back-stress nodes.** Each $\boldsymbol{X}_i$ with its own $(C_i, g_i, A_i, a_i)$ is a separate node summed by `LinearCombination` to produce $\boldsymbol{X} = \sum_i \boldsymbol{X}_i$. The DAG composition mechanism *is* the multi-mechanism modeling capability.

### 2.6 Damage / ductile rupture

**Gurson cavitation kinematics** (`GursonCavitation.cxx:39–63`):
$$\dot\phi = (1-\phi)\,\mathrm{tr}(\dot{\boldsymbol{\varepsilon}}^p)$$

This is *exact* mass conservation for an incompressible matrix material containing voids of volume fraction $\phi$: matrix volume change is forbidden, so all volumetric plastic strain has to go into void growth (or shrinkage when $\mathrm{tr}\dot{\boldsymbol{\varepsilon}}^p<0$). Coupled with the GTN yield surface above, this closes the loop: hydrostatic tension → dilatational plastic flow → void growth → softer yield → more flow.

### 2.7 Thermally-activated viscoplasticity (Kocks–Mecking)

**Flow viscosity** (`KocksMeckingFlowViscosity.cxx:38–82`):
$$\eta(T) = e^{B}\,\mu(T)\,\dot\varepsilon_0^{\,-kTA / (\mu b^3)}$$

**Yield stress** (`KocksMeckingYieldStress.cxx:38–63`):
$$\sigma_y(T) = e^{C}\,\mu(T)$$

The exponents $A,B,C$ are slope-intercept fits to a deformation-mechanism map: $A$ encodes the activation enthalpy of dislocation glide normalized by $\mu b^3$; $B,C$ are the curve intercepts. Temperature enters through $\mu(T)$ (modulus softening) *and* through the exponent — this is what reproduces the *power-law breakdown* and the *deformation regime change* between athermal and thermally-activated flow on the same constitutive law.

### 2.8 Crystal plasticity

NEML2 implements **rate-dependent crystal plasticity** with full multiplicative kinematics.

**Slip system geometry** (`crystallography/CrystalGeometry.cxx:112–116, 199–263`):

For each slip system $\alpha$ with slip direction $\boldsymbol{s}^\alpha$ and plane normal $\boldsymbol{m}^\alpha$ (unit vectors, $\boldsymbol{s}^\alpha\cdot\boldsymbol{m}^\alpha=0$):
$$\boldsymbol{P}^\alpha = \boldsymbol{s}^\alpha \otimes \boldsymbol{m}^\alpha, \quad \boldsymbol{M}^\alpha = \mathrm{sym}\,\boldsymbol{P}^\alpha, \quad \boldsymbol{W}^\alpha = \mathrm{skew}\,\boldsymbol{P}^\alpha$$

Crystal symmetry operators (e.g. the 24 elements of point group `m3m` for cubic crystals, `crystallography.h:36–72`) are applied to expand a single Miller-index family into all crystallographically equivalent systems.

**Resolved shear stress** (`solid_mechanics/crystal_plasticity/ResolvedShear.cxx:74–80`):
$$\tau^\alpha = \boldsymbol{\sigma} : \mathrm{rot}_{\boldsymbol{Q}}(\boldsymbol{M}^\alpha) = \boldsymbol{\sigma} : (\boldsymbol{Q}\boldsymbol{M}^\alpha\boldsymbol{Q}^\top)$$

Note: the *Schmid tensor* is rotated to the sample frame, not the stress to the crystal frame — algebraically equivalent but more efficient when the stress is itself an output of an upstream node.

**Power-law slip rate** (`PowerLawSlipRule.cxx:65`):
$$\dot\gamma^\alpha = \dot\gamma_0\left|\frac{\tau^\alpha}{g^\alpha}\right|^{n-1}\frac{\tau^\alpha}{g^\alpha}$$

This is the Hutchinson form (a.k.a. Norton–Hoff in dislocation kinetics): $n$ is the inverse strain-rate sensitivity $1/m$, $g^\alpha$ is the slip resistance.

**Single-slip Voce hardening** (`VoceSingleSlipHardeningRule.cxx:61`):
$$\dot g = \theta_0\!\left(1 - \frac{g}{\tau_f}\right)\sum_{\beta} |\dot\gamma^\beta|$$

All slip systems share the same $g$ — i.e. the *latent-hardening matrix* is implicitly $h_{\alpha\beta}=1$. (Anisotropic latent hardening, the Asaro–Needleman matrix, is a natural extension that the framework is structured to accept.)

**Plastic deformation rate and spin** (`PlasticDeformationRate.cxx:75`, `PlasticVorticity.cxx:75`):
$$\boldsymbol{D}^p_{\text{crystal}} = \sum_\alpha \dot\gamma^\alpha\,\boldsymbol{M}^\alpha, \qquad \boldsymbol{W}^p_{\text{crystal}} = \sum_\alpha \dot\gamma^\alpha\,\boldsymbol{W}^\alpha$$

then rotated to the sample frame via $\boldsymbol{Q}\,(\cdot)\,\boldsymbol{Q}^\top$.

**Elastic strain rate** (`ElasticStrainRate.cxx:96`):
$$\dot{\boldsymbol{\varepsilon}}^e = \boldsymbol{D} - \boldsymbol{D}^p + (\boldsymbol{W}\,\boldsymbol{\varepsilon}^e - \boldsymbol{\varepsilon}^e\,\boldsymbol{W})$$

The bracket term is the **Jaumann co-rotational correction** — what makes the elastic strain rate objective in a finitely rotating sample frame.

**Orientation rate** (`OrientationRate.cxx:97`):
$$\boldsymbol{\Omega}^e = \dot{\boldsymbol{Q}}\boldsymbol{Q}^\top = \boldsymbol{W} - \boldsymbol{W}^p + (\boldsymbol{D}^p\,\boldsymbol{\varepsilon}^e - \boldsymbol{\varepsilon}^e\,\boldsymbol{D}^p)$$

The first two terms are the textbook *crystal lattice spin* = total spin minus plastic spin. The third (commutator) term is the elastic-distortion-induced lattice rotation that becomes important at large strains.

> **Key Physical Insight #4 — Texture evolution emerges, it isn't programmed.** Once $\dot{\boldsymbol{Q}}$ is integrated by `BackwardEulerTimeIntegration`, repeated application over a polycrystal sample gives the full Taylor / VPSC-style texture without any extra orientation-update machinery.

### 2.9 Phase-field fracture

NEML2 implements the variational fracture framework of Bourdin–Francfort–Marigo with the Miehe energy decomposition.

**Total functional** (assembled in test input file `tests/regression/phase_field_fracture/elastic_brittle_fracture/small_deformation.i`):
$$\Pi[\boldsymbol{u}, d] = \int_\Omega \Bigl[g(d)\,\psi^+(\boldsymbol{\varepsilon}) + \psi^-(\boldsymbol{\varepsilon})\Bigr]\,dV \;+\; \frac{G_c}{c_0}\int_\Omega\Bigl[\frac{\alpha(d)}{\ell_c} + \ell_c\,|\nabla d|^2\Bigr]\,dV$$

**Degradation function** — power form (`PowerDegradationFunction.cxx:59`):
$$g(d) = (1-d)^p\,(1-\eta) + \eta$$

— rational form (`RationalDegradationFunction.cxx:66–69`):
$$g(d) = \frac{(1-d)^p}{(1-d)^p + b_1 d\,(1 + b_2 d + b_2 b_3 d^2)}\,(1-\eta) + \eta$$

The residual stiffness $\eta$ is what keeps the tangent matrix non-singular at fully broken material — a small but indispensable numerical regularization.

**Crack geometric function** — AT1 (`CrackGeometricFunctionAT1.cxx:51`): $\alpha(d) = d$. AT2 (`CrackGeometricFunctionAT2.cxx:52`): $\alpha(d) = d^2$.

The choice has *physical content*: AT1 has an elastic phase before damage initiates (no damage until a critical stress); AT2 begins damaging from the elastic origin. This affects whether the model can reproduce a true *crack initiation toughness* or only *propagation toughness*.

**Energy decomposition** (`LinearIsotropicStrainEnergyDensity.cxx:89–184`) — three options:

*No decomposition*:  $\psi^+ = \tfrac{1}{2}K(\mathrm{tr}\,\boldsymbol{\varepsilon})^2 + G\,\mathrm{dev}\boldsymbol{\varepsilon}:\mathrm{dev}\boldsymbol{\varepsilon}$, $\psi^- = 0$.

*Volumetric–deviatoric* (Amor):
$$\psi^+ = \tfrac{1}{2}K\langle\mathrm{tr}\,\boldsymbol{\varepsilon}\rangle_+^2 + G\,\mathrm{dev}\boldsymbol{\varepsilon}:\mathrm{dev}\boldsymbol{\varepsilon}, \quad \psi^- = \tfrac{1}{2}K\langle-\mathrm{tr}\,\boldsymbol{\varepsilon}\rangle_+^2$$

*Spectral* (Miehe):
$$\psi^\pm = \tfrac{1}{2}\lambda\,\langle\pm\mathrm{tr}\,\boldsymbol{\varepsilon}\rangle_+^2 + G\,\boldsymbol{\varepsilon}^\pm:\boldsymbol{\varepsilon}^\pm$$

with $\boldsymbol{\varepsilon}^\pm = \sum_i \langle\pm\lambda_i\rangle_+\,\boldsymbol{n}_i\otimes\boldsymbol{n}_i$ from the eigendecomposition (`eigh`, `ieigh`).

**Coupled stress** (`Normality` model on `psi`):
$$\boldsymbol{\sigma} = \frac{\partial \psi}{\partial \boldsymbol{\varepsilon}} = g(d)\,\frac{\partial \psi^+}{\partial \boldsymbol{\varepsilon}} + \frac{\partial \psi^-}{\partial \boldsymbol{\varepsilon}}$$

**Phase-field equation with KKT irreversibility:**
$$F := \frac{\partial \psi}{\partial d} \ge 0, \qquad \dot d \ge 0, \qquad F\,\dot d = 0$$

NEML2 enforces this by the *Fischer–Burmeister* smooth complementarity function ($\phi(a,b)=\sqrt{a^2+b^2}-(a+b)=0 \Leftrightarrow a,b\ge 0,\ ab=0$), which makes the KKT conditions tractable to a Newton solve without active-set switching.

> **Key Physical Insight #5 — Three orthogonal modeling choices.** $g(d)$ controls *softening shape*; $\alpha(d)$ controls *initiation behavior*; the energy split controls *unilateral contact*. Because each is a separate node, you can mix and match (Power+AT2+Spectral, Rational+AT1+VolDev, etc.) without rewriting any equation.

### 2.10 Porous flow

The porous-flow module supplies the *constitutive closures* needed by an external mass-balance solver — it does not advance the mass-balance PDE itself. The closures are:

**Effective saturation** (`EffectiveSaturation.cxx:65`):
$$S_e = \frac{\phi/\phi_{\max} - S_r}{1 - S_r}$$

**Porosity-permeability laws** (the $k\,k_{r}/\mu$ block in Darcy's law):

| Law | Equation | File |
|---|---|---|
| Power-law | $K = K_0 (\varphi/\varphi_0)^p$ | `PowerLawPermeability.cxx:58–62` |
| Kozeny–Carman | $K = K_0\,\frac{\varphi^n (1-\varphi_0)^m}{\varphi_0^m (1-\varphi)^n}$ | `KozenyCarmanPermeability.cxx:62–66` |
| Exponential | $K = K_0 \exp[-a(\varphi-\varphi_0)]$ | `ExponentialLawPermeability.cxx:59–62` |

Kozeny–Carman is the physically motivated choice (derived from capillary-bundle theory), the others are empirical fits for compacted soils and creep-densified powders.

**Capillary pressure curves** (`CapillaryPressure.cxx:76–94` for the log-extension trick):

van Genuchten (`VanGenuchtenCapillaryPressure.cxx:65–69`):
$$P_c = a\bigl(S_e^{-1/m} - 1\bigr)^{1-m}$$

Brooks–Corey (`BrooksCoreyCapillaryPressure.cxx:62–64`):
$$P_c = P_t\,S_e^{-1/p}$$

Both are blended below a transition saturation $S_p$ to a logarithmic extension to keep $P_c$ finite and the Newton tangent well-conditioned at the dry limit.

**Poromechanical coupling** (`AdvectiveStress.cxx:79–95`):
$$p_s = -\frac{c}{3J}\,J_s^{-5/3}\,J_t^{-2/3}\,\boldsymbol{P} : \boldsymbol{F}$$

This is the contribution of pore-pressure-driven volume change to the effective Cauchy stress, separating swelling Jacobian $J_s$, thermal Jacobian $J_t$, and total Jacobian $J = \det\boldsymbol{F}$ — the Biot decomposition done multiplicatively rather than additively, appropriate for finite-strain poromechanics.

### 2.11 Chemical reaction kinetics

The chemistry module implements **conversion-degree kinetics** — ODEs in the dimensionless conversion $a \in [0,1]$ — rather than detailed elementary reactions.

**Avrami–Erofeev nucleation-growth** (`AvramiErofeevNucleation.cxx:59`):
$$\dot a = k(T)\,(1-a)\,(-\ln(1-a))^n$$

The standard form for solid-state phase transformations (e.g. recrystallization, devitrification): the $(-\ln(1-a))^n$ kernel encodes the Johnson–Mehl–Avrami sigmoid, which is what gives the characteristic S-shaped conversion curve.

**Contracting geometry** (`ContractingGeometry.cxx:62`):
$$\dot a = k(T)\,(1-a)^n$$

Power-law decay; the *contracting cylinder* limit gives $n=1/2$ and the *contracting sphere* gives $n=2/3$.

**Diffusion-limited shrinking-core reaction** (`DiffusionLimitedReaction.cxx:72–78`):
$$r = \frac{2 D R_l R_s}{\omega}\,\frac{r_o}{r_o - r_i + \delta}$$

with $r_i, r_o$ the inner/outer radii of the reaction product layer. As the layer thickens ($r_i \to r_o$ small), the reaction rate $\to 0$ — diffusion-controlled passivation (think of metal oxidation that slows once a tight oxide layer forms).

**Geometric feedback** (`CylindricalChannelGeometry.cxx:63–76`):
$$r_i = \sqrt{1 - \phi_s - \phi_p}, \qquad r_o = \sqrt{1 - \phi_s}$$

closes the loop between volume fractions of solid reactant $\phi_s$ and solid product $\phi_p$ and the geometric quantities $r_i, r_o$ that drive the diffusion-limited rate.

**Volume balance** (`EffectiveVolume.cxx:87–95`):
$$V = \frac{M}{1 - \phi_o}\sum_i \frac{\omega_i}{\rho_i}$$

assembles the total volume from species mass fractions $\omega_i$ and densities $\rho_i$ accounting for an open (interconnected) porosity $\phi_o$.

Temperature enters through the rate constants $k(T)$, which are declared as `allow_nonlinear=true` — meaning $k$ itself can be the output of an upstream Arrhenius node like
$$k(T) = A_0\,\exp(-E_a/RT).$$

The module is therefore *Arrhenius-ready* without containing an Arrhenius law as a hardcoded part of any kinetic equation.

### 2.12 Finite-volume transport

`finite_volume/` discretizes 1-D scalar advection. The basic operators are:

**Cell-centered gradient** (`FiniteVolumeGradient.cxx:70–73`):
$$(\nabla u)_i \approx -\beta\,\frac{u_{i+1} - u_i}{\Delta x_{i+1/2}}$$

**Edge interpolation** (`LinearlyInterpolateToCellEdges.cxx:99–103`):
$$u_{i+1/2} = w_L\,u_i + w_R\,u_{i+1}, \qquad w_L + w_R = 1$$

**Upwind advective flux** (`FiniteVolumeUpwindedAdvectiveFlux.cxx:65–72`):
$$F_{i+1/2} = v^+\,u_i + v^-\,u_{i+1}, \quad v^\pm = \tfrac{1}{2}(v \pm |v|)$$

This is the canonical first-order monotone upwind scheme. It is *diffusive* (numerical diffusivity $\sim \tfrac{1}{2}|v|\Delta x$) but *unconditionally non-oscillatory* — the right starting point for nonlinear advection-reaction problems where stability matters more than sharpness.

The semi-discrete cell update is then
$$\frac{du_i}{dt} = -\frac{F_{i+1/2} - F_{i-1/2}}{\Delta x_i} + S_i$$
with $S_i$ supplied by an upstream reaction node. Time integration is delegated to `BackwardEulerTimeIntegration` (Section 3 below).

---

## 3. Implementation Insights

### 3.1 The implicit-update pattern

The pivot from "constitutive equations" to "solvable system" is `ImplicitUpdate.cxx`, which packages a model that *defines* a residual into a model that *solves* the residual.

For an internal-variable rate ODE $\dot{\boldsymbol{q}} = \boldsymbol{r}(\boldsymbol{q},\boldsymbol{\varepsilon},T)$, **backward Euler** (`BackwardEulerTimeIntegration.cxx:66`) writes the residual:
$$\boldsymbol{R}(\boldsymbol{q}_{n+1}) = \boldsymbol{q}_{n+1} - \boldsymbol{q}_n - \Delta t\,\boldsymbol{r}(\boldsymbol{q}_{n+1}, \boldsymbol{\varepsilon}_{n+1}, T_{n+1})$$

with analytic Jacobian
$$\frac{\partial \boldsymbol{R}}{\partial \boldsymbol{q}_{n+1}} = \boldsymbol{I} - \Delta t\,\frac{\partial \boldsymbol{r}}{\partial \boldsymbol{q}_{n+1}}.$$

`ImplicitUpdate` then drives a Newton solve over $\boldsymbol{q}_{n+1}$ until $\|\boldsymbol{R}\|<\text{tol}$. *Crucially*, $\partial \boldsymbol{r}/\partial \boldsymbol{q}$ is whatever the chained-derivative product across the constituent nodes evaluates to — so a J2 plasticity step uses the analytic flow-rule Jacobian, the Voce hardening Jacobian, and the elasticity tangent without any numerical perturbation.

### 3.2 Why Mandel + multiplicative + analytic?

The combination of **Mandel stress** + **multiplicative split** + **analytically differentiated residuals** is what makes NEML2 a fully implicit, batched, *quadratically convergent* solver on the GPU:

- *Mandel* makes flow rules frame-indifferent without per-step rotation algebra (no Hughes–Winget update needed).
- *Multiplicative* makes the elastic constitutive law strain-energy-based and exactly recoverable, even at large rotations and stretches.
- *Analytic Jacobians* propagate exactly through the DAG (chain rule), so the consistent algorithmic tangent the host code receives is the *exact* derivative of its converged stress with respect to its input strain — not a finite-difference probe. This is the property that gives quadratic FE convergence at the structural level.

### 3.3 The rate-form vs. integrated-form duality

NEML2 supports both formulations of any evolution law via two complementary nodes:

| | Rate form | Integrated form |
|---|---|---|
| Variable | $\dot{\boldsymbol{q}}$ given by constitutive node | $\boldsymbol{q}$ given by constitutive node |
| Time integration | `BackwardEulerTimeIntegration`: $\boldsymbol{R} = \boldsymbol{q}_{n+1} - \boldsymbol{q}_n - \Delta t\,\dot{\boldsymbol{q}}_{n+1}$ | `VariableRate`: $\dot{\boldsymbol{q}} = (\boldsymbol{q}_{n+1} - \boldsymbol{q}_n)/\Delta t$ |
| Use case | Internal variables driven by physics (stress, hardening, $\phi$, $d$) | Inputs driven by FE/loading (strain, temperature) |

This duality is what lets the same `YieldFunction` node serve both an inviscid plasticity model (rate-independent) and a viscoplastic model (rate-dependent), depending on whether $\dot\gamma$ is solved from a consistency condition $f=0$ or from a viscous flow rule like Perzyna.

### 3.4 Regularization choices and where they hide

| Regularizer | Where | Why |
|---|---|---|
| Machine-$\varepsilon$ in tensor norms | `AssociativeJ2FlowDirection.cxx:59`, `ChabochePlasticHardening.cxx:65` | Flow direction $\partial \bar\sigma/\partial \boldsymbol{M}$ is undefined at $\bar\sigma=0$. Adding $\varepsilon$ to the denominator makes the tangent finite without changing the direction anywhere it's physically meaningful. |
| Heaviside on Perzyna | `PerzynaPlasticFlowRate.cxx:63` | Smooth onset of viscoplastic flow; preserves $C^0$ continuity of the residual. |
| $\eta>0$ residual stiffness in $g(d)$ | `PowerDegradationFunction.cxx:59` | Keeps tangent invertible at $d=1$ (fully cracked) — without this the Newton iteration breaks at full damage. |
| $\delta$ in shrinking-core denominator | `DiffusionLimitedReaction.cxx:72–78` | Avoids $0/0$ when $r_i \to r_o$ at vanishing product layer thickness. |
| Logarithmic extension below $S_p$ | `CapillaryPressure.cxx:76–94` | $P_c \to \infty$ at $S_e=0$; replaces the singular branch with a finite log to keep the Newton step bounded near the dry limit. |
| Macaulay smoothing in spectral split | `LinearIsotropicStrainEnergyDensity.cxx:128–139` | Tensile-only energy must use $\langle\cdot\rangle_+$, but this is non-differentiable; the code uses the standard $\max(\cdot,0)$ + autograd subgradient, which is exact almost everywhere. |

These are the entire set of "numerical fudges" in the constitutive layer. Each one is mathematically transparent and isolated to a single line.

### 3.5 Solver-batched constitutive update

The Newton solve inside `ImplicitUpdate` is *batched over all material points simultaneously* — every quadrature point in a finite-element mesh, or every voxel in a PFM grid, advances together. The convergence test uses a global norm $\|\boldsymbol{R}\|$ over the whole batch (`Newton.cxx:75–105` in the architecture layer). The price is that all batch elements have to converge before the iteration stops; the benefit is that a single Newton step is one fused GPU kernel rather than an outer Python loop.

---

## 4. From Equation to Code (Mapping Reference)

A consolidated index of "equation ↔ source location" for the equations derived above. Use this as a quick lookup when reading the codebase.

### 4.1 Continuum mechanics

| Equation | LaTeX | File:Line |
|---|---|---|
| Linear isotropic stress | $\boldsymbol{\sigma}=K\,\mathrm{tr}(\boldsymbol{\varepsilon}_e)\boldsymbol{I}+2G\,\mathrm{dev}\,\boldsymbol{\varepsilon}_e$ | `LinearIsotropicElasticity.cxx:62` |
| Isotropic stiffness | $\mathbb{C}=3K\mathbb{I}_v+2G\mathbb{I}_d$ | `IsotropicElasticityTensor.cxx:60` |
| Cubic stiffness | $\mathbb{C}=C_1\mathbb{I}_1+C_2\mathbb{I}_2+C_3\mathbb{I}_3$ | `CubicElasticityTensor.cxx:62` |
| Green–Lagrange strain | $\boldsymbol{E}=\tfrac{1}{2}(\boldsymbol{F}^\top\boldsymbol{F}-\boldsymbol{I})$ | `GreenLagrangeStrain.cxx:61` |
| Mandel = Cauchy (small strain) | $\boldsymbol{M}=\boldsymbol{\sigma}$ | `IsotropicMandelStress.cxx:49` |

### 4.2 Plasticity (J2 and GTN)

| Equation | LaTeX | File:Line |
|---|---|---|
| J2 yield | $f=\sqrt{2/3}\,(\bar\sigma-\sigma_y-h)$ | `YieldFunction.cxx:68–71` |
| GTN yield | $f=(\bar\sigma/\sigma_e)^2 + 2q_1\phi\cosh(\frac{q_2}{2}\frac{3\sigma_h-\sigma_s}{\sigma_e})-(q_3\phi^2+1)$ | `GTNYieldFunction.cxx:89–91` |
| Associative flow | $\dot{\boldsymbol{\varepsilon}}^p=\dot\gamma\,\partial f/\partial\boldsymbol{M}$ | `AssociativePlasticFlow.cxx:62–67` |
| J2 flow direction | $\boldsymbol{N}_M=\tfrac{3}{2}\,\mathrm{dev}\boldsymbol{M}/\bar\sigma$ | `AssociativeJ2FlowDirection.cxx:60–62` |
| Perzyna VP | $\dot\gamma=\langle f/\eta\rangle^n$ | `PerzynaPlasticFlowRate.cxx:62–66` |
| Gurson cavitation | $\dot\phi=(1-\phi)\,\mathrm{tr}\dot{\boldsymbol{\varepsilon}}^p$ | `GursonCavitation.cxx:63` |

### 4.3 Hardening laws

| Equation | LaTeX | File:Line |
|---|---|---|
| Linear iso. | $h=K\bar\varepsilon^p$ | `LinearIsotropicHardening.cxx:57` |
| Voce iso. | $h=R[1-\exp(-d\bar\varepsilon^p)]$ | `VoceIsotropicHardening.cxx:61` |
| Voce iso. (rate) | $\dot h=\theta_0(1-h/R)\dot\gamma$ | `SlopeSaturationVoceIsotropicHardening.cxx:67` |
| Linear kin. | $\boldsymbol{X}=H\boldsymbol{K}_p$ | `LinearKinematicHardening.cxx:59` |
| Frederick–Armstrong | $\dot{\boldsymbol{X}}=(\tfrac{2}{3}C\boldsymbol{N}_M-g\boldsymbol{X})\dot\gamma$ | `FredrickArmstrongPlasticHardening.cxx:74` |
| Chaboche | $\dot{\boldsymbol{X}}=(\tfrac{2}{3}C\boldsymbol{N}_M-g\boldsymbol{X})\dot\gamma-A\|\boldsymbol{X}\|^{a-1}\boldsymbol{X}$ | `ChabochePlasticHardening.cxx:74` |
| Kocks–Mecking $\eta(T)$ | $\eta=e^B\mu\dot\varepsilon_0^{-kTA/\mu b^3}$ | `KocksMeckingFlowViscosity.cxx:82` |
| Kocks–Mecking $\sigma_y(T)$ | $\sigma_y=e^C\mu$ | `KocksMeckingYieldStress.cxx:63` |

### 4.4 Crystal plasticity

| Equation | LaTeX | File:Line |
|---|---|---|
| Schmid tensor | $\boldsymbol{P}^\alpha=\boldsymbol{s}^\alpha\otimes\boldsymbol{m}^\alpha$ | `CrystalGeometry.cxx:112–116` |
| Resolved shear | $\tau^\alpha=\boldsymbol{\sigma}:(\boldsymbol{Q}\boldsymbol{M}^\alpha\boldsymbol{Q}^\top)$ | `ResolvedShear.cxx:80` |
| Power-law slip | $\dot\gamma^\alpha=\dot\gamma_0\|\tau^\alpha/g^\alpha\|^{n-1}\tau^\alpha/g^\alpha$ | `PowerLawSlipRule.cxx:65` |
| Voce slip hardening | $\dot g=\theta_0(1-g/\tau_f)\sum\|\dot\gamma^\beta\|$ | `VoceSingleSlipHardeningRule.cxx:61` |
| Plastic def. rate | $\boldsymbol{D}^p=\boldsymbol{Q}(\sum\dot\gamma^\alpha\boldsymbol{M}^\alpha)\boldsymbol{Q}^\top$ | `PlasticDeformationRate.cxx:75` |
| Plastic spin | $\boldsymbol{W}^p=\boldsymbol{Q}(\sum\dot\gamma^\alpha\boldsymbol{W}^\alpha)\boldsymbol{Q}^\top$ | `PlasticVorticity.cxx:75` |
| Elastic strain rate | $\dot{\boldsymbol{\varepsilon}}^e=\boldsymbol{D}-\boldsymbol{D}^p+(\boldsymbol{W}\boldsymbol{\varepsilon}^e-\boldsymbol{\varepsilon}^e\boldsymbol{W})$ | `ElasticStrainRate.cxx:96` |
| Orientation rate | $\boldsymbol{\Omega}^e=\boldsymbol{W}-\boldsymbol{W}^p+[\boldsymbol{D}^p,\boldsymbol{\varepsilon}^e]$ | `OrientationRate.cxx:97` |

### 4.5 Phase-field fracture

| Equation | LaTeX | File:Line |
|---|---|---|
| Power degradation | $g(d)=(1-d)^p(1-\eta)+\eta$ | `PowerDegradationFunction.cxx:59` |
| Rational degradation | $g(d)=(1-d)^p/[(1-d)^p+Q(d)](1-\eta)+\eta$ | `RationalDegradationFunction.cxx:66–69` |
| AT1 | $\alpha(d)=d$ | `CrackGeometricFunctionAT1.cxx:51` |
| AT2 | $\alpha(d)=d^2$ | `CrackGeometricFunctionAT2.cxx:52` |
| Energy (no split) | $\psi^+ = \tfrac{1}{2}K(\mathrm{tr}\boldsymbol{\varepsilon})^2+G\,\mathrm{dev}:\mathrm{dev}$ | `LinearIsotropicStrainEnergyDensity.cxx:100–104` |
| Energy (vol-dev) | $\psi^\pm$ with Macaulay on $\mathrm{tr}\boldsymbol{\varepsilon}$ | `LinearIsotropicStrainEnergyDensity.cxx:128–139` |
| Energy (spectral) | $\psi^\pm$ from $\langle\lambda_i\rangle_\pm$ eigenvalues | `LinearIsotropicStrainEnergyDensity.cxx:162–175` |
| KKT (FB form) | $\sqrt{F^2+\dot d^2}-(F+\dot d)=0$ | `small_deformation.i:127–133` |

### 4.6 Porous flow & chemistry

| Equation | LaTeX | File:Line |
|---|---|---|
| Effective saturation | $S_e=(\phi/\phi_{\max}-S_r)/(1-S_r)$ | `EffectiveSaturation.cxx:65` |
| Power-law $K(\varphi)$ | $K=K_0(\varphi/\varphi_0)^p$ | `PowerLawPermeability.cxx:58–62` |
| Kozeny–Carman | $K=K_0\varphi^n(1-\varphi_0)^m/[\varphi_0^m(1-\varphi)^n]$ | `KozenyCarmanPermeability.cxx:62–66` |
| Exponential $K$ | $K=K_0\exp[-a(\varphi-\varphi_0)]$ | `ExponentialLawPermeability.cxx:59–62` |
| van Genuchten $P_c$ | $P_c=a(S_e^{-1/m}-1)^{1-m}$ | `VanGenuchtenCapillaryPressure.cxx:65–69` |
| Brooks–Corey $P_c$ | $P_c=P_t\,S_e^{-1/p}$ | `BrooksCoreyCapillaryPressure.cxx:62–64` |
| Advective stress | $p_s=-c/(3J)\,J_s^{-5/3}J_t^{-2/3}\boldsymbol{P}:\boldsymbol{F}$ | `AdvectiveStress.cxx:79–95` |
| Avrami–Erofeev | $\dot a=k(1-a)(-\ln(1-a))^n$ | `AvramiErofeevNucleation.cxx:59` |
| Contracting geometry | $\dot a=k(1-a)^n$ | `ContractingGeometry.cxx:62` |
| Diffusion-limited | $r=2DR_lR_s/\omega\cdot r_o/(r_o-r_i+\delta)$ | `DiffusionLimitedReaction.cxx:72–78` |
| Channel geometry | $r_i=\sqrt{1-\phi_s-\phi_p},\;r_o=\sqrt{1-\phi_s}$ | `CylindricalChannelGeometry.cxx:63–76` |
| Volume balance | $V=M/(1-\phi_o)\sum\omega_i/\rho_i$ | `EffectiveVolume.cxx:87–95` |

### 4.7 Generic time integration & finite volume

| Equation | LaTeX | File:Line |
|---|---|---|
| Backward Euler residual | $R=\boldsymbol{q}_{n+1}-\boldsymbol{q}_n-\Delta t\,\dot{\boldsymbol{q}}_{n+1}$ | `BackwardEulerTimeIntegration.cxx:66` |
| Forward Euler | $\boldsymbol{q}=\boldsymbol{q}_n+\Delta t\,\dot{\boldsymbol{q}}_n$ | `ForwardEulerTimeIntegration.cxx:69` |
| FD rate | $\dot f=(f-f_n)/(t-t_n)$ | `VariableRate.cxx:69` |
| Von Mises invariant | $\sigma_{\text{vm}}=\sqrt{3/2}\,\|\mathrm{dev}\boldsymbol{A}\|$ | `SR2Invariant.cxx:116` |
| Linear combination | $\boldsymbol{u}=\sum w_i\boldsymbol{v}_i+\boldsymbol{b}$ | `LinearCombination.cxx:135–138` |
| Upwind FV flux | $F_{i+1/2}=v^+u_i+v^-u_{i+1}$ | `FiniteVolumeUpwindedAdvectiveFlux.cxx:65–72` |
| FV gradient | $(\nabla u)_i\approx-\beta(u_{i+1}-u_i)/\Delta x$ | `FiniteVolumeGradient.cxx:70–73` |

---

## 5. Closing Remarks: The Physics Worldview

The library can be read as encoding a particular philosophical commitment about *how the constitutive response of materials should be formulated for computation*:

1. **Pointwise, not field.** Every equation in `models/` is local — a function of the strain, internal variables, and temperature *at a single material point*. Spatial PDE coupling lives outside the library. This makes the library composable across mesh-based, mesh-free, voxel-based, and integration-point-based hosts.

2. **Internal variables, not internal state.** The plastic strain, back-stresses, void fraction, slip strengths, and damage are all *primary unknowns* with their own evolution equations — not bookkeeping derived from the strain history. This is what makes implicit time integration of the constitutive law tractable: every internal variable contributes a residual and a row in the algorithmic Jacobian.

3. **Frame-indifference through Mandel + multiplicative decomposition.** Rather than complicating each flow rule with objective rates, the library lifts the *whole* plastic update into the intermediate Mandel-conjugate configuration. This is invisible to the model author writing a flow rule, but it is the architecture that makes large-rotation crystal plasticity work uniformly with small-rotation J2 plasticity.

4. **Variational fracture, not stress-based.** The phase-field module follows the Bourdin–Francfort–Marigo energy minimization paradigm: fracture is the minimization of $\Pi[\boldsymbol{u},d]$ subject to $\dot d \ge 0$. There is no separate "crack-tip" or "stress-intensity-factor" machinery — fracture mechanics emerges from the variational formulation, including initiation, branching, and coalescence.

5. **Differentiable through-and-through.** Every equation comes with its analytic Jacobian (no autograd-only paths in the constitutive layer). This is what permits the library to claim, simultaneously: implicit time integration with quadratic convergence; full sensitivity analysis for inverse problems / Bayesian calibration; gradient-based optimization for material design; and exact algorithmic tangents for the host FE code.

> **Final Insight — NEML2 is fundamentally a continuum-thermodynamics implementation, not a simulation engine.** What it implements is the *closure relations* between observable kinematic quantities ($\boldsymbol{\varepsilon}, T, t, \nabla d, \ldots$) and conjugate driving forces ($\boldsymbol{\sigma}, \dot\phi, \dot d, \ldots$). The host code provides the balance laws (momentum, mass, energy) and the geometry; NEML2 provides the *physical constitutive content* that closes them. Because every closure is an analytically differentiable, batched node in a DAG, the resulting system is solvable by a single batched Newton iteration on any modern accelerator.
