// SPDX-License-Identifier: GPL-3.0-or-later
//
// ftsim - Finite Theory N-body charge simulator
// Copyright (C) 2026 Phil Bouchard <phil@fornux.com>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// ---------------------------------------------------------------------------
//  ftsim - N random point charges under Finite Theory's time-dilation factor
//
//  f(r) = eta*r / (eta*r + S)        S = |q_source|   (Eq. 93, s -> inf)
//  v^2  = (K q_i q_j /(m_i r)) f^2   (dilated orbital speed squared)
//
//  Two ways to turn v^2 into an acceleration; they are NOT equivalent once
//  f != 1, so both are implemented and switchable at runtime:
//
//    CENTRIPETAL   a = v^2 / r          (Eq.103 derivation)
//                    = C * eta^2 / (eta r + S)^2
//                  sign-definite: no crossover, plateau at C eta^2/S^2
//
//    POTENTIAL     a = -(1/m) d(U f^2)/dr   (Eq.95 route, 1/2 restored)
//                    = -C * eta^2 (S - eta r) / (eta r + S)^3
//                  changes sign at r = S/eta: like charges ATTRACT inside it
//
//  where C = K q_i q_j / m_i.  Both are finite at r = 0, so no softening
//  length is needed - the f factor regularises the Coulomb singularity.
//
//  GRAVITOELECTRIC term: identical kernel, with
//      C_g = -G m_j        (acceleration, already per unit m_i)
//      S_g = m_j           (mass is the scale inside f, Eq.11 analogue)
//      eta_g = c^2/G       -> ETA_G = eta_g*a0/m_e = 7.82254e46
//  In these units G = G m_e^2/(K e^2) = 2.400610e-43, i.e. gravity between
//  two electrons is 43 orders below their electrostatic force.  A GAIN
//  multiplier is provided purely so the term can be SEEN; gain != 1 is a
//  visualisation choice, not physics.
//
//  Units (so the numbers stay O(1)):
//      length   a0      = 5.29177210903e-11 m
//      speed    v_Bohr  = 2.18769126e6 m/s
//      time     a0/v_B  = 2.4189e-17 s
//      charge   e,  mass m_e
//
//  Species:  electron  q = -1, m = 1
//            proton    q = +1, m = 1836.152673
//            neutron   q =  0, m = 1838.683662
//  A neutron is electrically inert automatically: as a SOURCE it gives
//  S = |q| = 0 and ftKernel returns 0; as a TARGET it gives C = 0.  With
//  gravity off it therefore drifts in a straight line - correct, not a bug.
//
//  DISCRETE STATIONARY STATES (optional).  Imposing the dilated standing-wave
//  condition L = n hbar f on the FT orbit balance makes f^2 cancel EXACTLY:
//
//      r_n = n^2 hbar^2/(m_e K q^2 Z) = n^2 a0 / Z      (no FT correction)
//      v_n = Z/n  v_Bohr,   omega_n = v_n/r_n = Z^2/n^3
//
//  f survives only in the ENERGY:
//
//      f_n = ETA n^2 / (ETA n^2 + Z^2)
//      E_n = -(Z^2 / 2n^2) * f_n^2   Hartree      (1 Ha = 27.211386 eV)
//
//  so the FT shift goes as 1/n^4.
//
//  RETROFIT OF lambda_e.  For S states the FT 1/r^2 term is exactly degenerate
//  with the proton finite-size term (both go as 1/n^3), so in electronic
//  hydrogen it appears as a LARGER apparent proton radius:
//      r_p^2(apparent) = r_p^2 + k a r_c,   r_c = e/lambda_e
//  k = 6.00 (1S-3S), 6.33 (2S-4P), 4.00 (2S-2P Lamb), 6.02 (2S-8D) after the
//  Rydberg is eliminated with 1S-2S.  Muonic hydrogen (a 186x smaller) is
//  immune.  A joint fit of r_p and r_c to the muonic radius 0.84087(39) fm
//  and five electronic measurements (Beyer 2017, Fleurbaey 2018, Bezginov
//  2019, Grinin 2020, Brandt 2022) gives
//      r_c = (5.1 +/- 1.5)e-23 m   ->  lambda_e = 3.2e3 C/m,  ETA = 1.04e12
//      95% limit r_c < 8.0e-23 m   ->  lambda_e > 2.0e3 C/m,  ETA > 6.6e11
//  The 3.4 sigma is the proton-radius-puzzle remnant and should be read as
//  an upper limit, not a detection: chi2/dof stays 8-12/3 whichever
//  experiment is dropped.  This supersedes the earlier "10 Hz" floor of
//  6.2e14, which wrongly compared 1S-2S precision against theory while
//  R_inf and r_p are themselves extracted from hydrogen.  The manuscript's
//  lambda_e = c^2/sqrt(KG) (ETA 3.8e25) is 13 orders above the retrofit.  The quantisation is IMPOSED here,
//  not derived: the FT force alone has a degree-1 numerator and therefore
//  exactly one root, so it cannot generate a series of radii.
//
//  DISCRETE CROSSOVER (a fourth route).  Instead of one crossover at S/eta,
//  the crossover is ROUNDED to the nearest stationary shell:
//      n     = round(sqrt(r*S))         since r_n = n^2/Z and S = Z
//      eta_L = S^2 / n^2                so that S/eta_L = n^2/Z = r_n
//  then the Potential formula is applied with eta_L. The force therefore has
//  a STABLE zero at every shell r_n = n^2 a0/Z and an unstable rim between
//  shells at (n+1/2)^2. Inside the ground shell it is repulsive, so nothing
//  collapses onto the nucleus.
//
//  Two things this route does NOT do, and they matter:
//    - it does not recover Coulomb. Because the crossover tracks r, the force
//      stays within a few percent of zero everywhere: |a|/|a_Coulomb| runs
//      about 0.03 - 0.15 rather than 1. A law with a zero at every shell
//      cannot also be inverse-square.
//    - it does not derive quantisation. The shells come from L = n hbar, and
//      this route hard-codes them into the force law. It is the "ditches"
//      idea expressed analytically, not a mechanism.
//
//  STRONG NUCLEAR FORCE (optional).  Malfliet-Tjon V, a standard NN
//  parametrisation, acting between NUCLEONS ONLY (protons and neutrons):
//
//      V(r) = [ V_R exp(-mu_R r) - V_A exp(-mu_A r) ] / r     MeV, r in fm
//      V_R = 1438.7200 MeV fm,  mu_R = 3.11 fm^-1   (repulsive core)
//      V_A =  626.8850 MeV fm,  mu_A = 1.55 fm^-1   (attraction)
//      F   = -dV/dr
//        = V_R e^{-mu_R r}(mu_R/r + 1/r^2) - V_A e^{-mu_A r}(mu_A/r + 1/r^2)
//
//  Verified properties: well depth -77.5 MeV at 0.819 fm, V = 0 at 0.5325 fm,
//  F = 0 at 0.8188 fm.  Conversions: r_fm = r_sim * a0/fm = r_sim * 52917.721,
//  and 1 MeV/fm = 1.944690e9 simulator force units.
//
//  It is CHARGE-INDEPENDENT by construction: pp, pn and nn get the identical
//  potential.  That is the physical point - it is why the deuteron (pn) binds
//  while the FT crossover, which couples to charge, gives pn exactly zero.
//  Electrons and composite nuclei are excluded (nucleon flag).
//
//  WARNING: the force is ~10^11 simulator units at 1 fm, so a proton sees
//  ~8e7 in acceleration.  Use dt <= 1e-9 or the integrator will explode.
//
//  MAGNETIC TERMS (optional).  The simulator's units are atomic units
//  (hbar = m_e = e = K = 1), so c = 1/alpha = 137.035999084 and the magnetic
//  constant mu0/4pi = K/c^2 = alpha^2 = 5.32513545e-5.  The motion is planar,
//  so every magnetic field here points along z and reduces to a scalar B_z.
//
//  "moving-charge field" - Biot-Savart for each moving point charge,
//      B_z(i) = alpha^2 q_j (v_j x r)_z / r^3,   r = r_i - r_j
//  acting on charge i through the Lorentz force a = q_i (v_i x B_z z)/m_i.
//  Relative to Coulomb this is O(v^2/c^2) ~ 5e-5 for a ground-state orbit.
//  Magnetic forces between point charges do NOT obey Newton's third law -
//  the missing momentum is carried by the field - so total momentum drifts.
//
//  "spin" - each particle carries a fixed intrinsic moment mu along +/- z:
//      electron  |mu| = (g_e/2) mu_B = 0.500579826 a.u., antiparallel to spin
//      proton    +2.79284734 mu_N,  neutron  -1.91304273 mu_N
//      nuclei    He-4: 0,  Li-7: +3.256427 mu_N,  Be-9: -1.17749 mu_N
//  (mu_B = 1/2, mu_N = mu_B m_e/m_p = 2.7230851e-4 a.u.)
//  Three effects follow, all through F = mu_z grad(B_z):
//    dipole-dipole   F_r = 3 alpha^2 mu_i mu_j / r^4   (parallel repel in-plane)
//    spin-orbit      moment i sees the field of charge j moving at v_j - v_i;
//                    the electron's own-motion part carries the Thomas factor
//                    (g_e - 1)/g_e = 0.500579
//    dipole field    B_z = -alpha^2 mu_j / r^3 acting on moving charges
//  A z-moment in a z-field feels no torque (mu x B = 0), so freezing the spin
//  direction is EXACT in this planar geometry, not an approximation.
//
//  GRAVITOMAGNETIC FIELD (optional).  The mass-current analogue of the
//  moving-charge field: every moving mass m_j produces
//      B_g,z(i) = (G/c^2) m_j (v_j x r)_z / r^3 * f_g(r)^2
//  and mass i responds with
//      a_i = -k (v_i x B_g z)                (independent of m_i)
//  The minus sign is the gravity sign flip: parallel mass currents REPEL,
//  where parallel like electric currents attract.  Relative to Newton the
//  force is -k v_i v_j / c^2.  In these units G/c^2 = G_UNITS * alpha^2 =
//  1.27836e-47, so at gain 1 it is ~5e-47 of the electric force.
//
//  The coupling k depends on the convention, and the choice is exposed:
//      k = 4    linearised GR (h_0i = -4 G m v_i / c^3 r, Mashhoon GEM)
//      k = 1.5  the manuscript's mu_g = 6 pi G/c^2 (Eq. 19)
//      k = 1    naive EM analogy, mu_g = 4 pi G/c^2
//  This is the gravitomagnetic piece only, not the full 1PN EIH dynamics,
//  which also carries v^2 corrections to the gravitoelectric term.
//  f_g^2 is applied for consistency with the FT gravitoelectric kernel; it
//  differs from 1 only within ~1e-47 a0 of a mass and is 1 on any visible scale.
//
//  QUARKS AND OTHER LEPTONS - electric charge only.  They interact through
//  the electric FT kernel (with S = |q|, so fractional charges get fractional
//  crossovers), the moving-charge and spin magnetic terms, and gravity.  They
//  are NOT nucleons, so the strong force does not act on them, and there is
//  no colour, no confinement and no weak interaction.  Free quarks are never
//  observed; this is a charge-only toy, not QCD.
//      species   q      m (m_e)       mu per unit spin (a.u.)
//      muon     -1      206.7682830   -2.420985e-3   (g_mu/2 = 1.00116592)
//      tau      -1      3477.23       -1.437926e-4   (g = 2)
//      nu_e/mu/tau 0    1.956951e-7   0              (0.1 eV PLACEHOLDER:
//                                                     only upper bounds known)
//      up      +2/3     4.227015      +7.885786e-2   | PDG 2022 central
//      down    -1/3     9.138962      -1.823694e-2   | current-quark masses,
//      strange -1/3     182.7792      -9.118468e-4   | MS-bar; Dirac moments
//      charm   +2/3     2485.328      +1.341205e-4   | q/(2m), g = 2.
//      bottom  -1/3     8180.056      -2.037476e-5   | Constituent u/d masses
//      top     +2/3     337945.9      +9.863512e-7   | are ~336 MeV instead.
//  Initial speeds follow equipartition but are capped at 0.5c (68.518 v_B):
//  neutrinos would otherwise start superluminal, and the integrator is
//  non-relativistic.  In stationary-state mode muons and taus bind with the
//  REDUCED mass mu = m M/(m+M): r_n = n^2/(Z mu), E_n = -Z^2 mu/(2n^2) f^2.
//  Muonic hydrogen (mu = 185.8408) sits 185.84x closer than ordinary hydrogen
//  and binds at -2528.5 eV.  (The level table panel stays at infinite nuclear
//  mass for the electron, so its numbers match the earlier verification.)
//
//  QCD BINDING (optional).  Two separate pieces, because "binding energy"
//  does two jobs in a real hadron:
//
//  1. CONFINEMENT - Cornell potential between quarks of the SAME hadron,
//     using the baryon "1/2 rule" (qq in a colour-singlet baryon carries half
//     the q-qbar strength):
//         V_qq(r) = -(2/3) alpha_s hbar c / r + (sigma/2) r
//         F_qq(r) = -[ (2/3) alpha_s hbar c / r^2 + sigma/2 ]   (attractive)
//     alpha_s = 0.35, sigma = 0.18 GeV^2 = 912.19 MeV/fm, so the long-range
//     pull is a CONSTANT 456.10 MeV/fm = 7.307e4 N per pair - the flux tube.
//     It acts only within a hadron (same colour singlet).  Between hadrons
//     the linear term would pull at any distance, which real string breaking
//     prevents, so it is switched off there.  Loose quarks feel charge only.
//     The absolute V carries an arbitrary constant in Cornell fits, so the E
//     readout is meaningful for changes, not as a mass.
//
//     KINEMATICS.  DEFAULTS: Potential force route with FT 1/2 kinematics
//     (p = m v/(1 - v^2/2c^2)).  Alternative: the positional FT dilation
//     f^2 times the mass of
//     the kinetic energy,  p = m v / f^2,  m = m0 + K/c^2,
//         f_i = eta / (eta + sum_j |q_j|/r_ij)     (single source: eta r/(eta r+S))
//     K is taken as the work actually done, dK = v dp.  Integrating
//     dm = dK/c^2 with p = m v gives dm/m = v dv/(c^2 - v^2), i.e.
//         m = m0 / sqrt(1 - v^2/c^2)        EXACTLY
//     (numerically 2.294150 vs 2.294157 at 0.9c).  So the rule m = m0 + K/c^2,
//     applied consistently, IS the special-relativistic mass, and the
//     default momentum is p = gamma m0 v / f^2.  Truncating to K = m0 v^2/2
//     ("FT f^2 first order") removes the speed limit: muon g-2 would orbit at
//     3.71c with a 40.16 ns period, against 149.1 ns measured.
//     With f^2 in the momentum, the matching force is plain COULOMB (route
//     default changed to it): a circular orbit then has v = sqrt(Kq^2/mr) f,
//     the Eq. 103 / worksheet result, exactly.  Putting f^2 in the force too
//     (Potential, Centripetal, Discrete routes) counts it twice.
//     WHAT IS INTEGRATED: (gamma m0 / f^2) dv/dt = F, with f re-evaluated at
//     each new position (reading A).  It reproduces Eq. 103 circular orbits
//     exactly, but conserves NEITHER energy (+7.9e-3 on an eccentric orbit at
//     ETA = 2) NOR, with more than two bodies, momentum (|dP| = 572 in the
//     default scene over 20 time units) - both dt-independent, so law-level.
//     No reading of "inertia m0/f^2 + Coulomb" keeps all three properties:
//       A  as above             Eq.103 yes  momentum no   energy no
//       B  dp/dt = F, p stored  Eq.103 yes  momentum yes  energy no - runs
//                               away: p gained where m0/f^2 is large flings
//                               the particle out where it is small
//       C  Lagrangian           Eq.103 NO (v = sqrt(Kq^2/mr) sqrt f)
//          L = m0 v^2/2f^2 - U  momentum yes energy yes - but m0 r^2/f^2 ->
//                               m0/eta^2 as r -> 0 erases the centrifugal
//                               barrier: an orbit 10% below circular falls
//                               into the nucleus in finite time
//
//     The velocity-only alternatives remain selectable.  FT with the 1/2 factor: Hypothesis 2 normalised by the rest
//     energy m c^2 (not m c^2/2), so Eq. 38 becomes 1 - v^2/(2c^2) and
//         p = m v / (1 - v^2/(2c^2)),   K = m c^2 [2/D - 2 + ln D]
//     This agrees with SR through v^4 (K = m v^2/2 + 3 m v^4/8c^2) and moves
//     the GPS cancellation altitude from 6378 km (Sec. 2.9.2) to R/2 =
//     3185.5 km, matching GR.  But D = 0 at v = sqrt(2) c, so the speed limit
//     becomes 1.414c: in the muon g-2 ring (p = 29.28 m c) it predicts
//     v = 1.380c and a 107.97 ns cyclotron period against 149.1 ns measured.
//     The published Eq. 41 form is kept selectable as "FT Eq. 41":  FT's factor is the square of SR's, so momentum
//     grows faster with speed: a uud triplet binds tighter under FT (2.07 fm,
//     0.972c) than under SR (2.22 fm, 0.9987c), energy conserved to 7e-8.
//     The conserved FT kinetic energy is m c^2 [G - 1 - ln(G)/2], NOT the
//     m c^2 (G - 1) of Eq. 52: that one is m v^2 at low speed, twice the
//     Newtonian value, so Eq. 52 is inconsistent with Eq. 41.
//
//     A NON-Newtonian factor is required inside hadrons: released from
//     rest at proton size, the flux tube hands the quarks ~1993 MeV against
//     1012 MeV of constituent mass, i.e. gamma ~ 3, v ~ 0.94c.  A Newtonian
//     integrator then sends them past c and the energy error explodes.  With
//     gamma-correct kicks a uud triplet stays within 2.22 fm and conserves
//     energy to ~1e-7 (dt = 1e-12).  Light-quark motion in real hadrons is
//     relativistic too, so this is physics, not a numerical patch.
//
//  2. MASS - in a classical simulation potential energy cannot become inertia,
//     so the 99% of the proton mass that is QCD field energy is represented
//     the standard way, by CONSTITUENT quark masses (Griffiths):
//         u 336, d 340, s 486, c 1550, b 4730 MeV   (top does not hadronise,
//         stays at its current mass).  uud = 1012 MeV vs 938.27, the residual
//         ~74 MeV being the colour-magnetic hyperfine term not modelled here.
//     The same masses reproduce mu_p = (4 mu_u - mu_d)/3 = 2.789 mu_N
//     (measured 2.793) - that is why they are the standard choice.
//
//  HYDROGEN-LIKE IONS.  Presets place a nucleus (Z protons, N neutrons) plus
//  ONE electron on a circular orbit:
//        r = a0/Z    ->  1/Z in these units
//        v = Z*alpha*c/n  ->  Z/n in these units (v_Bohr = alpha*c)
//    ion    Z  N   r (a0)    v (v_Bohr)
//    H      1  0   1.0000     1
//    He+    2  2   0.5000     2
//    Li2+   3  4   0.3333     3
//    Be3+   4  5   0.2500     4
//  The nucleus is ONE composite particle by default (q = +Z, m = Z*m_p+N*m_n)
//  because this simulator has no strong force.  "resolve nucleons" places the
//  nucleons individually instead - they then fly apart, which is the correct
//  outcome of Coulomb repulsion with no binding term, not a bug.  The one
//  exception is the Potential route with ETA low enough that the crossover
//  S/ETA reaches nuclear range (~1.9e-5 a0 = 1 fm needs ETA ~ 5.3e4): there
//  like charges attract at short range and the cluster can hold together.
//  In these units K = 1 exactly, and eta is entered as ETA = eta*a0/e.
//  eta must be a CHARGE PER LENGTH (C/m) so that eta*x has the units of q.
//  The manuscript value is lambda_e = lambda_g/lambda_eg = c^2/sqrt(K G)
//  = 1.160427e17 C/m = Planck charge / Planck length, i.e. ETA = 3.8327e25:
//  its crossover is 2.6e-26 a0 and nothing visible happens.  (c*eps0, used
//  in earlier worksheets, is 1/Z0 = 2.654e-3 SIEMENS, not C/m - it is not a
//  valid eta and has been removed.)  Lower ETA to see FT structure.
// ---------------------------------------------------------------------------

#include <QApplication>
#include <QWidget>
#include <QPainter>
#include <QPainterPath>
#include <QVector>
#include <QRectF>
#include <QFontMetrics>
#include <QTimer>
#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QFormLayout>
#include <QLabel>
#include <QSlider>
#include <QSpinBox>
#include <QDoubleSpinBox>
#include <QComboBox>
#include <QCheckBox>
#include <QPushButton>
#include <QGroupBox>
#include <QScrollArea>
#include <QScrollBar>
#include <QSplitter>
#include <QPointF>
#include <QRandomGenerator>
#include <QWheelEvent>
#include <QMouseEvent>
#include <QToolTip>
#include <QTabBar>
#include <QGridLayout>
#include <QDial>
#include <QAction>
#include <QMenu>
#include <QMenuBar>
#include <QRegularExpression>
#include <QValidator>
#include <QFont>
#include <QStringList>
#include <vector>
#include <cmath>
#include <deque>

// ---------------------------------------------------------------------------
// physics
// ---------------------------------------------------------------------------

enum class Route  { Centripetal, Potential, Coulomb, DiscreteCrossover };

// floor-quantisation, applied on top of whichever Route is selected:
//   None    a is used as computed
//   Radius  a is evaluated at r_q = floor(r/D)*D  (space comes in lumps)
//   Force   a_q = floor(a/D)*D                     (the force comes in lumps)
enum class Quant  { None, Radius, Force };

// how the initial nucleon positions are drawn
//   Random         uniform in the box, no regard for the crossover
//   ProtonCluster  protons packed inside the attractive region r < S/eta
//   MixedCluster   protons AND neutrons packed inside it
// The attractive region exists only on the Potential route: there like
// charges attract for r < S/eta.  On Centripetal and Coulomb the clusters
// are simply tight and will blow apart.
enum class Layout { Random, ProtonCluster, MixedCluster };

// masses in units of m_e
static constexpr double M_ELECTRON = 1.0;
static constexpr double M_PROTON   = 1836.152673;
static constexpr double M_NEUTRON  = 1838.683662;

// charge, mass (m_e) and intrinsic moment per unit spin (atomic units)
struct Species {
    const char *name;
    const char *qText;
    double      q, m, mu0;
    const char *color;
    bool        nucleon;
};
enum SpeciesId { SP_E = 0, SP_P, SP_N, SP_MU, SP_TAU, SP_NUE, SP_NUMU, SP_NUTAU,
                 SP_U, SP_D, SP_S, SP_C, SP_B, SP_T, SP_STAR, SP_PLANET, SP_COUNT };
static const Species SPECIES[SP_COUNT] = {
    {"electron", "-1",   -1.0,       1.0,          -5.005798e-01, "#3FB4CF", false},
    {"proton",   "+1",   +1.0,       1836.152673,  +7.605161e-04, "#C56A00", true },
    {"neutron",  " 0",    0.0,       1838.683662,  -5.209378e-04, "#8892A4", true },
    {"muon",     "-1",   -1.0,       206.7682830,  -2.420985e-03, "#2E8FD0", false},
    {"tau",      "-1",   -1.0,       3477.23,      -1.437926e-04, "#1F5FA8", false},
    {"nu_e",     " 0",    0.0,       1.956951e-07,  0.0,          "#C9D3E0", false},
    {"nu_mu",    " 0",    0.0,       1.956951e-07,  0.0,          "#B3BFD0", false},
    {"nu_tau",   " 0",    0.0,       1.956951e-07,  0.0,          "#9DABC0", false},
    {"up",       "+2/3",  2.0 / 3.0, 4.227015,     +7.885786e-02, "#F2C14E", false},
    {"down",     "-1/3", -1.0 / 3.0, 9.138962,     -1.823694e-02, "#C77DDB", false},
    {"strange",  "-1/3", -1.0 / 3.0, 182.7792,     -9.118468e-04, "#A35FC0", false},
    {"charm",    "+2/3",  2.0 / 3.0, 2485.328,     +1.341205e-04, "#E89B2F", false},
    {"bottom",   "-1/3", -1.0 / 3.0, 8180.056,     -2.037476e-05, "#7E45A0", false},
    {"top",      "+2/3",  2.0 / 3.0, 337945.9,     +9.863512e-07, "#D9772A", false},
    // neutral astronomical bodies (mass only): 1 M_sun and 1 M_earth in m_e
    {"star",     " 0",    0.0,       2.1833749e60,  0.0,          "#FFD966", false},
    {"planet",   " 0",    0.0,       6.5560967e54,  0.0,          "#5DADE2", false},
};
static constexpr double V_CAP = 68.518;   // 0.5 c in v_Bohr

// default starting scene used by Reset
static constexpr int DEFAULT_NE = 20, DEFAULT_NP = 20, DEFAULT_NN = 20;

// constituent quark masses (m_e), Griffiths: u 336, d 340, s 486, c 1550,
// b 4730 MeV; top kept at its current mass (no hadrons form around it)
static double constituentMass(int sp)
{
    switch (sp) {
    case SP_U: return 657.5356;
    case SP_D: return 665.3634;
    case SP_S: return 951.0783;
    case SP_C: return 3033.274;
    case SP_B: return 9256.379;
    default:   return SPECIES[sp].m;
    }
}
static bool isQuark(int sp) { return sp >= SP_U && sp <= SP_T; }

// on-screen symbol for each species, plus an optional subscript.  Standard
// particle symbols rather than literal first letters, because the literal
// letters collide: neutron/neutrino would both be 'n', tau/top both 't'.
//   e p n  mu tau  nu_e nu_mu nu_tau  u d s c b t
static const char *SYM[SP_COUNT] = {
    "e", "p", "n", "\xCE\xBC", "\xCF\x84",
    "\xCE\xBD", "\xCE\xBD", "\xCE\xBD",
    "u", "d", "s", "c", "b", "t", "S", "P"};
static const char *SUB[SP_COUNT] = {
    "", "", "", "", "",
    "e", "\xCE\xBC", "\xCF\x84",
    "", "", "", "", "", "", "", ""};

struct Particle {
    QPointF pos, vel, acc;
    double  q = 1.0;      // charge, units of e
    double  m = 1.0;      // mass,   units of m_e
    bool    nucleon = false;  // proton or neutron: eligible for the strong force
    int     species = -1;     // index into SPECIES, -1 = composite nucleus
    int     hadron  = -1;     // colour-singlet id; quarks sharing it are confined
    int     spin    = 0;      // +1 / -1 along z, 0 = no spin assigned
    double  mu      = 0.0;    // intrinsic magnetic moment along z, atomic units
    double  fpos    = 1.0;    // positional FT dilation f at this particle
    int     level  = 0;   // principal quantum number, 0 = not quantised
    int     anchor = -1;  // index of the nucleus it is bound to
    double  phase  = 0.0; // orbital phase, radians
    std::deque<QPointF> trail;
};

class Sim {
public:
    static constexpr double G_UNITS = 2.400610e-43;  // G m_e^2/(K e^2)

    // magnetic constants in atomic units
    static constexpr double ALPHA2   = 5.32513545e-5;          // 1/c^2
    static constexpr double G_E      = 2.00231930436;
    static constexpr double THOMAS_E = (G_E - 1.0) / G_E;       // 0.500579
    static constexpr double MU_B     = 0.5;
    static constexpr double MU_N     = 0.5 / 1836.152673;       // 2.7230851e-4
    static constexpr double MU_E     = 0.5 * G_E * MU_B;        // 0.500579826
    static constexpr double MU_P     = 2.79284734 * MU_N;
    static constexpr double MU_NEU   = -1.91304273 * MU_N;

    bool   magnetic = true;     // Biot-Savart field of moving charges
    bool   spinOn   = false;    // intrinsic moments
    int    spinMode = 0;        // 0 random, 1 all up, 2 alternating
    double magGain  = 1.0;      // display amplifier, 1 = true strength

    bool   gravMag  = true;     // gravitomagnetic field of moving masses

    // EXTERNAL UNIFORM FIELDS (atomic units).  The motion is planar, so:
    //   E, g   in-plane vectors (magnitude + in-plane angle from +x)
    //   B, B_g only their z-component acts (magnitude * cos(tilt from +z));
    //          an in-plane B or B_g would push out of the plane, which is
    //          not simulated
    //   a  = q/m (E + v x B)      E in 5.14220675e11 V/m, B in 2.35051757e5 T
    //   a += g                    g in 9.0e22 m/s^2 (acts on every mass)
    //   a += v x B_g              B_g in 4.1341e16 s^-1; equivalent to a
    //                             frame rotating at Omega = B_g/2 (Coriolis)
    double extEx = 0, extEy = 0, extBz = 0;
    double extGx = 0, extGy = 0, extBgz = 0;
    bool hasExternal() const
    { return extEx || extEy || extBz || extGx || extGy || extBgz; }

    void externalPass()
    {
        for (auto &a : p) {
            if (a.m <= 0.0) continue;
            const double qm = a.q / a.m;
            a.acc += QPointF(qm * (extEx + a.vel.y() * extBz),
                             qm * (extEy - a.vel.x() * extBz));
            a.acc += QPointF(extGx + a.vel.y() * extBgz,
                             extGy - a.vel.x() * extBgz);
        }
    }
    double gmK      = 4.0;      // coupling: 4 GR, 1.5 manuscript Eq.19, 1 naive

    std::vector<Particle> p;
    double eta     = 1.0;         // ETA   = eta_e*a0/e
    double etaG    = 7.822540e46; // ETA_G = eta_g*a0/m_e  (eta_g = c^2/G)
    double gravGain = 1.0;        // display amplifier, 1 = true strength
    bool   gravity  = true;
    Route  route    = Route::Potential; // default restored: FT force route
    Layout layout   = Layout::MixedCluster;
    double clusterFill = 0.90;   // pack so every pair sits at <= fill * r_c
    double t        = 0.0;

    // kinematic time-dilation factor used in the momentum update
    //   FTHalf     p = m v / (1 - v^2/2c^2)   FT with the 1/2 factor (default):
    //              Hypothesis 2 with KE normalised by m c^2, not m c^2/2
    //   FT         p = m v / (1 - v^2/c^2)    manuscript Eq. 41 as published
    //   SR         p = m v / sqrt(1 - v^2/c^2)
    //   Newtonian  p = m v
    //   F2SC       p = m v / f^2,  m = m0 + K/c^2 with K = work done
    //              -> m = m0/sqrt(1 - v^2/c^2) exactly, so p = gamma m0 v / f^2
    //   F2Lin      p = m v / f^2,  m = m0 (1 + v^2/2c^2)   (K = m0 v^2/2 only)
    //   F2Pos      p = m0 v / f^2  (positional dilation only)
    enum class Kin { Newtonian, FT, SR, FTHalf, F2SC, F2Lin, F2Pos };
    Kin    kin = Kin::FTHalf;   // default restored: FT with the 1/2 factor
    bool positionalKin() const
    { return kin == Kin::F2SC || kin == Kin::F2Lin || kin == Kin::F2Pos; }
    bool   qcd         = false;  // Cornell confinement within each hadron
    bool   constituent = false;  // constituent quark masses
    double qcdGain     = 1.0;    // display damper
    int    nextHadron  = 0;

    double quarkMass(int sp) const
    { return constituent ? constituentMass(sp) : SPECIES[sp].m; }

    // Cornell, baryon 1/2 rule; r in a0 -> simulator units
    static double qcdForce(double r_sim)
    {
        static constexpr double A0_OVER_FM   = 52917.721090;
        static constexpr double MEVFM_TO_SIM = 1.944690e9;
        static constexpr double K_C   = 46.04296;       // (2/3) a_s hbar c, MeV fm
        static constexpr double HALFS = 456.0958;       // sigma/2, MeV/fm
        double r = r_sim * A0_OVER_FM;
        if (r < 0.01) r = 0.01;
        return -(K_C / (r * r) + HALFS) * MEVFM_TO_SIM; // + = repulsive
    }
    static double qcdPotential(double r_sim)
    {
        static constexpr double A0_OVER_FM = 52917.721090;
        static constexpr double MEV_TO_SIM = 36749.32217565;
        static constexpr double K_C = 46.04296, HALFS = 456.0958;
        double r = r_sim * A0_OVER_FM;
        if (r < 0.01) r = 0.01;
        return (-K_C / r + HALFS * r) * MEV_TO_SIM;
    }

    // switch every existing quark between current and constituent mass
    void applyQuarkMasses()
    {
        for (auto &a : p)
            if (isQuark(a.species)) {
                a.m = quarkMass(a.species);
                if (a.spin != 0) a.mu = momentFor(a) * a.spin;
            }
        computeAcc();
    }

    bool   strong     = false;   // Malfliet-Tjon V between nucleons
    double strongGain = 1.0;     // display damper; 1 = true strength

    Quant  quant     = Quant::None;
    double quantStep = 0.25;    // D, in a0 for Radius, force units for Force

    bool   quantised = false;   // hold electrons on stationary states
    int    nMax      = 3;       // levels handed out round-robin, 1..nMax

    // crossover radius for a unit charge source (protons), in a0
    double crossover() const { return (eta > 0) ? 1.0 / eta : 0.0; }

    // --- stationary-state relations (f cancels in r and v, survives in E) ---
    // m = orbiting particle mass in m_e (1 = electron, 206.77 = muon, ...)
    // a0 scales as 1/m, speeds do not, energies scale as m
    static double r_n(int n, double Z, double m = 1.0) { return n * n / (Z * m); }
    static double v_n(int n, double Z)                 { return Z / double(n); }
    static double w_n(int n, double Z, double m = 1.0) { return Z * Z * m / double(n*n*n); }
    double f_n(int n, double Z, double m = 1.0) const
    {
        const double num = eta * n * n;           // f = eta r_n/(eta r_n + Z)
        return num / (num + Z * Z * m);
    }
    double E_n(int n, double Z, double m = 1.0) const                   // Hartree
    {
        const double f = f_n(n, Z, m);
        return -(Z * Z * m) / (2.0 * n * n) * f * f;
    }
    double E_bohr(int n, double Z, double m = 1.0) const
    { return -(Z * Z * m) / (2.0 * n * n); }

    // bind every electron to its nearest positive charge and hand out levels
    void assignLevels()
    {
        int next = 1;
        for (auto &e : p) {
            e.level = 0; e.anchor = -1;
            const bool lepton = (e.species == SP_E || e.species == SP_MU ||
                                 e.species == SP_TAU);
            if (!lepton) continue;
            double best = 1e300; int bi = -1;
            for (int j = 0; j < int(p.size()); ++j) {
                if (p[j].q <= 0.0) continue;
                double r = std::hypot(e.pos.x() - p[j].pos.x(),
                                      e.pos.y() - p[j].pos.y());
                if (r < best) { best = r; bi = j; }
            }
            if (bi < 0) continue;
            e.anchor = bi;
            e.level  = next;
            next = (next % nMax) + 1;
            e.phase  = std::atan2(e.pos.y() - p[bi].pos.y(),
                                  e.pos.x() - p[bi].pos.x());
        }
    }

    // advance the quantised electrons kinematically: no radial equation, so
    // they neither radiate nor spiral - the Bohr postulate, made explicit
    void stepQuantised(double dt)
    {
        for (auto &e : p) {
            if (e.level <= 0 || e.anchor < 0 ||
                e.anchor >= int(p.size())) continue;
            const Particle &nuc = p[e.anchor];
            const double Z = nuc.q;
            if (Z <= 0.0) continue;
            const double mu = e.m * nuc.m / (e.m + nuc.m);    // reduced mass
            const double r = r_n(e.level, Z, mu);
            const double w = w_n(e.level, Z, mu);
            e.phase += w * dt;
            const double cx = std::cos(e.phase), sy = std::sin(e.phase);
            e.pos = nuc.pos + QPointF(r * cx, r * sy);
            e.vel = nuc.vel + QPointF(-r * w * sy, r * w * cx);
        }
    }

    // total quantised energy, in Hartree, and the shift against pure Bohr
    void quantisedEnergy(double &E, double &shift) const
    {
        E = shift = 0.0;
        for (const auto &e : p) {
            if (e.level <= 0 || e.anchor < 0) continue;
            const double Z = p[e.anchor].q;
            if (Z <= 0.0) continue;
            const double M  = p[e.anchor].m;
            const double mu = e.m * M / (e.m + M);            // reduced mass
            E     += E_n(e.level, Z, mu);
            shift += E_n(e.level, Z, mu) - E_bohr(e.level, Z, mu);
        }
    }

    // diagnostic: how many charged pairs are inside the attractive region
    void attractivePairs(int &inside, int &total) const
    {
        inside = total = 0;
        const int n = static_cast<int>(p.size());
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                if (p[i].q == 0.0 || p[j].q == 0.0) continue;
                ++total;
                double r = std::hypot(p[i].pos.x() - p[j].pos.x(),
                                      p[i].pos.y() - p[j].pos.y());
                if (r < pairScaleQ(p[i].q, p[j].q) / eta) ++inside;
            }
    }

    void reset(int nE, int nP, int nN,
               double boxHalf, double vScale, quint32 seed)
    {
        QRandomGenerator rng(seed);
        p.clear();
        nextHadron = 0;
        p.reserve(nE + nP + nN);

        // uniform point in a disk of radius R (sqrt keeps the density flat)
        auto inDisk = [&](double R) {
            double u  = rng.generateDouble();
            double th = 2.0 * M_PI * rng.generateDouble();
            double rr = R * std::sqrt(u);
            return QPointF(rr * std::cos(th), rr * std::sin(th));
        };
        auto inBox = [&](double H) {
            return QPointF(rng.bounded(2.0 * H) - H, rng.bounded(2.0 * H) - H);
        };

        // every pair inside a disk of radius R is separated by at most 2R,
        // so R = fill * r_c / 2 guarantees all pairs sit within fill * r_c
        // floor at 1 fm: with lambda_e the crossover is 2.6e-26 a0, and
        // packing nucleons inside it would start them 1e-36 m apart
        const double R = std::max(0.5 * clusterFill * crossover(), 1.8897e-5);
        const bool clusterP = (layout != Layout::Random);
        const bool clusterN = (layout == Layout::MixedCluster);

        auto add = [&](int n, int sp, bool cluster) {
            const Species &S = SPECIES[sp];
            const double vs = std::min(vScale / std::sqrt(S.m), V_CAP);
            for (int i = 0; i < n; ++i) {
                Particle a;
                a.pos = cluster ? inDisk(R) : inBox(boxHalf);
                a.vel = QPointF((rng.bounded(2.0) - 1.0) * vs,
                                (rng.bounded(2.0) - 1.0) * vs);
                a.q = S.q;  a.m = S.m;
                a.nucleon = S.nucleon;
                a.species = sp;
                p.push_back(a);
            }
        };
        add(nE, SP_E, false);                      // electrons never clustered
        add(nP, SP_P, clusterP);
        add(nN, SP_N, clusterN);

        t = 0.0;
        computeAcc();
    }

    void clear() { p.clear(); nextHadron = 0; t = 0.0; }

    double hadronRadiusFm = 0.8409;   // triangle radius for uud/udd; proton rms charge radius

    // astronomical units in simulator units (a0, m_e, a0/v_Bohr)
    static constexpr double AU_A0   = 2.8269900e21;
    static constexpr double KPC_A0  = 5.8311279e29;
    static constexpr double MSUN_ME = 2.1833749e60;
    static constexpr double YEAR_T  = 1.3046345e24;

    Particle body(int sp, double m, QPointF at, QPointF v = QPointF(0, 0))
    {
        Particle a; a.pos = at; a.vel = v; a.q = 0.0; a.m = m;
        a.species = sp; a.nucleon = false; return a;
    }

    // give every body except 'center' the circular speed the CURRENT force
    // law implies at its radius (counter-clockwise), so a scene starts on
    // orbits whatever route, eta_g or gain is selected
    void circularize(int center)
    {
        computeAcc();
        const QPointF c = p[center].pos;
        for (int i = 0; i < int(p.size()); ++i) {
            if (i == center) continue;
            QPointF d = p[i].pos - c; double r = std::hypot(d.x(), d.y());
            if (r <= 0) continue;
            QPointF u = d / r;
            double ar = -(p[i].acc.x() * u.x() + p[i].acc.y() * u.y());  // inward
            double v  = (ar > 0) ? std::sqrt(ar * r) : 0.0;
            p[i].vel = p[center].vel + QPointF(-u.y() * v, u.x() * v);
        }
    }

    // bulge + N tracer stars (1e3 M_sun) on a disk, 0.5 .. R_kpc
    void buildGalaxy(double Mbulge_sun, int nStars, double R_kpc, quint32 seed)
    {
        clear();
        QRandomGenerator rng(seed);
        p.push_back(body(SP_STAR, Mbulge_sun * MSUN_ME, QPointF(0, 0)));
        for (int i = 0; i < nStars; ++i) {
            double r  = KPC_A0 * (0.5 + (R_kpc - 0.5) * std::sqrt(rng.generateDouble()));
            double th = 2.0 * M_PI * rng.generateDouble();
            // 1e3 M_sun tracers: with 1e6 M_sun, tracer-tracer encounters
            // at dt = 1e30 scatter ~1/5 of the disk within 0.5 Gyr
            p.push_back(body(SP_STAR, 1e3 * MSUN_ME, QPointF(r * std::cos(th), r * std::sin(th))));
        }
        circularize(0);
    }

    // Sun + the eight planets at their mean distances, phases random
    void buildSolar(quint32 seed)
    {
        clear();
        QRandomGenerator rng(seed);
        p.push_back(body(SP_STAR, MSUN_ME, QPointF(0, 0)));
        static const double au[8] = {0.387, 0.723, 1.000, 1.524, 5.203, 9.537, 19.19, 30.07};
        static const double me[8] = {0.0553, 0.815, 1.0, 0.107, 317.8, 95.2, 14.5, 17.1};  // Earth masses
        for (int k = 0; k < 8; ++k) {
            double th = 2.0 * M_PI * rng.generateDouble();
            p.push_back(body(SP_PLANET, me[k] * SPECIES[SP_PLANET].m,
                             QPointF(au[k] * AU_A0 * std::cos(th), au[k] * AU_A0 * std::sin(th))));
        }
        circularize(0);
    }

    // append n particles of any species.  Protons land in the cluster disk on
    // the proton- and mixed-cluster layouts, neutrons on the mixed one; all
    // else goes uniformly into the box.
    void addSpecies(int sp, int n, double boxHalf, double vScale, quint32 seed)
    {
        if (sp < 0 || sp >= SP_COUNT || n <= 0) return;
        QRandomGenerator rng(seed);
        const Species &S = SPECIES[sp];
        const double mass = isQuark(sp) ? quarkMass(sp) : S.m;
        const double vs = std::min(vScale / std::sqrt(mass), V_CAP);
        const bool inCluster = (sp == SP_P && layout != Layout::Random) ||
                               (sp == SP_N && layout == Layout::MixedCluster);
        const double R = std::max(0.5 * clusterFill * crossover(), 1.8897e-5);
        for (int i = 0; i < n; ++i) {
            Particle a;
            if (inCluster) {
                const double rr = R * std::sqrt(rng.generateDouble());
                const double th = 2.0 * M_PI * rng.generateDouble();
                a.pos = QPointF(rr * std::cos(th), rr * std::sin(th));
            } else {
                a.pos = QPointF(rng.bounded(2.0 * boxHalf) - boxHalf,
                                rng.bounded(2.0 * boxHalf) - boxHalf);
            }
            a.vel = QPointF((rng.bounded(2.0) - 1.0) * vs,
                            (rng.bounded(2.0) - 1.0) * vs);
            a.q = S.q;  a.m = mass;
            a.nucleon = S.nucleon;
            a.species = sp;
            p.push_back(a);
        }
        computeAcc();
    }

    // three valence quarks on an equilateral triangle: uud (proton) or udd
    // (neutron).  Triangle radius = proton rms charge radius 0.8409 fm.
    // Charge only - nothing confines them, see the header note.
    void addHadron(bool proton, int n, double boxHalf, double vScale,
                   quint32 seed)
    {
        if (n <= 0) return;
        static constexpr double A0_OVER_FM = 52917.721090;
        const double R = hadronRadiusFm / A0_OVER_FM;    // in a0
        const int q3[3] = { SP_U, proton ? SP_U : SP_D, SP_D };

        QRandomGenerator rng(seed);
        double M = 0.0;
        for (int sp : q3) M += quarkMass(sp);
        const double vs = std::min(vScale / std::sqrt(M), V_CAP);

        for (int h = 0; h < n; ++h) {
            const QPointF c(rng.bounded(2.0 * boxHalf) - boxHalf,
                            rng.bounded(2.0 * boxHalf) - boxHalf);
            const QPointF v((rng.bounded(2.0) - 1.0) * vs,
                            (rng.bounded(2.0) - 1.0) * vs);
            const double th0 = 2.0 * M_PI * rng.generateDouble();
            const int hid = nextHadron++;
            for (int k = 0; k < 3; ++k) {
                const double th = th0 + 2.0 * M_PI * k / 3.0;
                const double rr = R * (0.95 + 0.10 * rng.generateDouble());
                Particle a;                              // +/-5% radial jitter:
                a.pos = c + QPointF(rr * std::cos(th),   // an exact triangle
                                    rr * std::sin(th));  // collapses to a point
                a.vel = v;                               // hadron moves as one
                const Species &S = SPECIES[q3[k]];
                a.q = S.q;  a.m = quarkMass(q3[k]);  a.nucleon = false;
                a.species = q3[k];  a.hadron = hid;
                p.push_back(a);
            }
        }
        computeAcc();
    }

    // hydrogen-like ion: Z protons, N neutrons, one electron at r = a0/Z
    void addIon(int Z, int N, QPointF at, bool resolveNucleons,
                QRandomGenerator &rng)
    {
        if (Z < 1) return;
        const double rOrb = 1.0 / Z;                 // a0/Z
        const double vOrb = static_cast<double>(Z);  // Z*v_Bohr  (n = 1)
        const double Mnuc = Z * M_PROTON + N * M_NEUTRON;

        // electron: circular orbit, counter-clockwise
        Particle e;
        e.pos = at + QPointF(rOrb, 0.0);
        e.vel = QPointF(0.0, vOrb);
        e.q = -1.0; e.m = M_ELECTRON; e.species = SP_E;

        // nucleus recoils so the ion's net momentum is zero
        const QPointF recoil = -e.vel * (M_ELECTRON / Mnuc);

        if (!resolveNucleons) {
            Particle nuc;
            nuc.pos = at;
            nuc.vel = recoil;
            nuc.q = static_cast<double>(Z);
            nuc.m = Mnuc;
            p.push_back(nuc);
        } else {
            // ~1 fm cluster: 1 fm = 1.8897e-5 a0
            const double R = 1.8897e-5;
            for (int i = 0; i < Z + N; ++i) {
                Particle nucleon;
                double th = 2.0 * M_PI * i / (Z + N);
                double rr = R * (0.5 + 0.5 * rng.bounded(1.0));
                nucleon.pos = at + QPointF(rr * std::cos(th), rr * std::sin(th));
                nucleon.vel = recoil;
                nucleon.q = (i < Z) ? 1.0 : 0.0;
                nucleon.m = (i < Z) ? M_PROTON : M_NEUTRON;
                nucleon.nucleon = true;
                nucleon.species = (i < Z) ? SP_P : SP_N;
                p.push_back(nucleon);
            }
        }
        p.push_back(e);
        computeAcc();
    }

    void loadIons(int which, bool resolveNucleons, quint32 seed)
    {
        QRandomGenerator rng(seed);
        p.clear();
        t = 0.0;
        // Z, N for H, He-4, Li-7, Be-9
        static const int ZN[4][2] = {{1,0},{2,2},{3,4},{4,5}};
        if (which >= 0 && which < 4) {
            addIon(ZN[which][0], ZN[which][1], QPointF(0,0), resolveNucleons, rng);
        } else {
            for (int k = 0; k < 4; ++k)
                addIon(ZN[k][0], ZN[k][1],
                       QPointF(-6.0 + 4.0 * k, 0.0), resolveNucleons, rng);
        }
        computeAcc();
    }

    // SYMMETRIC pair scale inside f.  The FT force routes used S = |q_j|, the
    // OTHER body's charge, so for |q_i| != |q_j| the two bodies saw different
    // f and F_ij != -F_ji (electron / up quark at 2 a0, Potential route:
    // -0.046875 vs -0.024691 - the pair accelerated itself).  A pair force
    // must come from one symmetric pair potential U_ij f_ij^2, so S must be
    // symmetric in i and j:
    //   Potential, Centripetal  S_ij = sqrt(|q_i q_j|)   unit-independent;
    //                           equals |q| whenever |q_i| = |q_j|
    //   Discrete crossover      S_ij = |q_i q_j| (in e)  the Bohr scaling:
    //                           shells at n^2 a0 / (|q_i q_j|), so n^2 a0/Z
    //                           for an electron around Z, as before
    //   gravity                 S_ij = sqrt(m_i m_j)
    // Pairs that contain a unit charge are unchanged on the Discrete route;
    // on Potential/Centripetal an electron around Z now crosses over at
    // sqrt(Z)/eta instead of Z/eta.  The POSITIONAL clock factor f_i used by
    // the momentum law (computeDilation) is a field at a point, not a pair
    // force, and still sums |q_j| of the sources.
    double pairScaleQ(double qi, double qj) const
    {
        const double a = std::fabs(qi * qj);
        return (route == Route::DiscreteCrossover) ? a : std::sqrt(a);
    }
    static double pairScaleM(double mi, double mj) { return std::sqrt(mi * mj); }

    // one FT kernel, used by BOTH the electric and gravitoelectric terms.
    //   C  = coupling / m_i      S = source scale inside f      et = eta
    // positive result = outward (repulsive) along j -> i
    double ftKernel(double r, double C, double S, double et) const
    {
        if (S <= 0.0 || et <= 0.0) return 0.0;
        const double d = et * r + S;

        switch (route) {
        case Route::Centripetal:
            //  C * f^2 / r^2  =  C * et^2 / (et r + S)^2
            return C * et * et / (d * d);

        case Route::Potential:
            // -C * et^2 (S - et r) / (et r + S)^3
            return -C * et * et * (S - et * r) / (d * d * d);

        case Route::DiscreteCrossover: {
            // snap the crossover to the nearest shell r_n = n^2/S
            int n = static_cast<int>(std::lround(std::sqrt(r * S)));
            if (n < 1) n = 1;
            const double etaL = S * S / double(n * n);
            const double dL   = etaL * r + S;
            return -C * etaL * etaL * (S - etaL * r) / (dL * dL * dL);
        }

        case Route::Coulomb:
        default: {
            const double r2 = r * r;
            return (r2 > 1e-24) ? C / r2 : C / 1e-24;
        }
        }
    }

    // Malfliet-Tjon V force, r in simulator length units -> simulator force
    static double strongForce(double r_sim)
    {
        static constexpr double A0_OVER_FM = 52917.721090;
        static constexpr double MEVFM_TO_SIM = 1.944690e9;
        static constexpr double V_R = 1438.7200, MU_R = 3.11;
        static constexpr double V_A =  626.8850, MU_A = 1.55;

        double r = r_sim * A0_OVER_FM;               // fm
        if (r < 1e-4) r = 1e-4;                      // guard the 1/r^2
        const double t1 = V_R * std::exp(-MU_R * r) * (MU_R / r + 1.0 / (r * r));
        const double t2 = V_A * std::exp(-MU_A * r) * (MU_A / r + 1.0 / (r * r));
        return (t1 - t2) * MEVFM_TO_SIM;             // + = repulsive
    }
    static double strongPotential(double r_sim)      // simulator energy units
    {
        static constexpr double A0_OVER_FM = 52917.721090;
        static constexpr double MEV_TO_SIM = 1.602176634e-13 / 4.3597447e-18;
        static constexpr double V_R = 1438.7200, MU_R = 3.11;
        static constexpr double V_A =  626.8850, MU_A = 1.55;
        double r = r_sim * A0_OVER_FM;
        if (r < 1e-4) r = 1e-4;
        return (V_R * std::exp(-MU_R * r) - V_A * std::exp(-MU_A * r)) / r
               * MEV_TO_SIM;
    }

    // total radial acceleration on i due to j (electric + gravitoelectric),
    // with optional floor-quantisation of the radius or of the force
    double radialAcc(double r, double qi, double qj, double mi, double mj,
                     bool nuci, bool nucj) const
    {
        if (mi == 0.0) return 0.0;
        const double D = quantStep;

        double rr = r;
        if (quant == Quant::Radius && D > 0.0) {
            rr = std::floor(r / D) * D;
            if (rr < 0.5 * D) rr = 0.5 * D;      // clamp bin 0, avoids 1/0
        }

        double a = ftKernel(rr, qi * qj / mi, pairScaleQ(qi, qj), eta);
        if (gravity)
            a += ftKernel(rr, -G_UNITS * gravGain * mj, pairScaleM(mi, mj), etaG);

        // strong force: nucleons only, charge-independent
        if (strong && nuci && nucj)
            a += strongGain * strongForce(rr) / mi;

        if (quant == Quant::Force && D > 0.0)
            a = std::floor(a / D) * D;

        return a;
    }

    // intrinsic moment per unit spin, by species
    static double momentFor(const Particle &a)
    {
        if (isQuark(a.species) && a.m > 0.0)            // Dirac q/(2m), g = 2
            return SPECIES[a.species].q * 0.5 / a.m;
        if (a.species >= 0 && a.species < SP_COUNT) return SPECIES[a.species].mu0;
        if (a.q < 0.0)              return -MU_E;          // electron
        if (a.q > 1.5) {                                   // composite nucleus
            const int Z = int(std::lround(a.q));
            if (Z == 2) return 0.0;                        // He-4
            if (Z == 3) return  3.256427 * MU_N;           // Li-7
            if (Z == 4) return -1.17749  * MU_N;           // Be-9
            return 0.0;
        }
        if (a.q > 0.0)              return MU_P;           // proton / H nucleus
        if (a.nucleon)              return MU_NEU;         // neutron
        return 0.0;
    }

    void assignSpins(quint32 seed)
    {
        QRandomGenerator rng(seed);
        int k = 0;
        for (auto &a : p) {
            int s;
            if (spinMode == 1)      s = +1;
            else if (spinMode == 2) s = (k++ % 2) ? -1 : +1;
            else                    s = rng.bounded(2) ? +1 : -1;
            a.spin = s;
            a.mu   = momentFor(a) * s;
        }
    }
    void clearSpins() { for (auto &a : p) { a.spin = 0; a.mu = 0.0; } }

    // velocity- and moment-dependent forces; these are NOT central, so they
    // are accumulated as vectors after the radial pass
    void magneticPass()
    {
        const int n = static_cast<int>(p.size());
        const double k = ALPHA2 * magGain;

        for (int i = 0; i < n; ++i) {
            Particle &a = p[i];
            const double tau = (a.q < 0.0) ? THOMAS_E : 1.0;
            double Bz = 0.0, gx = 0.0, gy = 0.0;
            double Bg = 0.0;                   // gravitomagnetic, /(G/c^2)

            for (int j = 0; j < n; ++j) {
                if (j == i) continue;
                const Particle &b = p[j];
                const double rx = a.pos.x() - b.pos.x();
                const double ry = a.pos.y() - b.pos.y();
                const double r2 = rx * rx + ry * ry;
                if (r2 < 1e-30) continue;
                const double r  = std::sqrt(r2);
                const double r3 = r2 * r, r5 = r3 * r2;

                // Within one hadron the velocity-dependent point-charge terms
                // are skipped: quarks move at ~0.94c, where instantaneous
                // Biot-Savart/Lorentz is invalid and non-reciprocal (it made
                // an isolated uud drift 3.01 fm).  Bound-quark magnetism is
                // colour-magnetic and inside the Cornell fit.  The static,
                // reciprocal spin dipole-dipole term is kept.
                const bool inside = (a.hadron >= 0 && a.hadron == b.hadron);

                // gravitomagnetic field of moving mass j at i
                if (gravMag && !inside && b.m > 0.0) {
                    double fg2 = 1.0;
                    if (route != Route::Coulomb && etaG > 0.0) {
                        const double d = etaG * r + b.m;
                        fg2 = (etaG * r / d) * (etaG * r / d);
                    }
                    Bg += b.m * (b.vel.x() * ry - b.vel.y() * rx) / r3 * fg2;
                }

                // field acting on charge i's motion (lab frame)
                if (magnetic && !inside && b.q != 0.0)
                    Bz += b.q * (b.vel.x() * ry - b.vel.y() * rx) / r3;
                if (spinOn && !inside && b.mu != 0.0)
                    Bz += -b.mu / r3;

                // gradient of the field moment i sees (its rest frame)
                if (spinOn && a.mu != 0.0) {
                    if (!inside && b.q != 0.0) {
                        const double wx = b.vel.x() - tau * a.vel.x();
                        const double wy = b.vel.y() - tau * a.vel.y();
                        const double h  = wx * ry - wy * rx;
                        gx += b.q * (-wy / r3 - 3.0 * h * rx / r5);
                        gy += b.q * ( wx / r3 - 3.0 * h * ry / r5);
                    }
                    if (b.mu != 0.0) {
                        gx += 3.0 * b.mu * rx / r5;
                        gy += 3.0 * b.mu * ry / r5;
                    }
                }
            }

            if (a.m <= 0.0) continue;
            const double ax = a.q * ( a.vel.y() * Bz) + a.mu * gx;
            const double ay = a.q * (-a.vel.x() * Bz) + a.mu * gy;
            a.acc += QPointF(ax, ay) * (k / a.m);

            if (gravMag && Bg != 0.0) {
                const double kg = gmK * G_UNITS * ALPHA2 * gravGain;
                a.acc += QPointF(-a.vel.y() * Bg, a.vel.x() * Bg) * kg;
            }
        }
    }

    // dipole-dipole potential energy, for the E readout (velocity-dependent
    // terms have no potential and are left out, so expect a tiny drift)
    double dipolePotential() const
    {
        if (!spinOn) return 0.0;
        double u = 0.0;
        const int n = static_cast<int>(p.size());
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                if (p[i].mu == 0.0 || p[j].mu == 0.0) continue;
                const double r = std::hypot(p[i].pos.x() - p[j].pos.x(),
                                            p[i].pos.y() - p[j].pos.y());
                if (r < 1e-15) continue;
                u += ALPHA2 * magGain * p[i].mu * p[j].mu / (r * r * r);
            }
        return u;
    }

    // positional FT dilation at each particle, superposed as in Eq. 34:
    //   f_i = eta / (eta + sum_j |q_j| / r_ij)      electric
    //       x eta_g / (eta_g + sum_j m_j / r_ij)    gravitational, if on
    // A single source gives eta r/(eta r + S), the f of the worksheets.
    void computeDilation()
    {
        const int n = static_cast<int>(p.size());
        std::vector<double> se(n, 0.0), sg(n, 0.0);
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                // a hadron's own quarks do not dilate each other: the hadron
                // is one clock, set only by charges OUTSIDE it
                if (p[i].hadron >= 0 && p[i].hadron == p[j].hadron) continue;
                const double r = std::hypot(p[i].pos.x() - p[j].pos.x(),
                                            p[i].pos.y() - p[j].pos.y());
                if (r < 1e-300) continue;
                se[i] += std::fabs(p[j].q) / r;  se[j] += std::fabs(p[i].q) / r;
                sg[i] += p[j].m / r;             sg[j] += p[i].m / r;
            }
        for (int i = 0; i < n; ++i) {
            double f = (eta > 0.0) ? eta / (eta + se[i]) : 1.0;
            if (gravity && etaG > 0.0) {           // gain amplifies G, so
                const double eg = etaG / std::max(1.0, gravGain);   // eta_g ~ c^2/G
                f *= eg / (eg + sg[i]);
            }
            p[i].fpos = f;
        }
        // share ONE f across each hadron (mean over its quarks).  With a
        // common f the internal Cornell/Coulomb forces change the quarks'
        // momenta by equal and opposite amounts, so an isolated hadron
        // conserves momentum to round-off (2.4e-5 vs 1.7e5 with per-quark f).
        if (nextHadron > 0) {
            std::vector<double> fs(nextHadron, 0.0);
            std::vector<int>    nc(nextHadron, 0);
            for (const auto &a : p)
                if (a.hadron >= 0 && a.hadron < nextHadron) { fs[a.hadron] += a.fpos; ++nc[a.hadron]; }
            for (auto &a : p)
                if (a.hadron >= 0 && a.hadron < nextHadron && nc[a.hadron] > 0)
                    a.fpos = fs[a.hadron] / nc[a.hadron];
        }
    }

    void computeAcc()
    {
        const int n = static_cast<int>(p.size());
        for (int i = 0; i < n; ++i) p[i].acc = QPointF(0, 0);
        if (positionalKin()) computeDilation();

        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                QPointF d = p[i].pos - p[j].pos;
                double  r = std::hypot(d.x(), d.y());
                if (r < 1e-15) continue;
                QPointF u = d / r;                    // unit vector j -> i

                double ai = radialAcc(r, p[i].q, p[j].q, p[i].m, p[j].m,
                                      p[i].nucleon, p[j].nucleon);
                double aj = radialAcc(r, p[j].q, p[i].q, p[j].m, p[i].m,
                                      p[j].nucleon, p[i].nucleon);
                if (qcd && p[i].hadron >= 0 && p[i].hadron == p[j].hadron) {
                    const double F = qcdGain * qcdForce(r);
                    ai += F / p[i].m;
                    aj += F / p[j].m;
                }
                p[i].acc += u * ai;
                p[j].acc -= u * aj;
            }
        }
        if (magnetic || spinOn || gravMag) magneticPass();
        if (hasExternal()) externalPass();
    }

    // half kick.  The force F = m a changes the momentum p = G(v) m v, and v
    // is recovered from p, so for FT and SR |v| approaches c but never reaches
    // it.  With P = |p|/(m c) the inversions are closed form:
    //   SR  beta = P / sqrt(1 + P^2)
    //   FT  P = beta/(1 - beta^2)  ->  beta = 2P / (1 + sqrt(1 + 4P^2))
    //       (rationalised so small P does not cancel)
    void kick(Particle &a, double h) const
    {
        if (kin == Kin::Newtonian) { a.vel += a.acc * h; return; }
        static constexpr double C = 137.035999084;

        if (positionalKin()) {
            // p = (m0 / f^2) g(v) v ;  the force changes p ;  invert for v
            const double f2 = std::max(1e-300, a.fpos * a.fpos);
            const double v2c2 = (a.vel.x() * a.vel.x() + a.vel.y() * a.vel.y()) / (C * C);
            double g = 1.0;
            if (kin == Kin::F2SC)  g = 1.0 / std::sqrt(std::max(1e-15, 1.0 - v2c2));
            if (kin == Kin::F2Lin) g = 1.0 + 0.5 * v2c2;
            const double M  = a.m / f2;
            const double px = M * g * a.vel.x() + h * a.m * a.acc.x();
            const double py = M * g * a.vel.y() + h * a.m * a.acc.y();
            const double pm = std::hypot(px, py);
            if (pm == 0.0) { a.vel = QPointF(0, 0); return; }
            const double P = pm / (M * C);              // = g(beta) beta
            double beta;
            if (kin == Kin::F2SC) {
                beta = P / std::sqrt(1.0 + P * P);
            } else if (kin == Kin::F2Lin) {             // beta + beta^3/2 = P
                const double d = std::sqrt(P * P + 8.0 / 27.0);
                beta = std::cbrt(P + d) + std::cbrt(P - d);
                for (int k = 0; k < 3; ++k)             // polish the small-P case
                    beta -= (0.5 * beta * beta * beta + beta - P) / (1.5 * beta * beta + 1.0);
            } else {
                beta = P;
            }
            a.vel = QPointF(px, py) * (beta * C / pm);
            return;
        }
        const double v2c2 = (a.vel.x() * a.vel.x() + a.vel.y() * a.vel.y()) / (C * C);
        double G;
        switch (kin) {
        case Kin::FTHalf: G = 1.0 / std::max(1e-15, 1.0 - 0.5 * v2c2); break;
        case Kin::FT:     G = 1.0 / std::max(1e-15, 1.0 - v2c2);       break;
        default:          G = 1.0 / std::sqrt(std::max(1e-15, 1.0 - v2c2)); break;
        }
        const double px = G * a.m * a.vel.x() + h * a.m * a.acc.x();
        const double py = G * a.m * a.vel.y() + h * a.m * a.acc.y();
        const double pm = std::hypot(px, py);
        if (pm == 0.0) { a.vel = QPointF(0, 0); return; }
        const double P = pm / (a.m * C);
        double beta;
        switch (kin) {
        case Kin::FTHalf: beta = 2.0 * P / (1.0 + std::sqrt(1.0 + 2.0 * P * P)); break;
        case Kin::FT:     beta = 2.0 * P / (1.0 + std::sqrt(1.0 + 4.0 * P * P)); break;
        default:          beta = P / std::sqrt(1.0 + P * P);                     break;
        }
        a.vel = QPointF(px, py) * (beta * C / pm);
    }

    void step(double dt)
    {
        const bool q = quantised;
        auto free = [&](const Particle &a) { return !q || a.level <= 0; };

        for (auto &a : p) {
            if (!free(a)) continue;
            kick(a, 0.5 * dt);
            a.pos += a.vel * dt;
        }
        computeAcc();
        for (auto &a : p) {
            if (!free(a)) continue;
            kick(a, 0.5 * dt);
        }
        if (q) stepQuantised(dt);
        t += dt;
    }

    // kinetic energy, the one that is CONSERVED with the momentum law above,
    // i.e. K = integral of v dp:
    //   Newtonian  m v^2 / 2
    //   SR         (gamma - 1) m c^2
    //   FT         m c^2 [ G - 1 - ln(G)/2 ],   G = 1/(1 - v^2/c^2)
    //              -> m v^2/2 + 3 m v^4/(4 c^2) at low speed
    //   FTHalf     m c^2 [ 2/D - 2 + ln D ],    D = 1 - v^2/(2c^2)
    //              -> m v^2/2 + 3 m v^4/(8 c^2): identical to SR through v^4
    // NOT m c^2 (G - 1), the energy implied by manuscript Eq. 52: that is
    // m v^2 at low speed, twice the Newtonian value, and is not the work done
    // by a force under Eq. 41's momentum.
    double kinetic() const
    {
        static constexpr double C2 = 137.035999084 * 137.035999084;
        double k = 0;
        for (const auto &a : p) {
            const double v2 = a.vel.x() * a.vel.x() + a.vel.y() * a.vel.y();
            const double b2 = std::min(1.0 - 1e-15, v2 / C2);   // FT/SR only
            switch (kin) {
            case Kin::Newtonian: k += 0.5 * a.m * v2; break;
            case Kin::SR:        k += (1.0 / std::sqrt(1.0 - b2) - 1.0) * a.m * C2; break;
            case Kin::FT: {
                const double G = 1.0 / (1.0 - b2);
                k += (G - 1.0 - 0.5 * std::log(G)) * a.m * C2;
                break;
            }
            case Kin::FTHalf: {
                const double D = std::max(1e-15, 1.0 - 0.5 * v2 / C2);
                k += (2.0 / D - 2.0 + std::log(D)) * a.m * C2;
                break;
            }
            // positional kinds: integral v dp at fixed position, over f^2.
            // NOT conserved along a non-circular orbit - see the header.
            case Kin::F2SC: {
                const double f2 = std::max(1e-300, a.fpos * a.fpos);
                k += (1.0 / std::sqrt(1.0 - b2) - 1.0) * a.m * C2 / f2;
                break;
            }
            case Kin::F2Lin: {
                const double f2 = std::max(1e-300, a.fpos * a.fpos);
                k += a.m * (0.5 * v2 + 0.375 * v2 * v2 / C2) / f2;
                break;
            }
            case Kin::F2Pos: {
                const double f2 = std::max(1e-300, a.fpos * a.fpos);
                k += 0.5 * a.m * v2 / f2;
                break;
            }
            }
        }
        return k;
    }

    double potential() const
    {
        double u = 0;
        const int n = static_cast<int>(p.size());
        for (int i = 0; i < n; ++i)
            for (int j = i + 1; j < n; ++j) {
                QPointF d = p[i].pos - p[j].pos;
                double  r = std::hypot(d.x(), d.y());
                double  S = pairScaleQ(p[i].q, p[j].q);
                if (r < 1e-15) continue;
                const double Sg = pairScaleM(p[i].m, p[j].m);
                const double Cg = -G_UNITS * gravGain * p[i].m * p[j].m;
                if (route == Route::Potential) {
                    // U f^2 = C eta^2 r / (eta r + S)^2
                    double d2 = eta * r + S;
                    u += p[i].q * p[j].q * eta * eta * r / (d2 * d2);
                    if (gravity && Sg > 0) {
                        double dg = etaG * r + Sg;
                        u += Cg * etaG * etaG * r / (dg * dg);
                    }
                } else {
                    u += p[i].q * p[j].q / r;
                    if (gravity) u += Cg / r;
                }
                if (strong && p[i].nucleon && p[j].nucleon)
                    u += strongGain * strongPotential(r);
                if (qcd && p[i].hadron >= 0 && p[i].hadron == p[j].hadron)
                    u += qcdGain * qcdPotential(r);
            }
        // uniform E and g have potentials -q E.r and -m g.r; B and B_g do no work
        for (const auto &a : p)
            u -= a.q * (extEx * a.pos.x() + extEy * a.pos.y())
               + a.m * (extGx * a.pos.x() + extGy * a.pos.y());
        return u + dipolePotential();
    }
};

// ---------------------------------------------------------------------------
// view
// ---------------------------------------------------------------------------

class View : public QWidget {
    Q_OBJECT
public:
    Sim    *sim;
    double  scale = 40.0;            // pixels per a0
    QPointF pan   = QPointF(0, 0);   // world point shown at the widget centre
    bool    showTrails = true;
    bool    showLabels = true;

    explicit View(Sim *s, QWidget *parent = nullptr) : QWidget(parent), sim(s)
    {
        setMinimumSize(320, 240);
        setAutoFillBackground(true);
        setMouseTracking(true);
        setCursor(Qt::OpenHandCursor);
    }

    QPointF centrePx() const { return QPointF(width() / 2.0, height() / 2.0); }
    QPointF toScreen(const QPointF &w) const
    {
        QPointF d = w - pan;
        return centrePx() + QPointF(d.x() * scale, -d.y() * scale);
    }
    QPointF toWorld(const QPointF &s) const
    {
        QPointF d = s - centrePx();
        return pan + QPointF(d.x() / scale, -d.y() / scale);
    }

    double homeScale = 40.0;           // what double-click returns to
    void resetView()
    {
        pan = QPointF(0, 0);
        scale = homeScale;
        emit viewChanged(scale);
        update();
    }

signals:
    void viewChanged(double newScale);

protected:
    // wheel zooms about the cursor: the world point under it stays put
    void wheelEvent(QWheelEvent *ev) override
    {
        const QPointF sp = ev->position();
        const QPointF anchor = toWorld(sp);

        const double steps = ev->angleDelta().y() / 120.0;
        if (steps == 0.0) return;
        double next = scale * std::pow(1.15, steps);
        next = qBound(1e-45, next, 1e12);
        scale = next;

        // put the anchor back under the cursor
        QPointF d = sp - centrePx();
        pan = anchor - QPointF(d.x() / scale, -d.y() / scale);

        emit viewChanged(scale);
        update();
        ev->accept();
    }

    void mousePressEvent(QMouseEvent *ev) override
    {
        if (ev->button() == Qt::LeftButton) {
            dragging = true;
            dragFrom = ev->position();
            panFrom  = pan;
            setCursor(Qt::ClosedHandCursor);
        }
    }
    void mouseMoveEvent(QMouseEvent *ev) override
    {
        if (!dragging) return;
        QPointF d = ev->position() - dragFrom;
        pan = panFrom - QPointF(d.x() / scale, -d.y() / scale);
        update();
    }
    void mouseReleaseEvent(QMouseEvent *) override
    {
        dragging = false;
        setCursor(Qt::OpenHandCursor);
    }
    void mouseDoubleClickEvent(QMouseEvent *) override { resetView(); }

private:
    bool    dragging = false;
    QPointF dragFrom, panFrom;

protected:
    void paintEvent(QPaintEvent *) override
    {
        QPainter g(this);
        g.setRenderHint(QPainter::Antialiasing, true);
        g.fillRect(rect(), QColor("#0C1A33"));

        const QPointF c = toScreen(QPointF(0, 0));   // world origin on screen

        // crossover radius S/eta for unit charge, drawn as a guide circle
        if (sim->route == Route::Potential && sim->eta > 0) {
            double rc = 1.0 / sim->eta * scale;
            if (rc > 2 && rc < std::max(width(), height())) {
                g.setPen(QPen(QColor("#C56A00"), 1, Qt::DashLine));
                g.drawEllipse(c, rc, rc);
            }
        }

        // grid
        g.setPen(QPen(QColor("#1E2D4A"), 1));
        {   // grid rings at a decade that is actually visible
            double span = std::max(width(), height()) / scale;   // a0 across
            double unit = std::pow(10.0, std::floor(std::log10(span))) / 2.0;
            for (int k = 1; k <= 12; ++k) {
                double r = k * unit * scale;
                if (r > 1.5 * std::max(width(), height())) break;
                g.drawEllipse(c, r, r);
            }
            g.setPen(QColor("#5D6675"));
            g.drawText(int(c.x()) + 4, int(c.y()) - 4,
                       QString("ring = %1 a0").arg(unit, 0, 'g', 3));
            g.setPen(QPen(QColor("#1E2D4A"), 1));
        }

        // trails
        if (showTrails) {
            for (const auto &a : sim->p) {
                if (a.trail.size() < 2) continue;
                QColor col = colorOf(a);
                col.setAlpha(90);
                g.setPen(QPen(col, 1));
                QVector<QPointF> pts;
                pts.reserve(static_cast<int>(a.trail.size()));
                for (const auto &t : a.trail) pts.append(toScreen(t));
                g.drawPolyline(pts.constData(), pts.size());
            }
        }

        // stationary-state shells (also drawn for the discrete-crossover route,
        // whose stable zeros sit on exactly the same radii)
        if (sim->quantised || sim->route == Route::DiscreteCrossover) {
            for (const auto &nuc : sim->p) {
                if (nuc.q <= 0.0) continue;
                const int nTop = (sim->route == Route::DiscreteCrossover)
                                     ? qMax(sim->nMax, 6) : sim->nMax;
                for (int n = 1; n <= nTop; ++n) {
                    double rr = Sim::r_n(n, nuc.q) * scale;
                    if (rr < 3 || rr > 3.0 * std::max(width(), height())) continue;
                    g.setBrush(Qt::NoBrush);
                    g.setPen(QPen(QColor(63, 180, 207, 70), 1, Qt::DotLine));
                    g.drawEllipse(toScreen(nuc.pos), rr, rr);
                }
            }
        }

        // particles
        for (const auto &a : sim->p) {
            g.setBrush(colorOf(a));
            g.setPen(Qt::NoPen);
            double rad = (a.species < 0)  ? 6.0          // composite nucleus
                       : (a.m > 1000.0)   ? 5.0
                       : (a.m > 100.0)    ? 4.0
                       : (a.m < 1e-3)     ? 2.0          // neutrinos
                                          : 3.0;
            if (showLabels)                                // room for the letter
                rad = (a.species < 0) ? 10.0 : std::max(rad + 4.0, 7.5);
            g.drawEllipse(toScreen(a.pos), rad, rad);
            if (showLabels) drawSymbol(g, a, toScreen(a.pos), rad);
            if (sim->spinOn && a.spin != 0) {
                // spin arrow: up = +z (drawn toward screen top), down = -z
                const QPointF c0 = toScreen(a.pos);
                const double  L  = rad + 6.0;
                const double  sgn = (a.spin > 0) ? -1.0 : 1.0;   // screen y is down
                g.setPen(QPen(QColor("#E8EDF5"), 1.2));
                g.drawLine(c0, c0 + QPointF(0, sgn * L));
                g.drawLine(c0 + QPointF(0, sgn * L), c0 + QPointF(-2.5, sgn * (L - 3.5)));
                g.drawLine(c0 + QPointF(0, sgn * L), c0 + QPointF( 2.5, sgn * (L - 3.5)));
            }
            if (sim->quantised && a.level > 0) {
                g.setPen(QColor("#9AA6B8"));
                g.drawText(toScreen(a.pos) + QPointF(rad + 3, -rad + 2),
                           QString("n=%1").arg(a.level));
            }
        }

        // hud
        g.setPen(QColor("#E8EDF5"));
        double K = sim->kinetic(), U = sim->potential();
        g.drawText(10, 72, QString("zoom %1 px/a0   centre (%2, %3)   "
                                   "[wheel = zoom at cursor, drag = pan, "
                                   "double-click = reset]")
                               .arg(scale, 0, 'g', 4)
                               .arg(pan.x(), 0, 'g', 3).arg(pan.y(), 0, 'g', 3));
        g.drawText(10, 18, QString("t = %1 (%2)   N = %3   ETA = %4")
                               .arg(sim->t, 0, 'g', 4).arg(siTime(sim->t))
                               .arg(sim->p.size())
                               .arg(sim->eta, 0, 'g', 4));
        g.drawText(10, 36, QString("KE = %1   PE = %2   E = %3")
                               .arg(K, 0, 'g', 6).arg(U, 0, 'g', 6)
                               .arg(K + U, 0, 'g', 6));

        int inAtt = 0, totPair = 0;
        sim->attractivePairs(inAtt, totPair);
        if (sim->route == Route::Potential || sim->route == Route::DiscreteCrossover)
            g.drawText(10, 54, QString("r_c = S/ETA = %1 a0   charged pairs inside: %2 / %3")
                               .arg(sim->crossover(), 0, 'g', 4)
                               .arg(inAtt).arg(totPair));

        if (sim->quantised) drawLevels(g);
        drawLegend(g);
        drawFields(g);
    }

    // level table, top-right: E_n with f^2, and the shift against pure Bohr
    void drawLevels(QPainter &g)
    {
        double Z = 1.0;
        for (const auto &a : sim->p) if (a.q > 0.0) { Z = a.q; break; }
        const double HA = 27.211386245988;      // eV per Hartree
        const double HZ = 6.579683920502e15;    // Hz per Hartree

        QStringList L;
        L << QString("stationary states, Z = %1, ETA = %2")
                 .arg(Z, 0, 'g', 3).arg(sim->eta, 0, 'g', 4)
          << QString(" n    r_n (a0)     E_n (eV)        shift (eV)");
        for (int n = 1; n <= sim->nMax; ++n)
            L << QString(" %1  %2  %3  %4")
                     .arg(n)
                     .arg(Sim::r_n(n, Z), 10, 'f', 6)
                     .arg(sim->E_n(n, Z) * HA, 14, 'f', 8)
                     .arg((sim->E_n(n, Z) - sim->E_bohr(n, Z)) * HA, 12, 'e', 4);
        if (sim->nMax >= 2) {
            double d = (sim->E_n(2, Z) - sim->E_bohr(2, Z))
                     - (sim->E_n(1, Z) - sim->E_bohr(1, Z));
            L << QString("1s-2s FT shift: %1 Hz   (measured to ~10 Hz)")
                     .arg(std::fabs(d) * HZ, 0, 'e', 4);
        }

        QFontMetrics fm(g.font());
        int w = 0; for (const auto &t : L) w = qMax(w, fm.horizontalAdvance(t));
        const int pad = 8, lh = 16;
        const int bw = w + 2 * pad, bh = L.size() * lh + 2 * pad;
        const int x0 = width() - bw - 10, y0 = 10;
        g.setPen(QPen(QColor("#1E2D4A"), 1));
        g.setBrush(QColor(12, 26, 51, 215));
        g.drawRoundedRect(QRectF(x0, y0, bw, bh), 4, 4);
        g.setPen(QColor("#E8EDF5"));
        for (int i = 0; i < L.size(); ++i)
            g.drawText(x0 + pad, y0 + pad + lh * (i + 1) - 4, L[i]);
    }

    // ---- ruler helpers (drawn inside the legend) --------------------------
    // largest 1-2-5 x 10^k that fits under the target pixel length
    static double niceBelow(double L)
    {
        const double e = std::floor(std::log10(L));
        const double b = std::pow(10.0, e);
        const double m = L / b;
        return ((m >= 5.0) ? 5.0 : (m >= 2.0) ? 2.0 : 1.0) * b;
    }
    static QString siTime(double tsim)
    {
        const double s = tsim * 2.4188843e-17;
        if (s >= 3.15576e7) return QString("%1 yr").arg(s / 3.15576e7, 0, 'g', 4);
        if (s >= 86400)     return QString("%1 d").arg(s / 86400, 0, 'g', 4);
        if (s >= 1.0)       return QString("%1 s").arg(s, 0, 'g', 4);
        if (s >= 1e-9)      return QString("%1 ns").arg(s * 1e9, 0, 'g', 4);
        if (s >= 1e-15)     return QString("%1 fs").arg(s * 1e15, 0, 'g', 4);
        return QString("%1 as").arg(s * 1e18, 0, 'g', 4);
    }
    static QString siLabel(double metres)
    {
        struct P { double f; const char *s; };
        static const P pre[] = {
            {3.0857e19, "kpc"}, {3.0857e16, "pc"}, {9.4607e15, "ly"},
            {1.495979e11, "AU"}, {1e3, "km"},
            {1.0, "m"}, {1e-3, "mm"}, {1e-6, "\u00B5m"}, {1e-9, "nm"},
            {1e-10, "\u00C5"}, {1e-12, "pm"}, {1e-15, "fm"},
            {1e-18, "am"}, {1e-21, "zm"}, {1e-24, "ym"}};
        for (const auto &p : pre)
            if (metres >= p.f * 0.999)
                return QString("%1 %2").arg(metres / p.f, 0, 'g', 4)
                                       .arg(QString::fromUtf8(p.s));
        return QString("%1 m").arg(metres, 0, 'e', 3);
    }

    // particle symbol centred on its disc; text colour picked for contrast
    void drawSymbol(QPainter &g, const Particle &a, const QPointF &c, double rad)
    {
        g.save();                                      // font/pen restored below
        QString sym, sub;
        if (a.species >= 0 && a.species < SP_COUNT) {
            sym = QString::fromUtf8(SYM[a.species]);
            sub = QString::fromUtf8(SUB[a.species]);
        } else {                                       // composite nucleus by Z
            static const char *EL[] = {"?", "H", "He", "Li", "Be", "B", "C",
                                       "N", "O", "F", "Ne"};
            const int Z = int(std::lround(a.q));
            sym = (Z >= 1 && Z <= 10) ? EL[Z] : QString("Z%1").arg(Z);
        }

        const QColor fill = colorOf(a);
        const double lum = 0.299 * fill.redF() + 0.587 * fill.greenF()
                         + 0.114 * fill.blueF();
        const QColor ink = (lum > 0.55) ? QColor("#0C1A33") : QColor("#FFFFFF");

        QFont f = g.font();
        f.setBold(true);
        f.setPixelSize(int(std::max(7.0, rad * (sym.size() > 1 ? 0.95 : 1.25))));
        g.setFont(f);
        g.setPen(ink);

        if (sub.isEmpty()) {
            g.drawText(QRectF(c.x() - rad, c.y() - rad, 2 * rad, 2 * rad),
                       Qt::AlignCenter, sym);
        } else {
            // main symbol shifted left, subscript small and lower right
            g.drawText(QRectF(c.x() - rad - 2, c.y() - rad - 1, 2 * rad, 2 * rad),
                       Qt::AlignCenter, sym);
            QFont fs = f;
            fs.setPixelSize(std::max(6, int(f.pixelSize() * 0.62)));
            g.setFont(fs);
            g.drawText(QPointF(c.x() + rad * 0.15, c.y() + rad * 0.62), sub);
        }
        g.restore();
    }

    static QColor colorOf(const Particle &a)
    {
        if (a.species >= 0 && a.species < SP_COUNT)
            return QColor(SPECIES[a.species].color);
        return speciesColor(a.q);
    }
    static QColor speciesColor(double q)
    {
        if (q > 1.5) return QColor("#E08A2E");  // composite nucleus, Z >= 2
        if (q > 0)   return QColor("#C56A00");  // proton / H nucleus
        if (q < 0)   return QColor("#3FB4CF");  // electron
        return QColor("#8892A4");               // neutron
    }

    // ---- legend, bottom-left -------------------------------------------
    // compass for the external fields, top-right under the level table
    void drawFields(QPainter &g)
    {
        if (!sim->hasExternal()) return;
        struct F { const char *n; QColor c; double x, y; bool axial; };
        const F fs[4] = {
            {"E",   QColor("#3FB4CF"), sim->extEx, sim->extEy, false},
            {"B",   QColor("#B48CFF"), 0, sim->extBz, true},
            {"g",   QColor("#7FD4A0"), sim->extGx, sim->extGy, false},
            {"B_g", QColor("#E0A050"), 0, sim->extBgz, true}};
        const int R = 18, W = 4 * (2 * R + 14) + 10, H = 2 * R + 34;
        const int x0 = width() - W - 10, y0 = height() - H - 10;
        g.save();
        g.setPen(QPen(QColor("#1E2D4A"), 1));
        g.setBrush(QColor(12, 26, 51, 215));
        g.drawRoundedRect(QRectF(x0, y0, W, H), 4, 4);
        for (int k = 0; k < 4; ++k) {
            const QPointF c(x0 + 10 + R + k * (2 * R + 14), y0 + 10 + R);
            g.setBrush(Qt::NoBrush);
            g.setPen(QPen(QColor("#2A3B5A"), 1));
            g.drawEllipse(c, R, R);
            g.setPen(QPen(fs[k].c, 2));
            if (!fs[k].axial) {
                const double m = std::hypot(fs[k].x, fs[k].y);
                if (m > 0) {
                    const QPointF tip = c + QPointF(fs[k].x / m * (R - 3), -fs[k].y / m * (R - 3));
                    g.drawLine(c, tip);
                    const QPointF u = (tip - c) / (R - 3), n(-u.y(), u.x());
                    g.drawLine(tip, tip - u * 6 + n * 4);
                    g.drawLine(tip, tip - u * 6 - n * 4);
                }
            } else if (fs[k].y > 0) {           // out of the screen: dot
                g.setBrush(fs[k].c); g.drawEllipse(c, 3, 3);
            } else if (fs[k].y < 0) {           // into the screen: cross
                g.drawLine(c + QPointF(-6, -6), c + QPointF(6, 6));
                g.drawLine(c + QPointF(-6, 6), c + QPointF(6, -6));
            }
            g.setPen(QColor("#E8EDF5"));
            g.drawText(QRectF(c.x() - R - 6, c.y() + R + 2, 2 * R + 12, 14),
                       Qt::AlignHCenter, fs[k].n);
        }
        g.restore();
    }

    void drawLegend(QPainter &g)
    {
        int cnt[SP_COUNT] = {0};
        double mShown[SP_COUNT] = {0};
        int nNuc = 0;
        for (const auto &a : sim->p) {
            if (a.species >= 0 && a.species < SP_COUNT) {
                ++cnt[a.species]; mShown[a.species] = a.m;
            } else ++nNuc;
        }

        struct Row { QColor col; int shape; QString text; };
        //  shape: 0 = filled dot, 1 = dashed circle, 2 = solid circle
        QVector<Row> rows;
        if (nNuc)
            rows.append({QColor("#E08A2E"), 0,
                         QString("nucleus   q = +Z,  composite        (%1)").arg(nNuc)});
        for (int k = 0; k < SP_COUNT; ++k) {
            if (!cnt[k]) continue;
            const Species &S = SPECIES[k];
            rows.append({QColor(S.color), 0,
                         QString("%1 q = %2  m = %3   (%4)")
                             .arg(QString(S.name), -9)
                             .arg(QString(S.qText), -4)
                             .arg(mShown[k], -11, 'g', 7)
                             .arg(cnt[k])});
        }
        if (sim->route == Route::Potential && sim->eta > 0)
            rows.append({QColor("#C56A00"), 1,
                         QString("crossover  r = S/ETA = %1 a0")
                             .arg(1.0 / sim->eta, 0, 'g', 4)});
        rows.append({QColor("#1E2D4A"), 2, QString("grid rings, 1 a0 apart")});
        if (sim->gravMag)
            rows.append({QColor("#7FD4A0"), -1,
                         QString("+ gravitomagnetic, k = %1  x gravity gain %2")
                             .arg(sim->gmK, 0, 'g', 3)
                             .arg(sim->gravGain, 0, 'g', 3)});
        if (sim->magnetic)
            rows.append({QColor("#B48CFF"), -1,
                         QString("+ moving-charge field (Biot-Savart)  x gain %1")
                             .arg(sim->magGain, 0, 'g', 3)});
        if (sim->spinOn) {
            int up = 0, dn = 0;
            for (const auto &a : sim->p) { if (a.spin > 0) ++up; else if (a.spin < 0) ++dn; }
            rows.append({QColor("#B48CFF"), -1,
                         QString("+ spin moments, g_e = 2.0023   (up %1 / down %2)")
                             .arg(up).arg(dn)});
        }
        if (sim->positionalKin() && sim->route != Route::Coulomb)
            rows.append({QColor("#E05050"), -1,
                         QString("! f^2 in momentum AND force route - counted twice")});
        if (sim->qcd || sim->constituent)
            rows.append({QColor("#F2C14E"), -1,
                         QString("+ QCD: %1%2")
                             .arg(sim->qcd ? "Cornell confinement within hadrons" : "")
                             .arg(sim->constituent
                                  ? (sim->qcd ? ", constituent masses"
                                              : "constituent quark masses")
                                  : "")});
        if (sim->strong)
            rows.append({QColor("#D46A8A"), -1,
                         QString("+ strong (Malfliet-Tjon V)  x gain %1")
                             .arg(sim->strongGain, 0, 'g', 3)});
        if (sim->gravity)
            rows.append({QColor("#7FD4A0"), -1,
                         QString("+ gravity  G = 2.4006e-43  x gain %1   ETA_G = %2")
                             .arg(sim->gravGain, 0, 'g', 3)
                             .arg(sim->etaG, 0, 'g', 4)});

        QString routeName =
            (sim->route == Route::Potential)   ? "force  -d(U f^2)/dr   [crossover]"
          : (sim->route == Route::Centripetal) ? "force  v^2 / r        [no crossover]"
                                               : "force  Coulomb, f = 1 [reference]";
        rows.append({QColor("#E8EDF5"), -1, routeName});

        const int pad = 8, lh = 18, sw = 22;
        QFontMetrics fm(g.font());
        int wMax = 0;
        for (const auto &r : rows) wMax = qMax(wMax, fm.horizontalAdvance(r.text));

        // ruler: largest 1-2-5 x 10^k fitting in 150 px, labelled a0 + SI
        const double A0   = 5.29177210903e-11;
        const double rL   = niceBelow(150.0 / scale);          // a0
        const double rPx  = rL * scale;                         // px
        const double mant = rL / std::pow(10.0, std::floor(std::log10(rL)));
        const int    nSub = (mant > 4.5) ? 5 : (mant > 1.5) ? 4 : 5;
        const QString rTxt = QString("%1 a0  =  %2")
                                 .arg(rL, 0, 'g', 3).arg(siLabel(rL * A0));
        const int rulerH = 30;

        const int innerW = qMax(sw + wMax + pad,
                                qMax(int(std::ceil(rPx)),
                                     fm.horizontalAdvance(rTxt)));
        const int boxW = innerW + 2 * pad;
        const int boxH = rows.size() * lh + rulerH + 2 * pad;
        const int x0 = 10, y0 = height() - boxH - 10;

        g.setPen(QPen(QColor("#1E2D4A"), 1));
        g.setBrush(QColor(12, 26, 51, 210));
        g.drawRoundedRect(QRectF(x0, y0, boxW, boxH), 4, 4);

        for (int i = 0; i < rows.size(); ++i) {
            const Row &r = rows[i];
            const QPointF sym(x0 + pad + sw / 2.0, y0 + pad + lh * (i + 0.5));
            switch (r.shape) {
            case 0:
                g.setBrush(r.col); g.setPen(Qt::NoPen);
                g.drawEllipse(sym, 4, 4);
                break;
            case 1:
                g.setBrush(Qt::NoBrush);
                g.setPen(QPen(r.col, 1, Qt::DashLine));
                g.drawEllipse(sym, 6, 6);
                break;
            case 2:
                g.setBrush(Qt::NoBrush);
                g.setPen(QPen(r.col, 1));
                g.drawEllipse(sym, 6, 6);
                break;
            default:
                break;                       // text-only row
            }
            g.setPen(QColor("#E8EDF5"));
            g.drawText(QPointF(x0 + pad + sw + pad,
                               y0 + pad + lh * (i + 1) - 5), r.text);
        }

        // ---- ruler row ----
        const double ry = y0 + pad + rows.size() * lh + 8;   // bar baseline
        const double bx = x0 + pad;

        g.setPen(QPen(QColor("#1E2D4A"), 1));
        g.drawLine(QPointF(x0 + pad, ry - 6), QPointF(x0 + boxW - pad, ry - 6));

        g.setPen(QPen(QColor("#E8EDF5"), 2));
        g.drawLine(QPointF(bx, ry + 4),       QPointF(bx + rPx, ry + 4));
        g.drawLine(QPointF(bx, ry - 1),       QPointF(bx, ry + 9));
        g.drawLine(QPointF(bx + rPx, ry - 1), QPointF(bx + rPx, ry + 9));
        g.setPen(QPen(QColor("#9AA6B8"), 1));
        for (int k = 1; k < nSub; ++k) {
            const double x = bx + rPx * k / nSub;
            g.drawLine(QPointF(x, ry + 1), QPointF(x, ry + 7));
        }
        g.setPen(QColor("#E8EDF5"));
        g.drawText(QPointF(bx, ry + 24), rTxt);
    }
};

// ---------------------------------------------------------------------------
// main window
// ---------------------------------------------------------------------------

// QDoubleSpinBox that shows and accepts scientific notation, and steps
// logarithmically - needed for values like ETA = 3.8327e25 or ETA_G = 7.8e46,
// which a fixed-decimal spin box cannot even display
class SciSpinBox : public QDoubleSpinBox {
public:
    using QDoubleSpinBox::QDoubleSpinBox;
    QString textFromValue(double v) const override
    { return QString::number(v, 'g', 7); }
    double valueFromText(const QString &t) const override
    { return t.trimmed().toDouble(); }
    QValidator::State validate(QString &t, int &) const override
    {
        bool ok = false;
        const double v = t.trimmed().toDouble(&ok);
        if (ok && v >= minimum() && v <= maximum()) return QValidator::Acceptable;
        static const QRegularExpression partial(
            "^\\s*[+-]?[0-9]*\\.?[0-9]*([eE][+-]?[0-9]*)?\\s*$");
        return partial.match(t).hasMatch() ? QValidator::Intermediate
                                           : QValidator::Invalid;
    }
    void stepBy(int steps) override               // x10^(0.1) per step
    { setValue(qBound(minimum(), value() * std::pow(10.0, 0.1 * steps), maximum())); }
};

class Window : public QWidget {
    Q_OBJECT

    // add a form row whose label AND field share one tooltip, so hovering
    // either side explains the control
    // wrap in <pre> so the alignment and line breaks survive Qt's rich-text
    // detection; escape the few characters that would be read as markup
    static QString tipHtml(const QString &s)
    {
        QString e = s.toHtmlEscaped();
        return QStringLiteral("<pre style='font-family:monospace'>%1</pre>").arg(e);
    }
    static void addRow(QFormLayout *f, const QString &text,
                       QWidget *field, const QString &tip)
    {
        auto *lab = new QLabel(text);
        const QString h = tipHtml(tip);
        lab->setToolTip(h);
        field->setToolTip(h);
        f->addRow(lab, field);
    }
    static void addRow(QFormLayout *f, QWidget *field, const QString &tip)
    {
        field->setToolTip(tipHtml(tip));
        f->addRow(field);
    }

public:
    Sim    sim;

    // one frozen simulation + its settings per scale tab; switching tabs
    // parks the current one and resumes the other where it was left
    struct TabState {
        bool    init = false;
        Sim     sim;
        double  eta = 1, etaG = 7.8e46, dt = 2e-4, box = 5, scale = 40, home = 40;
        QPointF pan;
        int     etaPreset = 4, route = 0, kin = 3, layout = 2;
        bool    qcd = false, constituent = false;
    } tabState[4];
    int curTab = -1;
    View  *view;
    QTimer timer;

    QSpinBox       *sbSpecN;
    QPushButton    *btClear, *btRestart;
    QSlider        *slField[4];
    QDial          *dlField[4];
    QLabel         *lbField[4];
    QDoubleSpinBox *sbHadR;
    QComboBox      *cbSpecies;
    QPushButton    *btAddSpecies;
    QDoubleSpinBox *sbEta, *sbDt, *sbBox, *sbV, *sbEtaG, *sbGain, *sbFill;
    QComboBox      *cbEtaPreset;
    QTabBar        *tabs;
    QComboBox      *cbRoute, *cbIon, *cbLayout, *cbQuantMode;
    QDoubleSpinBox *sbStep;
    QCheckBox      *chTrails, *chGravity, *chResolve, *chQuant, *chStrong;
    QCheckBox      *chMagnetic, *chSpin, *chGravMag, *chLabels;
    QCheckBox      *chQCD, *chConstituent;
    QComboBox      *cbKin;
    QDoubleSpinBox *sbQcdGain;
    QComboBox      *cbGmK;
    QComboBox      *cbSpinMode;
    QDoubleSpinBox *sbMagGain;
    QDoubleSpinBox *sbSGain;
    QSpinBox       *sbNmax;
    QPushButton    *btIon;
    QSlider        *slZoom;
    QPushButton    *btReset, *btPause;
    int             stepsPerFrame = 4;
    quint32         seed = 12345;

    Window()
    {
        view = new View(&sim, this);


        cbSpecies = new QComboBox;
        for (int k = 0; k < SP_COUNT; ++k)
            cbSpecies->addItem(QString("%1   q = %2   m = %3 m_e")
                                   .arg(SPECIES[k].name)
                                   .arg(SPECIES[k].qText)
                                   .arg(SPECIES[k].m, 0, 'g', 7));
        // group entries, appended after the individual species
        cbSpecies->insertSeparator(cbSpecies->count());
        cbSpecies->addItem("all leptons   e, mu, tau, nu_e, nu_mu, nu_tau", 1000);
        cbSpecies->addItem("all quarks    u, d, s, c, b, t",                1001);
        cbSpecies->addItem("all leptons + all quarks",                      1002);
        cbSpecies->insertSeparator(cbSpecies->count());
        cbSpecies->addItem("proton   u + u + d   (q = +1)",                 1003);
        cbSpecies->addItem("neutron  u + d + d   (q =  0)",                 1004);
        cbSpecies->setCurrentIndex(SP_E);
        sbSpecN = new QSpinBox; sbSpecN->setRange(1, 400); sbSpecN->setValue(5);
        btAddSpecies = new QPushButton("Add to scene");
        sbEta = new SciSpinBox; sbEta->setDecimals(300);
        sbEta->setRange(1e-6, 1e40); sbEta->setValue(1.0);

        // presets, in simulator units ETA = eta * a0 / e   (eta in C/m)
        cbEtaPreset = new QComboBox;
        cbEtaPreset->addItem("lambda_e = c^2/sqrt(KG)  manuscript   3.8327e25", 3.832733e25);
        cbEtaPreset->addItem("retrofit: H e/mu proton radii       1.044e12",  1.043742e12);
        cbEtaPreset->addItem("spectroscopy 95% floor             6.606e11",  6.606457e11);
        cbEtaPreset->addItem("nuclear: r_c = 1 fm                5.2918e4",  52917.72109);
        cbEtaPreset->addItem("visible on screen: r_c = 1 a0      1",         1.0);
        cbEtaPreset->addItem("custom", -1.0);
        cbEtaPreset->setCurrentIndex(4);                 // matches sbEta = 1
        sbEta->setSingleStep(0.1);
        sbDt = new SciSpinBox; sbDt->setDecimals(300);
        sbDt->setRange(1e-15, 1e36);  sbDt->setValue(2e-4);
        sbBox = new SciSpinBox; sbBox->setDecimals(300);
        sbBox->setRange(1e-6, 1e36); sbBox->setValue(5.0);
        sbV   = new QDoubleSpinBox; sbV->setRange(0.0, 5.0);  sbV->setValue(0.3);
        sbV->setSingleStep(0.05);

        cbLayout = new QComboBox;
        cbLayout->addItems({
            "random  - uniform in the box",
            "proton cluster  - protons inside r < S/ETA",
            "mixed cluster   - protons + neutrons inside r < S/ETA"});
        cbLayout->setCurrentIndex(2);                 // default: mixed cluster
        sbFill = new QDoubleSpinBox; sbFill->setDecimals(3);
        sbFill->setRange(0.05, 2.0); sbFill->setValue(0.90);
        sbFill->setSingleStep(0.05);

        cbQuantMode = new QComboBox;
        cbQuantMode->addItems({"none",
                               "floor the RADIUS   a(floor(r/D)*D)",
                               "floor the FORCE    floor(a/D)*D"});
        sbStep = new QDoubleSpinBox; sbStep->setDecimals(6);
        sbStep->setRange(1e-6, 100.0); sbStep->setValue(0.25);
        sbStep->setSingleStep(0.05);

        cbRoute = new QComboBox;
        cbRoute->addItems({"Potential  -d(U f^2)/dr      [one crossover]",
                           "Centripetal  v^2/r           [no crossover]",
                           "Coulomb  (f = 1)             [reference]",
                           "Discrete crossover           [zero at every shell]"});
        cbRoute->setCurrentIndex(0);       // default: Potential

        sbEtaG = new SciSpinBox; sbEtaG->setDecimals(300);
        sbEtaG->setRange(1e-6, 1e60); sbEtaG->setValue(7.822540e46);
        sbGain = new SciSpinBox; sbGain->setDecimals(300);
        sbGain->setRange(1.0, 1e50);  sbGain->setValue(1.0);
        chMagnetic = new QCheckBox("moving-charge magnetic field");
        chMagnetic->setChecked(true);
        chSpin     = new QCheckBox("spin magnetic moments");
        cbSpinMode = new QComboBox;
        cbSpinMode->addItems({"random", "all up", "alternating"});
        sbMagGain  = new QDoubleSpinBox; sbMagGain->setDecimals(4);
        sbMagGain->setRange(1.0, 1e6); sbMagGain->setValue(1.0);

        chQCD = new QCheckBox("QCD confinement (Cornell, within hadrons)");
        chConstituent = new QCheckBox("constituent quark masses");
        cbKin = new QComboBox;
        cbKin->addItems({"FT f^2 x (m0 + K/c^2)   p = gamma m0 v / f^2",
                         "FT f^2 x first order    p = m0 v (1+v^2/2c^2) / f^2",
                         "FT f^2 positional only  p = m0 v / f^2",
                         "FT 1/2                  p = m v / (1 - v^2/2c^2)",
                         "FT Eq.41                p = m v / (1 - v^2/c^2)",
                         "SR                      p = m v / sqrt(1 - v^2/c^2)",
                         "Newtonian               p = m v"});
        cbKin->setCurrentIndex(3);              // default: FT 1/2
        sbQcdGain = new QDoubleSpinBox; sbQcdGain->setDecimals(6);
        sbQcdGain->setRange(1e-9, 10.0); sbQcdGain->setValue(1.0);

        chStrong = new QCheckBox("strong nuclear force (Malfliet-Tjon V)");
        sbSGain  = new QDoubleSpinBox; sbSGain->setDecimals(6);
        sbSGain->setRange(1e-9, 10.0); sbSGain->setValue(1.0);

        chGravity = new QCheckBox("gravitoelectric term");
        chGravity->setChecked(true);
        chGravMag = new QCheckBox("gravitomagnetic field (moving masses)");
        chGravMag->setChecked(true);
        cbGmK = new QComboBox;
        cbGmK->addItems({"k = 4     linearised GR",
                         "k = 1.5   manuscript Eq.19, mu_g = 6 pi G/c^2",
                         "k = 1     naive EM analogy, mu_g = 4 pi G/c^2"});

        chTrails = new QCheckBox("trails");  chTrails->setChecked(true);
        chLabels = new QCheckBox("particle symbols"); chLabels->setChecked(true);
        chQuant  = new QCheckBox("discrete stationary states");
        sbNmax   = new QSpinBox; sbNmax->setRange(1, 8); sbNmax->setValue(3);

        cbIon = new QComboBox;
        cbIon->addItems({"H     Z=1 N=0   r=1.0000 a0  v=1",
                         "He+   Z=2 N=2   r=0.5000 a0  v=2",
                         "Li2+  Z=3 N=4   r=0.3333 a0  v=3",
                         "Be3+  Z=4 N=5   r=0.2500 a0  v=4",
                         "all four, spaced 4 a0"});
        cbIon->setCurrentIndex(4);
        chResolve = new QCheckBox("resolve nucleons (unbound)");
        sbHadR = new SciSpinBox; sbHadR->setDecimals(300);
        sbHadR->setRange(1e-12, 100.0); sbHadR->setValue(0.8409);
        btIon = new QPushButton("Load ion(s)");

        // logarithmic: scale = 10^(v/100 - 1)  ->  v in [0,1300] = 0.1 .. 1e12
        slZoom = new QSlider(Qt::Horizontal);
        slZoom->setRange(-4400, 1300);   // 1e-45 .. 1e12 px per a0
        slZoom->setValue(int(100.0 * (std::log10(40.0) + 1.0)));

        btReset = new QPushButton("Reset");
        btReset->setToolTip(tipHtml(
            "Restore the current tab's defaults (eta presets, dt, box,\n"
            "zoom) and rebuild its scene. Use Restart to keep settings."));
        btRestart = new QPushButton("Restart");
        btRestart->setToolTip(tipHtml(
            "Rebuild the current tab's scene with a new random seed and\n"
            "keep every setting as it is - presets, dt, zoom, route."));
        btClear = new QPushButton("Clear");
        btClear->setToolTip(tipHtml(
            "Empty the scene, then build it from the 'add species' list."));
        btPause = new QPushButton("Pause");
        btPause->setToolTip(tipHtml(
            "Freeze the integrator. The view stays live, so you\n"
            "can still zoom and pan while paused."));

        auto *form = new QFormLayout;

        // ---- external uniform fields: magnitude slider + angle dial each ----
        // slider s in [-600, 40] -> magnitude 10^(s/10) in atomic units,
        // s = -600 means off.  Spans 1e-60 .. 1e4 so every tab has range.
        auto *fieldBox  = new QGroupBox("External uniform fields");
        auto *fieldGrid = new QGridLayout(fieldBox);
        static const char *fName[4] = {"E", "B", "g", "B_g"};
        static const char *fTip[4] = {
            "Uniform ELECTRIC field acting on every charge: a = q E / m.\n"
            "Magnitude in atomic units (1 = 5.142e11 V/m); log slider\n"
            "1e-60 .. 1e4, far left = off. Dial: in-plane direction,\n"
            "0 deg = +x, counter-clockwise. Potential -q E.r is in E.",
            "Uniform MAGNETIC field: a = q (v x B) / m.\n"
            "Magnitude in atomic units (1 = 2.3505e5 T). The motion is\n"
            "planar, so only B_z acts: dial = tilt from +z, B_z = |B| cos.\n"
            "0 deg = out of the screen, 180 = into it, 90/270 = in-plane\n"
            "(no in-plane effect). Charges gyrate at omega = q B_z / m.",
            "Uniform GRAVITOELECTRIC field g, acting on every mass\n"
            "including neutrons and neutrinos: a = g.\n"
            "Magnitude in atomic units (1 = 9.0e22 m/s^2; Earth's 9.81\n"
            "m/s^2 is 1.09e-22). Dial: in-plane direction from +x.",
            "Uniform GRAVITOMAGNETIC field: a = v x B_g, on every mass.\n"
            "Magnitude in atomic units (1 = 4.134e16 s^-1). Equivalent to\n"
            "a frame rotating at Omega = B_g/2 (Coriolis); Earth's frame\n"
            "dragging at the surface is ~1e-14 s^-1 = 2.4e-31 a.u.\n"
            "Dial = tilt from +z, as for B."};
        for (int k = 0; k < 4; ++k) {
            auto *nm = new QLabel(QString("<b>%1</b>").arg(fName[k]));
            slField[k] = new QSlider(Qt::Horizontal);
            slField[k]->setRange(-600, 40);  slField[k]->setValue(-600);
            dlField[k] = new QDial;
            dlField[k]->setRange(0, 359);    dlField[k]->setWrapping(true);
            dlField[k]->setNotchesVisible(true); dlField[k]->setNotchTarget(15);
            dlField[k]->setFixedSize(40, 40);
            lbField[k] = new QLabel("off");
            lbField[k]->setMinimumWidth(150);
            const QString tip = tipHtml(fTip[k]);
            for (QWidget *w : {(QWidget*)nm, (QWidget*)slField[k],
                               (QWidget*)dlField[k], (QWidget*)lbField[k]})
                w->setToolTip(tip);
            fieldGrid->addWidget(nm,          2 * k, 0, 2, 1);
            fieldGrid->addWidget(slField[k],  2 * k, 1);
            fieldGrid->addWidget(lbField[k],  2 * k + 1, 1);
            fieldGrid->addWidget(dlField[k],  2 * k, 2, 2, 1);
            connect(slField[k], &QSlider::valueChanged, this, [this] { applyFields(); });
            connect(dlField[k], &QDial::valueChanged,   this, [this] { applyFields(); });
        }
        auto *btFieldsOff = new QPushButton("Fields off");
        fieldGrid->addWidget(btFieldsOff, 8, 0, 1, 3);
        connect(btFieldsOff, &QPushButton::clicked, this, [this] {
            for (auto *sl : slField) { sl->blockSignals(true); sl->setValue(-600); sl->blockSignals(false); }
            applyFields();
        });
        form->addRow(fieldBox);

        addRow(form, "add species", cbSpecies,
               "Every particle type is added from this list.\n\n"
               "electron  q = -1 e   m = 1 m_e\n"
               "proton    q = +1 e   m = 1836.152673 m_e\n"
               "neutron   q =  0     m = 1838.683662 m_e - electrically\n"
               "          inert; with gravity off it drifts in a straight\n"
               "          line, which is correct, not a bug\n\n"
               "Protons go inside the cluster disk when the layout is a\n"
               "proton or mixed cluster, neutrons when it is mixed; every\n"
               "other species, and nucleons on the random layout, go\n"
               "uniformly into the box.\n\n"
               "Quarks and the other leptons, ELECTRIC CHARGE ONLY.\n\n"
               "They feel the FT electric kernel (S = |q|, so fractional\n"
               "charges get fractional crossovers), the magnetic and spin\n"
               "terms, and gravity. They are not nucleons, so the strong\n"
               "force skips them. No colour, no confinement, no weak force:\n"
               "free quarks are never observed - this is a charge-only toy.\n\n"
               "species   q      m (m_e)      note\n"
               "muon     -1      206.768283   g_mu/2 = 1.00116592\n"
               "tau      -1      3477.23\n"
               "nu x3     0      1.957e-7     0.1 eV PLACEHOLDER, only\n"
               "                              upper bounds are known\n"
               "up      +2/3     4.227015     PDG 2022 current-quark\n"
               "down    -1/3     9.138962     masses (MS-bar); constituent\n"
               "strange -1/3     182.7792     u/d masses are ~336 MeV,\n"
               "charm   +2/3     2485.328     i.e. ~658 m_e\n"
               "bottom  -1/3     8180.056\n"
               "top     +2/3     337945.9\n\n"
               "Initial speeds are capped at 0.5 c: neutrinos would otherwise\n"
               "start superluminal in a non-relativistic integrator.");
        addRow(form, "count", sbSpecN,
               "How many of the selected species to add, uniform in the box.\n"
               "For the group entries this count applies to EACH species:\n"
               "count 5 with 'all leptons + all quarks' adds 60 particles.\n\n"
               "For proton / neutron it is the number of HADRONS: each one is\n"
               "three quarks on a triangle of radius 0.8409 fm (the proton rms\n"
               "charge radius), moving together.\n\n"
               "  proton  uud  q = 2/3 + 2/3 - 1/3 = +1\n"
               "  neutron udd  q = 2/3 - 1/3 - 1/3 =  0\n\n"
               "With charge only, nothing holds them: the u-u pair repels,\n"
               "each u-d pair attracts, and the triplet does not stay a\n"
               "hadron. Turn on 'QCD confinement' to bind them, and\n"
               "'constituent quark masses' to give uud 1012 MeV instead of\n"
               "the 8.99 MeV of current masses (proton 938.27 MeV).\n"
               "Zoom to ~1e6 px/a0 to see them.");
        addRow(form, btAddSpecies,
               "Append to the current scene without clearing it.\n\n"
               "With stationary states on, muons and taus bind too, using the\n"
               "reduced mass: muonic hydrogen (mu = 185.84 m_e) sits 185.84x\n"
               "closer than ordinary hydrogen and binds at -2528.5 eV.\n\n"
               "The FT 1/r^2 term grows as mu^2 there. In QM it also splits\n"
               "2S from 2P: at ETA = 8.7673e5 it would add ~357 meV to the\n"
               "muonic Lamb shift, which is measured as 202.3706(23) meV and\n"
               "fully explained without it. That bounds ETA > ~1.4e11, close\n"
               "to the electronic/muonic proton-radius retrofit (ETA > 6.6e11).");
        addRow(form, "eta preset", cbEtaPreset,
               "eta is a charge per length (C/m), entered as ETA = eta*a0/e.\n\n"
               "lambda_e = c^2/sqrt(K G)   ETA 3.8327e25   eta 1.1604e17 C/m\n"
               "  the manuscript's value, lambda_g/lambda_eg (Sec. 2.5.5).\n"
               "  Equal to Planck charge / Planck length exactly (hbar\n"
               "  cancels). Crossover 2.6e-26 a0: FT invisible here.\n"
               "retrofit                   ETA 1.044e12   eta 3.2e3 C/m\n"
               "  joint fit of r_p and r_c = e/eta to the muonic proton\n"
               "  radius and five electronic H measurements: the FT 1/r^2\n"
               "  term mimics a larger proton in electronic H (k ~ 6 for\n"
               "  1S-3S, 4 for the Lamb shift) and not in muonic H.\n"
               "  r_c = (5.1 +/- 1.5)e-23 m, 3.4 sigma - the proton-radius\n"
               "  puzzle remnant; treat as an upper limit, not a detection.\n"
               "95% floor                  ETA 6.606e11   eta 2.0e3 C/m\n"
               "  r_c < 8.0e-23 m from the same fit. Supersedes the old\n"
               "  '10 Hz' floor, which ignored that R_inf and r_p are\n"
               "  themselves extracted from hydrogen.\n"
               "nuclear, r_c = 1 fm         ETA 5.2918e4   - excluded\n"
               "visible, r_c = 1 a0         ETA 1          - excluded\n\n"
               "c*eps0 is NOT offered: it is 1/Z0 = 2.654e-3 siemens, not\n"
               "C/m, so eta*x would not have the units of a charge.\n\n"
               "Cluster layouts are floored at 1 fm radius so a large eta\n"
               "does not pack nucleons 1e-36 m apart.");
        addRow(form, "ETA = eta*a0/e", sbEta,
               "Electric time-dilation scale, ETA = eta * a0 / e.\n"
               "f(r) = ETA*r / (ETA*r + S),  S = |q_source|.\n\n"
               "eta must be in C/m. Crossover sits at r_c = S / ETA.\n"
               "Use the preset list above, or type any value (scientific\n"
               "notation accepted, e.g. 3.8327e25). Wheel/arrows step x1.26.");
        addRow(form, "dt", sbDt,
               "Velocity-Verlet time step, in units of a0/v_Bohr = 2.4189e-17 s.\n\n"
               "2e-4 is fine for hydrogen (period 6.283).\n"
               "Use dt <= 1e-5 for Be3+ - its period is 16x shorter.\n"
               "Watch the E readout: on velocity-only kinematics (SR, FT 1/2,\n"
               "FT Eq.41, Newtonian) a drift means dt is too large. With the\n"
               "f^2 momentum factors energy is NOT conserved by the law\n"
               "itself, so drift there is expected and says nothing about dt.");
        addRow(form, "box half-width (a0)", sbBox,
               "Half-width of the square the random layout draws from.\n"
               "Ignored by the cluster layouts, which use r_c instead.");
        addRow(form, "initial |v| (v_Bohr)", sbV,
               "Initial speed scale, in units of v_Bohr = 2.187691e6 m/s.\n\n"
               "Applied with equipartition, v ~ 1/sqrt(m), so protons and\n"
               "neutrons start about 43x slower than the electrons.");

        addRow(form, "initial layout", cbLayout,
               "How nucleon positions are drawn.\n\n"
               "random         - uniform in the box\n"
               "proton cluster - protons packed inside r < r_c\n"
               "mixed cluster  - protons AND neutrons inside r < r_c\n\n"
               "The attractive region exists only on the Potential route.\n"
               "Even there an 8-proton cluster expands ~235x - slower than\n"
               "Coulomb's 1692x, but it does not bind.");
        addRow(form, "cluster fill (x r_c)", sbFill,
               "Cluster disk radius = fill * r_c / 2, so every pair starts\n"
               "within fill * r_c of every other.\n\n"
               "fill < 1 puts all pairs inside the attractive region.\n"
               "Below ~0.2 the cluster holds together long enough to film.");

        addRow(form, "acceleration", cbRoute,
               "Which acceleration the same f is turned into.\n"
               "These are NOT equivalent once f != 1:\n\n"
               "Potential    a = -(1/m) d(U f^2)/dr\n"
               "             = -C eta^2 (S - eta r)/(eta r + S)^3\n"
               "             changes sign at r_c: like charges ATTRACT inside\n\n"
               "Centripetal  a = v^2 / r  (the Eq.103 derivation)\n"
               "             = C eta^2/(eta r + S)^2, sign-definite, no crossover\n\n"
               "Coulomb      f = 1, plain inverse square, for reference.\n\n"
               "Discrete crossover\n"
               "             the crossover is ROUNDED to the nearest shell:\n"
               "               n     = round(sqrt(r*S))\n"
               "               eta_L = S^2/n^2   so S/eta_L = n^2/Z = r_n\n"
               "             giving a STABLE zero at every r_n = n^2 a0/Z and\n"
               "             an unstable rim between shells at (n+1/2)^2.\n"
               "             Repulsive inside the ground shell, so nothing\n"
               "             collapses onto the nucleus.\n\n"
               "             It does NOT recover Coulomb: |a|/|a_Coulomb| runs\n"
               "             0.03 - 0.15, not 1, because the crossover tracks r.\n"
               "             It does NOT derive quantisation either - the shells\n"
               "             come from L = n hbar and are hard-coded here.\n\n"
               "Both FT routes are finite at r = 0, so no softening is used.");

        addRow(form, "floor quantisation", cbQuantMode,
               "Applied on top of the selected acceleration route.\n\n"
               "floor the RADIUS   a is evaluated at r_q = floor(r/D)*D,\n"
               "                   so space comes in lumps of D. Bin 0 is\n"
               "                   clamped to D/2 to avoid 1/0.\n"
               "floor the FORCE    a_q = floor(a/D)*D, so the force itself\n"
               "                   takes only multiples of D.\n\n"
               "Both make the force a staircase and both BREAK momentum\n"
               "conservation slightly, since a_ij and a_ji get rounded\n"
               "independently. Watch the E readout.");
        addRow(form, "step D", sbStep,
               "Quantisation step.\n\n"
               "Radius mode: D is in a0. Visible staircase at D ~ 0.1 - 1.\n"
               "Force mode:  D is in force units (1 = Ke^2/a0^2 = 8.2387e-8 N).\n"
               "             Visible at D ~ 0.05 - 0.2.\n\n"
               "Note on the literal floor(f(x)/a0)*a0: f is DIMENSIONLESS, so\n"
               "f/a0 is about 9.4e9 and the staircase has ~1.9e10 steps across\n"
               "0 < f < 1 - far too fine to see. Use D of order the quantity\n"
               "being stepped, not a0.");
        addRow(form, chMagnetic,
               "Biot-Savart field of every moving point charge,\n"
               "  B_z = alpha^2 q_j (v_j x r)_z / r^3\n"
               "acting on the other charges through the Lorentz force.\n\n"
               "alpha^2 = 1/c^2 = 5.32513545e-5 in these atomic units, so the\n"
               "effect is O(v^2/c^2) ~ 5e-5 of Coulomb for a ground-state orbit.\n"
               "Two electrons moving side by side in parallel feel a magnetic\n"
               "ATTRACTION reducing their Coulomb repulsion by exactly v^2/c^2.\n\n"
               "Magnetic forces between point charges break Newton's third\n"
               "law (the field carries momentum), so total momentum drifts.");
        addRow(form, chSpin,
               "Fixed intrinsic moment along +/- z for every particle:\n"
               "  electron  0.500579826 a.u. = (g_e/2) mu_B, opposite its spin\n"
               "  proton    +2.79284734 mu_N    neutron  -1.91304273 mu_N\n"
               "  He-4 0,  Li-7 +3.256427 mu_N,  Be-9 -1.17749 mu_N\n\n"
               "Gives dipole-dipole forces, F_r = 3 alpha^2 mu_i mu_j / r^4,\n"
               "spin-orbit forces (moment in the field of a moving charge,\n"
               "with the Thomas factor (g_e-1)/g_e for the electron), and the\n"
               "dipole field acting on moving charges.\n\n"
               "Freezing the spin direction is EXACT here: a z-moment in a\n"
               "z-field feels no torque, mu x B = 0.");
        addRow(form, "spin assignment", cbSpinMode,
               "random       each spin +1 or -1 at random\n"
               "all up       every spin +1 (electron moments all point -z)\n"
               "alternating  +1, -1, +1, ... in creation order");
        addRow(form, "magnetic gain (display)", sbMagGain,
               "Multiplies both magnetic terms so they can be SEEN.\n\n"
               "At gain 1 they are ~1e-5 of Coulomb and invisible on screen.\n"
               "gain != 1 is a visualisation choice, NOT physics.");
        addRow(form, chQCD,
               "Cornell potential between quarks of the SAME hadron, with\n"
               "the baryon 1/2 rule:\n"
               "  V_qq = -(2/3) a_s hbar c / r + (sigma/2) r\n"
               "  a_s = 0.35,  sigma = 0.18 GeV^2 = 912.19 MeV/fm\n\n"
               "Beyond ~0.3 fm the pull is a CONSTANT 456.10 MeV/fm\n"
               "(7.307e4 N) per pair - the flux tube. It never weakens with\n"
               "distance, so a proton/neutron triplet stays together: this\n"
               "is what was missing from the charge-only toy.\n\n"
               "Each hadron shares ONE time-dilation factor, set by charges\n"
               "outside it, and the velocity-dependent magnetic terms are\n"
               "skipped between its own quarks (they move at ~0.94c, where\n"
               "point-charge Biot-Savart is invalid and non-reciprocal).\n"
               "An isolated hadron then conserves momentum to round-off.\n\n"
               "Only within a hadron. Between hadrons a constant pull at any\n"
               "range is unphysical - real flux tubes break by pair creation -\n"
               "so it is off there. Loose quarks feel charge only.\n\n"
               "WARNING: ~8.9e11 simulator force units. Use dt <= 1e-12 with\n"
               "constituent masses, dt <= 1e-13 with current masses.\n"
               "A non-Newtonian momentum factor is required; FT is the default.");
        addRow(form, chConstituent,
               "Replace current quark masses with CONSTITUENT masses, the\n"
               "standard way to put QCD field energy into inertia when the\n"
               "dynamics is classical:\n\n"
               "  quark   current (MeV)   constituent (MeV)   (m_e)\n"
               "  u          2.16            336              657.54\n"
               "  d          4.67            340              665.36\n"
               "  s         93.4             486              951.08\n"
               "  c       1270              1550             3033.27\n"
               "  b       4180              4730             9256.38\n"
               "  t      172690             (unchanged, no hadrons form)\n\n"
               "uud = 1012 MeV vs proton 938.27; udd = 1016 vs 939.57. The\n"
               "~74 MeV excess is the colour-magnetic hyperfine term, not\n"
               "modelled. The same masses give mu_p = (4mu_u - mu_d)/3 =\n"
               "2.789 mu_N against 2.793 measured.\n\n"
               "Applies to existing quarks immediately.");
        addRow(form, "momentum factor", cbKin,
               "How a force changes velocity. The force changes p; v is\n"
               "recovered from p.\n\n"
               "FT f^2 x (m0 + K/c^2)\n"
               "  p = m v / f^2,  f = eta/(eta + sum |q_j|/r_j) at the particle\n"
               "  m = m0 + K/c^2 with K the work done (dK = v dp). Solving\n"
               "  that exactly gives m = m0/sqrt(1 - v^2/c^2): the rule IS\n"
               "  the relativistic mass, so p = gamma m0 v / f^2 and |v| < c.\n\n"
               "FT f^2 x first order\n"
               "  K truncated to m0 v^2/2:  p = m0 v (1 + v^2/2c^2) / f^2.\n"
               "  No speed limit. Muon g-2 (p = 29.28 m c): v = 3.71c,\n"
               "  period 40.16 ns vs 149.1 ns measured.\n\n"
               "FT f^2 positional only   p = m0 v / f^2\n\n"
               "With any f^2 choice, use the COULOMB route: a circular orbit\n"
               "then has v = sqrt(Kq^2/mr) f exactly (Eq. 103, the v2 sheet).\n"
               "FT force routes put f^2 in the force as well - counted twice.\n"
               "Inertia m0/f^2 varies with position, so energy is NOT\n"
               "conserved off circular orbits; the E readout will drift.\n\n"
               "FT 1/2 (DEFAULT), FT Eq.41, SR, Newtonian: velocity-only\n"
               "factors (FT 1/2 limit sqrt2 c, Eq.41 limit c).\n"
               "QCD confinement needs a speed limit and switches the three\n"
               "unlimited choices back to the default.");
        addRow(form, "QCD gain (display)", sbQcdGain,
               "Scales the Cornell term. 1 = alpha_s 0.35, sigma 0.18 GeV^2.\n"
               "Reduce it to run at a larger dt; gain != 1 is not physics.");
        addRow(form, chStrong,
               "Malfliet-Tjon V, a standard NN parametrisation, between\n"
               "NUCLEONS ONLY (protons and neutrons):\n\n"
               "  V(r) = [1438.72 e^-3.11r - 626.885 e^-1.55r]/r  MeV, r in fm\n"
               "  well depth -77.5 MeV at 0.819 fm, V = 0 at 0.5325 fm\n\n"
               "CHARGE-INDEPENDENT: pp, pn and nn get the identical potential.\n"
               "That is the physical point - it is why the deuteron (pn) binds\n"
               "while the FT crossover, coupling to charge, gives pn exactly 0.\n\n"
               "Electrons and composite nuclei are excluded.\n\n"
               "WARNING: ~1.5e11 simulator force units at 1 fm, so a proton\n"
               "sees ~8e7 in acceleration. Set dt <= 1e-9 first or the\n"
               "integrator will explode.");
        addRow(form, "strong gain (display)", sbSGain,
               "Scales the strong term. 1 = true Malfliet-Tjon strength.\n\n"
               "Reduce it to keep a cluster on screen at a workable dt;\n"
               "gain != 1 is a visualisation choice, not physics.");
        addRow(form, chGravity,
               "Add the gravitoelectric term through the same kernel, with\n"
               "C = -G m_j, S = m_j, eta = eta_g = c^2/G.\n\n"
               "G = 2.400610e-43 in these units, so between two electrons it\n"
               "is 43 orders below the electric force and changes nothing.");
        addRow(form, chGravMag,
               "Mass-current analogue of the moving-charge field:\n"
               "  B_g = (G/c^2) m_j (v_j x r)/r^3 * f_g^2\n"
               "  a_i = -k (v_i x B_g)        independent of m_i\n\n"
               "Sign flipped relative to EM: parallel mass currents REPEL.\n"
               "Relative to Newton the force is -k v_i v_j / c^2.\n\n"
               "G/c^2 = 1.27836e-47 here, so at gravity gain 1 it is ~5e-47\n"
               "of the electric force. It uses the gravity gain.\n\n"
               "This is the gravitomagnetic piece only, not the complete 1PN\n"
               "(EIH) dynamics, which also corrects the gravitoelectric term.");
        addRow(form, "GEM coupling", cbGmK,
               "The coefficient depends on the convention:\n\n"
               "k = 4    linearised GR, h_0i = -4 G m v_i/(c^3 r)\n"
               "         (Mashhoon's GEM, Lense-Thirring, Gravity Probe B)\n"
               "k = 1.5  the manuscript's mu_g = 6 pi G/c^2 (Eq. 19)\n"
               "k = 1    naive EM analogy, mu_g = 4 pi G/c^2\n\n"
               "Gravity Probe B measured frame dragging at 37.2 +/- 7.2 mas/yr\n"
               "against the GR (k = 4) prediction of 39.2 mas/yr.");
        addRow(form, "ETA_G = eta_g*a0/m_e", sbEtaG,
               "Gravitational time-dilation scale.\n"
               "7.822540e46 corresponds to eta_g = c^2/G = 1.3466e27 kg/m.\n\n"
               "Its crossover for one electron mass is 1.28e-47 a0 = 6.76e-58 m,\n"
               "about 10^-19 of the Planck length.");
        addRow(form, "gravity gain (display)", sbGain,
               "Multiplies the gravitational term so it can be SEEN.\n\n"
               "gain != 1 is a visualisation choice, NOT physics.\n"
               "Gravity would need gain = 4.1656e42 to match the electric\n"
               "force between two electrons at 1 a0.");

        addRow(form, "zoom (log, px per a0)", slZoom,
               "Logarithmic: scale = 10^(v/100 - 1), spanning 0.1 to 1e12.\n\n"
               "Mouse wheel zooms about the cursor, left-drag pans,\n"
               "double-click resets. Grid rings relabel as you zoom.\n"
               "A 1 fm cluster needs about 1e6 px/a0 to be visible.");
        addRow(form, chQuant,
               "Hold every electron on a stationary state instead of\n"
               "integrating its radial equation.\n\n"
               "  r_n = n^2 a0 / Z      (f cancels EXACTLY)\n"
               "  v_n = Z/n v_Bohr,  omega_n = Z^2/n^3\n"
               "  f_n = ETA n^2/(ETA n^2 + Z^2)\n"
               "  E_n = -(Z^2/2n^2) f_n^2  Hartree\n\n"
               "So the dilation contributes NOTHING to the radii and shows up\n"
               "only in the energies, as a 1/n^4 shift. Electrons held this way\n"
               "neither radiate nor spiral - that is the Bohr postulate made\n"
               "explicit, IMPOSED here rather than derived: the FT force has a\n"
               "degree-1 numerator, so it has exactly one root and cannot\n"
               "generate a series of radii by itself.");
        addRow(form, "levels n = 1..", sbNmax,
               "Highest principal quantum number handed out. Electrons get\n"
               "levels round-robin from 1 to this value.\n\n"
               "Shell radii scale as n^2/Z, periods as 2 pi n^3/Z^2.");
        addRow(form, chLabels,
               "Draw each particle's symbol on its disc:\n"
               "  e  p  n   mu  tau   nu_e nu_mu nu_tau   u d s c b t\n"
               "  composite nuclei by element: H He Li Be\n\n"
               "Standard symbols are used rather than literal first letters,\n"
               "which would collide (neutron/neutrino 'n', tau/top 't').\n"
               "Discs grow to fit the letter while this is on.");
        addRow(form, chTrails,
               "Draw the last 600 positions of each particle.\n"
               "Turn off for a large N - trails dominate the frame time.");

        auto *ionHdr = new QLabel("<b>Hydrogen-like ions</b>");
        ionHdr->setToolTip(tipHtml(QString(
               "Presets place a nucleus (Z protons, N neutrons) and ONE\n"
               "electron on a circular orbit at r = a0/Z with v = Z v_Bohr.\n"
               "Verified circular to 9 digits over 3 periods.")));
        form->addRow(ionHdr);
        addRow(form, "preset", cbIon,
               "ion    Z  N   r (a0)   v (v_Bohr)   period\n"
               "H      1  0   1.0000   1            6.2832\n"
               "He+    2  2   0.5000   2            1.5708\n"
               "Li2+   3  4   0.3333   3            0.69813\n"
               "Be3+   4  5   0.2500   4            0.39270\n\n"
               "'all four' spaces them 4 a0 apart.");
        addRow(form, "hadron radius (fm)", sbHadR,
               "Radius of the quark triangle for the uud / udd entries and\n"
               "the Quark tab (quark-quark spacing = radius x sqrt 3).\n"
               "Default 0.8409 fm, the proton rms charge radius.\n\n"
               "For the FT attraction a like pair must sit inside\n"
               "r_c = sqrt(|q_i q_j|)/ETA: for u-u that is (2/3)/ETA a0, e.g.\n"
               "3.4e-23 m at the retrofit ETA, 0.67 fm at ETA = 5.29e4.\n"
               "Cornell (QCD on) and Coulomb both grow as 1/r^2 inside, so a\n"
               "much smaller radius needs a much smaller dt.");
        addRow(form, chResolve,
               "Place the nucleons individually instead of one composite\n"
               "particle of charge +Ze.\n\n"
               "They WILL fly apart - this simulator has no strong force,\n"
               "and that is the correct outcome, not a bug. The exception is\n"
               "the Potential route with ETA ~ 5.3e4, where r_c reaches 1 fm.");
        addRow(form, btIon,
               "Clear the scene and load the selected ion preset.\n"
               "Uses the current ETA, route and gravity settings.");

        auto *btns = new QHBoxLayout;
        btns->addWidget(btReset);
        btns->addWidget(btRestart);
        btns->addWidget(btClear);
        btns->addWidget(btPause);

        auto *note = new QLabel(
            "ETA = 3.8327e25 is the manuscript's lambda_e = c^2/sqrt(KG);\n"
            "its crossover is 2.6e-26 a0 and the run is plain Coulomb.\n"
            "Use ETA ~ 0.1 - 10 to make the FT structure visible.\n\n"
            "Gravity: G = 2.4006e-43 in these units, so between two electrons\n"
            "it is 43 orders below the electric force and changes nothing at\n"
            "gain = 1. ETA_G = 7.8225e46 puts its crossover at 1.28e-47 a0\n"
            "(6.76e-58 m), far below the Planck length. Raising the gain is a\n"
            "visualisation aid, not physics.\n\n"
            "Neutrons carry mass only. With gravity off they feel no force at\n"
            "all and drift in straight lines - that is correct, not a bug.\n"
            "Initial speeds follow equipartition, v ~ 1/sqrt(m), so the heavy\n"
            "species start ~43x slower than the electrons.\n\n"
            "Ions: the nucleus is one composite particle (q=+Ze) because there\n"
            "is no strong force here. 'Resolve nucleons' places them singly and\n"
            "they fly apart - correct, not a bug - unless ETA is low enough that\n"
            "the crossover reaches nuclear range (ETA ~ 5.3e4 puts it at 1 fm),\n"
            "where like charges attract at short range and the cluster holds.\n"
            "Use dt <= 1e-5 for Be3+: its orbital period is 16x shorter than H's.\n\n"
            "Cluster layouts pack nucleons in a disk of radius fill*r_c/2, so\n"
            "EVERY pair starts within fill*r_c and is inside the attractive\n"
            "region. This only attracts on the Potential route - on Centripetal\n"
            "and Coulomb like charges always repel and the cluster explodes.\n"
            "Neutrons in the mixed cluster feel nothing unless gravity is on.\n"
            "For a nuclear-scale r_c (1 fm = 1.8897e-5 a0) set ETA = 5.292e4,\n"
            "then zoom to ~1e6 px/a0 to see the cluster.\n\n"
            "View: mouse wheel zooms about the cursor, left-drag pans,\n"
            "double-click resets. The grid rings relabel themselves as you\n"
            "zoom, so the ring spacing is always a readable power of ten.\n"
            "The ruler at the foot of the legend snaps to 1-2-5 steps and\n"
            "reads in both a0 and SI - fm at nuclear scale, A at atomic.");
        note->setWordWrap(true);
        note->setStyleSheet("color:#9AA6B8; font-size:11px;");

        auto *box = new QGroupBox("Controls");
        auto *bl  = new QVBoxLayout(box);
        bl->addLayout(form);
        bl->addLayout(btns);
        bl->addWidget(note);
        bl->addStretch();

        // tooltips here are pre-formatted tables; keep the line breaks
        QToolTip::setFont(QFont("monospace", 9));

        // controls scroll vertically when the window is shorter than the form;
        // the splitter lets the panel be dragged wider or narrower
        auto *scroll = new QScrollArea;
        scroll->setWidget(box);
        scroll->setWidgetResizable(true);
        scroll->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
        scroll->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
        scroll->setFrameShape(QFrame::NoFrame);
        scroll->setMinimumWidth(box->sizeHint().width()
                                + scroll->verticalScrollBar()->sizeHint().width() + 4);

        auto *split = new QSplitter(Qt::Horizontal);
        split->addWidget(view);
        split->addWidget(scroll);
        split->setStretchFactor(0, 1);
        split->setStretchFactor(1, 0);
        split->setChildrenCollapsible(false);

        // ---- top-level scale tabs ----
        tabs = new QTabBar;
        tabs->addTab("Galactic");  tabs->addTab("Solar");
        tabs->addTab("Atomic");    tabs->addTab("Quark");
        tabs->setToolTip(tipHtml("Each tab keeps its own scene: switching away freezes it,\n"
                                 "switching back resumes it where it was."));
        tabs->setTabToolTip(0, tipHtml(
            "1e10 M_sun bulge + 200 tracer stars of 1e3 M_sun on a disk to\n"
            "15 kpc, each started on the circular orbit the current force\n"
            "law gives. The force route is left as the profile set it.\n"
            "CAUTION: with the symmetric pair scale S = sqrt(M m) the\n"
            "bulge/tracer crossover sits at 0.0016 kpc, not M/eta_g = 5 kpc,\n"
            "so the tracers orbit Keplerian (107 km/s at 3.7 kpc, 54 at 15)\n"
            "whatever eta_g. The paper's Eq. 103 uses S = M alone.\n"
            "ETA_G preset 'galactic fit': r_c = M/eta_g = 5 kpc, the per-galaxy\n"
            "fitted value of Sec. 4.4.1 (c^2/G puts r_c at 7e13 m, invisible).\n"
            "dt = 1e30 (~300 steps per orbit at 8 kpc). Ruler reads kpc."));
        tabs->setTabToolTip(1, tipHtml(
            "Sun + eight planets at mean distance, circular, random phases.\n"
            "eta_g = c^2/G: the FT correction is GM/(c^2 r) ~ 1e-8 at Earth,\n"
            "so this is Newtonian to the eye. dt = 1e21 (~300 steps per\n"
            "Mercury orbit). Ruler reads AU."));
        tabs->setTabToolTip(2, tipHtml(
            "20 electrons, 20 protons, 20 neutrons in the current layout.\n"
            "ETA preset 'visible' (ETA = 1, r_c = 1 a0) so the FT structure\n"
            "shows; dt = 0.02. Choose 'retrofit' (1.04e12, from the\n"
            "electronic/muonic proton-radius fit) for the bounded value,\n"
            "where f = 1 here to 1e-12."));
        tabs->setTabToolTip(3, tipHtml(
            "One resolved uud proton: constituent masses, Cornell confinement,\n"
            "SR kinematics (speed limit c), dt = 1e-12, zoomed to ~1 fm.\n"
            "ETA preset 'retrofit': at r_c = 5e-23 m the FT correction to\n"
            "quark binding at 1 fm is ~1e-7 relative, so QCD does the work."));

        auto *root = new QVBoxLayout(this);
        root->setContentsMargins(4, 4, 4, 4);
        root->addWidget(tabs);
        root->addWidget(split, 1);
        connect(tabs, &QTabBar::currentChanged, this, &Window::switchTab);

        // ---- top-level Presets menu ----
        auto *mb = new QMenuBar(this);
        root->setMenuBar(mb);
        QMenu *mPre = mb->addMenu("&Presets");
        QAction *aSR = mPre->addAction("SR at all levels");
        QAction *aFT = mPre->addAction("FT at all levels");
        aSR->setToolTip("Mainstream physics everywhere");
        aFT->setToolTip("Finite Theory everywhere");
        mPre->setToolTipsVisible(true);
        connect(aSR, &QAction::triggered, this, [this] { applyProfile(false); });
        connect(aFT, &QAction::triggered, this, [this] { applyProfile(true);  });

        connect(btReset, &QPushButton::clicked, this, &Window::doReset);
        connect(sbHadR, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) { sim.hadronRadiusFm = v; });
        connect(btRestart, &QPushButton::clicked, this, &Window::restartScene);
        connect(btClear, &QPushButton::clicked, this, [this] {
            syncSettings(); sim.clear(); view->update();
        });
        connect(btPause, &QPushButton::clicked, this, [this] {
            if (timer.isActive()) { timer.stop();  btPause->setText("Resume"); }
            else                  { timer.start(); btPause->setText("Pause");  }
        });
        connect(sbEta, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) {
            sim.eta = v; sim.computeAcc();
            // typed by hand: mark the preset list as custom unless it matches
            const double pv = cbEtaPreset->currentData().toDouble();
            if (pv <= 0.0 || std::fabs(v - pv) > 1e-6 * pv) {
                cbEtaPreset->blockSignals(true);
                cbEtaPreset->setCurrentIndex(cbEtaPreset->count() - 1);
                cbEtaPreset->blockSignals(false);
            }
            view->update();
        });
        connect(cbEtaPreset, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int) {
            const double v = cbEtaPreset->currentData().toDouble();
            if (v > 0.0) sbEta->setValue(v);          // 'custom' leaves it alone
        });
        connect(cbRoute, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int i) {
                    sim.route = (i == 0) ? Route::Potential
                              : (i == 1) ? Route::Centripetal
                              : (i == 2) ? Route::Coulomb
                                         : Route::DiscreteCrossover;
                    sim.computeAcc();
                });
        connect(btAddSpecies, &QPushButton::clicked, this, [this] {
            static const int LEPTONS[] = {SP_E, SP_MU, SP_TAU,
                                          SP_NUE, SP_NUMU, SP_NUTAU};
            static const int QUARKS[]  = {SP_U, SP_D, SP_S, SP_C, SP_B, SP_T};
            const int idx   = cbSpecies->currentIndex();
            const int group = cbSpecies->currentData().isValid()
                              ? cbSpecies->currentData().toInt() : -1;
            const int n = sbSpecN->value();
            const double box = sbBox->value(), vs = sbV->value();

            auto addSet = [&](const int *set, int len) {
                for (int k = 0; k < len; ++k)
                    sim.addSpecies(set[k], n, box, vs, seed++);
            };
            if (group == 1000 || group == 1002) addSet(LEPTONS, 6);
            if (group == 1001 || group == 1002) addSet(QUARKS, 6);
            if (group == 1003) sim.addHadron(true,  n, box, vs, seed++);
            if (group == 1004) sim.addHadron(false, n, box, vs, seed++);
            if (group < 1000 && idx >= 0 && idx < SP_COUNT)
                sim.addSpecies(idx, n, box, vs, seed++);
            if (sim.spinOn)    { sim.assignSpins(seed++); sim.computeAcc(); }
            if (sim.quantised)   sim.assignLevels();
            view->update();
        });
        connect(chGravMag, &QCheckBox::toggled, this, [this](bool b) {
            sim.gravMag = b; sim.computeAcc(); view->update();
        });
        connect(cbGmK, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int i) {
            sim.gmK = (i == 0) ? 4.0 : (i == 1) ? 1.5 : 1.0;
            sim.computeAcc(); view->update();
        });
        connect(chMagnetic, &QCheckBox::toggled, this, [this](bool b) {
            sim.magnetic = b; sim.computeAcc(); view->update();
        });
        connect(chSpin, &QCheckBox::toggled, this, [this](bool b) {
            sim.spinOn = b;
            if (b) sim.assignSpins(seed++); else sim.clearSpins();
            sim.computeAcc(); view->update();
        });
        connect(cbSpinMode, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int i) {
            sim.spinMode = i;
            if (sim.spinOn) { sim.assignSpins(seed++); sim.computeAcc(); }
            view->update();
        });
        connect(sbMagGain, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) { sim.magGain = v; sim.computeAcc();
                                         view->update(); });
        connect(chQCD, &QCheckBox::toggled, this, [this](bool b) {
            sim.qcd = b;
            // QCD needs a speed limit: first-order, positional-only and
            // Newtonian have none
            {
                const int k = cbKin->currentIndex();
                if (b && (k == 1 || k == 2 || k == 6)) cbKin->setCurrentIndex(3); // FT 1/2
            }
            sim.computeAcc(); view->update();
        });
        connect(cbKin, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int i) {
            static const Sim::Kin K[] = {Sim::Kin::F2SC, Sim::Kin::F2Lin,
                Sim::Kin::F2Pos, Sim::Kin::FTHalf, Sim::Kin::FT, Sim::Kin::SR,
                Sim::Kin::Newtonian};
            sim.kin = K[qBound(0, i, 6)];
            sim.computeAcc();
            view->update();
        });
        connect(chConstituent, &QCheckBox::toggled, this, [this](bool b) {
            sim.constituent = b; sim.applyQuarkMasses(); view->update();
        });
        connect(sbQcdGain, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) { sim.qcdGain = v; sim.computeAcc();
                                         view->update(); });
        connect(chStrong, &QCheckBox::toggled, this, [this](bool b) {
            sim.strong = b; sim.computeAcc(); view->update();
        });
        connect(sbSGain, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) { sim.strongGain = v; sim.computeAcc();
                                         view->update(); });
        connect(cbQuantMode, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int i) {
            sim.quant = (i == 0) ? Quant::None
                      : (i == 1) ? Quant::Radius : Quant::Force;
            sim.computeAcc(); view->update();
        });
        connect(sbStep, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) { sim.quantStep = v; sim.computeAcc();
                                         view->update(); });
        connect(chQuant, &QCheckBox::toggled, this, [this](bool b) {
            sim.quantised = b;
            if (b) sim.assignLevels();
            else for (auto &a : sim.p) { a.level = 0; a.anchor = -1; }
            view->update();
        });
        connect(sbNmax, QOverload<int>::of(&QSpinBox::valueChanged),
                this, [this](int v) {
            sim.nMax = v;
            if (sim.quantised) sim.assignLevels();
            view->update();
        });
        connect(cbLayout, QOverload<int>::of(&QComboBox::currentIndexChanged),
                this, [this](int) { restartScene(); });
        connect(sbFill, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) { sim.clusterFill = v; restartScene(); });
        connect(sbEtaG, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) { sim.etaG = v; sim.computeAcc(); });
        connect(sbGain, QOverload<double>::of(&QDoubleSpinBox::valueChanged),
                this, [this](double v) { sim.gravGain = v; sim.computeAcc(); });
        connect(chGravity, &QCheckBox::toggled, this,
                [this](bool b) { sim.gravity = b; sim.computeAcc(); view->update(); });
        connect(slZoom, &QSlider::valueChanged, this, [this](int v) {
            view->scale = std::pow(10.0, v / 100.0 - 1.0);
            view->update();
        });
        // wheel zoom writes back to the slider without re-triggering it
        connect(view, &View::viewChanged, this, [this](double sc) {
            slZoom->blockSignals(true);
            slZoom->setValue(int(100.0 * (std::log10(sc) + 1.0)));
            slZoom->blockSignals(false);
        });
        connect(chLabels, &QCheckBox::toggled, this,
                [this](bool b) { view->showLabels = b; view->update(); });
        connect(chTrails, &QCheckBox::toggled, this,
                [this](bool b) { view->showTrails = b; view->update(); });

        connect(btIon, &QPushButton::clicked, this, [this] {
            syncSettings();
            sim.loadIons(cbIon->currentIndex(), chResolve->isChecked(), seed++);
            if (sim.quantised) sim.assignLevels();
            if (sim.spinOn) { sim.assignSpins(seed++); sim.computeAcc(); }
            view->update();
        });

        connect(&timer, &QTimer::timeout, this, &Window::tick);
        timer.setInterval(16);

        tabs->setCurrentIndex(2);          // emits currentChanged -> switchTab(2)
        if (sim.p.empty()) switchTab(2);   // in case index 2 was already current
        timer.start();
        setWindowTitle("Finite Theory - N charge simulator");
    }

private slots:
    // copy every control into the simulation, so a reset never runs with a
    // stale value (earlier, layout and cluster fill were never applied here)
    void syncSettings()
    {
        sim.eta      = sbEta->value();
        sim.etaG     = sbEtaG->value();
        sim.gravGain = sbGain->value();
        sim.gravity  = chGravity->isChecked();
        sim.gravMag  = chGravMag->isChecked();
        sim.gmK      = (cbGmK->currentIndex() == 0) ? 4.0
                     : (cbGmK->currentIndex() == 1) ? 1.5 : 1.0;
        sim.route = (cbRoute->currentIndex() == 0) ? Route::Potential
                  : (cbRoute->currentIndex() == 1) ? Route::Centripetal
                  : (cbRoute->currentIndex() == 2) ? Route::Coulomb
                                                   : Route::DiscreteCrossover;
        sim.layout = (cbLayout->currentIndex() == 0) ? Layout::Random
                   : (cbLayout->currentIndex() == 1) ? Layout::ProtonCluster
                                                     : Layout::MixedCluster;
        sim.clusterFill = sbFill->value();
        sim.quant = (cbQuantMode->currentIndex() == 0) ? Quant::None
                  : (cbQuantMode->currentIndex() == 1) ? Quant::Radius
                                                       : Quant::Force;
        sim.quantStep   = sbStep->value();
        sim.magnetic    = chMagnetic->isChecked();
        sim.spinOn      = chSpin->isChecked();
        sim.spinMode    = cbSpinMode->currentIndex();
        sim.magGain     = sbMagGain->value();
        sim.strong      = chStrong->isChecked();
        sim.strongGain  = sbSGain->value();
        sim.hadronRadiusFm = sbHadR->value();
        sim.qcd         = chQCD->isChecked();
        {
            static const Sim::Kin K[] = {Sim::Kin::F2SC, Sim::Kin::F2Lin,
                Sim::Kin::F2Pos, Sim::Kin::FTHalf, Sim::Kin::FT, Sim::Kin::SR,
                Sim::Kin::Newtonian};
            sim.kin = K[qBound(0, cbKin->currentIndex(), 6)];
        }
        sim.constituent = chConstituent->isChecked();
        sim.qcdGain     = sbQcdGain->value();
        sim.nMax        = sbNmax->value();
        sim.quantised   = chQuant->isChecked();
    }

    // One-click profiles.  Only the SR-vs-FT choices are touched; species
    // counts, spin, strong force, QCD, dt and the view are left as they are.
    //
    //   level                   SR at all levels          FT at all levels
    //   electric force route    Coulomb (f = 1)           Potential -d(U f^2)/dr
    //   momentum factor         SR, gamma                 FT 1/2
    //   moving-charge field     on (Biot-Savart)          on
    //   gravitoelectric         on, Newtonian             on, FT kernel, eta_g = c^2/G
    //   gravitomagnetic         on, k = 4 (GR)            on, k = 1.5 (Eq. 19)
    //   eta                     lambda_e (f = 1 anyway)   visible preset, ETA = 1
    //   layout                  random                    mixed cluster
    //   floor quantisation      none                      none
    void applyProfile(bool ft)
    {
        cbRoute->setCurrentIndex(ft ? 0 : 2);        // Potential / Coulomb
        cbKin->setCurrentIndex(ft ? 3 : 5);          // FT 1/2   / SR
        chMagnetic->setChecked(true);
        chGravity->setChecked(true);
        chGravMag->setChecked(true);
        cbGmK->setCurrentIndex(ft ? 1 : 0);          // k = 1.5  / k = 4
        sbEtaG->setValue(7.822540e46);               // eta_g = c^2/G
        cbEtaPreset->setCurrentIndex(ft ? 4 : 0);    // visible  / lambda_e
        cbQuantMode->setCurrentIndex(0);
        cbLayout->blockSignals(true);                // one reset, below
        cbLayout->setCurrentIndex(ft ? 2 : 0);       // mixed cluster / random
        cbLayout->blockSignals(false);
        profileName = ft ? "FT at all levels" : "SR at all levels";
        for (auto &t : tabState) t.init = false;   // other tabs rebuild under the new profile
        applyScale(curTab);
    }

private:                                  // data and helpers: not slots
    static double sliderMag(int s) { return (s <= -600) ? 0.0 : std::pow(10.0, s / 10.0); }

    // read the four slider/dial pairs into the simulation and label them
    void applyFields()
    {
        static const double SI[4]  = {5.14220675e11, 2.35051757e5, 9.0e22, 4.1341e16};
        static const char  *U[4]   = {"V/m", "T", "m/s\u00B2", "s\u207B\u00B9"};
        double m[4], th[4];
        for (int k = 0; k < 4; ++k) {
            m[k]  = sliderMag(slField[k]->value());
            th[k] = dlField[k]->value() * M_PI / 180.0;
            lbField[k]->setText(m[k] == 0.0 ? QString("off")
                : QString("%1 a.u. = %2 %3   %4\u00B0")
                      .arg(m[k], 0, 'g', 3).arg(m[k] * SI[k], 0, 'g', 3)
                      .arg(QString::fromUtf8(U[k])).arg(dlField[k]->value()));
        }
        sim.extEx  = m[0] * std::cos(th[0]);  sim.extEy = m[0] * std::sin(th[0]);
        sim.extBz  = m[1] * std::cos(th[1]);
        sim.extGx  = m[2] * std::cos(th[2]);  sim.extGy = m[2] * std::sin(th[2]);
        sim.extBgz = m[3] * std::cos(th[3]);
        sim.computeAcc();
        view->update();
    }

    // put the widgets back to what a restored tab's simulation holds
    void fieldsFromSim()
    {
        auto set = [this](int k, double mag, double deg) {
            const int s = (mag <= 0.0) ? -600
                        : qBound(-600, int(std::lround(10.0 * std::log10(mag))), 40);
            slField[k]->blockSignals(true); slField[k]->setValue(s); slField[k]->blockSignals(false);
            int d = int(std::lround(deg)) % 360; if (d < 0) d += 360;
            dlField[k]->blockSignals(true); dlField[k]->setValue(d); dlField[k]->blockSignals(false);
        };
        auto deg = [](double y, double x) { return std::atan2(y, x) * 180.0 / M_PI; };
        set(0, std::hypot(sim.extEx, sim.extEy), deg(sim.extEy, sim.extEx));
        set(1, std::fabs(sim.extBz),  sim.extBz  < 0 ? 180.0 : 0.0);
        set(2, std::hypot(sim.extGx, sim.extGy), deg(sim.extGy, sim.extGx));
        set(3, std::fabs(sim.extBgz), sim.extBgz < 0 ? 180.0 : 0.0);
        applyFields();
    }

    QString profileName = "FT at all levels";

    void setZoom(double s)
    {
        view->homeScale = s; view->scale = s; view->pan = QPointF(0, 0);
        slZoom->blockSignals(true);
        slZoom->setValue(int(100.0 * (std::log10(s) + 1.0)));
        slZoom->blockSignals(false);
    }

    // B may be a base of W (QCheckBox::setChecked lives in QAbstractButton)
    template <class W, class B, class V> static void quiet(W *w, void (B::*set)(V), V v)
    { w->blockSignals(true); (w->*set)(v); w->blockSignals(false); }

private slots:

    // ---- scale tabs: scene, time step, zoom and the presets that fit ----
    // Physics choices (route, momentum factor, k) stay with the profile; a
    // tab only changes what the Presets menu does not.  The force route is
    // never touched: under FT it stays Potential at every scale.  Note for
    // galactic: the symmetric pair scale S = sqrt(M m) (Newton's-third-law
    // fix) puts the bulge/tracer crossover at sqrt(M m)/eta_g = 0.0016 kpc,
    // not at M/eta_g = 5 kpc, so tracers orbit Keplerian on either route.
    void parkTab(int k)
    {
        if (k < 0 || k > 3) return;
        TabState &t = tabState[k];
        t.init = true;          t.sim = sim;
        t.eta = sbEta->value(); t.etaG = sbEtaG->value();
        t.dt = sbDt->value();   t.box = sbBox->value();
        t.scale = view->scale;  t.home = view->homeScale;  t.pan = view->pan;
        t.etaPreset = cbEtaPreset->currentIndex();
        t.route = cbRoute->currentIndex();  t.kin = cbKin->currentIndex();
        t.layout = cbLayout->currentIndex();
        t.qcd = chQCD->isChecked();         t.constituent = chConstituent->isChecked();
    }

    // leave the running scene frozen in its tab; resume (or build) the other
    void switchTab(int k)
    {
        if (k == curTab) return;
        parkTab(curTab);
        curTab = k;
        TabState &t = tabState[k];
        if (!t.init) { applyScale(k); return; }
        sim = t.sim;
        quiet(sbEta,  &QDoubleSpinBox::setValue, t.eta);
        quiet(sbEtaG, &QDoubleSpinBox::setValue, t.etaG);
        quiet(sbDt,   &QDoubleSpinBox::setValue, t.dt);
        quiet(sbBox,  &QDoubleSpinBox::setValue, t.box);
        quiet(cbEtaPreset, &QComboBox::setCurrentIndex, t.etaPreset);
        quiet(cbRoute,     &QComboBox::setCurrentIndex, t.route);
        quiet(cbKin,       &QComboBox::setCurrentIndex, t.kin);
        quiet(cbLayout,    &QComboBox::setCurrentIndex, t.layout);
        quiet(chQCD,         &QCheckBox::setChecked, t.qcd);
        quiet(chConstituent, &QCheckBox::setChecked, t.constituent);
        fieldsFromSim();                  // fields travel with the tab's simulation
        view->homeScale = t.home; view->scale = t.scale; view->pan = t.pan;
        slZoom->blockSignals(true);
        slZoom->setValue(int(100.0 * (std::log10(t.scale) + 1.0)));
        slZoom->blockSignals(false);
        syncSettings();
        sim.computeAcc();
        static const char *names[4] = {"Galactic", "Solar", "Atomic", "Quark"};
        setWindowTitle(QString("Finite Theory simulator - %1 scale   [%2]")
                           .arg(names[k]).arg(profileName));
        view->update();
    }

    // Tab defaults: presets, dt, box, zoom.  Then the scene.
    void applyScale(int k)
    {
        syncSettings();
        switch (k) {
        case 0:  sbEtaG->setValue(7.4773e39);        // r_c = M_bulge/eta_g = 5 kpc
                 sbDt->setValue(1e30);  sbBox->setValue(15 * Sim::KPC_A0);
                 setZoom(250.0 / (15 * Sim::KPC_A0)); break;
        case 1:  sbEtaG->setValue(7.822540e46);      // c^2/G
                 sbDt->setValue(1e21);  sbBox->setValue(31 * Sim::AU_A0);
                 setZoom(250.0 / (6 * Sim::AU_A0)); break;   // inner system
        case 2:  cbEtaPreset->setCurrentIndex(4);    // visible on screen, ETA = 1
                 sbDt->setValue(0.02);  sbBox->setValue(5.0);
                 setZoom(40.0); break;
        default: cbEtaPreset->setCurrentIndex(1);    // retrofit
                 chConstituent->setChecked(true);
                 chQCD->setChecked(true);
                 if (cbKin->currentIndex() != 5) cbKin->setCurrentIndex(5);   // SR: limit c
                 sbDt->setValue(1e-12); sbBox->setValue(1e-4);
                 setZoom(150.0 / 1.6e-5); break;
        }
        restartScene();
    }

    // Scene only: rebuild the current tab's bodies with a new seed, keeping
    // every control, the time step and the view exactly as they are.
    void restartScene()
    {
        syncSettings();
        const int k = tabs->currentIndex();
        static const char *names[4] = {"Galactic", "Solar", "Atomic", "Quark"};
        switch (k) {
        case 0:  sim.buildGalaxy(1e10, 200, 15.0, seed++); break;
        case 1:  sim.buildSolar(seed++); break;
        case 2:  sim.reset(DEFAULT_NE, DEFAULT_NP, DEFAULT_NN, sbBox->value(), sbV->value(), seed++);
                 if (sim.quantised) sim.assignLevels(); break;
        default: sim.clear(); sim.addHadron(true, 1, 0.0, 0.0, seed++); break;
        }
        if (sim.spinOn) { sim.assignSpins(seed++); sim.computeAcc(); }
        setWindowTitle(QString("Finite Theory simulator - %1 scale   [%2]")
                           .arg(names[qBound(0, k, 3)]).arg(profileName));
        view->update();
    }

    void doReset() { applyScale(curTab); }   // defaults + scene for the current tab

    void tick()
    {
        const double dt = sbDt->value();
        for (int k = 0; k < stepsPerFrame; ++k) sim.step(dt);
        for (auto &a : sim.p) {
            a.trail.push_back(a.pos);
            if (a.trail.size() > 600) a.trail.pop_front();
        }
        view->update();
    }
};

#include "main.moc"

int main(int argc, char **argv)
{
    QApplication app(argc, argv);
    Window w;
    w.resize(1200, 760);
    w.show();
    return app.exec();
}
