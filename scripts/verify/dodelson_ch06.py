#!/usr/bin/env python
"""ตรวจขั้นพีชคณิตของคู่มืออ่าน Dodelson บทที่ 6 (reference/dodelson-ch06.html)

รัน:  python scripts/verify/dodelson_ch06.py
จุด (overdot) = d/d eta (conformal) · ตัวห้อย _t = d/dt · H = (da/dt)/a = a'/a^2
"""
import sympy as sp

ok = lambda name: print(f"OK  {name}")
zero = lambda e: sp.simplify(sp.expand(e)) == 0

eta = sp.symbols('eta')
k, G, Hc = sp.symbols('k G H', positive=True)

# ---------------------------------------------------------------- (6.4) a ∝ eta
rho = sp.symbols('rho', positive=True)
a_r = sp.symbols('a', positive=True)
# Friedmann ในเวลา conformal: (a'/a)^2 = 8πGρa^2/3 และ a'/a = 1/eta
sol = sp.solve(sp.Eq(1 / eta**2, 8 * sp.pi * G * rho * a_r**2 / 3), rho)[0]
assert zero(16 * sp.pi * G * sol * a_r**2 / 3 - 2 / eta**2)
ok("(6.3)-(6.4) 16πGρa²/3 = 2/eta² ในยุค radiation")

# ---------------------------------------------------------------- (6.6)-(6.9)
Phi = sp.Function('Phi')(eta)
fv, Th0, N0 = sp.symbols('f_nu Theta0 N0')
# d/deta ของ (6.6) โดยใช้ Theta0' = N0' = -Phi' และ Psi = -Phi
lhs = sp.diff(sp.diff(Phi, eta) * eta + Phi, eta)        # Phi'eta - Psi กับ Psi=-Phi
rhs = 2 * ((1 - fv) * (-sp.diff(Phi, eta)) + fv * (-sp.diff(Phi, eta)))
e68 = sp.expand(lhs - rhs)
assert zero(e68 - (sp.diff(Phi, eta, 2) * eta + 4 * sp.diff(Phi, eta)))
p = sp.symbols('p')
roots = sp.solve(sp.expand((eta**p * 0 + p * (p - 1) + 4 * p)), p)
assert set(roots) == {0, -3}
for pp in (0, -3):
    f = eta**pp
    assert zero(sp.diff(f, eta, 2) * eta + 4 * sp.diff(f, eta))
ok("(6.7)-(6.9) Phi''eta + 4Phi' = 0 ⇒ Phi ∝ eta^0, eta^-3")
# (6.10)-(6.12)
assert zero(2 * ((1 - fv) * Th0 + fv * Th0) - 2 * Th0)
ok("(6.10)-(6.12) Phi = 2 Theta0 เมื่อ Theta0 = N0")

# ---------------------------------------------------------------- (6.14)-(6.16)
eps = sp.symbols('epsilon')
d_, T0 = sp.symbols('delta Theta0')
r = sp.series((1 + eps * d_) / (1 + 3 * eps * T0), eps, 0, 2).removeO()
assert zero(r - (1 + eps * (d_ - 3 * T0)))
ok("(6.14)-(6.15) (1+δ)/(1+3Θ0) ≈ 1 + δ - 3Θ0")

# ---------------------------------------------------------------- (6.21)-(6.22)
t = sp.symbols('t')
a = sp.Function('a')(t)
rh, P = sp.symbols('rho P')
add = sp.symbols('addot')
H2 = 8 * sp.pi * G * rh / 3
add_sol = sp.solve(sp.Eq(add / a + H2 / 2, -4 * sp.pi * G * P), add)[0]
assert zero(add_sol / a + 4 * sp.pi * G / 3 * (rh + 3 * P))
ok("(6.21)-(6.22) a_tt/a = -(4πG/3)(ρ+3P)")

# ---------------------------------------------------------------- (6.25)-(6.27)
ph_t, V0 = sp.symbols('phidot_t V')
g = sp.diag(-1, a**2, a**2, a**2)
ginv = g.inv()
dphi = [ph_t, 0, 0, 0]
T = sp.Matrix(4, 4, lambda al, be: sum(ginv[al, n] * dphi[n] * dphi[be] for n in range(4))
              - (1 if al == be else 0) * (sp.Rational(1, 2) * sum(ginv[m, n] * dphi[m] * dphi[n]
                                                                  for m in range(4) for n in range(4)) + V0))
assert zero(-T[0, 0] - (ph_t**2 / 2 + V0))
assert zero(T[1, 1] - (ph_t**2 / 2 - V0))
ok("(6.24)-(6.27) ρ = φ_t²/2 + V และ P = φ_t²/2 - V")

# ---------------------------------------------------------------- (6.30)-(6.32)
phi = sp.Function('phi')(t)
Vf = sp.Function('V')
rho_phi = sp.diff(phi, t)**2 / 2 + Vf(phi)
P_phi = sp.diff(phi, t)**2 / 2 - Vf(phi)
H = sp.sqrt(8 * sp.pi * G * rho_phi / 3)
# continuity dρ/dt + 3H(ρ+P) = 0 เทียบเท่า φ_t (φ_tt + 3Hφ_t + V') = 0
cont = sp.diff(rho_phi, t) + 3 * H * (rho_phi + P_phi)
kg = sp.diff(phi, t, 2) + 3 * H * sp.diff(phi, t) + sp.diff(Vf(phi), phi)
assert zero(cont - sp.diff(phi, t) * kg)
ok("(6.30)-(6.32) φ_tt + 3Hφ_t + V' = 0 (เทียบเท่าการอนุรักษ์พลังงาน)")
# dH/dt = -4πG φ_t² (ใช้ต่อใน 6.76 และ Exercise 12)
dHdt = sp.diff(H, t).subs(sp.diff(phi, t, 2), sp.solve(kg, sp.diff(phi, t, 2))[0])
assert zero(dHdt + 4 * sp.pi * G * sp.diff(phi, t)**2)
ok("(6.21)+(6.32) dH/dt = -4πG φ_t²")

# ---------------------------------------------------------------- (6.33) conformal
ac = sp.Function('a')(eta)
ph = sp.Function('phi')(eta)
Vp = sp.Function('Vp')(eta)                  # V'(φ(eta))
Hc_ = sp.diff(ac, eta) / ac**2
ph_t_ = sp.diff(ph, eta) / ac
ph_tt_ = sp.diff(ph_t_, eta) / ac
e632 = ph_tt_ + 3 * Hc_ * ph_t_ + Vp
e633 = sp.diff(ph, eta, 2) + 2 * ac * Hc_ * sp.diff(ph, eta) + ac**2 * Vp
assert zero(e632 * ac**2 - e633)
ok("(6.32)->(6.33) Exercise 8")

# ---------------------------------------------------------------- (6.34)-(6.36), Exercise 12
aa = sp.symbols('a', positive=True)
ae = sp.symbols('a_e', positive=True)
eta_int = sp.integrate(1 / (Hc * aa**2), (aa, ae, sp.Symbol('A', positive=True)))
A_ = sp.Symbol('A', positive=True)
assert zero(eta_int - (1 / (Hc * ae) - 1 / (Hc * A_)))   # → -1/(aH) เมื่อ a_e >> a
ok("(6.34) eta = (1/H)(1/a_e - 1/a) ≈ -1/(aH)")
eps_c = -sp.diff(Hc_, eta) / (ac * Hc_**2)                 # (6.35)
assert zero(eps_c.subs(ac, eta).doit() - 2)                # radiation: a ∝ eta
assert zero(sp.simplify(eps_c.subs(ac, -1 / (Hc * eta)).doit()))   # de Sitter
assert zero(sp.diff(1 / (ac * Hc_), eta) - (eps_c - 1))   # Ex.12(a)
ok("(6.35) ε = 2 ในยุค radiation · ε = 0 ใน de Sitter · Ex.12(a) d(1/aH)/deta = ε-1")
delta_t = ph_tt_ / (Hc_ * ph_t_)                           # (6.36) บรรทัดแรก
delta_c = -1 / (ac * Hc_ * sp.diff(ph, eta)) * (ac * Hc_ * sp.diff(ph, eta) - sp.diff(ph, eta, 2))
assert zero(delta_t - delta_c)
delta_3 = delta_c.subs(sp.diff(ph, eta, 2), sp.solve(e633, sp.diff(ph, eta, 2))[0])
assert zero(delta_3 + 1 / (ac * Hc_ * sp.diff(ph, eta)) * (3 * ac * Hc_ * sp.diff(ph, eta) + ac**2 * Vp))
ok("(6.36) สามรูปของ δ เท่ากัน (ใช้ 6.33)")
# Ex.12(b),(c): ใช้ ε = 4πG φ_t²/H² (จาก dH/dt = -4πGφ_t²) แล้วหา dε/deta
eps_phi = 4 * sp.pi * G * ph_t_**2 / Hc_**2
assert zero(4 * sp.pi * G * sp.diff(ph, eta)**2 - eps_phi * ac**2 * Hc_**2)          # (6.112)
deps = sp.diff(eps_phi, eta)
plus = 2 * ac * Hc_ * eps_phi * (eps_c + delta_t)
# ใช้เงื่อนไข dH/dt = -4πGφ_t² เพื่อแทน ε ในรูป H ด้วยรูป φ
cond = sp.solve(sp.Eq(eps_c, eps_phi), sp.diff(ac, eta, 2))[0]
assert zero((deps - plus.subs(eps_c, eps_phi)).subs(sp.diff(ac, eta, 2), cond))
ok("Ex.12(c) dε/deta = +2aHε(ε+δ)  ⚠ หนังสือ (6.113) พิมพ์เครื่องหมายลบ — ดูหมายเหตุใน html")

# ---------------------------------------------------------------- (6.37)-(6.44)
w, tt = sp.symbols('omega t', positive=True)
v = sp.exp(-sp.I * w * tt) / sp.sqrt(2 * w)
assert zero(sp.diff(v, tt, 2) + w**2 * v)
comm = v * sp.conjugate(sp.diff(v, tt)) - sp.conjugate(v) * sp.diff(v, tt)   # [x,p] = v v̇* - v* v̇
assert zero(comm - sp.I)
assert zero(v * sp.conjugate(v) - 1 / (2 * w))
ok("(6.38)-(6.44) [x,p] = i เมื่อ v = e^{-iωt}/√(2ω) และ <|x|²> = 1/(2ω)")

# ---------------------------------------------------------------- (6.46)-(6.49)
ht = sp.Function('htilde')(eta)
h = sp.sqrt(16 * sp.pi * G) * ht / ac
e645 = sp.diff(h, eta, 2) + 2 * sp.diff(ac, eta) / ac * sp.diff(h, eta) + k**2 * h
e649 = (sp.diff(ht, eta, 2) + (k**2 - sp.diff(ac, eta, 2) / ac) * ht) / ac
assert zero(e645 / sp.sqrt(16 * sp.pi * G) - e649)
ok("(6.46)-(6.49) h̃ = ah/√(16πG) ทำให้ไม่มีพจน์หน่วง")

# ---------------------------------------------------------------- (6.55)-(6.59)
ads = -1 / (Hc * eta)                                        # de Sitter, eta < 0
assert zero(sp.diff(ads, eta, 2) / ads - 2 / eta**2)
ok("(6.55) a''/a = 2/eta²")
vv = sp.exp(-sp.I * k * eta) / sp.sqrt(2 * k) * (1 - sp.I / (k * eta))
assert zero(sp.diff(vv, eta, 2) + (k**2 - 2 / eta**2) * vv)
wr = vv * sp.conjugate(sp.diff(vv, eta)) - sp.conjugate(vv) * sp.diff(vv, eta)
etar = sp.symbols('eta_r', real=True)
assert zero(wr.subs(eta, etar) - sp.I)
assert sp.limit((vv * sp.sqrt(2 * k) * sp.exp(sp.I * k * eta)), eta, -sp.oo) == 1
ok("(6.56)-(6.57) v แก้สมการ · Wronskian = i (normalization แบบ 6.41) · → e^{-ik eta}/√(2k) เมื่อ |k eta|→∞")
v2 = sp.simplify(sp.expand(vv * sp.conjugate(vv)).subs(eta, etar))
small = sp.limit(v2 * etar**2, etar, 0)                       # |v|² ~ 1/(2k³eta²)
assert zero(small - 1 / (2 * k**3))
Ph = 16 * sp.pi * G * (1 / (2 * k**3 * etar**2)) / ads.subs(eta, etar)**2
assert zero(Ph - 8 * sp.pi * G * Hc**2 / k**3)
ok("(6.58)-(6.59) P_h = 16πG|v|²/a² → 8πGH²/k³")
Pdphi = (1 / (2 * k**3 * etar**2)) / ads.subs(eta, etar)**2
assert zero(Pdphi - Hc**2 / (2 * k**3))
ok("(6.72) P_δφ = H²/(2k³)")

# ---------------------------------------------------------------- (6.69)->(6.71)
dp = sp.Function('dphi')(eta)
Vpp = sp.Function('Vpp')(eta)                                # V''(φ)
# dV'/deta = V'' φ' : ใช้ Vp เป็นฟังก์ชันที่อนุพันธ์คือ Vpp*φ'
VpF = sp.Function('VpF')
subsV = {sp.Derivative(Vp, eta): Vpp * sp.diff(ph, eta)}
dT00 = -sp.diff(ph, eta) * sp.diff(dp, eta) / ac**2 - Vp * dp     # (6.67)
dTii = sp.diff(ph, eta) * sp.diff(dp, eta) / ac**2 - Vp * dp      # (6.68) ไม่รวมผลรวม
e669 = (sp.diff(dT00, eta) / ac + 3 * Hc_ * dT00) - k**2 / ac**3 * sp.diff(ph, eta) * dp - 3 * Hc_ * dTii
e669 = e669.subs(subsV)
e670 = (-sp.diff(ph, eta) * sp.diff(dp, eta, 2)
        + sp.diff(dp, eta) * (-sp.diff(ph, eta, 2) - 4 * ac * Hc_ * sp.diff(ph, eta) - ac**2 * Vp)
        + dp * (-ac**2 * Vpp * sp.diff(ph, eta) - k**2 * sp.diff(ph, eta)))
assert zero(e669 * ac**3 - e670)
ok("(6.69)->(6.70) คูณ a³ และ ∂V'/∂eta = V''φ'")
coef = (-sp.diff(ph, eta, 2) - 4 * ac * Hc_ * sp.diff(ph, eta) - ac**2 * Vp).subs(
    sp.diff(ph, eta, 2), sp.solve(e633, sp.diff(ph, eta, 2))[0])
assert zero(coef + 2 * ac * Hc_ * sp.diff(ph, eta))
e671 = (e670.subs(Vpp, 0).subs(sp.diff(ph, eta, 2), sp.solve(e633, sp.diff(ph, eta, 2))[0])) / (-sp.diff(ph, eta))
assert zero(e671 - (sp.diff(dp, eta, 2) + 2 * ac * Hc_ * sp.diff(dp, eta) + k**2 * dp))
ok("(6.70)->(6.71) ทิ้ง V'' แล้วได้สมการเดียวกับ (6.45)")

# ---------------------------------------------------------------- (6.76)
assert zero((ph_t_**2) / (3 * Hc_**2 / (8 * sp.pi * G)) - sp.Rational(2, 3) * eps_phi)
ok("(6.76) (P+ρ)/ρ = 2ε/3")

# ---------------------------------------------------------------- (6.77)-(6.79)
dph, phd, aHs, Th1, Ps, rr = sp.symbols('dphi phidot aH Theta1 Psi rho_r')
# scalar field: δT^0_i = g^00 φ_,0 φ_,i = -(phidot/a) i k_i δφ ⇒ i k_i δT^0_i = (phidot/a) k² δφ
a_s = sp.symbols('a_s', positive=True)
Hs = aHs / a_s
ikT = phd / a_s * k**2 * dph
zeta = -ikT * Hs / (k**2 * (phd / a_s)**2)                   # Psi ทิ้ง
assert zero(zeta + aHs * dph / phd)
ok("(6.78) ζ = -aH δφ/φ' ที่ horizon crossing")
ikT_r = 4 * a_s * k * rr * Th1
zeta_r = -ikT_r * Hs / (k**2 * sp.Rational(4, 3) * rr) - Ps
assert zero(zeta_r - (-3 * aHs * Th1 / k - Ps))
Phis = sp.symbols('Phi_s')
zeta_r2 = zeta_r.subs({Th1: -k * Phis / (6 * aHs), Ps: -Phis})
assert zero(zeta_r2 - sp.Rational(3, 2) * Phis)
ok("(6.79) ζ = -3aHΘ1/k - Ψ = -3Ψ/2 (ใช้ 6.16 และ Ψ = -Φ)")

# ---------------------------------------------------------------- (6.80)-(6.83), (6.95)-(6.97), (6.100)
eps_s, V_s, Vp_s = sp.symbols('epsilon V Vprime', positive=True)
Psi_post = sp.Rational(2, 3) * aHs * dph / phd
P_dphi = Hc**2 / (2 * k**3)
P_Psi = sp.Rational(4, 9) * (aHs / phd)**2 * P_dphi
ratio = 4 * sp.pi * G / eps_s                                 # (aH/φ')² = 4πG/ε  (Ex.12b)
P_Psi_e = P_Psi.subs((aHs / phd)**2, ratio)
P_Psi_e = sp.Rational(4, 9) * ratio * P_dphi
assert zero(P_Psi_e - 8 * sp.pi * G * Hc**2 / (9 * k**3 * eps_s))
ok("(6.80)-(6.82) P_Ψ = 8πGH²/(9k³ε)")
# Ex.12(b) ↔ ratio: φ' = aφ_t, ε = 4πGφ_t²/H²
assert zero((ac * Hc_ / (ac * ph_t_))**2 - 4 * sp.pi * G / eps_phi)
# (6.83) ด้วย slow-roll ε = V'²/(16πG V²), H² = 8πGV/3
P83 = (8 * sp.pi * G * Hc**2 / (9 * k**3 * eps_s)).subs(eps_s, Vp_s**2 / (16 * sp.pi * G * V_s**2))
assert zero(P83 - 128 * sp.pi**2 * G**2 / (9 * k**3) * Hc**2 * V_s**2 / Vp_s**2)
ok("(6.83) P_Φ = (128π²G²/9k³)(H²V²/V'²)")
P_zeta = ratio * P_dphi
assert zero(P_zeta - 2 * sp.pi * G * Hc**2 / (eps_s * k**3))
assert zero(sp.Rational(4, 9) * P_zeta - 8 * sp.pi * G * Hc**2 / (9 * eps_s * k**3))
ok("(6.95)-(6.97) P_ζ = 2πGH²/(εk³) และ P_Φ = 4P_ζ/9")
mpl = 1 / sp.sqrt(G)
assert zero(8 * sp.pi * Hc**2 / (9 * k**3 * eps_s * mpl**2) - 8 * sp.pi * G * Hc**2 / (9 * k**3 * eps_s))
assert zero(8 * sp.pi * Hc**2 / (k**3 * mpl**2) - 8 * sp.pi * G * Hc**2 / k**3)
ok("(6.100) รูปที่ใช้ m_Pl = G^{-1/2} ตรงกับ (6.82) และ (6.59)")
# tensor/scalar
assert zero((8 * sp.pi * G * Hc**2 / k**3) / (8 * sp.pi * G * Hc**2 / (9 * k**3 * eps_s)) - 9 * eps_s)
ok("P_h/P_Φ = 9ε (สัดส่วน ~ ε ตามที่หนังสือพูด)")

# ---------------------------------------------------------------- (6.101)-(6.106)
# ที่ horizon crossing k = aH ⇒ d ln k = d ln(aH) = aH(1-ε) d eta
e_, d2 = sp.symbols('e d')
dlnk_deta = aHs * (1 - e_)
dlnH_deta = -aHs * e_                         # Ḣ/H = -aHε
dlneps_deta = 2 * aHs * (e_ + d2)             # จาก Ex.12(c) ที่ sympy ได้
nT = 2 * sp.series(dlnH_deta / dlnk_deta, e_, 0, 2).removeO()
assert zero(sp.expand(nT).coeff(e_, 1) + 2)
n_minus_1 = sp.series((2 * dlnH_deta - dlneps_deta) / dlnk_deta, e_, 0, 2).removeO()
first = sp.expand(n_minus_1).subs({e_**2: 0, e_ * d2: 0})
assert zero(first - (-4 * e_ - 2 * d2))
ok("(6.104) n_T = -2ε · (6.106) n = 1 - 4ε - 2δ (อันดับหนึ่ง)")

# ---------------------------------------------------------------- Exercise 14
phs = sp.symbols('phi')
Vfun = sp.Function('V')(phs)
H2s = 8 * sp.pi * G * Vfun / 3                     # slow roll: ρ ≈ V
phidot_t = -sp.diff(Vfun, phs) / (3 * sp.sqrt(H2s))
eps14 = 4 * sp.pi * G * phidot_t**2 / H2s
assert zero(eps14 - sp.diff(Vfun, phs)**2 / (16 * sp.pi * G * Vfun**2))
# δ = φ_tt/(Hφ_t) โดยหา φ_tt จาก φ_t(φ) ด้วย chain rule
phi_tt = sp.diff(phidot_t, phs) * phidot_t
delta14 = phi_tt / (sp.sqrt(H2s) * phidot_t)
assert zero(delta14 - (eps14 - sp.diff(Vfun, phs, 2) / (8 * sp.pi * G * Vfun)))
ok("Exercise 14: ε = (V'/V)²/(16πG) · δ = ε - V''/(8πG V)")

# ---------------------------------------------------------------- (6.19), (6.29), e-folds
te, tb, ts = sp.symbols('t_e t_b s', real=True)
eta_prim = sp.integrate(1 / (ae * sp.exp(Hc * (ts - te))), (ts, tb, te))
assert zero(eta_prim - (sp.exp(Hc * (te - tb)) - 1) / (Hc * ae))
assert abs(float(sp.log(10**28)) - 64.47) < 0.01
ok("(6.29) eta_prim = (e^{H(te-tb)} - 1)/(H a_e) · ln(10^28) ≈ 64")

print("\nผ่านทั้งหมด — dodelson_ch06")
