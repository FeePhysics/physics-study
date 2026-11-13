#!/usr/bin/env python
"""ตรวจขั้นพีชคณิตของคู่มืออ่าน Dodelson บทที่ 4 (reference/dodelson-ch04.html)

แต่ละ check พิมพ์ OK หรือ assert ล้ม · รัน:  python scripts/verify/dodelson_ch04.py
สมการอ้างเลขของหนังสือ (Dodelson 2003) · ใช้ metric (4.9) signature (-+++)
"""
import sympy as sp

ok = lambda name: print(f"OK  {name}")

# ---------------------------------------------------------------- (4.8)
x, p, m, k = sp.symbols('x p m k', positive=True)
F = sp.Function('F')
E = p**2 / (2 * m) + k * x**2 / 2
expr = p / m * sp.diff(F(E), x) - k * x * sp.diff(F(E), p)
assert sp.simplify(expr) == 0
ok("(4.8) f = F(E) ทำให้ (p/m) df/dx - kx df/dp = 0")

# ---------------------------------------------------------------- (4.14)
eps = sp.symbols('epsilon')
Psi0, Phi0 = sp.symbols('Psi Phi')
P0 = p / sp.sqrt(1 + 2 * eps * Psi0)
assert sp.simplify(-(1 + 2 * eps * Psi0) * P0**2 + p**2) == 0
assert sp.simplify(sp.series(P0, eps, 0, 2).removeO() - p * (1 - eps * Psi0)) == 0
ok("(4.12)-(4.14) P^0 = p/sqrt(1+2Psi) ~ p(1-Psi)")

# ---------------------------------------------------------------- (4.18)-(4.20)
a = sp.symbols('a', positive=True)
C = p / (a * sp.sqrt(1 + 2 * eps * Phi0))
assert sp.simplify(a**2 * (1 + 2 * eps * Phi0) * C**2 - p**2) == 0
assert sp.simplify(sp.series(C, eps, 0, 2).removeO() - p * (1 - eps * Phi0) / a) == 0
ratio = (p * (1 - eps * Phi0) / a) / (p * (1 - eps * Psi0))       # P^i/P^0 ต่อ p^i
lin = sp.series(ratio, eps, 0, 2).removeO()
assert sp.simplify(lin - (1 + eps * Psi0 - eps * Phi0) / a) == 0
ok("(4.18)-(4.20) C = p(1-Phi)/a และ dx^i/dt = (p^i/a)(1+Psi-Phi)")

# ---------------------------------------------------------------- (4.22)-(4.32) และ (4.68)
# geodesic สำหรับอนุภาคมีมวล (massless = กรณี m->0) ใน metric (4.9)
t, X, Y, Z = sp.symbols('t x y z')
crd = [t, X, Y, Z]
af = sp.Function('a')(t)
Psi = sp.Function('Psi')(t, X, Y, Z)
Phi = sp.Function('Phi')(t, X, Y, Z)
g = sp.diag(-(1 + 2 * eps * Psi), *[af**2 * (1 + 2 * eps * Phi)] * 3)
ginv = sp.diag(1 / g[0, 0], *[1 / g[i, i] for i in (1, 2, 3)])   # metric ทแยง


def Gamma(mu, al, be):
    return sp.Rational(1, 2) * sum(
        ginv[mu, nu] * (sp.diff(g[nu, al], crd[be]) + sp.diff(g[nu, be], crd[al])
                        - sp.diff(g[al, be], crd[nu])) for nu in range(4))


n1, n2, n3, pm, mm = sp.symbols('n1 n2 n3 p m', positive=True)
nvec = [n1, n2, n3]
Eng = sp.sqrt(pm**2 + mm**2)
Pup = [Eng * (1 - eps * Psi)] + [pm * ni * (1 - eps * Phi) / af for ni in nvec]   # (4.66)


def first_order(e):
    return sp.series(sp.expand(e), eps, 0, 2).removeO()


# ตรวจว่า P^mu ที่เลือกสอดคล้องกับ constraint (4.64) ถึงอันดับหนึ่ง
norm = sum(g[i, i] * Pup[i]**2 for i in range(4))
norm1 = first_order(norm).subs(n3**2, 1 - n1**2 - n2**2)
assert sp.simplify(norm1 + mm**2) == 0
ok("(4.64)/(4.66) g_{mu nu}P^mu P^nu = -m^2 ถึงอันดับหนึ่ง")

GPP = sum(Gamma(0, i, j) * Pup[i] * Pup[j] for i in range(4) for j in range(4))
# d/dt [E(1-Psi)] = -Gamma^0 PP / P^0  (เพราะ d/dlambda = P^0 d/dt)
dxdt = [pm * ni / (af * Eng) for ni in nvec]          # อันดับศูนย์ พอแล้วเพราะคูณกับ eps
dPsidt = sp.diff(Psi, t) + sum(dxdt[i] * sp.diff(Psi, crd[i + 1]) for i in range(3))
dEdt = (-GPP / Pup[0] + eps * Eng * dPsidt) / (1 - eps * Psi)
dEdt1 = first_order(dEdt).subs(n3**2, 1 - n1**2 - n2**2)
H = sp.diff(af, t) / af
target = -(H * pm**2 / Eng + eps * pm**2 / Eng * sp.diff(Phi, t)
           + eps * pm / af * sum(nvec[i] * sp.diff(Psi, crd[i + 1]) for i in range(3)))
unit = {n3: sp.sqrt(1 - n1**2 - n2**2)}           # p-hat เป็นเวกเตอร์หนึ่งหน่วย
assert sp.simplify(sp.expand((dEdt1 - target).subs(unit))) == 0
ok("(4.68) dE/dt = -[H p^2/E + (p^2/E) dPhi/dt + (p^i p/a) dPsi/dx^i]")
# massless limit -> (4.32)
dpdt0 = sp.simplify(target.subs(mm, 0) / pm)
target32 = -(H + eps * sp.diff(Phi, t) + eps / af * sum(nvec[i] * sp.diff(Psi, crd[i + 1]) for i in range(3)))
assert sp.simplify(dpdt0 - target32) == 0
ok("(4.32) (1/p) dp/dt = -H - dPhi/dt - (p^i/a) dPsi/dx^i  (กรณี m = 0)")

# ---------------------------------------------------------------- (4.35)
T, Th = sp.symbols('T Theta', positive=True)
f0 = lambda TT: 1 / (sp.exp(p / TT) - 1)
lhs = sp.series(f0(T * (1 + eps * Th)), eps, 0, 2).removeO()
rhs = f0(T) - eps * p * sp.diff(f0(T), p) * Th
assert sp.simplify(lhs - rhs) == 0
assert sp.simplify(T * sp.diff(f0(T), T) + p * sp.diff(f0(T), p)) == 0
ok("(4.35) f ~ f0 - p (df0/dp) Theta และ T df0/dT = -p df0/dp")

# ---------------------------------------------------------------- (4.37)-(4.39)
Tt = sp.Function('T')(t)
z0 = sp.diff(f0(Tt), t) - H * p * sp.diff(f0(Tt), p)
sol = sp.simplify(z0.subs(Tt, 1 / af).doit())
assert sol == 0
ok("(4.37)-(4.39) T ∝ 1/a ทำให้สมการอันดับศูนย์เป็นศูนย์")

# ---------------------------------------------------------------- (4.40)->(4.42)
# ตรวจว่าพจน์ที่ 3 ของ (4.40) หักล้างกับชิ้นที่สองของ (4.41) เมื่อ dT/dt / T = -H
Thf = sp.Function('Theta')(t)
F0 = f0(Tt)
term1 = -p * sp.diff(sp.diff(F0, p) * Thf, t)
term3 = H * p * Thf * sp.diff(p * sp.diff(F0, p), p)
left = term1 + term3
want = -p * sp.diff(F0, p) * sp.diff(Thf, t)
diff = sp.simplify((left - want).subs(Tt, 1 / af).doit())
assert diff == 0
ok("(4.40)-(4.42) พจน์ H หักล้างกันพอดีเมื่อ T ∝ 1/a")

# ---------------------------------------------------------------- (4.47)
qv = sp.Matrix(sp.symbols('q1:4'))
pv = sp.Matrix(sp.symbols('p1:4'))
pp = sp.Matrix(sp.symbols('pp1:4'))
me = sp.symbols('m_e', positive=True)
dE = (qv.dot(qv) - (qv + pv - pp).dot(qv + pv - pp)) / (2 * me)
lead = (pp - pv).dot(qv) / me
rest = sp.expand(dE - lead)
assert sp.simplify(rest + (pv - pp).dot(pv - pp) / (2 * me)) == 0
ok("(4.47) Ee(q)-Ee(q+p-p') = (p'-p).q/me - |p-p'|^2/(2me)  (ตัวหลังเล็กเมื่อ q >> p)")

# ---------------------------------------------------------------- (4.52)/(4.97)/(4.99)
mu = sp.symbols('mu')
Leg = lambda l: sp.legendre(l, mu)
Th_mu = sp.Function('Theta')(mu)
moment = lambda l: sp.Integral(Leg(l) * Th_mu, (mu, -1, 1)) / 2 / (-sp.I)**l
assert sp.simplify(moment(1) - sp.I * sp.Integral(mu * Th_mu, (mu, -1, 1)) / 2) == 0
assert sp.expand(Leg(2) - (3 * mu**2 - 1) / 2) == 0
ok("(4.97)/(4.99) Theta_1 = i ∫dmu/2 mu Theta  และ P_2 = (3mu^2-1)/2")

# ---------------------------------------------------------------- (4.96)
Th0, vb = sp.symbols('Theta0 v_b')
assert sp.integrate(mu * (Th0 + vb * mu), (mu, -1, 1)) / 2 == vb / 3
# ∫ dp p^4 df/dp = -4 ∫ dp p^3 f  (ทดลองกับ f0 ของ Bose-Einstein ที่ T=1 ด้วยตัวเลข)
fb = 1 / (sp.exp(p) - 1)
I1 = sp.Integral(p**4 * sp.diff(fb, p), (p, 0, sp.oo)).evalf(30)
I2 = sp.Integral(p**3 * fb, (p, 0, sp.oo)).evalf(30)
assert abs(I1 + 4 * I2) < 1e-20
assert abs(I2 - sp.pi**4 / 15) < 1e-20
ok("(4.96) ∫dmu/2 mu(Theta0+vb mu) = vb/3  และ ∫p^4 f' = -4∫p^3 f")

# ---------------------------------------------------------------- (4.72)/(4.79)/(4.80)
fg = sp.exp(-p**2)
assert sp.simplify(sp.integrate(p**3 * sp.diff(fg, p), (p, 0, sp.oo))
                   + 3 * sp.integrate(p**2 * fg, (p, 0, sp.oo))) == 0
Ev = sp.sqrt(p**2 + m**2)
assert sp.simplify(sp.diff(p**4 / Ev, p) - (4 * p**3 / Ev - p**5 / Ev**3)) == 0
th, ph = sp.symbols('theta phi')
nn = [sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)]
for i in range(3):
    for j in range(3):
        val = sp.integrate(sp.integrate(nn[i] * nn[j] * sp.sin(th), (th, 0, sp.pi)), (ph, 0, 2 * sp.pi))
        assert sp.simplify(val - (4 * sp.pi / 3 if i == j else 0)) == 0
ok("(4.72) ∫p^3 f' = -3∫p^2 f · (4.79) d(p^4/E)/dp · (4.80) ∫dΩ p^i p^j = 4π/3 δ^ij")

# ---------------------------------------------------------------- (4.75)
n0 = sp.Function('n')(t)
sol = sp.dsolve(sp.Eq(sp.diff(n0, t) + 3 * sp.diff(af, t) / af * n0, 0), n0)
assert sp.simplify(sol.rhs * af**3).free_symbols <= {sp.Symbol('C1')}
ok("(4.74)-(4.75) n ∝ a^-3")

# ---------------------------------------------------------------- (4.56)->(4.63)
# ∂_t = (1/a)∂_eta · ∂_i -> i k_i · p^i k_i = k mu · -tau_dot = n_e sigma_T a
kk, a_, ne_sT, Th_, Phd, Psi_, Th0_, vb_ = sp.symbols('k a nsT Theta Phidot Psi Theta0 vb')
Thd = sp.symbols('Thetadot')
eq456 = Thd + sp.I * kk * mu * Th_ + Phd + sp.I * kk * mu * Psi_ - ne_sT * a_ * (Th0_ - Th_ + mu * vb_)
taud = -ne_sT * a_
eq463 = Thd + sp.I * kk * mu * Th_ + Phd + sp.I * kk * mu * Psi_ + taud * (Th0_ - Th_ + mu * vb_)
assert sp.simplify(eq456 - eq463) == 0
ok("(4.56)->(4.63) แทน d/dx^i -> ik^i และ tau_dot = -n_e sigma_T a")

# ---------------------------------------------------------------- (4.96)->(4.98)->(4.106)
# ใช้ ∫dmu/2 mu Theta(mu) = -i Theta_1 (จากนิยาม 4.97) และผลอินทิเกรต p = -4 rho_gamma (ตามหนังสือ)
Th1, rg, rb, nsT = sp.symbols('Theta1 rho_gamma rho_b nsT')
mu_int = 0 - (-sp.I * Th1) + vb / 3               # ∫dmu/2 mu [Theta0 - Theta(mu) + vb mu]
rhs_t = nsT / rb * (-4 * rg) * mu_int              # (4.96) ในเวลา t
rhs_eta = a_ * rhs_t                               # คูณ a เมื่อเปลี่ยนเป็น conformal time
taud = -nsT * a_
assert sp.simplify(rhs_eta - taud * 4 * rg / (3 * rb) * (3 * sp.I * Th1 + vb)) == 0
R = 3 * rb / (4 * rg)                              # (4.108)
assert sp.simplify(rhs_eta - taud / R * (vb + 3 * sp.I * Th1)) == 0
ok("(4.96)->(4.98)->(4.106) พจน์ Compton ของ baryon = (tau_dot/R)[vb + 3i Theta1]")

# ---------------------------------------------------------------- (4.77)/(4.82) -> (4.83)/(4.84)
# v^i = (k^i/k) v · d/dx^i -> i k_i · คูณ (4.82) ด้วย k^j/k
k1, k2, k3 = sp.symbols('k1:4', real=True)
kv = sp.Matrix([k1, k2, k3]); kn = sp.sqrt(k1**2 + k2**2 + k3**2)
v_, Ps_ = sp.symbols('v Psi')
div_v = sum(sp.I * kv[i] * kv[i] / kn * v_ for i in range(3))
assert sp.simplify(div_v - sp.I * kn * v_) == 0
grad_proj = sum(kv[j] / kn * sp.I * kv[j] * Ps_ for j in range(3))
assert sp.simplify(grad_proj - sp.I * kn * Ps_) == 0
ok("(4.83)/(4.84) div v -> ikv และ k-hat . grad Psi -> ik Psi")

print("\nผ่านทั้งหมด — dodelson_ch04")
