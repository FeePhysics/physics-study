"""คู่มืออ่าน Dodelson บทที่ 8 (Anisotropies) — ตรวจขั้นพีชคณิตที่คู่มือเติมเอง

    python scripts/verify/dodelson_ch08.py
"""
import sympy as sp
import mpmath as mp

ok = 0


def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("✅", name)


def is_zero(e):
    return sp.simplify(sp.together(sp.expand(e))) == 0


eta, k, t, w, m = sp.symbols('eta k t omega m', positive=True)

# ── ① กล่องหน้า 220: (8.2) → (8.3) ────────────────────────────────────
Aa, F0, kk = sp.symbols('A F0 kappa', positive=True)   # kappa = force constant (หนังสือใช้ k)
om = sp.sqrt(kk / m)
x = Aa * sp.cos(om * t) + F0 / (m * om**2)
check("(8.3) แก้ (8.2): ẍ + (k/m)x = F0/m", is_zero(sp.diff(x, t, 2) + kk / m * x - F0 / m))
check("(8.3): ẋ(0) = 0", sp.diff(x, t).subs(t, 0) == 0)
# ยอดคี่สูงกว่ายอดคู่: x² ที่ ωt = π (คี่) เทียบ ωt = 2π (คู่) เมื่อ A < 0 (เริ่มต้นตรงข้ามแรง) — ตัวอย่างตัวเลข
xv = lambda wt: (-1 * sp.cos(wt) + sp.Rational(1, 2))**2
check("ตัวอย่าง A=−1, F0/mω²=1/2: x²(π)=9/4 > x²(2π)=1/4", xv(sp.pi) == sp.Rational(9, 4) and xv(2 * sp.pi) == sp.Rational(1, 4))

# ── ② §8.2: (8.5)-(8.8) ────────────────────────────────────────────────
Pp = sp.symbols('Phi_p')
Phis = sp.Rational(9, 10) * Pp
Th0 = -Phis + sp.Rational(3, 2) * Pp
check("(8.5): Θ0(η*) = 2Φ(η*)/3", is_zero(Th0 - sp.Rational(2, 3) * Phis))
Psi = -Phis
check("(8.6): Θ0 + Ψ = Ψ/3", is_zero(Th0 + Psi - Psi / 3))
dl = sp.Rational(3, 2) * Pp - 3 * (Phis - Pp)
check("(8.7): δ(η*) = 2Φ(η*)", is_zero(dl - 2 * Phis))
check("(8.8): Θ0 + Ψ = −δ/6", is_zero(Th0 + Psi + dl / 6))

# ── ③ (8.9)-(8.10): ทวน Legendre recurrence ด้วยการอินทิเกรตจริง ───────
mu = sp.symbols('mu')
c = sp.symbols('c0:9')
# inverse ของ (4.99): Θ(μ) = Σ (−i)^l (2l+1) P_l(μ) Θ_l
Theta = sum((-sp.I)**l * (2 * l + 1) * sp.legendre(l, mu) * c[l] for l in range(9))
Th_l = lambda l: sp.expand(sp.integrate(sp.legendre(l, mu) * Theta, (mu, -1, 1)) / 2 / (-sp.I)**l)
check("(4.99) กลับด้านถูก: (1/(−i)^l)∫dμ/2 P_l Θ = Θ_l (l=0..4)", all(is_zero(Th_l(l) - c[l]) for l in range(5)))
for l in (3, 4, 5):
    term = sp.expand(k / (-sp.I)**(l + 1) * sp.integrate(mu * sp.legendre(l, mu) * Theta, (mu, -1, 1)) / 2)
    target = -k * l / (2 * l + 1) * c[l - 1] + k * (l + 1) / (2 * l + 1) * c[l + 1]
    check(f"(8.9) ⇒ (8.10) ที่ l={l}", is_zero(term - target))

# ── ④ (8.13)+(8.15) ⇒ (8.16) ⇒ (8.18) ⇒ (8.20) ───────────────────────
f = lambda n: sp.Function(n)(eta)
T0, T1, Ph, Ps, a, td = f('Theta0'), f('Theta1'), f('Phi'), f('Psi'), f('a'), f('taudot')
r0 = sp.symbols('r0', positive=True)
R = r0 * a                                   # R ∝ a (4.108: R = 3ρ_b/4ρ_γ)
H = a.diff(eta) / a
vb = -3 * sp.I * T1 + R / td * (-3 * sp.I * T1.diff(eta) - 3 * sp.I * H * T1 + sp.I * k * Ps)   # (8.15)
e813 = T1.diff(eta) - k * T0 / 3 - k * Ps / 3 - td * (T1 - sp.I * vb / 3)
e816 = T1.diff(eta) + H * R / (1 + R) * T1 - k / (3 * (1 + R)) * T0 - k * Ps / 3
check("(8.13)+(8.15) ⇒ (8.16) (คูณด้วย 1+R)", is_zero(e813 - (1 + R) * e816))
T1_of = -(T0.diff(eta) + Ph.diff(eta)) / k                      # (8.12)
T1dot = sp.solve(sp.Eq(e816, 0), T1.diff(eta))[0]
e817 = (T0.diff(eta, 2) + k * T1dot + Ph.diff(eta, 2)).subs(T1, T1_of)
cs2 = 1 / (3 * (1 + R))
e818 = T0.diff(eta, 2) + H * R / (1 + R) * T0.diff(eta) + k**2 * cs2 * T0 \
    - (-k**2 / 3 * Ps - H * R / (1 + R) * Ph.diff(eta) - Ph.diff(eta, 2))
check("(8.12)+(8.16) ⇒ (8.18) ด้วย c_s² = 1/[3(1+R)]", is_zero(e817 - e818))
X = T0 + Ph
Rd = R.diff(eta)
e820 = X.diff(eta, 2) + Rd / (1 + R) * X.diff(eta) + k**2 * cs2 * X - k**2 / 3 * (Ph / (1 + R) - Ps)
check("(8.18) ⇔ (8.20) (ใช้ Ṙ/R = ȧ/a)", is_zero(e818 - e820))

# ── ⑤ (8.21)-(8.24): Green's function ─────────────────────────────────
F = sp.Function('F')        # F = Φ − Ψ
etap = sp.symbols('etap')
cs = 1 / sp.sqrt(3)
sol = (k / sp.sqrt(3)) * sp.Integral(F(etap) * sp.sin(k * cs * (eta - etap)), (etap, 0, eta))
check("พจน์อินทิกรัลใน (8.24) แก้ X'' + k²X/3 = (k²/3)(Φ−Ψ)",
      sp.simplify((sp.diff(sol, eta, 2) + k**2 / 3 * sol).doit() - k**2 / 3 * F(eta)) == 0)
C = sp.symbols('C')
hom = C * sp.cos(k * cs * eta)
check("C cos(kr_s) แก้สมการ homogeneous และ Ẋ(0)=0", is_zero(sp.diff(hom, eta, 2) + k**2 / 3 * hom)
      and sp.diff(hom, eta).subs(eta, 0) == 0)
S1, S2 = sp.sin(k * cs * etap), sp.cos(k * cs * etap)
Wr = sp.simplify(S1 * sp.diff(S2, etap) - sp.diff(S1, etap) * S2)
check("ตัวส่วนใน (8.23) = −k c_s = −k/√3", is_zero(Wr + k / sp.sqrt(3)))
num = sp.sin(k * cs * etap) * sp.cos(k * cs * eta) - sp.sin(k * cs * eta) * sp.cos(k * cs * etap)
check("ตัวเศษใน (8.23) = −sin[k(r_s − r_s')]", is_zero(sp.expand_trig(num + sp.sin(k * cs * (eta - etap)))))
# (8.26) = −(1/k) d/dη ของ (8.24)
full = hom + sol
dip = -sp.diff(full, eta).doit() / k
dip_book = C / sp.sqrt(3) * sp.sin(k * cs * eta) - k / 3 * sp.Integral(F(etap) * sp.cos(k * cs * (eta - etap)), (etap, 0, eta))
check("(8.12) + (8.24) ⇒ (8.26)", sp.simplify(dip - dip_book.doit()) == 0)

# ── ⑥ §8.4 diffusion damping: (8.33)-(8.40) ──────────────────────────
wR, tdd, Rr, kk2 = sp.symbols('omega taudot R k', nonzero=True)
eps = sp.symbols('epsilon')     # ตัวนับกำลังของ 1/τ̇
exact = 1 / (1 - sp.I * wR * Rr * eps)
ser = sp.series(exact, eps, 0, 3).removeO()
check("(8.33): [1 − iωR/τ̇]⁻¹ ≈ 1 + iωR/τ̇ − (ωR/τ̇)²", is_zero(ser - (1 + sp.I * wR * Rr * eps - (wR * Rr * eps)**2)))
Th1 = sp.symbols('Theta1')
Th2 = -4 * kk2 / (9 * tdd) * Th1                         # (8.34)
check("(8.29) ไม่มี Θ̇2 ⇒ (8.34)", is_zero(-2 * kk2 / 5 * Th1 - sp.Rational(9, 10) * tdd * Th2))
Th0v = -kk2 * Th1 / (sp.I * wR)                          # (8.35)
ivb3 = Th1 * (1 + sp.I * wR * Rr / tdd - (wR * Rr / tdd)**2)
e828 = sp.I * wR * Th1 + kk2 * (sp.Rational(2, 3) * Th2 - Th0v / 3) - tdd * (Th1 - ivb3)
e836 = (sp.I * wR - 8 * kk2**2 / (27 * tdd) + kk2**2 / (3 * sp.I * wR)) - tdd * (1 - (1 + sp.I * wR * Rr / tdd - (wR * Rr / tdd)**2))
check("(8.28) ⇒ (8.36)", is_zero(e828 / Th1 - e836))
e837 = wR**2 * (1 + Rr) - kk2**2 / 3 + sp.I * wR / tdd * (wR**2 * Rr**2 + 8 * kk2**2 / 27)
check("(8.36) × (−iω) ⇒ (8.37)", is_zero(sp.expand(e836 * (-sp.I * wR)) - sp.expand(e837)))
# δω: แทน ω = kc_s + δω, เก็บอันดับ 1/τ̇
dw, z = sp.symbols('delta_omega z')
csv = 1 / sp.sqrt(3 * (1 + Rr))
e = e837.subs(tdd, 1 / z).subs(wR, kk2 * csv + dw)
lin = sp.series(e, z, 0, 2).removeO()
lin = sp.expand(lin.subs(dw, sp.Symbol('D') * z))
dsol = sp.solve(sp.expand(lin).coeff(z, 1), sp.Symbol('D'))[0] * z
book838 = -sp.I * kk2**2 / (2 * (1 + Rr)) * z * (csv**2 * Rr**2 + sp.Rational(8, 27))
check("(8.37) ⇒ (8.38) δω อันดับแรกใน 1/τ̇", is_zero(dsol - book838))
# (8.39)-(8.40): i∫δω dη = −k²/k_D² เมื่อ τ̇ = −n_eσ_T a
nsa = sp.symbols('n_e_sigma_a', positive=True)
integrand = sp.I * book838.subs(z, -1 / nsa)
kD2inv = 1 / (6 * (1 + Rr) * nsa) * (Rr**2 / (1 + Rr) + sp.Rational(8, 9))
check("(8.38) ⇒ (8.39)-(8.40): iδω = −k²·[integrand ของ 1/k_D²]", is_zero(integrand + kk2**2 * kD2inv))
check("แบบฝึกหัด 8 ตั้งต้น: R = 0 ⇒ 1/k_D² = ∫ dη (8/9)/(6 n_eσ_T a)", is_zero(kD2inv.subs(Rr, 0) - sp.Rational(8, 54) / nsa))

# ── ⑦ (8.47): ∫dμ/2 P_l e^{ixμ} = j_l(x)/(−i)^l ──────────────────────
mp.mp.dps = 20
for l in range(5):
    for xv in (0.7, 3.3, 11.0):
        lhs = mp.quad(lambda u: mp.legendre(l, u) * mp.exp(1j * xv * u), [-1, 1]) / 2
        rhs = mp.sqrt(mp.pi / (2 * xv)) * mp.besselj(l + 0.5, xv) / (-1j)**l
        assert abs(lhs - rhs) < 1e-12
check("(8.47) ตรวจตัวเลข l=0..4, x=0.7, 3.3, 11", True)
check("j_l(−x) = (−1)^l j_l(x) (ใช้ใน 8.51)", all(is_zero(sp.expand_func(sp.jn(l, -sp.Symbol('x'))) - (-1)**l * sp.expand_func(sp.jn(l, sp.Symbol('x')))) for l in range(5)))

# ── ⑧ (8.53): ∫ g dη = 1 ───────────────────────────────────────────
tau = sp.Function('tau')(eta)
g = -tau.diff(eta) * sp.exp(-tau)
check("(8.53): g = d(e^{−τ})/dη ⇒ ∫_0^{η0} g = e^{−τ(η0)} − e^{−τ(0)} = 1 − 0", is_zero(g - sp.diff(sp.exp(-tau), eta)))

# ── ⑨ (8.64) cosmic variance ──────────────────────────────────────────
Cl, l_ = sp.symbols('C_l l', positive=True)
N = 2 * l_ + 1
# Ĉ = (1/N)Σa², a ~ Gaussian var C: E[a⁴] = 3C²
var = (3 * Cl**2 - Cl**2) / N
check("(8.64): Var(Ĉ_l)/C_l² = 2/(2l+1)", is_zero(var / Cl**2 - 2 / N))
check("COBE: 41,253 deg²/(7°)² ≈ 840 ⇒ (l_max+1)² = 840 ⇒ l_max ≈ 28",
      abs(4 * mp.pi * (180 / mp.pi)**2 / 49 - 842) < 2 and int(mp.sqrt(840)) - 1 == 27)

# ── ⑩ Sachs-Wolfe (8.70)-(8.75) ─────────────────────────────────────
Om, H0, D1, d0, kq, dH, n, P = sp.symbols('Omega_m H0 D1 delta0 k delta_H n P', positive=True)
Phi0 = sp.Rational(3, 2) * Om * H0**2 * d0 / kq**2          # (7.7) ที่ a=1
check("(8.69)+(7.7) ⇒ (8.70)", is_zero(-Phi0 / (3 * D1) + Om * H0**2 / (2 * kq**2 * D1) * d0))
pref = sp.Rational(2, 1) / sp.pi * kq**2 * (Om * H0**2 / (2 * kq**2 * D1))**2
check("(8.68)+(8.70) ⇒ (8.71): ค่าหน้าอินทิกรัล = Ω_m²H0⁴/(2πD1²k²)", is_zero(pref - Om**2 * H0**4 / (2 * sp.pi * D1**2 * kq**2)))
Pk = 2 * sp.pi**2 * dH**2 * kq**n / H0**(n + 3)                # (7.9) ที่ a=1, T=1
check("(8.71)+(7.9) ⇒ (8.72)", is_zero(sp.simplify(pref * Pk - sp.pi * H0**(1 - n) * (Om / D1)**2 * dH**2 * kq**(n - 2))))
for l in (2, 5, 10):
    fI = lambda x: (mp.sqrt(mp.pi / (2 * x)) * mp.besselj(l + 0.5, x))**2 / x
    I1 = mp.quad(fI, [0, l, 2 * l]) + mp.quadosc(fI, [2 * l, mp.inf], period=mp.pi)
    check(f"n=1: ∫dx j_l²/x = 1/[2l(l+1)] ที่ l={l}", abs(I1 * 2 * l * (l + 1) - 1) < 1e-5)
nn = mp.mpf('0.6')
for l in (3, 8):
    fn = lambda x: x**(nn - 2) * (mp.sqrt(mp.pi / (2 * x)) * mp.besselj(l + 0.5, x))**2
    In = mp.quad(fn, [0, l, 2 * l]) + mp.quadosc(fn, [2 * l, mp.inf], period=mp.pi)
    G = 2**(nn - 4) * mp.pi * mp.gamma(l + nn / 2 - 0.5) * mp.gamma(3 - nn) / (mp.gamma(l + 2.5 - nn / 2) * mp.gamma(2 - nn / 2)**2)
    check(f"(8.74) สูตร gamma ตรงกับอินทิกรัลตัวเลข (n=0.6, l={l})", abs(In / G - 1) < 1e-5)
check("Γ(2)/Γ(3/2)² = 4/π", sp.gamma(2) / sp.gamma(sp.Rational(3, 2))**2 == 4 / sp.pi)
lpos = sp.symbols('l', positive=True, integer=True)
check("Γ(l)/Γ(l+2) = 1/[l(l+1)]", is_zero(sp.gammasimp(sp.gamma(lpos) / sp.gamma(lpos + 2)) - 1 / (lpos * (lpos + 1))))
ClSW = 2**(1 - 4) * sp.pi**2 * (Om / D1)**2 * dH**2 * sp.gamma(lpos) * sp.gamma(2) / (sp.gamma(lpos + 2) * sp.gamma(sp.Rational(3, 2))**2)
check("(8.74), n=1 ⇒ (8.75): l(l+1)C_l = (π/2)(Ω_m/D1)²δ_H²",
      is_zero(sp.gammasimp(lpos * (lpos + 1) * ClSW) - sp.pi / 2 * (Om / D1)**2 * dH**2))

# ── ⑪ ตัวเลขเล็ก ๆ ในเนื้อความ ────────────────────────────────────────
check("dipole 30% บวกแบบ incoherent: 1² + 0.3² ≈ 1.1", abs(1 + 0.3**2 - 1.09) < 1e-12)
check("ISW 5% บวกแบบ coherent: (1.05)² ≈ 1.1", abs(1.05**2 - 1.1025) < 1e-12)
Tb, Th, tau_ = sp.symbols('T Theta tau')
check("(8.79): T(1+Θ)e^{−τ} + T(1−e^{−τ}) = T(1 + Θe^{−τ})",
      is_zero(Tb * (1 + Th) * sp.exp(-tau_) + Tb * (1 - sp.exp(-tau_)) - Tb * (1 + Th * sp.exp(-tau_))))

print(f"\nผ่านทั้งหมด {ok} ข้อ")
