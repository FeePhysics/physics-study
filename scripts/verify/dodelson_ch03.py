"""คู่มืออ่าน Dodelson บทที่ 3 — ตรวจขั้นที่คู่มือเติมเอง (reference/dodelson-ch03.html)

รัน: python scripts/verify/dodelson_ch03.py
หน่วยธรรมชาติ ħ = c = k_B = 1 · แปลงเป็นวินาที/เซนติเมตรด้วยค่าคงตัวที่กำกับว่า "นอกหนังสือ"
"""
import sys
import sympy as sp
import mpmath as mp

sys.stdout.reconfigure(encoding="utf-8")
ok = 0


def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("✅", name)


# ค่าคงตัว (นอกหนังสือ เว้นแต่ระบุ)
hbar = 6.5821e-16          # eV s
hbarc = 1.97327e-5         # eV cm
mpl = 1.2209e28            # eV · G = 1/mpl² (หนังสือใช้ 1.2e19 GeV ในโจทย์ 2.1)
me, mprot = 0.511e6, 938.27e6   # eV
T0 = 2.725 / 11605         # eV (หนังสือ: 11605 K = 1 eV)
rcr_eV4 = 8.098e-11        # h² eV⁴ (หนังสือ 2.70)

E, p, m, T = sp.symbols('E p m T', positive=True)

# ① (3.2): ∫dE δ(E² − p² − m²) = 1/(2E) ที่ E = √(p²+m²)
Ep = sp.sqrt(p**2 + m**2)
fprime = sp.diff(E**2 - p**2 - m**2, E).subs(E, Ep)
check("δ(f(E)) = δ(E−E₀)/|f'(E₀)| กับ f' = 2E ⇒ ได้ 1/(2E)", sp.simplify(1 / fprime - 1 / (2 * Ep)) == 0)

# ② (3.4): f = e^{μ/T} e^{−E/T} และ E1+E2 = E3+E4
E1, E2, E3, mu1, mu2, mu3, mu4 = sp.symbols('E1 E2 E3 mu1 mu2 mu3 mu4', real=True)
E4 = E1 + E2 - E3
f = lambda mu, e: sp.exp(mu / T) * sp.exp(-e / T)
lhs = f(mu3, E3) * f(mu4, E4) - f(mu1, E1) * f(mu2, E2)
rhs = sp.exp(-(E1 + E2) / T) * (sp.exp((mu3 + mu4) / T) - sp.exp((mu1 + mu2) / T))
check("f3f4 − f1f2 = e^{−(E1+E2)/T}{e^{(μ3+μ4)/T} − e^{(μ1+μ2)/T}}  (3.4)", sp.simplify(lhs - rhs) == 0)

# ③ (3.6) สองลิมิต
g = sp.symbols('g', positive=True)
n_nr = g * 4 * sp.pi / (2 * sp.pi)**3 * sp.integrate(p**2 * sp.exp(-(m + p**2 / (2 * m)) / T), (p, 0, sp.oo))
check("m ≫ T: g∫d³p/(2π)³ e^{−(m+p²/2m)/T} = g(mT/2π)^{3/2} e^{−m/T}",
      sp.simplify(n_nr - g * (m * T / (2 * sp.pi))**sp.Rational(3, 2) * sp.exp(-m / T)) == 0)
n_r = g * 4 * sp.pi / (2 * sp.pi)**3 * sp.integrate(p**2 * sp.exp(-p / T), (p, 0, sp.oo))
check("m ≪ T: g∫d³p/(2π)³ e^{−p/T} = g T³/π²", sp.simplify(n_r - g * T**3 / sp.pi**2) == 0)

# ④ (3.15)→(3.16) deuterium
mn, mp_, mD, B = sp.symbols('m_n m_p m_D B', positive=True)
n0 = lambda gg, mm: gg * (mm * T / (2 * sp.pi))**sp.Rational(3, 2) * sp.exp(-mm / T)
r315 = sp.simplify(n0(3, mD) / (n0(2, mn) * n0(2, mp_)))
book315 = sp.Rational(3, 4) * (2 * sp.pi * mD / (mn * mp_ * T))**sp.Rational(3, 2) * sp.exp((mn + mp_ - mD) / T)
check("n_D⁰/(n_n⁰n_p⁰) = (3/4)(2πm_D/(m_n m_p T))^{3/2} e^{(m_n+m_p−m_D)/T}  (3.15)",
      sp.simplify(r315 / book315) == 1)
pref = (2 * sp.pi * mD / (mn * mp_ * T)).subs({mD: 2 * mp_, mn: mp_})
check("m_D = 2m_p = 2m_n ⇒ (4π/(m_p T))^{3/2}  (3.16)", sp.simplify(pref - 4 * sp.pi / (mp_ * T)) == 0)
check("B_D = 1877.84 − 1875.62 = 2.22 MeV  (3.12)", abs((1877.84 - 1875.62) - 2.22) < 1e-9)

# ⑤ (3.18)–(3.21)
check("(m_p/m_n)^{3/2} = 0.998 ≈ 1 (นอกหนังสือ: m_n = 939.57 MeV)", abs((938.27 / 939.57)**1.5 - 0.998) < 1e-3)
Q = sp.symbols('Q', positive=True)
XEQ = 1 / (1 + sp.exp(Q / T))
check("X_n,EQ = 1/(1+e^{Q/T}) → 1/2 เมื่อ T ≫ Q", sp.limit(XEQ, T, sp.oo) == sp.Rational(1, 2))

# ⑥ (3.22)→(3.23)→(3.27)
t = sp.symbols('t')
lam, X, x, H1 = sp.symbols('lambda X x H1', positive=True)
rhs323 = lam * ((1 - X) * sp.exp(-Q / T) - X)
check("(3.22) หาร (n_n+n_p): λ{ n_p e^{−Q/T}/(n_n+n_p) − X } = λ{(1−X)e^{−Q/T} − X}",
      sp.simplify(lam * ((1 - X) * sp.exp(-Q / T) - X) - rhs323) == 0)
# dX/dx = (dX/dt)/(dx/dt), dx/dt = xH, H = H1/x²
dXdx = rhs323.subs(T, Q / x) / (x * H1 / x**2)
book327 = x * lam / H1 * (sp.exp(-x) - X * (1 + sp.exp(-x)))
check("dx/dt = xH และ H = H(x=1)/x² ⇒ (3.27)", sp.simplify(dXdx - book327) == 0)
# H(T) ในยุค radiation จาก (2.39)+(3.26)
Gn, gs = sp.symbols('G g_*', positive=True)
HT = sp.sqrt(8 * sp.pi * Gn / 3 * gs * sp.pi**2 / 30 * T**4)
check("H = √(8πGρ/3), ρ = g*π²T⁴/30 ⇒ H = √(4π³Gg*/45) T²",
      sp.simplify(HT - sp.sqrt(4 * sp.pi**3 * Gn * gs / 45) * T**2) == 0)
Qv = 1.293e6
H1v = mp.sqrt(4 * mp.pi**3 / 45 * 10.75) * Qv**2 / mpl / hbar
check(f"H(x=1) = {float(H1v):.3f} s⁻¹ ≈ 1.13  (3.28)", abs(H1v - 1.13) < 0.01)
check("g* = 2 + (7/8)(6+2+2) = 10.75", 2 + 7 / 8 * (6 + 2 + 2) == 10.75)
check(f"λ_np(x=1) = (255/886.7)·19 = {255/886.7*19:.2f} ≈ 5.5 s⁻¹  (3.29)", abs(255 / 886.7 * 19 - 5.5) < 0.05)
# (3.30) หลัง e± หมด: g* = 3.36 (นี่คือ 2 + (7/8)·6·(4/11)^{4/3}) · RD ⇒ t = 1/(2H)
gs2 = 2 + 7 / 8 * 6 * (4 / 11)**(4 / 3)
check(f"g* หลัง e± หมด = 2 + (7/8)·6·(4/11)^(4/3) = {gs2:.2f}", abs(gs2 - 3.36) < 0.01)
Tv = 0.1e6
tv = 1 / (2 * mp.sqrt(4 * mp.pi**3 / 45 * 3.36) * Tv**2 / mpl) * hbar
check(f"t(T=0.1 MeV) = 1/(2H) = {float(tv):.0f} s ≈ 132  (3.30)", abs(tv - 132) < 2)
dec = mp.e**(-(132 / 886.7) * (0.1 / 0.07)**2)
check(f"exp[−(132/886.7)(0.1/0.07)²] = {float(dec):.3f} ≈ 0.74", abs(dec - 0.74) < 0.005)
check("X_n(T_nuc) = 0.15 × 0.74 ≈ 0.11  (3.31)", abs(0.15 * 0.74 - 0.111) < 1e-3)
check("X_4 = 2 X_n = 0.22  (3.33)", abs(2 * 0.11 - 0.22) < 1e-12)

# ⑦ (3.11) η_b และ (3.32) T_nuc และ (3.34)
nb = 0.020 * 1.878e-29 / 1.6726e-24          # cm⁻³ ที่ Ω_b h² = 0.02 (m_p = 1.6726e-24 g นอกหนังสือ)
ng = float(2 * mp.zeta(3) / mp.pi**2 * (T0 / hbarc)**3)
check(f"η_b = n_b/n_γ = {nb/ng:.2e} ≈ 5.5e-10 ที่ Ω_b h² = 0.020  (3.11)", abs(nb / ng / 5.5e-10 - 1) < 0.02)
Tn = mp.findroot(lambda TT: mp.log(5.5e-10) + 1.5 * mp.log(TT / 938.27) + 2.22 / TT, (0.04, 0.1), solver='bisect')
check(f"(3.32) แก้ ~ เป็น = ได้ T_nuc = {float(Tn):.3f} MeV ≈ 0.07 (ระดับเดียวกัน)", 0.05 < Tn < 0.08)
Yp = 0.2262 + 0.0135 * mp.log(5.5)
check(f"Y_p(η=5.5e-10) = {float(Yp):.3f}  (3.34)", abs(Yp - 0.249) < 0.001)

# ⑧ Saha (3.35)–(3.37)
me_, mH, Xe, ntot = sp.symbols('m_e m_H X_e n_tot', positive=True)
saha = sp.simplify(n0(2, me_) * n0(2, mp_) / n0(4, mH))
book337 = (me_ * mp_ / mH * T / (2 * sp.pi))**sp.Rational(3, 2) * sp.exp(-(me_ + mp_ - mH) / T)
check("g_e g_p / g_H = 2·2/4 = 1 ⇒ n_e⁰n_p⁰/n_H⁰ = (m_e m_p T/(2π m_H))^{3/2} e^{−ε₀/T}",
      sp.simplify(saha / book337) == 1)
ne = Xe * ntot
nH = (1 - Xe) * ntot
check("n_e n_p/n_H = X_e² n_tot/(1−X_e)  (ซ้ายมือของ 3.37)", sp.simplify(ne * ne / nH - Xe**2 * ntot / (1 - Xe)) == 0)
eps0 = 13.6
big = 1e9 * (me / (2 * mp.pi * eps0))**1.5
check(f"ที่ T = ε₀: 10⁹(m_e/2πT)^{{3/2}} = {float(big):.1e} ~ 10^15 (ระดับขนาด)", 1e14 < big < 1e16)
# (3.38)→(3.39): หารด้วย n_b ทั้งสองข้าง
nbS, al2, beta = sp.symbols('n_b alpha2 beta', positive=True)
r338 = nbS * ((1 - Xe) * beta - Xe**2 * nbS * al2)
check("(3.38)/n_b ⇒ dX_e/dt = (1−X_e)β − X_e² n_b α⁽²⁾  (3.39)",
      sp.simplify(r338 / nbS - ((1 - Xe) * beta - Xe**2 * nbS * al2)) == 0)

# ⑨ decoupling (3.43)–(3.47)
nb1 = 1.878e-29 / 1.67e-24          # ใช้ m_p ≈ 1.67e-24 g
coef = nb1 * 0.665e-24
check(f"n_b σ_T ที่ Ω_b h²=1: {coef:.3e} cm⁻¹ ≈ 7.477e-30  (3.44)", abs(coef / 7.477e-30 - 1) < 2e-3)
H0cm = 100e5 / 3.0856e24 / 2.998e10       # cm⁻¹ ต่อ h
check(f"7.477e-30 / (H0/c) = {7.477e-30/H0cm:.4f} ≈ 0.0692  (3.45)", abs(7.477e-30 / H0cm - 0.0692) < 5e-4)
z_ = sp.symbols('z', positive=True)
Obh2, Omh2 = sp.symbols('Obh2 Omh2', positive=True)
# H/H0 = Ωm^{1/2} a^{-3/2} [1 + a_eq/a]^{1/2}, a_eq = 4.15e-5/Ωmh²
c346 = 0.0692 * 0.02 * 0.15**-0.5 * 1000**1.5
check(f"สัมประสิทธิ์ของ (3.46) = 0.0692·0.02·0.15^(-1/2)·1000^(3/2) = {c346:.1f} ≈ 113", abs(c346 - 113) < 0.5)
check(f"a_eq/a = (1+z)·4.15e-5/Ω_mh² เทียบ (1+z)/3600 · 0.15/Ω_mh²: {0.15/3600:.3e} ≈ 4.15e-5",
      abs(0.15 / 3600 / 4.15e-5 - 1) < 0.01)
zd = 1000 * 113**(-2 / 3)
check(f"X_e = 1, ตัดวงเล็บ a_eq: 1+z = 1000·113^(−2/3) = {zd:.1f} ≈ 43  (3.47)", abs(zd - 43) < 0.5)
zdec = sp.solve(sp.Eq(113 * (Obh2 / sp.Rational(2, 100)) * (sp.Rational(15, 100) / Omh2)**sp.Rational(1, 2)
                      * (z_ / 1000)**sp.Rational(3, 2), 1), z_)[0]
check("กำลังของ (3.47): (0.02/Ω_bh²)^{2/3} (Ω_mh²/0.15)^{1/3}",
      sp.simplify(zdec / (Obh2**sp.Rational(-2, 3) * Omh2**sp.Rational(1, 3))).free_symbols == set())

# ⑩ dark matter (3.48)–(3.59)
Y, YEQ, sv, Hm = sp.symbols('Y Y_EQ sigma_v H_m', positive=True)
dYdt = T**3 * sv * (YEQ**2 - Y**2)                              # (3.50)
dYdx = (dYdt / (x * Hm / x**2)).subs(T, m / x)                  # dx/dt = Hx, H = H(m)/x²
lam_ = m**3 * sv / Hm                                           # (3.53)
check("(3.50) + dx/dt = Hx + H = H(m)/x² ⇒ dY/dx = −(λ/x²)(Y² − Y_EQ²)  (3.52)",
      sp.simplify(dYdx - (-lam_ / x**2 * (Y**2 - YEQ**2))) == 0)
xf, Yf, Yinf = sp.symbols('x_f Y_f Y_inf', positive=True)
Yx = sp.Function('Y')
sol = sp.dsolve(sp.Eq(Yx(x).diff(x), -lam * Yx(x)**2 / x**2), Yx(x), ics={Yx(xf): Yf})
Yinf_exact = sp.limit(sol.rhs, x, sp.oo)
check("อินทิเกรต (3.54) จาก x_f ถึง ∞ ⇒ 1/Y∞ − 1/Y_f = λ/x_f  (3.55)",
      sp.simplify(1 / Yinf_exact - 1 / Yf - lam / xf) == 0)
check("Y_f ≫ Y∞ ⇒ Y∞ ≈ x_f/λ  (3.56)", sp.limit(Yinf_exact, Yf, sp.oo) == xf / lam)
T0s, rc = sp.symbols('T_0 rho_cr', positive=True)
Om_X = xf * m * T0s**3 / (lam_ * 30 * rc)
check("Ω_X = x_f m T0³/(30λρcr) = H(m) x_f T0³/(30 m² ⟨σv⟩ ρcr)  (3.58)",
      sp.simplify(Om_X - Hm * xf * T0s**3 / (30 * m**2 * sv * rc)) == 0)
Om359 = Om_X.subs(Hm, sp.sqrt(4 * sp.pi**3 * Gn * gs / 45) * m**2)
check("แทน H(m) = √(4π³Gg*/45) m² ⇒ m หายไป  (3.59)", m not in sp.simplify(Om359).free_symbols)
# (3.60) แทนตัวเลขเอง
cm2 = 1 / hbarc**2                    # 1 cm² ในหน่วย eV⁻²
def OmX(sigv_cm2):
    return 10 * mp.sqrt(4 * mp.pi**3 * 100 / 45) / mpl * T0**3 / (30 * sigv_cm2 * cm2 * rcr_eV4)
v39 = OmX(1e-39)
v37 = OmX(1e-37)
print(f"   (3.59) ที่ x_f=10, g*=100: ⟨σv⟩=1e-39 cm² ⇒ Ω_X h² = {float(v39):.1f} · ⟨σv⟩=1e-37 cm² ⇒ {float(v37):.3f}")
check("⚠️ (3.60) ตามที่พิมพ์ (0.3 h⁻² ที่ 1e-39 cm²) ไม่ตรงกับ (3.59) — ห่างราว 100 เท่า", 50 < v39 / 0.3 < 150)
check("ถ้าใช้ ⟨σv⟩ = 1e-37 cm² จะได้ Ω_X h² ≈ 0.3", abs(v37 - 0.3) < 0.05)

print(f"\nผ่านทั้งหมด {ok} ข้อ")
