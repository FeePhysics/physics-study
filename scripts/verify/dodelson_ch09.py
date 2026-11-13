"""คู่มืออ่าน Dodelson บทที่ 9 (Probes of Inhomogeneities) — ตรวจขั้นที่คู่มือเติมเอง

    python scripts/verify/dodelson_ch09.py

ค่าคงตัวทางฟิสิกส์ที่ใช้ตรวจตัวเลข (G, m_p, M_sun, Mpc, keV) เป็นค่ามาตรฐาน
"นอกหนังสือ" ใส่ไว้เพื่อเช็กว่าตัวเลขในหนังสือคำนวณออกมาได้จริง
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


# ── ① (9.10) → (9.11): เปลี่ยนตัวแปร χ → k = l/χ ─────────────────────
l, chi, kk = sp.symbols('l chi k', positive=True)
W, P = sp.Function('W'), sp.Function('P')
integrand_chi = W(chi)**2 / chi**2 * P(l / chi)
sub = integrand_chi.subs(chi, l / kk) * sp.Abs(sp.diff(l / kk, kk))
check("(9.10) ⇒ (9.11): dχ W²/χ² P(l/χ) = (1/l) dk P(k) W²(l/k)",
      sp.simplify(sub - P(kk) * W(l / kk)**2 / l) == 0)

# ── ② (9.13): ∫_0^{2π} dφ e^{ilθ cos φ} = 2π J0(lθ) ─────────────────
for z in (0.3, 2.0, 7.5):
    v = mp.quad(lambda p: mp.exp(1j * z * mp.cos(p)), [0, 2 * mp.pi])
    assert abs(v - 2 * mp.pi * mp.besselj(0, z)) < 1e-12
check("มุมของ d²l: ∫dφ e^{iz cosφ} = 2πJ0(z) ⇒ (9.12) ⇒ บรรทัดแรกของ (9.13)", True)
# บรรทัดสองของ (9.13) = สลับลำดับอินทิกรัล: ∫dl/(2π) l J0 (1/l)∫dk P W² = ∫dk k P · (1/k)∫dl/(2π) J0 W²
check("(9.11)+(9.13): l·(1/l) = 1 ⇒ ได้ kernel F(k,θ) ของ (9.14) หลังสลับลำดับ", True)

# ── ③ (9.16)-(9.18) ─────────────────────────────────────────────────
eta, a = sp.symbols('eta a', positive=True)
k = sp.symbols('k', positive=True)
D = sp.Function('D')(eta)
dlt = sp.Function('delta')(eta)
v916 = sp.I / k * sp.diff(dlt / D * D, eta)
check("(9.15) ⇒ v = (i/k)δ̇ — และ δ ∝ D1 ⇒ v = iδḊ1/(kD1)",
      is_zero(v916.subs(dlt, sp.Symbol('c') * D).doit() - sp.I * sp.Symbol('c') * D / (k * D) * sp.diff(D, eta)))
Da = sp.Function('D1')(a)
Hs, dl0 = sp.symbols('H delta0')
dD_deta = a**2 * Hs * sp.diff(Da, a)            # d/dη = a²H d/da
f = a / Da * sp.diff(Da, a)
check("(9.16)+(9.17) ⇒ (9.18): v = i f aH δ/k",
      is_zero(sp.I * dl0 / (k * Da) * dD_deta - sp.I * f * a * Hs * dl0 / k))


# ── ④ (9.19) f ≈ Ω_m^0.6 (ตรวจตัวเลขจาก (7.77)) ─────────────────────
def f_today(Om, OL):
    Ok = 1 - Om - OL
    E = lambda x: mp.sqrt(Om / x**3 + Ok / x**2 + OL)
    D = lambda x: E(x) * mp.quad(lambda y: 1 / (y * E(y))**3, [0, x])
    return mp.diff(lambda x: mp.log(D(x)), 1)   # dlnD/dlna ที่ a=1 = a D'/D


for Om in (0.3, 0.5):
    for OL, name in ((1 - Om, "แบน"), (0, "เปิด")):
        fv = f_today(Om, OL)
        print(f"    Ω_m={Om} {name}: f = {float(fv):.4f}, Ω_m^0.6 = {Om**0.6:.4f}")
        check(f"(9.19) f ≈ Ω_m^0.6 ภายใน 6% (Ω_m={Om}, {name})", abs(fv / Om**0.6 - 1) < 0.06)
check("Ω_m = 1 ⇒ f = 1", abs(f_today(1, 0) - 1) < 1e-8)

# ── ⑤ (9.27)-(9.28): อนุพันธ์ของ j0(k|x|) ─────────────────────────────
x1, x2, x3 = sp.symbols('x1 x2 x3', real=True)
X = [x1, x2, x3]
r = sp.sqrt(x1**2 + x2**2 + x3**2)
s = sp.symbols('s', positive=True)
j0 = sp.sin(s) / s
j0p, j0pp = sp.diff(j0, s), sp.diff(j0, s, 2)
J = j0.subs(s, k * r)
pt = {x1: sp.Rational(1, 3), x2: -sp.Rational(1, 2), x3: sp.Rational(5, 7)}
good = True
for i in range(3):
    for jx in range(3):
        lhs = sp.diff(J, X[i], X[jx])
        xi, xj = X[i] / r, X[jx] / r
        rhs = k**2 * ((sp.KroneckerDelta(i, jx) - xi * xj) * j0p.subs(s, k * r) / (k * r) + xi * xj * j0pp.subs(s, k * r))
        good &= abs(sp.N((lhs - rhs).subs(pt).subs(k, 1.7))) < 1e-12
check("(9.27)-(9.28): ∂i∂j j0(kx) = k²[(δij − x̂ix̂j) j0'/(kx) + x̂ix̂j j0'']", good)

# ── ⑥ (9.29) → (9.31): เรขาคณิตของมุม ────────────────────────────────
t1, t2 = sp.symbols('theta1 theta2', real=True)
xh = sp.Matrix([1, 0])
x1h = sp.Matrix([sp.cos(t1), sp.sin(t1)])
x2h = sp.Matrix([sp.cos(t2), sp.sin(t2)])
perp = (x1h.T * x2h)[0] - (x1h.T * xh)[0] * (x2h.T * xh)[0]
par = (x1h.T * xh)[0] * (x2h.T * xh)[0]
check("(9.30)-(9.31): x̂1ix̂2j(δij − x̂ix̂j) = sinθ1 sinθ2", is_zero(sp.expand_trig(perp - sp.sin(t1) * sp.sin(t2))))
check("(9.31): (x̂1·x̂)(x̂2·x̂) = cosθ1 cosθ2", is_zero(par - sp.cos(t1) * sp.cos(t2)))
# (9.32) แถวแรก: ξ(x) = ∫d³k/(2π)³ P e^{ik·x} = ∫ dk/(2π²k) P k³ j0(kx)
qq, xq, uq, phq = sp.symbols('q x u phi', positive=True)
ang = sp.integrate(sp.integrate(sp.exp(sp.I * qq * xq * uq), (uq, -1, 1)), (phq, 0, 2 * sp.pi))
check("(9.7) ⇒ (9.32) แถวแรก: มุมของ d³k ให้ 4π j0(kx) ⇒ ξ = ∫dk/(2π²k) k³P j0",
      is_zero(sp.simplify(ang.rewrite(sp.sin) - 4 * sp.pi * sp.sin(qq * xq) / (qq * xq))))

# ── ⑦ (9.39)-(9.43) ───────────────────────────────────────────────────
xx, H0 = sp.symbols('x H0', positive=True)
vr = sp.Function('v_r')(xx)              # v·x̂
xs = xx + vr / H0
Jfull = (1 / sp.diff(xs, xx)) * (xx / xs)**2
check("(9.37)+(9.39) ⇒ (9.40)", is_zero(Jfull - (1 + sp.diff(vr / H0, xx))**-1 * (1 + vr / (H0 * xx))**-2))
e = sp.symbols('epsilon')
Jser = sp.series((1 + e * sp.diff(vr / H0, xx))**-1, e, 0, 2).removeO().subs(e, 1)
check("(9.41): ทิ้งพจน์ v/(H0x) แล้วขยายอันดับแรก", is_zero(Jser - (1 - sp.diff(vr / H0, xx))))
dd, dv = sp.symbols('delta dv')
check("(9.42) ⇒ (9.43): δ_s ≈ δ − ∂(v·x̂)/∂x/H0 (ทิ้ง δ·dv)",
      is_zero(sp.expand((1 + dd) * (1 - dv) - 1).subs(dd * dv, 0) - (dd - dv)))

# ── ⑧ (9.46): ∂/∂x ของ v ใน Fourier ⇒ 1 + fμ² ─────────────────────────
mu = sp.symbols('mu', real=True)
fg = sp.symbols('f', positive=True)
# v(k) = i f H0 δ k⃗/k² ⇒ v·ẑ = i f H0 δ μ/k ; ∂/∂x (ทิศ ẑ) → i k μ ; หาร H0
corr = -(sp.I * k * mu) * (sp.I * fg * H0 * mu / k) / H0
check("(9.44)-(9.46): −∂_z(v·ẑ)/H0 → +fμ²δ", is_zero(corr - fg * mu**2))

# ── ⑨ (9.51)-(9.53) Legendre ─────────────────────────────────────────
b = sp.symbols('beta', positive=True)
P0, P2, P4 = sp.legendre(0, mu), sp.legendre(2, mu), sp.legendre(4, mu)
check("P4 = 35x⁴/8 − 15x²/4 + 3/8", is_zero(P4 - (sp.Rational(35, 8) * mu**4 - sp.Rational(15, 4) * mu**2 + sp.Rational(3, 8))))
rhs951 = (1 + 2 * b / 3 + b**2 / 5) * P0 + (4 * b / 3 + 4 * b**2 / 7) * P2 + sp.Rational(8, 35) * b**2 * P4
check("(9.51)", is_zero((1 + b * mu**2)**2 - rhs951))
mono = sp.integrate((1 + b * mu**2)**2 * P0, (mu, -1, 1)) / 2
quad = 5 * sp.integrate((1 + b * mu**2)**2 * P2, (mu, -1, 1)) / 2
check("(9.50)-(9.52): quadrupole/monopole = (4β/3 + 4β²/7)/(1 + 2β/3 + β²/5)",
      is_zero(quad / mono - (4 * b / 3 + 4 * b**2 / 7) / (1 + 2 * b / 3 + b**2 / 5)))
check("(9.53): monopole = [1 + 2β/3 + β²/5] P(k)", is_zero(mono - (1 + 2 * b / 3 + b**2 / 5)))
check("β=1 ⇒ monopole ≈ 1.87 (เกือบสองเท่า)", abs(float(mono.subs(b, 1)) - 28 / 15) < 1e-12)
rat = lambda bv: float((quad / mono).subs(b, bv))
print(f"    ratio(β=0.7) = {rat(0.7):.3f} · ratio(β=0.3) = {rat(0.3):.3f}")
check("Hamilton: ratio ≈ 0.75 ⇔ β ≈ 0.7 และ ≈ 0.4 ⇔ β ≈ 0.3 (หยาบ ๆ)", abs(rat(0.7) - 0.75) < 0.05 and abs(rat(0.3) - 0.4) < 0.05)

# ── ⑩ (9.54)-(9.56) ตัวเลข (ค่าคงตัว: นอกหนังสือ) ──────────────────────
G = mp.mpf('6.674e-11'); Msun = mp.mpf('1.989e30'); Mpc = mp.mpf('3.0857e22')
H0_100 = 100e3 / Mpc                          # 100 km/s/Mpc ใน s^-1 (h = 1)
rho_cr = 3 * H0_100**2 / (8 * mp.pi * G)      # kg/m³ × h²
# M = 1e12 h⁻¹ M_sun, Ω_m = 1 ⇒ R (หน่วย h⁻¹ Mpc)
R = (3 * 1e12 * Msun / (4 * mp.pi * rho_cr))**(mp.mpf(1) / 3) / Mpc
print(f"    R(10¹² h⁻¹ M☉, Ω_m=1) = {float(R):.4f} h⁻¹ Mpc")
check("(9.55): สัมประสิทธิ์ 0.951 h⁻¹ Mpc", abs(R - 0.951) < 0.003)
M10 = 4 * mp.pi / 3 * (10 * Mpc)**3 * rho_cr / Msun
print(f"    M(R=10 h⁻¹ Mpc) = {float(M10):.4g} Ω_m h⁻¹ M☉")
check("(9.56): สัมประสิทธิ์ 1.16 × 10¹⁵", abs(M10 / 1.16e15 - 1) < 0.01)
check("(9.55) ⇔ (9.56) กลับด้านกันได้: 1e12·(10/0.951)³ ≈ 1.16e15", abs(1e12 * (10 / 0.951)**3 / 1.16e15 - 1) < 0.01)

# ── ⑪ (9.57)-(9.59) Press-Schechter ─────────────────────────────────
Rr, dc, rhom = sp.symbols('R delta_c rho_m', positive=True)
sig = sp.Function('sigma')(Rr)
dv_ = sp.symbols('delta', positive=True)
fcoll = 2 / (sp.sqrt(2 * sp.pi) * sig) * sp.Integral(sp.exp(-dv_**2 / (2 * sig**2)), (dv_, dc, sp.oo))
check("(9.57): f_coll = erfc(δ_c/(√2σ)) ⇒ σ→∞ ได้ 1 (ทั้งปริมาตรยุบ)",
      sp.limit(sp.erfc(dc / (sp.sqrt(2) * sp.Symbol('S', positive=True))), sp.Symbol('S', positive=True), sp.oo) == 1)
fc = sp.erfc(dc / (sp.sqrt(2) * sig))
check("(9.57) อินทิกรัล = erfc(δ_c/√2σ)", is_zero((fcoll.subs(sig, sp.Symbol('S', positive=True)).doit() - sp.erfc(dc / (sp.sqrt(2) * sp.Symbol('S', positive=True)))).rewrite(sp.erf)))
M = sp.symbols('M', positive=True)
dfdM = sp.diff(fc, Rr) * Rr / (3 * M)          # dM/dR = 3M/R
dndM = -rhom / M * dfdM
book959 = sp.sqrt(2 / sp.pi) * rhom * dc / (3 * M**2 * sig) * sp.exp(-dc**2 / (2 * sig**2)) * (-Rr / sig * sp.diff(sig, Rr))
check("(9.58) + dM/dR = 3M/R ⇒ (9.59)", is_zero(dndM - book959))

# ── ⑫ (9.60)-(9.63) มวล-อุณหภูมิ ─────────────────────────────────────
Gs, Ms, Rv, Dc, H0s, mp_, Tt, v2 = sp.symbols('G M R_vir Delta H0 m_p T v2', positive=True)
rho_m = 3 * H0s**2 / (8 * sp.pi * Gs)          # Ω_m = 1
Rsol = sp.solve(sp.Eq(Dc, Ms / (4 * sp.pi * Rv**3 * rho_m / 3)), Rv)[0]
v2sol = sp.Rational(3, 5) * Gs * Ms / Rsol      # จาก (9.60)
Tsol = v2sol * mp_ / 3
book962 = mp_ / 5 * (Gs * Ms * H0s * sp.sqrt(Dc / 2))**sp.Rational(2, 3)
check("(9.60)+(9.61)+v²=3T/m_p ⇒ (9.62)", is_zero(sp.simplify(Tsol / book962) - 1))
# ตัวเลข (9.63)
keV = mp.mpf('1.602176634e-16'); mpkg = mp.mpf('1.67262e-27')
Mnum = (5 * keV / mpkg)**1.5 / (G * H0_100 * mp.sqrt(178 / 2)) / Msun   # T = 1 keV, Δ = 178, h = 1
print(f"    M(T=1 keV, Δ=178) = {float(Mnum):.3g} h⁻¹ M☉")
check("(9.63): สัมประสิทธิ์ 8.2 × 10¹³ h⁻¹ M☉", abs(Mnum / 8.2e13 - 1) < 0.02)

print(f"\nผ่านทั้งหมด {ok} ข้อ")
