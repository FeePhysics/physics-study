"""คู่มืออ่าน Dodelson บทที่ 7 (Inhomogeneities) — ตรวจขั้นพีชคณิตที่คู่มือเติมเอง

    python scripts/verify/dodelson_ch07.py

ทุกข้อเริ่มจากสมการที่พิมพ์ในหนังสือ (เปิดดูภาพหน้าจริงแล้ว) แล้วให้ sympy
ทำขั้นที่หนังสือข้าม ไม่มีการพิมพ์ผลลัพธ์จากความจำ
"""
import sympy as sp
import mpmath as mp

ok = 0


def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("✅", name)


def is_zero(expr):
    return sp.simplify(sp.together(sp.expand(expr))) == 0


y, eta, k, a = sp.symbols('y eta k a', positive=True)
Phi = sp.Function('Phi')
P = Phi(y)
dP, ddP = P.diff(y), P.diff(y, 2)

# ── ① (7.20) → (7.23): เปลี่ยนตัวแปรเวลาเป็น y = a/a_eq ─────────────────────
# (7.20): 3(ȧ/a)(Φ̇ + (ȧ/a)Φ) = 4πG a² ρ_dm δ [1 + 4/(3y)]
# ใช้ d/dη = aH y d/dy  (7.22) · ȧ/a = aH · 8πGρ_dm/3 = H² y/(y+1)
aH, H, delta = sp.symbols('aH H delta', positive=True)
lhs = 3 * aH * (aH * y * dP + aH * P)
rhs = sp.Rational(3, 2) * aH**2 * y / (y + 1) * delta * (1 + sp.Rational(4, 3) / y)  # 4πGρ_dm a² = (3/2)(aH)² y/(y+1)
eq723 = sp.Eq(y * dP + P, (3 * y + 4) / (6 * (y + 1)) * delta)
check("(7.20)+(7.22) ⇒ (7.23): yΦ'+Φ = (3y+4)δ/[6(y+1)]",
      is_zero((lhs - rhs) / (3 * aH**2) - (eq723.lhs - eq723.rhs)))

# ── ② (7.24) → (7.25): หาอนุพันธ์แล้วจัดรูป ─────────────────────────────
delta_y = 6 * (y + 1) / (3 * y + 4) * (y * dP + P)          # δ จาก (7.23)
ode = sp.diff(delta_y, y) + 3 * dP                          # δ' = -3Φ'  (7.18)
book725 = ddP + (21 * y**2 + 54 * y + 32) / (2 * y * (y + 1) * (3 * y + 4)) * dP \
    + P / (y * (y + 1) * (3 * y + 4))
ratio = sp.simplify(ode / book725)
check("(7.24) ⇒ (7.25) (ต่างกันแค่ตัวคูณที่ไม่มี Φ)", not ratio.has(Phi) and ratio != 0)
print("    ตัวคูณ =", sp.factor(ratio))

# ── ③ (7.26) → (7.27): u = y³Φ/√(1+y) ไม่มีพจน์ u เปล่า ────────────────────
u = sp.Function('u')(y)
Phi_u = u * sp.sqrt(1 + y) / y**3
ode_u = book725.subs(P, Phi_u).doit()
ode_u = ode_u.subs({sp.Derivative(Phi(y), (y, 2)): sp.diff(Phi_u, y, 2),
                    sp.Derivative(Phi(y), y): sp.diff(Phi_u, y)})
ode_u = sp.simplify(ode_u.subs(Phi(y), Phi_u))
book727 = u.diff(y, 2) + u.diff(y) * (-2 / y + sp.Rational(3, 2) / (1 + y) - 3 / (3 * y + 4))
r2 = sp.simplify(ode_u / book727)
check("(7.25) ⇒ (7.27) (ไม่มีพจน์ u)", not r2.has(u) and not r2.has(sp.Derivative))

# ── ④ (7.28)-(7.30): u' = A y²(3y+4)/(1+y)^{3/2} ─────────────────────────
up = y**2 * (3 * y + 4) / (1 + y)**sp.Rational(3, 2)
check("(7.30) แก้ (7.27)", is_zero(sp.diff(up, y) + up * (-2 / y + sp.Rational(3, 2) / (1 + y) - 3 / (3 * y + 4))))

# ── ⑤ (7.31) → (7.32): อินทิเกรตด้วย x = √(1+y) (แบบฝึกหัด 4c) ─────────────
yp = sp.symbols('yp', positive=True)
Phi0 = sp.symbols('Phi0')
A = 3 * Phi0 / 4
integral = sp.integrate(yp**2 * (3 * yp + 4) / (1 + yp)**sp.Rational(3, 2), (yp, 0, y))
Phi_731 = A * sp.sqrt(1 + y) / y**3 * integral
Phi_732 = Phi0 / 10 / y**3 * (16 * sp.sqrt(1 + y) + 9 * y**3 + 2 * y**2 - 8 * y - 16)
check("(7.31) กับ A=3Φ(0)/4 ⇒ (7.32)", is_zero(Phi_731 - Phi_732))
check("(7.32) แก้ (7.25)", is_zero(book725.subs(P, Phi_732).doit().subs(
    {sp.Derivative(Phi(y), (y, 2)): sp.diff(Phi_732, y, 2), sp.Derivative(Phi(y), y): sp.diff(Phi_732, y)})
    .subs(Phi(y), Phi_732)))
check("(7.32): y→0 ได้ Φ(0)", sp.limit(Phi_732, y, 0) == Phi0)
check("(7.32): y→∞ ได้ (9/10)Φ(0)", sp.limit(Phi_732, y, sp.oo) == sp.Rational(9, 10) * Phi0)
check("small y: integrand ≈ 4y'² ⇒ Φ → 4A/3", sp.limit(integral / y**3, y, 0) * sp.Rational(1, 1) == sp.Rational(4, 3))

# ── ⑥ §7.2.2: ในยุค matter Φ = ค่าคงที่ เป็นคำตอบของ (7.33)-(7.35) ─────────
# เส้นทางของคู่มือเอง (นอกหนังสือ): matter era a ∝ η² ⇒ aH = 2/η แล้วแก้ตรง ๆ
Pc = sp.symbols('Phi_c')
v = sp.I * k * Pc * eta / 3                     # particular solution ของ v̇ + (2/η)v = ikΦ
check("v = ikΦη/3 แก้ (7.34) เมื่อ aH = 2/η", is_zero(sp.diff(v, eta) + 2 / eta * v - sp.I * k * Pc))
aHm = 2 / eta
d = sp.solve(sp.Eq(k**2 * Pc, sp.Rational(3, 2) * aHm**2 * (sp.Symbol('dd') + 3 * aHm * sp.I * v / k)),
             sp.Symbol('dd'))[0]                 # δ จาก (7.35)
check("δ จาก (7.35) = 2Φ + k²η²Φ/6", is_zero(d - (2 * Pc + k**2 * eta**2 * Pc / 6)))
check("δ, v นี้แก้ (7.33): δ̇ + ikv = 0", is_zero(sp.diff(d, eta) + sp.I * k * v))
check("matter era: d(aH)/dη = -(aH)²/2 และ d(aH)⁻¹/dη = 1/2",
      is_zero(sp.diff(aHm, eta) + aHm**2 / 2) and is_zero(sp.diff(1 / aHm, eta) - sp.Rational(1, 2)))

# ── ⑦ §7.3.1: (7.42),(7.43) ⇒ (7.44) ⇒ (7.45) ⇒ (7.46) ─────────────────
Ph = sp.Function('Phi')(eta)
Th1 = -k * eta**2 / 6 * (Ph.diff(eta) + Ph / eta)       # Θ_r,1 จาก (7.43)
e742 = Th1.diff(eta) + Th1 / eta + k / 3 * Ph * (1 - k**2 * eta**2 / 6)
book744 = Ph.diff(eta, 2) + 4 / eta * Ph.diff(eta) + k**2 / 3 * Ph
r3 = sp.simplify(e742 / book744)
check("(7.42)+(7.43) ⇒ (7.44)", not r3.has(Ph) and not r3.has(sp.Derivative))
print("    ตัวคูณ =", r3)
uu = sp.Function('u')(eta)
e = book744.subs(Ph, uu / eta).doit()
book745 = uu.diff(eta, 2) + 2 / eta * uu.diff(eta) + (k**2 / 3 - 2 / eta**2) * uu
check("u = Φη ⇒ (7.45)", is_zero(sp.simplify(e * eta) - book745))
Pp = sp.symbols('Phi_p')
xx = k * eta / sp.sqrt(3)
Phi746 = 3 * Pp * (sp.sin(xx) - xx * sp.cos(xx)) / xx**3
check("(7.46) แก้ (7.44)", is_zero(sp.simplify(book744.subs(Ph, Phi746).doit())))
check("(7.46): η→0 ได้ Φ_p (วงเล็บ → 1/3)", sp.limit(Phi746, eta, 0) == Pp)

# ── ⑧ (7.47)-(7.49): สมการ δ ในยุครังสี ─────────────────────────────────
dl, vv = sp.Function('delta')(eta), sp.Function('v')(eta)
# ยุครังสี ȧ/a = 1/η · (7.14): v̇ = -v/η + ikΦ · (7.13): v = (−3Φ̇ − δ̇)/(ik)
v_of = (-3 * Ph.diff(eta) - dl.diff(eta)) / (sp.I * k)
lhs747 = dl.diff(eta, 2) + sp.I * k * (-v_of / eta + sp.I * k * Ph) + 3 * Ph.diff(eta, 2)
S = -3 * Ph.diff(eta, 2) + k**2 * Ph - 3 / eta * Ph.diff(eta)
check("(7.47)+(7.13) ⇒ (7.48),(7.49)",
      is_zero(sp.expand(lhs747) - sp.expand(dl.diff(eta, 2) + dl.diff(eta) / eta - S)))
for sol in (sp.Integer(1), sp.log(eta)):
    check(f"homogeneous ของ (7.48): δ = {sol}", is_zero(sp.diff(sol, eta, 2) + sp.diff(sol, eta) / eta))
# Green's function ใน (7.50): ∫ S η' (ln kη' − ln kη) ให้ particular solution
f = sp.Function('S')
etap = sp.symbols('etap', positive=True)
part = -sp.Integral(f(etap) * etap * (sp.log(k * etap) - sp.log(k * eta)), (etap, 0, eta))
res = sp.diff(part, eta, 2) + sp.diff(part, eta) / eta
check("(7.50): อินทิกรัลเป็น particular solution ของ δ̈ + δ̇/η = S", sp.simplify(res.doit() - f(eta)) == 0)

# ── ⑨ Meszaros (7.58)-(7.62) ─────────────────────────────────────────
mes = lambda D: sp.diff(D, y, 2) + (2 + 3 * y) / (2 * y * (y + 1)) * sp.diff(D, y) - 3 / (2 * y * (y + 1)) * D
D1 = y + sp.Rational(2, 3)
D2 = D1 * sp.log((sp.sqrt(1 + y) + 1) / (sp.sqrt(1 + y) - 1)) - 2 * sp.sqrt(1 + y)
check("(7.59) D1 = y + 2/3 แก้ Meszaros (7.58)", is_zero(mes(D1)))
check("(7.62) D2 แก้ Meszaros (7.58)", is_zero(mes(D2)))
check("D1'/D1 = 3/(2+3y)", is_zero(sp.diff(D1, y) / D1 - 3 / (2 + 3 * y)))
W = sp.simplify(D1 * sp.diff(D2, y) - sp.diff(D1, y) * D2)
check("Wronskian D1D2' − D1'D2 = −(4/9) y⁻¹ (y+1)^{-1/2}",
      is_zero(W + sp.Rational(4, 9) / y / sp.sqrt(y + 1)))
check("y≫1: D2 ≈ (8/45) y^{-3/2}", sp.limit(D2 * y**sp.Rational(3, 2), y, sp.oo) == sp.Rational(8, 45))
# small y
mp.mp.dps = 30
for yv in (1e-6, 1e-8):
    D2n = sp.N(D2.subs(y, yv), 25)
    approx = sp.Rational(2, 3) * sp.log(4 / sp.Float(yv, 30)) - 2
    check(f"y→0: D2 → (2/3)ln(4/y) − 2 (ที่ y={yv})", abs(D2n - approx) < 1e-4)
    D2p = sp.N(sp.diff(D2, y).subs(y, yv), 25)
    check(f"y→0: D2' → −2/(3y) (ที่ y={yv})", abs(D2p * yv + sp.Rational(2, 3)) < 1e-4)
# u = δ/D1 → (7.60)
uu = sp.Function('u')(y)
e760 = sp.expand(mes(D1 * uu) * sp.Rational(3, 2))
fixed760 = (1 + sp.Rational(3, 2) * y) * uu.diff(y, 2) + uu.diff(y) / (y * (y + 1)) * (sp.Rational(21, 4) * y**2 + 6 * y + 1)
printed760 = (1 + sp.Rational(3, 2) * y) * uu.diff(y, 2) + uu.diff(y) / (y * (y + 1)) * (sp.Rational(21, 4) * y**2 + 3 * y + 1)
check("u = δ/D1 ⇒ (7.60) ไม่มีพจน์ u · สัมประสิทธิ์ที่ได้คือ (21/4)y² + 6y + 1", is_zero(e760 - fixed760))
check("⚠️ ที่พิมพ์ในหนังสือ '(21/4)y² + 3y + 1' ไม่ตรง (น่าจะพิมพ์ผิด)", not is_zero(e760 - printed760))
book760 = fixed760
up761 = (y + sp.Rational(2, 3))**-2 / y / sp.sqrt(y + 1)
check("(7.61) u' แก้ (7.60) ฉบับแก้ 6y — ยืนยันว่า 3y เป็นที่พิมพ์ผิด", is_zero(book760.subs(uu, sp.Function('U')(y)).subs(
    sp.Derivative(sp.Function('U')(y), (y, 2)), sp.diff(up761, y)).subs(sp.Derivative(sp.Function('U')(y), y), up761)))

# ── ⑩ (7.65) → (7.66) → (7.67) ───────────────────────────────────────
Am, Bm, ym, yH, Pp = sp.symbols('A B y_m y_H Phi_p', positive=True)
D2s = sp.Rational(2, 3) * sp.log(4 / ym) - 2
D2ps = -sp.Rational(2, 3) / ym
den = -sp.Rational(4, 9) / ym
C1 = (D2ps * Am * sp.log(Bm * ym / yH) - D2s * Am / ym) / den * Pp
C1_book = -sp.Rational(9, 4) * Am * Pp * (-sp.Rational(2, 3) * sp.log(Bm * ym / yH) - sp.Rational(2, 3) * sp.log(4 / ym) + 2)
check("(7.65) ⇒ (7.66) ในลิมิต y_m ≪ 1", is_zero(C1 - C1_book))
C1_simpl = sp.simplify(sp.expand_log(C1_book, force=True))
check("(7.66) ไม่ขึ้นกับ y_m", sp.diff(C1_simpl, ym) == 0)
# C1 D1(y) ที่ y ≫ 1: D1 ≈ y = a/a_eq ⇒ (7.67)
aeq, aHs, D1a = sp.symbols('a_eq a_H D1a', positive=True)
delta767 = sp.Rational(3, 2) * Am * Pp * sp.log(4 * Bm * sp.exp(-3) * aeq / aHs) * D1a
check("(7.66) × D1 ⇒ (7.67) [D1 = a/a_eq, y_H = a_H/a_eq]",
      is_zero(sp.expand_log(C1_simpl.subs(yH, aHs / aeq) * D1a / aeq - delta767 / aeq, force=True)))

# (7.67) = (7.8) ⇒ (7.68) — ใช้ได้เมื่อ "D1" ใน (7.67) คือ D1(y) ≈ a/a_eq ของ Meszaros (7.59)
Om_, H0_, k_, T_, a__ = sp.symbols('Omega_m H0 k T a', positive=True)
L = sp.Symbol('L')  # ln[4Be^{-3} a_eq/a_H]
Tsol = sp.solve(sp.Eq(sp.Rational(3, 5) * k_**2 / (Om_ * H0_**2) * Pp * T_ * a__,
                      sp.Rational(3, 2) * Am * Pp * L * a__ / aeq), T_)[0]
check("(7.8) = (7.67) [D1(y)=a/a_eq] ⇒ (7.68): T = 5AΩ_mH0² ln[..]/(2k²a_eq)",
      is_zero(Tsol - 5 * Am * Om_ * H0_**2 * L / (2 * k_**2 * aeq)))

# ── ⑪ (7.68) → (7.69): ตัวเลข A = 9.6, B = 0.44 ────────────────────────
# k_eq = a_eq H(a_eq) = √(2Ω_m) H0 a_eq^{-1/2} ⇒ 5AΩ_m H0²/(2k² a_eq) = 5A k_eq²/(4k²)
Om, H0, kk, keq = sp.symbols('Omega_m H0 k k_eq', positive=True)
pref = 5 * 9.6 * Om * H0**2 / (2 * kk**2 * aeq)
pref_keq = pref.subs(aeq, 2 * Om * H0**2 / keq**2)
check("prefactor (7.68) = 12 k_eq²/k² (ถึงสองหลัก)", abs(sp.N(pref_keq * kk**2 / keq**2) - 12) < 0.01)
arg = 4 * 0.44 * mp.e**-3 * mp.sqrt(2)
check("4Be⁻³√2 ≈ 1/8 ⇒ ln[k/(8k_eq)] (ถึงสองหลัก)", abs(1 / arg - 8) < 0.1)
print("    1/(4Be⁻³√2) =", float(1 / arg))

# ── ⑫ (7.5)-(7.9): จาก Poisson ถึง P(k) ───────────────────────────────
T, Dl, Dl1, Ppr, dH, n, a_ = sp.symbols('T D_1 D_1one Phi_p delta_H n a', positive=True)
Phi75 = sp.Rational(9, 10) * Ppr * T * Dl / a_
delta77 = kk**2 * Phi75 * a_ / (sp.Rational(3, 2) * Om * H0**2)
check("(7.5)+(7.7) ⇒ (7.8): δ = (3/5)(k²/Ω_mH0²)Φ_p T D1",
      is_zero(delta77 - sp.Rational(3, 5) * kk**2 / (Om * H0**2) * Ppr * T * Dl))
PPhi = 50 * sp.pi**2 / (9 * kk**3) * (kk / H0)**(n - 1) * dH**2 * (Om / Dl1)**2
Pk = (sp.Rational(3, 5) * kk**2 / (Om * H0**2) * T * Dl)**2 * PPhi
check("(7.8)+(6.100) ⇒ (7.9)",
      is_zero(sp.simplify(Pk - 2 * sp.pi**2 * dH**2 * kk**n / H0**(n + 3) * T**2 * (Dl / Dl1)**2)))
check("(7.6) ⇒ (7.7): ρ_m = Ω_mρ_cr/a³, 4πGρ_cr = (3/2)H0²",
      is_zero(sp.solve(sp.Eq(Phi75, sp.Rational(3, 2) * H0**2 * Om / a_**3 * a_**2 * sp.Symbol('dd') / kk**2),
                       sp.Symbol('dd'))[0] - delta77))

# ── ⑬ §7.5 growth function ───────────────────────────────────────────
OL = sp.symbols('Omega_L', positive=True)
Hm = H0 * sp.sqrt(Om / a**3)
grow = lambda D, Hf: sp.diff(D, a, 2) + (sp.diff(sp.log(Hf), a) + 3 / a) * sp.diff(D, a) \
    - 3 * Om * H0**2 / (2 * a**5 * Hf**2) * D
check("matter ล้วน: δ ∝ H = a^{-3/2} แก้ (7.73)", is_zero(grow(Hm, Hm)))
terms = [sp.simplify(t.subs(Om, 1) / a**sp.Rational(-7, 2) / H0) for t in
         (sp.diff(Hm, a, 2), sp.diff(sp.log(Hm), a) * sp.diff(Hm, a) + 3 / a * sp.diff(Hm, a),
          -3 * H0**2 / (2 * a**5 * Hm**2) * Hm)]
check("สามพจน์มีสัมประสิทธิ์ 15/4, −9/4, −3/2 (หนังสือเขียน 3/2 = ขนาด)",
      terms == [sp.Rational(15, 4), sp.Rational(-9, 4), sp.Rational(-3, 2)])
HL = H0 * sp.sqrt(Om / a**3 + 1 - Om)
check("แบบฝึกหัด 7: δ ∝ H แก้ (7.73) เมื่อแบน + Λ", is_zero(grow(HL, HL)))
HO = H0 * sp.sqrt(Om / a**3 + (1 - Om) / a**2)
check("แบบฝึกหัด 7: δ ∝ H แก้ (7.73) เมื่อเปิด (Ω_Λ=0)", is_zero(grow(HO, HO)))
# (7.74)-(7.76): D = H·I(a), I' = (aH)^{-3} แก้ (7.73) สำหรับ H ทั่วไปที่ H² = H0²(Ω_m a⁻³ + Ω_Λ)
I = sp.Function('I')(a)
for Hf, name in ((HL, "แบน+Λ"), (HO, "เปิด")):
    D = Hf * I
    ex = grow(D, Hf).subs(sp.Derivative(I, (a, 2)), sp.diff((a * Hf)**-3, a)).subs(sp.Derivative(I, a), (a * Hf)**-3)
    check(f"(7.76) D1 = H∫da/(aH)³ แก้ (7.73) [{name}]", is_zero(ex))
# (7.77) normalisation: ยุคต้น D1 → a
for Omv, OLv in ((sp.Rational(3, 10), sp.Rational(7, 10)), (1, 0)):
    Hn = lambda av: mp.sqrt(Omv / av**3 + OLv)
    D1n = lambda av: 5 * Omv / 2 * Hn(av) * mp.quad(lambda x: 1 / (x * Hn(x))**3, [0, av])
    check(f"(7.77) D1(a) → a ที่ a = 10⁻³ (Ω_m={Omv})", abs(D1n(mp.mpf('1e-3')) / mp.mpf('1e-3') - 1) < 1e-3)
check("(7.77) Ω_m = 1 ⇒ D1(a) = a ทุก a", abs(5 / 2 * mp.quad(lambda x: x**1.5, [0, 0.7]) * 0.7**-1.5 - 0.7) < 1e-12)
D1L = lambda av: 5 * 0.3 / 2 * mp.sqrt(0.3 / av**3 + 0.7) * mp.quad(lambda x: 1 / (x * mp.sqrt(0.3 / x**3 + 0.7))**3, [0, av])
print("    D1(a=1) ใน Ω_m=0.3 แบน =", float(D1L(1)), " (< 1 = การโตถูกกดตอนปลาย ตามรูป 7.12)")
check("Λ กดการโต: D1(1) < 1 เมื่อ Ω_m = 0.3", D1L(1) < 1)

# ── ⑭ (7.10): Δ² = k³P/(2π²) จาก d³k/(2π)³ ─────────────────────────────
check("(7.10): 4πk²dk/(2π)³ = (dk/k)·k³/(2π²)", is_zero(4 * sp.pi * kk**2 / (2 * sp.pi)**3 - kk**3 / (2 * sp.pi**2) / kk))

print(f"\nผ่านทั้งหมด {ok} ข้อ")
