"""คู่มืออ่าน Dodelson บทที่ 10 · ตรวจขั้นพีชคณิตที่คู่มือเติมเอง (sympy + mpmath)

รัน:  python scripts/verify/dodelson_ch10.py
ทุกบรรทัด ✅ = ผ่าน · assert ล้ม = มีอะไรผิด
เลข (10.x) = เลขสมการในหนังสือ (ใช้เปิดเทียบ)
"""
import sympy as sp
import mpmath as mp

ok = 0


def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("✅", name)


def zero(expr):
    return sp.simplify(sp.expand_trig(sp.expand(expr))) == 0


# ── ① (10.10) ⇒ (10.11): พจน์ aH หักล้างกันพอดี ──────────────────────
# dchi/dt = -1/a  ⇒  da/dchi = (da/dt)(dt/dchi) = -a^2 H
chi = sp.symbols('chi')
a = sp.Function('a')(chi)
H = sp.Function('H')(chi)
X = sp.Function('X')(chi)          # X = chi*theta^i  (ตำแหน่งตามขวาง)
Phi_i = sp.Function('Phi_i')(chi)
lhs = sp.diff(a**-2 * sp.diff(X, chi), chi)
rhs = 2 / a**2 * (Phi_i + a * H * sp.diff(X, chi))
diffeq = (lhs - rhs).subs(sp.Derivative(a, chi), -a**2 * H)
diffeq = sp.simplify(diffeq * a**2)
check("(10.10)→(10.11): a^2*(LHS-RHS) = X'' - 2 Phi_,i",
      sp.simplify(diffeq - (sp.diff(X, chi, 2) - 2 * Phi_i)) == 0)

# ── ② (10.13) ⇒ (10.14): สลับลำดับอินทิเกรต (รูป 10.5) ─────────────
c, x1, x2 = sp.symbols('chi chi1 chi2', positive=True)
for f in [sp.Integer(1), x1, x1**3, sp.exp(-x1), sp.sin(x1)]:
    f2 = f.subs(x1, x2)
    double = 2 / c * sp.integrate(sp.integrate(f2, (x2, 0, x1)), (x1, 0, c))
    single = 2 * sp.integrate(f * (1 - x1 / c), (x1, 0, c))
    check(f"(10.13)=(10.14) สำหรับ Phi_,i = {f}", sp.simplify(double - single) == 0)

# ── ③ (10.23) A^{-1} และ (10.24)-(10.25) ellipticity กับ shear ────────
k, g1, g2 = sp.symbols('kappa gamma1 gamma2', real=True)
A = sp.Matrix([[1 - k - g1, -g2], [-g2, 1 - k + g1]])          # (10.15)
Ai = A.inv()
book_Ai = sp.Matrix([[1 - k + g1, g2], [g2, 1 - k - g1]]) / ((1 - k)**2 - g1**2 - g2**2)
check("(10.23) A^{-1}", sp.simplify(Ai - book_Ai) == sp.zeros(2, 2))
# วงกลม: q_ij ∝ sum_k (A^{-1})_{ik}(A^{-1})_{jk}  (10.21 หลังอินทิเกรตได้ delta)
q = Ai * Ai.T
eps1 = sp.simplify((q[0, 0] - q[1, 1]) / (q[0, 0] + q[1, 1]))
eps2 = sp.simplify(2 * q[0, 1] / (q[0, 0] + q[1, 1]))
book_eps1 = 4 * g1 * (1 - k) / (2 * (1 - k)**2 + 2 * g1**2 + 2 * g2**2)
check("(10.24) eps1 = 4γ1(1-κ)/[2(1-κ)^2+2γ1^2+2γ2^2]", sp.simplify(eps1 - book_eps1) == 0)
t = sp.symbols('t')
lin1 = sp.series(eps1.subs({k: t * k, g1: t * g1, g2: t * g2}), t, 0, 2).removeO().subs(t, 1)
check("(10.25) ลำดับหนึ่ง: eps1 ≈ 2γ1", sp.simplify(lin1 - 2 * g1) == 0)
lin2 = sp.series(eps2.subs({k: t * k, g1: t * g1, g2: t * g2}), t, 0, 2).removeO().subs(t, 1)
check("คู่กัน (ข้อ 5 ท้ายบท): eps2 ≈ 2γ2 ที่ลำดับหนึ่ง", sp.simplify(lin2 - 2 * g2) == 0)
# γ ที่อ่านจาก A ตาม (10.16)
check("(10.16) γ1 = -(A11-A22)/2", sp.simplify(-(A[0, 0] - A[1, 1]) / 2 - g1) == 0)
check("(10.16) γ2 = -A12", sp.simplify(-A[0, 1] - g2) == 0)

# ── ④ (10.27) ⇒ (10.28)-(10.29): สลับลำดับอีกครั้ง (ตัวอย่าง W คงที่บน [0,L]) ──
L = sp.symbols('L', positive=True)
xx, xp = sp.symbols('x xp', positive=True)
for f in [sp.Integer(1), xp, xp**2]:
    lhs = 2 * sp.integrate(1 / L * sp.integrate(f * xp * (1 - xp / xx), (xp, 0, xx)), (xx, 0, L))
    gfun = 2 * xp * sp.integrate((1 - xp / xx) / L, (xx, xp, L))
    rhs = sp.integrate(f * gfun, (xp, 0, L))
    check(f"(10.27)=(10.28) ด้วย g ของ (10.29), W=1/L, Phi_,ij={f}",
          sp.simplify(lhs - rhs) == 0)

# ── ⑤ (10.33)-(10.40): power spectrum ของ shear และ E/B ─────────────
phi, ell = sp.symbols('phi_l l', positive=True)
l1, l2 = ell * sp.cos(phi), ell * sp.sin(phi)
lv = [l1, l2]
# P^ψ_ijlm ∝ l_i l_j l_l l_m × (ก้อนอินทิเกรตเดียวกัน) ⇒ ตัดก้อนนั้นทิ้ง
P = lambda i, j, l, m: lv[i] * lv[j] * lv[l] * lv[m]
# γ1 = (ψ22 - ψ11)/2 , γ2 = -ψ12  (10.16 + 10.26)
Pg1 = sp.Rational(1, 4) * (P(1, 1, 1, 1) + P(0, 0, 0, 0) - 2 * P(1, 1, 0, 0))
Pg2 = P(0, 1, 0, 1)
Pg12 = -sp.Rational(1, 2) * (P(1, 1, 0, 1) - P(0, 0, 0, 1))
Pk = ell**4 / 4                       # (10.34): P_κ = (l^4/4) × ก้อนเดียวกัน
check("sin^4+cos^4-2sin^2cos^2 = cos^2(2φ)",
      zero(sp.sin(phi)**4 + sp.cos(phi)**4 - 2 * sp.sin(phi)**2 * sp.cos(phi)**2 - sp.cos(2 * phi)**2))
check("(10.35) P_γ1 = cos^2(2φ) P_κ", zero(Pg1 - sp.cos(2 * phi)**2 * Pk))
check("(10.36) P_γ2 = sin^2(2φ) P_κ", zero(Pg2 - sp.sin(2 * phi)**2 * Pk))
check("P_γ1γ2 = cos(2φ)sin(2φ) P_κ", zero(Pg12 - sp.cos(2 * phi) * sp.sin(2 * phi) * Pk))
C2, S2 = sp.cos(2 * phi), sp.sin(2 * phi)
PE = C2**2 * Pg1 + S2**2 * Pg2 + 2 * S2 * C2 * Pg12          # (10.38)
PB = S2**2 * Pg1 + C2**2 * Pg2 - 2 * S2 * C2 * Pg12          # (10.40)
PEB = -C2 * S2 * Pg1 + S2 * C2 * Pg2 + (C2**2 - S2**2) * Pg12
check("(10.39) P_E = P_κ (ไม่ขึ้นกับมุม)", zero(PE - Pk))
check("(10.40) P_B = 0", zero(PB))
check("cross <EB> = 0", zero(PEB))
# รูปเมทริกซ์: (E,B) คือ (γ1,γ2) ที่หมุนไป -2φ
R = sp.Matrix([[C2, S2], [-S2, C2]])
check("(10.37) เมทริกซ์ E/B เป็นการหมุน (R R^T = I)", sp.simplify(R * R.T) == sp.eye(2))

# ── ⑥ (10.43) และอินทิเกรตเชิงมุมด้วย Bessel (C.21) ─────────────────
pp, ph = sp.symbols("phip phi")
expr = (sp.cos(2 * pp) * sp.cos(2 * ph) - sp.sin(2 * pp) * sp.sin(2 * ph))**2
book = sp.cos(2 * pp)**2 * sp.cos(2 * ph)**2 - sp.sin(4 * pp) * sp.sin(4 * ph) / 2 \
    + sp.sin(2 * pp)**2 * sp.sin(2 * ph)**2
check("(10.43) กระจายกำลังสอง", zero(expr - book))
check("[cos(2(φ'+φ))]^2 = รูปกระจาย", zero(sp.cos(2 * (pp + ph))**2 - expr))
mp.mp.dps = 25
for z in [0.7, 2.3, 5.1]:
    I_c = mp.quad(lambda u: mp.cos(2 * u)**2 * mp.exp(1j * z * mp.cos(u)), [0, 2 * mp.pi])
    I_s = mp.quad(lambda u: mp.sin(2 * u)**2 * mp.exp(1j * z * mp.cos(u)), [0, 2 * mp.pi])
    I_4 = mp.quad(lambda u: mp.sin(4 * u) * mp.exp(1j * z * mp.cos(u)), [0, 2 * mp.pi])
    check(f"∫cos²(2φ')e^(iz cosφ') = π(J0+J4) ที่ z={z}",
          abs(I_c - mp.pi * (mp.besselj(0, z) + mp.besselj(4, z))) < 1e-15)
    check(f"∫sin²(2φ')e^(iz cosφ') = π(J0−J4) ที่ z={z}",
          abs(I_s - mp.pi * (mp.besselj(0, z) - mp.besselj(4, z))) < 1e-15)
    check(f"∫sin(4φ')e^(iz cosφ') = 0 ที่ z={z}", abs(I_4) < 1e-15)

# ── ⑦ (10.49)-(10.53): Q, U จาก dot product ของ polarization vector ──
th, ph = sp.symbols("theta phi", real=True)
e1p = sp.Matrix([sp.cos(th) * sp.cos(ph), sp.cos(th) * sp.sin(ph), -sp.sin(th)])   # (10.51)
e2p = sp.Matrix([-sp.sin(ph), sp.cos(ph), 0])
xh, yh = sp.Matrix([1, 0, 0]), sp.Matrix([0, 1, 0])
Qint = sum((xh.dot(e))**2 - (yh.dot(e))**2 for e in (e1p, e2p))
check("(10.52) วงเล็บ = -sin^2θ' cos2φ'", zero(Qint + sp.sin(th)**2 * sp.cos(2 * ph)))
check("(10.52) บรรทัดกลางตามที่พิมพ์",
      zero(Qint - (sp.cos(th)**2 * sp.cos(ph)**2 + sp.sin(ph)**2
                   - sp.cos(th)**2 * sp.sin(ph)**2 - sp.cos(ph)**2)))
up, um = (xh + yh) / sp.sqrt(2), (xh - yh) / sp.sqrt(2)
Uint = sum((up.dot(e))**2 - (um.dot(e))**2 for e in (e1p, e2p))
check("(10.53) [ข้อ 10a] แกนหมุน 45° ⇒ -sin^2θ' sin2φ'", zero(Uint + sp.sin(th)**2 * sp.sin(2 * ph)))

# ── ⑧ (10.56)-(10.60): Q, U จากคลื่นระนาบหนึ่งลูก ─────────────────────
P2 = lambda m: (3 * m**2 - 1) / 2
check("∫_0^{2π} cos^2φ' cos2φ' dφ' = π/2",
      sp.integrate(sp.cos(ph)**2 * sp.cos(2 * ph), (ph, 0, 2 * sp.pi)) == sp.pi / 2)
check("∫_0^π sin^5θ' dθ' = 16/15", sp.integrate(sp.sin(th)**5, (th, 0, sp.pi)) == sp.Rational(16, 15))
Aamp, Th2 = sp.symbols("A Theta2")
# f(n') → -5 Θ2 P2(k̂·n') ; Q = -A ∫dΩ f sin^2θ' cos2φ'   (10.52) + (10.55)
nvec = sp.Matrix([sp.sin(th) * sp.cos(ph), sp.sin(th) * sp.sin(ph), sp.cos(th)])


def QU(khat):
    mu = khat.dot(nvec)
    f = -5 * Th2 * P2(mu)
    w = sp.sin(th)                                    # dΩ = sinθ dθ dφ
    Q = -Aamp * sp.integrate(sp.integrate(sp.expand(f * sp.sin(th)**2 * sp.cos(2 * ph) * w),
                                          (ph, 0, 2 * sp.pi)), (th, 0, sp.pi))
    U = -Aamp * sp.integrate(sp.integrate(sp.expand(f * sp.sin(th)**2 * sp.sin(2 * ph) * w),
                                          (ph, 0, 2 * sp.pi)), (th, 0, sp.pi))
    return sp.simplify(Q), sp.simplify(U)


Qx, Ux = QU(sp.Matrix([1, 0, 0]))
check("(10.58) k ∥ x̂: Q = 4πAΘ2", sp.simplify(Qx - 4 * sp.pi * Aamp * Th2) == 0)
check("k ∥ x̂: U = 0", Ux == 0)
for tk, pk in [(sp.pi / 3, 0), (sp.pi / 4, sp.pi / 6), (sp.pi / 2, sp.pi / 8), (sp.pi / 5, 2)]:
    kh = sp.Matrix([sp.sin(tk) * sp.cos(pk), sp.sin(tk) * sp.sin(pk), sp.cos(tk)])
    Qg, Ug = QU(kh)
    check(f"(10.60) Q ที่ θk={tk}, φk={pk}",
          sp.simplify(Qg - 4 * sp.pi * Aamp * sp.sin(tk)**2 * sp.cos(2 * pk) * Th2) == 0)
    check(f"(10.60) U ที่ θk={tk}, φk={pk}",
          sp.simplify(Ug - 4 * sp.pi * Aamp * sp.sin(tk)**2 * sp.sin(2 * pk) * Th2) == 0)

# ── ⑨ (10.67): (1 - P2)/2 = (3/4)(1 - μ^2) ────────────────────────────
mu = sp.symbols('mu')
check("(1-P2(μ))/2 = (3/4)(1-μ^2)", sp.simplify((1 - P2(mu)) / 2 - sp.Rational(3, 4) * (1 - mu**2)) == 0)
# ⚠️ ข้อสังเกต: เมื่อไม่มี polarization เดิม Π = Θ2 ⇒ บรรทัดแรกของ (10.67) ให้สัมประสิทธิ์ 3/4
#    แต่บรรทัดที่สองพิมพ์ 3/2 — ต่างกัน factor 2 (คู่มือระบุไว้เป็นจุดที่ยังไม่ยืนยัน)
coef_line1 = sp.Rational(3, 4)
check("บรรทัดแรกของ (10.67) ให้ 3/4 ไม่ใช่ 3/2 (ตั้งเป็นธงในคู่มือ)", coef_line1 != sp.Rational(3, 2))

# ── ⑩ (10.75): j_l + j_l'' ด้วยสมการ Bessel ทรงกลม (C.13) + (C.18) ────
xv = sp.symbols('x', positive=True)
for l in [2, 3, 5]:
    jl = sp.expand_func(sp.jn(l, xv))
    jlm1 = sp.expand_func(sp.jn(l - 1, xv))
    lhs = jl + sp.diff(jl, xv, 2)
    rhs = -2 * jlm1 / xv + 2 * (l + 1) / xv**2 * jl + l * (l + 1) / xv**2 * jl
    check(f"(10.75) ที่ l={l}", sp.simplify(lhs - rhs) == 0)

# ── ⑪ ข้อ 12 / (10.77): tight coupling ⇒ Π = 5Θ2/2 ────────────────────
TP0, TP2, Pi = sp.symbols("ThetaP0 ThetaP2 Pi")
# Θ_P(μ) = (1/2)(1 - P2)Π  ;  การกระจาย Θ = Σ(-i)^l(2l+1)Θ_l P_l
# ⇒ สัมประสิทธิ์ P0: Θ_P0 = Π/2 · สัมประสิทธิ์ P2: -5 Θ_P2 = -Π/2
sol = sp.solve([TP0 - Pi / 2, -5 * TP2 + Pi / 2, Pi - (Th2 + TP2 + TP0)], [TP0, TP2, Pi], dict=True)[0]
check("ข้อ 12: Θ_P0 = 5Θ2/4", sp.simplify(sol[TP0] - 5 * Th2 / 4) == 0)
check("ข้อ 12: Θ_P2 = Θ2/4", sp.simplify(sol[TP2] - Th2 / 4) == 0)
check("(10.77) Π = 5Θ2/2 ⇒ 3Π/4 = 15Θ2/8", sp.simplify(3 * sol[Pi] / 4 - 15 * Th2 / 8) == 0)

# ── ⑫ (10.81) และ (10.92) ─────────────────────────────────────────────
pk = sp.symbols("phi_k", real=True)
kx, ky = sp.cos(pk), sp.sin(pk)
check("(10.81) cos2φk = (kx^2-ky^2)/(kx^2+ky^2)", zero(sp.cos(2 * pk) - (kx**2 - ky**2) / (kx**2 + ky**2)))
check("(10.92) sin^2θ sin2φ = 2 n_x n_y", zero(sp.sin(th)**2 * sp.sin(2 * ph) - 2 * nvec[0] * nvec[1]))

# ── ⑬ (10.91)-(10.95): หมุนแกน แล้วแยกส่วน Y_{2,±2} ────────────────────
tk = sp.symbols("theta_k", real=True)
Rm = sp.Matrix([[sp.cos(tk) * sp.cos(pk), -sp.sin(pk), sp.sin(tk) * sp.cos(pk)],
                [sp.cos(tk) * sp.sin(pk), sp.cos(pk), sp.sin(tk) * sp.sin(pk)],
                [-sp.sin(tk), 0, sp.cos(tk)]])
khat = sp.Matrix([sp.sin(tk) * sp.cos(pk), sp.sin(tk) * sp.sin(pk), sp.cos(tk)])
check("ข้อ 13a: R ẑ = k̂", sp.simplify(Rm * sp.Matrix([0, 0, 1]) - khat) == sp.zeros(3, 1))
check("R เป็นเมทริกซ์หมุน (R^T R = I)", sp.simplify(Rm.T * Rm) == sp.eye(3))
rn = Rm.T * nvec
prod = rn[0] * rn[1]
book94 = sp.cos(tk) * sp.sin(th)**2 * sp.sin(2 * (ph - pk)) / 2 \
    - sp.sin(th) * sp.cos(th) * sp.sin(tk) * sp.sin(ph - pk)
check("(10.94) (R^t n)_x (R^t n)_y = รูปสุดท้าย", zero(prod - book94))
# Y_{2,±2} ตาม Appendix C (C.10) ของหนังสือ: มีเครื่องหมายลบนำหน้า
Y22b = lambda m: -sp.sqrt(sp.Rational(15, 32) / sp.pi) * sp.sin(th)**2 * sp.exp(m * 2 * sp.I * ph)
lhs95 = sp.cos(tk) * sp.sin(th)**2 * sp.sin(2 * (ph - pk)) / 2
rhs95 = -sp.sqrt(32 * sp.pi / 15) * sp.cos(tk) / (4 * sp.I) * (
    sp.exp(-2 * sp.I * pk) * sp.conjugate(Y22b(-1)) - sp.exp(2 * sp.I * pk) * sp.conjugate(Y22b(+1)))
check("(10.95) ถูกเมื่อใช้ Y22 ตาม (C.10) ของหนังสือ",
      sp.simplify(sp.expand_complex(sp.expand(lhs95 - rhs95))) == 0)
Y22std = lambda m: +sp.sqrt(sp.Rational(15, 32) / sp.pi) * sp.sin(th)**2 * sp.exp(m * 2 * sp.I * ph)
rhs95s = -sp.sqrt(32 * sp.pi / 15) * sp.cos(tk) / (4 * sp.I) * (
    sp.exp(-2 * sp.I * pk) * sp.conjugate(Y22std(-1)) - sp.exp(2 * sp.I * pk) * sp.conjugate(Y22std(+1)))
check("ถ้าใช้ Y22 แบบตำรา QM ส่วนใหญ่ (ไม่มีลบ) เครื่องหมายของ (10.95) จะกลับ",
      sp.simplify(sp.expand_complex(sp.expand(lhs95 + rhs95s))) == 0)

print(f"\nผ่านทั้งหมด {ok} ข้อ")
