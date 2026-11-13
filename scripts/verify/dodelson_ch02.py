"""คู่มืออ่าน Dodelson บทที่ 2 — ตรวจขั้นที่คู่มือเติมเอง (reference/dodelson-ch02.html)

รัน: python scripts/verify/dodelson_ch02.py
Christoffel และ Ricci คำนวณจากสูตรของหนังสือ (2.19) และ (2.31) ตรง ๆ ด้วย sympy
— ไม่ได้พิมพ์ผลจากความจำ แล้วเทียบกับ (2.22) (2.23) (2.34) (2.35) (2.37) (2.39)
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


def christoffel(g, X):
    """(2.19): Γ^μ_{αβ} = ½ g^{μν} (∂_β g_{αν} + ∂_α g_{βν} − ∂_ν g_{αβ})"""
    n = len(X)
    gi = g.inv()
    return [[[sp.simplify(sum(gi[m, v] * (sp.diff(g[al, v], X[be]) + sp.diff(g[be, v], X[al])
                                          - sp.diff(g[al, be], X[v])) for v in range(n)) / 2)
              for be in range(n)] for al in range(n)] for m in range(n)]


def ricci(G, X):
    """(2.31): R_{μν} = Γ^α_{μν,α} − Γ^α_{μα,ν} + Γ^α_{βα}Γ^β_{μν} − Γ^α_{βν}Γ^β_{μα}"""
    n = len(X)
    R = sp.zeros(n)
    for m in range(n):
        for v in range(n):
            R[m, v] = sp.simplify(
                sum(sp.diff(G[al][m][v], X[al]) for al in range(n))
                - sum(sp.diff(G[al][m][al], X[v]) for al in range(n))
                + sum(G[al][be][al] * G[be][m][v] for al in range(n) for be in range(n))
                - sum(G[al][be][v] * G[be][m][al] for al in range(n) for be in range(n)))
    return R


# ───────── ① อุ่นเครื่อง: ระนาบในพิกัดเชิงขั้ว (โจทย์ข้อ 2 · หนังสือให้คำตอบ (2.92) ไว้แล้ว) ─────────
r, th, tt = sp.symbols('r theta t', positive=True)
gpol = sp.diag(1, r**2)
Gp = christoffel(gpol, [r, th])
check("ขั้ว: Γ^r_θθ = −r", sp.simplify(Gp[0][1][1] + r) == 0)
check("ขั้ว: Γ^θ_rθ = Γ^θ_θr = 1/r", sp.simplify(Gp[1][0][1] - 1 / r) == 0 and sp.simplify(Gp[1][1][0] - 1 / r) == 0)
nz = [(m, a, b) for m in range(2) for a in range(2) for b in range(2) if Gp[m][a][b] != 0]
check("ขั้ว: ตัวอื่นเป็นศูนย์หมด (เหลือ 3 ตัว)", len(nz) == 3)
# ได้ Γ ตัวเดียวกันจากวิธีที่สองคือวงเล็บใน (2.17): (∂x/∂x')^{-1} ∂²x/∂x'∂x'
x1, x2 = sp.symbols('x1 x2')
cart = sp.Matrix([x1 * sp.cos(x2), x1 * sp.sin(x2)])
J = cart.jacobian([x1, x2])
Jinv = sp.simplify(J.inv())
G217 = [[[sp.simplify(sum(Jinv[l, i] * sp.diff(cart[i], [x1, x2][j], [x1, x2][k]) for i in range(2)))
          for k in range(2)] for j in range(2)] for l in range(2)]
same = all(sp.simplify(G217[l][j][k] - Gp[l][j][k].subs({r: x1, th: x2})) == 0
           for l in range(2) for j in range(2) for k in range(2))
check("ขั้ว: วงเล็บใน (2.17) ให้ Γ ชุดเดียวกับ (2.19)", same)
# สมการเคลื่อนที่: r'' − r θ'^2 = 0 , θ'' + (2/r) r' θ' = 0 ใช้ได้กับเส้นตรงจริง
s = sp.symbols('s')
b0, v0 = sp.symbols('b v', positive=True)
rs = sp.sqrt(b0**2 + (v0 * s)**2)            # เส้นตรง y = b เดินด้วยความเร็ว v
ths = sp.atan2(b0, v0 * s)
eq_r = sp.diff(rs, s, 2) + Gp[0][1][1].subs(r, rs) * sp.diff(ths, s)**2
eq_t = sp.diff(ths, s, 2) + 2 * Gp[1][0][1].subs(r, rs) * sp.diff(rs, s) * sp.diff(ths, s)
check("เส้นตรง y=b เป็นคำตอบของ geodesic ในพิกัดขั้ว",
      sp.simplify(eq_r) == 0 and sp.simplify(eq_t) == 0)

# ───────── ② FRW แบน (2.4) ─────────
t, x, y, z = sp.symbols('t x y z')
X = [t, x, y, z]
a = sp.Function('a')(t)
ad, add = sp.diff(a, t), sp.diff(a, t, 2)
g = sp.diag(-1, a**2, a**2, a**2)
check("g^{μν}: ส่วน space เป็น 1/a^2", g.inv() == sp.diag(-1, a**-2, a**-2, a**-2))
G = christoffel(g, X)
check("Γ^0_00 = 0 · Γ^0_0i = 0  (2.22)", G[0][0][0] == 0 and all(G[0][0][i] == 0 for i in (1, 2, 3)))
check("Γ^0_ij = δ_ij a ȧ  (2.22)",
      all(sp.simplify(G[0][i][j] - (a * ad if i == j else 0)) == 0 for i in (1, 2, 3) for j in (1, 2, 3)))
check("Γ^i_0j = Γ^i_j0 = δ_ij ȧ/a  (2.23)",
      all(sp.simplify(G[i][0][j] - (ad / a if i == j else 0)) == 0 and sp.simplify(G[i][j][0] - G[i][0][j]) == 0
          for i in (1, 2, 3) for j in (1, 2, 3)))
check("Γ^i_00 = Γ^i_jk = 0  (2.23 'ตัวอื่นเป็นศูนย์')",
      all(G[i][0][0] == 0 for i in (1, 2, 3)) and
      all(G[i][j][k] == 0 for i in (1, 2, 3) for j in (1, 2, 3) for k in (1, 2, 3)))

R = ricci(G, X)
check("R_00 = −3 ä/a  (2.34)", sp.simplify(R[0, 0] + 3 * add / a) == 0)
check("R_ij = δ_ij (2ȧ² + aä)  (2.35)",
      all(sp.simplify(R[i, j] - ((2 * ad**2 + a * add) if i == j else 0)) == 0 for i in (1, 2, 3) for j in (1, 2, 3)))
check("R_0i = 0 (โจทย์ 5b)", all(sp.simplify(R[0, i]) == 0 for i in (1, 2, 3)))
Rs = sp.simplify(sum(g.inv()[m, v] * R[m, v] for m in range(4) for v in range(4)))
check("Ricci scalar = 6(ä/a + ȧ²/a²)  (2.37)", sp.simplify(Rs - 6 * (add / a + ad**2 / a**2)) == 0)
G00 = sp.simplify(R[0, 0] - g[0, 0] * Rs / 2)
check("G_00 = R_00 + R/2 = 3ȧ²/a²  (2.38 → 2.39)", sp.simplify(G00 - 3 * ad**2 / a**2) == 0)
Gn, rho, P, H0, rcr = sp.symbols('G_N rho P H0 rho_cr', positive=True)
Hsq = sp.solve(sp.Eq(3 * sp.Symbol('Hs'), 8 * sp.pi * Gn * rho), sp.Symbol('Hs'))[0]
check("(ȧ/a)^2 = 8πGρ/3  (2.39)", sp.simplify(Hsq - 8 * sp.pi * Gn * rho / 3) == 0)
check("หารด้วย H0² และใช้ ρcr = 3H0²/8πG ⇒ H²/H0² = ρ/ρcr  (2.40)",
      sp.simplify((Hsq / H0**2).subs(rho, rho) - rho / (3 * H0**2 / (8 * sp.pi * Gn))) == 0)

# ───────── ③ geodesic ของอนุภาคไร้มวล (2.26)–(2.29) ─────────
E = sp.Function('E')(t)
P1, P2, P3 = sp.symbols('P1 P2 P3')
Pv = [P1, P2, P3]
rhs = -sum(G[0][i + 1][j + 1] * Pv[i] * Pv[j] for i in range(3) for j in range(3))
check("−Γ^0_ij P^i P^j = −a ȧ δ_ij P^i P^j", sp.simplify(rhs + a * ad * (P1**2 + P2**2 + P3**2)) == 0)
# (2.27): a² δ_ij P^i P^j = E²  ⇒ RHS = −(ȧ/a) E²
rhs_E = sp.simplify(rhs.subs(P1, sp.sqrt(E**2 / a**2 - P2**2 - P3**2)))
check("แทน (2.27): RHS = −(ȧ/a)E²", sp.simplify(rhs_E + ad / a * E**2) == 0)
sol = sp.dsolve(sp.Eq(sp.diff(E, t) + ad / a * E, 0), E)
check("dE/dt + (ȧ/a)E = 0 ⇒ E ∝ 1/a  (2.28→2.29)", sp.simplify(sol.rhs * a).free_symbols <= {sp.Symbol('C1')})

# ───────── ④ conservation (2.52)–(2.58) ─────────
rh = sp.Function('rho')(t)
Pp = sp.Function('P')(t)
Tud = sp.diag(-rh, Pp, Pp, Pp)        # T^μ_ν (2.51)
nu = 0
div = sum(sp.diff(Tud[m, nu], X[m]) for m in range(4)) \
    + sum(G[m][al][m] * Tud[al, nu] for m in range(4) for al in range(4)) \
    - sum(G[al][nu][m] * Tud[m, al] for m in range(4) for al in range(4))
check("T^μ_0;μ = 0  ⇔  ρ̇ + (ȧ/a)(3ρ + 3P) = 0  (2.55)",
      sp.simplify(-div - (sp.diff(rh, t) + ad / a * (3 * rh + 3 * Pp))) == 0)
lhs256 = a**-3 * sp.diff(rh * a**3, t)
check("(2.55) ⇔ a^-3 ∂(ρa³)/∂t = −3(ȧ/a)P  (2.56)",
      sp.simplify(lhs256 - (-3 * ad / a * Pp) - (sp.diff(rh, t) + ad / a * (3 * rh + 3 * Pp))) == 0)
check("radiation P=ρ/3: ρ̇ + 4(ȧ/a)ρ = a^-4 ∂(ρa⁴)/∂t  (2.58)",
      sp.simplify(sp.diff(rh, t) + 4 * ad / a * rh - a**-4 * sp.diff(rh * a**4, t)) == 0)
# (2.85) ที่ w คงที่: ρ ∝ a^{-3(1+w)}
A, w = sp.symbols('A w')
aa = sp.symbols('a', positive=True)
rho_w = sp.exp(-3 * sp.integrate((1 + w) / aa, aa))
for wv, pw in ((sp.Rational(1, 3), -4), (0, -3), (-1, 0)):
    check(f"(2.85) w={wv}: ρ ∝ a^{pw}", sp.simplify(rho_w.subs(w, wv) - aa**pw) == 0)
# ตรวจว่า (2.85) แก้ (2.55) จริง (ทั่วไป w คงที่)
rw = aa**(-3 * (1 + w))
check("ρ ∝ a^{-3(1+w)} แก้ dρ/da + (3/a)(1+w)ρ = 0", sp.simplify(sp.diff(rw, aa) + 3 / aa * (1 + w) * rw) == 0)

# ───────── ⑤ สถิติ (2.59)–(2.69) ─────────
xx = sp.symbols('x', positive=True)
# sympy อินทิเกรตตัวนี้ตรง ๆ ไม่ออก ⇒ ตรวจสองชั้น: 6ζ(4) = π⁴/15 แบบสัญลักษณ์ + อินทิกรัลเชิงตัวเลข
check("6ζ(4) = π⁴/15 (สัญลักษณ์)", sp.simplify(6 * sp.zeta(4) - sp.pi**4 / 15) == 0)
check("∫₀^∞ x³/(eˣ−1) dx = π⁴/15 (ตัวเลข)",
      abs(mp.quad(lambda u: u**3 / (mp.e**u - 1), [0, mp.inf]) - mp.pi**4 / 15) < 1e-12)
I_BE = sp.pi**4 / 15
T = sp.symbols('T', positive=True)
rho_g = 8 * sp.pi * T**4 / (2 * sp.pi)**3 * I_BE      # (2.68)
check("ρ_γ = π²T⁴/15  (2.69)", sp.simplify(rho_g - sp.pi**2 * T**4 / 15) == 0)
I_FD = sp.Rational(7, 8) * sp.pi**4 / 15
check("Fermi–Dirac: ∫x³/(eˣ+1) = (7/8)π⁴/15 (ตรวจเชิงตัวเลข)",
      abs(mp.quad(lambda u: u**3 / (mp.e**u + 1), [0, mp.inf]) - float(I_FD)) < 1e-12)
# P = ρ/3 สำหรับไร้มวล: E = p ⇒ p²/(3E) = E/3 (2.62 เทียบ 2.59)
pp = sp.symbols('p', positive=True)
check("ไร้มวล: p²/(3E) = E/3 ⇒ P = ρ/3", sp.simplify(pp**2 / (3 * pp) - pp / 3) == 0)
# (2.63) กับ P = ρ/3, ρ = cT⁴
c = sp.symbols('c', positive=True)
check("ρ=cT⁴, P=ρ/3 ⇒ ∂P/∂T = (ρ+P)/T  (2.63)",
      sp.simplify(sp.diff(c * T**4 / 3, T) - (c * T**4 + c * T**4 / 3) / T) == 0)

# ───────── ⑥ entropy (2.64)–(2.66): เอกลักษณ์ที่ใช้ใน (2.65) ─────────
Tf = sp.Function('T')(t)
S = sp.Function('S')(t)      # S = ρ + P
lhs = a**-3 * sp.diff(S * a**3, t) - sp.diff(Tf, t) * S / Tf
rhs = a**-3 * Tf * sp.diff(S * a**3 / Tf, t)
check("a⁻³∂[(ρ+P)a³]/∂t − (Ṫ/T)(ρ+P) = a⁻³T ∂[(ρ+P)a³/T]/∂t  (2.65)", sp.simplify(lhs - rhs) == 0)

# ───────── ⑦ ตัวเลข: ρcr ในหน่วย eV⁴ และ (2.70) ─────────
eV_per_g = 5.6096e32          # นอกหนังสือ: 1 g c² = 5.6096e32 eV
hbarc = 1.97327e-5            # นอกหนังสือ: ħc = 1.97327e-5 eV cm
rcr_eV4 = 1.878e-29 * eV_per_g * hbarc**3
check(f"ρcr = {rcr_eV4:.4e} h² eV⁴ ≈ 8.098e-11", abs(rcr_eV4 / 8.098e-11 - 1) < 3e-3)
T0 = 2.725 / 11605
Og = float(mp.pi**2 / 15 * T0**4 / 8.098e-11)
check(f"Ω_γ h² = {Og:.3e} ≈ 2.47e-5  (2.70)", abs(Og / 2.47e-5 - 1) < 5e-3)

# ───────── ⑧ neutrino (2.73)–(2.78), (2.86)–(2.88) ─────────
g1 = 2 + sp.Rational(7, 8) * (2 + 2 + 3 + 3)
check("2 + (7/8)·10 = 43/4", g1 == sp.Rational(43, 4))
check("s(a1) = (2π²/45)(43/4)T³ = (43π²/90)T³  (2.73)",
      sp.simplify(2 * sp.pi**2 / 45 * g1 - 43 * sp.pi**2 / 90) == 0)
Tg, Tn = sp.symbols('T_gamma T_nu', positive=True)
# (2.74): s(a2) = 2π²/45 [2Tγ³ + (7/8)·6 Tν³] ; ตั้ง (a1T1)³ = (a2Tν)³ ⇒ 43/4 Tν³ = 2Tγ³ + 21/4 Tν³
ratio = sp.solve(sp.Eq(g1 * Tn**3, 2 * Tg**3 + sp.Rational(7, 8) * 6 * Tn**3), Tn)[0] / Tg
check("Tν/Tγ = (4/11)^{1/3}  (2.76)", sp.simplify(ratio - sp.Rational(4, 11)**sp.Rational(1, 3)) == 0)
check("(2.75): 43/2 (a1T1)³ = 4[(Tγ/Tν)³ + 21/8](Tν a2)³ คือ (2.73)=(2.74) คูณ 45/π²",
      sp.simplify(sp.Rational(43, 2) - 4 * (sp.Rational(11, 4) + sp.Rational(21, 8))) == 0)
fac = 3 * sp.Rational(7, 8) * sp.Rational(4, 11)**sp.Rational(4, 3)
check(f"ρν/ργ = 3·(7/8)(4/11)^(4/3) = {float(fac):.4f} ; ต่อชนิด {float(fac/3):.3f} ≈ 0.23",
      abs(float(fac / 3) - 0.227) < 1e-3)
check(f"Ων h² = {float(fac)*2.47e-5:.3e} ≈ 1.68e-5  (2.78)", abs(float(fac) * 2.47e-5 / 1.68e-5 - 1) < 3e-3)
check("Ω_r h² = 2.47e-5 + 1.68e-5 = 4.15e-5  (2.86)", abs(2.47e-5 + 1.68e-5 - 4.15e-5) < 1e-9)
Om, hh = sp.symbols('Omega_m h', positive=True)
aeq = sp.solve(sp.Eq(4.15e-5 / (hh**2 * aa**4), Om / aa**3), aa)[0]
check("Ω_m a⁻³ = 4.15e-5/(h²a⁴) ⇒ a_eq = 4.15e-5/(Ω_m h²)  (2.87)",
      sp.simplify(aeq - 4.15e-5 / (Om * hh**2)) == 0)
check(f"1 + z_eq = 1/a_eq = {1/4.15e-5:.3e} Ω_m h² ≈ 2.4e4  (2.88)", abs(1 / 4.15e-5 / 2.4e4 - 1) < 0.01)
# (2.80): n_ν = (3/11) n_γ, n_γ = 2ζ(3)T³/π²
ng = float(2 * mp.zeta(3) / mp.pi**2 * T0**3 / hbarc**3)
check(f"n_γ = {ng:.0f} cm⁻³ (≈ 411 ที่หนังสือใช้ในเฉลย)", abs(ng - 411) < 2)
nnu = 3 / 11 * ng
denom = 1.878e-29 * eV_per_g / nnu
check(f"Ων = mν/({denom:.1f} h² eV) ≈ mν/(94 h² eV)  (2.80)", abs(denom - 94) < 1)

# ───────── ⑨ ระยะทาง (2.41)–(2.50), (2.82)–(2.83) ─────────
ap = sp.symbols("a'", positive=True)
chi_MD = sp.integrate(1 / (ap**2 * H0 * ap**sp.Rational(-3, 2)), (ap, aa, 1))
check("MD: χ(a) = (2/H0)(1 − a^{1/2})", sp.simplify(chi_MD - 2 / H0 * (1 - sp.sqrt(aa))) == 0)
zz = sp.symbols('z', positive=True)
chi_z = sp.simplify(chi_MD.subs(aa, 1 / (1 + zz)))
check("แทน a = 1/(1+z) ⇒ (2/H0)[1 − 1/√(1+z)]  (2.43)",
      sp.simplify(chi_z - 2 / H0 * (1 - 1 / sp.sqrt(1 + zz))) == 0)
check("z เล็ก: χ ≈ z/H0", sp.simplify(sp.series(chi_z, zz, 0, 2).removeO() - zz / H0) == 0)
check("z → ∞: χ → 2/H0", sp.limit(chi_z, zz, sp.oo) == 2 / H0)
eta_MD = sp.integrate(1 / (ap**2 * H0 * ap**sp.Rational(-3, 2)), (ap, 0, aa))
eta_RD = sp.integrate(1 / (ap**2 * H0 * ap**-2), (ap, 0, aa))
check("η ∝ a^{1/2} (matter) · η ∝ a (radiation)  (หน้า 34)",
      sp.simplify(eta_MD / sp.sqrt(aa)).free_symbols == {H0} and sp.simplify(eta_RD / aa).free_symbols == {H0})
Ok, chi = sp.symbols('Omega_k chi', positive=True)
dA_open = aa / (H0 * sp.sqrt(Ok)) * sp.sinh(sp.sqrt(Ok) * H0 * chi)
check("(2.46) Ωk → 0 กลับเป็น aχ", sp.limit(dA_open, Ok, 0) == aa * chi)
L, F = sp.symbols('L F', positive=True)
dL = sp.sqrt(L / (4 * sp.pi * (L * aa**2 / (4 * sp.pi * chi**2))))
check("F = La²/(4πχ²) และ F = L/(4π d_L²) ⇒ d_L = χ/a  (2.49→2.50)", sp.simplify(dL - chi / aa) == 0)
check("d_L/d_A = 1/a² = (1+z)²", sp.simplify((chi / aa) / (aa * chi) - aa**-2) == 0)
v = 10**((24.32 - 16.08) / 5) * 0.026
check(f"H0 d_L(0.83) = 10^((24.32−16.08)/5)·0.026 = {v:.3f} ≈ 1.16  (2.83)", abs(v - 1.16) < 0.01)
zs = 0.83
dl_md = 2 * (1 + zs) * (1 - 1 / mp.sqrt(1 + zs))
check(f"Ωm=1: H0 d_L(0.83) = {float(dl_md):.3f} ≈ 0.95", abs(dl_md - 0.95) < 0.01)
dl_L = (1 + zs) * mp.quad(lambda q: 1 / mp.sqrt(0.3 * (1 + q)**3 + 0.7), [0, zs])
check(f"Ωm=0.3, ΩΛ=0.7: H0 d_L(0.83) = {float(dl_L):.3f} ≈ 1.23", abs(dl_L - 1.23) < 0.01)

print(f"\nผ่านทั้งหมด {ok} ข้อ")
