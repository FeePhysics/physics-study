"""คู่มืออ่าน Dodelson บทที่ 11 · ตรวจขั้นพีชคณิตที่คู่มือเติมเอง (sympy + mpmath)

รัน:  python scripts/verify/dodelson_ch11.py
ทุกบรรทัด ✅ = ผ่าน · assert ล้ม = มีอะไรผิด
เลข (11.x) = เลขสมการในหนังสือ
"""
import sympy as sp
import mpmath as mp

ok = 0


def check(name, cond):
    global ok
    assert cond, name
    ok += 1
    print("✅", name)


# ── ① §11.1.1 ชั่งน้ำหนัก: (11.6)-(11.13) ────────────────────────────
N = 5                                            # ใช้ N จำกัดให้ sympy กระจายได้
d = sp.symbols(f"d1:{N + 1}", real=True)
w = sp.symbols("w", real=True)
s2 = sp.symbols("sigma2", positive=True)         # σ_w^2
lnL = -sp.Rational(N, 2) * sp.log(2 * sp.pi * s2) - sum((di - w)**2 for di in d) / (2 * s2)
wbar = sum(d) / N
check("(11.8) ∂lnL/∂w = 0 ⇒ w = ค่าเฉลี่ย", sp.simplify(sp.solve(sp.diff(lnL, w), w)[0] - wbar) == 0)
s2hat = sp.solve(sp.diff(lnL, s2).subs(w, wbar), s2)[0]
check("(11.10) σ_w^2 ที่ดีที่สุด = Σ(d-w̄)^2/N",
      sp.simplify(s2hat - sum((di - wbar)**2 for di in d) / N) == 0)
curv = sp.diff(lnL, w, 2)
check("(11.11) ∂²lnL/∂w² = -N/σ_w² ⇒ error = σ_w/√N", sp.simplify(curv + N / s2) == 0)
# ln L เป็นกำลังสองใน w เป๊ะ ⇒ Taylor ถึงอันดับ 2 คือทั้งหมด
check("(11.11) lnL(w) = lnL(w̄) - N(w-w̄)²/(2σ²) เป๊ะ",
      sp.simplify(lnL - (lnL.subs(w, wbar) - N * (w - wbar)**2 / (2 * s2))) == 0)

# ── ② (11.16)-(11.19): convolution ของ Gaussian สองตัว ────────────────
s, D = sp.symbols("s Delta", real=True)
CS, CN = sp.symbols("C_S C_N", positive=True)
integrand = sp.exp(-s**2 / (2 * CS)) / sp.sqrt(2 * sp.pi * CS) * \
    sp.exp(-(D - s)**2 / (2 * CN)) / sp.sqrt(2 * sp.pi * CN)
Lconv = sp.simplify(sp.integrate(integrand, (s, -sp.oo, sp.oo)))
target = sp.exp(-D**2 / (2 * (CS + CN))) / sp.sqrt(2 * sp.pi * (CS + CN))
check("(11.19) ∫ds N(s;0,C_S) N(Δ;s,C_N) = N(Δ;0,C_S+C_N)", sp.simplify(Lconv / target) == 1)
# (11.17) การทำกำลังสองสมบูรณ์
C = CS + CN
lhs17 = -s**2 * C / (2 * CS * CN) + D * s / CN - D**2 / (2 * CN)
rhs17 = -C / (2 * CS * CN) * (s - CS * D / C)**2 + CS * D**2 / (2 * C * CN) - D**2 / (2 * CN)
check("(11.17) ทำกำลังสองสมบูรณ์", sp.simplify(lhs17 - rhs17) == 0)

# ── ③ (11.21)-(11.26): Np pixel, C = (C_S + C_N) I ───────────────────
Np, S2sum = sp.symbols("N_p Sigma", positive=True)      # Sigma = Σ Δ_i^2
lnL2 = -Np / 2 * sp.log(CS + CN) - S2sum / (2 * (CS + CN))
CShat = sp.solve(sp.diff(lnL2, CS), CS)[0]
check("(11.23) C_S ที่ดีที่สุด = ΣΔ²/Np - C_N", sp.simplify(CShat - (S2sum / Np - CN)) == 0)
d2 = sp.diff(lnL2, CS, 2)
check("(11.25) อนุพันธ์อันดับสอง",
      sp.simplify(d2 - (Np / 2 / (CS + CN)**2 - S2sum / (CS + CN)**3)) == 0)
sig = sp.simplify((-d2.subs(S2sum, Np * (CS + CN)))**sp.Rational(-1, 2))
check("(11.26) σ_CS = √(2/Np)(C_S+C_N)", sp.simplify(sig - sp.sqrt(2 / Np) * (CS + CN)) == 0)

# ── ④ (11.45)-(11.51): window function ของ beam ───────────────────────
x, y, lx, ly = sp.symbols("x y l_x l_y", real=True)
sg = sp.symbols("sigma", positive=True)
B = sp.exp(-(x**2 + y**2) / (2 * sg**2)) / (2 * sp.pi * sg**2)
Bt = sp.integrate(sp.integrate(B * sp.exp(-sp.I * (lx * x + ly * y)), (x, -sp.oo, sp.oo)), (y, -sp.oo, sp.oo))
check("(11.45) FT ของ Gaussian beam = exp(-l²σ²/2)",
      sp.simplify(Bt - sp.exp(-(lx**2 + ly**2) * sg**2 / 2)) == 0)
fwhm_factor = 1 / sp.sqrt(8 * sp.log(2))
check("σ = FWHM/√(8 ln2) ≈ 0.4247 FWHM (หนังสือพิมพ์ 0.4245)",
      abs(float(fwhm_factor) - 0.42466) < 1e-4)
x0 = sp.symbols("x0", positive=True)
# (11.47) chop สองจุด: FT = e^{-i lx x0} - e^{+i lx x0} = -2i sin(lx x0) ; |B|^2 = 4 sin^2
Bchop = sp.exp(-sp.I * lx * x0) - sp.exp(sp.I * lx * x0)
check("(11.48) |B̃|² = 4 sin²(l_x x0) (เครื่องหมายหน้า 2i ไม่มีผล)",
      sp.simplify(sp.expand_complex(Bchop * sp.conjugate(Bchop)) - 4 * sp.sin(lx * x0)**2) == 0)
mp.mp.dps = 20
for l, xv in [(30, 0.02), (120, 0.01), (300, 0.004)]:
    Wnum = mp.quad(lambda p: 4 * mp.sin(l * xv * mp.cos(p))**2, [0, 2 * mp.pi]) / (2 * mp.pi)
    check(f"(11.49) W = 2(1-J0(2 l x0)) ที่ l={l}", abs(Wnum - 2 * (1 - mp.besselj(0, 2 * l * xv))) < 1e-12)
    # (11.41) เป็นการประมาณ l ใหญ่: คลาดราว 1/l (ถ้าใช้ l+½ แทน l จะแม่นกว่ามาก — นอกหนังสือ)
    check(f"P_l(cos 2x0) ≈ J0(2 l x0) คลาดไม่เกิน 2/l ที่ l={l}",
          abs(mp.legendre(l, mp.cos(2 * xv)) - mp.besselj(0, 2 * l * xv)) < 2 / l)
    check(f"P_l(cos 2x0) ≈ J0(2(l+½)x0) คลาด < 1e-3 ที่ l={l}",
          abs(mp.legendre(l, mp.cos(2 * xv)) - mp.besselj(0, 2 * (l + 0.5) * xv)) < 1e-3)

# ── ⑤ (11.57)-(11.59): survey ทรงกลมรัศมี R ───────────────────────────
q, R, r = sp.symbols("q R r", positive=True)
mu = sp.symbols("mu")
ft = sp.integrate(sp.integrate(2 * sp.pi * r**2 * sp.exp(sp.I * q * r * mu), (mu, -1, 1)), (r, 0, R))
check("(11.57) ∫_{|x|<R} e^{iq·x} d³x = 4π(sin qR - qR cos qR)/q³",
      sp.simplify(sp.expand_complex(ft) - 4 * sp.pi * (sp.sin(q * R) - q * R * sp.cos(q * R)) / q**3) == 0)
yv = sp.symbols("y", positive=True)
j1 = (sp.sin(yv) - yv * sp.cos(yv)) / yv**2
check("(sin y - y cos y)²/y⁶ = j1(y)²/y²", sp.simplify((sp.sin(yv) - yv * sp.cos(yv))**2 / yv**6 - j1**2 / yv**2) == 0)
# แปลง dμ → dy: y² = (k²+ki²-2 k ki μ)R²  ⇒ |dμ| = y dy/(k ki R²) ; V = 4πR³/3
kk, ki = sp.symbols("k k_i", positive=True)
V = 4 * sp.pi * R**3 / 3
pref = 8 * sp.pi**2 * R**6 / V**2 / (kk * ki * R**2)     # (8π²R⁶/V²)·(1/(k ki R²)) × ∫ y dy j1²/y²
check("(11.58)→(11.59) สัมประสิทธิ์ = 9/(2 k k_i R²)", sp.simplify(pref - 9 / (2 * kk * ki * R**2)) == 0)
mp.mp.dps = 20
kR, kiR = 2.3, 3.1      # ตรวจตัวเลข: เฉลี่ยมุมตรง ๆ เทียบกับรูป (11.59)
j1n = lambda t: (mp.sin(t) - t * mp.cos(t)) / t**2
lhs = mp.quad(lambda m: 0.5 * (j1n(mp.sqrt(kR**2 + kiR**2 - 2 * kR * kiR * m)))**2
              / (kR**2 + kiR**2 - 2 * kR * kiR * m), [-1, 1]) * (4 * mp.pi)**2 / (4 * mp.pi / 3)**2
rhs = 9 / (2 * kR * kiR) * mp.quad(lambda t: j1n(t)**2 / t, [abs(kR - kiR), kR + kiR])
check("(11.58)=(11.59) ตรวจตัวเลข (R=1)", abs(lhs - rhs) < 1e-12)

# ── ⑥ Karhunen-Loève ตัวอย่างสอง pixel (11.74)-(11.82) ─────────────────
eps, th = sp.symbols("epsilon theta", real=True)
sn, ss = sp.symbols("sigma_n sigma_s", positive=True)
Rth = sp.Matrix([[sp.cos(th), sp.sin(th)], [-sp.sin(th), sp.cos(th)]])
M = sp.Matrix([[1, eps], [eps, 1]])
RMRt = sp.simplify(Rth * M * Rth.T)
check("(11.79) RMR^T = [[1+ε sin2θ, ε cos2θ],[ε cos2θ, 1-ε sin2θ]]",
      sp.simplify(RMRt - sp.Matrix([[1 + eps * sp.sin(2 * th), eps * sp.cos(2 * th)],
                                    [eps * sp.cos(2 * th), 1 - eps * sp.sin(2 * th)]])) == sp.zeros(2, 2))
CSm = ss**2 * M
CNm = sn**2 * sp.eye(2)
R2 = sp.eye(2) / sn
R3 = Rth.subs(th, sp.pi / 4)
CSp = sp.simplify(R3 * R2 * CSm * R2.T * R3.T)
CNp = sp.simplify(R3 * R2 * CNm * R2.T * R3.T)
check("(11.80) C_S' = (σs²/σn²) diag(1+ε, 1-ε)",
      sp.simplify(CSp - ss**2 / sn**2 * sp.diag(1 + eps, 1 - eps)) == sp.zeros(2, 2))
check("หลังหมุน C_N' ยังเป็น I", CNp == sp.eye(2))
D1, D2 = sp.symbols("Delta1 Delta2")
modes = R3 * R2 * sp.Matrix([D1, D2])
check("(11.82) โหมดใหม่ = (Δ1+Δ2)/(√2σn), (-Δ1+Δ2)/(√2σn)",
      sp.simplify(modes - sp.Matrix([(D1 + D2), (-D1 + D2)]) / (sp.sqrt(2) * sn)) == sp.zeros(2, 1))

# ── ⑦ (11.89), (11.93): อนุพันธ์ของ lnL แบบ Gaussian หลายมิติ ─────────
# ใช้เมทริกซ์ตรรกยะ 3×3 ที่บวกแน่นอน ⇒ ตรวจแบบเป๊ะ ไม่ใช่ตัวเลขทศนิยม
lamS = sp.symbols("lambda", real=True)
C0 = sp.Matrix([[4, 1, 0], [1, 3, 1], [0, 1, 5]])
C1 = sp.Matrix([[2, 1, 1], [1, 2, 0], [1, 0, 1]])
Cl = C0 + lamS * C1
Dv = sp.Matrix([1, -2, 3])
lnLf = -sp.log(Cl.det()) / 2 - (Dv.T * Cl.inv() * Dv)[0] / 2
Ci = Cl.inv()
ana1 = (Dv.T * Ci * C1 * Ci * Dv)[0] / 2 - (Ci * C1).trace() / 2
check("(11.89) (lnL),λ = ½ΔC⁻¹C,λC⁻¹Δ - ½Tr[C⁻¹C,λ]", sp.simplify(sp.diff(lnLf, lamS) - ana1) == 0)
ana2 = -(Dv.T * Ci * C1 * Ci * C1 * Ci * Dv)[0] + (Ci * C1 * Ci * C1).trace() / 2   # C,λλ = 0
check("(11.90) (lnL),λλ (C เชิงเส้นใน λ ⇒ C,λλ = 0)", sp.simplify(sp.diff(lnLf, lamS, 2) - ana2) == 0)
# ค่าเฉลี่ยเมื่อ <ΔΔ^T> = C : <Δ^T M Δ> = Tr[M C]
meanQ = (Ci * C1 * Ci * C1 * Ci * Cl).trace()
F = -(-meanQ + (Ci * C1 * Ci * C1).trace() / 2)
check("(11.93) F = <-(lnL),λλ> = ½Tr[C,λC⁻¹C,λC⁻¹]", sp.simplify(F - (C1 * Ci * C1 * Ci).trace() / 2) == 0)
check("(C⁻¹),λ = -C⁻¹C,λC⁻¹", sp.simplify(sp.diff(Ci, lamS) + Ci * C1 * Ci) == sp.zeros(3, 3))
check("∂ ln det C/∂λ = Tr[C⁻¹C,λ]", sp.simplify(sp.diff(sp.log(Cl.det()), lamS) - (Ci * C1).trace()) == 0)

# ── ⑧ Fisher ตัวอย่างง่าย: ชั่งน้ำหนักซ้ำ ใช้ (11.95) ให้ผลเดียวกับ §11.1.1 ──
# ข้อมูล N ตัว ค่าเฉลี่ย w ความแปรปรวน σ² : Fisher ของ w = N/σ² ⇒ error σ/√N
Fw = -sp.diff(lnL, w, 2)
check("Fisher ของ w ในตัวอย่างชั่งน้ำหนัก = N/σ²", sp.simplify(Fw - N / s2) == 0)
# Fisher ของ σ² (ข้อ 2 ท้ายบท ทำแบบนี้ได้): <-∂²lnL/∂(σ²)²> โดย <(d-w)²> = σ²
SS = sp.symbols("SS", positive=True)                          # SS = Σ(d-w)²
Fs = -sp.diff(-sp.Rational(N, 2) * sp.log(s2) - SS / (2 * s2), s2, 2)
check("-∂²lnL/∂(σ²)² = -N/(2σ⁴) + Σ(d-w)²/σ⁶", sp.simplify(Fs - (-sp.Rational(N, 2) / s2**2 + SS / s2**3)) == 0)
Fs_direct = Fs.subs(SS, N * s2)
check("Fisher ของ σ² = N/(2σ⁴) (ใช้ <Σ(d-w)²> = Nσ²)", sp.simplify(Fs_direct - sp.Rational(N, 2) / s2**2) == 0)
check("ตัวอย่าง Fisher ด้วยสูตร (11.95) ½Tr[...] : C = σ² I ⇒ F_σ²σ² = N/(2σ⁴)",
      sp.simplify(sp.Rational(1, 2) * sp.trace((sp.eye(N) / s2) * (sp.eye(N) / s2)) - sp.Rational(N, 2) / s2**2) == 0)

# ── ⑨ (11.109)-(11.110): Fisher ของ C_l ทั้งท้องฟ้า ────────────────────
l = 7
lam, ww, sb = sp.symbols("lambda_l w sigma_b", positive=True)
bl = sp.exp(-l**2 * sb**2)
# บล็อก (l,m) ของ C เป็นเส้นทแยง: λ b + 1/w ; อนุพันธ์ = b
Fll = sp.Rational(1, 2) * sum(bl**2 / (lam * bl + 1 / ww)**2 for _ in range(2 * l + 1))
check("(11.109) F_ll = (2l+1)/2 · e^{-2l²σ²}[C_l e^{-l²σ²} + 1/w]^{-2}",
      sp.simplify(Fll - sp.Rational(2 * l + 1, 2) * sp.exp(-2 * l**2 * sb**2) / (lam * bl + 1 / ww)**2) == 0)
check("(11.110) δC_l = √(2/(2l+1)) [C_l + e^{l²σ²}/w]",
      sp.simplify(1 / sp.sqrt(Fll) - sp.sqrt(sp.Rational(2, 2 * l + 1)) * (lam + sp.exp(l**2 * sb**2) / ww)) == 0)

# ── ⑩ (11.126)-(11.127): marginalize ⇒ error = √(F⁻¹)₁₁ ────────────────
F11, F22, F12 = sp.symbols("F11 F22 F12", real=True)
l1, l2 = sp.symbols("lambda1 lambda2", real=True)
Fm = sp.Matrix([[F11, F12], [F12, F22]])
expo = -sp.Rational(1, 2) * (sp.Matrix([[l1, l2]]) * Fm * sp.Matrix([l1, l2]))[0]
# ทำกำลังสองสมบูรณ์ใน λ2
compl = -F22 / 2 * (l2 + F12 * l1 / F22)**2 - l1**2 / 2 * (F11 * F22 - F12**2) / F22
check("(11.127) ทำกำลังสองสมบูรณ์ใน λ2", sp.simplify(sp.expand(expo - compl)) == 0)
check("(F11F22-F12²)/F22 = 1/(F⁻¹)₁₁", sp.simplify((F11 * F22 - F12**2) / F22 - 1 / Fm.inv()[0, 0]) == 0)
# ตัวเลข: F = [[4,2],[2,3]] ⇒ ไม่ marginalize 1/√4 = 0.5 ; marginalize √(3/8)
Fn = sp.Matrix([[4, 2], [2, 3]])
check("ตัวอย่าง F=[[4,2],[2,3]]: σ1 (λ2 รู้ค่า) = 1/2", 1 / sp.sqrt(Fn[0, 0]) == sp.Rational(1, 2))
check("ตัวอย่าง F=[[4,2],[2,3]]: σ1 (marginalize) = √(3/8)", sp.sqrt(Fn.inv()[0, 0]) == sp.sqrt(sp.Rational(3, 8)))
ll1, ll2 = sp.symbols("u v", real=True)
marg = sp.integrate(sp.exp(-sp.Rational(1, 2) * (4 * ll1**2 + 4 * ll1 * ll2 + 3 * ll2**2)), (ll2, -sp.oo, sp.oo))
check("อินทิเกรต λ2 จริง ⇒ exp(-λ1²/(2·3/8))",
      sp.simplify(sp.log(marg).diff(ll1, 2) + 1 / sp.Rational(3, 8)) == 0)

# ── ⑪ (11.122)-(11.124): χ² ⇒ Fisher ─────────────────────────────────
a_, b_ = sp.symbols("a b", real=True)
ls = [2, 3, 4]
Cl_model = lambda A: [A * (lv + 1) for lv in ls]                # โมเดลของเล่น C_l = A(l+1)
dC = [sp.Rational(1, 2), 1, 2]
Ao = sp.Integer(3)
obs = Cl_model(Ao)                                               # ข้อมูลตรงกับโมเดลพอดี ⇒ เทอมที่สองหาย
chi2 = sum((cm - co)**2 / dc**2 for cm, co, dc in zip(Cl_model(a_), obs, dC))
Fcal = sp.diff(chi2, a_, 2) / 2
F124 = sum((sp.diff(cm, a_))**2 / dc**2 for cm, dc in zip(Cl_model(a_), dC))
check("(11.122)=(11.124) เมื่อ C_l^obs ตรงโมเดล", sp.simplify(Fcal - F124) == 0)

# ── ⑫ mapmaking (11.129)-(11.137) ────────────────────────────────────
P = sp.Matrix([[1, 0], [1, 0], [0, 1], [1, 0], [0, 1]])          # 5 ช่วงเวลา, 2 pixel
dt = sp.Matrix(sp.symbols("d1:6"))
Nn = sp.symbols("N", positive=True)
Ninv = sp.eye(5) / Nn
CNmap = (P.T * Ninv * P).inv()
Dmap = CNmap * P.T * Ninv * dt
check("(11.136) C_N = N/π_i (π = 3, 2 ครั้ง)", sp.simplify(CNmap - sp.diag(Nn / 3, Nn / 2)) == sp.zeros(2, 2))
check("(11.137) map = ค่าเฉลี่ยของข้อมูลที่ชี้ pixel นั้น",
      sp.simplify(Dmap - sp.Matrix([(dt[0] + dt[1] + dt[3]) / 3, (dt[2] + dt[4]) / 2])) == sp.zeros(2, 1))
sv = sp.Matrix(sp.symbols("s1:3"))
chi2m = ((dt - P * sv).T * Ninv * (dt - P * sv))[0]
sol = sp.solve([sp.diff(chi2m, si) for si in sv], list(sv), dict=True)[0]
check("(11.130)-(11.134) ∂χ²/∂s = 0 ให้ map เดียวกัน",
      sp.simplify(sp.Matrix([sol[si] for si in sv]) - Dmap) == sp.zeros(2, 1))
# noise ของ map: Δ - s = C_N P^T N⁻¹ η ⇒ cov = C_N P^T N⁻¹ N N⁻¹ P C_N = C_N
Ncov = CNmap * P.T * Ninv * (Nn * sp.eye(5)) * Ninv * P * CNmap
check("(11.135) covariance ของ noise ใน map = (P^T N⁻¹ P)⁻¹", sp.simplify(Ncov - CNmap) == sp.zeros(2, 2))

# ── ⑬ foreground สองความถี่ (11.147)-(11.153) ─────────────────────────
sn2 = sp.symbols("sigma_n2", positive=True)
dd1, dd2 = sp.symbols("d1 d2")
Nf = sn2 * sp.eye(2)
W0 = sp.Matrix([[1], [1]])
c0 = (W0.T * Nf.inv() * W0)[0]
check("(11.148) WN⁻¹W = 2/σn²", sp.simplify(c0 - 2 / sn2) == 0)
est0 = sp.simplify((W0.T * Nf.inv() * sp.Matrix([dd1, dd2]))[0] / c0)
check("(11.149) ไม่มี foreground ⇒ (d1+d2)/2", sp.simplify(est0 - (dd1 + dd2) / 2) == 0)
Wm = sp.Matrix([[1, 1], [1, sp.Rational(1, 2)]])                  # คอลัมน์ = CMB, synchrotron
CNinv = Wm.T * Nf.inv() * Wm
check("(11.151) C_N⁻¹ = (1/σn²)[[2,3/2],[3/2,5/4]]",
      sp.simplify(CNinv - sp.Matrix([[2, sp.Rational(3, 2)], [sp.Rational(3, 2), sp.Rational(5, 4)]]) / sn2) == sp.zeros(2, 2))
CNf = CNinv.inv()
check("(11.152) C_N = 4σn²[[5/4,-3/2],[-3/2,2]]",
      sp.simplify(CNf - 4 * sn2 * sp.Matrix([[sp.Rational(5, 4), -sp.Rational(3, 2)], [-sp.Rational(3, 2), 2]])) == sp.zeros(2, 2))
est = sp.simplify((CNf * Wm.T * Nf.inv() * sp.Matrix([dd1, dd2]))[0])
check("(11.153) ตัวประมาณ CMB = -d1 + 2d2", sp.simplify(est - (-dd1 + 2 * dd2)) == 0)
check("ตัวประมาณไม่เห็น foreground: ใส่ d = (1, 1/2) ได้ 0", est.subs({dd1: 1, dd2: sp.Rational(1, 2)}) == 0)
check("FDF = √(5σ² / (σ²/2)) = √10", sp.sqrt(CNf[0, 0] / (1 / CNinv[0, 0])) == sp.sqrt(10))

# ── ⑭ constraint matrix (11.154)-(11.158) ─────────────────────────────
kap, aa = sp.symbols("kappa a", positive=True)                    # a = σn² + σs²
Cc = sp.Matrix([[aa + kap, kap], [kap, aa + kap]])
check("(11.156) det C = 2κa + a² → 2κa", sp.simplify(Cc.det() - (2 * kap * aa + aa**2)) == 0)
Cinv_lim = Cc.inv().applyfunc(lambda e: sp.limit(e, kap, sp.oo))
check("(11.157) C⁻¹ → [[1,-1],[-1,1]]/(2a)", sp.simplify(Cinv_lim - sp.Matrix([[1, -1], [-1, 1]]) / (2 * aa)) == sp.zeros(2, 2))
quad = sp.simplify((sp.Matrix([[D1, D2]]) * Cinv_lim * sp.Matrix([D1, D2]))[0] / 2)
check("½ΔᵀC⁻¹Δ → (Δ1-Δ2)²/(4a)  [หนังสือพิมพ์ /(2a) ใน (11.158) — ตั้งเป็นธง]",
      sp.simplify(quad - (D1 - D2)**2 / (4 * aa)) == 0)
check("ความแปรปรวนของ Δ1-Δ2 เมื่อสอง pixel อิสระ = 2a", sp.simplify(sp.Matrix([[1, -1]]) * (aa * sp.eye(2)) * sp.Matrix([1, -1]))[0] == 2 * aa)

print(f"\nผ่านทั้งหมด {ok} ข้อ")
