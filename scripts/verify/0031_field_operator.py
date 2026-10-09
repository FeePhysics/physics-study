"""ตรวจทุกข้ออ้างของบทที่ 19 — สนามเขียนเป็นผลรวมของปุ่ม a_k, a_k†

1+1 มิติ · g = diag(-1, +1) · hbar = c = 1 · กล่องยาว L ต่อปลายเป็นวง ⇒ k = 2 pi n / L
phi(t, x) = sum_k 1/sqrt(2 w_k L) [ a_k e^{i(kx - w_k t)} + a_k† e^{-i(kx - w_k t)} ]
pi = d phi / dt

ตัวดำเนินการใช้ sympy noncommutative symbols แล้วจัดลำดับเองด้วย [a_k, a_q†] = delta
⇒ ไม่ได้พิมพ์ผลลัพธ์ไว้ก่อน (ผลของ H กับ [phi, pi] ออกมาจากการคูณจริง)
"""
import itertools
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


t, x, y, L, m = sp.symbols('t x y L m', real=True)
L = sp.Symbol('L', positive=True)
m = sp.Symbol('m', positive=True)

# ── ② แต่ละรูปร่างแก้สมการสนามเมื่อ w^2 = k^2 + m^2 ─────────────
k, w = sp.symbols('k w', real=True)
mode = sp.exp(sp.I*(k*x - w*t))
kg = sp.diff(mode, t, 2) - sp.diff(mode, x, 2) + m**2*mode
chk("phi_tt - phi_xx + m^2 phi = (-w^2 + k^2 + m^2) phi", sp.simplify(kg/mode), -w**2 + k**2 + m**2)
chk("แทน w = sqrt(k^2 + m^2) แล้วเป็นศูนย์", sp.simplify((kg/mode).subs(w, sp.sqrt(k**2 + m**2))), 0)
modec = sp.conjugate(mode).subs({sp.conjugate(k): k, sp.conjugate(w): w})
chk("พจน์ a† (สังยุค) ก็แก้สมการเดียวกัน",
    sp.simplify((sp.diff(modec, t, 2) - sp.diff(modec, x, 2) + m**2*modec)/modec), -w**2 + k**2 + m**2)
ok.append(("ตัวควบคุม: w ผิด (w = k) แล้วไม่เป็นศูนย์", sp.simplify((kg/mode).subs(w, k)) != 0, None, None))

# ── k ที่ใส่ในกล่องได้ + ตั้งฉากกัน ───────────────────────────
n1, n2 = sp.symbols('n1 n2', integer=True)
for a_, b_ in [(3, 3), (3, -3), (2, 5), (-4, -4), (0, 0), (1, 0)]:
    val = sp.integrate(sp.exp(sp.I*2*sp.pi*(a_ - b_)*x/L), (x, 0, L))
    chk(f"∫_0^L e^(i(k_{a_}-k_{b_})x) dx = {'L' if a_ == b_ else '0'}", sp.simplify(val), L if a_ == b_ else 0)
chk("e^{ikx} กลับมาค่าเดิมที่ x = L เมื่อ k = 2 pi n/L", sp.simplify(sp.exp(sp.I*2*sp.pi*n1/L*L)), 1)

# ── ตัวดำเนินการ: เขียน phi, pi เป็นผลรวมของพจน์ (สัมประสิทธิ์, ตัวดำเนินการ, k, เครื่องหมาย) ──
N = 2
ns = list(range(-N, N + 1))
A = {n: sp.Symbol(f'a{n}', commutative=False) for n in ns}
B = {n: sp.Symbol(f'b{n}', commutative=False) for n in ns}       # b = a†
kv = lambda n: 2*sp.pi*n/L
wv = lambda n: sp.sqrt(kv(n)**2 + m**2)
norm = lambda n: 1/sp.sqrt(2*wv(n)*L)


def field_terms(tt):
    """phi = sum c * op * e^{i s k x}  คืน [(c, op, n, s)] · s = +1 ของ a · -1 ของ a†"""
    out = []
    for n in ns:
        out.append((norm(n)*sp.exp(-sp.I*wv(n)*tt), A[n], n, +1))
        out.append((norm(n)*sp.exp(+sp.I*wv(n)*tt), B[n], n, -1))
    return out


def dt(terms, tt):
    return [(sp.diff(c, tt), op, n, s) for c, op, n, s in terms]


def dx(terms):
    return [(c*sp.I*s*kv(n), op, n, s) for c, op, n, s in terms]


def integrate_product(T1, T2):
    """∫_0^L T1 T2 dx — ใช้ความตั้งฉาก: เหลือคู่ที่ s1 k1 + s2 k2 = 0 คูณ L"""
    tot = 0
    for (c1, o1, n1_, s1), (c2, o2, n2_, s2) in itertools.product(T1, T2):
        if s1*n1_ + s2*n2_ == 0:
            tot += c1*c2*L*o1*o2
    return tot


def normal_order(expr):
    """จัดทุกพจน์สองตัวให้ b มาก่อน a ด้วย a_n b_n = b_n a_n + 1 · คนละ n สลับได้เลย"""
    expr = sp.expand(expr)
    res = 0
    for term in sp.Add.make_args(expr):
        c, nc = term.args_cnc()
        coef = sp.Mul(*c)
        if len(nc) == 2:
            o1, o2 = nc
            n_1, n_2 = int(o1.name[1:]), int(o2.name[1:])
            if o1.name[0] == 'a' and o2.name[0] == 'b':
                res += coef*o2*o1 + (coef if n_1 == n_2 else 0)
                continue
            if o1.name[0] == o2.name[0] and n_1 > n_2:        # a a หรือ b b คนละ n สลับเรียงได้
                res += coef*o2*o1
                continue
        res += term
    return sp.expand(res)


phi = field_terms(t)
pi_ = dt(phi, t)
H = normal_order(sp.Rational(1, 2)*integrate_product(pi_, pi_)
                 + sp.Rational(1, 2)*integrate_product(dx(phi), dx(phi))
                 + sp.Rational(1, 2)*m**2*integrate_product(phi, phi))
H_want = sum(wv(n)*(B[n]*A[n] + sp.Rational(1, 2)) for n in ns)
ok.append(("H = ∫(½pi² + ½phi_x² + ½m²phi²)dx = Σ w_k (a_k†a_k + ½)  [N=2 · t ทั่วไป]",
           sp.simplify(sp.expand(H - H_want)) == 0, None, None))
# พจน์ a_k a_{-k} ต้องหักกันหมด — ตัวควบคุม: เอาเฉพาะพลังงานจลน์ ยังเหลือพจน์คู่
Hk = normal_order(sp.Rational(1, 2)*integrate_product(pi_, pi_))
ok.append(("ตัวควบคุม: ½pi² อย่างเดียวยังมีพจน์ a_k a_{-k} ค้าง (ต้องรวมสามก้อนถึงหักกัน)",
           any(str(f).startswith('a') and 'a' in str(f)[1:] for f in sp.Add.make_args(Hk)
               for f in [sp.Mul(*f.args_cnc()[1])]), None, None))
P = normal_order(-integrate_product(pi_, dx(phi)))
P_want = sum(kv(n)*B[n]*A[n] for n in ns)
ok.append(("P = -∫ pi phi_x dx = Σ k a_k†a_k  (ค่าคงที่หักกันเพราะ k กับ -k)",
           sp.simplify(sp.expand(P - P_want)) == 0, sp.simplify(sp.expand(P - P_want)), 0))

# ── ③ กับสุญญากาศ: <0|phi|0> = 0 · phi|0> มีแต่สถานะหนึ่งอนุภาค ─────
# ทำบน Fock space จริง: แต่ละรูปร่างตัดที่ 3 ขั้น
D = 3
a1 = sp.zeros(D, D)
for j in range(1, D):
    a1[j - 1, j] = sp.sqrt(j)
modes = [-1, 0, 1]
I_ = sp.eye(D)


def kron(*Ms):
    out = Ms[0]
    for M in Ms[1:]:
        out = sp.kronecker_product(out, M)
    return out


aop = {q: kron(*[a1 if r == q else I_ for r in modes]) for q in modes}
vac = sp.zeros(D**len(modes), 1); vac[0] = 1
xx = sp.Symbol('xx', real=True)
phi0 = sum((norm(q)*(aop[q]*sp.exp(sp.I*kv(q)*xx) + aop[q].H*sp.exp(-sp.I*kv(q)*xx)) for q in modes),
           sp.zeros(D**len(modes)))
st = phi0*vac
chk("<0|phi(x)|0> = 0", (vac.T*st)[0], 0)
one = lambda q: aop[q].H*vac
for q in modes:
    chk(f"<k={q}|phi(x)|0> = e^(-ikx)/sqrt(2wL)", (one(q).T*st)[0], norm(q)*sp.exp(-sp.I*kv(q)*xx))
num = sum((aop[q].H*aop[q] for q in modes), sp.zeros(D**len(modes)))
ok.append(("phi(x)|0> มีหนึ่งอนุภาคพอดี (ตัวนับให้ 1)", sp.simplify(num*st - st) == sp.zeros(D**len(modes), 1),
           None, None))

# ── ④ [phi(x), pi(y)] = (i/L) Σ_k e^{ik(x-y)} ────────────────────
def commute_terms(T1, T2):
    """[Σ c1 o1, Σ c2 o2] ใช้แค่ [a_n, b_n] = 1 และ [b_n, a_n] = -1"""
    tot = 0
    for (c1, o1, n1_, s1), (c2, o2, n2_, s2) in itertools.product(T1, T2):
        if n1_ != n2_:
            continue
        if o1.name[0] == 'a' and o2.name[0] == 'b':
            tot += c1*c2
        elif o1.name[0] == 'b' and o2.name[0] == 'a':
            tot -= c1*c2
    return tot


def at(terms, pos):
    return [(c*sp.exp(sp.I*s*kv(n)*pos), op, n, s) for c, op, n, s in terms]


comm = commute_terms(at(phi, x), at(pi_, y))
want = sp.I/L*sum(sp.exp(sp.I*kv(n)*(x - y)) for n in ns)
chk("[phi(x), pi(y)] = (i/L) Σ_k e^(ik(x-y))  [N=2]", sp.simplify((comm - want).rewrite(sp.cos)), 0)
chk("ที่ x = y ได้ i(2N+1)/L", sp.simplify(comm.subs(y, x)), sp.I*(2*N + 1)/L)
comm_pp = commute_terms(at(phi, x), at(phi, y))
chk("[phi(x), phi(y)] = 0 (เวลาเดียวกัน)", sp.simplify(comm_pp.rewrite(sp.cos)), 0)
# ผลรวม (1/L)Σ e^{ik(x-y)} อินทิเกรตเทียบ x ได้ 1 ไม่ว่า N เท่าไร ⇒ เดลต้า
for NN in (2, 5):
    ker = sum(sp.exp(sp.I*2*sp.pi*n*(x - y)/L) for n in range(-NN, NN + 1))/L
    chk(f"∫_0^L (1/L)Σ e^(ik(x-y)) dx = 1  (N={NN})", sp.simplify(sp.integrate(ker, (x, 0, L))), 1)

# ── ⑤ กล่อง → ต่อเนื่อง: a_cont = sqrt(L) a_box ⇒ [a_cont, a_cont†] = L δ ↔ 2 pi δ(k-q) ──
# Σ_k (1/L) f(k) → ∫ dk/2pi f(k) ตรวจด้วยผลรวมของ f = e^{-k^2} บนกล่องใหญ่
Lb = 200
s_box = sum(sp.exp(-(2*sp.pi*n/Lb)**2) for n in range(-200, 201))/Lb
ok.append(("(1/L) Σ_k e^{-k²} ≈ ∫ dk/2π e^{-k²} = 1/(2√π) ที่ L = 200",
           abs(float(s_box) - float(1/(2*sp.sqrt(sp.pi)))) < 1e-9, float(s_box), float(1/(2*sp.sqrt(sp.pi)))))

# ลืม 2 ใต้ราก ⇒ ได้ 2i δ (ข้ออ้างในข้อ ④)
comm_bad = commute_terms([(c*sp.sqrt(2), o, nn, s) for c, o, nn, s in at(phi, x)],
                         [(c*sp.sqrt(2), o, nn, s) for c, o, nn, s in at(pi_, y)])
chk("ลืม 2 ใต้ราก ⇒ [phi, pi] เป็นสองเท่า", sp.simplify((comm_bad - 2*want).rewrite(sp.cos)), 0)
# ตาราง paper: (1/L) Σ (1/sqrt(2w)) (sqrt(L) a) e^{ikx} = Σ 1/sqrt(2wL) a e^{ikx}
aa, ww, kk_ = sp.symbols('aa ww kk', positive=True)
chk("(1/L)·(1/sqrt(2w))·sqrt(L) = 1/sqrt(2wL)", (1/L)*(1/sp.sqrt(2*ww))*sp.sqrt(L), 1/sp.sqrt(2*ww*L))
chk("[sqrt(L)a, sqrt(L)a†] = L ⇒ (1/L)Σ_q L δ_kq = 1 ↔ ∫dq/2π 2π δ(k-q) = 1",
    (1/L)*L, sp.integrate(sp.DiracDelta(kk_ - sp.Symbol('q', real=True))*2*sp.pi/(2*sp.pi),
                          (sp.Symbol('q', real=True), -sp.oo, sp.oo)))

# ── ตัวเลขในข้อฝึก ──────────────────────────────────────────
kbox = lambda nn, LL: 2*sp.pi*nn/LL
chk("ข้อ 1: L = 4π, n = 3 ⇒ k = 3/2", kbox(3, 4*sp.pi), sp.Rational(3, 2))
chk("ข้อ 1: L = 4π, n = -4 ⇒ k = -2", kbox(-4, 4*sp.pi), -2)
chk("ข้อ 1: L = 2π, n = 5 ⇒ k = 5", kbox(5, 2*sp.pi), 5)
chk("ข้อ 2: m=4 k=3 ⇒ w = 5 ⇒ phi_tt = -25 phi", -(sp.sqrt(9 + 16))**2, -25)
chk("ข้อ 2: phi_xx = -9 phi", -(3**2), -9)
chk("ข้อ 2: -25 + 9 + 16 = 0", -25 + 9 + 16, 0)
chk("ข้อ 3: แอมพลิจูด ∝ 1/sqrt(w) ⇒ w=5 กับ w=20 ต่างกัน 2 เท่า",
    (1/sp.sqrt(2*5*L))/(1/sp.sqrt(2*20*L)), 2)
chk("ข้อ 4: N = 3 ⇒ 7 รูปร่าง", len(range(-3, 4)), 7)
chk("ข้อ 4: L = 2 ⇒ [phi, pi] ที่จุดเดียวกัน = i*7/2", sp.I*7/2, sp.I*sp.Rational(7, 2))
chk("ข้อ 4: N = 10 ⇒ i*21/2", sp.I*len(range(-10, 11))/2, sp.I*sp.Rational(21, 2))
chk("ข้อ 5: L = 10π ⇒ ระยะห่างของ k = 2π/L = 1/5", 2*sp.pi/(10*sp.pi), sp.Rational(1, 5))
chk("ข้อ 5: ช่วง k จาก 0 ถึง 3 (ไม่รวม) มี 15 รูปร่าง", len([nn for nn in range(0, 1000) if kbox(nn, 10*sp.pi) < 3]), 15)
chk("ข้อ 5: L = 100π ⇒ 150 รูปร่าง", len([nn for nn in range(0, 5000) if kbox(nn, 100*sp.pi) < 3]), 150)
chk("ข้อ 6: L = 2π m = 4: a_3†a_0†|0> ⇒ P = 3 + 0 = 3", kbox(3, 2*sp.pi) + kbox(0, 2*sp.pi), 3)
ok.append(("ข้อ 6 ตัวควบคุม: คูณมวลได้ 12 ≠ 3", 3*4 != 3, 12, 3))
chk("ข้อ 6: E = w_3 + w_0 = 5 + 4 = 9", sp.sqrt(9 + 16) + sp.sqrt(0 + 16), 9)

bad_ = [r for r in ok if not r[1]]
for name, good, got, want_ in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want_}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
