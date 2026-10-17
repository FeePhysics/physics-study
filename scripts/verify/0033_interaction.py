"""ตรวจทุกข้ออ้างของบทที่ 21 — อันตรกิริยา lambda phi^4 และ Feynman diagram รูปแรก

1+1 มิติ · g = diag(-1, +1) · hbar = c = 1
L = ½phi_t² - ½phi_x² - ½m²phi² - (lambda/4!) phi^4
"""
import itertools
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


t, x = sp.symbols('t x', real=True)
m, lam = sp.symbols('m lambda', positive=True)
f = sp.Function('phi')(t, x)

# ── ① สมการการเคลื่อนที่มีพจน์ phi^3 ─────────────────────────────
Lag = sp.Rational(1, 2)*sp.diff(f, t)**2 - sp.Rational(1, 2)*sp.diff(f, x)**2 \
    - sp.Rational(1, 2)*m**2*f**2 - lam/sp.factorial(4)*f**4
EL = sp.euler_equations(Lag, f, [t, x])[0].lhs
want = -(sp.diff(f, t, 2) - sp.diff(f, x, 2) + m**2*f + lam/6*f**3)
chk("EL: phi_tt - phi_xx + m²phi + (λ/6)phi³ = 0", EL - want, 0)
p = sp.Symbol('p')
chk("∂/∂phi ของ -(λ/4!)phi^4 = -(λ/6)phi³", sp.diff(-lam/24*p**4, p), -lam/6*p**3)
chk("λ = 12 ⇒ สัมประสิทธิ์ λ/6 = 2", (lam/6).subs(lam, 12), 2)
chk("λ = 12 ⇒ λ/4! = 1/2", (lam/24).subs(lam, 12), sp.Rational(1, 2))

# ซ้อนคลื่นไม่ได้แล้ว: u, v แก้สมการเชิงเส้น ⇒ u + v แก้ได้ · พจน์ phi^3 ทำให้เหลือเศษ
A, B = sp.symbols('A B')
lin = lambda g: sp.diff(g, t, 2) - sp.diff(g, x, 2) + m**2*g
nonlin = lambda g: lin(g) + lam/6*g**3
u = A*sp.cos(sp.sqrt(9 + m**2)*t - 3*x)
v = B*sp.cos(sp.sqrt(4 + m**2)*t + 2*x)
chk("สมการเชิงเส้น: u + v ยังแก้ได้", sp.simplify(lin(u + v)), 0)
leftover = sp.expand(nonlin(u + v) - nonlin(u) - nonlin(v))
chk("พจน์ phi³: เศษจากการซ้อน = (λ/6)[(u+v)³ - u³ - v³]", leftover - sp.expand(lam/6*((u + v)**3 - u**3 - v**3)), 0)
ok.append(("ตัวควบคุม: เศษนั้นไม่เป็นศูนย์ (มีพจน์ u²v)", sp.simplify(leftover) != 0, None, None))

# ── ② (a + a†)^4 กางเป็นสตริงจริง 16 ตัว · นับตามจำนวน a† ─────────
a_, b_ = sp.symbols('a b', commutative=False)   # b = a†
words = list(itertools.product([a_, b_], repeat=4))
chk("(a + a†)^4 มี 2^4 = 16 สตริง", len(words), 16)
expanded = sp.expand((a_ + b_)**4)
chk("sympy กาง (a + a†)^4 ได้ 16 พจน์ (ลำดับสำคัญ ไม่รวมกัน)", len(sp.Add.make_args(expanded)), 16)
counts = [sum(1 for w in words if w.count(b_) == j) for j in range(5)]
ok.append(("นับตามจำนวน a†: 1, 4, 6, 4, 1 (ทวินาม)", counts == [1, 4, 6, 4, 1], counts, [1, 4, 6, 4, 1]))

# บน Fock space ของสปริงตัวเดียว: สตริงที่มี a† สองตัว a สองตัว พา |2> ไป |2> ได้
D = 6
A1 = sp.zeros(D, D)
for j in range(1, D):
    A1[j - 1, j] = sp.sqrt(j)
B1 = A1.H
ket = lambda n: sp.Matrix([1 if i == n else 0 for i in range(D)])
mat = {a_: A1, b_: B1}


def word_mat(w):
    M = sp.eye(D)
    for s in w:
        M = M*mat[s]
    return M


chk("สตริงที่มี a† 4 ตัว (a†a†a†a†) พา |0> ไป |4>: สัมประสิทธิ์ sqrt(24)", (ket(4).T*word_mat((b_,)*4)*ket(0))[0], sp.sqrt(24))
nonzero_22 = [w for w in words if (ket(2).T*word_mat(w)*ket(2))[0] != 0]
ok.append(("สตริงที่พา |2> → |2> ได้ คือ 6 สตริงที่มี a† สองตัวพอดี",
           len(nonzero_22) == 6 and all(w.count(b_) == 2 for w in nonzero_22), len(nonzero_22), 6))

# ── ③ 2 → 2 ด้วยสี่รูปร่างต่างกัน: <k3 k4| phi^4 |k1 k2> = 4! × c1c2c3c4 ──
Dm = 3
a1 = sp.zeros(Dm, Dm)
for j in range(1, Dm):
    a1[j - 1, j] = sp.sqrt(j)
modes = [0, 1, 2, 3]


def kron(*Ms):
    out = Ms[0]
    for M in Ms[1:]:
        out = sp.kronecker_product(out, M)
    return out


aop = {q: kron(*[a1 if r == q else sp.eye(Dm) for r in modes]) for q in modes}
c = {0: 2, 1: 3, 2: 5, 3: 7}                  # น้ำหนักต่างกันเพื่อดูว่าคูณกันครบสี่ตัว
phi0 = sum((c[q]*(aop[q] + aop[q].H) for q in modes), sp.zeros(Dm**4))
vac = sp.zeros(Dm**4, 1); vac[0] = 1
ket_in = aop[0].H*aop[1].H*vac                 # |k1 k2>
ket_out = aop[2].H*aop[3].H*vac                # |k3 k4>
phi4 = phi0**4
amp = (ket_out.T*phi4*ket_in)[0]
chk("<k3 k4| phi^4 |k1 k2> = 4! × c1c2c3c4 = 24 × 210", amp, 24*2*3*5*7)
chk("คูณ λ/4! แล้วเหลือ λ × c1c2c3c4 (4! ตัดกันพอดี)", lam/24*amp, lam*210)
chk("จำนวนวิธีจับอนุภาคภายนอก 4 ตัวกับ phi 4 ตัว = 4! = 24",
    len(list(itertools.permutations(range(4)))), 24)
# ตัวควบคุม: ถ้าใช้ phi^3 แทน ไม่มีทางพา 2 ไป 2 (จำนวนปุ่มคี่)
amp3 = (ket_out.T*phi0**3*ket_in)[0]
chk("ตัวควบคุม: phi^3 พา 2 อนุภาคไป 2 อนุภาคไม่ได้ (ได้ 0)", amp3, 0)

# โมเมนตัมอนุรักษ์: ∫_0^L e^{i(k1+k2-k3-k4)x} dx = L ถ้าผลรวมเป็น 0 ไม่งั้น 0
L = sp.Symbol('L', positive=True)
for ns, want_ in [((1, 2, 2, 1), L), ((1, 2, 3, 0), L), ((1, 2, 3, 1), 0)]:
    n1, n2, n3, n4 = ns
    val = sp.integrate(sp.exp(sp.I*2*sp.pi*(n1 + n2 - n3 - n4)*x/L), (x, 0, L))
    chk(f"∫ e^(i(k1+k2-k3-k4)x) dx ที่ n = {ns}", sp.simplify(val), want_)

# 1+1 มิติ มวลเท่ากัน: อนุรักษ์ทั้ง P และ E ⇒ ขาออกเป็นชุดเดียวกับขาเข้า
k3, k4 = sp.symbols('k3 k4', real=True)
mm = 12
sol = sp.solve([k3 + k4 - (5 + (-5)), sp.sqrt(k3**2 + mm**2) + sp.sqrt(k4**2 + mm**2) - 26], [k3, k4], dict=True)
ok.append(("1+1 มิติ: เข้า (5, -5) ออกได้แค่ {5, -5}", sorted({tuple(sorted((s[k3], s[k4]))) for s in sol}) == [(-5, 5)],
           sol, "{(5,-5)}"))
sol2 = sp.solve([k3 + k4 - (5 + (-9)), sp.sqrt(k3**2 + mm**2) + sp.sqrt(k4**2 + mm**2) - 28], [k3, k4], dict=True)
ok.append(("1+1 มิติ: เข้า (5, -9) ออกได้แค่ {5, -9}", sorted({tuple(sorted((s[k3], s[k4]))) for s in sol2}) == [(-9, 5)],
           sol2, "{(-9,5)}"))

# ── ④ นับเส้นใน diagram ────────────────────────────────────────
def internal_lines(V, E):
    """จุดยอดละ 4 ขา · ขาภายนอก E · ที่เหลือจับคู่เป็นเส้นภายใน"""
    return sp.Rational(4*V - E, 2)


chk("1 จุดยอด 4 ขาภายนอก ⇒ เส้นภายใน 0", internal_lines(1, 4), 0)
chk("2 จุดยอด 4 ขาภายนอก ⇒ ปลายขา 8 · เส้นภายใน 2", internal_lines(2, 4), 2)
chk("ห่วง (loop) = เส้นภายใน - จุดยอด + 1: 2 - 2 + 1 = 1", internal_lines(2, 4) - 2 + 1, 1)
chk("3 จุดยอด 4 ขาภายนอก ⇒ เส้นภายใน 4", internal_lines(3, 4), 4)

# ── ตัวเลขในข้อฝึก ──────────────────────────────────────────
chk("ข้อ 1: λ = 12 ⇒ λ/4! = 1/2", sp.Rational(12, 24), sp.Rational(1, 2))
chk("ข้อ 1: 4 × (1/2) = 2 (อนุพันธ์ของ phi^4 ได้ 4phi^3)", 4*sp.Rational(1, 2), 2)
ok.append(("ข้อ 2: 16 · 6 · 1", (len(words), counts[2], counts[4]) == (16, 6, 1), None, None))
chk("ข้อ 3: 4! = 24", sp.factorial(4), 24)
chk("ข้อ 3: λ = 6 ⇒ (λ/4!)·4! = 6", sp.Rational(6, 24)*24, 6)
chk("ข้อ 4: P เข้า = 5 + (-5) = 0", 5 + (-5), 0)
chk("ข้อ 4: k4 = 5 ⇒ k3 = -5", 0 - 5, -5)
chk("ข้อ 4: E = 13 + 13 = 26", 2*sp.sqrt(25 + 144), 26)
ok.append(("ข้อ 5: ปลายขา 8 · ภายนอก 4 · ภายใน 2", (4*2, 4, internal_lines(2, 4)) == (8, 4, 2), None, None))
chk("ข้อ 6: 3 จุดยอด ⇒ λ^3", 3, 3)
chk("ข้อ 6: λ = 1/10 ⇒ λ²/λ = 1/10", sp.Rational(1, 10)**2/sp.Rational(1, 10), sp.Rational(1, 10))
chk("ข้อ 6: λ³ = 1/1000", sp.Rational(1, 10)**3, sp.Rational(1, 1000))

bad_ = [r_ for r_ in ok if not r_[1]]
for name, good, got, want_ in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want_}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
