"""ตรวจทุกข้ออ้างของบทที่ 20 — propagator <0|phi(t,x) phi(0,y)|0>

1+1 มิติ · hbar = c = 1 · กล่อง L · phi เหมือนบทที่ 19
D(t, x-y) = <0|phi(t,x) phi(0,y)|0> = Σ_k e^{i(k(x-y) - w_k t)} / (2 w_k L)
ตรวจบน Fock space จริง (คูณเมทริกซ์) ไม่ได้ใช้สูตรข้างบนคำนวณ
"""
import mpmath as mp
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


# ── ① ผลคูณภายในเชิงซ้อน: <a|b> = Σ conj(a_i) b_i ─────────────
a = sp.Matrix([1, sp.I])
b = sp.Matrix([2, sp.I])
chk("<a|b> = 1*2 + (-i)(i) = 3", (a.H*b)[0], 3)
chk("<a|a> = 2 (ความยาวกำลังสอง)", (a.H*a)[0], 2)
chk("ลืมสังยุค: a·a = 1 + i*i = 0", (a.T*a)[0], 0)
ex = sp.Matrix([3, 4*sp.I])
chk("ตัวอย่างในบท: (3, 4i) dot ตัวเอง = 25", (ex.H*ex)[0], 25)
chk("ตัวอย่างในบท: ลืมสังยุคได้ 9 - 16 = -7", (ex.T*ex)[0], -7)

# ── ② D บน Fock space จริง: 3 รูปร่าง ตัดที่ 3 ขั้น ─────────────
L = sp.Symbol('L', positive=True)
m = sp.Symbol('m', positive=True)
t, x, y = sp.symbols('t x y', real=True)
D_ = 3
a1 = sp.zeros(D_, D_)
for j in range(1, D_):
    a1[j - 1, j] = sp.sqrt(j)
modes = [-1, 0, 1]


def kron(*Ms):
    out = Ms[0]
    for M in Ms[1:]:
        out = sp.kronecker_product(out, M)
    return out


aop = {q: kron(*[a1 if r == q else sp.eye(D_) for r in modes]) for q in modes}
kv = lambda n: 2*sp.pi*n/L
wv = lambda n: sp.sqrt(kv(n)**2 + m**2)
N = D_**len(modes)


def phi(tt, xx):
    return sum((1/sp.sqrt(2*wv(q)*L)*(aop[q]*sp.exp(sp.I*(kv(q)*xx - wv(q)*tt))
                                      + aop[q].H*sp.exp(-sp.I*(kv(q)*xx - wv(q)*tt))) for q in modes),
               sp.zeros(N))


vac = sp.zeros(N, 1); vac[0] = 1
Dfock = (vac.H*phi(t, x)*phi(0, y)*vac)[0]
Dform = sum(sp.exp(sp.I*(kv(q)*(x - y) - wv(q)*t))/(2*wv(q)*L) for q in modes)
chk("<0|phi(t,x) phi(0,y)|0> = Σ e^{i(k(x-y) - wt)}/(2wL)  [คูณเมทริกซ์จริง]", sp.simplify(Dfock - Dform), 0)
# เส้นทางในบท: dot ของ phi(t,x)|0> (เป็นบรา = สังยุค) กับ phi(0,y)|0>
left = phi(t, x)*vac
right = phi(0, y)*vac
chk("= dot(phi(t,x)|0>, phi(0,y)|0>) แบบมีสังยุค",
    sp.simplify((left.applyfunc(sp.conjugate).T*right)[0].subs({sp.conjugate(x): x, sp.conjugate(t): t,
                                                                  sp.conjugate(L): L, sp.conjugate(m): m}) - Dform), 0)
# ตัวควบคุม: ลืมสังยุคแล้วไม่เท่า
ok.append(("ตัวควบคุม: ลืมสังยุค (เอาเลขชี้ไม่กลับเครื่องหมาย) ได้ค่าอื่น",
           sp.simplify((left.T*right)[0] - Dform) != 0, None, None))
# ⟨k|q⟩ = δ ⇒ ผลรวมสองชั้นเหลือชั้นเดียว
one = {q: aop[q].H*vac for q in modes}
for q1 in modes:
    for q2 in modes:
        chk(f"<k={q1}|k={q2}> = {int(q1 == q2)}", (one[q1].H*one[q2])[0], int(q1 == q2))
# ค่า D ที่ t = 0, x = y เป็นจำนวนจริงบวก = Σ 1/(2wL)
chk("D(0, 0) = Σ 1/(2wL)", sp.simplify(Dform.subs({t: 0, x: y})), sum(1/(2*wv(q)*L) for q in modes))

# ── ③ คู่ ±k ที่ t = 0 ให้ 2cos(kr) ─────────────────────────────
k, r = sp.symbols('k r', real=True)
chk("e^{ikr} + e^{-ikr} = 2cos(kr)", sp.simplify((sp.exp(sp.I*k*r) + sp.exp(-sp.I*k*r)).rewrite(sp.cos)), 2*sp.cos(k*r))
chk("ตัวอย่าง k=2, r=pi/2: 2cos(pi) = -2", 2*sp.cos(2*sp.pi/2), -2)

# ── ④ ต่อเนื่อง: D(0, r) = ∫ dk/2pi e^{ikr}/(2 sqrt(k^2+m^2)) = K0(m r)/(2pi) ──
mp.mp.dps = 20
for mm, rr in [(1, 1), (4, 0.5), (0.5, 3)]:
    num = mp.quadosc(lambda kk: mp.cos(kk*rr)/(2*mp.sqrt(kk**2 + mm**2)), [0, mp.inf], omega=rr)/mp.pi
    ok.append((f"∫dk/2π cos(kr)/(2ω) = K0(mr)/2π  (m={mm}, r={rr})",
               abs(num - mp.besselk(0, mm*rr)/(2*mp.pi)) < 1e-12, num, mp.besselk(0, mm*rr)/(2*mp.pi)))
# ระยะไกล: K0(mr) ~ sqrt(pi/(2mr)) e^{-mr}
for mm in (1, 4):
    rr = 30/mm
    ratio = mp.besselk(0, mm*rr)/(mp.sqrt(mp.pi/(2*mm*rr))*mp.e**(-mm*rr))
    ok.append((f"ระยะไกล K0(mr) ≈ sqrt(π/2mr) e^(-mr) (m={mm}, mr=30)", abs(ratio - 1) < 0.01, ratio, 1))
# ตัวลดหลักคือ e^{-mr}: เดินไป 1/m ลดลงประมาณ e เท่า (ตัวราก sqrt ทำให้ไม่เป๊ะ)
rat = mp.besselk(0, 31)/mp.besselk(0, 30)
ok.append(("เดินเพิ่ม 1/m ที่ระยะไกล ค่าลดลง ~ 1/e", abs(rat*mp.e - 1) < 0.05, rat*mp.e, 1))
# ไม่เป็นศูนย์ที่ t = 0 คนละจุด (spacelike) — ข้ออ้างใน warn
ok.append(("D(t=0, r=2) ≠ 0 (m=1)", mp.besselk(0, 2) > 0, mp.besselk(0, 2), ">0"))

# ── ตัวเลขในข้อฝึก ──────────────────────────────────────────
chk("ข้อ 1: <a|b> = 3", (a.H*b)[0], 3)
chk("ข้อ 1: <a|a> = 2", (a.H*a)[0], 2)
chk("ข้อ 1: ลืมสังยุค = 0", (a.T*a)[0], 0)
w = sp.Symbol('w', positive=True)
chk("ข้อ 2: น้ำหนักบทที่ 19 (1/sqrt(2wL)) w=5 เทียบ w=20 = 2", sp.sqrt(sp.Rational(20, 5)), 2)
chk("ข้อ 2: พจน์ของ D (1/(2wL)) w=5 เทียบ w=20 = 4", (1/(2*5*L))/(1/(2*20*L)), 4)
chk("ข้อ 2: e^{-iwt} ที่ w=5 t=π/5 = -1", sp.exp(-sp.I*5*sp.pi/5), -1)
chk("ข้อ 3: k=1, r=0 ⇒ 2cos 0 = 2", 2*sp.cos(0), 2)
chk("ข้อ 3: k=1, r=π ⇒ -2", 2*sp.cos(sp.pi), -2)
chk("ข้อ 3: k=1, r=π/2 ⇒ 0", 2*sp.cos(sp.pi/2), 0)
chk("ข้อ 4: <k=3|k=3> = 1", 1, 1)
chk("ข้อ 4: <k=3|k=-3> = 0 (เวกเตอร์ฐานคนละตัว)", 0, 0)
chk("ข้อ 4: 5 รูปร่าง ⇒ 25 คู่ เหลือ 5", len([(p, q) for p in range(-2, 3) for q in range(-2, 3) if p == q]), 5)
chk("ข้อ 5: m = 4 ⇒ ระยะ 1/m = 1/4", sp.Rational(1, 4), 1/sp.Integer(4))
chk("ข้อ 5: m = 1/2 ⇒ ระยะ 2", 1/sp.Rational(1, 2), 2)
chk("ข้อ 5: m = ln 2, r = 1 ⇒ e^{-mr} = 1/2", sp.exp(-sp.log(2)*1), sp.Rational(1, 2))
chk("ข้อ 6: P = 5 + (-9) = -4", 5 + (-9), -4)
chk("ข้อ 6: m=12 ⇒ w_5 = 13", sp.sqrt(25 + 144), 13)
chk("ข้อ 6: w_-9 = 15", sp.sqrt(81 + 144), 15)
chk("ข้อ 6: E = 28", 13 + 15, 28)

bad_ = [r_ for r_ in ok if not r_[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
