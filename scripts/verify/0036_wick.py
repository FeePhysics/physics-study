"""ตรวจทุกข้ออ้างของบทที่ 24 — Wick's theorem ด้วยการจับคู่

① ⟨0|(a+a†)^n|0⟩ จากการคูณเมทริกซ์ปุ่มจริง = จำนวนวิธีจับคู่ (n-1)!!
② สตริงที่รอดใน n = 4 และค่าของแต่ละสตริง
③ x = (a+a†)/sqrt(2ω) ⇒ ⟨x^n⟩ = (n-1)!! (1/2ω)^{n/2} · เทียบกับ Gaussian ของ ψ0
④ สนามสองรูปร่าง: ⟨φ1φ2φ3φ4⟩ = D12 D34 + D13 D24 + D14 D23 (คงลำดับในคู่)
"""
import itertools
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


def ladder(D):
    a = sp.zeros(D, D)
    for j in range(1, D):
        a[j - 1, j] = sp.sqrt(j)
    return a


def pairings(items):
    """ทุกวิธีแบ่ง items เป็นคู่ (คงลำดับเดิมภายในคู่)"""
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for i, p in enumerate(rest):
        for tail in pairings(rest[:i] + rest[i + 1:]):
            yield [(first, p)] + tail


def dfact(n):
    return sp.Integer(1) if n <= 0 else n*dfact(n - 2)


# ── ① สปริงตัวเดียว ──────────────────────────────────────────────
D = 12
a = ladder(D); ad = a.T
vac = sp.Matrix([1] + [0]*(D - 1))
X = a + ad
for n in range(0, 9):
    got = (vac.T*X**n*vac)[0]
    want = dfact(n - 1) if n % 2 == 0 else 0
    chk(f"⟨0|(a+a†)^{n}|0⟩ = {want}", got, want)
    chk(f"  = จำนวนวิธีจับคู่ {n} สิ่ง", len(list(pairings(list(range(n))))) if n % 2 == 0 else 0, want)

# ตัวอย่างในบท n = 2: สี่สตริง มีแค่ a a† ที่รอด
ops = {'a': a, 'A': ad}
for s, want in [('aA', 1), ('Aa', 0), ('aa', 0), ('AA', 0)]:
    M = sp.eye(D)
    for c in s:
        M = M*ops[c]
    chk(f"n=2 สตริง {s.replace('A', 'a†')} บน |0⟩ เหลือ |0⟩ = {want}", (vac.T*M*vac)[0], want)

# ── ② n = 4: 16 สตริง รอดสอง ─────────────────────────────────────
alive = {}
for s in itertools.product('aA', repeat=4):
    M = sp.eye(D)
    for c in s:
        M = M*ops[c]
    v = (vac.T*M*vac)[0]
    if v != 0:
        alive[''.join(s)] = v
chk("n=4 สตริงที่รอดมี 2 สตริง", len(alive), 2)
chk("a a a† a† |0⟩ → 2|0⟩", alive.get('aaAA', 0), 2)
chk("a a† a a† |0⟩ → 1|0⟩", alive.get('aAaA', 0), 1)
chk("ผลรวม 2 + 1 = 3 = วิธีจับคู่ 4 สิ่ง", sum(alive.values()), 3)
# ทีละปุ่มของ a a a† a† (อ่านจากขวา): √1 · √2 · √2 · √1
chk("ตัวคูณทีละปุ่ม √1·√2·√2·√1 = 2", sp.sqrt(1)*sp.sqrt(2)*sp.sqrt(2)*sp.sqrt(1), 2)

# การนับคู่แบบ "คู่ของตัวแรกมี n-1 ทาง"
chk("6 สิ่ง: 5 × 3 × 1 = 15", 5*3*1, len(list(pairings(list(range(6))))))
chk("4 สิ่ง: 3 × 1 = 3", 3, len(list(pairings(list(range(4))))))
chk("8 สิ่ง: 7!! = 105", dfact(7), len(list(pairings(list(range(8))))))

# ── ③ x = (a+a†)/sqrt(2ω) ───────────────────────────────────────
w = sp.Symbol('omega', positive=True)
x = X/sp.sqrt(2*w)
chk("⟨x²⟩ = 1/(2ω)", (vac.T*x**2*vac)[0], 1/(2*w))
chk("⟨x⁴⟩ = 3/(2ω)²", (vac.T*x**4*vac)[0], 3/(2*w)**2)
chk("⟨x⁶⟩ = 15/(2ω)³", (vac.T*x**6*vac)[0], 15/(2*w)**3)
chk("ω = 1: ⟨x⁴⟩ = 3/4 (ตัวอย่างในบท)", ((vac.T*x**4*vac)[0]).subs(w, 1), sp.Rational(3, 4))
chk("ω = 2: ⟨x²⟩ = 1/4", ((vac.T*x**2*vac)[0]).subs(w, 2), sp.Rational(1, 4))
chk("ω = 2: ⟨x⁴⟩ = 3/16", ((vac.T*x**4*vac)[0]).subs(w, 2), sp.Rational(3, 16))
chk("ตัวควบคุม: ลืมยกกำลังสอง 3 × 1/4 = 3/4 ≠ 3/16", 3*sp.Rational(1, 4) - sp.Rational(3, 16) != 0, True)
chk("⟨x³⟩ = 0", (vac.T*x**3*vac)[0], 0)

# ภาพคุ้นเคย: |ψ0|² เป็น Gaussian ความแปรปรวน 1/(2ω) (Griffiths §2.3) ⇒ โมเมนต์เดียวกัน
q = sp.Symbol('q', real=True)
rho = sp.sqrt(w/sp.pi)*sp.exp(-w*q**2)
chk("∫|ψ0|² = 1", sp.integrate(rho, (q, -sp.oo, sp.oo)), 1)
chk("Gaussian: E[q²] = 1/(2ω)", sp.integrate(q**2*rho, (q, -sp.oo, sp.oo)), 1/(2*w))
chk("Gaussian: E[q⁴] = 3/(2ω)²", sp.integrate(q**4*rho, (q, -sp.oo, sp.oo)), 3/(2*w)**2)

# ── ④ สนาม: สองรูปร่าง, φ_i = Σ_k (c_ik a_k + conj(c_ik) a_k†) ───────
N = 5  # ระดับต่อรูปร่าง พอสำหรับปุ่มสี่ตัว
a1 = sp.kronecker_product(ladder(N), sp.eye(N))
a2 = sp.kronecker_product(sp.eye(N), ladder(N))
V = sp.Matrix([1] + [0]*(N*N - 1))
coef = [(1 + 2*sp.I, sp.Rational(1, 2) - sp.I), (3 - sp.I, 2 + sp.I),
        (-1 + sp.I, 1), (2, -3*sp.I)]
phi = [c1*a1 + sp.conjugate(c1)*a1.T + c2*a2 + sp.conjugate(c2)*a2.T for c1, c2 in coef]
Dm = [[(V.T*phi[i]*phi[j]*V)[0] for j in range(4)] for i in range(4)]
four = (V.T*phi[0]*phi[1]*phi[2]*phi[3]*V)[0]
wick = Dm[0][1]*Dm[2][3] + Dm[0][2]*Dm[1][3] + Dm[0][3]*Dm[1][2]
chk("⟨φ1φ2φ3φ4⟩ = D12D34 + D13D24 + D14D23 (คูณเมทริกซ์จริง)", four, wick)
wrong = Dm[0][1]*Dm[2][3] + Dm[0][2]*Dm[1][3] + Dm[3][0]*Dm[2][1]
chk("ตัวควบคุม: กลับลำดับในคู่ (D41 D32) ให้ค่าผิด", sp.simplify(four - wrong) != 0, True)
chk("ตัวควบคุม: D21 ≠ D12 (ลำดับสำคัญ)", sp.simplify(Dm[1][0] - Dm[0][1]) != 0, True)
chk("D21 = สังยุคของ D12", sp.simplify(Dm[1][0] - sp.conjugate(Dm[0][1])), 0)
chk("⟨φ1φ2φ3⟩ = 0", (V.T*phi[0]*phi[1]*phi[2]*V)[0], 0)

# ตัวเลขในข้อฝึก
chk("ข้อฝึก: 2·3 + 1·5 + 4·½ = 13", 2*3 + 1*5 + 4*sp.Rational(1, 2), 13)
chk("หกจุด: 15 พจน์ แต่ละพจน์มี D สามตัว", len(list(pairings(list(range(6))))), 15)

# ฉบับแก้ (ผล 4/6)
chk("ω=1 ต่อวิธี ½×½ = ¼ · รวม 3×¼ = ¾", 3*sp.Rational(1, 4), sp.Rational(3, 4))
chk("สามจุด: คู่ที่หยิบได้ C(3,2) = 3", sp.binomial(3, 2), 3)
chk("สามจุด: วิธีจับคู่ทุกจุด = 0", sum(1 for p in pairings([0, 1, 2]) if len(p)*2 == 3), 0)
chk("หกจุด: ทุกวิธีมี 3 คู่", min(len(p) for p in pairings(list(range(6)))) == max(len(p) for p in pairings(list(range(6)))) == 3, True)

bad = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad else '❌ ตก'} — {len(ok)-len(bad)}/{len(ok)}")
raise SystemExit(1 if bad else 0)
