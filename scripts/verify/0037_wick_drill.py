"""ตรวจทุกตัวเลขของชุดฝึก Wick (0036-wick-drill.html)

สปริงสองตัว x (ω = 1) กับ y (ω = 2) บน Fock space ผลคูณเทนเซอร์ — คูณเมทริกซ์จริงทุกข้อ
"""
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(got - want) == 0, got, want))


def ladder(D):
    a = sp.zeros(D, D)
    for j in range(1, D):
        a[j - 1, j] = sp.sqrt(j)
    return a


N = 8
a = ladder(N)
I = sp.eye(N)
ax, ay = sp.kronecker_product(a, I), sp.kronecker_product(I, a)
x = (ax + ax.T)/sp.sqrt(2*1)      # ω = 1
y = (ay + ay.T)/sp.sqrt(2*2)      # ω = 2
V = sp.Matrix([1] + [0]*(N*N - 1))
ev = lambda *ops: (V.T*sp.prod(ops, start=sp.eye(N*N))*V)[0] if ops else 1


def E(*ops):
    M = sp.eye(N*N)
    for o in ops:
        M = M*o
    return (V.T*M*V)[0]


# ── ค่าของหนึ่งคู่ ───────────────────────────────────────────────
chk("⟨x²⟩ = 1/2 (ω = 1)", E(x, x), sp.Rational(1, 2))
chk("⟨y²⟩ = 1/4 (ω = 2)", E(y, y), sp.Rational(1, 4))
chk("⟨xy⟩ = 0 (สปริงคนละตัว)", E(x, y), 0)
# ตัวอย่างในเนื้อบท: x y |0⟩ ไปจบที่ |1⟩_x|1⟩_y
v = x*y*V
chk("x y|0⟩ ไม่มีส่วน |0⟩", v[0], 0)
chk("x y|0⟩ มีส่วน |1,1⟩ = (1/√2)(1/√4)", v[1*N + 1], 1/sp.sqrt(2)/sp.sqrt(4))

# ── ตารางในเนื้อบท: x⁴ ที่ ω = 1 ─────────────────────────────────
chk("x⁴ ω=1: แต่ละวิธี ½·½ = ¼ · บวกสามวิธี = ¾", E(x, x, x, x), sp.Rational(3, 4))

# ── ข้อ 1: ค่า D ไม่เท่ากัน (ตัวเลขกำหนดเอง) ────────────────────────
D12, D34, D13, D24, D14, D23 = 1, 2, 3, 0, 1, 4
chk("ข้อ 1: 1·2 + 3·0 + 1·4 = 6", D12*D34 + D13*D24 + D14*D23, 6)

# ── ข้อ 2: x x y y ──────────────────────────────────────────────
chk("ข้อ 2: ⟨x x y y⟩ = ½·¼ = 1/8", E(x, x, y, y), sp.Rational(1, 8))
# ── ข้อ 3: x y x y (ลำดับต่าง ชื่อเดิม) ──────────────────────────
chk("ข้อ 3: ⟨x y x y⟩ = 1/8 (เหลือวิธี {13, 24})", E(x, y, x, y), sp.Rational(1, 8))
chk("ข้อ 3 ตัวควบคุม: ⟨x y y x⟩ ก็ 1/8 (เหลือวิธี {14, 23})", E(x, y, y, x), sp.Rational(1, 8))
chk("⟨x x y⟩ = 0 (สามจุด)", E(x, x, y), 0)

# ── ข้อ 4: x⁶ ที่ ω = 1 ─────────────────────────────────────────
chk("ข้อ 4: ⟨x⁶⟩ = 15 · (½)³ = 15/8", E(x, x, x, x, x, x), sp.Rational(15, 8))

# ── ข้อ 5: ย้อนกลับจาก ⟨z⁴⟩ = 3/16 หา ω ───────────────────────────
w = sp.Symbol('w', positive=True)
sol = sp.solve(sp.Eq(3*(1/(2*w))**2, sp.Rational(3, 16)), w)
chk("ข้อ 5: 3/16 ÷ 3 = 1/16 ⇒ ⟨z²⟩ = 1/4 ⇒ ω = 2", sol[0], 2)

# ── ข้อ 6: ⟨y⁴⟩ − ⟨y²⟩⟨y²⟩ ที่ ω = 2 ──────────────────────────────
chk("ข้อ 6: ⟨y⁴⟩ = 3/16", E(y, y, y, y), sp.Rational(3, 16))
chk("ข้อ 6: ⟨y²⟩·⟨y²⟩ = 1/16 = ค่าของวิธี {12, 34}", E(y, y)*E(y, y), sp.Rational(1, 16))
chk("ข้อ 6: ผลต่าง = 2/16 = 1/8 = สองวิธีที่เหลือ", E(y, y, y, y) - E(y, y)**2, sp.Rational(1, 8))

bad = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad else '❌ ตก'} — {len(ok)-len(bad)}/{len(ok)}")
raise SystemExit(1 if bad else 0)
