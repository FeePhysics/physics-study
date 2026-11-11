"""ตรวจทุกตัวเลขของบทที่ 25 — ใส่จุดยอด φ_v⁴ แล้วจับคู่ (Wick) · แยกตามรูป

ไล่ทุกวิธีจับคู่ของตำแหน่งจริง แล้วจัดกลุ่มตามรูปกราฟ (ใครต่อกับใคร) ไม่ใช่นับด้วยสูตร
"""
import itertools
from fractions import Fraction
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, got == want, got, want))


def pairings(items):
    if not items:
        yield []
        return
    first, rest = items[0], items[1:]
    for i, p in enumerate(rest):
        for tail in pairings(rest[:i] + rest[i + 1:]):
            yield [(first, p)] + tail


def shape(pairing, label):
    """รูปของการจับคู่ = multiset ของเส้น (ชื่อปลายสองข้าง) — ขาของจุดยอดทุกขาชื่อ 'v'"""
    return tuple(sorted(tuple(sorted((label[a], label[b]))) for a, b in pairing))


def classify(ext):
    label = {i: e for i, e in enumerate(ext)}
    n = len(ext)
    for k in range(4):
        label[n + k] = 'v'
    counts = {}
    for p in pairings(list(range(n + 4))):
        s = shape(p, label)
        counts[s] = counts.get(s, 0) + 1
    return counts


# ── สองจุด ⟨φ1 φ2 φv⁴⟩ (ตัวอย่างในเนื้อบท) ───────────────────────────
c2 = classify(['1', '2'])
chk("สองจุด: ทั้งหมด 5!! = 15 วิธี", sum(c2.values()), 15)
conn2 = c2[(('1', 'v'), ('2', 'v'), ('v', 'v'))]
disc2 = c2[(('1', '2'), ('v', 'v'), ('v', 'v'))]
chk("สองจุด: แบบต่อถึงกัน (1–v–2 + ห่วงที่ v) = 12", conn2, 12)
chk("สองจุด: แบบแยก (1–2 + ห่วงคู่ที่ v) = 3", disc2, 3)
chk("สองจุด: มีแค่สองรูป", len(c2), 2)
chk("ตัวคูณของรูปต่อถึงกัน 12/4! = 1/2", Fraction(conn2, 24), Fraction(1, 2))
chk("ตัวคูณของรูปแยก 3/4! = 1/8", Fraction(disc2, 24), Fraction(1, 8))
chk("วิธีนับในบท: 4 × 3 × 1 = 12", 4*3*1, 12)

# ── สี่จุด ⟨φ1 φ2 φ3 φ4 φv⁴⟩ (ข้อฝึก) ─────────────────────────────────
c4 = classify(['1', '2', '3', '4'])
chk("สี่จุด: ทั้งหมด 7!! = 105", sum(c4.values()), 105)
chk("7 × 5 × 3 × 1 = 105", 7*5*3*1, 105)
X = c4[(('1', 'v'), ('2', 'v'), ('3', 'v'), ('4', 'v'))]
chk("รูป X (ทุกขาภายนอกต่อที่ v) = 4! = 24", X, 24)
chk("X: λ/4! × 24 = λ (ตัวคูณ 1)", Fraction(X, 24), 1)
disc = sum(n for s, n in c4.items() if s.count(('v', 'v')) == 2)
chk("แบบแยกทั้งหมด (ภายนอกจับกันเอง × ห่วงคู่ที่ v) = 3 × 3 = 9", disc, 9)
mid = sum(n for s, n in c4.items() if s.count(('v', 'v')) == 1)
chk("แบบกลาง (สองขาเข้า v · สองขาจับกัน · ห่วงหนึ่งที่ v) = 72", mid, 72)
chk("C(4,2) × 4 × 3 = 72", sp.binomial(4, 2)*4*3, 72)
chk("24 + 72 + 9 = 105", X + mid + disc, 105)
# แบบกลางแยกตามว่าคู่ไหนเข้า v: 6 รูป รูปละ 12
mids = [n for s, n in c4.items() if s.count(('v', 'v')) == 1]
chk("แบบกลางมี 6 รูป รูปละ 12", (len(mids), set(mids)), (6, {12}))
chk("แบบแยก: 3 รูป รูปละ 3", sorted(n for s, n in c4.items() if s.count(('v', 'v')) == 2), [3, 3, 3])

# ── ตรวจไขว้ด้วยเมทริกซ์จริง: ถ้าทุกเส้นมีค่าเท่ากัน (สปริงเดียว) ผลรวม = ⟨x⁸⟩ ──
D = 12
a = sp.zeros(D, D)
for j in range(1, D):
    a[j - 1, j] = sp.sqrt(j)
X_ = a + a.T
vac = sp.Matrix([1] + [0]*(D - 1))
chk("ตรวจไขว้: ⟨0|(a+a†)⁸|0⟩ = 105", (vac.T*X_**8*vac)[0], 105)
chk("ตรวจไขว้: ⟨0|(a+a†)⁶|0⟩ = 15", (vac.T*X_**6*vac)[0], 15)

bad = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad else '❌ ตก'} — {len(ok)-len(bad)}/{len(ok)}")
raise SystemExit(1 if bad else 0)
