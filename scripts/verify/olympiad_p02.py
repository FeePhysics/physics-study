"""ตรวจโจทย์บท P2 (รังจากเศษที่จับคู่ r กับ -r) + ใบงาน ด้วย brute force บนเศษ
ทุกข้อตรวจสองฝั่ง: จำนวนที่โจทย์ให้บังคับได้จริง และน้อยกว่าหนึ่งตัวไม่บังคับ"""
from itertools import combinations, combinations_with_replacement as cwr

def forced_mod(m, k, prop):
    """จำนวนเต็ม k ตัวใด ๆ มีคู่ที่ prop(a,b,m) จริงเสมอ — พอตรวจบนเศษ mod m เพราะ prop ขึ้นกับเศษเท่านั้น"""
    return all(any(prop(a, b, m) for a, b in combinations(R, 2)) for R in cwr(range(m), k))

sum_or_diff = lambda a, b, m: (a - b) % m == 0 or (a + b) % m == 0
sq_diff = lambda a, b, m: (a * a - b * b) % m == 0
nests = lambda m: m // 2 + 1          # {0}, {±1}, ..., (และ {m/2} ถ้า m คู่)

checks = []
def ok(name, c): checks.append((name, c))

# จำนวนรังจากการนับจริง (ไล่จับคู่ r กับ m-r) ตรงกับสูตร
for m in (7, 10, 12, 14, 100, 1000):
    groups = {frozenset({r % m, (-r) % m}) for r in range(m)}
    ok(f"จำนวนรัง mod {m} = {len(groups)} = m//2+1", len(groups) == nests(m))

for m, k in ((10, 7), (7, 5), (12, 8)):
    ok(f"ผลบวกหรือผลต่าง mod {m}: {k} ตัวบังคับ", forced_mod(m, k, sum_or_diff))
    ok(f"ผลบวกหรือผลต่าง mod {m}: {k-1} ตัวไม่บังคับ", not forced_mod(m, k - 1, sum_or_diff))
ok("ผลต่างอย่างเดียว mod 12: 13 ตัวบังคับ / 12 ไม่", forced_mod(12, 13, lambda a,b,m:(a-b)%m==0) and not forced_mod(12, 12, lambda a,b,m:(a-b)%m==0))

# ข้อสอบข้อ 2: 0..500 ไม่มีคู่ (ขอบเขต 502 แน่น)
ok("ข้อสอบ 2: {0..500} ไม่มีคู่", not any(sum_or_diff(a, b, 1000) for a, b in combinations(range(501), 2)))
ok("ตย. m=10: {0..5} ไม่มีคู่", not any(sum_or_diff(a, b, 10) for a, b in combinations(range(6), 2)))
ok("1 + 999 ≡ 0 mod 1000", (1 + 999) % 1000 == 0)

# ข้อสอบข้อ 4 แบบย่อส่วน: รังเดียวกันแล้ว a^2-b^2 หาร m ลงตัวเสมอ
for m in (10, 14, 100):
    ok(f"a≡±b mod {m} ⇒ a²−b² ≡ 0", all(sq_diff(a, b, m) for a in range(m) for b in (a, (-a) % m)))
ok("ใบงาน: 9 ตัว ⇒ a²−b² หาร 14 ลงตัว", forced_mod(14, 9, sq_diff))

# ค่ายแบบฝึก 18: S = {1,4,...,100} · 19 ตัวบังคับผลบวก 104 · 18 ไม่
S = list(range(1, 101, 3))
pairs = [(a, 104 - a) for a in S if a < 52 and 104 - a in S]
ok("ค่าย 18: |S| = 34", len(S) == 34)
ok("ค่าย 18: คู่ที่บวกได้ 104 มี 16 คู่", len(pairs) == 16)
ok("ค่าย 18: 16 คู่ + {1} + {52} ครอบ S ครบ (18 รัง)", sorted([x for p in pairs for x in p] + [1, 52]) == S)
bad18 = [1, 52] + [a for a, _ in pairs]
ok("ค่าย 18: เซต 18 ตัว (1, 52, ตัวเล็กของทุกคู่) ไม่มีคู่", len(bad18) == 18 and not any(a + b == 104 for a, b in combinations(bad18, 2)))

bad = [n for n, c in checks if not c]
for n, c in checks: print("✅" if c else "❌", n)
raise SystemExit(1 if bad else 0)
