"""ตรวจโจทย์บท P3 (จับคู่แล้วเหลือเศษ · parity) + ใบงาน ด้วย brute force"""
from itertools import combinations, permutations
from math import prod

checks = []
def ok(name, c): checks.append((name, c))

# ── ข้อสอบข้อ 3 (ย่อส่วน): วางหมาก k ตัวบนกระดาน n×n แบบสมมาตรตามเส้นทแยง (ช่องไม่ซ้ำ)
def sym_sets(n, k):
    cells = [(i, j) for i in range(n) for j in range(n)]
    for S in combinations(cells, k):
        s = set(S)
        if all((j, i) in s for i, j in s):
            yield s
for n, k in ((3, 3), (4, 3), (5, 5), (4, 5)):
    ok(f"กระดาน {n}×{n} หมาก {k} ตัว (คี่) สมมาตร ⇒ มีตัวบนเส้นทแยงเสมอ",
       all(any(i == j for i, j in s) for s in sym_sets(n, k)))
ok("กระดาน 4×4 หมาก 4 ตัว (คู่) สมมาตร ⇒ ไม่จำเป็นต้องมีบนเส้นทแยง",
   any(not any(i == j for i, j in s) for s in sym_sets(4, 4)))
ok("ตัวอย่างค้านข้อ 6 ในบท: (i,i+1) กับ (9,1) ไม่มีตัวบนเส้นทแยง", not any(i == j for i, j in [(i, i + 1) for i in range(1, 9)] + [(9, 1)]))
ok("ช่อง (3,7) สะท้อนไป (7,3)", (7, 3) == (3, 7)[::-1])

# ── ข้อสอบข้อ 6 (ย่อส่วน): A ขนาดคี่ ⇒ ผลคูณ (a_i - b_i) เป็นคู่เสมอ · ขนาดคู่มีตัวอย่างค้าน
for A in ([3, 8, -1, 10, 7], [2, 5, 0, 9, 4, 11, 6]):
    ok(f"|A| = {len(A)} (คี่): ทุกการเรียงสลับให้ผลคูณคู่",
       all(prod(a - b for a, b in zip(A, p)) % 2 == 0 for p in permutations(A)))
for A in ([1, 2], [3, 8, -1, 10]):
    ok(f"|A| = {len(A)} (คู่): มีการเรียงสลับที่ให้ผลคูณคี่",
       any(prod(a - b for a, b in zip(A, p)) % 2 for p in permutations(A)))
ok("ตัวอย่างค้าน a=(1,2) b=(2,1): ผลคูณ = -1", (1 - 2) * (2 - 1) == -1)
ok("ใบงานข้อ 6: a=(1,2,3,4) b=(2,1,4,3) ได้ผลคูณ 1 (คี่)", prod(x - y for x, y in zip((1,2,3,4),(2,1,4,3))) == 1)
ok("ผลบวกของ 13 จำนวนคี่ เป็นคี่ (ไล่ 1..13)", sum(range(1, 27, 2)) % 2 == 1)

# ── จับมือ: k-regular graph บน v จุด มีได้ก็ต่อเมื่อ v·k คู่ (ตรวจ v = 5, k = 3 · v = 4, k = 3 · v = 6, k = 3)
def regular_exists(v, k):
    E = list(combinations(range(v), 2))
    for mask in range(1 << len(E)):
        deg = [0] * v
        for t, (a, b) in enumerate(E):
            if mask >> t & 1: deg[a] += 1; deg[b] += 1
        if all(d == k for d in deg): return True
    return False
ok("5 คน จับมือคนละ 3 คน เป็นไปไม่ได้", not regular_exists(5, 3))
ok("4 คน จับมือคนละ 3 คน เป็นไปได้", regular_exists(4, 3))
ok("6 คน จับมือคนละ 3 คน เป็นไปได้", regular_exists(6, 3))
ok("15 × 3 = 45 เป็นคี่ · 7 × 3 = 21 เป็นคี่", 45 % 2 == 1 and 21 % 2 == 1)

# ── โดมิโน: กระดาน 3×3 ปูด้วยโดมิโน 1×2 ไม่ได้ (ตรวจด้วยการค้นหา) · 5×5 มี 25 ช่อง (คี่)
def tileable(r, c):
    cells = frozenset((i, j) for i in range(r) for j in range(c))
    def go(left):
        if not left: return True
        i, j = min(left)
        for d in ((0, 1), (1, 0)):
            q = (i + d[0], j + d[1])
            if q in left and go(left - {(i, j), q}): return True
        return False
    return go(cells)
ok("3×3 ปูโดมิโนไม่ได้", not tileable(3, 3))
ok("2×3 ปูได้ (ตัวตรวจทำงาน)", tileable(2, 3))

# ── ทางเลือก (รังนก) ของข้อ 6: ชนิดที่มีอย่างน้อย 7 ตัว · 7 + 7 = 14 > 13
ok("13 = คี่ + คู่ ⇒ ชนิดหนึ่งมี ≥ 7", all(max(k, 13 - k) >= 7 for k in range(14)))

bad = [n for n, c in checks if not c]
for n, c in checks: print("✅" if c else "❌", n)
raise SystemExit(1 if bad else 0)
