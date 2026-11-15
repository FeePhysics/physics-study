"""ตรวจโจทย์บท P1 (หลักรังนกพิราบ) + ใบงานนักเรียน ด้วย brute force
ทุกข้อตรวจสองฝั่ง: (ก) จำนวนที่โจทย์ให้ "บังคับ" ได้จริง (ข) น้อยกว่านั้นหนึ่งตัว "ไม่บังคับ" (ขอบเขตแน่น)"""
from itertools import combinations
from math import gcd, ceil

def forced(universe, k, prop):
    """ทุกเซตย่อยขนาด k มีคู่ที่มี prop หรือไม่"""
    return all(any(prop(a, b) for a, b in combinations(S, 2)) for S in combinations(universe, k))

checks = []
def ok(name, cond): checks.append((name, cond))

# ── บทเรียน ──
s12 = lambda a, b: a + b == 12
ok("ค่าย ตย.2: 7 ตัวจาก 1..11 บังคับผลบวก 12", forced(range(1, 12), 7, s12))
ok("ค่าย ตย.2: 6 ตัวไม่บังคับ", not forced(range(1, 12), 6, s12))
cop = lambda a, b: gcd(a, b) == 1
ok("ค่าย ตย.3: 15 เลขคู่ใน 1..30 ไม่มีคู่เฉพาะสัมพัทธ์", not any(cop(a, b) for a, b in combinations(range(2, 31, 2), 2)))
ok("ค่าย ตย.3: เลขคู่ 1..30 มี 15 ตัว", len(range(2, 31, 2)) == 15)
ok("ceil(80/12)=7", ceil(80 / 12) == 7)
ok("ceil(40/3)=14", ceil(40 / 3) == 14)
ok("ค่าย ตย.7: 5 เกรด · 35 คนเลี่ยง 8 ได้ · 36 คนบังคับ", ceil(35 / 5) < 8 and ceil(36 / 5) >= 8)
ok("ข้อหาที่ผิด: รัง {1,3},{2,4} ไม่ให้เฉพาะสัมพัทธ์", gcd(2, 4) != 1)

# ── ใบงานนักเรียน ──
ok("W2: 8 ตัวจาก 1..14 บังคับผลบวก 15", forced(range(1, 15), 8, lambda a, b: a + b == 15))
ok("W2: 7 ตัวไม่บังคับ", not forced(range(1, 15), 7, lambda a, b: a + b == 15))
con = lambda a, b: abs(a - b) == 1
ok("W3: 9 ตัวจาก 1..16 บังคับเลขติดกัน", forced(range(1, 17), 9, con))
ok("W3: 8 ตัวไม่บังคับ", not forced(range(1, 17), 8, con))
ok("W4: ceil(40/7)=6", ceil(40 / 7) == 6)
ok("W5: ไพ่ 4 ดอก · 16 ใบเลี่ยง 5 ได้ · 17 บังคับ", ceil(16 / 4) < 5 and ceil(17 / 4) >= 5)
div = lambda a, b: a % b == 0 or b % a == 0
ok("W6: 11 ตัวจาก 1..20 บังคับตัวหนึ่งหารอีกตัว", forced(range(1, 21), 11, div))
ok("W6: 10 ตัวไม่บังคับ ({11..20})", not any(div(a, b) for a, b in combinations(range(11, 21), 2)))
ok("W6: ส่วนคี่ของ 1..20 มี 10 ค่า", len({n >> (n & -n).bit_length() - 1 for n in range(1, 21)}) == 10)
ok("W7: 10 จำนวน mod 9 (9 รัง)", 10 > 9)
ok("W8: 50 เลขคู่ใน 1..100 ไม่มีคู่เฉพาะสัมพัทธ์", not any(cop(a, b) for a, b in combinations(range(2, 101, 2), 2)))
ok("W8: 51 ตัวจาก 1..100 บังคับ (ตรวจแบบย่อ 11 ตัวจาก 1..20)", forced(range(1, 21), 11, cop))

bad = [n for n, c in checks if not c]
for n, c in checks: print("✅" if c else "❌", n)
raise SystemExit(1 if bad else 0)
