"""ตรวจโจทย์บท P4 (หลักสุดขีด + รังเรขาคณิต) + ใบงาน
ข้อสอบข้อ 1 ตรวจด้วยการสุ่มชุดจุดที่ผ่านเงื่อนไข แล้วดูว่ามีวงกลมรัศมี 1 รอบจุดหนึ่งที่ครอบ ≥ 13 จุด"""
import random, math
from math import ceil, dist, sqrt
from itertools import combinations

checks = []
def ok(name, c): checks.append((name, c))

def cond(P, d):        # สามจุดใด ๆ มีสองจุดห่างไม่เกิน d
    return all(min(dist(a, b), dist(a, c), dist(b, c)) <= d + 1e-12 for a, b, c in combinations(P, 3))
def best_disk(P, r):   # วงกลมรัศมี r ที่มีจุดศูนย์กลางเป็นจุดในชุด ครอบได้มากสุดกี่จุด
    return max(sum(dist(p, q) <= r + 1e-12 for q in P) for p in P)

random.seed(1)
def cluster(c, n, rad):
    out = []
    while len(out) < n:
        x, y = random.uniform(-rad, rad), random.uniform(-rad, rad)
        if x * x + y * y <= rad * rad: out.append((c[0] + x, c[1] + y))
    return out
trials_ok = True; tested = 0
for _ in range(150):
    k = random.randint(0, 25)
    far = random.uniform(1.2, 3)
    P = cluster((0, 0), k, 0.5) + cluster((far, random.uniform(-1, 1)), 25 - k, 0.5)
    if not cond(P, 1): continue
    tested += 1
    trials_ok &= best_disk(P, 1) >= 13
ok(f"ข้อสอบ 1: สุ่ม {tested} ชุดที่ผ่านเงื่อนไข ทุกชุดมีวงกลมรัศมี 1 ครอบ ≥ 13 จุด", trials_ok and tested > 50)
# ข้อพิสูจน์: ถ้า AB > 1 ทุกจุด C ห่าง A หรือ B ไม่เกิน 1 — ตรวจบนชุดสุ่ม
P = cluster((0, 0), 12, 0.5) + cluster((3, 0), 13, 0.5)
ok("ชุด 12 + 13 จุด (สองกลุ่มห่างกัน) ผ่านเงื่อนไข และครอบได้มากสุด 13 พอดี (13 แน่น)", cond(P, 1) and best_disk(P, 1) == 13)
A, B = P[0], P[-1]
ok("ชุดเดียวกัน: AB > 1 และทุกจุดอยู่ในวงกลมรอบ A หรือ B", dist(A, B) > 1 and all(dist(C, A) <= 1 or dist(C, B) <= 1 for C in P))
ok("ceil(25/2) = 13 · ceil(9/2) = 5 · ceil(7/2) = 4", (ceil(25/2), ceil(9/2), ceil(7/2)) == (13, 5, 4))

# รังแบบแบ่งพื้นที่
ok("สามเหลี่ยมด้านเท่า 1 แบ่งเป็น 4 รูปด้าน 1/2 (พื้นที่ 1 = 4 × 1/4)", abs(1 - 4 * (0.5 ** 2)) < 1e-12)
ok("ceil(5/4) = 2", ceil(5 / 4) == 2)
ok("สี่เหลี่ยมด้าน 4 แบ่งเป็น 16 ช่องด้าน 1 · เส้นทแยง = √2", 4 * 4 == 16 and abs(dist((0, 0), (1, 1)) - sqrt(2)) < 1e-12)
ok("ค่ายแบบฝึก 11: สี่เหลี่ยมด้าน 2 แบ่ง 4 ช่องด้าน 1 · 5 บ้าน · ห่างไม่เกิน √2", ceil(5 / 4) == 2)
ok("ค่ายแบบฝึก 8: สามเหลี่ยมด้าน 3 แบ่ง 9 รูปด้าน 1 · 10 จุด", 9 * 1 == (3 ** 2) and ceil(10 / 9) == 2)
# ค่ายตัวอย่าง 5: บ้าน 5 หลังบนเส้นยาว 3 · แบ่ง [0,1],(1,2],(2,3] · ceil(5/3) = 2
ok("ค่ายตัวอย่าง 5: ceil(5/3) = 2", ceil(5 / 3) == 2)
# ใบงาน: 9 จุด เงื่อนไข d = 2 ⇒ วงกลมรัศมี 2 ครอบ ≥ 5
good = True; t = 0
for _ in range(100):
    k = random.randint(0, 9)
    P = cluster((0, 0), k, 1) + cluster((random.uniform(2.5, 6), 0), 9 - k, 1)
    if cond(P, 2): t += 1; good &= best_disk(P, 2) >= 5
ok(f"ใบงาน: 9 จุด (d = 2) สุ่ม {t} ชุด ทุกชุดมีวงกลมรัศมี 2 ครอบ ≥ 5", good and t > 30)

# ใบงานข้อ 9: ส่วนสี่ของวงกลมรัศมี 1 — ระยะไกลสุดระหว่างจุดในชิ้น = √2 (ตรวจด้วยการสุ่มจุดในชิ้น)
pts = [(r * math.cos(t), r * math.sin(t)) for r, t in ((math.sqrt(random.random()), random.uniform(0, math.pi / 2)) for _ in range(400))] + [(1, 0), (0, 1), (0, 0)]
ok("ใบงาน 9: ระยะไกลสุดในส่วนสี่ของวงกลมรัศมี 1 = √2", abs(max(dist(a, b) for a, b in combinations(pts, 2)) - sqrt(2)) < 1e-9)
ok("ใบงาน 10: ช่องด้าน 1/3 · เส้นทแยง √2/3 · ceil(10/9) = 2", abs(dist((0, 0), (1/3, 1/3)) - sqrt(2) / 3) < 1e-12 and ceil(10 / 9) == 2)
# ใบงาน 11: สามเหลี่ยมที่จุดยอดอยู่ในสี่เหลี่ยมด้าน 1/2 มีพื้นที่ ≤ 1/8 (สุ่มตรวจ)
def area(a, b, c): return abs((b[0]-a[0])*(c[1]-a[1]) - (c[0]-a[0])*(b[1]-a[1])) / 2
mx = max(area(*[(random.uniform(0, .5), random.uniform(0, .5)) for _ in range(3)]) for _ in range(20000))
ok(f"ใบงาน 11: สุ่มสามเหลี่ยมในช่องด้าน 1/2 พื้นที่มากสุด {mx:.4f} ≤ 1/8 · ceil(9/4) = 3", mx <= 0.125 and ceil(9 / 4) == 3)
ok("ใบงาน 11: สามเหลี่ยมใช้มุมสามมุมของช่อง ได้ 1/8 พอดี", abs(area((0, 0), (.5, 0), (0, .5)) - 0.125) < 1e-12)

bad = [n for n, c in checks if not c]
for n, c in checks: print("✅" if c else "❌", n)
raise SystemExit(1 if bad else 0)
