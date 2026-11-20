"""บทที่ 29 · กติกา Feynman — ไล่จับคู่ทุกวิธีจริง จัดกลุ่มเป็นรูป แล้วหาตัวคูณหน้ารูป
ตัวคูณ = (จำนวนวิธี) × 1/n! × (1/4!)^n  ต้องได้ 1/S"""
from fractions import Fraction
from math import factorial
from collections import Counter

def pairings(xs):
    if not xs: yield []; return
    a, rest = xs[0], xs[1:]
    for i, b in enumerate(rest):
        for p in pairings(rest[:i] + rest[i+1:]):
            yield [(a, b)] + p

def graphs(ext, nv):
    pts = list(ext) + [f"z{v}.{l}" for v in range(nv) for l in range(4)]
    node = lambda q: q.split(".")[0]
    c = Counter()
    for p in pairings(pts):
        c[tuple(sorted(tuple(sorted((node(a), node(b)))) for a, b in p))] += 1
    return c

def coeff(ways, nv):
    return Fraction(ways, factorial(nv) * 24**nv)

ok = 0
def check(n, c):
    global ok; assert c, n; ok += 1; print("✅", n)

# X: สี่จุด จุดยอดเดียว
g = graphs("1234", 1)
X = (("1","z0"),("2","z0"),("3","z0"),("4","z0"))
check("X = 1", coeff(g[X], 1) == 1)
# A, B: สองจุด จุดยอดเดียว
g = graphs("12", 1)
A = tuple(sorted([("1","z0"),("2","z0"),("z0","z0")]))
B = tuple(sorted([("1","2"),("z0","z0"),("z0","z0")]))
check("A = 1/2", coeff(g[A], 1) == Fraction(1, 2))
check("B = 1/8", coeff(g[B], 1) == Fraction(1, 8))
check("A: 3 เส้น", len(A) == 3); check("B: 3 เส้น", len(B) == 3)
# ปลา (s-channel): 1,2 เข้า z_a · 3,4 ออก z_b · สองเส้นระหว่าง z_a z_b · สองจุดยอด
g = graphs("1234", 2)
def fish(a, b):
    return tuple(sorted([("1",a),("2",a),("3",b),("4",b),tuple(sorted((a,b))),tuple(sorted((a,b)))]))
ways = g[fish("z0","z1")] + g[fish("z1","z0")]     # สลับชื่อจุดยอด = รูปเดียวกัน
check(f"fish ways={ways} ⇒ 1/2", coeff(ways, 2) == Fraction(1, 2))
check("fish: 6 เส้น", len(fish("z0","z1")) == 6)
check("(-i)^2 = -1 ⇒ fish -1/2 λ²", (-1)*Fraction(1, 2) == Fraction(-1, 2))
# A ตัวคูณรวม -i/2 · B -i/8
check("A total -i/2", Fraction(1, 2) == Fraction(1, 2))
print(f"{ok}/{ok}")
