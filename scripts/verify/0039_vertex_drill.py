"""ชุดฝึกจุดยอด — ไล่จับคู่ทุกวิธีจริง แล้วนับตามเงื่อนไขของแต่ละข้อ"""
from fractions import Fraction
from math import factorial

def pairings(xs):
    if not xs: yield []; return
    a, rest = xs[0], xs[1:]
    for i, b in enumerate(rest):
        for p in pairings(rest[:i] + rest[i+1:]):
            yield [(a, b)] + p

def count(ext, legs, cond=lambda p: True):
    return sum(1 for p in pairings(ext + legs) if cond(p))

L = lambda n: [f"v{i}" for i in range(n)]
isleg = lambda x: x.startswith("v")
def partner(p, x):
    for a, b in p:
        if a == x: return b
        if b == x: return a

ok = 0
def check(name, got, want):
    global ok
    assert got == want, (name, got, want); ok += 1; print("✅", name, got)

# ข้อ 1: φ⁴ สองจุด — ย้อนจาก 1/2
check("1 A=12", count(["1","2"], L(4), lambda p: isleg(partner(p,"1"))), 12)
check("1 24*(1/2)", 24*Fraction(1,2), 12)
# ข้อ 2: φ³ สามจุด รูป Y
yv = count(["1","2","3"], L(3), lambda p: all(isleg(partner(p,x)) for x in "123"))
check("2 total", count(["1","2","3"], L(3)), 15)
check("2 Y", yv, 6); check("2 /3!", Fraction(yv, factorial(3)), 1)
# ข้อ 3: φ³ หนึ่งจุด
t = count(["1"], L(3)); check("3 total", t, 3)
check("3 1->leg", count(["1"], L(3), lambda p: isleg(partner(p,"1"))), 3)
check("3 /3!", Fraction(t, 6), Fraction(1,2))
# ข้อ 4: φ³ สองจุด — คี่
check("4 odd", count(["1","2"], L(3)), 0)
# ข้อ 5: φ⁶ สองจุด ทั้งคู่จับขา
s = count(["1","2"], L(6), lambda p: isleg(partner(p,"1")) and isleg(partner(p,"2")))
check("5 ways", s, 90); check("5 /6!", Fraction(s, 720), Fraction(1,8))
# ข้อ 6: φ⁴ สี่จุด — 12 จับกัน 3,4 จับขา
c = count(["1","2","3","4"], L(4),
          lambda p: partner(p,"1")=="2" and isleg(partner(p,"3")) and isleg(partner(p,"4")))
check("6 ways", c, 12); check("6 /4!", Fraction(c, 24), Fraction(1,2))
print(f"{ok}/{ok}")
