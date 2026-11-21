"""บทที่ 31 · เช็คพอยต์อ่าน Tong บทที่ 3 — ตรวจทุกตัวเลขที่ถามในบท"""
import sympy as sp
from itertools import combinations
from math import comb
t, t1, t2 = sp.symbols('t t1 t2', real=True)
ok = 0
def check(n, c):
    global ok; assert c, n; ok += 1; print("✅", n)
# (3.18) vs (3.23): (-i)^2/2 × จัตุรัส = (-i)^2 × สามเหลี่ยม  (กรณีสลับได้ ทดสอบด้วย H(t)=cos t, e^t)
for H in [sp.cos(t), sp.exp(t), t**2 + 3]:
    h = lambda s: H.subs(t, s)
    sq = sp.integrate(h(t1), (t1, 0, t))**2
    tri = sp.integrate(sp.integrate(h(t1)*h(t2), (t2, 0, t1)), (t1, 0, t))
    check(f"(3.18)=(3.23) for {H}", sp.simplify((-sp.I)**2/2*sq - (-sp.I)**2*tri) == 0)
check("(-i)^2/2 = -1/2", (-sp.I)**2/2 == sp.Rational(-1, 2))
check("(-i)^2 = -1", (-sp.I)**2 == -1)
# (3.43): T(φ1φ2φ3φ4) — นับพจน์
check("contraction เดียว = C(4,2) = 6", comb(4, 2) == 6)
full = {frozenset(map(frozenset, p)) for p in [((1,2),(3,4)),((1,3),(2,4)),((1,4),(2,3))]}
check("จับครบ = 3", len(full) == 3)
check("รวม 1 + 6 + 3 = 10", 1 + 6 + 3 == 10)
# (3.74)-(3.75): λ=24 · 4! contraction · iA = -iλ ⇒ A = -λ
check("H_int coef +1, L coef -1", (24/24, -24/24) == (1.0, -1.0))
check("4!/4! = 1", sp.factorial(4)/24 == 1)
check("iA=-iλ ⇒ A=-λ", sp.simplify(-sp.I/sp.I) == -1)
# §3.5.4 connected: X, A, ปลา ต่อกับเส้นภายนอกทุกชิ้น · B มีห่วงลอยไม่ต่อ
def connected(edges, ext):
    nodes = {n for e in edges for n in e}; seen = set(); stack = [next(iter(ext))]
    while stack:
        n = stack.pop()
        if n in seen: continue
        seen.add(n); stack += [b for a, b in edges if a == n] + [a for a, b in edges if b == n]
    return seen == nodes
X = [("1","z"),("2","z"),("3","z"),("4","z")]
A = [("1","z"),("2","z"),("z","z")]
B = [("1","2"),("z","z"),("z","z")]
F = [("1","a"),("2","a"),("3","b"),("4","b"),("a","b"),("a","b")]
res = [connected(g, "1") for g in (X, A, B, F)]
check(f"connected X,A,B,fish = {res}", res == [True, True, False, True])
check("connected count = 3", sum(res) == 3)
print(f"{ok}/{ok}")
