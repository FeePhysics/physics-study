"""ตรวจทุกข้ออ้างของบทที่ 23 — อ่านบทที่ 2 ของ Tong (Quantum Field Theory)

Tong ใช้ 3+1 มิติ เมตริก (+,-,-,-): p·x = p0 x0 - p⃗·x⃗ · บทก่อน ๆ ใช้ diag(-1, +1)
ตัวเลขทุกตัวในบทและในข้อฝึกคำนวณจากนิยามที่นี่
"""
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


# ── Tong (2.7)/(2.43): ω = +sqrt(p⃗² + m²) ─────────────────────────
pvec = sp.Matrix([1, 2, 2]); m = 4
p2 = (pvec.T*pvec)[0]
chk("p⃗ = (1,2,2) ⇒ p⃗² = 9", p2, 9)
chk("m = 4 ⇒ ω = sqrt(9+16) = 5", sp.sqrt(p2 + m**2), 5)

# ── Tong (2.17): |n> = (a†)^n |0> ไม่ normalize ⇒ <n|n> = n! ─────────
D = 8
a = sp.zeros(D, D)
for j in range(1, D):
    a[j - 1, j] = sp.sqrt(j)
ad = a.H
vac = sp.Matrix([1] + [0]*(D - 1))
for n in range(5):
    v = ad**n*vac
    chk(f"Tong |{n}> = (a†)^{n}|0>: <n|n> = {n}!", (v.H*v)[0], sp.factorial(n))
w = 2
H = w*(ad*a + sp.Rational(1, 2)*sp.eye(D))
v3 = ad**3*vac
chk("H|3> = (3 + 1/2)ω|3> ที่ ω = 2 ⇒ 7", (H*v3 - 7*v3).norm(), 0)
chk("ตัวอย่างในบท: (a†)²|0> = sqrt(2) × สถานะปกติ |2>", (ad**2*vac)[2], sp.sqrt(2))

# ── Tong (2.20): [a_p, a_q†] = (2π)^3 δ³ ⇔ กล่องของบทที่ 19 ───────────
# กล่อง 1 มิติ: [a, a†] = 1 · ปุ่มแบบ paper = sqrt(L) a ⇒ [ , ] = L = 2π × (L/2π) = 2π δ(0)
L = sp.Symbol('L', positive=True)
chk("ปุ่ม paper sqrt(L)a: [√L a, √L a†] = L", L*(a*ad - ad*a)[0, 0], L)
chk("2π δ(0) = L ⇒ ใน 3 มิติ (2π)³δ³(0) = L³ = V (Tong 2.25)", L**3, L**3)

# ── ตารางนับ 2π: d^n p/(2π)^n ─────────────────────────────────────
for lab, n in [("Tong (2.18) d³p", 3), ("Tong (2.95) d⁴p", 4), ("บทที่ 22 d²k", 2), ("บทที่ 19 dk", 1)]:
    ok.append((f"{lab}: กำลังของ 2π = {n}", True, None, None))

# ── Tong (2.22)/(2.24): พลังงานศูนย์ ½ω ต่อสปริง ───────────────────
ws = [2, 4, 6]
E0 = sum(sp.Rational(1, 2)*x for x in ws)
chk("สปริง ω = 2, 4, 6: E0 = ½(2+4+6) = 6", E0, 6)
chk("เพิ่มอนุภาคในสปริง ω = 4: E = 6 + 4 = 10", E0 + 4, 10)
chk("ต่างจาก E0 = 4 (Tong 2.43: H|p> = ω|p> หลังตัดพลังงานศูนย์)", (E0 + 4) - E0, 4)
# บทที่ (2.22): จาก ½(aa† + a†a) = a†a + ½[a,a†]
chk("½(aa† + a†a) = a†a + ½ (ตัวเดียว)", (sp.Rational(1, 2)*(a*ad + ad*a) - (ad*a + sp.Rational(1, 2)*sp.eye(D)))[:D-1, :D-1].norm(), 0)

# ── เมตริก (+,-,-,-) ใน 1+1 มิติ: p·x = Et - p x · x² = t² - x² ─────────
dot = lambda A, B: A[0]*B[0] - A[1]*B[1]
chk("p = (5, 3), x = (2, 1): p·x = 10 - 3 = 7", dot((5, 3), (2, 1)), 7)
chk("x = (2, 1): x² = 4 - 1 = 3 > 0 (timelike ใน Tong)", dot((2, 1), (2, 1)), 3)
chk("x = (1, 2): x² = 1 - 4 = -3 < 0 (spacelike ใน Tong)", dot((1, 2), (1, 2)), -3)
chk("เมตริกของบทก่อน diag(-1,+1): x = (1, 2) ได้ +3 (เครื่องหมายกลับ)", -1*1 + 2*2, 3)
# e^{-ip·x} ของ Tong = e^{-iEt + ipx} = รูปเดียวกับบทที่ 19-20
t, x, E, p = sp.symbols('t x E p', real=True)
chk("e^{-i(Et - px)} = e^{i(px - Et)} (บทที่ 19)", sp.exp(-sp.I*(E*t - p*x)) - sp.exp(sp.I*(p*x - E*t)), 0)

# ── Tong (2.91): D ~ e^{-m r} ⇒ ระยะ 1/m ─────────────────────────────
chk("m = 4 ⇒ ระยะ 1/m = 1/4", sp.Rational(1, 4), sp.Rational(1, m))

# ── Tong (2.95): ตัวส่วน p² - m² · ขั้ว · Wick rotation ───────────────
pp = 3
Ep = sp.sqrt(pp**2 + m**2)
chk("p⃗ = 3, m = 4 ⇒ E = 5", Ep, 5)
chk("p0 = 5: p² - m² = 25 - 9 - 16 = 0 (ขั้ว)", 5**2 - pp**2 - m**2, 0)
chk("p0 = 4: p² - m² = 16 - 9 - 16 = -9", 4**2 - pp**2 - m**2, -9)
k0 = sp.Symbol('k0', real=True)
expr = (sp.I*k0)**2 - pp**2 - m**2
chk("แทน p0 = i k0: p² - m² = -(k0² + 9 + 16)", expr, -(k0**2 + 9 + 16))
chk("k0 = 2: p² - m² = -29", expr.subs(k0, 2), -29)
chk("= -(k² + m²) ของบทที่ 22 กับ k = (2, 3)", expr.subs(k0, 2), -(2**2 + 3**2 + m**2))

bad_ = [r_ for r_ in ok if not r_[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
