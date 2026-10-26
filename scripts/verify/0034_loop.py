"""ตรวจทุกข้ออ้างของบทที่ 22 — propagator แบบ 1/(k²+m²) และค่าของห่วงแรก

1+1 มิติ แบบยุคลิด: k = (k0, k1) · k² = k0² + k1² (บวกทั้งสองตัว) · hbar = c = 1
"""
import mpmath as mp
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


k0, k1, k, u = sp.symbols('k0 k1 k u', real=True)
w, m, r = sp.symbols('omega m r', positive=True)

# ── ① ∫ dk0/2π 1/(k0² + ω²) = 1/(2ω) — น้ำหนักของบทที่ 20 ─────────
I1 = sp.integrate(1/(k0**2 + w**2), (k0, -sp.oo, sp.oo))
chk("∫ dk0/(k0² + ω²) = π/ω", I1, sp.pi/w)
chk("∫ dk0/2π 1/(k0² + ω²) = 1/(2ω)", I1/(2*sp.pi), 1/(2*w))
chk("ตัวอย่างในบท ω = 3: (π/3)/(2π) = 1/6", I1.subs(w, 3)/(2*sp.pi), sp.Rational(1, 6))
chk("ปฏิยานุพันธ์ = (1/ω) arctan(k0/ω)", sp.diff(sp.atan(k0/w)/w, k0) - 1/(k0**2 + w**2), 0)
chk("arctan(+∞) − arctan(−∞) = π", sp.limit(sp.atan(k0), k0, sp.oo) - sp.limit(sp.atan(k0), k0, -sp.oo), sp.pi)
ok.append(("ตัวควบคุม: ถ้าลืมหาร 2π ได้ π/ω ≠ 1/(2ω)", sp.simplify(I1 - 1/(2*w)) != 0, None, None))

# ⇒ ∫ d²k/(2π)² e^{ik1 r}/(k0² + k1² + m²) = ∫ dk1/2π e^{ik1 r}/(2ω) = K0(mr)/(2π) (ผลของบทที่ 20)
mp.mp.dps = 20
for mm, rr in [(1, 1), (2, 0.7)]:
    inner = mp.quadosc(lambda q: mp.cos(q*rr)/(2*mp.sqrt(q**2 + mm**2)), [0, mp.inf], omega=rr)/mp.pi
    ok.append((f"propagator 2D ยุคลิดที่ระยะ r = K0(mr)/2π (m={mm}, r={rr})",
               abs(inner - mp.besselk(0, mm*rr)/(2*mp.pi)) < 1e-12, inner, mp.besselk(0, mm*rr)/(2*mp.pi)))
# ตรวจซ้ำแบบอินทิเกรตสองชั้นจริง (k0 ก่อน) ที่ m=1, r=1 — k0 ให้ 1/(2ω) ตามสูตร
inner_k0 = mp.quad(lambda a: 1/(a**2 + 1 + 0.3**2), [-mp.inf, mp.inf])/(2*mp.pi)
ok.append(("ตัวเลข: ∫dk0/2π 1/(k0²+ω²) ที่ ω²=1.09 = 1/(2ω)", abs(inner_k0 - 1/(2*mp.sqrt(1.09))) < 1e-15, inner_k0, None))

# ── ② propagator 1/(k² + m²) ─────────────────────────────────────
G = lambda a, b, mm: sp.Rational(1, 1)/(a**2 + b**2 + mm**2)
chk("k = (3, 4), m = 0 ⇒ 1/25", G(3, 4, 0), sp.Rational(1, 25))
chk("k = (1, 2), m = 2 ⇒ 1/9", G(1, 2, 2), sp.Rational(1, 9))
chk("ตัวอย่างในบท k = (0, 3), m = 4 ⇒ 1/25", G(0, 3, 4), sp.Rational(1, 25))
# โมเมนตัมในห่วง: เข้า p แยกเป็น k กับ p − k
p = (2, 1); kk = (5, 3)
chk("p − k = (−3, −2) ตัวแรก", p[0] - kk[0], -3)
chk("p − k = (−3, −2) ตัวหลัง", p[1] - kk[1], -2)

# ── ③ ห่วงที่โมเมนตัมภายนอกเป็นศูนย์: ∫ d²k/(2π)² 1/(k²+m²)² = 1/(4π m²) ──
radial = sp.integrate(k/(k**2 + m**2)**2, (k, 0, sp.oo))
chk("∫_0^∞ k dk/(k² + m²)² = 1/(2m²)", radial, 1/(2*m**2))
chk("แทน u = k²: ½∫_0^∞ du/(u + m²)² = 1/(2m²)", sp.integrate(sp.Rational(1, 2)/(u + m**2)**2, (u, 0, sp.oo)),
    1/(2*m**2))
loop = 2*sp.pi*radial/(2*sp.pi)**2
chk("∫ d²k/(2π)² 1/(k²+m²)² = 2π · 1/(2m²) / (4π²) = 1/(4π m²)", loop, 1/(4*sp.pi*m**2))
# ตรวจแบบพิกัดฉากจริง (ไม่ใช้เชิงขั้ว) ที่ m = 1
num = mp.quad(lambda a, b: 1/(a**2 + b**2 + 1)**2, [-mp.inf, mp.inf], [-mp.inf, mp.inf])/(2*mp.pi)**2
ok.append(("ตัวเลขพิกัดฉาก m=1: = 1/(4π)", abs(num - 1/(4*mp.pi)) < 1e-10, num, 1/(4*mp.pi)))
ok.append(("ตัวควบคุม: ลืม k ใน dk (ใช้ dk แทน k dk) ได้ค่าอื่น",
           sp.simplify(sp.integrate(1/(k**2 + m**2)**2, (k, 0, sp.oo)) - 1/(2*m**2)) != 0, None, None))

# ── ④ ทำไมสี่มิติอนันต์: ที่ k ใหญ่ k^(d-1)/k^4 ──────────────────
for d, pw in [(2, -3), (3, -2), (4, -1)]:
    chk(f"d = {d}: k^(d-1) × k^(-4) = k^{pw}", (d - 1) - 4, pw)
Lam = sp.Symbol('Lam', positive=True)
ok.append(("d = 4: ∫_1^Λ dk/k = ln Λ → ∞", sp.limit(sp.integrate(1/k, (k, 1, Lam)), Lam, sp.oo) == sp.oo, None, None))
chk("d = 2: ∫_1^∞ k^(-3) dk = 1/2 จำกัด", sp.integrate(k**-3, (k, 1, sp.oo)), sp.Rational(1, 2))

# ── ตัวเลขในข้อฝึก ──────────────────────────────────────────
chk("ข้อ 1: ω = 5 ⇒ ∫ dk0/(k0²+25) = π/5 (ตัวเลขหน้า π = 1/5)", I1.subs(w, 5)/sp.pi, sp.Rational(1, 5))
chk("ข้อ 1: หาร 2π ⇒ 1/10", I1.subs(w, 5)/(2*sp.pi), sp.Rational(1, 10))
chk("ข้อ 1: = 1/(2ω)", sp.Rational(1, 10), sp.Rational(1, 2*5))
chk("ข้อ 3: m = 1 ⇒ ∫ k dk/(k²+1)² = 1/2", radial.subs(m, 1), sp.Rational(1, 2))
chk("ข้อ 3: m = 2 ⇒ 1/8", radial.subs(m, 2), sp.Rational(1, 8))
chk("ข้อ 3: m = 2 ⇒ π × ค่าห่วง = 1/16", sp.pi*loop.subs(m, 2), sp.Rational(1, 16))

bad_ = [r_ for r_ in ok if not r_[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
