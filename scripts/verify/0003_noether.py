"""ทุกข้ออ้างในบทที่ 2 — คำนวณจริงด้วย sympy ไม่ได้คัดลอกจากตำรา"""
import sympy as sp

t = sp.symbols('t')
m, k, g, l, w = sp.symbols('m k g l omega', positive=True)

def EL(L, q):
    return sp.simplify(sp.diff(sp.diff(L, sp.diff(q,t)), t) - sp.diff(L, q))

def partial_t(L, q):
    """∂L/∂t = อนุพันธ์เทียบ t ที่โผล่ **ตรง ๆ** เท่านั้น

    ⚠️ sp.diff(L, t) ให้อนุพันธ์ *รวม* — sympy รู้ว่า q(t) ขึ้นกับ t จึงไล่เข้าไป
    ในทุกพจน์ · ต้องตรึง q และ q̇ เป็นสัญลักษณ์อิสระก่อน แล้วค่อย diff
    (subs q̇ ก่อน q เสมอ ไม่งั้น q(t) ข้างใน Derivative โดนแทนไปก่อน)
    """
    Q, Qd = sp.symbols('Q Qd')
    Ls = L.subs(sp.diff(q, t), Qd).subs(q, Q)
    return sp.diff(Ls, t).subs({Q: q, Qd: sp.diff(q, t)})


def energy(L, qs):
    """H = Σ q̇·∂L/∂q̇ − L  (พลังงานที่มาจากสมมาตรการเลื่อนเวลา)"""
    return sp.simplify(sum(sp.diff(q,t)*sp.diff(L, sp.diff(q,t)) for q in qs) - L)

print("=== ① เลื่อนที่ได้ (L ไม่มี x) ⇒ โมเมนตัมอนุรักษ์ ===")
x, y = sp.Function('x')(t), sp.Function('y')(t)
L1 = sp.Rational(1,2)*m*(sp.diff(x,t)**2 + sp.diff(y,t)**2) - m*g*y   # แรงโน้มถ่วงแนวดิ่ง
print("  ∂L/∂x =", sp.diff(L1, x), "⇒ x เป็นพิกัดวัฏจักร")
px = sp.diff(L1, sp.diff(x,t))
print("  p_x =", px, " · d/dt(p_x) =", sp.simplify(sp.diff(px,t).subs(sp.diff(x,t,2), sp.solve(EL(L1,x), sp.diff(x,t,2))[0])))
print("  แต่แนวดิ่งไม่อนุรักษ์: ∂L/∂y =", sp.diff(L1, y), "≠ 0 ✓ (มีแรงโน้มถ่วง)")

print("\n=== ② หมุนได้ (V ขึ้นกับ r เท่านั้น) ⇒ โมเมนตัมเชิงมุมอนุรักษ์ ===")
r, th = sp.Function('r')(t), sp.Function('theta')(t)
V = sp.Function('V')(r)
L2 = sp.Rational(1,2)*m*(sp.diff(r,t)**2 + r**2*sp.diff(th,t)**2) - V
print("  ∂L/∂θ =", sp.diff(L2, th), "⇒ θ เป็นพิกัดวัฏจักร")
print("  ปริมาณอนุรักษ์ p_θ =", sp.diff(L2, sp.diff(th,t)), " ← โมเมนตัมเชิงมุม")
print("  ยืนยัน EL ของ θ ยุบเป็น d/dt(p_θ) = 0 :", sp.simplify(EL(L2,th) - sp.diff(sp.diff(L2, sp.diff(th,t)), t)) == 0)

print("\n=== ③ เลื่อนเวลาได้ (L ไม่มี t โผล่ตรง ๆ) ⇒ พลังงานอนุรักษ์ ===")
q = sp.Function('q')(t)
L3 = sp.Rational(1,2)*m*sp.diff(q,t)**2 - sp.Rational(1,2)*k*q**2          # สปริงปกติ
H3 = energy(L3, [q])
dH3 = sp.simplify(sp.diff(H3,t).subs(sp.diff(q,t,2), sp.solve(EL(L3,q), sp.diff(q,t,2))[0]))
print("  H =", H3, "  ← T + V พอดี")
print("  dH/dt เมื่อใช้สมการการเคลื่อนที่ =", dH3, " ⇒ อนุรักษ์" if dH3 == 0 else " ⇒ ไม่อนุรักษ์")

print("\n=== ④ เคสที่พัง: สปริงที่ค่าคงที่เปลี่ยนตามเวลา (L มี t โผล่ตรง ๆ) ===")
L4 = sp.Rational(1,2)*m*sp.diff(q,t)**2 - sp.Rational(1,2)*k*sp.exp(w*t)*q**2
H4 = energy(L4, [q])
dH4 = sp.simplify(sp.diff(H4,t).subs(sp.diff(q,t,2), sp.solve(EL(L4,q), sp.diff(q,t,2))[0]))
pt = sp.simplify(partial_t(L4, q))
print("  ∂L/∂t =", pt, "≠ 0")
print("  dH/dt =", dH4, " ⇒ พลังงานไม่อนุรักษ์ · สมมาตรหายไป กฎอนุรักษ์หายตาม")
print("  ตรงกับสูตร dH/dt = −∂L/∂t :", sp.simplify(dH4 + pt) == 0)
print("\n  (เทียบ: กรณีสปริงปกติ ∂L/∂t =", sp.simplify(partial_t(L3, q)), "⇒ พลังงานอนุรักษ์)")
