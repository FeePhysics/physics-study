"""ตัวเลขทุกตัวในชุดฝึกบทที่ 1 — คำนวณจริง ไม่ได้ประมาณ"""
import sympy as sp
t = sp.symbols('t', positive=True)
m = 1

print("=== อนุภาคอิสระ: จาก x=0 ถึง x=10 ใน 5 วินาที (V=0) ===")
paths = {
    "ความเร็วคงที่  x=2t":      2*t,
    "เร่งสม่ำเสมอ   x=0.4t²":   sp.Rational(2,5)*t**2,
    "ช้าก่อนแล้วพุ่ง x=10(t/5)³": 10*(t/5)**3,
}
for name, x in paths.items():
    v = sp.diff(x, t)
    S = sp.integrate(sp.Rational(1,2)*m*v**2, (t, 0, 5))
    print(f"  {name:26} x(5)={x.subs(t,5)}  S={S} = {float(S):.2f}")

print("\n=== ลูกตุ้ม: หา 'พิกัดวัฏจักร' — L ไม่มี φ ⇒ อนุรักษ์โมเมนตัมเชิงมุม ===")
th, ph = sp.Function('theta')(t), sp.Function('phi')(t)
M, g, l = sp.symbols('m g l', positive=True)
# ลูกตุ้มทรงกลม: L = ½ml²(θ̇² + sin²θ·φ̇²) + mgl·cosθ
L = sp.Rational(1,2)*M*l**2*(sp.diff(th,t)**2 + sp.sin(th)**2*sp.diff(ph,t)**2) + M*g*l*sp.cos(th)
print("  ∂L/∂φ =", sp.diff(L, ph), " ⇒ φ เป็นพิกัดวัฏจักร")
p_phi = sp.simplify(sp.diff(L, sp.diff(ph,t)))
print("  ปริมาณที่อนุรักษ์ p_φ =", p_phi)

print("\n=== ตรวจว่า S ของเส้นทางจริงเป็น 'ต่ำสุด' ไม่ใช่ 'ศูนย์' ===")
y = 5*t*(2-t)
S_true = sp.integrate(sp.Rational(1,2)*sp.diff(y,t)**2 - 10*y, (t,0,2))
print(f"  S ของเส้นทางจริงในบทที่ 1 = {S_true} = {float(S_true):.2f}  ← ไม่ใช่ 0")
