import sympy as sp

t, m, g, l = sp.symbols('t m g l', positive=True)

def euler_lagrange(L, q):
    """สมการออยเลอร์-ลากรองจ์: d/dt(∂L/∂q̇) − ∂L/∂q = 0"""
    qd = sp.diff(q, t)
    return sp.simplify(sp.diff(sp.diff(L, qd), t) - sp.diff(L, q))

print("=== 1. วัตถุตกอิสระ: L = ½mẏ² − mgy ===")
y = sp.Function('y')(t)
L1 = sp.Rational(1,2)*m*sp.diff(y,t)**2 - m*g*y
eq1 = euler_lagrange(L1, y)
print("  EL:", sp.Eq(eq1, 0))
print("  แก้ได้:", sp.solve(eq1, sp.diff(y,t,2)), " → ตรงกับ F = ma:", sp.solve(eq1, sp.diff(y,t,2)) == [-g])

print("\n=== 2. ลูกตุ้ม: L = ½ml²θ̇² + mgl·cos θ ===")
th = sp.Function('theta')(t)
L2 = sp.Rational(1,2)*m*l**2*sp.diff(th,t)**2 + m*g*l*sp.cos(th)
sol2 = sp.solve(euler_lagrange(L2, th), sp.diff(th,t,2))
print("  θ̈ =", sol2[0], " → ตรงกับสมการลูกตุ้มมาตรฐาน:", sp.simplify(sol2[0] + g*sp.sin(th)/l) == 0)

print("\n=== 3. ทำไมต้องเป็น T − V ไม่ใช่ T + V ===")
L3 = sp.Rational(1,2)*m*sp.diff(y,t)**2 + m*g*y      # ลองใช้ T + V
print("  ใช้ T+V ได้:", sp.solve(euler_lagrange(L3, y), sp.diff(y,t,2)), "← เครื่องหมายกลับด้าน ผิด")
import sympy as sp
t, a = sp.symbols('t alpha', real=True)
m, g, T = 1, 10, 2                      # kg, m/s², s

y_true = sp.Rational(1,2)*g*t*(T - t)   # เส้นทางจริง: ขึ้นแล้วตกกลับที่เดิม
bump   = sp.sin(sp.pi*t/T)              # การรบกวน: เป็นศูนย์ที่ปลายทั้งสองข้าง
y      = y_true + a*bump

L = sp.Rational(1,2)*m*sp.diff(y,t)**2 - m*g*y
S = sp.integrate(L, (t, 0, T))
S = sp.simplify(sp.expand(S))

print("เส้นทางจริง y(t) =", y_true, " → จุดสูงสุด", y_true.subs(t,1), "m ที่ t = 1 s")
print("\nS(α) =", S)
print("S(0) =", S.subs(a,0), "J·s")
print("dS/dα ที่ α=0 :", sp.diff(S,a).subs(a,0), "  ← ศูนย์ = จุดวิกฤต (stationary)")
print("d²S/dα²      :", sp.diff(S,a,2), "  ← เป็นบวก = ค่าต่ำสุดจริง ไม่ใช่ค่าสูงสุด")
print()
for A in [-2,-1,-sp.Rational(1,2),0,sp.Rational(1,2),1,2]:
    v = sp.nsimplify(S.subs(a,A)); print(f"  α = {float(A):+5.2f} → S = {float(v):9.4f}")

# ตรวจซ้ำด้วยการอินทิเกรตเชิงตัวเลข (คนละวิธีกับ symbolic)
import numpy as np
def S_num(al, n=200001):
    tt = np.linspace(0, 2, n)
    yy = 5*tt*(2-tt) + al*np.sin(np.pi*tt/2)
    vv = np.gradient(yy, tt)
    return np.trapezoid(0.5*vv**2 - 10*yy, tt)
print("\nตรวจซ้ำเชิงตัวเลข: S(0) =", round(S_num(0),4), " S(1) =", round(S_num(1),4),
      " S(-1) =", round(S_num(-1),4))
