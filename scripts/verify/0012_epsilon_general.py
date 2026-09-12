"""ที่มาของสัมประสิทธิ์ ε — ตรวจว่าสูตรทั่วไปมี m อยู่ใน "ทั้งสอง" พจน์
และตรวจว่าสูตรที่พิมพ์ไว้ในบท 0009 ตกตัว m ไปหนึ่งที่จริงหรือไม่
"""
import sympy as sp

t, eps, m, g = sp.symbols('t varepsilon m g', positive=True)

def eps_coeff(y0, eta, V, a, b, mass):
    """กางแอ็กชันจริงแล้วหยิบสัมประสิทธิ์ของ eps — ไม่ได้ใช้สูตรลัดเลย"""
    y = y0 + eps*eta
    L = sp.Rational(1, 2)*mass*sp.diff(y, t)**2 - V(y)
    S = sp.integrate(sp.expand(L), (t, a, b))
    return sp.nsimplify(sp.expand(S).coeff(eps, 1))

# ---- 1) สูตรทั่วไป
y0f, etaf = sp.Function('y0')(t), sp.Function('eta')(t)
lin = sp.expand(sp.Rational(1,2)*m*sp.diff(y0f + eps*etaf, t)**2
                - m*g*(y0f + eps*etaf)).coeff(eps, 1)
print('พจน์เชิงเส้นใน ε (V = mgy):', lin)
print('  ⇒  m·ẏ0·η̇ − m·g·η   — m อยู่ทั้งสองพจน์\n')

# ---- 2) โจทย์ B ของบท 0011 : m=5, y0=2t, eta=t(1-t), [0,1], ไม่มี V
b = eps_coeff(2*t, t*(1-t), lambda Y: 0, 0, 1, 5)
short = sp.integrate(5*sp.diff(2*t, t)*sp.diff(t*(1-t), t), (t, 0, 1))
print('โจทย์ B  สัมประสิทธิ์ ε =', b, ' · สูตรลัด ∫m ẏ0 η̇ dt =', short)
assert b == 0 == short

# ---- 3) เคสที่ "ตกตัว m" ให้คำตอบผิดจริง — ต้องเลือก η ที่ ∫η̇ ≠ 0 ไม่ได้
#     (η ปลายตรึงเสมอ ⇒ ∫η̇ = 0) จึงต้องใช้ y0 ที่ ẏ0 ไม่คงที่
y0c, etac, M, G = t**2, t*(1-t), 5, 10
true = eps_coeff(y0c, etac, lambda Y: M*G*Y, 0, 1, M)
withm = (sp.integrate(M*sp.diff(y0c,t)*sp.diff(etac,t), (t,0,1))
         - sp.integrate(M*G*etac, (t,0,1)))
nom   = (sp.integrate(  sp.diff(y0c,t)*sp.diff(etac,t), (t,0,1))
         - sp.integrate(M*G*etac, (t,0,1)))
print('\nเคส y0=t², m=5, V=mgy')
print('  ค่าจริง (กางแอ็กชัน) =', true)
print('  สูตรที่มี m ครบ      =', withm, ' ตรง:', withm == true)
print('  สูตรที่ตกตัว m       =', nom,   ' ตรง:', nom == true, '  ← ผิด')

# ---- 4) ตัวอย่างของบท 0009 เอง (ตาราง T = ½ẏ² ⇒ m = 1, V = 10y)
y09, eta9 = 8*t - 4*t**2, t*(2-t)
true9 = eps_coeff(y09, eta9, lambda Y: 10*Y, 0, 2, 1)
printed9 = (sp.integrate(sp.diff(y09,t)*sp.diff(eta9,t), (t,0,2))
            - sp.integrate(m*10*eta9, (t,0,2)))       # <-- ตามที่พิมพ์ไว้จริง
print('\nบท 0009 (m=1):')
print('  ค่าจริง            =', true9)
print('  สูตรที่พิมพ์ไว้ในบท =', printed9, '  ← มี m ค้างอยู่ในพจน์เดียว')
print('  เท่ากันเมื่อ m =', sp.solve(sp.Eq(printed9, true9), m))
print('\n⇒ สูตรในกล่องของบท 0009 สอดคล้องเฉพาะตอน m = 1 เท่านั้น')
print('  ถ้า m ≠ 1 มันผิด และนั่นคือสิ่งที่โจทย์ B (m=5) ขอ')
