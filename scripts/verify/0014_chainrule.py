"""ตรวจทุกตัวเลขของบทซ่อม "กฎที่เคยพูดผ่าน ๆ"

หัวใจ: สัมประสิทธิ์ของ ε มาจาก **กฎลูกโซ่** บรรทัดเดียว ไม่ใช่สูตรที่ต้องท่อง
    สัมประสิทธิ์ = ∫ [ (∂L/∂y)·η + (∂L/∂ẏ)·η̇ ] dt
และเครื่องหมายลบของ V โผล่เอง จาก ∂L/∂y = −V'(y)
"""
import sympy as sp
t, eps = sp.symbols('t varepsilon')
ok = lambda n, got, want: print(f"{'✓' if sp.simplify(got-want)==0 else '✗'} {n:46s} = {got}")

def coeff_direct(Lfun, y0, eta, a, b):
    """กางแอ็กชันจริงแล้วหยิบสัมประสิทธิ์ ε — เส้นทางอ้างอิง"""
    y = y0 + eps*eta
    return sp.expand(sp.integrate(sp.expand(Lfun(y, sp.diff(y,t))), (t,a,b))).coeff(eps,1)

def coeff_chain(Lfun, y0, eta, a, b):
    """ใช้กฎลูกโซ่ — เส้นทางที่บทนี้สอน (ต้องได้เท่ากัน)"""
    Y, V_ = sp.symbols('Y Vd')
    L = Lfun(Y, V_)
    dLdy  = sp.diff(L, Y ).subs({Y: y0, V_: sp.diff(y0,t)})
    dLdyd = sp.diff(L, V_).subs({Y: y0, V_: sp.diff(y0,t)})
    return sp.integrate(sp.expand(dLdy*eta + dLdyd*sp.diff(eta,t)), (t,a,b))

print('── ข้อ 1 · อ่าน η ออกจากนิพจน์ ' + '─'*30)
ye = t**3 + eps*t*(2-t)
eta1 = sp.expand(sp.diff(ye, eps))          # η = สัมประสิทธิ์ของ ε
ok('η', eta1, 2*t - t**2)
ok('η̇', sp.diff(eta1,t), 2 - 2*t)
print('   ปลายตรึงจริงไหม:', sp.simplify(eta1.subs(t,0)), sp.simplify(eta1.subs(t,2)))

print('\n── ข้อ 2 · โจทย์เดิมที่พลาด · m=3, V=mgy, g=10 ' + '─'*13)
m, g = 3, 10
Lg = lambda Y, Yd: sp.Rational(1,2)*m*Yd**2 - m*g*Y
y0, eta = t**2, t*(1-t)
Y, Vd = sp.symbols('Y Vd')
ok('∂L/∂ẏ  แทน y0 แล้ว', sp.diff(Lg(Y,Vd),Vd).subs({Y:y0,Vd:sp.diff(y0,t)}), 6*t)
ok('∂L/∂y', sp.diff(Lg(Y,Vd),Y), -30)
i_y  = sp.integrate(-30*eta, (t,0,1))
i_yd = sp.integrate(6*t*sp.diff(eta,t), (t,0,1))
ok('∫ (∂L/∂y)·η dt', i_y, -5)
ok('∫ (∂L/∂ẏ)·η̇ dt', i_yd, -1)
ok('รวม', i_y+i_yd, -6)
ok('กฎลูกโซ่ = กางแอ็กชันจริง',
   coeff_chain(Lg,y0,eta,0,1), coeff_direct(Lg,y0,eta,0,1))
print('   ⇒ −5 ที่เคยตอบคือค่าที่ถูกของพจน์นี้ เมื่อเขียนในรูปกฎลูกโซ่')

print('\n── ข้อ 3 · V คนละแบบ (สปริง) m=2, k=6 ' + '─'*22)
Ls = lambda Y, Yd: sp.Rational(1,2)*2*Yd**2 - sp.Rational(1,2)*6*Y**2
y0b, etab = t, t*(1-t)
ok('∂L/∂y = −ky แทน y0', sp.diff(Ls(Y,Vd),Y).subs({Y:y0b}), -6*t)
ok('∂L/∂ẏ แทน y0', sp.diff(Ls(Y,Vd),Vd).subs({Vd:sp.diff(y0b,t)}), 2)
c3 = coeff_chain(Ls, y0b, etab, 0, 1)
ok('สัมประสิทธิ์ ε', c3, sp.Rational(-1,2))
ok('  กางแอ็กชันจริงได้เท่ากัน', coeff_direct(Ls,y0b,etab,0,1), c3)

print('\n── ข้อ 4-5 · ตรวจสมมาตรก่อนตอบ ' + '─'*29)
mm,k,K,C = sp.symbols('m k K C', positive=True)
def check(L, xs, dq, label):
    s = sp.simplify(sum(sp.diff(L,x)*d for x,d in zip(xs,dq)))
    Q = sp.simplify(sum(sp.diff(L,sp.diff(x,t))*d for x,d in zip(xs,dq)))
    print(f'   {label:34s} Σ(∂L/∂x)δx = {s}' + ('   ⇒ สมมาตร ✓' if s==0 else '   ⇒ ไม่ใช่สมมาตร'))
    return s, Q

x1,x2,x3 = (sp.Function(f'x{i}')(t) for i in (1,2,3))
T3 = sp.Rational(1,2)*mm*sum(sp.diff(x,t)**2 for x in (x1,x2,x3))
# ข้อ 4 — ลูกโซ่ล้วน: เลื่อนทั้งสามตัวเป็นสมมาตรจริง
L4 = T3 - sp.Rational(1,2)*k*(x1-x2)**2 - sp.Rational(1,2)*K*(x2-x3)**2
s4,Q4 = check(L4, (x1,x2,x3), (1,1,1), 'ลูกโซ่ล้วน · เลื่อนทั้งสาม')
assert s4 == 0; print(f'      Q = {Q4}  ⇒ m*v1+m*v2+m*v3')
# ข้อ 5 — มีสปริงยึดพื้นที่ตัว 1: ไม่มีสมมาตรการเลื่อนเลย
L5 = (sp.Rational(1,2)*mm*(sp.diff(x1,t)**2+sp.diff(x2,t)**2)
      - sp.Rational(1,2)*k*(x1-x2)**2 - sp.Rational(1,2)*C*x1**2)
for dq, lbl in [((1,1),'ตัว 1 กับ 2 พร้อมกัน'), ((1,0),'เฉพาะตัว 1'), ((0,1),'เฉพาะตัว 2')]:
    s,_ = check(L5, (x1,x2), dq, 'ยึดพื้นที่ตัว 1 · '+lbl); assert s != 0
print('   ⇒ ปริมาณอนุรักษ์จากสมมาตรการเลื่อน = 0 ตัว')

print('\n── ข้อ 6-8 · เมตริกเป็นตัวเลข ' + '─'*31)
gm = sp.diag(-1,1,1,1)
dot = lambda a,b: (sp.Matrix(a).T*gm*sp.Matrix(b))[0]
ok('v^μ=(1,2,3,4) · v·v', dot([1,2,3,4],[1,2,3,4]), 28)
print('      บวก ⇒ ช่วงแบบอวกาศ (spacelike)')
ok('v^μ=(3,1,2,2) · v·v', dot([3,1,2,2],[3,1,2,2]), 0)
print('      ศูนย์ ⇒ ช่วงแบบแสง (null) — แสงเดินบนเส้นนี้')
ok('A^μ=(2,1,0,0) · B^μ=(3,4,0,0) · A_μB^μ', dot([2,1,0,0],[3,4,0,0]), -2)

print('\n── ข้อ 9 · ปฏิสมมาตรใน D มิติ ' + '─'*31)
import itertools
for D in (3,4,5):
    M = sp.zeros(D,D)
    for a_,b_ in itertools.combinations(range(D),2):
        s = sp.Symbol(f'f{a_}{b_}'); M[a_,b_], M[b_,a_] = s, -s
    assert M.T == -M and len(M.free_symbols) == D*(D-1)//2
    print(f'   D={D}: นับจากเมทริกซ์จริงได้ {len(M.free_symbols)} = D(D-1)/2')
