"""ตรวจทุกข้อของชุดทดสอบรอบสอง — คำนวณจากศูนย์ ไม่ได้เชื่อเลขที่พิมพ์ไว้

กฎของโฟลเดอร์: ข้ออ้างทุกข้อในบทเรียนต้องวัดก่อนเขียน
และ **ห้ามออกข้อที่ใช้กฎซึ่งไม่เคยสอน** (พลาดมาแล้วสามครั้ง)
"""
import sympy as sp

t = sp.symbols('t')
ok = lambda name, got, want: print(
    f"{'✓' if sp.simplify(got - want) == 0 else '✗'} {name:44s} = {got}")

print('── ข้อ 1 · EL บนระบบใหม่ (รอกชั่งสองข้าง) ' + '─'*20)
# L = 2 xdot^2 + 20 x   (m1=3, m2=1, g=10 ⇒ ½(m1+m2)=2 · (m1-m2)g=20)
v, a = sp.symbols('v a')
L1 = 2*v**2 + 20*sp.Symbol('x')
p1 = sp.diff(L1, v)                      # ∂L/∂v
f1 = sp.diff(L1, sp.Symbol('x'))         # ∂L/∂x
ok('∂L/∂v', p1, 4*v)
ok('∂L/∂x', f1, 20)
acc = sp.solve(sp.Eq(4*a, f1), a)[0]     # d/dt(4v) = 4a = ∂L/∂x
ok('ความเร่ง', acc, 5)
# เทียบกับสูตรรอกมาตรฐาน (m1-m2)g/(m1+m2) — คนละทางกับที่ไล่ข้างบน
ok('สูตรรอก (m1-m2)g/(m1+m2)', sp.Rational(3-1,3+1)*10, 5)

print('\n── ข้อ 2 · ประกอบ S เมื่อ m ≠ 1 ' + '─'*28)
m2_, y2 = 4, 3*t - t**2
S2 = sp.integrate(sp.Rational(1,2)*m2_*sp.diff(y2, t)**2, (t, 0, 1))
ok('S = ∫ ½m ẏ² dt  บน [0,1]', S2, sp.Rational(26,3))
print(f'   (= {float(S2):.4f} ⇒ ตอบ 26/3 หรือ 8.667 ก็ได้)')

print('\n── ข้อ 3 · สัมประสิทธิ์ ε เมื่อ m ≠ 1 และมี V ' + '─'*13)
m3, g3, y0, eta = 3, 10, t**2, t*(1-t)
eps = sp.Symbol('varepsilon')
# กางแอ็กชันจริง ไม่ใช้สูตรลัด
y = y0 + eps*eta
S3 = sp.integrate(sp.expand(sp.Rational(1,2)*m3*sp.diff(y,t)**2 - m3*g3*y), (t,0,1))
coef = sp.expand(S3).coeff(eps, 1)
i1 = sp.integrate(m3*sp.diff(y0,t)*sp.diff(eta,t), (t,0,1))
i2 = sp.integrate(m3*g3*eta, (t,0,1))
ok('η̇', sp.diff(eta,t), 1-2*t)
ok('∫ m ẏ0 η̇ dt', i1, -1)
ok('∫ m g η dt', i2, 5)
ok('สัมประสิทธิ์ ε (กางแอ็กชันจริง)', coef, -6)
ok('   สูตรลัดให้ค่าเดียวกัน', i1 - i2, coef)
print(f'   ≠ 0 ⇒ y0 = t² ไม่ใช่เส้นทางจริง')
# ⚠ ถ้าตกตัว m ในพจน์แรกจะได้เท่าไร — ต้องต่างจากคำตอบจริง ไม่งั้นข้อนี้ไม่ได้วัดอะไร
bad = sp.integrate(sp.diff(y0,t)*sp.diff(eta,t), (t,0,1)) - i2
print(f'   สูตรที่ตกตัว m ให้ {bad}  ← ต่างจาก {coef} ⇒ ข้อนี้วัดตัว m ได้จริง')
assert bad != coef

print('\n── ข้อ 4 · ปริมาณอนุรักษ์เมื่อ δq ไม่เท่ากันทุกตัว ' + '─'*8)
mm, k, K = sp.symbols('m k K', positive=True)
x1,x2,x3 = (sp.Function(f'x{i}')(t) for i in (1,2,3))
L4 = (sp.Rational(1,2)*mm*sum(sp.diff(x,t)**2 for x in (x1,x2,x3))
      - sp.Rational(1,2)*k*(x1-x2)**2 - sp.Rational(1,2)*K*x3**2)
# สมมาตร: เลื่อนเฉพาะตัว 1 กับ 2  ⇒ δx1 = δx2 = 1 · δx3 = 0
dq = {x1: 1, x2: 1, x3: 0}
test = sum(sp.diff(L4, x)*d for x, d in dq.items())
ok('Σ (∂L/∂x_i)·δx_i  (ต้องเป็น 0 จึงมีสมมาตร)', sp.simplify(test), 0)
Q = sum(sp.diff(L4, sp.diff(x,t))*d for x, d in dq.items())
print(f'   Q = {sp.simplify(Q)}   ⇒ m*v1 + m*v2')
# ยืนยันว่าคงที่จริงด้วยสมการการเคลื่อนที่
eom = {sp.diff(x,t,2): sp.solve(sp.Eq(sp.diff(sp.diff(L4,sp.diff(x,t)),t),
                                      sp.diff(L4,x)), sp.diff(x,t,2))[0]
       for x in (x1,x2,x3)}
ok('dQ/dt เมื่อใช้สมการการเคลื่อนที่', sp.simplify(sp.diff(Q,t).subs(eom)), 0)
# และตัวลวงหลัก "บวกทั้งสามตัว" ต้องไม่คงที่ ไม่งั้นข้อนี้ไม่ได้วัดอะไร
Qall = sum(sp.diff(L4, sp.diff(x,t)) for x in (x1,x2,x3))
dall = sp.simplify(sp.diff(Qall,t).subs(eom))
print(f'   ตัวลวง m(v1+v2+v3): dQ/dt = {dall}  ← ไม่เป็นศูนย์ ⇒ ข้อนี้วัดได้จริง')
assert dall != 0

print('\n── ข้อ 5 · หดดัชนีด้วยเมตริก diag(-1,1,1,1) ' + '─'*14)
gm = sp.diag(-1, 1, 1, 1)                 # ธรรมเนียมเดียวกับบทที่ 3 เป๊ะ
vup = sp.Matrix([5, 1, 2, 2])
vdn = gm*vup
print(f'   v^μ = {list(vup)}  ⇒  v_μ = {list(vdn)}')
ok('v^μ v_μ', (vup.T*gm*vup)[0], -16)
print('   ติดลบ ⇒ ช่วงแบบเวลา (timelike) ตามตาราง §③')

print('\n── ข้อ 8 · องค์ประกอบอิสระของ F_{μν} ใน 4 มิติ ' + '─'*12)
import itertools
F = sp.zeros(4,4)
syms = {}
for a_, b_ in itertools.combinations(range(4), 2):
    s = sp.Symbol(f'F{a_}{b_}')
    syms[(a_,b_)] = s
    F[a_,b_], F[b_,a_] = s, -s          # ปฏิสมมาตร + แนวทแยงศูนย์
print(f'   นับจากการสร้างจริง: {len(set(F.free_symbols))} ตัว')
assert len(F.free_symbols) == 6 and F.T == -F
ok('จำนวนองค์ประกอบอิสระ', len(F.free_symbols), 6)
print('   = E สามตัว + B สามตัว ตรงกับที่บทที่ 5 แจกแจง (F10=Ex · F12=Bz)')

print('\n── ตรวจ "ข้ออ้างในเฉลย" ที่ประตูจับไม่ได้ ' + '─'*18)
yy = 3*t - t**2
print('ข้อ 2 · ลืมคูณ ½m ทั้งก้อน =', sp.integrate(sp.diff(yy,t)**2, (t,0,1)),
      ' · เผลอใช้ m=1 =', sp.integrate(sp.Rational(1,2)*sp.diff(yy,t)**2, (t,0,1)))
assert sp.integrate(sp.diff(yy,t)**2,(t,0,1)) == sp.Rational(13,3)
assert sp.integrate(sp.Rational(1,2)*sp.diff(yy,t)**2,(t,0,1)) == sp.Rational(13,6)

# ข้อ 4 · "มีปริมาณอนุรักษ์ตัวเดียว" เป็นเท็จ — พลังงานก็อนุรักษ์ ⇒ คำถามต้องระบุสมมาตร
H = sum(sp.diff(L4, sp.diff(x,t))*sp.diff(x,t) for x in (x1,x2,x3)) - L4
ok('ข้อ 4 · dH/dt (พลังงานก็อนุรักษ์)', sp.simplify(sp.diff(H,t).subs(eom)), 0)
for lbl, d in [('เลื่อนเฉพาะตัว 1', {x1:1,x2:0,x3:0}),
               ('เลื่อนเฉพาะตัว 3', {x1:0,x2:0,x3:1}),
               ('เลื่อนทั้งสามตัว',  {x1:1,x2:1,x3:1})]:
    val = sp.simplify(sum(sp.diff(L4,x)*v for x,v in d.items()))
    print(f'   Σ(∂L/∂x_i)δx_i [{lbl}] = {val}  ← ไม่เป็นศูนย์ ⇒ ไม่ใช่สมมาตร')
    assert val != 0
print('   ⇒ ชุด δq ที่ใช้ได้มีแบบเดียวจริง: เลื่อนตัว 1 กับ 2 พร้อมกัน')
