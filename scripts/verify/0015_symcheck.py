"""ตรวจทุกข้อของบทฝึก "ตรวจสมมาตรทีละตัว"

หัวใจ: ตัวตรวจ Σ(∂L/∂q_i)δq_i ต้องถูกลงมือ **ทีละอนุภาค แล้วค่อยบวก**
ผลรอบก่อนชี้ว่าอนุพันธ์รายตัวถูก แต่ขั้นที่หายไปคือการบวก
"""
import sympy as sp
t = sp.symbols('t'); m = sp.Symbol('m', positive=True)
X = [sp.Function(f'x{i}')(t) for i in range(1, 5)]
V_ = lambda *args: None

def kinetic(n):  return sp.Rational(1,2)*m*sum(sp.diff(x,t)**2 for x in X[:n])
def report(label, L, xs, dq):
    print(f'\n{label}')
    parts = []
    for i, x in enumerate(xs, 1):
        d = sp.expand(sp.diff(L, x)); parts.append(d)
        print(f'   ∂L/∂x{i} = {d}')
    s = sp.simplify(sum(p*q for p, q in zip(parts, dq)))
    Q = sp.simplify(sum(sp.diff(L, sp.diff(x,t))*q for x, q in zip(xs, dq)))
    print(f'   Σ(∂L/∂x_i)·δx_i [δ={dq}] = {s}' + ('   ⇒ สมมาตร ✓' if s == 0 else '   ⇒ ไม่ใช่'))
    if s == 0: print(f'   Q = {Q}')
    return s, Q

x1,x2,x3,x4 = X
# A · สองอนุภาค ไม่มีตัวยึด (k=4)
LA = kinetic(2) - sp.Rational(1,2)*4*(x1-x2)**2
sA,QA = report('A · N=2 ไม่มีตัวยึด (k=4)', LA, [x1,x2], (1,1)); assert sA == 0
# B · สองอนุภาค ยึดที่ x1 (k=4, C=6) — ข้อที่พลาดรอบก่อน
LB = LA - sp.Rational(1,2)*6*x1**2
sB,_ = report('B · N=2 ยึดที่ x1 (k=4, C=6)', LB, [x1,x2], (1,1))
assert sB == -6*x1
for dq in [(1,0),(0,1)]:
    s,_ = report(f'   B · ลองชุด δ={dq}', LB, [x1,x2], dq); assert s != 0
print('   ⇒ ไม่มีชุดไหนผ่าน ⇒ ปริมาณอนุรักษ์จากการเลื่อน = 0 ตัว')
# C · สามอนุภาคลูกโซ่ (k=4, K=2) — ข้อที่ขั้น 1 พลาด
LC = kinetic(3) - sp.Rational(1,2)*4*(x1-x2)**2 - sp.Rational(1,2)*2*(x2-x3)**2
sC,QC = report('C · N=3 ลูกโซ่ (k=4, K=2)', LC, [x1,x2,x3], (1,1,1)); assert sC == 0
# D · สามอนุภาค ยึดที่ x3 (k=4, C=6)
LD = kinetic(3) - sp.Rational(1,2)*4*(x1-x2)**2 - sp.Rational(1,2)*6*x3**2
sD,_  = report('D · N=3 ยึดที่ x3 · ลองเลื่อนทั้งสาม', LD, [x1,x2,x3], (1,1,1)); assert sD != 0
sD2,QD = report('D · N=3 ยึดที่ x3 · เลื่อนเฉพาะ 1 กับ 2', LD, [x1,x2,x3], (1,1,0)); assert sD2 == 0
# E · อนุภาคเดียวใน 2 มิติ · V ขึ้นกับ y เท่านั้น
xs, ys = sp.Function('x')(t), sp.Function('y')(t)
Vy = sp.Function('V')(ys)
LE = sp.Rational(1,2)*m*(sp.diff(xs,t)**2 + sp.diff(ys,t)**2) - Vy
sE,QE = report('E · 2 มิติ · V(y) เท่านั้น · δx=1, δy=0', LE, [xs,ys], (1,0)); assert sE == 0
sE2,_ = report('E · ลองเลื่อนแนว y แทน (δx=0, δy=1)', LE, [xs,ys], (0,1))
print('      ⇒ ขึ้นกับ V จึงไม่เป็นศูนย์ทั่วไป ⇒ ทิศ y ไม่ใช่สมมาตร')
# F · สามมิติ · V ขึ้นกับ x และ y (ไม่ขึ้นกับ z) ⇒ นับได้กี่ตัว
zx,zy,zz = (sp.Function(c)(t) for c in 'XYZ')
LF = sp.Rational(1,2)*m*sum(sp.diff(c,t)**2 for c in (zx,zy,zz)) - sp.Function('V')(zx,zy)
n = 0
for lbl, dq in [('x',(1,0,0)), ('y',(0,1,0)), ('z',(0,0,1))]:
    s = sp.simplify(sum(sp.diff(LF,c)*q for c,q in zip((zx,zy,zz), dq)))
    good = (s == 0); n += good
    print(f'   F · เลื่อนแนว {lbl}: Σ = {s}' + ('   ✓' if good else ''))
print(f'   ⇒ ปริมาณอนุรักษ์จากการเลื่อน = {n} ตัว'); assert n == 1
# G · สี่อนุภาคลูกโซ่ ยึดปลายหนึ่งข้าง
LG = (kinetic(4) - sp.Rational(1,2)*4*((x1-x2)**2+(x2-x3)**2+(x3-x4)**2)
      - sp.Rational(1,2)*6*x1**2)
cnt = 0
import itertools
for dq in itertools.product((0,1), repeat=4):
    if dq == (0,0,0,0): continue
    if sp.simplify(sum(sp.diff(LG,x)*q for x,q in zip(X[:4], dq))) == 0: cnt += 1
print(f'\nG · N=4 ลูกโซ่ ยึดที่ x1 · ชุด δ ที่ผ่านทั้งหมด (ไล่ครบ 15 ชุด): {cnt}')
assert cnt == 0
print('   ⇒ 0 ตัว')
