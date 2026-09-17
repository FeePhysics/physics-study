"""ตรวจทุกข้ออ้างของบทที่ 6 · แฮมิลโทเนียน

กฎบ้าน: ห้ามเชื่อความจำ — ทุกสูตรต้องให้ sympy ยืนยัน
และทุกกฎที่เขียนลงบทต้องมีข้อฝึกคู่กัน (บทเรียนจาก 0014)
"""
import sympy as sp
t = sp.symbols('t')
ok = lambda n,g,w: print(f"{'✓' if sp.simplify(g-w)==0 else '✗'} {n:44s} = {sp.simplify(g)}")

# ── เครื่องมือกลาง ─────────────────────────────────────────
def legendre(L, qs, vs):
    """H = Σ p·q̇ − L  โดยกำจัด q̇ ให้หมด · คืน (H(q,p), สูตร p, สูตร q̇)"""
    ps = sp.symbols(' '.join(f'p{i}' for i in range(len(qs))))
    ps = (ps,) if not isinstance(ps, tuple) else ps
    pdef = [sp.diff(L, v) for v in vs]
    sol = sp.solve([sp.Eq(p, d) for p, d in zip(ps, pdef)], list(vs), dict=True)[0]
    H = sp.simplify((sum(p*sol[v] for p, v in zip(ps, vs)) - L).subs(sol))
    return sp.simplify(H), pdef, sol, ps

def pb(A, B, qs, ps):
    return sp.simplify(sum(sp.diff(A,q)*sp.diff(B,p) - sp.diff(A,p)*sp.diff(B,q)
                           for q, p in zip(qs, ps)))

print('── ① เลอฌ็องดร์ · ออสซิลเลเตอร์ m=2, k=8 ' + '─'*18)
y, v = sp.symbols('y v')
L1 = sp.Rational(1,2)*2*v**2 - sp.Rational(1,2)*8*y**2
H1, pdef, sol, (p,) = legendre(L1, (y,), (v,))
ok('p = ∂L/∂v', pdef[0], 2*v)
ok('ẏ เขียนกลับด้วย p', sol[v], p/2)
ok('H = p·ẏ − L', H1, p**2/4 + 4*y**2)
print('   ⇒ H เป็นฟังก์ชันของ (y, p) ล้วน ไม่มี v เหลือ:', v not in H1.free_symbols)
# H ต้องเท่ากับ T+V สำหรับระบบแบบนี้
ok('H เท่ากับ T+V ไหม', H1 - (p**2/4 + 4*y**2), 0)

print('\n── ② สมการแฮมิลตัน ' + '─'*40)
ok('ẏ = ∂H/∂p', sp.diff(H1, p), p/2)
ok('ṗ = −∂H/∂y', -sp.diff(H1, y), -8*y)
# ต้องให้สมการการเคลื่อนที่เดียวกับ EL
a = sp.Function('Y')(t); Pt = 2*sp.diff(a,t)
el  = sp.diff(sp.diff(L1,v).subs(v, sp.diff(a,t)), t) - sp.diff(L1,y).subs(y,a)
ham = sp.diff(Pt,t) - (-8*a)
ok('EL กับ แฮมิลตัน ให้สมการเดียวกัน', sp.expand(el-ham), 0)

print('\n── ③ พิกัดเชิงขั้ว · θ เป็นพิกัดวัฏจักร ' + '─'*20)
m, r, w, vr = sp.symbols('m r w vr', positive=True)
th = sp.Symbol('theta')
L3 = sp.Rational(1,2)*m*(vr**2 + r**2*w**2)
ok('p_θ = ∂L/∂θ̇', sp.diff(L3, w), m*r**2*w)
pr, pth = sp.symbols('pr pth')
H3 = sp.simplify(pr**2/(2*m) + pth**2/(2*m*r**2))
print(f'   H = {H3}')
ok('∂H/∂θ (θ ไม่โผล่ใน H เลย)', sp.diff(H3, th), 0)
ok('ṗ_θ = −∂H/∂θ', -sp.diff(H3, th), 0)
# ยืนยันด้วยเลอฌ็องดร์จริง ไม่ใช่พิมพ์ H ไว้เอง
H3b, _, _, _ = legendre(L3, (r, th), (vr, w))
print('   เลอฌ็องดร์จริงให้:', sp.simplify(H3b))

print('\n── ④ วงเล็บปัวซง ' + '─'*44)
q = sp.Symbol('q')
ok('{y, H}  (ต้องเท่า ∂H/∂p)', pb(y, H1, (y,), (p,)), p/2)
ok('{p, H}  (ต้องเท่า −∂H/∂y)', pb(p, H1, (y,), (p,)), -8*y)
ok('{q, p}', pb(q, sp.Symbol('pq'), (q,), (sp.Symbol('pq'),)), 1)
ok('{H, H}  ⇒ พลังงานอนุรักษ์', pb(H1, H1, (y,), (p,)), 0)
Hf = p**2/4 + 3*q
ok('{p, H} เมื่อ H = p²/4 + 3q', pb(p, Hf, (q,), (p,)), -3)
print('   ⇒ ไม่เป็นศูนย์ ⇒ โมเมนตัมไม่อนุรักษ์ (มีแรงคงที่ −3)')

print('\n── ⑤ {A,H} = 0 ⟺ อนุรักษ์ — ตรงกับเนอเธอร์ ' + '─'*15)
# สองอนุภาคเชื่อมสปริง: โมเมนตัมรวมอนุรักษ์ (เคยได้จากเนอเธอร์)
q1,q2,P1,P2,k = sp.symbols('q1 q2 P1 P2 k', positive=False)
Hs = (P1**2 + P2**2)/(2*sp.Symbol('m', positive=True)) + sp.Rational(1,2)*k*(q1-q2)**2
ok('{P1+P2, H} สองอนุภาคเชื่อมสปริง', pb(P1+P2, Hs, (q1,q2), (P1,P2)), 0)
ok('{P1, H} ตัวเดียวไม่อนุรักษ์', pb(P1, Hs, (q1,q2), (P1,P2)), -k*(q1-q2))
print('   ⇒ ผลเดียวกับที่ได้จากเนอเธอร์ทุกประการ · คนละภาษา ฟิสิกส์เดียวกัน')
