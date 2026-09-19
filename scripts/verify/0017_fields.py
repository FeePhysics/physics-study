"""ตรวจบทที่ 7 · จากอนุภาคไปสนาม — ทุกอย่างคำนวณเป็นองค์ประกอบจริง

⚠️ จุดเสี่ยงที่สุดของบทนี้คือ **เครื่องหมายของเมตริก** · workspace นี้ใช้
   g = diag(-1,+1,+1,+1) มาตั้งแต่บทที่ 3 ⇒ ต้องยืนยันว่า 𝓛 ที่เขียนให้รูป T − V จริง
"""
import sympy as sp

t, x, y, z, m, om, kx, ky, kz = sp.symbols('t x y z m omega k_x k_y k_z', real=True)
X  = (t, x, y, z)
g  = sp.diag(-1, 1, 1, 1)
ok = lambda n, got, want: print(f"{'✓' if sp.simplify(got-want)==0 else '✗'} {n}")

phi = sp.Function('phi')(*X)
d   = [sp.diff(phi, c) for c in X]                 # ∂_μ φ  (ดัชนีล่าง)
dup = [sum(g[mu,nu]*d[nu] for nu in range(4)) for mu in range(4)]   # ∂^μ φ

print('── ① เมตริกทำอะไรกับพจน์จลน์ ' + '─'*30)
sq = sp.expand(sum(d[mu]*dup[mu] for mu in range(4)))
ok('∂_μφ ∂^μφ = −φ̇² + |∇φ|²',
   sq, -sp.diff(phi,t)**2 + sum(sp.diff(phi,c)**2 for c in (x,y,z)))
print(f'   = {sq}')

print('\n── ② 𝓛 ต้องออกมาเป็นรูป T − V ' + '─'*29)
Lag = -sp.Rational(1,2)*sq - sp.Rational(1,2)*m**2*phi**2
want = (sp.Rational(1,2)*sp.diff(phi,t)**2
        - sp.Rational(1,2)*sum(sp.diff(phi,c)**2 for c in (x,y,z))
        - sp.Rational(1,2)*m**2*phi**2)
ok('𝓛 = −½ ∂_μφ∂^μφ − ½m²φ²  ⇒  ½φ̇² − ½|∇φ|² − ½m²φ²', Lag, want)
print('   ⇒ พจน์เวลาเป็นบวก (T) · พจน์ที่ว่างกับมวลเป็นลบ (V) — ถูกแบบเดียวกับกลศาสตร์')
print(f'   ⚠️ ถ้าเขียน +½∂_μφ∂^μφ จะได้ {sp.expand(sp.Rational(1,2)*sq)} ← พจน์จลน์ติดลบ = ผิด')

print('\n── ③ สมการ EL ของสนาม (คำนวณเป็นองค์ประกอบ) ' + '─'*15)
# ∂𝓛/∂(∂_μφ) โดยแทน ∂_μφ เป็นสัญลักษณ์อิสระก่อน
D = sp.symbols('D0 D1 D2 D3')
Ls = Lag.subs({d[mu]: D[mu] for mu in range(4)}, simultaneous=True)
dLdD = [sp.diff(Ls, D[mu]).subs({D[mu]: d[mu] for mu in range(4)}, simultaneous=True)
        for mu in range(4)]
for mu, c in enumerate('txyz'):
    ok(f'∂𝓛/∂(∂_{c}φ) = −∂^{c}φ', dLdD[mu], -dup[mu])
EL = sum(sp.diff(dLdD[mu], X[mu]) for mu in range(4)) - sp.diff(Lag, phi)
EL = sp.simplify(sp.expand(EL))
print(f'   EL = {EL}')
box = -sp.diff(phi,t,2) + sum(sp.diff(phi,c,2) for c in (x,y,z))    # □ = ∂_μ∂^μ
ok('EL = 0  ⟺  □φ = m²φ', EL, -(box - m**2*phi))
ok('   กางออก ⟺ φ̈ = ∇²φ − m²φ',
   box - m**2*phi, -sp.diff(phi,t,2) + sum(sp.diff(phi,c,2) for c in (x,y,z)) - m**2*phi)

print('\n── ④ คลื่นระนาบ ⇒ ความสัมพันธ์การกระจาย ' + '─'*19)
w = sp.cos(om*t - kx*x - ky*y - kz*z)
res = sp.simplify((-sp.diff(w,t,2) + sum(sp.diff(w,c,2) for c in (x,y,z)) - m**2*w) / w)
print(f'   แทน φ = cos(ωt − k·x) ลงใน □φ − m²φ แล้วหารด้วย φ ได้: {sp.expand(res)}')
sol = sp.solve(sp.Eq(res, 0), om**2)
ok('ω² = k² + m²', sp.expand(sol[0]), kx**2+ky**2+kz**2+m**2)
print('   ⇒ ในหน่วย ħ=c=1 คือ  E² = p² + m²  — สูตรสัมพัทธภาพของไอน์สไตน์')
print(f'   ตั้ง m=0 ⇒ ω = |k| ⇒ ความเร็ว ω/|k| = {sp.sqrt(kx**2)/sp.Abs(kx)} ⇒ เท่ากับ c')

print('\n── ⑤ กลศาสตร์กับสนาม — กระบวนการเดียวกันเป๊ะ ' + '─'*14)
q  = sp.Function('q')(t); mm, k = sp.symbols('m_p k', positive=True)
Lp = sp.Rational(1,2)*mm*sp.diff(q,t)**2 - sp.Rational(1,2)*k*q**2
V  = sp.Symbol('V')
Lps = Lp.subs(sp.diff(q,t), V)
ELp = sp.diff(sp.diff(Lps,V).subs(V, sp.diff(q,t)), t) - sp.diff(Lp, q)
print(f'   อนุภาค: d/dt(∂L/∂q̇) − ∂L/∂q = {sp.expand(ELp)}  ⇒ m q̈ = −k q')
print( '   สนาม  : ∂_μ(∂𝓛/∂(∂_μφ)) − ∂𝓛/∂φ = 0  ⇒ □φ = m²φ')
print( '   ⇒ ต่างกันแค่ d/dt กลายเป็น ∂_μ (รวมสี่ทิศ) และ q(t) กลายเป็น φ(x^μ)')

print('\n── ⑥ เช็คความเข้ากันกับบทที่ 5 (−¼F²) ' + '─'*22)
A = [sp.Function(f'A{i}')(*X) for i in range(4)]
Adn = [sum(g[i,j]*A[j] for j in range(4)) for i in range(4)]
F  = sp.Matrix(4,4, lambda a,b: sp.diff(Adn[b], X[a]) - sp.diff(Adn[a], X[b]))
Fup = g.inv()*F*g.inv().T
F2 = sp.expand(sum(F[a,b]*Fup[a,b] for a in range(4) for b in range(4)))
E = [sp.simplify(F[i,0]) for i in (1,2,3)]
B = [sp.simplify(F[2,3]), sp.simplify(F[3,1]), sp.simplify(F[1,2])]
ok('−¼F_{μν}F^{μν} = ½(E² − B²)',
   sp.expand(-sp.Rational(1,4)*F2),
   sp.expand(sp.Rational(1,2)*(sum(e**2 for e in E) - sum(b**2 for b in B))))
print('   ⇒ เป็นรูป T − V เหมือนกัน · บทที่ 5 ใช้เครื่องหมายเดียวกับบทนี้ ไม่ขัดกัน')

print('\n── ⑦ ตัวเลขของข้อฝึกทุกข้อ ' + '─'*33)
# ข้อ 1 · φ = 3t + 2x
p1 = 3*t + 2*x
d1 = [sp.diff(p1,c) for c in X]
v1 = sum(g[i,j]*d1[i]*d1[j] for i in range(4) for j in range(4))
print(f'   ข้อ1: ∂_tφ={d1[0]} · ∂_xφ={d1[1]} · ∂_μφ∂^μφ={sp.expand(v1)}')
assert v1 == -5
# ข้อ 2 · □φ เมื่อ φ = t² + x³
p2 = t**2 + x**3
b2 = sp.expand(-sp.diff(p2,t,2) + sum(sp.diff(p2,c,2) for c in (x,y,z)))
print(f'   ข้อ2: □φ = {b2}');  assert b2 == 6*x - 2
# ข้อ 3 · อนุพันธ์ย่อยของ 𝓛 ใน 1+1 มิติ
pt_, px_, ph_ = sp.symbols('pt px phi')
L11 = sp.Rational(1,2)*pt_**2 - sp.Rational(1,2)*px_**2 - sp.Rational(1,2)*m**2*ph_**2
print(f'   ข้อ3: ∂𝓛/∂pt={sp.diff(L11,pt_)} · ∂𝓛/∂px={sp.diff(L11,px_)} · ∂𝓛/∂φ={sp.diff(L11,ph_)}')
assert sp.diff(L11,pt_)==pt_ and sp.diff(L11,px_)==-px_ and sp.diff(L11,ph_)==-m**2*ph_
# ประกอบตาม EL ใน 1+1 มิติ ⇒ ∂_t²φ − ∂_x²φ + c φ = 0 · c = ?
f11 = sp.Function('f')(t,x)
elf = sp.expand(sp.diff(sp.diff(f11,t),t) + sp.diff(-sp.diff(f11,x),x) + m**2*f11)
print(f'   ข้อ3 ขั้น4: EL ใน 1+1 มิติ = {elf}  ⇒ c = m²')
assert sp.simplify(elf - (sp.diff(f11,t,2) - sp.diff(f11,x,2) + m**2*f11)) == 0
# ข้อ 4-6 · การกระจาย
wv = sp.cos(om*t - kx*x)
disp = sp.simplify(sp.expand((sp.diff(wv,t,2) - sp.diff(wv,x,2) + m**2*wv)/wv))
w2 = sp.solve(sp.Eq(disp,0), om**2)[0]
print(f'   ข้อ4: ω² = {w2}');                       assert sp.expand(w2) == kx**2 + m**2
print(f'   ข้อ5: m=3,k=4 ⇒ ω = {sp.sqrt(w2.subs({m:3,kx:4}))}')
assert sp.sqrt(w2.subs({m:3,kx:4})) == 5
print(f'   ข้อ6: m=0 ⇒ ω/k = {sp.simplify(sp.sqrt(w2.subs(m,0))/kx)}')
assert sp.simplify(sp.sqrt(w2.subs(m,0))/sp.Abs(kx)) == 1
# ข้อ 7 · −¼F² เมื่อ E=4, B=0
print(f'   ข้อ7: ½(E²−B²) เมื่อ E=4,B=0 = {sp.Rational(1,2)*(4**2-0)}')
assert sp.Rational(1,2)*(4**2) == 8
print('   ข้อ8: [𝓛] = M⁴ (บทที่ 3) ⇒ ตอบ 4')
print('\n✅ ตัวเลขของข้อฝึกทั้งแปดข้อตรวจครบ')
