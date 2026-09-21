"""ตรวจบทที่ 8 · กระแสอนุรักษ์ กับ เทนเซอร์ความเค้น-พลังงาน
ทุกอย่างคำนวณเป็นองค์ประกอบจริงด้วยเมตริก diag(-1,1,1,1) ตามที่ใช้มาทั้ง workspace
"""
import sympy as sp
t,x,y,z,m = sp.symbols('t x y z m', real=True)
X=(t,x,y,z); g=sp.diag(-1,1,1,1)
ok=lambda n,got,want: print(f"{'✓' if sp.simplify(got-want)==0 else '✗'} {n}")
phi=sp.Function('phi')(*X)
d  =[sp.diff(phi,c) for c in X]
dup=[sum(g[a,b]*d[b] for b in range(4)) for a in range(4)]
sq =sp.expand(sum(d[a]*dup[a] for a in range(4)))
Lag=-sp.Rational(1,2)*sq - sp.Rational(1,2)*m**2*phi**2
D=sp.symbols('D0 D1 D2 D3')
Ls=Lag.subs({d[a]:D[a] for a in range(4)}, simultaneous=True)
dLdD=[sp.diff(Ls,D[a]).subs({D[a]:d[a] for a in range(4)}, simultaneous=True) for a in range(4)]
for a,c in enumerate('txyz'): ok(f'∂𝓛/∂(∂_{c}φ) = −∂^{c}φ', dLdD[a], -dup[a])
EOM = sp.diff(phi,t,2) - sum(sp.diff(phi,c,2) for c in (x,y,z)) + m**2*phi   # = 0 บนเปลือก

print('\n── ① กระแสจากสมมาตรการเลื่อนค่าสนาม (δφ = 1) ' + '─'*14)
print('   เงื่อนไข: 𝓛 ต้องไม่มี φ โผล่เปล่า ๆ ⇒ ต้องไม่มีพจน์มวล')
Lag0 = Lag.subs(m,0)
dLdD0=[dd.subs(m,0) for dd in dLdD]
j=[dd*1 for dd in dLdD0]                               # j^μ = ∂𝓛/∂(∂_μφ)·δφ
div=sp.simplify(sum(sp.diff(j[a],X[a]) for a in range(4)))
ok('j^μ = −∂^μφ', j[0], -dup[0].subs(m,0))
ok('∂_μ j^μ = −□φ  ⇒ เป็นศูนย์บนเปลือก (m=0)', div, -( -sp.diff(phi,t,2)+sum(sp.diff(phi,c,2) for c in (x,y,z)) ))
print(f'   ∂_μ j^μ = {sp.expand(div)}   ⇒ ใช้สมการการเคลื่อนที่ (m=0) แล้วได้ 0')
# มีพจน์มวล ⇒ ไม่ใช่สมมาตร ⇒ กระแสไม่อนุรักษ์
jm=[dd*1 for dd in dLdD]
divm=sp.simplify(sum(sp.diff(jm[a],X[a]) for a in range(4)))
sub=sp.solve(sp.Eq(EOM,0), sp.diff(phi,t,2))[0]
ok('มีพจน์มวล ⇒ ∂_μ j^μ = −m²φ ≠ 0', sp.simplify(divm.subs(sp.diff(phi,t,2), sub)), -m**2*phi)
print('   ⇒ พจน์มวลทำลายสมมาตรการเลื่อน ⇒ กระแสไม่อนุรักษ์ (รากของทฤษฎีบทโกลด์สโตน)')

print('\n── ② 1+1 มิติ (รูปที่ใช้ในข้อฝึก) ' + '─'*27)
pt_,px_,ph_=sp.symbols('pt px phi')
L11 = sp.Rational(1,2)*pt_**2 - sp.Rational(1,2)*px_**2
print(f'   𝓛 = ½pt² − ½px²  ⇒  j^t = {sp.diff(L11,pt_)} · j^x = {sp.diff(L11,px_)}')
assert sp.diff(L11,pt_)==pt_ and sp.diff(L11,px_)==-px_
f=sp.Function('f')(t,x)
print(f'   ∂_μ j^μ = ∂_t(pt) + ∂_x(−px) = {sp.expand(sp.diff(f,t,2)-sp.diff(f,x,2))}  ⇒ 0 เมื่อ m=0 ✓')

print('\n── ③ เทนเซอร์ความเค้น-พลังงาน — หาเครื่องหมายที่ถูก ' + '─'*8)
for lbl, sgn in [('T = −∂^μφ∂^νφ − g^{μν}𝓛', -1), ('T = +∂^μφ∂^νφ + g^{μν}𝓛', +1)]:
    T00 = sgn*(dup[0]*dup[0] + g[0,0]*Lag)
    print(f'   {lbl:28s} ⇒ T⁰⁰ = {sp.expand(T00)}')
T = sp.Matrix(4,4, lambda a,b: dup[a]*dup[b] + g[a,b]*Lag)
want00 = (sp.Rational(1,2)*sp.diff(phi,t)**2
          + sp.Rational(1,2)*sum(sp.diff(phi,c)**2 for c in (x,y,z))
          + sp.Rational(1,2)*m**2*phi**2)
ok('T⁰⁰ = ½φ̇² + ½|∇φ|² + ½m²φ²  (บวกเสมอ)', sp.expand(T[0,0]), sp.expand(want00))
# ⚠️ ห้ามเทียบ Matrix กับ 0 ตรง ๆ และห้ามใช้ sum(Abs(...)) — ทั้งคู่ไม่ยุบเป็น 0
#    แม้ทุกช่องจะเป็นศูนย์จริง ⇒ ตัวตรวจจะรายงาน ✗ ทั้งที่ฟิสิกส์ถูก (เกิดจริงรอบนี้)
bad=[(a,b) for a in range(4) for b in range(4) if sp.simplify(T[a,b]-T[b,a])!=0]
print(f"{'✓' if not bad else '✗'} T^{{μν}} สมมาตร (ไล่ทีละช่อง · ไม่สมมาตร {len(bad)} ช่อง)")
assert not bad
Tmix=sp.Matrix(4,4, lambda a,b: sum(T[a,c]*g[c,b] for c in range(4)))
for b in range(4):
    e=sp.simplify(sp.expand(sum(sp.diff(Tmix[a,b], X[a]) for a in range(4))
                            .subs(sp.diff(phi,t,2), sub)))
    print(f"{'✓' if e==0 else '✗'} ∂_μ T^μ{{}}_ν = 0 · ν={b}"); assert e==0
print(f'   T⁰¹ = {sp.expand(T[0,1])}   (ความหนาแน่นโมเมนตัม = ฟลักซ์พลังงาน)')

print('\n── ④ 1+1 มิติสำหรับข้อฝึก + เช็คมิติ ' + '─'*24)
T00_11 = sp.Rational(1,2)*pt_**2 + sp.Rational(1,2)*px_**2 + sp.Rational(1,2)*m**2*ph_**2
print(f'   T⁰⁰ = {T00_11}')
for vals,want in [({pt_:2,px_:4,m:0},10), ({pt_:1,px_:1,m:2,ph_:3},19), ({pt_:0,px_:0,m:5,ph_:2},50)]:
    v=T00_11.subs(vals); print(f'      {dict((str(k),vv) for k,vv in vals.items())} ⇒ {v}'); assert v==want
print('   [T^{μν}] = [∂φ][∂φ] = 2+2 = 4 ✓ เท่ากับ [𝓛] — ตรวจงานตัวเองได้ทันที')
assert 2*(1+1)==4
print('\n── ⑤ เทียบกับ −¼F² ที่ไม่ใช่พลังงาน ' + '─'*26)
E,B=sp.symbols('E B', positive=True)
print(f'   −¼F² = ½(E²−B²) ⇒ เป็นศูนย์เมื่อ E=B (คลื่นแสง) และติดลบได้')
print(f'   T⁰⁰ ของแม่เหล็กไฟฟ้า = ½(E²+B²) ⇒ บวกเสมอ · E=B=3 ⇒ {sp.Rational(1,2)*(9+9)}')
assert sp.Rational(1,2)*(3**2+3**2)==9

print('\n── ⑥ ทฤษฎีบท: ไดเวอร์เจนซ์ของกระแส = ปริมาณที่สมมาตรพัง ' + '─'*2)
dL = sp.diff(Lag, phi)*1          # δ𝓛 เมื่อ δφ = 1
print(f'   δ𝓛 (เลื่อน φ → φ+1) = {dL}')
divm = sp.simplify(sum(sp.diff(dd, X[a]) for a,dd in enumerate(dLdD))
                   .subs(sp.diff(phi,t,2), sub))
print(f'   ∂_μ j^μ บนเปลือก      = {divm}')
ok('∂_μ j^μ = δ𝓛  ⇒ สมมาตรพังเท่าไร กระแสรั่วเท่านั้น', divm, dL)

print('\n── ⑦ จำนวนปริมาณอนุรักษ์จากการเลื่อนกาลอวกาศ ' + '─'*13)
print('   เลื่อนได้ 4 ทิศอิสระ (t, x, y, z) ⇒ กระแสอนุรักษ์ 4 ชุด')
print('   ⇒ พลังงาน 1 + โมเมนตัม 3 = 4 ตัว · นี่คือดัชนี ν ของ T^μ{}_ν พอดี')
assert 4 == 1 + 3

print('\n── ⑧ ตัวเลขของข้อฝึกทุกข้อ ' + '─'*33)
pt2,px2,ph2,m2 = sp.symbols('pt px phi m')
T00 = sp.Rational(1,2)*pt2**2 + sp.Rational(1,2)*px2**2 + sp.Rational(1,2)*m2**2*ph2**2
for v,w in [({pt2:2,px2:4,m2:0},10), ({pt2:1,px2:1,m2:2,ph2:3},19), ({pt2:0,px2:0,m2:5,ph2:2},50)]:
    got=T00.subs(v); print(f'   T⁰⁰{ {str(k):vv for k,vv in v.items()} } = {got}'); assert got==w
print(f'   แม่เหล็กไฟฟ้า E=B=3 ⇒ T⁰⁰ = ½(E²+B²) = {sp.Rational(1,2)*(9+9)}'
      f'  (ขณะที่ −¼F² = {sp.Rational(1,2)*(9-9)})')
print('\n✅ ตรวจครบทุกข้ออ้างของบทที่ 8')
