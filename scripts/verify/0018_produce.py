"""ตรวจบทซ่อม "ผลิตเอง ไม่ใช่รับมา"

สามข้อที่พลาดมีรูปร่างเดียวกัน: เป็นข้อเท็จจริงที่ถูก **ให้ไว้** ซ้ำ ๆ
แต่ไม่เคยถูก **ถามเอา** — วัดได้ว่า [𝓛] = M⁴ ให้ไว้ 6 ครั้ง ถามเอา 0 ครั้ง
"""
import sympy as sp, itertools

t,x,y,z = sp.symbols('t x y z', real=True)
X = (t,x,y,z); g = sp.diag(-1,1,1,1)
ok = lambda n,got,want: print(f"{'✓' if sp.simplify(got-want)==0 else '✗'} {n:46s} = {sp.simplify(got)}")

print('── ① เครื่องหมายของพจน์เวลา vs พจน์ที่ว่าง ' + '─'*18)
pt,px,ph,m = sp.symbols('pt px phi m')
for lbl, L in [('𝓛 = 3/2 pt² − 5/2 px² − 2φ²',
                sp.Rational(3,2)*pt**2 - sp.Rational(5,2)*px**2 - 2*ph**2),
               ('𝓛 = 2pt² − 7px²', 2*pt**2 - 7*px**2)]:
    print(f'   {lbl}')
    print(f'      ∂/∂pt = {sp.diff(L,pt)} · ∂/∂px = {sp.diff(L,px)}'
          + (f' · ∂/∂φ = {sp.diff(L,ph)}' if ph in L.free_symbols else ''))
assert sp.diff(sp.Rational(3,2)*pt**2,pt) == 3*pt
assert sp.diff(-sp.Rational(5,2)*px**2,px) == -5*px
assert sp.diff(-2*ph**2,ph) == -4*ph
assert sp.diff(2*pt**2,pt) == 4*pt and sp.diff(-7*px**2,px) == -14*px
print('   ⚠️ สัมประสิทธิ์ต่างกันโดยตั้งใจ ⇒ ลอกรูปจากขั้นก่อนแล้วได้เลขผิดทันที')

print('\n── ② ที่มาของ −¼F² = ½(E²−B²) — กางจริงทั้ง 16 ช่อง ' + '─'*8)
E = sp.symbols('E1 E2 E3'); B = sp.symbols('B1 B2 B3')
F = sp.zeros(4,4)
for i in range(3):
    F[0,i+1], F[i+1,0] = E[i], -E[i]          # F_{0i} = E_i
F[1,2], F[2,1] =  B[2], -B[2]
F[2,3], F[3,2] =  B[0], -B[0]
F[3,1], F[1,3] =  B[1], -B[1]
assert F.T == -F and len(F.free_symbols) == 6
print(f'   F ปฏิสมมาตร ✓ · องค์ประกอบอิสระ {len(F.free_symbols)} ตัว = E สาม + B สาม')
Fup = g.inv()*F*g.inv().T
print(f'   ยกดัชนีแล้ว: F^0i / F_0i = {sp.simplify(Fup[0,1]/F[0,1])}'
      f'   ·   F^12 / F_12 = {sp.simplify(Fup[1,2]/F[1,2])}')
assert sp.simplify(Fup[0,1]/F[0,1]) == -1 and sp.simplify(Fup[1,2]/F[1,2]) == 1
F2 = sp.expand(sum(F[a,b]*Fup[a,b] for a in range(4) for b in range(4)))
E2 = sum(e**2 for e in E); B2 = sum(b**2 for b in B)
ok('F_{μν}F^{μν} = 2(B² − E²)', F2, 2*(B2 - E2))
ok('−¼F² = ½(E² − B²)', -sp.Rational(1,4)*F2, sp.Rational(1,2)*(E2 - B2))
for Ev,Bv,want in [(3,1,4),(2,4,-6),(4,0,8),(0,5,sp.Rational(-25,2))]:
    v = (-sp.Rational(1,4)*F2).subs({E[0]:Ev,E[1]:0,E[2]:0,B[0]:Bv,B[1]:0,B[2]:0})
    print(f'   E={Ev}, B={Bv} ⇒ F² = {F2.subs({E[0]:Ev,E[1]:0,E[2]:0,B[0]:Bv,B[1]:0,B[2]:0})}'
          f' · −¼F² = {v}');  assert v == want

print('\n── ③ บันไดมิติเชิงมวล — ผลิตทีละขั้น ' + '─'*24)
S, d4x, Lg, phi_, dmu = sp.symbols('S d4x Lg phi dmu')
print('   [S]    = 0   (เพราะสิ่งที่โผล่จริงคือ e^{iS} และ ħ=1)')
print('   [d⁴x]  = -4  (ความยาว = 1/มวล ยกกำลังสี่)')
print('   [𝓛]    = 4   ← ต้องหักล้าง [d⁴x] ให้ [S] = 0')
print('   [φ]    = 1   (จาก m²φ² ต้องได้ M⁴)')
print('   [∂]    = 1   (∂ ~ 1/ความยาว = มวล)')
assert 0 - (-4) == 4
# ตรวจด้วยพจน์จลน์: [(∂φ)²] ต้องได้ 4 เช่นกัน
assert 2*(1+1) == 4
print('   ตรวจไขว้: [(∂φ)²] = 2([∂]+[φ]) = 2(1+1) = 4 ✓ ตรงกับ [m²φ²] = 2+2 = 4 ✓')

def coupling(term_dims):
    """[g] = 4 − ผลรวมมิติของสิ่งที่เหลือในพจน์"""
    return 4 - sum(term_dims)
print(f'\n   g·φ⁴          ⇒ [g] = {coupling([1]*4)}')
print(f'   g·φ⁶          ⇒ [g] = {coupling([1]*6)}')
print(f'   g·φ·∂_μφ∂^μφ  ⇒ [g] = {coupling([1, 1,1, 1,1])}')
print(f'   g·φⁿ ไร้มิติ  ⇒ n = {sp.solve(sp.Eq(4 - sp.Symbol("n"), 0), sp.Symbol("n"))[0]}')
assert coupling([1]*4)==0 and coupling([1]*6)==-2 and coupling([1,1,1,1,1])==-1

print('\n── ④ ยืนยันว่ากับดักของแต่ละข้อให้คำตอบคนละค่ากับเฉลย ' + '─'*6)
traps = [('ลอกเครื่องหมายจากขั้นก่อน', '5*px', '-5*px'),
         ('ตอบมิติของ S แทนของ 𝓛', 0, 4),
         ('ลืมว่ายกดัชนีเวลาแล้วติดลบ', 2*(E2+B2), F2)]
for lbl, bad, good in traps:
    same = sp.simplify(sp.sympify(bad) - sp.sympify(good)) == 0
    print(f'   {"✗ ชนกัน!" if same else "✓ ต่างกันจริง"}  {lbl}')
    assert not same
