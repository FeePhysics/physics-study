"""ตรวจบทฝึก "ประกอบมิติ" — ข้อ 6 ที่พลาดเป็นพจน์แรกในประวัติที่มีอนุพันธ์

วัดแล้ว: [∂] ให้ไว้ 4 ครั้ง ถามเอา 0 ครั้ง
และข้อ [g] ทุกข้อก่อนหน้า (λφ⁴ · gφ³ · gφ⁶) ไม่มีอนุพันธ์เลยสักข้อ
⇒ การประกอบ [(∂φ)²] = 2([∂]+[φ]) ไม่เคยถูกใช้จริง
"""
import sympy as sp
D_PHI, D_D, D_L = 1, 1, 4          # [φ] = [∂] = M · [𝓛] = M⁴
def dim(*factors):  return sum(factors)
def coupling(rest): return D_L - rest
ok = lambda n,g,w: print(f"{'✓' if g==w else '✗'} {n:52s} = {g}")

print('── ฐานที่ต้องผลิตเองให้ได้ ' + '─'*33)
ok('[∂]', D_D, 1)
ok('[φ]', D_PHI, 1)
ok('[∂_μφ ∂^μφ]  = 2([∂]+[φ])', dim(D_D,D_PHI,D_D,D_PHI), 4)
print('   ⇒ พจน์จลน์มีมิติ M⁴ พอดี ⇒ สัมประสิทธิ์ของมันไร้มิติ (นั่นคือ ½)')

print('\n── หา [g] ของแต่ละพจน์ ' + '─'*37)
terms = {
 'g φ ∂_μφ∂^μφ      (ข้อที่พลาด)': dim(D_PHI, D_D,D_PHI, D_D,D_PHI),
 'g φ² ∂_μφ∂^μφ'                 : dim(D_PHI,D_PHI, D_D,D_PHI, D_D,D_PHI),
 'g (∂_μφ∂^μφ)²'                 : dim(*([D_D,D_PHI]*4)),
 'g φ²              (= พจน์มวล)' : dim(D_PHI,D_PHI),
 'g φ⁴'                          : dim(*([D_PHI]*4)),
 'g φ⁶'                          : dim(*([D_PHI]*6)),
}
for lbl, rest in terms.items():
    print(f'   {lbl:34s} ที่เหลือ = M^{rest}  ⇒  [g] = {coupling(rest)}')
assert coupling(terms['g φ ∂_μφ∂^μφ      (ข้อที่พลาด)']) == -1
assert coupling(terms['g φ² ∂_μφ∂^μφ'])  == -2
assert coupling(terms['g (∂_μφ∂^μφ)²'])  == -4
assert coupling(terms['g φ²              (= พจน์มวล)']) == 2
assert coupling(terms['g φ⁴']) == 0 and coupling(terms['g φ⁶']) == -2
print('   ⚠️ g φ² ให้ [g] = 2 ⇒ ตรงกับ m² พอดี — พจน์มวลคือค่าคงตัวคู่ควบที่มีมิติ 2')

print('\n── [A_μ] จาก −¼F² (เช็คย้อนบทที่ 5) ' + '─'*22)
A = sp.Symbol('A')
# [F] = [∂] + [A] · [F²] = 2[F] = 4
F_dim = sp.solve(sp.Eq(2*(D_D + A), D_L), A)
print(f'   2([∂] + [A]) = 4  ⇒  [A] = {F_dim[0]}')
assert F_dim[0] == 1
print('   ⇒ สนามเกจมีมิติเท่ากับสนามสเกลาร์พอดี (M ทั้งคู่)')
print(f'   ตรวจซ้ำ: [F] = [∂]+[A] = {D_D + F_dim[0]} ⇒ [F²] = {2*(D_D+F_dim[0])} ✓')

print('\n── กลับทาง: g φ^a ∂_μφ∂^μφ ให้ [g] = 0 ⇒ a = ? ' + '─'*11)
a = sp.Symbol('a')
sol = sp.solve(sp.Eq(D_L - (a*D_PHI + 4), 0), a)
print(f'   4 − (a + 4) = 0  ⇒  a = {sol[0]}')
assert sol[0] == 0
print('   ⇒ พจน์จลน์เปล่า ๆ มีมิติครบ 4 อยู่แล้ว ไม่ต้องมี φ มาเติม')

print('\n── ยืนยันว่าคำตอบที่พลาดต่างจากเฉลยจริง ' + '─'*20)
print(f'   ตอบ 4 = [𝓛] ซึ่งเป็นมิติของ **ทั้งพจน์** ไม่ใช่ของ g เดี่ยว ·'
      f' เฉลย {coupling(5)} ≠ 4')
assert coupling(5) != 4
