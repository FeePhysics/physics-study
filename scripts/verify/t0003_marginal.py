"""Python · sympy — ตรวจทุกข้ออ้างของ T3 (สองพารามิเตอร์ · marginalise)
ตัวอย่างหลัก (Trotta หน้า 74-75): แสงรวม T = s + b · s = แหล่งที่สนใจ · b = พื้นหลัง (nuisance)
likelihood ของแสงรวม: exp(-(s+b-T)²/2σT²) · ความรู้เรื่องพื้นหลัง: exp(-(b-m_b)²/2σb²) · prior บน s แบน"""
import sympy as sp
from fractions import Fraction as Fr
s, b = sp.symbols("s b", real=True)
T, mb = sp.symbols("T m_b", real=True)
sT, sb = sp.symbols("sigma_T sigma_b", positive=True)
n = 0
def ok(c, name):
    global n
    assert bool(c), name; n += 1

joint = sp.exp(-(s + b - T)**2 / (2*sT**2)) * sp.exp(-(b - mb)**2 / (2*sb**2))
# (1) marginal: อินทิเกรต b ทิ้ง
marg = sp.simplify(sp.integrate(joint, (b, -sp.oo, sp.oo)))
# ต้องเป็น Gaussian ใน s จุดกลาง T - m_b ความแปรปรวน σT² + σb²
ratio = sp.simplify(marg / sp.exp(-(s - (T - mb))**2 / (2*(sT**2 + sb**2))))
ok(sp.simplify(sp.diff(ratio, s)) == 0, "(2) marginal ∝ N(T−m_b, σT²+σb²)")
# ถ้ารู้ b แน่นอน (= m_b): ตัด b = m_b ⇒ กว้าง σT
cond = joint.subs(b, mb)
ok(sp.simplify(sp.diff(sp.simplify(cond / sp.exp(-(s - (T - mb))**2/(2*sT**2))), s)) == 0, "slice ที่ b = m_b กว้าง σT")
# profile (maximise b ทีละ s) ให้รูปร่างเดียวกับ marginal สำหรับ Gaussian (Trotta §5.3.2)
lnj = sp.log(joint)
bstar = sp.solve(sp.diff(lnj, b), b)[0]
prof = sp.simplify(joint.subs(b, bstar))
ok(sp.simplify(sp.diff(sp.simplify(prof / sp.exp(-(s-(T-mb))**2/(2*(sT**2+sb**2)))), s)) == 0, "profile = marginal (Gaussian)")
# สหสัมพันธ์ของ s กับ b เป็นลบ: ρ = −σb/√(σT²+σb²)
# precision matrix จาก -2 ln joint
Q = sp.expand(-2*lnj)
P = sp.Matrix([[sp.diff(Q, s, 2), sp.diff(Q, s, b)], [sp.diff(Q, b, s), sp.diff(Q, b, 2)]]) / 2
C = sp.simplify(P.inv())
rho = sp.simplify(C[0, 1] / sp.sqrt(C[0, 0]*C[1, 1]))
ok(sp.simplify(rho + sb/sp.sqrt(sT**2 + sb**2)) == 0, "ρ = −σb/√(σT²+σb²)")
ok(sp.simplify(C[0, 0] - (sT**2 + sb**2)) == 0, "ความแปรปรวนของ s จาก covariance = σT²+σb²")
# ตัวอย่างในเนื้อบท: T = 10 ± 1 · b = 3 ± 2
ok((sp.sqrt(sp.Integer(1) + 4)) == sp.sqrt(5), "เนื้อบท: s = 7 ± √5")
# ข้อ 3: T = 12 ± 3 · b = 5 ± 4 ⇒ s = 7 ± 5 · ถ้ารู้ b: ± 3
ok(12 - 5 == 7 and sp.sqrt(9 + 16) == 5, "ข้อ3: 7 ± 5 (3-4-5)")

# ── ตารางแบบนับ (T0) ──
# ข้อทาย: แถว s = 1,2,3 · คอลัมน์ b = 0,1,2 · รวม 100 รอบ
W = {1: [0, 0, 28], 2: [14, 14, 14], 3: [10, 10, 10]}
ok(sum(map(sum, W.values())) == 100, "ข้อทาย รวม 100")
prof_pick = max(W, key=lambda k: max(W[k])); marg_pick = max(W, key=lambda k: sum(W[k]))
ok(prof_pick == 1 and marg_pick == 2, "ข้อทาย: ช่องสูงสุดชี้ s=1 · รวมแถวชี้ s=2")
# ข้อ 1: แถว s = 1,2,3 · คอลัมน์ b = 0,1
E1 = {1: [10, 15], 2: [25, 20], 3: [5, 25]}
ok(sum(map(sum, E1.values())) == 100, "ข้อ1 รวม 100")
ok(Fr(sum(E1[2]), 100) == Fr(9, 20), "ข้อ1: P(s=2) = 45/100")
ok(Fr(E1[1][0] + E1[2][0] + E1[3][0], 100) == Fr(2, 5), "ข้อ1: P(b=0) = 40/100")
ok(max(E1, key=lambda k: sum(E1[k])) == 2, "ข้อ1: marginal ชี้ s = 2")
# ข้อ 2: ตารางที่ maximise กับ marginalise ไม่ตรงกัน
E2 = {1: [5, 5, 5, 5], 2: [0, 0, 32, 0], 3: [12, 12, 12, 12]}
ok(sum(map(sum, E2.values())) == 100, "ข้อ2 รวม 100")
ok(max(E2, key=lambda k: max(E2[k])) == 2, "ข้อ2: ช่องสูงสุดชี้ s = 2")
ok(max(E2, key=lambda k: sum(E2[k])) == 3, "ข้อ2: รวมแถวชี้ s = 3")
ok(Fr(sum(E2[3]), 100) == Fr(12, 25), "ข้อ2: P(s=3) = 48/100")
print(f"{n}/{n}")
