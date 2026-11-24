"""Python · fractions — ตรวจคำตอบข้อฝึกและตัวอย่างใน T1 ฉบับระฆัง (กฎ 2 ข้อ + สมการ 3)
กฎข้อ 1: precision รวม = 1/Σ² + N/σ² · กฎข้อ 2: μ̄ = (μ/Σ² + m·N/σ²) / precision รวม
ตรวจอีกทางด้วยรูปของ Trotta (9)(10) · ส่วนการเติมกำลังสองอยู่ใน t0001_bayes_gauss.py"""
from fractions import Fraction as F
n = 0
def rule(mu, S2, m, s2, N=1):
    P, D = 1/F(S2), F(N)/F(s2)
    return P, D, P + D, mu*P + m*D, (mu*P + m*D)/(P + D)
def trotta(mu, S2, m, s2):                     # Trotta (9)(10) · r = (σ/Σ)²
    r = F(s2)/F(S2); return (m + r*mu)/(1 + r), F(s2)/(1 + r)
def eq(a, b, name):
    global n; assert a == b, (name, a, b); n += 1
# ตัวอย่างในเนื้อบท + ข้อทาย: 0 ± 3 กับ 10 ± 1 → 9
P, D, T, num_, mb = rule(0, 9, 10, 1)
eq(T, F(10, 9), "ตัวอย่าง precision"); eq(1/T, F(9, 10), "ตัวอย่าง τ²"); eq(mb, 9, "ตัวอย่าง μ̄ = ข้อทาย")
eq(trotta(0, 9, 10, 1), (F(9), F(9, 10)), "ตัวอย่าง = Trotta")
# ข้อ 1
P, D, T, num_, mb = rule(0, 4, 5, 1)
eq(P, F(1, 4), "ข้อ1 ข1"); eq(T, F(5, 4), "ข้อ1 ข2"); eq(1/T, F(4, 5), "ข้อ1 ข3"); eq(num_, 5, "ข้อ1 ข4"); eq(mb, 4, "ข้อ1 ข5")
eq(trotta(0, 4, 5, 1)[0], 4, "ข้อ1 = Trotta")
# ข้อ 2
P, D, T, num_, mb = rule(10, 1, 5, 1)
eq(P, 1, "ข้อ2 ข1"); eq(D, 1, "ข้อ2 ข2"); eq(T, 2, "ข้อ2 ข3"); assert 1/T != T; eq(num_, 15, "ข้อ2 ข4"); eq(mb, F(15, 2), "ข้อ2 ข5"); eq(abs(mb - 4), F(7, 2), "ข้อ2 ห่าง Alice")
# ข้อ 3
P, D, T, num_, mb = rule(10, 1, 5, 1, 99)
eq(D, 99, "ข้อ3 ข1"); eq(T, 100, "ข้อ3 ข2"); eq(10*P, 10, "ข้อ3 ข3"); eq(5*D, 495, "ข้อ3 ข4"); eq(num_, 505, "ข้อ3 ข5"); eq(mb, F(101, 20), "ข้อ3 ข6")
alice = rule(0, 4, 5, 1, 99)[4]
assert abs(mb - alice) < F(7, 100); n += 1      # "เหลือไม่ถึง 0.07"
assert round(float(alice), 2) == 4.99; n += 1
# ข้อ 4: σ → ใหญ่ ⇒ μ̄ → μ
assert abs(rule(3, 1, 10, 10**6)[4] - 3) < F(1, 10**4); n += 1
# กฎข้อ 1: posterior แคบกว่าทั้งคู่เสมอ (สุ่มค่าเป็นเศษส่วน)
import random; random.seed(1)
for _ in range(200):
    S2, s2 = F(random.randint(1, 50), 4), F(random.randint(1, 50), 4)
    tau2 = 1/rule(0, S2, 1, s2)[2]; assert tau2 < S2 and tau2 < s2
n += 1
print(f"{n}/{n}")
