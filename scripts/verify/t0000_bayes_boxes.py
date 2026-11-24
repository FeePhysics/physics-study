"""Python · fractions — ตรวจทุกตัวเลขใน T0 (lessons/t0000-bayes-boxes.html)
นับแบบ 100 รอบ แล้วเทียบกับสูตร Bayes ตรง ๆ อีกวิธี ต้องตรงกันทุกข้อ"""
from fractions import Fraction as F
ok = 0
def check(name, got, want):
    global ok
    assert got == want, f"{name}: {got} != {want}"
    ok += 1

def rounds(pB, rA, rB, N=100):
    """นับรอบ: คืน (รอบ A, รอบ B, แดงจาก A, แดงจาก B)"""
    nB = N * pB; nA = N - nB
    return nA, nB, nA * rA, nB * rB

def bayes(pB, lA, lB):
    """สูตรตรง: P(B|ข้อมูล) = lB pB / (lB pB + lA pA)"""
    return lB * pB / (lB * pB + lA * (1 - pB))

# เนื้อบท: prior ½ · แดง 0.3 / 0.7
nA, nB, a, b = rounds(F(1,2), F(3,10), F(7,10))
check("เนื้อบท แดงจาก A", a, 15); check("เนื้อบท แดงจาก B", b, 35); check("evidence", a + b, 50)
check("posterior นับ", b / (a + b), F(7,10)); check("posterior สูตร", bayes(F(1,2), F(3,10), F(7,10)), F(7,10))
check("likelihood×prior/evidence", F(7,10) * F(1,2) / F(1,2), F(7,10))
# เนื้อบท ④: ลูกที่สองแดงอีก — prior ใหม่ 0.7
nA, nB, a, b = rounds(F(7,10), F(3,10), F(7,10))
check("④ รอบ B", nB, 70); check("④ แดงจาก B", b, 49); check("④ แดงจาก A", a, 9)
check("④ posterior", b / (a + b), F(49,58))
check("④ = สองลูกพร้อมกัน", bayes(F(1,2), F(9,100), F(49,100)), F(49,58))
# ข้อ 1: prior ½ · แดง 0.2 / 0.6
nA, nB, a, b = rounds(F(1,2), F(1,5), F(3,5))
check("ข้อ1 ขั้น1", a, 10); check("ข้อ1 ขั้น2", b, 30); check("ข้อ1 ขั้น3", a + b, 40)
check("ข้อ1 ขั้น4", b / (a + b), F(3,4)); check("ข้อ1 สูตร", bayes(F(1,2), F(1,5), F(3,5)), F(3,4))
# ข้อ 2: prior B = 0.2 · แดง 0.25 / 0.75
nA, nB, a, b = rounds(F(1,5), F(1,4), F(3,4))
check("ข้อ2 ขั้น1", nA, 80); check("ข้อ2 ขั้น2", a, 20); check("ข้อ2 ขั้น3", b, 15)
check("ข้อ2 ขั้น4", a + b, 35); check("ข้อ2 ขั้น5", b / (a + b), F(3,7))
check("ข้อ2 สูตร", bayes(F(1,5), F(1,4), F(3,4)), F(3,7))
assert F(3,7) < F(1,2)          # B แดงเยอะกว่า แต่ยังน่าจะเป็น A มากกว่า — ประเด็นของข้อ
# ข้อ 4: ต่อจากข้อ 1 · prior ใหม่ ¾ · แดงอีกลูก
nA, nB, a, b = rounds(F(3,4), F(1,5), F(3,5))
check("ข้อ4 ขั้น1", nB, 75); check("ข้อ4 ขั้น2", b, 45); check("ข้อ4 ขั้น3", a, 5)
check("ข้อ4 ขั้น4", a + b, 50); check("ข้อ4 ขั้น5", b / (a + b), F(9,10))
check("ข้อ4 = สองลูกพร้อมกัน", bayes(F(1,2), F(1,25), F(9,25)), F(9,10))
# ข้อ 5: เริ่มใหม่จากข้อ 1 แต่ได้ "ขาว" · ขาว = 1 − แดง
nA, nB, a, b = rounds(F(1,2), 1 - F(1,5), 1 - F(3,5))
check("ข้อ5 ขั้น1", a, 40); check("ข้อ5 ขั้น2", b, 20); check("ข้อ5 ขั้น3", a + b, 60)
check("ข้อ5 ขั้น4", a / (a + b), F(2,3))
# ตัวควบคุม: คำตอบผิดที่น่าจะเจอต้องไม่บังเอิญเท่าคำตอบ
assert F(3,10) != F(3,4)        # ข้อ1: ลืมหารด้วย evidence (0.3)
assert F(3,4) != F(3,7)         # ข้อ2: ลืม prior (ได้ 0.75 เหมือน likelihood)
assert F(3,5) != F(9,10)        # ข้อ4: ใช้ prior ½ ซ้ำ ได้ 0.75 ไม่ใช่ 0.9
# เดิมพันของ Laplace (ส่วน ⓪)
for L, M, d in [(3512, F(3497898, 1000), F(40, 100)), (1070, F(10473486, 10000), F(216, 100))]:
    dev = (F(L) / M - 1) * 100
    assert abs(dev - d) < F(1, 100), (L, float(dev))
    ok += 1
print(f"{ok}/{ok}")
