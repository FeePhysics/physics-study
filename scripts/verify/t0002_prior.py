"""Python · sympy — ตรวจทุกข้ออ้างของ T2 (prior และการเปลี่ยนตัวแปร)
ตัวอย่างหลัก: สี่เหลี่ยมจัตุรัส ด้าน L ∈ [1,3] · พื้นที่ A = L² ∈ [1,9]
กฎ: p_A(A) = p_L(L) / |dA/dL|  (ความน่าจะเป็นในแถบเล็ก ๆ ต้องเท่ากันทั้งสองแกน)
Jeffreys: ถ้า ln A แบน ⇒ p(A) ∝ 1/A"""
import sympy as sp
L, A, u = sp.symbols("L A u", positive=True)
n = 0
def ok(cond, name):
    global n
    assert bool(cond), name; n += 1

pL = sp.Rational(1, 2)                                   # แบนบน [1,3]
ok(sp.integrate(pL, (L, 1, 3)) == 1, "p_L normalize")
# ⓪ ข้อทาย: P(L<2)
ok(sp.integrate(pL, (L, 1, 2)) == sp.Rational(1, 2), "แบนบน L: P(L<2) = 1/2")
pA_flat = sp.Rational(1, 8)                              # แบนบน [1,9]
ok(sp.integrate(pA_flat, (A, 1, 4)) == sp.Rational(3, 8), "แบนบน A: P(L<2)=P(A<4) = 3/8")
# ② กฎเปลี่ยนตัวแปร: p_A = p_L / (dA/dL) แล้วเขียนด้วย A
dAdL = sp.diff(L**2, L)
pA = sp.simplify((pL / dAdL).subs(L, sp.sqrt(A)))
ok(sp.simplify(pA - 1/(4*sp.sqrt(A))) == 0, "p_A = 1/(4√A)")
ok(sp.integrate(pA, (A, 1, 9)) == 1, "p_A normalize — กฎถูก")
ok(sp.integrate(pA, (A, 1, 4)) == sp.Rational(1, 2), "แถบเดียวกัน ความน่าจะเป็นเท่ากันทั้งสองแกน")
# การนับในเนื้อบท: แถบ L∈[1,1.2] กับ [2.8,3] ได้ความน่าจะเป็นเท่ากัน แต่กว้างบนแกน A ไม่เท่ากัน
w1 = sp.Rational(6, 5)**2 - 1; w2 = 9 - sp.Rational(14, 5)**2
ok(w1 == sp.Rational(11, 25) and w2 == sp.Rational(29, 25), "ความกว้างแถบบนแกน A: 0.44 กับ 1.16")
# ข้อ 2
ok(dAdL.subs(L, 1) == 2 and dAdL.subs(L, 3) == 6 and dAdL.subs(L, 2) == 4, "dA/dL ที่ L = 1, 2, 3")
ok(pA.subs(A, 1) == sp.Rational(1, 4) and pA.subs(A, 9) == sp.Rational(1, 12), "p_A(1) = 1/4 · p_A(9) = 1/12")
ok(pA.subs(A, 4) == sp.Rational(1, 8), "p_A(4) = 1/8")
# ข้อ 3: amplitude A ∈ [1,100]
ok(sp.Rational(10 - 1, 100 - 1) == sp.Rational(1, 11), "แบนบน A: P(A<10) = 9/99 = 1/11")
ok(sp.Rational(1 - 0, 2 - 0) == sp.Rational(1, 2), "แบนบน log10 A ∈ [0,2]: P(A<10) = 1/2")
# เนื้อบท: [1e-10, 1e-8] แบน ⇒ P(A < 1e-9) = 1/11 เช่นกัน
e = sp.Rational(1, 10**10)
ok((10*e - e) / (100*e - e) == sp.Rational(1, 11), "ตัวอย่าง A_s แบนบน [1e-10, 1e-8]: P(A < 1e-9) = 1/11")
# ④ Jeffreys: u = ln A แบน ⇒ p_A = c · du/dA = c/A
c = sp.symbols("c", positive=True)
pJ = c * sp.diff(sp.log(A), A)
ok(sp.simplify(pJ - c/A) == 0, "Jeffreys ∝ 1/A")
ok(sp.simplify(pJ.subs(A, 10) / pJ.subs(A, 1)) == sp.Rational(1, 10), "ข้อ 4: อัตราส่วน 1/10")
# Jeffreys: ทุกช่วงสิบเท่ามีความน่าจะเป็นเท่ากัน
ok(sp.simplify(sp.integrate(1/A, (A, 1, 10)) - sp.integrate(1/A, (A, 10, 100))) == 0, "ทุก decade เท่ากัน")
# ข้อ 1: เขียว (แบนบน A) · L < 2.5
a = sp.Rational(5, 2)**2
ok(a == sp.Rational(25, 4) and a - 1 == sp.Rational(21, 4), "ข้อ1 ขั้น1-2")
ok((a - 1) / 8 == sp.Rational(21, 32), "ข้อ1 ขั้น4 = 21/32")
ok(sp.integrate(pA, (A, 1, a)) == sp.Rational(3, 4), "แดง (แบนบน L) ได้ 0.75 — ตัวควบคุมใน fb")
ok(sp.Rational(21, 32) != sp.Rational(3, 4), "สองคำตอบไม่บังเอิญเท่ากัน")
# ④ รูปของ Trotta (11): p(ψ) = p(θ)|dθ/dψ| · θ = L, ψ = A
dLdA = sp.diff(sp.sqrt(A), A)
ok(dLdA.subs(A, 4) == sp.Rational(1, 4), "|dL/dA| ที่ A = 4 = 1/4")
ok(sp.simplify(pL * dLdA - pA) == 0, "Trotta (11) = สมการ (1) ของบท ทุกจุด")
# ⑤ / ข้อ 5: Cobaya รูป 9 · logamplitude แบนบน [-2, -0.6]
lo, hi = sp.Integer(-2), sp.Rational(-3, 5)
ok(sp.Integer(10)**lo == sp.Rational(1, 100), "ข้อ5 ขั้น1")
ok(hi - lo == sp.Rational(7, 5), "ข้อ5 ขั้น2")
ok(sp.log(sp.Rational(1, 10), 10) == -1, "ข้อ5 ขั้น3")
ok((-1 - lo) / (hi - lo) == sp.Rational(5, 7), "ข้อ5 ขั้น4 = 5/7")
flat = (sp.Rational(1, 10) - sp.Rational(1, 100)) / (10**hi - sp.Rational(1, 100))
ok(abs(float(flat) - 0.37) < 0.01, "ตัวควบคุมใน fb: แบนบน amplitude ได้ราว 0.37")
# รอบสอง: ขั้นใหม่ของข้อ 3 และตัวควบคุมคำตอบผิดจริง
ok(sp.Integer(100) - 1 == 99, "ข้อ3 ขั้น2 = 99")
ok(sp.Rational(1, 12) != sp.Rational(1, 6) and sp.Rational(1, 11) != sp.Rational(1, 10), "คำตอบผิดจริงไม่บังเอิญเท่าคีย์")
print(f"{n}/{n}")
