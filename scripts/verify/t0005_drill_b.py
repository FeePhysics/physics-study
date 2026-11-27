"""ชุดฝึก T0–T4 จาก [B] Trotta, Bayesian Methods in Cosmology · ตรวจทุกคำตอบด้วย sympy และเทียบกับเฉลยใน [B]

รัน: python3 scripts/verify/t0005_drill_b.py
"""
import sympy as sp
from sympy import Rational as R

ok = 0
def check(name, got, want, tol=None):
    global ok
    if tol is None:
        assert sp.simplify(got - want) == 0, f"{name}: ได้ {got} ต้องการ {want}"
    else:
        assert abs(float(got) - float(want)) <= tol, f"{name}: ได้ {float(got)} ต้องการ {want}"
    ok += 1
    print("✓", name, "=", got if tol is None else round(float(got), 6))

# ── ข้อ 1 · [B] 3.8.1 (i) ตรวจโรค BB · เฉลย [B] น.50: P(+) = 0.059 · P(BB|+) = 0.16 ──
pB, sens, fp = R(1, 100), R(95, 100), R(5, 100)
Pplus = sens * pB + fp * (1 - pB)
check("ข้อ1 ขั้น1 P(+)", Pplus, R(59, 1000))
post = sens * pB / Pplus
check("ข้อ1 ขั้น2 P(ป่วย|+) ≈ 0.161", post, 0.161, tol=5e-4)
# ขั้น 3: โรค 1 ในล้าน · ตรวจไม่พลาดเลย (P(+|ป่วย)=1) · fp 5% ⇒ ในคนผลบวก 100,000 คน ป่วยจริงกี่คน
p6 = R(1, 10**6)
post6 = 1 * p6 / (1 * p6 + fp * (1 - p6))
check("ข้อ1 ขั้น3 ต่อ 100,000 คนผลบวก ≈ 2", post6 * 100000, 2, tol=1e-3)

# ── ข้อ 2 · [B] 3.8.1 (ii) ประตูสามบาน · เฉลย [B] น.51-53 ──
pr = R(1, 3)
openC = {"A": R(1, 2), "B": 1, "C": 0}          # พิธีกรรู้ · เราเลือก A
PC = sum(openC[d] * pr for d in "ABC")
check("ข้อ2 ขั้น1 P(พิธีกรเปิด C)", PC, R(1, 2))
check("ข้อ2 ขั้น2 P(รางวัลหลัง A | เปิด C)", openC["A"] * pr / PC, R(1, 3))
check("ข้อ2 ขั้น3 P(รางวัลหลัง B | เปิด C)", openC["B"] * pr / PC, R(2, 3))
openC2 = {"A": R(1, 2), "B": R(1, 2), "C": 0}   # เราสุ่มเปิดเอง แล้วเจอว่าว่าง
PC2 = sum(openC2[d] * pr for d in "ABC")
check("ข้อ2 ขั้น4 แบบสุ่มเปิดเอง P(B | C ว่าง)", openC2["B"] * pr / PC2, R(1, 2))

# ── ข้อ 3 · [B] 3.8.1 (iii) นักการเมือง · เฉลย [B] น.53: 4/7 ──
AT, BT = R(4, 5), R(1, 4)
PBsays = BT * AT + (1 - BT) * (1 - AT)
check("ข้อ3 ขั้น1 P(B บอกว่าจริง)", PBsays, R(7, 20))
check("ข้อ3 ขั้น2 P(S จริง | B บอกว่าจริง)", BT * AT / PBsays, R(4, 7))

# ── ข้อ 4 · [B] 3.8.2 (v) Gaussian prior N(0, Σ²) + N ครั้ง · ตัวเลขที่เลือก: N=4 σ=2 Σ=3 x̄=10 ──
N, sig, Sig, xbar = 4, 2, 3, 10
check("ข้อ4 ขั้น1 ความกว้างของ likelihood σ/√N", sig / sp.sqrt(N), 1)
mean = xbar * R(Sig**2) / (Sig**2 + R(sig**2, N))       # สูตร [B] น.44
var = 1 / (R(1, Sig**2) + R(N, sig**2))
check("ข้อ4 ขั้น2 จุดกลาง posterior", mean, 9)
check("ข้อ4 ขั้น3 ความแปรปรวน posterior", var, R(9, 10))
# ตรวจซ้ำอีกวิธี: ค่าเฉลี่ยถ่วง precision แบบ T1
check("ข้อ4 ตรวจซ้ำด้วยสูตร T1", (0 * R(1, 9) + xbar * 1) / (R(1, 9) + 1), 9)
Nn = sp.symbols("N", positive=True)
lim = sp.limit(xbar * Sig**2 / (Sig**2 + sig**2 / Nn), Nn, sp.oo)
check("ข้อ4 ขั้น4 N → ∞", lim, 10)

# ── ข้อ 5 · [B] 2.7 (i) อุณหภูมิดาวอังคาร · เฉลย [B] น.19: 199.7 · 1.58 · N = 25 ──
T = [R(x) for x in "191.9 201.6 206.1 200.4 203.2 201.6 196.5 199.5 194.1 202.4".split()]
Tml = sum(T) / len(T)
check("ข้อ5 ขั้น1 T_ML", Tml, R(199730, 1000))
S = 5 / sp.sqrt(10)
check("ข้อ5 ขั้น2 ความคลาดของค่าเฉลี่ย", S, 1.5811, tol=1e-4)
check("ข้อ5 ขั้น3 ครึ่งความกว้าง 2σ", 2 * S, 3.1623, tol=1e-4)
check("ข้อ5 ขั้น4 N ที่ทำให้ 5/√N = 1", sp.solve(sp.Eq(5 / sp.sqrt(Nn), 1), Nn)[0], 25)
# ตรวจด้วย Δχ²: χ²(T) เป็นพาราโบลา ห่าง 1σ ⇒ Δχ² = 1
t = sp.symbols("t")
chi2 = sum(((x - t) / 5) ** 2 for x in T)
check("ข้อ5 ตรวจ Δχ² = 1 ที่ห่าง Σ", sp.expand(chi2.subs(t, Tml + S) - chi2.subs(t, Tml)), 1)

# ── ข้อ 6 · [B] 2.7 (vii) เลเซอร์สี่เครื่อง · เฉลย [B] น.25: 29.2 ± 1.5 · ค่าเฉลี่ยธรรมดา 30.3 ──
F = [R("34.7"), R("28.9"), R("27.1"), R("30.6")]; s = [5, 2, 3, 4]
prec = sum(R(1, x**2) for x in s)
check("ข้อ6 ขั้น1 precision รวม", prec, R(1, 25) + R(1, 4) + R(1, 9) + R(1, 16))
print("   precision รวม ≈", float(prec))
Fml = sum(f / x**2 for f, x in zip(F, s)) / prec
check("ข้อ6 ขั้น2 F_ML ≈ 29.2", Fml, 29.198, tol=1e-3)
check("ข้อ6 ขั้น3 ความคลาด ≈ 1.47", 1 / sp.sqrt(prec), 1.4687, tol=1e-4)
check("ข้อ6 ขั้น4 ค่าเฉลี่ยธรรมดา", sum(F) / 4, R("30.325"))

# ── ข้อ 7 · [B] 2.7 (vi) ดาว + พื้นหลัง · เฉลย [B] น.24: λs = 3 ──
rt, tt, lb = 30, 2, 12
check("ข้อ7 ขั้น1 อัตรารวมที่วัดได้", R(rt, tt), 15)
check("ข้อ7 ขั้น2 อัตราของดาว", R(rt, tt) - lb, 3)
check("ข้อ7 ขั้น3 ถ้านับได้ 20 ใน 2 นาที", R(20, 2) - lb, -2)

print(f"\nผ่าน {ok} ข้อ")
