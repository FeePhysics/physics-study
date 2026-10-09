"""ตรวจทุกข้ออ้างของบทที่ 18 — ปุ่มขึ้นลงบันได a, a† · นับด้วย a†a · [a, a†] = 1 · หลายโหมด

หน่วย hbar = c = 1 · สถานะ |n> = เวกเตอร์ฐานตัวที่ n
⚠️ เมทริกซ์ตัดที่ขนาด D จะผิดที่แถวสุดท้ายเสมอ ([a, a†] ที่มุมล่างขวาไม่ใช่ 1)
   ⇒ ตรวจเฉพาะสถานะ |n> ที่ n + 2 < D
"""
import sympy as sp

ok = []


def chk(name, got, want):
    ok.append((name, sp.simplify(sp.expand(got - want)) == 0, got, want))


D = 8
ket = lambda n: sp.Matrix([1 if i == n else 0 for i in range(D)])
a = sp.zeros(D, D)
for n in range(1, D):
    a[n - 1, n] = sp.sqrt(n)          # a|n> = sqrt(n)|n-1>
ad = a.H                              # a† = ทรานสโพสสังยุค
ok.append(("a† ที่ได้จาก dagger มี sqrt(n+1) ใต้แนวทแยง",
           all(ad[n + 1, n] == sp.sqrt(n + 1) for n in range(D - 1)), None, None))

# ── ① ขึ้นลงบันได ───────────────────────────────────────────
for n in range(5):
    ok.append((f"a†|{n}> = sqrt({n+1}) |{n+1}>", ad*ket(n) == sp.sqrt(n + 1)*ket(n + 1), None, None))
    if n:
        ok.append((f"a|{n}> = sqrt({n}) |{n-1}>", a*ket(n) == sp.sqrt(n)*ket(n - 1), None, None))
ok.append(("a|0> = เวกเตอร์ศูนย์", a*ket(0) == sp.zeros(D, 1), None, None))
# 4x4 ที่โชว์ในบท
a4 = sp.Matrix([[0, 1, 0, 0], [0, 0, sp.sqrt(2), 0], [0, 0, 0, sp.sqrt(3)], [0, 0, 0, 0]])
ok.append(("เมทริกซ์ 4x4 ในบท = มุมบนซ้ายของ a", a[:4, :4] == a4, None, None))
# ตัวอย่างเมทริกซ์คูณเวกเตอร์ในข้อ ① (ฉบับแก้): a† 4x4 คูณคอลัมน์ |2> ได้ sqrt3 |3>
ad4 = sp.Matrix([[0, 0, 0, 0], [1, 0, 0, 0], [0, sp.sqrt(2), 0, 0], [0, 0, sp.sqrt(3), 0]])
ok.append(("a† 4x4 ที่พิมพ์ในบท = ทรานสโพสสังยุคของ a 4x4", ad4 == a4.H, None, None))
ok.append(("ตัวอย่าง: a† [0,0,1,0] = [0,0,0,sqrt3] = sqrt3 |3>",
           ad4*sp.Matrix([0, 0, 1, 0]) == sp.Matrix([0, 0, 0, sp.sqrt(3)]), None, None))
ok.append(("ตัวควบคุม: ถ้าใช้ a แทน a† ได้ sqrt2 |1> ไม่ใช่ |3>",
           a4*sp.Matrix([0, 0, 1, 0]) == sp.Matrix([0, sp.sqrt(2), 0, 0]), None, None))
# ลำดับการกด: ตัวติดเค็ตทำก่อน — a†a|3> (ลงก่อน) ต่างจาก a a†|3> (ขึ้นก่อน)
ok.append(("a†a|3> = 3|3> แต่ a a†|3> = 4|3>", ad*a*ket(3) == 3*ket(3) and a*ad*ket(3) == 4*ket(3), None, None))

# ── ② นับ: a†a = diag(n) ────────────────────────────────────
N = ad*a
ok.append(("a†a = diag(0, 1, ..., D-1)", N == sp.diag(*range(D)), N, None))
H = lambda w: w*(N + sp.eye(D)/2)
w = sp.symbols('w', positive=True)
for n in range(5):
    chk(f"H|{n}> = w({n} + 1/2)|{n}>  (เท่าบทที่ 17)", (H(w)*ket(n))[n], w*(n + sp.Rational(1, 2)))

# ── ③ [a, a†] = 1 (ยกเว้นแถวสุดท้ายของเมทริกซ์ตัด) ───────────
C = a*ad - ad*a
ok.append(("[a, a†] = 1 บนสถานะ |0>..|D-2>", all(C[i, i] == 1 for i in range(D - 1))
           and all(C[i, j] == 0 for i in range(D) for j in range(D) if i != j), C, None))
ok.append(("ตัวควบคุม: มุมล่างขวาของเมทริกซ์ตัด ≠ 1 (เพราะตัด ไม่ใช่ฟิสิกส์)", C[D - 1, D - 1] != 1, C[D - 1, D - 1], "≠1"))
# H = 1/2 p^2 + 1/2 w^2 x^2 กับ x = (a + a†)/sqrt(2w), p = i sqrt(w/2)(a† - a) ⇒ w(a†a + 1/2) (นอกแถวสุดท้าย)
X = (a + ad)/sp.sqrt(2*w)
P = sp.I*sp.sqrt(w/2)*(ad - a)
Hxp = sp.simplify(P*P/2 + w**2*X*X/2)
ok.append(("1/2 p^2 + 1/2 w^2 x^2 = w(a†a + 1/2) บน |0>..|D-2>",
           all(sp.simplify(Hxp[i, j] - H(w)[i, j]) == 0 for i in range(D - 1) for j in range(D - 1)), None, None))

# ── ตัวเลขในข้อฝึก ──────────────────────────────────────────
coef = lambda v, n: v[n]                      # สัมประสิทธิ์หน้า |n>
chk("ข้อ 1: a†|2> = c|3>, c^2 = 3", coef(ad*ket(2), 3)**2, 3)
chk("ข้อ 1: a|2> = c|1>, c^2 = 2", coef(a*ket(2), 1)**2, 2)
chk("ข้อ 1: ขนาดของ a|0> = 0", (a*ket(0)).norm(), 0)
chk("ข้อ 2: a|3> = c|2>, c^2 = 3", coef(a*ket(3), 2)**2, 3)
chk("ข้อ 2: a†a|3> = 3|3>", coef(N*ket(3), 3), 3)
chk("ข้อ 2: เทรซของ a†a ขนาด 4x4 = 6", (ad[:4, :4]*a[:4, :4]).trace(), 6)
chk("ข้อ 3: a a†|2> = 3|2>", coef(a*ad*ket(2), 2), 3)
chk("ข้อ 3: a†a|2> = 2|2>", coef(N*ket(2), 2), 2)
chk("ข้อ 3: ต่างกัน 1", coef(C*ket(2), 2), 1)
chk("ข้อ 4: w=4, |3> ⇒ E = 14", coef(H(4)*ket(3), 3), 14)
chk("ข้อ 4: w=4, |0> ⇒ E = 2", coef(H(4)*ket(0), 0), 2)
chk("ข้อ 4: เทียบสุญญากาศ = 12 = w a†a", 14 - 2, 4*3)
Enw = lambda ww, n: ww*(n + sp.Rational(1, 2))
chk("ข้อ 5: ทั้งหมด 2*3.5 + 6*1.5 = 16", Enw(2, 3) + Enw(6, 1), 16)
chk("ข้อ 5: สุญญากาศ = 4", Enw(2, 0) + Enw(6, 0), 4)
chk("ข้อ 5: เทียบสุญญากาศ = 12 = 3*2 + 1*6", 16 - 4, 3*2 + 1*6)
wk = lambda kk, mm: sp.sqrt(kk**2 + mm**2)
chk("ข้อ 6: m=4, k=±3 ⇒ w = 5 ทั้งคู่ ⇒ E = 10", wk(3, 4) + wk(-3, 4), 10)
chk("ข้อ 6: P = 3 + (-3) = 0", 3 + (-3), 0)
chk("ข้อ 7: a†a†|0> = sqrt(2)|2> ⇒ c^2 = 2", coef(ad*ad*ket(0), 2)**2, 2)
chk("ข้อ 7: E = 2*5 = 10", 2*wk(3, 4), 10)
chk("ข้อ 7: P = 2*3 = 6", 2*3, 6)
# กล่อง "k คือโมเมนตัมแล้ว ไม่ต้องคูณมวล": m γ v = k ทั่วไป และตัวเลข m=4 k=3
kk, mm = sp.symbols('k m', positive=True)
vv = kk/sp.sqrt(kk**2 + mm**2)
chk("m γ v = k ทั่วไป (v = k/ω)", sp.simplify(mm*vv/sp.sqrt(1 - vv**2)), kk)
chk("ตัวอย่าง: v = 3/5", 3/wk(3, 4), sp.Rational(3, 5))
chk("ตัวอย่าง: γ = 1.25", 1/sp.sqrt(1 - sp.Rational(9, 25)), sp.Rational(5, 4))
chk("ตัวอย่าง: m γ v = 4 * 1.25 * 0.6 = 3", 4*sp.Rational(5, 4)*sp.Rational(3, 5), 3)
ok.append(("ตัวควบคุม: คำตอบ 24 = (2*3)*m คือคูณมวลซ้ำ ≠ P", 2*3*4 != 6, 24, 6))

bad_ = [r for r in ok if not r[1]]
for name, good, got, want in ok:
    print(("  ✓ " if good else "  ✗ ") + name + ("" if good else f"\n      ได้ {got}\n      ควรเป็น {want}"))
print(f"\n{'✅ ผ่าน' if not bad_ else '❌ ตก'} — {len(ok)-len(bad_)}/{len(ok)}")
raise SystemExit(1 if bad_ else 0)
