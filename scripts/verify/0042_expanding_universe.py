"""บทที่ 28 · เอกภพขยายตัว — scale factor ⇒ กฎ Hubble ⇒ redshift (Dodelson §1.1–1.2)"""
import sympy as sp
ok = 0
def check(n, c):
    global ok; assert c, n; ok += 1; print("OK", n)

t = sp.symbols('t', positive=True)
x = sp.symbols('x', positive=True)          # comoving distance — ไม่ขึ้นกับเวลา
a = sp.Function('a', positive=True)(t)

# (1)+(2)+(3): d = a x, H = a'/a  ⇒  d' = H d  สำหรับ a(t) ใด ๆ
d = a*x
H = sp.diff(a, t)/a
check("v = d' = H d (a ใด ๆ, x คงที่)", sp.simplify(sp.diff(d, t) - H*d) == 0)

# ตัวอย่างในเนื้อบท: a = t^(2/3) ⇒ H = 2/(3t), H t = 2/3
am = t**sp.Rational(2, 3)
Hm = sp.simplify(sp.diff(am, t)/am)
check("a=t^(2/3) => H = 2/(3t)", sp.simplify(Hm - sp.Rational(2, 3)/t) == 0)
check("H t = 2/3", sp.simplify(Hm*t) == sp.Rational(2, 3))
# a = t^(1/2) (ยุครังสี — กล่าวถึงในเนื้อบท) ⇒ H = 1/(2t)
ar = sp.sqrt(t)
check("a=t^(1/2) => H = 1/(2t)", sp.simplify(sp.diff(ar, t)/ar - 1/(2*t)) == 0)

# ข้อ 1: x = 6 ⇒ d(a=1) = 6, d(a=1/3) = 2
check("ข้อ1 d วันนี้", 1*6 == 6)
check("ข้อ1 d ตอน a=1/3", sp.Rational(1, 3)*6 == 2)
# ข้อ 2: t = 8
check("ข้อ2 a(8) = 4", am.subs(t, 8) == 4)
check("ข้อ2 a'(8) = 1/3", sp.diff(am, t).subs(t, 8) == sp.Rational(1, 3))
check("ข้อ2 H(8) = 1/12", Hm.subs(t, 8) == sp.Rational(1, 12))
check("ข้อ2 H t = 2/3", Hm.subs(t, 8)*8 == sp.Rational(2, 3))
# ข้อ 3: a' = 3, x = 2, a = 6
ad, xx, aa = 3, 2, 6
check("ข้อ3 v = a' x = 6", ad*xx == 6)
check("ข้อ3 d = a x = 12", aa*xx == 12)
check("ข้อ3 H = 1/2", sp.Rational(ad, aa) == sp.Rational(1, 2))
check("ข้อ3 H d = 6 = v", sp.Rational(ad, aa)*aa*xx == ad*xx)
# ข้อ 4: ความชัน
check("ข้อ4 Hubble 1929: 1000/2 = 500", 1000/2 == 500)
check("ข้อ4 v = 72 * 50 = 3600", 72*50 == 3600)
# ข้อ 5: 1 + z = 1/a
z_of = lambda A: 1/sp.Rational(A) - 1
check("ข้อ5 z=1 => a=1/2", z_of(sp.Rational(1, 2)) == 1)
check("ข้อ5 a=1/4 => z=3", z_of(sp.Rational(1, 4)) == 3)
check("ข้อ5 1216*(1+3) = 4864", 1216*4 == 4864)
check("ข้อ5 x=8, z=3 => d=2", sp.Rational(1, 4)*8 == 2)
# ข้อ 6: [t] = M^-1 ⇒ [H] = [a'/a] = M^{+1}
check("ข้อ6 [H] = M^1", -(-1) == 1)
print(f"\n{ok}/{ok} ผ่าน")
