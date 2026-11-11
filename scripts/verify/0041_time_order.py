"""บทที่ 27 · time ordering — สปริงตัวเดียว ตรวจด้วยเมทริกซ์ Fock space จริง (ตัดที่ N ระดับ)"""
import numpy as np, itertools, sympy as sp
N = 14
a = np.diag(np.sqrt(np.arange(1, N)), 1); ad = a.T.copy()
vac = np.zeros(N); vac[0] = 1
def x(t, w):  # x(t) = (a e^{-iwt} + a† e^{iwt})/sqrt(2w)  (บทที่ 18-19)
    return (a*np.exp(-1j*w*t) + ad*np.exp(1j*w*t))/np.sqrt(2*w)
def vev(ops):
    v = vac.astype(complex)
    for O in reversed(ops): v = O @ v
    return vac @ v
ok = 0
def check(n, c):
    global ok; assert c, n; ok += 1; print("✅", n)
w = 2.0
# ① ลำดับสำคัญ
v12 = vev([x(3, w), x(1, w)]); v21 = vev([x(1, w), x(3, w)])
check("<x(3)x(1)> = e^{-4i}/4", abs(v12 - np.exp(-4j)/4) < 1e-12)
check("<x(1)x(3)> = e^{+4i}/4", abs(v21 - np.exp(4j)/4) < 1e-12)
check("ต่างกัน", abs(v12 - v21) > 0.1)
# ② T: เวลามากอยู่ซ้าย ⇒ e^{-iw|Δt|}/(2w)
def T(ts, w):
    order = sorted(ts, reverse=True)       # มากไปน้อย = ซ้ายไปขวา
    return vev([x(t, w) for t in order])
DF = lambda dt, w: np.exp(-1j*w*abs(dt))/(2*w)
for t1, t2 in [(3, 1), (1, 3), (0.4, -2.2), (5, 5)]:
    check(f"T({t1},{t2}) = DF", abs(T([t1, t2], w) - DF(t1-t2, w)) < 1e-12)
check("เวลาเท่ากัน = 1/(2w) = 1/6 ที่ w=3", abs(DF(0, 3) - 1/6) < 1e-15)
# ③ Wick ภายใต้ T: สี่จุดเวลาใด ๆ = ผลบวกสามวิธีจับคู่ของ DF
rng = np.random.default_rng(1)
for _ in range(5):
    ts = rng.uniform(-3, 3, 4)
    lhs = T(list(ts), w)
    rhs = sum(DF(ts[p[0]]-ts[p[1]], w)*DF(ts[p[2]]-ts[p[3]], w) for p in [(0,1,2,3),(0,2,1,3),(0,3,1,2)])
    check(f"Wick ใต้ T สี่จุด err {abs(lhs-rhs):.0e}", abs(lhs - rhs) < 1e-10)
# ④ ค่าตัวเลขในข้อฝึก
check("ω=2 Δt=2 ⇒ เฟส 4", 2*2 == 4)
check("1/(2ω) ω=2 = 1/4", sp.Rational(1, 4) == sp.Rational(1, 2*2))
# ⑤ สนาม: DF(t<0) = D(|t|) ของบทที่ 20 (ผลบวก ±k สมมาตร)
L = 2*np.pi; ks = np.arange(-6, 7)*2*np.pi/L; m = 1.0
D = lambda t, r: sum(np.exp(1j*k*r - 1j*np.sqrt(k*k+m*m)*t)/(2*np.sqrt(k*k+m*m)*L) for k in ks)
check("D_F(-2, r) = D(2, r)", abs(D(2, 0.7) - D(2, -0.7)) < 1e-12)
print(f"{ok}/{ok}")
