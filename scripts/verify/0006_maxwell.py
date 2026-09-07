"""บทที่ 5 — แมกซ์เวลล์จาก Lagrangian · ตรวจทุกขั้นด้วย sympy"""
import sympy as sp
import itertools

t, x, y, z = sp.symbols('t x y z', real=True)
X = [t, x, y, z]
g = sp.diag(1, -1, -1, -1)          # (+−−−) แบบฟิสิกส์อนุภาค

# ศักย์สี่ A_μ(t,x,y,z) — ยังไม่กำหนดรูป
A = [sp.Function(f'A{i}')(t, x, y, z) for i in range(4)]

def d(f, mu):                       # ∂_μ
    return sp.diff(f, X[mu])

print("=== ① F_μν = ∂_μA_ν − ∂_νA_μ เป็นปฏิสมมาตรโดยโครงสร้าง ===")
F = sp.Matrix(4, 4, lambda m, n: d(A[n], m) - d(A[m], n))
print("  F_μν + F_νμ = 0 ทุกคู่ :", sp.simplify(F + F.T) == sp.zeros(4, 4))
print("  แนวทแยงเป็นศูนย์เสมอ  :", all(sp.simplify(F[i, i]) == 0 for i in range(4)))

print("\n=== ② องค์ประกอบของ F คือสนามไฟฟ้าและสนามแม่เหล็กจริง ===")
# ในหน่วย c=1: E_i = F_{i0} = ∂_i A_0 − ∂_0 A_i  (เครื่องหมายตาม convention)
Ex = sp.simplify(F[1, 0]); print("  F₁₀ =", Ex, "  ← Eₓ (ศักย์ไฟฟ้า + การเปลี่ยนของ A)")
Bz = sp.simplify(F[1, 2]); print("  F₁₂ =", Bz, "  ← −B_z (เคิร์ลของ A)")

print("\n=== ③ EL ของ ℒ = −¼F_μν F^μν ให้สมการแมกซ์เวลล์ ∂_μF^μν = 0 ===")
g_inv = g.inv()
F_up = g_inv * F * g_inv.T                              # ยกดัชนีทั้งสอง
Lag = -sp.Rational(1,4) * sum(F[m,n]*F_up[m,n] for m in range(4) for n in range(4))
Lag = sp.expand(Lag)

for nu in range(2):                                      # ตรวจสองสมการก็พอ (ที่เหลือรูปเดียวกัน)
    el = -sum(d(sp.diff(Lag, d(A[nu], mu)), mu) for mu in range(4))   # ∂_μ(∂ℒ/∂(∂_μA_ν))
    target = sum(d(F_up[mu, nu], mu) for mu in range(4))              # ∂_μ F^{μν}
    print(f"  ν={nu}: EL ตรงกับ ∂_μF^μν :", sp.simplify(el - target) == 0)

print("\n=== ④ เอกลักษณ์เบียงกี ∂_[λ F_μν] = 0 ⇒ กฎฟาราเดย์ + ไม่มีขั้วแม่เหล็กเดี่ยว ===")
ok = all(sp.simplify(d(F[m,n], l) + d(F[n,l], m) + d(F[l,m], n)) == 0
         for l, m, n in itertools.combinations(range(4), 3))
print("  จริงทุกชุดดัชนี :", ok, " ← เป็นจริงโดยอัตโนมัติเพราะ F มาจาก A")

print("\n=== ⑤ สมมาตรเกจ: A_μ → A_μ + ∂_μΛ แล้ว F ไม่เปลี่ยน ===")
Lam = sp.Function('Lambda')(t, x, y, z)
A2 = [A[i] + d(Lam, i) for i in range(4)]
F2 = sp.Matrix(4, 4, lambda m, n: d(A2[n], m) - d(A2[m], n))
print("  F เดิม = F ใหม่ :", sp.simplify(F2 - F) == sp.zeros(4, 4))
print("  ⇒ A มีอิสระที่วัดไม่ได้ · ฟิสิกส์อยู่ที่ F ไม่ใช่ A")

print("\n=== ⑥ มิติเชิงมวลของสนามเกจ (ħ=c=1) ===")
print("  [ℒ]=M⁴ · ℒ มี (∂A)² ⇒ M²·[A]² = M⁴ ⇒ [A] = M¹  ← เท่ากับสนามสเกลาร์")
print("  พจน์มวล m²A² จะมีมิติ M⁴ พอดีเช่นกัน — แต่มันทำลายสมมาตรเกจ (ข้อ ⑤)")
mA = sp.Matrix(4, 4, lambda m, n: 0)   # ตรวจว่าพจน์มวลไม่คงตัวภายใต้เกจ
lhs = sum(g_inv[i,i]*A[i]**2 for i in range(4))
rhs = sum(g_inv[i,i]*A2[i]**2 for i in range(4))
print("  A_μA^μ คงตัวภายใต้เกจไหม :", sp.simplify(rhs - lhs) == 0, " ⇒ โฟตอนไร้มวลโดยบังคับ")
