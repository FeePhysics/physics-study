"""บทที่ 31 · spin ½ — เมทริกซ์ Pauli และการหมุน 2π ได้ −1 (Tong (4.17), (4.26)–(4.28))"""
import sympy as sp
I=sp.I; th=sp.symbols('theta', real=True)
s1=sp.Matrix([[0,1],[1,0]]); s2=sp.Matrix([[0,-I],[I,0]]); s3=sp.Matrix([[1,0],[0,-1]]); E=sp.eye(2)
ok=0
def check(n,c):
    global ok; assert c,n; ok+=1; print("✅",n)
S=[s1,s2,s3]
for a in range(3):
    check(f"σ{a+1} Hermitian", S[a].H==S[a]); check(f"tr σ{a+1}=0", S[a].trace()==0)
    check(f"(σ{a+1})²=1", S[a]*S[a]==E)
check("σ1σ2 = iσ3", s1*s2==I*s3); check("σ2σ3 = iσ1", s2*s3==I*s1); check("σ3σ1 = iσ2", s3*s1==I*s2)
for a in range(3):
    for b in range(3):
        check(f"{{σ{a+1},σ{b+1}}}=2δ", S[a]*S[b]+S[b]*S[a]==(2*E if a==b else sp.zeros(2)))
check("(σ2)² top-left = (-i)(i) = 1", (0*0+(-I)*I)==1)
check("σ1σ2 top-left = i", (s1*s2)[0,0]==I); check("σ2σ1 top-left = -i", (s2*s1)[0,0]==-I)
check("eig σ3 = ±1", set(s3.eigenvals())=={1,-1})
check("spin z = ±1/2", set((s3/2).eigenvals())=={sp.Rational(1,2),-sp.Rational(1,2)})
U=lambda t:(I*t*s3/2).exp()
check("U(θ)=diag(e^{iθ/2},e^{-iθ/2})", sp.simplify(U(th)-sp.diag(sp.exp(I*th/2),sp.exp(-I*th/2)))==sp.zeros(2))
check("U(π) top-left = i", sp.simplify(U(sp.pi)[0,0])==I)
check("U(2π) = -1", sp.simplify(U(2*sp.pi))==-E)
check("U(4π) = +1", sp.simplify(U(4*sp.pi))==E)
R=lambda t:sp.Matrix([[sp.cos(t),-sp.sin(t)],[sp.sin(t),sp.cos(t)]])
check("vector R(2π)=1", R(2*sp.pi)==E)
# สูตรทั่วไป: (n·σ)²=1 ⇒ e^{iθ n·σ/2} = cos(θ/2) + i sin(θ/2) n·σ
n1,n2,n3=sp.Rational(2,3),sp.Rational(1,3),sp.Rational(2,3)   # |n|=1
N=n1*s1+n2*s2+n3*s3
check("(n·σ)²=1", sp.simplify(N*N)==E)
# sympy ทำ exp ของเมทริกซ์ไม่ใช่ทแยงให้ simplify เป็นศูนย์ไม่ได้ ⇒ ตรวจเชิงตัวเลขด้วย scipy
import numpy as np; from scipy.linalg import expm
Nn=np.array(N.evalf(),dtype=complex)
for t in [np.pi/3, 2*np.pi, 1.234]:
    err=np.abs(expm(1j*t*Nn/2)-(np.cos(t/2)*np.eye(2)+1j*np.sin(t/2)*Nn)).max()
    check(f"e^(iθn·σ/2) formula θ={t:.3f} err {err:.0e}", err<1e-12)
check("(σ1)^4 = 1", s1**4==E)
check("cos π = -1, sin π = 0", (sp.cos(sp.pi),sp.sin(sp.pi))==(-1,0))
Nq=sp.Rational(3,5)*s1+sp.Rational(4,5)*s3
check("ข้อ 6: n=(0.6,0,0.8) ⇒ (n·σ)²=1", Nq*Nq==E)
print(f"{ok}/{ok}")
