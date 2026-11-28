"""รูปประกอบตำรา Bayes ฉบับภาษาไทย (ตาม Trotta [B] บทที่ 2–3) · วาดเองจากตัวอย่างของตำรา ไม่ลอกรูปต้นฉบับ

รัน: python3 scripts/verify/trotta_b_figs.py  → เขียน SVG ทับระหว่าง <!--FIG:x--> ... <!--/FIG:x-->

  coin      likelihood ของเหรียญ 5 ครั้ง เมื่อได้หัว 1, 2, 4 ครั้ง (แกนนอน θ · แกนตั้ง L)
  converge  prior สองแบบ → หลังวัด 1 ครั้ง → หลังวัด 100 ครั้ง (ตัวเลขตรงกับ trotta_b_ch2_3.py)
"""
import math
import re
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
PAGE = ROOT / "reference" / "trotta-b-ch2-3.html"
BLUE, ORANGE, RED, GREEN = "#3b6fb6", "#c8692a", "#d62728", "#2ca02c"
FONT = 'font-family="Sarabun, sans-serif" fill="currentColor"'


def path(pts):
    return "M" + " L".join(f"{x:.1f},{y:.1f}" for x, y in pts)


def fig_coin():
    W, H, L, B, T, R = 400, 270, 50, 48, 14, 385
    X = lambda t: L + t * (R - L)
    Y = lambda v: H - B - v / 0.45 * (H - B - T)
    o = [f'<line x1="{L}" y1="{T}" x2="{L}" y2="{H-B}" stroke="currentColor" opacity=".5"/>',
         f'<line x1="{L}" y1="{H-B}" x2="{R}" y2="{H-B}" stroke="currentColor" opacity=".5"/>']
    for t in (0, 0.2, 0.4, 0.6, 0.8, 1):
        o.append(f'<text x="{X(t):.1f}" y="{H-B+16}" font-size="13" text-anchor="middle" {FONT}>{t:g}</text>')
    for v in (0, 0.1, 0.2, 0.3, 0.4):
        o.append(f'<text x="{L-6}" y="{Y(v)+4:.1f}" font-size="13" text-anchor="end" {FONT}>{v:g}</text>')
    o.append(f'<text x="{(L+R)/2}" y="{H-8}" font-size="14" text-anchor="middle" {FONT}>θ = โอกาสออกหัวในหนึ่งครั้ง</text>')
    o.append(f'<text transform="translate(15,{(T+H-B)/2}) rotate(-90)" font-size="14" text-anchor="middle" {FONT}>likelihood L(θ)</text>')
    for r, col, dash in ((1, RED, ""), (2, GREEN, ' stroke-dasharray="6 4"'), (4, BLUE, ' stroke-dasharray="2 3"')):
        c = math.comb(5, r)
        pts = [(X(t), Y(c * t**r * (1 - t) ** (5 - r))) for t in [i / 200 for i in range(201)]]
        o.append(f'<path d="{path(pts)}" fill="none" stroke="{col}" stroke-width="3"{dash}/>')
        o.append(f'<text x="{X(r/5):.1f}" y="{Y(c*(r/5)**r*(1-r/5)**(5-r))-8:.1f}" font-size="13" text-anchor="middle" '
                 f'fill="{col}" font-family="Sarabun, sans-serif">หัว {r} ครั้ง</text>')
    return (f'<figure><svg role="img" aria-label="likelihood ของเหรียญห้าครั้ง ได้หัวหนึ่งครั้งยอดอยู่ที่ 0.2 สองครั้งยอดที่ 0.4 สี่ครั้งยอดที่ 0.8" '
            f'viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" style="max-width:{W}px">{"".join(o)}</svg></figure>')


def gauss(x, m, s):
    return math.exp(-((x - m) ** 2) / (2 * s * s)) / s


def upd(m0, s0, xbar, s, n):
    p = 1 / s0**2 + n / s**2
    return (m0 / s0**2 + n * xbar / s**2) / p, 1 / math.sqrt(p)


def fig_converge():
    W, PH, L, R = 400, 110, 50, 385
    panels = [("prior ก่อนวัด", (750, 100), (900, 20)),
              ("หลังวัด 1 ครั้ง (ได้ 850 ± 50)", upd(750, 100, 850, 50, 1), upd(900, 20, 850, 50, 1)),
              ("หลังวัด 100 ครั้ง (เฉลี่ย 850)", upd(750, 100, 850, 50, 100), upd(900, 20, 850, 50, 100))]
    X = lambda v: L + (v - 600) / 400 * (R - L)
    o = []
    top = 8
    for k, (title, A, Bp) in enumerate(panels):
        y0 = top + k * (PH + 14)
        base = y0 + PH - 22
        hmax = max(gauss(A[0], *A), gauss(Bp[0], *Bp))
        o.append(f'<text x="{L}" y="{y0+12}" font-size="14" font-weight="600" {FONT}>{title}</text>')
        for (m, s), col, dash in ((A, BLUE, ""), (Bp, ORANGE, ' stroke-dasharray="6 4"')):
            pts = [(X(v), base - (PH - 42) * gauss(v, m, s) / hmax) for v in [600 + i * 2 for i in range(201)]]
            o.append(f'<path d="{path(pts)}" fill="none" stroke="{col}" stroke-width="2.5"{dash}/>')
        o.append(f'<line x1="{L}" y1="{base}" x2="{R}" y2="{base}" stroke="currentColor" opacity=".5"/>')
        for v in (600, 700, 800, 900, 1000):
            o.append(f'<text x="{X(v):.1f}" y="{base+15}" font-size="12" text-anchor="middle" {FONT}>{v}</text>')
    yb = top + 3 * (PH + 14) + 4
    o.append(f'<text x="{(L+R)/2}" y="{yb}" font-size="14" text-anchor="middle" {FONT}>ค่าที่อยากรู้ θ (แกนตั้งของแต่ละช่อง = ความเชื่อ)</text>')
    o.append(f'<line x1="{L}" y1="{yb+16}" x2="{L+24}" y2="{yb+16}" stroke="{BLUE}" stroke-width="3"/>'
             f'<text x="{L+30}" y="{yb+20}" font-size="13" {FONT}>นักวิจัย A (เริ่ม 750 ± 100)</text>')
    o.append(f'<line x1="{L+200}" y1="{yb+16}" x2="{L+224}" y2="{yb+16}" stroke="{ORANGE}" stroke-width="3" stroke-dasharray="6 4"/>'
             f'<text x="{L+230}" y="{yb+20}" font-size="13" {FONT}>B (เริ่ม 900 ± 20)</text>')
    H = yb + 30
    return (f'<figure><svg role="img" aria-label="สามช่อง ก่อนวัด prior ของ A กว้างที่ 750 ของ B แคบที่ 900 หลังวัดหนึ่งครั้งยังต่างกัน '
            f'หลังวัดร้อยครั้งทับกันที่ราว 850" viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" style="max-width:{W}px">'
            f'{"".join(o)}</svg></figure>')


FIGS = {"coin": fig_coin, "converge": fig_converge}

if __name__ == "__main__":
    html = PAGE.read_text(encoding="utf-8")
    for k, f in FIGS.items():
        pat = re.compile(rf"(<!--FIG:{k}-->).*?(<!--/FIG:{k}-->)", re.S)
        assert pat.search(html), f"ไม่เจอ marker FIG:{k}"
        s = f()
        assert "$" not in s
        html = pat.sub(lambda m: m.group(1) + s + m.group(2), html)
    PAGE.write_text(html, encoding="utf-8")
    print("ฝัง", ", ".join(FIGS))
