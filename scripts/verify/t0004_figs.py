"""รูปประกอบ T4 (SVG เขียนมือ · ไม่มี $ ข้างใน · ข้อความใช้ currentColor)

รัน: python3 scripts/verify/t0004_figs.py  → เขียน SVG ทับระหว่าง <!--FIG:x--> ... <!--/FIG:x--> ในบท

  ex4     ตารางของข้อ 4 เป็นจุดบนกราฟ (ไม่มีเส้น Δχ²)
  warmup  ข้อมูลสี่จุด ± σ (ชุดเดียวกับของเล่น) + เส้น A (H0 = 70) และเส้น B (H0 = 90)
"""
import re
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[2]
LESSON = ROOT / "lessons" / "t0004-chi2.html"
BLUE, ORANGE = "#3b6fb6", "#c8692a"
FONT = 'font-family="Sarabun, sans-serif" fill="currentColor"'

XS = [10, 20, 30, 40]; YS = [700, 1500, 2100, 3600]; SS = [100, 100, 100, 600]


def fig_warmup():
    W, H, L, B, T, R = 400, 300, 58, 48, 12, 385
    X = lambda d: L + d / 45 * (R - L)
    Y = lambda v: H - B - v / 4400 * (H - B - T)
    o = []
    o.append(f'<line x1="{L}" y1="{T}" x2="{L}" y2="{H-B}" stroke="currentColor" opacity=".5"/>')
    o.append(f'<line x1="{L}" y1="{H-B}" x2="{R}" y2="{H-B}" stroke="currentColor" opacity=".5"/>')
    for d in (0, 10, 20, 30, 40):
        o.append(f'<text x="{X(d):.1f}" y="{H-B+16}" font-size="13" text-anchor="middle" {FONT}>{d}</text>')
    for v in (0, 1000, 2000, 3000, 4000):
        o.append(f'<text x="{L-5}" y="{Y(v)+4:.1f}" font-size="13" text-anchor="end" {FONT}>{v}</text>')
    o.append(f'<text x="{(L+R)/2}" y="{H-8}" font-size="14" text-anchor="middle" {FONT}>ระยะทาง d (Mpc)</text>')
    o.append(f'<text transform="translate(14,{(T+H-B)/2}) rotate(-90)" font-size="14" text-anchor="middle" {FONT}>ความเร็ว v (km/s)</text>')
    for k, col, dash, name in ((70, BLUE, "", "A"), (90, ORANGE, ' stroke-dasharray="7 5"', "B")):
        o.append(f'<line x1="{X(0)}" y1="{Y(0)}" x2="{X(45):.1f}" y2="{Y(k*45):.1f}" stroke="{col}" stroke-width="3"{dash}/>')
        lx = 25 if name == "A" else 14            # ป้ายวางกลางเส้น ห่างจากแถบของจุดที่สี่
        o.append(f'<text x="{X(lx):.1f}" y="{Y(k*lx) + (22 if name=="A" else -10):.1f}" font-size="16" font-weight="600" '
                 f'text-anchor="{"start" if name=="A" else "end"}" fill="{col}" font-family="Sarabun, sans-serif">เส้น {name}</text>')
    for x, y, s in zip(XS, YS, SS):
        o.append(f'<path d="M{X(x):.1f},{Y(y-s):.1f} V{Y(y+s):.1f} M{X(x)-5:.1f},{Y(y-s):.1f} h10 M{X(x)-5:.1f},{Y(y+s):.1f} h10" '
                 f'stroke="currentColor" stroke-width="1.6" fill="none"/>')
        o.append(f'<circle cx="{X(x):.1f}" cy="{Y(y):.1f}" r="4.5" fill="currentColor"/>')
    body = "".join(o)
    return (f'<figure><svg role="img" aria-label="ข้อมูลสี่จุดพร้อมแถบความคลาด สามจุดแรกแถบสั้น จุดที่สี่แถบยาวมาก '
            f'เส้น A ผ่านใกล้สามจุดแรก เส้น B ผ่านจุดที่สี่พอดีแต่ห่างสามจุดแรก" viewBox="0 0 {W} {H}" '
            f'xmlns="http://www.w3.org/2000/svg" style="max-width:{W}px">{body}</svg></figure>')


def fig_ex4():
    """ตารางข้อ 4 เป็นจุดบนกราฟ (ไม่ลากเส้น Δχ² ให้ · ไม่บอกคำตอบ)"""
    W, H, L, B, T, R = 400, 230, 50, 46, 12, 385
    X = lambda h: L + (h - 64) / 12 * (R - L)
    Y = lambda c: H - B - c / 8 * (H - B - T)
    o = [f'<line x1="{L}" y1="{T}" x2="{L}" y2="{H-B}" stroke="currentColor" opacity=".5"/>',
         f'<line x1="{L}" y1="{H-B}" x2="{R}" y2="{H-B}" stroke="currentColor" opacity=".5"/>']
    for h in (66, 68, 70, 72, 74):
        o.append(f'<text x="{X(h):.1f}" y="{H-B+16}" font-size="13" text-anchor="middle" {FONT}>{h}</text>')
    for c in (0, 2, 4, 6, 8):
        o.append(f'<text x="{L-6}" y="{Y(c)+4:.1f}" font-size="13" text-anchor="end" {FONT}>{c}</text>')
        o.append(f'<line x1="{L}" y1="{Y(c):.1f}" x2="{R}" y2="{Y(c):.1f}" stroke="currentColor" opacity=".12"/>')
    o.append(f'<text x="{(L+R)/2}" y="{H-8}" font-size="14" text-anchor="middle" {FONT}>H₀ (km/s/Mpc)</text>')
    o.append(f'<text transform="translate(16,{(T+H-B)/2}) rotate(-90)" font-size="14" text-anchor="middle" {FONT}>χ²</text>')
    for h, c in ((66, 7), (68, 4), (70, 3), (72, 4), (74, 7)):
        o.append(f'<circle cx="{X(h):.1f}" cy="{Y(c):.1f}" r="5" fill="{BLUE}"/>')
    return (f'<figure><svg role="img" aria-label="จุดห้าจุดจากตาราง H0 66 68 70 72 74 มี chi-square 7 4 3 4 7 เป็นรูปพาราโบลา" '
            f'viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg" style="max-width:{W}px">{"".join(o)}</svg></figure>')


FIGS = {"warmup": fig_warmup, "ex4": fig_ex4}

if __name__ == "__main__":
    html = LESSON.read_text(encoding="utf-8")
    for k, f in FIGS.items():
        pat = re.compile(rf"(<!--FIG:{k}-->).*?(<!--/FIG:{k}-->)", re.S)
        assert pat.search(html), f"ไม่เจอ marker FIG:{k}"
        s = f()
        assert "$" not in s
        html = pat.sub(lambda m: m.group(1) + s + m.group(2), html)
    LESSON.write_text(html, encoding="utf-8")
    print("ฝัง", ", ".join(FIGS))
