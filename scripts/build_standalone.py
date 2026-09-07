#!/usr/bin/env python3
"""แปลงบทเรียนเป็นไฟล์เดียวจบ — เปิดที่ไหนก็ได้ ไม่ต้องมีโฟลเดอร์ข้าง ๆ

    python3 scripts/build_standalone.py lessons/0001-action-principle.html
    python3 scripts/build_standalone.py --all          # ทุกบท + reference

ทำไมต้องมี: ไฟล์ในรีโปลิงก์ `../assets/lesson.css` แบบ relative ⇒ ส่งไฟล์เดียว
ให้ใครหรืออัปโหลดขึ้นเว็บ จะไม่มีสไตล์และ quiz ไม่ทำงาน · ตัวนี้ฝัง CSS + JS
ลงในไฟล์ และ**แปลง LaTeX เป็น MathML** ให้ด้วย ⇒ สมการขึ้นโดยไม่ต้องโหลด KaTeX
จาก CDN เลย (บางที่บล็อก CDN — ดู CLAUDE.md)

ผลลัพธ์ไปที่ build/ ซึ่งไม่เข้า git · **ห้ามแก้ไฟล์ใน build/ โดยตรง**
มันถูกสร้างใหม่ทับทุกครั้ง — แก้ที่ lessons/ แล้วรันใหม่
"""
import pathlib, re, sys

try:
    from latex2mathml.converter import convert as tex2mml
except ImportError:
    sys.exit("ต้องติดตั้งก่อน: pip install latex2mathml")

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "build"

FONTS = ('<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>\n'
         '<link rel="stylesheet" href="https://fonts.googleapis.com/css2?'
         'family=Sarabun:wght@400;600&family=Source+Serif+4:opsz,wght@8..60,400;8..60,600'
         '&family=IBM+Plex+Mono:wght@400&display=swap">')


def to_mathml(html: str) -> str:
    """แทน $$...$$ และ $...$ ด้วย MathML

    ทำ $$ ก่อนเสมอ ไม่งั้นตัวจับ $ เดี่ยวจะกินครึ่งหนึ่งของ $$ ไปก่อน
    ⚠️ ห้ามแตะข้างใน <svg> — MathML ในนั้นไม่เรนเดอร์ (กฎเดียวกับที่ check_lessons บังคับ)
    """
    svgs = []
    def stash(m):
        svgs.append(m.group(0)); return f"\x00SVG{len(svgs)-1}\x00"
    html = re.sub(r"<svg.*?</svg>", stash, html, flags=re.S)

    def block(m):
        return f'<math xmlns="http://www.w3.org/1998/Math/MathML" display="block">' \
               + tex2mml(m.group(1).strip()).split(">", 1)[1]
    html = re.sub(r"\$\$(.+?)\$\$", block, html, flags=re.S)
    html = re.sub(r"(?<!\\)\$([^$]+?)\$", lambda m: tex2mml(m.group(1)), html)

    for i, s in enumerate(svgs):
        html = html.replace(f"\x00SVG{i}\x00", s)
    return html


def build(src: pathlib.Path) -> pathlib.Path:
    html = src.read_text(encoding="utf-8")
    css = (ROOT / "assets" / "lesson.css").read_text(encoding="utf-8")
    js = (ROOT / "assets" / "quiz.js").read_text(encoding="utf-8")

    html = to_mathml(html)

    # ตัด <link>/<script> ที่ชี้ไฟล์ข้าง ๆ และ CDN ออกให้หมด แล้วฝังของจริงแทน
    html = re.sub(r'\s*<link rel="stylesheet" href="\.\./assets/lesson\.css">', "", html)
    html = re.sub(r'\s*<link[^>]*katex[^>]*>', "", html)
    html = re.sub(r'\s*<script[^>]*katex[^>]*>\s*</script>', "", html, flags=re.S)
    html = re.sub(r'\s*<script[^>]*quiz\.js[^>]*>\s*</script>', "", html)

    html = html.replace("</head>", f"{FONTS}\n<style>\n{css}\n</style>\n"
                                   f"<script>\n{js}\n</script>\n</head>", 1)

    # ลิงก์ไปไฟล์อื่นในรีโปใช้ไม่ได้เมื่อไฟล์ยืนอยู่เดี่ยว ๆ ⇒ ตัดเหลือข้อความ
    html = re.sub(r'<a href="\.\./[^"]+">(.*?)</a>', r"\1", html, flags=re.S)
    html = re.sub(r'<a href="\./">(.*?)</a>', r"\1", html, flags=re.S)

    OUT.mkdir(exist_ok=True)
    dst = OUT / src.name
    dst.write_text(html, encoding="utf-8")

    # ฉบับสำหรับ publish เป็น Artifact — ตัวห่อ <html>/<head>/<body> ถูกใส่ให้ตอน
    # publish ⇒ ไฟล์ต้องมีแต่เนื้อหา ไม่งั้นแท็กซ้อนกัน
    inner = re.search(r"<body[^>]*>(.*)</body>", html, flags=re.S).group(1)
    head = "\n".join(re.findall(r"<title>.*?</title>|<link rel=\"stylesheet\"[^>]*>|"
                                r"<link rel=\"preconnect\"[^>]*>|<style>.*?</style>|"
                                r"<script>.*?</script>", html, flags=re.S))
    (OUT / ("artifact-" + src.name)).write_text(head + "\n" + inner.strip(), encoding="utf-8")
    return dst


if __name__ == "__main__":
    args = sys.argv[1:]
    if args == ["--all"]:
        srcs = sorted((ROOT/"lessons").glob("*.html")) + sorted((ROOT/"reference").glob("*.html"))
    elif args:
        srcs = [ROOT / a for a in args]
    else:
        sys.exit(__doc__)

    for s in srcs:
        d = build(s)
        kb = d.stat().st_size / 1024
        left = d.read_text(encoding="utf-8").count("$")
        print(f"✅ {d.relative_to(ROOT)}  {kb:.0f} KB  · $ ที่เหลือ: {left}")
