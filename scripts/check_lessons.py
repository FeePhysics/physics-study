#!/usr/bin/env python3
"""ประตูตรวจบทเรียน — รันก่อน commit ทุกครั้ง

ทำไมต้องมี: กฎบ้านคือ "ตรวจบนเบราว์เซอร์จริง KaTeX error 0" แต่ Claude Code
บนคลาวด์ทำแบบนั้นไม่ได้ — egress proxy บล็อก CDN ทุกเจ้า (jsdelivr, cdnjs)
⇒ สมการไม่เรนเดอร์ในคลาวด์เสมอ ไม่ว่าบทเรียนจะถูกหรือผิด

ตัวนี้จึงตรวจสิ่งที่ตรวจได้โดยไม่ต้องเรนเดอร์ · ไม่ได้แทนการเปิดดูด้วยตา
บนเครื่อง แต่จับความพลาดที่เกิดบ่อยที่สุดได้ก่อนถึงมือคนเรียน

    python3 scripts/check_lessons.py
"""
import re, sys, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
problems = []


def check(f: pathlib.Path, html: str):
    def bad(msg): problems.append(f"{f.relative_to(ROOT)}: {msg}")

    # ── $ ต้องเป็นคู่ ──────────────────────────────────────
    # $ เดี่ยวทำให้ KaTeX กินข้อความที่เหลือทั้งย่อหน้าเข้าไปในโหมดคณิต
    # นับเฉพาะ $ ที่ไม่ได้ถูก escape และไม่ได้อยู่ใน $$
    body = re.sub(r'\$\$.*?\$\$', '', html, flags=re.S)      # ตัดบล็อกคู่ออกก่อน
    if len(re.findall(r'(?<!\\)\$', body)) % 2:
        bad("$ ไม่เป็นคู่ — สมการจะกินข้อความที่เหลือ")
    if html.count('$$') % 2:
        bad("$$ ไม่เป็นคู่")

    # ── $ ห้ามอยู่ใน <svg> ────────────────────────────────
    # renderMathInElement ไม่เดินเข้าไปใน SVG ⇒ $...$ ตรงนั้นขึ้นเป็น LaTeX ดิบบนจอ
    # โดยไม่มี error อะไรฟ้อง · ใช้ยูนิโคด (α β θ − ⁻¹) แทน
    for svg in re.findall(r'<svg.*?</svg>', html, flags=re.S):
        if '$' in svg:
            bad("มี $ อยู่ใน <svg> — KaTeX ไม่เรนเดอร์ในนั้น จะขึ้น LaTeX ดิบ")

    # ── ไฟล์ที่ลิงก์ต้องมีอยู่จริง ─────────────────────────
    for rel in re.findall(r'(?:href|src)="((?!https?:|#|mailto:)[^"]+)"', html):
        target = (f.parent / rel.split('#')[0]).resolve()
        if rel.endswith('/'):
            continue                                          # ลิงก์ไปโฟลเดอร์ ปล่อยผ่าน
        if not target.exists():
            bad(f"ลิงก์ตาย: {rel}")

    # ── โครงที่ทุกบทต้องมี ────────────────────────────────
    if 'assets/lesson.css' not in html:
        bad("ไม่ได้ link assets/lesson.css — บทนี้จะหน้าตาไม่เหมือนบทอื่น")
    if 'class="ask"' not in html:
        bad('ไม่มีกล่อง .ask — สกิลกำหนดว่าทุกบทต้องเตือนให้ผู้เรียนถามครูได้')
    if not re.search(r'<h1[ >]', html):
        bad("ไม่มี <h1>")

    # ── แบบฝึก ────────────────────────────────────────────
    for m in re.finditer(r'<div class="quiz"([^>]*)>(.*?)</div>\s*(?=<h|<div class="ask"|</body)',
                         html, flags=re.S):
        attrs, block = m.group(1), m.group(2)
        am = re.search(r'data-answer="(\d+)"', attrs)
        n = len(re.findall(r'<button', block))
        if not am:
            bad("quiz ไม่มี data-answer — คลิกแล้วไม่มีอะไรเกิดขึ้น")
            continue
        if not 0 <= int(am.group(1)) < n:
            bad(f"quiz: data-answer={am.group(1)} แต่มี {n} ตัวเลือก")
        if 'class="fb"' not in block:
            bad("quiz ไม่มี .fb — ตอบผิดแล้วไม่รู้ว่าผิดตรงไหน")

        # ตัวเลือกต้องยาวพอ ๆ กัน — ตัวที่ยาวกว่าเพื่อนคือเฉลยที่มองเห็นได้
        lens = [len(re.sub(r'<[^>]+>', '', c).strip())
                for c in re.findall(r'<button[^>]*>(.*?)</button>', block, flags=re.S)]
        if lens and max(lens) > 1.6 * min(lens):
            bad(f"quiz: ตัวเลือกยาวไม่เท่ากัน {lens} — ความยาวกลายเป็นคำใบ้")


targets = sorted((ROOT / 'lessons').glob('*.html')) + sorted((ROOT / 'reference').glob('*.html'))
tpl = ROOT / 'assets' / '_template.html'
if tpl.exists():
    targets.append(tpl)

for f in targets:
    html = f.read_text(encoding='utf-8')
    if f == tpl:
        html = html.replace('href="../', 'href="' + str(ROOT) + '/')   # template อยู่คนละชั้น
        html = html.replace('src="../', 'src="' + str(ROOT) + '/')
    check(f, html)

if problems:
    print(f"❌ เจอ {len(problems)} จุด\n")
    for p in problems:
        print("  •", p)
    sys.exit(1)

print(f"✅ ผ่าน — ตรวจ {len(targets)} ไฟล์")
