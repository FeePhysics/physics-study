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


_NUMRE = re.compile(r"^[+-]?\d*\.?\d+$")


def _is_num(v: str) -> bool:
    """คำตอบนี้เป็นตัวเลขล้วนไหม (กฎเดียวกับ num() ใน quiz.js · รับเศษส่วน a/b)"""
    t = re.sub(r"[\s*·]", "", str(v).strip().lower()).replace("\u2212", "-")
    if re.match(r"^[+-]?\d*\.?\d+/[+-]?\d*\.?\d+$", t):
        return True
    return bool(_NUMRE.match(t))


def _placeholder_declares_number(ph: str) -> bool:
    """ช่องพิมพ์บอกไหมว่าขอ "ตัวเลข" — บอกตรง ๆ หรือยกตัวอย่างที่เป็นตัวเลขก็ได้"""
    if re.search(r"ตัวเลข|จำนวน", ph):
        return True
    m = re.match(r"\s*เช่น\s+(.+)$", ph)
    if m and _is_num(m.group(1)):
        return True
    # "1 หรือ 2" — ไล่ค่าที่เป็นไปได้ทั้งหมด ชัดกว่าคำว่า "ตัวเลข" ด้วยซ้ำ
    toks = [t for t in re.split(r"[\s,/]+|หรือ", ph) if t]
    return bool(toks) and all(_is_num(t) for t in toks)


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

    # ── บทที่มีแบบฝึกต้องมี <div id="score"> ──────────────
    # quiz.js สร้างปุ่ม "📋 คัดลอกผล" ไว้ข้างใน #score เท่านั้น (board())
    # ไม่มี element นี้ = ปุ่มไม่โผล่ · หน้ายังใช้ได้ปกติ ไม่มี error อะไรฟ้อง
    # แต่ผู้เรียนส่งผลกลับมาไม่ได้ ซึ่งเป็นทั้งหมดของ loop นี้ (เกิดจริงในบทที่ 9)
    if 'class="quiz"' in html and 'id="score"' not in html:
        bad('มีแบบฝึกแต่ไม่มี <div id="score"> — ปุ่ม 📋 คัดลอกผล จะไม่โผล่')

    # ── แบบฝึก — มีสองแบบ: ปรนัย กับ พิมพ์คำตอบเอง ───────
    # ⚠️ เดิมจับด้วย regex ตัวเดียว `<div class="quiz"...>(.*?)</div>` + lookahead
    #    ซึ่ง **ยุบทุกข้อในหน้าเป็นก้อนเดียว** (วัดแล้ว: 7 ข้อ → 1 block)
    #    เพราะ </div> ปิดข้อหนึ่งตามด้วย <div class="quiz"> ของข้อถัดไป ไม่ใช่ <h
    #    ⇒ ด่านรายข้อทุกด่านในนี้ตรวจก้อนรวมมาตลอด และบางด่านยิงไม่ออกเลย
    #    หั่นด้วยตำแหน่งเริ่มของแต่ละข้อแทน
    starts = [mm.start() for mm in re.finditer(r'<div class="quiz"', html)]
    stop = min([i for i in (html.find('<div class="ask"'), html.find('</body')) if i != -1]
               or [len(html)])
    for k, st in enumerate(starts):
        en = starts[k + 1] if k + 1 < len(starts) else stop
        chunk = html[st:en]
        hm = re.match(r'<div class="quiz"([^>]*)>', chunk)
        if not hm:
            continue
        attrs, block = hm.group(1), chunk[hm.end():]

        # ── ขั้นที่คำตอบเป็น "ตัวเลข" ช่องพิมพ์ต้องบอกว่าขอตัวเลข ──
        # ⚠️ ไม่บอก = ผู้เรียนตอบเป็น *นิพจน์* ซึ่งถูกตามที่คำถามเขียน แต่ถูกนับว่าผิด
        #    (เกิดจริง 2026-09-24 · บทที่ 10 ข้อ 4: ถาม "ค่า off-shell" ตอบ ptt-pxx
        #     ซึ่งเป็นคำตอบของข้อ 1-2 ในบทเดียวกันเป๊ะ · ข้อที่ผ่านคือข้อที่ช่องพิมพ์
        #     เขียนว่า "พิมพ์ตัวเลข")
        for em in re.finditer(r'<(?:li|div)\b[^>]*data-answer="([^"]*)"[^>]*>', block):
            a = em.group(1)
            if not _is_num(a):
                continue
            nx = re.search(r'<(?:li|div)\b[^>]*data-answer="', block[em.end():])
            seg = block[em.end():em.end() + (nx.start() if nx else len(block))]
            ph = re.search(r'placeholder="([^"]*)"', seg)
            if ph and not _placeholder_declares_number(ph.group(1)):
                bad(f'คำตอบ {a!r} เป็นตัวเลข แต่ช่องพิมพ์เขียน {ph.group(1)!r}'
                    ' — ต้องบอกว่าขอตัวเลข ไม่งั้นผู้เรียนตอบเป็นนิพจน์')

        if 'data-steps' in attrs:
            # ── ต้องมี <p class="fb"> สรุปท้ายข้อ (หลัง </ol>) ────────
            # check_layout.py บังคับข้อนี้อยู่แล้ว แต่ต้องเปิดเบราว์เซอร์ ~3 นาที
            # ⇒ ตรวจจากตัวไฟล์ตรงนี้ด้วย จะได้รู้ใน 1 วินาที (ลืมมาแล้วสองรอบ)
            if not re.search(r'</ol>\s*<p class="fb"', block):
                bad("quiz[data-steps] ไม่มี <p class=\"fb\"> สรุปท้ายข้อหลัง </ol>")

            # โจทย์ไล่ขั้น: เฉลยอยู่ที่ <li> แต่ละขั้น ไม่ใช่ที่ตัว .quiz
            steps = re.findall(r'<li ([^>]*)>(.*?)</li>', block, flags=re.S)

            if len(steps) < 2:
                bad("quiz[data-steps] มีขั้นเดียว — ใช้ data-open แทน")
            for j, (sa, sb) in enumerate(steps, 1):
                if not re.search(r'data-answer="[^"]+"', sa):
                    bad(f"quiz[data-steps] ขั้นที่ {j} ไม่มี data-answer")
                if '<input' not in sb or '<button' not in sb:
                    bad(f"quiz[data-steps] ขั้นที่ {j} ไม่มีช่องกรอกหรือปุ่มตรวจ")
                if 'class="fb"' not in sb:
                    bad(f"quiz[data-steps] ขั้นที่ {j} ไม่มี .fb — ผิดแล้วไม่รู้ว่าผิดตรงไหน")
            if 'class="choices"' in block:
                bad("quiz[data-steps] มี .choices ปนอยู่ด้วย")
            continue

        am = re.search(r'data-answer="([^"]*)"', attrs)
        if not am:
            bad("quiz ไม่มี data-answer — ตอบแล้วไม่มีอะไรเกิดขึ้น")
            continue
        if 'class="fb"' not in block:
            bad("quiz ไม่มี .fb — ตอบผิดแล้วไม่รู้ว่าผิดตรงไหน")

        if 'data-open' in attrs:
            # แบบพิมพ์เอง: ต้องมีช่องกรอกกับปุ่มตรวจ และเฉลยต้องไม่ว่าง
            if '<input' not in block:
                bad("quiz[data-open] ไม่มีช่องให้พิมพ์คำตอบ")
            if '<button' not in block:
                bad("quiz[data-open] ไม่มีปุ่มตรวจ")
            if not am.group(1).strip():
                bad("quiz[data-open] เฉลยว่าง")
            # ห้ามมี .choices ปนมา — จะกลายเป็นสองกลไกในข้อเดียว
            if 'class="choices"' in block:
                bad("quiz[data-open] มี .choices ปนอยู่ด้วย")
            continue

        # ปรนัย: data-answer ต้องเป็นดัชนีที่อยู่ในช่วง
        n = len(re.findall(r'<button', block))
        if not am.group(1).isdigit():
            bad(f"quiz ปรนัย: data-answer=\"{am.group(1)}\" ต้องเป็นดัชนีตัวเลข"
                " (ถ้าตั้งใจให้พิมพ์คำตอบเอง ต้องใส่ data-open)")
        elif not 0 <= int(am.group(1)) < n:
            bad(f"quiz: data-answer={am.group(1)} แต่มี {n} ตัวเลือก")

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
