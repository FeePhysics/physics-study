# physics-study — Teaching Workspace ฟิสิกส์ทฤษฎี

> **Context file สำหรับ Claude Code** — โหลดอัตโนมัติทุก session ในโฟลเดอร์นี้

## นี่คืออะไร

**Teaching workspace** ไม่ใช่โปรเจกต์ซอฟต์แวร์ — ของที่ผลิตออกมาคือ *บทเรียน*
ไม่ใช่ฟีเจอร์ · ขับด้วยสกิล **`/teach`** (`.claude/skills/teach/`) ซึ่งแบต้อง
พิมพ์เรียกเอง (`disable-model-invocation: true` — Claude เรียกให้ไม่ได้)

- **เจ้าของ:** แบฟี — ครูคณิตศาสตร์ ม.ปลาย · ติว สอวน. คอมพิวเตอร์ · กำลังจะเรียนต่อ **ป.เอก ฟิสิกส์ทฤษฎี**
- **แยกออกมาจาก `rpg-quiz` โดยตั้งใจ** — repo นั้นเป็นเว็บแอปคณิต มี `CLAUDE.md`
  4,000+ บรรทัดที่ไม่เกี่ยวกับการเรียนฟิสิกส์เลย · ปนกันแล้วทั้งคอนเท็กซ์และ
  git history จะสับสน

## โครงไฟล์ (สกิล `teach` เป็นคนกำหนด)

```
MISSION.md            เหตุผลที่อยากเรียน — ground ทุกบทเรียน (สกิลถามก่อนสอนบทแรก)
RESOURCES.md          แหล่งอ้างอิงที่เชื่อถือได้ + คุณภาพของแต่ละแหล่ง
NOTES.md              สไตล์การสอนที่แบชอบ · ข้อควรจำ
lessons/NNNN-*.html   บทเรียนทีละชิ้น — สั้น จบในหนึ่งนั่ง มีชัยชนะหนึ่งอย่าง
reference/*.html      ชีทสรุป/glossary — ของที่กลับมาเปิดซ้ำ (บทเรียนไม่ค่อยถูกเปิดซ้ำ)
learning-records/     บันทึกว่าเรียนอะไรไปแล้ว → ใช้คำนวณ zone of proximal development
assets/               stylesheet + quiz widget ที่ทุกบทเรียนใช้ร่วมกัน
```

## กฎของ workspace นี้

1. **ภาษาไทยเสมอ** — ศัพท์เทคนิคใช้อังกฤษตามเดิม (Lagrangian, gauge, tensor,
   eigenstate) เพราะแบต้องอ่านตำราและ paper ภาษาอังกฤษอยู่แล้ว การแปลศัพท์
   เป็นไทยจะทำให้จำคำที่ใช้จริงไม่ได้
2. **กระชับ ตรงประเด็น** — ไม่ต้องคำนำ ไม่ต้องสรุปซ้ำ
3. ⚠️ **ห้ามเชื่อความรู้ในหัวตัวเอง** (กฎของสกิล: *"Never trust your parametric
   knowledge"*) — ทุกข้ออ้างในบทเรียนต้องมี citation ไปแหล่งจริง · ฟิสิกส์ทฤษฎี
   เป็นสาขาที่ผิดแล้วดูเนียนที่สุด (เครื่องหมาย convention ต่างตำรา · factor ของ
   $2\pi$ · metric signature $(-+++)$ กับ $(+---)$)
4. **คณิตต้องคำนวณจริง ไม่ใช่พิมพ์ผลลัพธ์จากความจำ** — ใช้ `sympy` ตรวจการอนุพันธ์
   /อินทิเกรต/พีชคณิตเทนเซอร์ก่อนเขียนลงบทเรียนเสมอ · แบเป็นครูคณิต จับผิด
   พีชคณิตได้ทันที แต่นั่นไม่ใช่หน้าที่ของคนเรียน
5. **แบเป็นครูคณิต ไม่ใช่มือใหม่คณิต** — ข้ามการปูพื้นแคลคูลัส/พีชคณิตเชิงเส้น/
   ความน่าจะเป็นได้เลย · สิ่งที่ต้องปูคือ **ฟิสิกส์** และ **คณิตของฟิสิกส์**
   (variational calculus, differential geometry, group theory) ที่คณิต ม.ปลาย
   ไม่มี
6. **LaTeX ในบทเรียนใช้ KaTeX ผ่าน CDN** delimiter `$...$` / `$$...$$` (แบชินกับ
   ชุดนี้อยู่แล้วจาก rpg-quiz)

## Environment

```bash
pip install -r requirements.txt   # sympy numpy scipy matplotlib
python -c "import sympy; print(sympy.__version__)"
```

บทเรียนเป็นไฟล์ HTML ล้วน เปิดด้วยเบราว์เซอร์ได้ตรง ๆ ไม่ต้องมี server

## สั่งงานจาก iPad/มือถือ

`claude.ai/code` → เลือก repo `FeePhysics/physics-study` → พิมพ์ `/teach`
· session รันบนคลาวด์ แก้เสร็จเปิด PR ให้ · merge จากแอป GitHub ได้
· ⚠️ **cloud session เห็นเฉพาะของที่ push แล้ว** — ก่อนออกจากคอมให้ commit+push

## Assets ที่มีแล้ว — ใช้ซ้ำก่อนเขียนใหม่เสมอ

| ไฟล์ | ใช้ทำอะไร |
|---|---|
| `assets/lesson.css` | stylesheet ร่วมของทุกบท — Tufte-ish · dark mode · **print styles** (บทเรียนถูกปรินต์จริง) |
| `assets/quiz.js` | แบบฝึกที่ให้ feedback ทันที · แค่วาง markup ไม่ต้องเขียน JS เพิ่ม |
| `assets/_template.html` | โครงบทเรียน — ก๊อปไปเป็น `lessons/NNNN-*.html` แล้วเติมเนื้อ |

**สร้างของใหม่ที่บทถัดไปอาจใช้ซ้ำ → เขียนเป็น component ใน `assets/`** อย่า inline
ลงบทเรียนเดียว (กฎของสกิล: *"Reuse is the default, not the exception"*)

## ประตูตรวจ — รันก่อน commit ทุกครั้ง

```bash
python3 scripts/check_lessons.py
```

ตรวจ: `$`/`$$` เป็นคู่ · ลิงก์ไม่ตาย · ทุกบท link stylesheet + มีกล่อง `.ask` ·
`data-answer` อยู่ในช่วง · quiz มี `.fb` · **ตัวเลือกยาวไม่ต่างกันเกิน 1.6 เท่า**
(ตัวที่ยาวกว่าเพื่อน = เฉลยที่มองเห็นได้โดยไม่ต้องคิด)

ซ้อมทุบแล้วยิงครบ 5/5 เคส และไม่ยิงตอนไฟล์ปกติ

## ⚠️ ข้อจำกัดของ Claude Code บนคลาวด์ (วัดแล้ว 2026-09-06)

**egress proxy บล็อก CDN ทุกเจ้า** — `cdn.jsdelivr.net` และ `cdnjs.cloudflare.com`
ตอบ `connect_rejected` ⇒ **สมการไม่เรนเดอร์ในคลาวด์เสมอ ไม่ว่าบทเรียนจะถูกหรือผิด**

⇒ กฎบ้าน "ตรวจบนเบราว์เซอร์จริง KaTeX error 0" ทำในคลาวด์ไม่ได้ · ที่ทำได้คือ
`check_lessons.py` แล้ว **ให้แบเปิดไฟล์บนเครื่องตัวเองยืนยันอีกชั้น** ·
ห้ามรายงานว่า "ตรวจแล้วบนเบราว์เซอร์จริง" ถ้าตรวจจากคลาวด์

ทดสอบ `quiz.js` ในคลาวด์ยังทำได้ (ไม่พึ่ง CDN) — Chromium ติดตั้งไว้แล้วที่
`/opt/pw-browsers/chromium` ใช้ผ่าน `pip install playwright` (อย่ารัน `playwright install`)

เว็บที่บล็อกด้วย: `goodtheorist.science`, `webspace.science.uu.nl` (ดู `RESOURCES.md`)
