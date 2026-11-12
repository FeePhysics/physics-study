# RESOURCES

แหล่งอ้างอิงสำหรับ ground การสอน · **ทุกข้ออ้างในบทเรียนต้อง cite กลับมาที่นี่**

> สถานะ: ตั้งต้น — ยังไม่ได้ประเมินเชิงลึก · จะโตขึ้นเมื่อเริ่มสอนแต่ละหัวข้อ
> และต้องระบุ *คุณภาพ* กับ *เหมาะกับใคร* ไม่ใช่แค่ชื่อกับลิงก์

## ⚠️ เข้าถึงไม่ได้จาก Claude Code บนคลาวด์

egress proxy บล็อกไว้ — ถ้าจะให้ Claude อ่าน ต้องรันจากเครื่องแบ หรือ save
เนื้อหามาวางใน repo

| แหล่ง | สถานะ |
|---|---|
| `goodtheorist.science` | ❌ blocked (ตรวจ 2026-09-06) |
| `webspace.science.uu.nl/~hooft101/` | ❌ blocked (ตรวจ 2026-09-06) |

## แผนที่ทั้งเส้นทาง (อ่านก่อนตัดสินใจว่าจะเรียงลำดับยังไง)

### How to become a GOOD Theoretical Physicist — Gerard 't Hooft
- https://www.goodtheorist.science/ · mirror: https://webspace.science.uu.nl/~hooft101/theorist.html
- **คืออะไร:** โนเบล 1999 เรียงลำดับหัวข้อทั้งหมดที่ต้องผ่าน ตั้งแต่คณิตพื้นฐาน
  ไปจนถึง QFT/GR/string พร้อมลิงก์ lecture notes ฟรีของแต่ละหัวข้อ
- **ทำไมสำคัญ:** เป็น *แผนที่* ไม่ใช่ตำรา — ใช้ตัดสินว่า MISSION ของแบต้องผ่าน
  อะไรบ้างและลำดับไหน · เขียนสำหรับคนที่เรียนเองนอกระบบพอดี
- ⚠️ ยังไม่ได้เปิดอ่านจริง (proxy บล็อก) — **ต้องเปิดอ่านก่อนวางหลักสูตร**

## Paper ของ thesis (อาจารย์ที่ปรึกษาให้มา · 2026-11-11)

### Trotta — *Bayes in the sky* (Contemp. Phys. 49 (2008) 71-104)
- doi:10.1080/00107510802066753 · 34 หน้า · บทความทบทวนระดับบัณฑิตศึกษา
- **คุณภาพ:** สูง — ผู้เขียนทำงานด้านนี้โดยตรง ไล่จากหลักการถึงงานวิจัย มีตัวอย่าง Gaussian ที่คำนวณมือได้ทุกหัวข้อ
- **จุดอ่อน:** ปี 2008 ⇒ ข้อมูล WMAP3 ไม่ใช่ Planck · ตาราง 3-5 เป็นตัวเลขเก่า ใช้เรียนวิธีคิด อย่าใช้ตัวเลข
- ⚠️ ไฟล์ PDF แบแนบใน session 2026-11-11 · ไม่ได้ commit (หน้าแรกของ Taylor & Francis ห้ามแจกจ่าย) · session ใหม่ต้องแนบใหม่

### Torrado & Lewis — *Cobaya* (arXiv:2005.05290v2, 2021)
- paper ของโปรแกรม Cobaya · เอกสาร https://cobaya.readthedocs.io · โค้ด https://github.com/CobayaSampler/cobaya
- **คุณภาพ:** สูง — Lewis เขียน CAMB/CosmoMC/GetDist เอง · แต่เป็น paper ของซอฟต์แวร์ ⇒ สมมติว่ารู้สถิติแล้ว
- ใช้ตัวอย่างวงแหวน (รูป 2-4) เป็นแบบฝึกที่รันจริงได้ · `pip install cobaya` ในคลาวด์ยังไม่ได้ลอง

## วิดีโอ / คอร์ส

### The Theoretical Minimum — Leonard Susskind (Stanford)
- https://theoreticalminimum.com/
- 6 คอร์สต่อกัน: classical mechanics (Lagrange/Hamilton) → relativity → QM →
  cosmology/black holes · ฟรีทั้งหมด (YouTube + iTunes U)
- **จุดแข็ง:** สอนแบบ "ขั้นต่ำที่ต้องรู้จริง ๆ" ไม่อ้อมค้อม เหมาะกับคนที่คณิตแน่นแล้ว
- **จุดอ่อน:** เป็นวิดีโอ ⇒ ค้นย้อนกลับยาก · ใช้เป็นตัวปูความเข้าใจ แล้วไปเอา
  ความแม่นยำจาก lecture notes

### MIT OpenCourseWare 8.323 — Relativistic Quantum Field Theory I
- https://ocw.mit.edu/courses/8-323-relativistic-quantum-field-theory-i-spring-2023/
- รุ่น 2023 มี **วิดีโอ + problem sets + เฉลย** · รุ่น 2008 มี lecture notes แยก
  (https://ocw.mit.edu/courses/physics/8-323-relativistic-quantum-field-theory-i-spring-2008/lecture-notes/)
- ฟรีภายใต้ Creative Commons
- **มีเฉลย = มี feedback loop** ⇒ ใช้เป็นแหล่งโจทย์ฝึกได้ ไม่ใช่แค่อ่าน

## ยังไม่ได้หา — คิวถัดไป

- [ ] ชุมชนสำหรับถามและทดสอบความเข้าใจ (Physics Stack Exchange · r/AskPhysics ·
      กลุ่มไทย) — สกิลบอกว่า *wisdom* มาจากการโต้ตอบกับคนจริง ไม่ใช่จากบทเรียน
- [ ] ตำราหลักของแต่ละวิชา + ฉบับที่มีเฉลย
- [ ] แหล่งคณิตของฟิสิกส์โดยเฉพาะ (differential geometry, group theory,
      variational calculus) — ส่วนที่คณิต ม.ปลายไม่มี
- [ ] arXiv: หมวดที่ตรงกับหัวข้อวิจัยที่แบสนใจ (รอ MISSION.md)
