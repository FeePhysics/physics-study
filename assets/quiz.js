/* quiz.js — แบบฝึกที่ให้ feedback ทันที
   สกิลกำหนดว่า skill acquisition ต้องมี feedback loop ที่แน่นที่สุดเท่าที่ทำได้
   ⇒ ตอบปุ๊บรู้ผลปั๊บ ไม่มีปุ่ม "ส่งคำตอบ" ไม่ต้องรอ

   วิธีใช้ — ไม่ต้องเขียน JS เพิ่มในบทเรียน แค่วาง markup นี้:

     <div class="quiz" data-answer="1">
       <h3>คำถาม…</h3>
       <div class="choices">
         <button>ตัวเลือกแรก</button>
         <button>ตัวเลือกที่สอง</button>
       </div>
       <p class="fb">อธิบายว่าทำไมข้อนี้ถูก และทำไมข้ออื่นผิด</p>
     </div>

   data-answer = ดัชนีของตัวเลือกที่ถูก เริ่มที่ 0

   ⚠️ กฎเนื้อหาที่ CSS/JS บังคับให้ไม่ได้ (สกิลกำหนด แต่คนเขียนต้องคุมเอง):
      ทุกตัวเลือกต้องยาวเท่ากัน — จำนวนคำและจำนวนอักขระ
      ตัวเลือกที่ยาวกว่าเพื่อน = เฉลยที่มองเห็นได้โดยไม่ต้องคิด */

document.addEventListener('click', function (e) {
  var btn = e.target.closest('.quiz .choices button');
  if (!btn) return;

  var quiz = btn.closest('.quiz');
  var buttons = Array.prototype.slice.call(quiz.querySelectorAll('.choices button'));
  var answer = parseInt(quiz.dataset.answer, 10);
  var picked = buttons.indexOf(btn);

  buttons.forEach(function (b) { b.disabled = true; });
  buttons[answer].classList.add('right');
  if (picked !== answer) btn.classList.add('wrong');

  var fb = quiz.querySelector('.fb');
  if (fb) fb.classList.add('show');

  quiz.dataset.result = (picked === answer) ? 'right' : 'wrong';
});
