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

   ตัวนับคะแนน (ไม่บังคับ) — วาง <div id="score"></div> ไว้ที่ไหนก็ได้ในหน้า
   แล้วมันจะอัปเดตเอง พร้อมปุ่มเริ่มใหม่เมื่อทำครบ

   ⚠️ กฎเนื้อหาที่ CSS/JS บังคับให้ไม่ได้ (สกิลกำหนด แต่คนเขียนต้องคุมเอง):
      ทุกตัวเลือกต้องยาวเท่ากัน — จำนวนคำและจำนวนอักขระ
      ตัวเลือกที่ยาวกว่าเพื่อน = เฉลยที่มองเห็นได้โดยไม่ต้องคิด */

(function () {
  function board() {
    var box = document.getElementById('score');
    if (!box) return;
    var all = document.querySelectorAll('.quiz');
    var done = document.querySelectorAll('.quiz[data-result]');
    var ok = document.querySelectorAll('.quiz[data-result="right"]');
    var pct = all.length ? Math.round(done.length / all.length * 100) : 0;

    box.innerHTML =
      '<div class="score-bar"><span style="width:' + pct + '%"></span></div>' +
      '<p class="score-text">ทำแล้ว <b>' + done.length + '</b> จาก ' + all.length +
      ' ข้อ · ถูก <b>' + ok.length + '</b></p>' +
      (done.length === all.length && all.length
        ? '<button type="button" id="again">↻ เริ่มใหม่ทั้งชุด</button>' : '');
  }

  document.addEventListener('click', function (e) {
    var btn = e.target.closest('.quiz .choices button');
    if (btn) {
      var quiz = btn.closest('.quiz');
      if (quiz.dataset.result) return;                 // ตอบไปแล้ว
      var buttons = Array.prototype.slice.call(quiz.querySelectorAll('.choices button'));
      var answer = parseInt(quiz.dataset.answer, 10);
      var picked = buttons.indexOf(btn);

      buttons.forEach(function (b) { b.disabled = true; });
      buttons[answer].classList.add('right');
      if (picked !== answer) btn.classList.add('wrong');

      var fb = quiz.querySelector('.fb');
      if (fb) fb.classList.add('show');

      quiz.dataset.result = (picked === answer) ? 'right' : 'wrong';
      board();
      return;
    }

    if (e.target.id === 'again') {                     // ล้างทั้งชุดเพื่อทำซ้ำ
      document.querySelectorAll('.quiz').forEach(function (q) {
        delete q.dataset.result;
        q.querySelectorAll('.choices button').forEach(function (b) {
          b.disabled = false; b.classList.remove('right', 'wrong');
        });
        var fb = q.querySelector('.fb');
        if (fb) fb.classList.remove('show');
      });
      board();
      window.scrollTo({ top: 0, behavior: 'smooth' });
    }
  });

  document.addEventListener('DOMContentLoaded', board);
  if (document.readyState !== 'loading') board();
})();
