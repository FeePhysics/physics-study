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

   ── แบบพิมพ์คำตอบเอง (ไม่มีตัวเลือกให้เดา) ──────────────

     <div class="quiz" data-open data-answer="-k*y/m" data-alt="-(k/m)y">
       <h3>คำถาม…</h3>
       <div class="ans"><input type="text" placeholder="พิมพ์คำตอบ"><button>ตรวจ</button></div>
       <p class="fb">อธิบาย…</p>
     </div>

   ปรนัยเดาถูกได้ 33–50% ⇒ คะแนนเต็มไม่ได้พิสูจน์ว่าเข้าใจ (ผู้เรียนรายงานเองว่า
   "ตอบถูกเพราะเดา") · แบบพิมพ์เองตัดช่องนั้นทิ้งทั้งหมด

   การตรวจยืดหยุ่นเรื่อง *รูปแบบ* แต่ไม่ยืดหยุ่นเรื่อง *ค่า*:
   ตัวเลขเทียบเชิงค่า (รับ 4/15 · 0.2667 · .267) · นิพจน์ตัด * และช่องว่างก่อนเทียบ
   · `data-alt` ใส่รูปแบบอื่นที่ยอมรับได้ คั่นด้วย |
   · พิมพ์อะไรที่อ่านไม่ออกเลย = **ไม่นับผิด** บอกให้ลองใหม่

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
        ? '<div class="score-btns">' +
          '<button type="button" id="copy">📋 คัดลอกผล</button>' +
          '<button type="button" id="again">↻ เริ่มใหม่ทั้งชุด</button></div>' : '');
  }

  /* ── ตรวจคำตอบที่พิมพ์เอง ───────────────────────────── */
  function norm(v) {
    return String(v).trim().toLowerCase()
      .replace(/[−–—]/g, '-')          // ขีดยูนิโคด → ลบธรรมดา
      .replace(/\*\*/g, '^')
      .replace(/[\s*·]/g, '')          // ช่องว่างและเครื่องหมายคูณไม่สำคัญ
      .replace(/[()]/g, '');
  }
  function num(v) {                     // อ่านเป็นตัวเลข รับเศษส่วน a/b
    var s = norm(v), m = s.match(/^(-?\d*\.?\d+)\/(-?\d*\.?\d+)$/);
    if (m) return parseFloat(m[1]) / parseFloat(m[2]);
    return /^-?\d*\.?\d+$/.test(s) ? parseFloat(s) : null;
  }
  function judge(quiz, typed) {
    if (!typed.trim()) return 'empty';
    var keys = [quiz.dataset.answer].concat((quiz.dataset.alt || '').split('|').filter(Boolean));
    var got = num(typed);
    for (var i = 0; i < keys.length; i++) {
      var want = num(keys[i]);
      if (got !== null && want !== null) {
        // ยอมให้ปัดเศษได้ 0.5% — โจทย์บอกว่าตอบเป็นทศนิยมก็ได้ แล้ว 4/15 → 0.2666
        // ต่างจากค่าจริง 6.7e-5 ซึ่งเกิน tolerance แบบสัมบูรณ์ที่แคบไป
        var tol = Math.max(1e-9, Math.abs(want) * 5e-3);
        if (Math.abs(got - want) <= tol) return 'right';
      }
      else if (norm(typed) === norm(keys[i])) return 'right';
    }
    return 'wrong';
  }

  document.addEventListener('click', function (e) {
    var ansBtn = e.target.closest('.quiz[data-open] .ans button');
    if (ansBtn) {
      var q = ansBtn.closest('.quiz');
      if (q.dataset.result) return;
      var input = q.querySelector('.ans input');
      var verdict = judge(q, input.value);

      if (verdict === 'empty') {        // ยังไม่ได้ตอบ ≠ ตอบผิด
        input.focus();
        input.placeholder = 'พิมพ์คำตอบก่อนกดตรวจ';
        return;
      }
      input.readOnly = true;
      input.classList.add(verdict);
      ansBtn.disabled = true;
      var fb = q.querySelector('.fb');
      if (fb) fb.classList.add('show');
      if (verdict === 'wrong') {
        var key = document.createElement('p');
        key.className = 'ans-key';
        key.textContent = 'คำตอบ: ' + q.dataset.answer;
        fb.parentNode.insertBefore(key, fb);
      }
      q.dataset.result = verdict;
      board();
      return;
    }

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

    /* คัดลอกผลเป็นข้อความสั้น ๆ ให้วางกลับไปบอกครูได้ในคลิกเดียว
       — ไม่งั้นต้องไล่ถามว่าผิดข้อไหน ซึ่งเป็นงานที่ไม่ควรตกเป็นของผู้เรียน */
    if (e.target.id === 'copy') {
      var all = [].slice.call(document.querySelectorAll('.quiz'));
      var wrong = [];
      all.forEach(function (q, i) { if (q.dataset.result === 'wrong') wrong.push(i + 1); });
      var text = 'ถูก ' + (all.length - wrong.length) + '/' + all.length +
                 (wrong.length ? ' · ผิดข้อ ' + wrong.join(', ') : ' · ถูกหมด');

      var done = function () {
        e.target.textContent = '✓ คัดลอกแล้ว';
        setTimeout(function () { e.target.textContent = '📋 คัดลอกผล'; }, 2000);
      };
      // navigator.clipboard ใช้ไม่ได้บน http และในบาง iframe ⇒ ต้องมีทางสำรอง
      if (navigator.clipboard && window.isSecureContext) {
        navigator.clipboard.writeText(text).then(done, function () { fallback(text, done); });
      } else { fallback(text, done); }
      return;
    }

    if (e.target.id === 'again') {                     // ล้างทั้งชุดเพื่อทำซ้ำ
      document.querySelectorAll('.ans-key').forEach(function (n) { n.remove(); });
      document.querySelectorAll('.quiz[data-open]').forEach(function (q) {
        var i = q.querySelector('.ans input');
        i.value = ''; i.readOnly = false; i.classList.remove('right', 'wrong');
        q.querySelector('.ans button').disabled = false;
      });
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

  function fallback(text, ok) {
    var ta = document.createElement('textarea');
    ta.value = text; ta.style.position = 'fixed'; ta.style.opacity = '0';
    document.body.appendChild(ta); ta.select();
    try { document.execCommand('copy'); ok(); } catch (err) { window.prompt('คัดลอกข้อความนี้:', text); }
    document.body.removeChild(ta);
  }

  document.addEventListener('keydown', function (e) {   // กด Enter = กดตรวจ
    if (e.key === 'Enter' && e.target.matches('.quiz[data-open] .ans input')) {
      e.preventDefault();
      e.target.parentNode.querySelector('button').click();
    }
  });

  document.addEventListener('DOMContentLoaded', board);
  if (document.readyState !== 'loading') board();
})();
