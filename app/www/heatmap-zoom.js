/* Pan/zoom navigation for the Heatmap tab image.
 *
 * Mirrors the interactive report's controls on the (static) heatmap
 * image: scroll pans up/down and Shift+scroll left/right (the wrapper's
 * native overflow scrolling), Cmd/Ctrl+scroll zooms at the cursor, a
 * drag selects a region to zoom into, and double-click or Escape
 * resets. One set of delegated listeners: the image is re-rendered by
 * Shiny whenever the figure changes, so per-element bindings would not
 * survive.
 *
 * Zooming works by widening the image (native scrolling then provides
 * the panning); the zoom factor lives in img.dataset.rdZoom.
 */
(function () {
  'use strict';

  var MAX_ZOOM = 12;
  var drag = null; // {wrap, img, startX, startY, box}

  function wrapOf(target) {
    return target && target.closest && target.closest('.rd-heatmap-wrap');
  }

  function imgOf(wrap) {
    return wrap ? wrap.querySelector('img') : null;
  }

  function zoomOf(img) {
    var z = parseFloat(img.dataset.rdZoom || '1');
    return isFinite(z) && z > 0 ? z : 1;
  }

  function applyZoom(wrap, img, newZoom, anchorContentX, anchorContentY, viewX, viewY) {
    var oldZoom = zoomOf(img);
    newZoom = Math.min(MAX_ZOOM, Math.max(1, newZoom));
    var base = img.getBoundingClientRect().width / oldZoom;
    if (newZoom === 1) {
      img.style.width = '';
      img.style.maxWidth = '';
    } else {
      img.style.maxWidth = 'none';
      img.style.width = base * newZoom + 'px';
    }
    img.dataset.rdZoom = String(newZoom);
    var scale = newZoom / oldZoom;
    wrap.scrollLeft = anchorContentX * scale - viewX;
    wrap.scrollTop = anchorContentY * scale - viewY;
  }

  function reset(wrap, img) {
    img.style.width = '';
    img.style.maxWidth = '';
    img.dataset.rdZoom = '1';
    wrap.scrollLeft = 0;
    wrap.scrollTop = 0;
  }

  document.addEventListener(
    'wheel',
    function (ev) {
      var wrap = wrapOf(ev.target);
      if (!wrap || !(ev.ctrlKey || ev.metaKey)) {
        return; // plain / Shift+scroll: the wrapper's native panning
      }
      var img = imgOf(wrap);
      if (!img) return;
      ev.preventDefault();
      var rect = wrap.getBoundingClientRect();
      var viewX = ev.clientX - rect.left;
      var viewY = ev.clientY - rect.top;
      var factor = Math.exp(-ev.deltaY * 0.0015);
      applyZoom(
        wrap,
        img,
        zoomOf(img) * factor,
        wrap.scrollLeft + viewX,
        wrap.scrollTop + viewY,
        viewX,
        viewY
      );
    },
    { passive: false }
  );

  document.addEventListener('mousedown', function (ev) {
    var wrap = wrapOf(ev.target);
    if (!wrap || ev.button !== 0) return;
    var img = imgOf(wrap);
    if (!img) return;
    ev.preventDefault(); // no native image ghost-drag
    var rect = wrap.getBoundingClientRect();
    drag = {
      wrap: wrap,
      img: img,
      startX: wrap.scrollLeft + ev.clientX - rect.left,
      startY: wrap.scrollTop + ev.clientY - rect.top,
      box: null,
    };
  });

  document.addEventListener('mousemove', function (ev) {
    if (!drag) return;
    var rect = drag.wrap.getBoundingClientRect();
    var x = drag.wrap.scrollLeft + ev.clientX - rect.left;
    var y = drag.wrap.scrollTop + ev.clientY - rect.top;
    if (!drag.box) {
      if (Math.abs(x - drag.startX) < 4 && Math.abs(y - drag.startY) < 4) {
        return; // ignore jitter until it reads as a drag
      }
      drag.box = document.createElement('div');
      drag.box.className = 'rd-hm-selbox';
      drag.wrap.appendChild(drag.box);
    }
    var left = Math.min(x, drag.startX);
    var top = Math.min(y, drag.startY);
    drag.box.style.left = left + 'px';
    drag.box.style.top = top + 'px';
    drag.box.style.width = Math.abs(x - drag.startX) + 'px';
    drag.box.style.height = Math.abs(y - drag.startY) + 'px';
  });

  document.addEventListener('mouseup', function (ev) {
    if (!drag) return;
    var d = drag;
    drag = null;
    if (!d.box) return;
    var selLeft = parseFloat(d.box.style.left);
    var selTop = parseFloat(d.box.style.top);
    var selW = parseFloat(d.box.style.width);
    var selH = parseFloat(d.box.style.height);
    d.box.remove();
    if (selW < 12 || selH < 12) return; // too small to mean "zoom here"
    var factor = Math.min(
      d.wrap.clientWidth / selW,
      d.wrap.clientHeight / selH
    );
    var newZoom = Math.min(MAX_ZOOM, zoomOf(d.img) * factor);
    // Anchor the selection centre to the viewport centre.
    applyZoom(
      d.wrap,
      d.img,
      newZoom,
      selLeft + selW / 2,
      selTop + selH / 2,
      d.wrap.clientWidth / 2,
      d.wrap.clientHeight / 2
    );
  });

  document.addEventListener('dblclick', function (ev) {
    var wrap = wrapOf(ev.target);
    var img = imgOf(wrap);
    if (wrap && img) reset(wrap, img);
  });

  document.addEventListener('keydown', function (ev) {
    if (ev.key !== 'Escape') return;
    document.querySelectorAll('.rd-heatmap-wrap').forEach(function (wrap) {
      var img = imgOf(wrap);
      if (img && zoomOf(img) !== 1) reset(wrap, img);
    });
  });
})();
