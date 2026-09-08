/* Row-click selection for the cluster-assignment table.
 *
 * One delegated listener (the table is re-rendered with each new result,
 * so per-row bindings would not survive).  Clicking a row toggles its
 * whole cluster; the selected cluster names go to Shiny as the
 * 'cluster_select' input, and the server forwards them to the report
 * iframe (bridge.js, 'rd_highlight_clusters') to dim non-members.
 *
 * bridge.js dispatches the 'rd-cluster-clear' DOM event when the server
 * clears the highlight (new result / new assignments), so the table's
 * row styling cannot go stale.
 */
(function () {
  'use strict';

  var selected = {};

  function applyRowClasses() {
    document.querySelectorAll('tr.rd-cluster-row').forEach(function (row) {
      row.classList.toggle(
        'rd-cluster-selected',
        !!selected[row.getAttribute('data-cluster')]
      );
    });
  }

  function announce() {
    if (window.Shiny && typeof window.Shiny.setInputValue === 'function') {
      window.Shiny.setInputValue(
        'cluster_select',
        { clusters: Object.keys(selected) },
        { priority: 'event' }
      );
    }
  }

  document.addEventListener('click', function (ev) {
    var row =
      ev.target && ev.target.closest && ev.target.closest('tr.rd-cluster-row');
    if (!row) return;
    var name = row.getAttribute('data-cluster');
    if (!name) return;
    if (selected[name]) {
      delete selected[name];
    } else {
      selected[name] = true;
    }
    applyRowClasses();
    announce();
  });

  document.addEventListener('rd-cluster-clear', function () {
    selected = {};
    applyRowClasses();
  });
})();
