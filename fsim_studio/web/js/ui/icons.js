// Hand-authored icon set: 24-unit grid, 1.5px stroke, round joins, no fills except where noted.
// Drawn for this product: breadboard holes, beam, instrument traces, projector screen.

const P = {
  overview: '<rect x="3.5" y="4" width="7" height="7" rx="1"/><rect x="13.5" y="4" width="7" height="4.5" rx="1"/><rect x="13.5" y="11.5" width="7" height="8.5" rx="1"/><rect x="3.5" y="14" width="7" height="6" rx="1"/>',
  design: '<circle cx="5" cy="6" r=".9"/><circle cx="12" cy="6" r=".9"/><circle cx="19" cy="6" r=".9"/><circle cx="5" cy="18" r=".9"/><circle cx="19" cy="18" r=".9"/><path d="M2.5 12h4.2M17.3 12h4.2"/><rect x="6.7" y="9" width="4" height="6" rx=".6"/><rect x="13.3" y="9" width="4" height="6" rx=".6"/><path d="M10.7 12h2.6"/>',
  compare: '<path d="M3 18c3-0.2 4.5-9 7.5-9s4 6 10.5 6.5"/><path d="M3 14.5c3.5 0 5-7.5 8.5-7.5S17 11 21 11" stroke-dasharray="2 2.2"/><path d="M3 20.5h18"/>',
  results: '<path d="M4 3.5v17h16.5"/><circle cx="8.5" cy="15.5" r="1.2"/><circle cx="12" cy="11" r="1.2"/><circle cx="15.5" cy="13.5" r="1.2"/><circle cx="18.5" cy="7" r="1.2"/>',
  explain: '<circle cx="10.5" cy="10.5" r="6"/><path d="M15 15l5.5 5.5"/><path d="M7.5 11.5c1-2.5 2-2.5 3 0s2 2.5 3 0"/>',
  library: '<path d="M6 3.8h11.5a.9.9 0 0 1 .9.9v12.9H7.1A1.6 1.6 0 0 0 5.5 19.2"/><path d="M5.5 19.2V5.4A1.6 1.6 0 0 1 7.1 3.8"/><path d="M5.5 19.2a1.6 1.6 0 0 0 1.6 1.6h11.3v-3.2"/><path d="M9.5 3.8v6.4l1.6-1.2 1.6 1.2V3.8"/>',
  story: '<rect x="3" y="4" width="18" height="11.5" rx="1"/><path d="M12 15.5v3.5M8 20.5h8"/><path d="M7 12l3-3 2.5 2 4-4"/>',
  sun: '<circle cx="12" cy="12" r="4"/><path d="M12 2.8v2.4M12 18.8v2.4M2.8 12h2.4M18.8 12h2.4M5.5 5.5l1.7 1.7M16.8 16.8l1.7 1.7M5.5 18.5l1.7-1.7M16.8 7.2l1.7-1.7"/>',
  moon: '<path d="M19.5 14.5A8 8 0 0 1 9.5 4.5a8 8 0 1 0 10 10z"/>',
  present: '<rect x="2.5" y="4" width="19" height="12.5" rx="1"/><path d="M12 16.5V20M8 21h8"/><path d="M10 8l4 2.3-4 2.3z"/>',
  close: '<path d="M6 6l12 12M18 6L6 18"/>',
  run: '<path d="M7.5 5l11 7-11 7z"/>',
  stop: '<rect x="6.5" y="6.5" width="11" height="11" rx="1"/>',
  pin: '<path d="M9 3.5h6l-1 5 3 3H7l3-3z"/><path d="M12 11.5v9"/>',
  table: '<rect x="3.5" y="5" width="17" height="14" rx="1"/><path d="M3.5 9.5h17M3.5 14h17M9.5 9.5V19"/>',
  chart: '<path d="M4 4v16h16"/><path d="M6.5 16c3-1 4-8 7-8s3 5 5.5 5"/>',
  download: '<path d="M12 4v11M7.5 10.5L12 15l4.5-4.5"/><path d="M4.5 19.5h15"/>',
  chevronDown: '<path d="M6.5 9.5l5.5 5.5 5.5-5.5"/>',
  chevronRight: '<path d="M9.5 6.5l5.5 5.5-5.5 5.5"/>',
  info: '<circle cx="12" cy="12" r="8.5"/><path d="M12 11v5.5M12 7.6v.4"/>',
  range: '<path d="M7 6H4.5v12H7M17 6h2.5v12H17"/><path d="M9 12h6"/>',
  warn: '<path d="M12 4l9 15.5H3z"/><path d="M12 10v4.5M12 17v.3"/>',
  check: '<path d="M5 12.5l4.5 4.5L19 7.5"/>',
  cross: '<path d="M6.5 6.5l11 11M17.5 6.5l-11 11"/>',
  arrowUp: '<path d="M12 19V5M6.5 10.5L12 5l5.5 5.5"/>',
  arrowDown: '<path d="M12 5v14M6.5 13.5L12 19l5.5-5.5"/>',
  minus: '<path d="M6 12h12"/>',
  retry: '<path d="M19 12a7 7 0 1 1-2.1-5"/><path d="M19 4.5V8h-3.5"/>',
  folder: '<path d="M3.5 6.5h6l2 2h9v10h-17z"/>',
  menu: '<path d="M4 7h16M4 12h16M4 17h16"/>',
};

export function icon(name, { size = 20, label = null, cls = "" } = {}) {
  const span = document.createElement("span");
  span.className = `ico ${cls}`.trim();
  const aria = label ? `role="img" aria-label="${label}"` : 'aria-hidden="true"';
  span.innerHTML = `<svg viewBox="0 0 24 24" width="${size}" height="${size}" fill="none" stroke="currentColor" stroke-width="1.5" stroke-linecap="round" stroke-linejoin="round" ${aria}>${P[name] || ""}</svg>`;
  return span;
}
