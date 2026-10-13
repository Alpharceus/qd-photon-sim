// Screenshot assertions for the browser verifiers (R2 audit 14): a capture must not be blank.
//
//   const { pngStats, assertPng } = require("./png_check.cjs");
//   const st = assertPng(file);   // { ok, reason, width, height, variance, topColourShare, distinctColours }
//
// Dependency-free PNG reader (8-bit, non-interlaced, colour types 0/2/3/4/6, which is what Chromium writes).
// Non-blank means: luminance standard deviation >= minStd AND the most common colour covers <= maxShare of the pixels
// AND at least minColours distinct colours (a flat fallback plate, an empty chart or a failed render fails).
const fs = require("fs");
const zlib = require("zlib");

function decode(buf) {
  if (buf.length < 33 || buf.readUInt32BE(0) !== 0x89504e47) throw new Error("not a PNG");
  let pos = 8, w = 0, h = 0, depth = 8, ctype = 2, interlace = 0, plte = null;
  const idat = [];
  while (pos < buf.length) {
    const len = buf.readUInt32BE(pos), type = buf.toString("ascii", pos + 4, pos + 8), data = buf.subarray(pos + 8, pos + 8 + len);
    if (type === "IHDR") { w = data.readUInt32BE(0); h = data.readUInt32BE(4); depth = data[8]; ctype = data[9]; interlace = data[12]; }
    else if (type === "PLTE") plte = data;
    else if (type === "IDAT") idat.push(data);
    else if (type === "IEND") break;
    pos += 12 + len;
  }
  if (depth !== 8 || interlace !== 0) throw new Error(`unsupported PNG (depth ${depth}, interlace ${interlace})`);
  const ch = { 0: 1, 2: 3, 3: 1, 4: 2, 6: 4 }[ctype];
  if (!ch) throw new Error(`unsupported colour type ${ctype}`);
  const raw = zlib.inflateSync(Buffer.concat(idat));
  const stride = w * ch, out = Buffer.alloc(h * stride);
  for (let y = 0; y < h; y++) {
    const f = raw[y * (stride + 1)], src = y * (stride + 1) + 1, dst = y * stride;
    for (let x = 0; x < stride; x++) {
      const a = x >= ch ? out[dst + x - ch] : 0, b = y ? out[dst - stride + x] : 0, c = x >= ch && y ? out[dst - stride + x - ch] : 0;
      let v = raw[src + x];
      if (f === 1) v += a; else if (f === 2) v += b; else if (f === 3) v += (a + b) >> 1;
      else if (f === 4) { const p = a + b - c, pa = Math.abs(p - a), pb = Math.abs(p - b), pc = Math.abs(p - c); v += pa <= pb && pa <= pc ? a : pb <= pc ? b : c; }
      out[dst + x] = v & 255;
    }
  }
  return { w, h, ch, ctype, px: out, plte };
}

function pngStats(file, { step = 3 } = {}) {
  const { w, h, ch, ctype, px, plte } = decode(fs.readFileSync(file));
  let n = 0, sum = 0, sum2 = 0;
  const counts = new Map();
  for (let y = 0; y < h; y += step) {
    for (let x = 0; x < w; x += step) {
      const i = (y * w + x) * ch;
      let r, g, b;
      if (ctype === 0 || ctype === 4) r = g = b = px[i];
      else if (ctype === 3) { r = plte[px[i] * 3]; g = plte[px[i] * 3 + 1]; b = plte[px[i] * 3 + 2]; }
      else { r = px[i]; g = px[i + 1]; b = px[i + 2]; }
      const l = 0.2126 * r + 0.7152 * g + 0.0722 * b;
      sum += l; sum2 += l * l; n++;
      const key = ((r >> 2) << 12) | ((g >> 2) << 6) | (b >> 2);   // 6 bits per channel: ignores anti-alias noise
      counts.set(key, (counts.get(key) || 0) + 1);
    }
  }
  const mean = sum / n, variance = Math.max(0, sum2 / n - mean * mean);
  let top = 0;
  for (const c of counts.values()) if (c > top) top = c;
  return { width: w, height: h, variance, std: Math.sqrt(variance), topColourShare: top / n, distinctColours: counts.size };
}

function assertPng(file, { minStd = 4, maxShare = 0.97, minColours = 12 } = {}) {
  if (!fs.existsSync(file)) return { ok: false, reason: "file missing" };
  let st;
  try { st = pngStats(file); } catch (e) { return { ok: false, reason: `unreadable PNG: ${e.message}` }; }
  const bad = [];
  if (st.std < minStd) bad.push(`luminance std ${st.std.toFixed(2)} < ${minStd}`);
  if (st.topColourShare > maxShare) bad.push(`one colour covers ${(100 * st.topColourShare).toFixed(1)} % > ${100 * maxShare} %`);
  if (st.distinctColours < minColours) bad.push(`only ${st.distinctColours} distinct colours`);
  return { ok: bad.length === 0, reason: bad.join("; "), ...st };
}

module.exports = { pngStats, assertPng };

if (require.main === module) {
  let bad = 0;
  for (const f of process.argv.slice(2)) {
    const r = assertPng(f);
    console.log(`${r.ok ? "PASS" : "FAIL"}  ${f}  ${r.ok ? `std ${r.std.toFixed(1)}, ${r.distinctColours} colours` : r.reason}`);
    if (!r.ok) bad++;
  }
  process.exit(bad ? 1 : 0);
}
