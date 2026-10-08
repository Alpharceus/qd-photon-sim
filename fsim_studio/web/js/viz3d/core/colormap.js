// Colour mapping for the 3D views: perceptual viridis for data, wavelength
// -> sRGB for emitted light (NIR rendered as a flagged false colour), and
// role tints for layers. No physics here: inputs are SceneSpec numbers.

// Viridis (matplotlib, CC0) sampled at 17 stops; linear interpolation
// between stops in sRGB stays perceptually monotone at this density.
const VIRIDIS = [
  [68, 1, 84], [72, 26, 108], [71, 47, 125], [65, 68, 135], [57, 86, 140],
  [49, 104, 142], [42, 120, 142], [35, 136, 142], [31, 152, 139], [34, 168, 132],
  [53, 183, 121], [84, 197, 104], [122, 209, 81], [165, 219, 54], [210, 226, 27],
  [236, 229, 27], [253, 231, 37],
];

export function viridis(t) {
  if (!Number.isFinite(t)) return [0.35, 0.36, 0.38];
  const x = Math.min(1, Math.max(0, t)) * (VIRIDIS.length - 1);
  const i = Math.min(VIRIDIS.length - 2, Math.floor(x));
  const f = x - i;
  const a = VIRIDIS[i], b = VIRIDIS[i + 1];
  return [0, 1, 2].map((k) => (a[k] + (b[k] - a[k]) * f) / 255);
}

export function viridisCss(t) {
  const [r, g, b] = viridis(t);
  return `rgb(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)})`;
}

export function viridisGradientCss(dir = 'to right') {
  const stops = VIRIDIS.map((c, i) => `rgb(${c[0]},${c[1]},${c[2]}) ${(100 * i / (VIRIDIS.length - 1)).toFixed(1)}%`);
  return `linear-gradient(${dir}, ${stops.join(', ')})`;
}

// Wavelength (nm) -> linear-ish sRGB triple in 0..1 plus a false-colour flag.
// Visible range uses the Bruton piecewise approximation (display only);
// above 700 nm the light is invisible, so a fixed deep red is used and
// flagged so the UI can print "NIR, false colour".
export function lambdaToRGB(nm) {
  if (!Number.isFinite(nm)) return { rgb: [0.85, 0.85, 0.85], falseColour: true };
  if (nm > 700) return { rgb: [0.78, 0.06, 0.09], falseColour: true };
  if (nm < 380) return { rgb: [0.45, 0.2, 0.75], falseColour: true };
  let r = 0, g = 0, b = 0;
  if (nm < 440) { r = -(nm - 440) / 60; b = 1; }
  else if (nm < 490) { g = (nm - 440) / 50; b = 1; }
  else if (nm < 510) { g = 1; b = -(nm - 510) / 20; }
  else if (nm < 580) { r = (nm - 510) / 70; g = 1; }
  else if (nm < 645) { r = 1; g = -(nm - 645) / 65; }
  else { r = 1; }
  let f = 1;
  if (nm > 680) f = 0.3 + 0.7 * (700 - nm) / 20;
  else if (nm < 420) f = 0.3 + 0.7 * (nm - 380) / 40;
  return { rgb: [r * f, g * f, b * f], falseColour: false };
}

export function rgbCss([r, g, b], a = 1) {
  return `rgba(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)},${a})`;
}

// Layer role -> base tint (machined, desaturated; the beam ink is never used).
export const ROLE_TINT = {
  substrate: 0x1c1f25,
  cladding: 0x8794a6,
  core: 0x34587f,
  active: 0xd9a845,
  mirror: 0xa9b6c6,
  contact: 0x9a7a52,
  oxide: 0x56626e,
  doped_p: 0x7d4c6d,
  doped_n: 0x3c6d66,
  intrinsic: 0x666d77,
  thermal: 0x6b7480,
};

// Fixed enclosure palette (the viewport is a dark window in both themes).
export const ENCLOSURE = {
  bg: 0x0b0d10,
  floor: 0x111419,
  grid: 0x1b1f26,
  ink: '#e9e8e3',
  ink2: '#b8b7ae',
  muted: '#7f7d77',
};
