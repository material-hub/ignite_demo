/* ignite similarity-only planner, in the browser (his word 21 Oct 2026, CdSe0.50Te0.50: "can not it use similarity
 * model to get synthesis route, impurity etc?").
 *
 * A port of the parts of plan.context() that need no Materials Project, NIST-JANAF, kinetics, ReactCA or phonon data:
 * the learned precursor ranker (learned.py) over the candidate generators of planner.py, the learned firing
 * temperature (temperature.py) with its conformal range, prediction set and abstention (uncertainty.py), the route
 * share, batch, off-gas, atmosphere model, crucible, heating schedule and warnings (practical.py), the likely
 * impurities (impurity.py), the most similar known syntheses, the robot-lab outcomes and the measured pairwise
 * reactions. The knowledge is ignite_lite.json(.gz), written by synth/agentic/lite.py from the files plan.py reads.
 * The text matches synth.agentic.lite.reference() line for line (tests/lite_parity.js); numbers are formatted as
 * Python formats them (half-even on the exact binary value).
 *
 *   IGNITE_LITE.load(url) -> Promise   fetch the knowledge once (gzip, DecompressionStream)
 *   IGNITE_LITE.plan(formula)          {target, ranked_sets, fire_C, temperature, sections}, or {error}
 *   IGNITE_LITE.stats()                {targets, recipes} of the stored solid-state knowledge
 */
(function (root) {
  'use strict';
  const ENV = new Set(['O', 'H', 'C', 'N']);
  const ANIONS = new Set(['F', 'Cl', 'Br', 'I', 'S', 'Se', 'Te']);
  const GASES = ['CO2', 'H2O', 'N2', 'O2'];
  const GAS_OF = [['C', 'CO2'], ['H', 'H2O'], ['N', 'N2'], ['O', 'O2']];  // scorecard.GASES, in its order
  const VAPOURS = ['S', 'Se', 'Te'];
  const PER_ELEMENT = 6, BEAM = 300;
  const CLASSES = ['air', 'O2', 'inert', 'reducing'];
  const SEVERE = new Set(['H300', 'H301', 'H310', 'H311', 'H330', 'H331', 'H340', 'H341', 'H350', 'H351', 'H360', 'H361', 'H372']);
  let B = null, K = null;

  // ---------- Python number formatting and sums ----------
  const DV = new DataView(new ArrayBuffer(8));
  function bits(a) {  // a >= 0 finite: [m, e] with a = m * 2^e exactly
    DV.setFloat64(0, a);
    const hi = DV.getUint32(0), lo = DV.getUint32(4);
    const E = (hi >>> 20) & 0x7ff;
    let m = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo);
    if (E === 0) return [m, -1074];
    return [m | (1n << 52n), E - 1075];
  }
  function roundScaled(a, k) {  // round_half_even(a * 10^k) as BigInt, exactly
    const [m, e] = bits(a);
    let num = m, den = 1n;
    if (e >= 0) num <<= BigInt(e); else den <<= BigInt(-e);
    if (k >= 0) num *= 10n ** BigInt(k); else den *= 10n ** BigInt(-k);
    let q = num / den;
    const r2 = (num % den) * 2n;
    if (r2 > den || (r2 === den && (q & 1n) === 1n)) q += 1n;
    return q;
  }
  function special(x) { return Number.isNaN(x) ? 'nan' : x > 0 ? 'inf' : '-inf'; }
  function fixed(x, d) {  // f"{x:.{d}f}"
    if (!Number.isFinite(x)) return special(x);
    const neg = x < 0 || Object.is(x, -0);
    let s = roundScaled(Math.abs(x), d).toString();
    if (d > 0) {
      s = s.padStart(d + 1, '0');
      s = s.slice(0, s.length - d) + '.' + s.slice(s.length - d);
    }
    return (neg ? '-' : '') + s;
  }
  const plus = (x, d) => (x < 0 || Object.is(x, -0) ? '' : '+') + fixed(x, d);  // f"{x:+.{d}f}"
  const pct = (x, d = 0) => fixed(x * 100, d) + '%';  // f"{x:.{d}%}"
  function gfmt(x, p = 6) {  // f"{x:.{p}g}"
    if (!Number.isFinite(x)) return special(x);
    if (p === 0) p = 1;
    if (x === 0) return (Object.is(x, -0) ? '-' : '') + '0';
    const neg = x < 0, a = Math.abs(x);
    let E = Math.floor(Math.log10(a)), D;
    const top = 10n ** BigInt(p), bot = 10n ** BigInt(p - 1);
    for (;;) {
      D = roundScaled(a, p - 1 - E);
      if (D >= top) { E += 1; continue; }
      if (D < bot) { E -= 1; continue; }
      break;
    }
    let ds = D.toString(), s;
    if (E >= -4 && E < p) {
      const dec = p - 1 - E;
      if (dec > 0) {
        ds = ds.padStart(dec + 1, '0');
        s = (ds.slice(0, ds.length - dec) + '.' + ds.slice(ds.length - dec)).replace(/0+$/, '').replace(/\.$/, '');
      } else s = ds;
    } else {
      let mant = ds.length > 1 ? (ds[0] + '.' + ds.slice(1)).replace(/0+$/, '').replace(/\.$/, '') : ds;
      s = mant + 'e' + (E < 0 ? '-' : '+') + String(Math.abs(E)).padStart(2, '0');
    }
    return (neg ? '-' : '') + s;
  }
  const pyround = (x, n) => Number(fixed(x, n));
  function pysum(a) {  // Python 3.12 sum() of floats (Neumaier)
    let f = 0, c = 0;
    for (const x of a) {
      const t = f + x;
      if (Math.abs(f) >= Math.abs(x)) c += (f - t) + x; else c += (x - t) + f;
      f = t;
    }
    if (c !== 0 && Number.isFinite(c)) f += c;
    return f;
  }
  function pw(a, lo, n) {  // numpy pairwise_sum
    if (n < 8) { let r = 0; for (let i = 0; i < n; i++) r += a[lo + i]; return r; }
    if (n <= 128) {
      const r = a.slice(lo, lo + 8);
      let i;
      for (i = 8; i < n - (n % 8); i += 8) for (let j = 0; j < 8; j++) r[j] += a[lo + i + j];
      let res = ((r[0] + r[1]) + (r[2] + r[3])) + ((r[4] + r[5]) + (r[6] + r[7]));
      for (; i < n; i++) res += a[lo + i];
      return res;
    }
    let n2 = Math.floor(n / 2);
    n2 -= n2 % 8;
    return pw(a, lo, n2) + pw(a, lo + n2, n - n2);
  }
  const npsum = a => 0 + pw(a, 0, a.length);
  const npmean = a => npsum(a) / a.length;
  const npavg = (a, w) => npsum(a.map((x, i) => x * w[i])) / npsum(w);
  function npmedian(a) {
    const s = [...a].sort((x, y) => x - y), n = s.length, h = n >> 1;
    return n % 2 ? s[h] : (s[h - 1] + s[h]) / 2;
  }
  function searchsorted(arr, v) { let i = 0; while (i < arr.length && arr[i] < v) i++; return i; }
  function cumsum(a) { const o = []; let s = 0; for (const x of a) { s += x; o.push(s); } return o; }
  function wmedian(v, w) {
    const o = v.map((_, i) => i).sort((i, j) => v[i] - v[j]);
    const vs = o.map(i => v[i]), ws = o.map(i => w[i]);
    return vs[searchsorted(cumsum(ws), npsum(ws) / 2)];
  }
  const sortDesc = (arr, f) => arr.map((x, i) => [x, i]).sort((a, b) => f(b[0]) - f(a[0]) || a[1] - b[1]).map(x => x[0]);

  // ---------- formulas (pymatgen Composition) ----------
  const GROUP = /([A-Z][a-z]*)\s*([-*.e\d]*)/g, PAREN = /\(([^()]+)\)\s*([.e\d]*)/;
  function pyfloat(s) {
    if (!/^[+-]?(\d+\.?\d*|\.\d+)(e[+-]?\d+)?$/.test(s)) throw new Error(`could not convert string to float: '${s}'`);
    return Number(s);
  }
  function symDict(form, factor) {
    const d = new Map();
    for (const m of [...form.matchAll(GROUP)]) {
      const amt = m[2].trim() !== '' ? pyfloat(m[2]) : 1.0;
      d.set(m[1], (d.get(m[1]) || 0) + amt * factor);
      form = form.replace(m[0], '');
    }
    if (form.trim()) throw new Error(`${form} is an invalid formula!`);
    return d;
  }
  const PARSE = new Map();
  function comp(formula) {  // Composition(formula).get_el_amt_dict(), as a Map in pymatgen's order
    if (PARSE.has(formula)) return PARSE.get(formula);
    if (/^[\s\d.*/]*$/.test(formula)) throw new Error(`Invalid formula=${formula}`);
    let f = formula.replace(/@/g, '').replace(/[[{]/g, '(').replace(/[\]}]/g, ')');
    let m;
    while ((m = PAREN.exec(f))) {
      const factor = m[2] !== '' ? pyfloat(m[2]) : 1.0;
      const u = symDict(m[1], factor);
      f = f.replace(m[0], [...u].map(([e, a]) => e + pyrepr(a)).join(''));
    }
    const out = new Map();
    for (const [e, a] of symDict(f, 1)) {
      if (!Object.hasOwn(B.elements, e)) throw new Error(`${e} is not a valid element`);
      if (a < -1e-8) throw new Error('Amounts in Composition cannot be negative!');
      if (Math.abs(a) >= 1e-8) out.set(e, a);
    }
    PARSE.set(formula, out);
    return out;
  }
  function pyrepr(x) {  // repr(float), for the expanded parenthesis text that is parsed again
    if (Number.isInteger(x) && Math.abs(x) < 1e16) return x.toFixed(1);
    let s = String(x);
    if (/e/.test(s)) { const [a, b] = s.split('e'); s = a + 'e' + (b[0] === '-' ? '-' : '+') + b.replace(/^[+-]/, '').padStart(2, '0'); }
    return s;
  }
  const amt = (f, e) => comp(f).get(e) || 0;
  const nAtoms = f => { let s = 0; for (const a of comp(f).values()) s += Math.abs(a); return s; };
  function weight(f) { let s = 0; for (const [e, a] of comp(f)) s += a * B.elements[e][2]; return s; }
  function fracVec(f) { const n = nAtoms(f), o = new Map(); for (const [e, a] of comp(f)) o.set(e, a / n); return o; }
  const MAIN = new Map();
  function main(f) {
    if (!MAIN.has(f)) MAIN.set(f, [...comp(f).keys()].filter(e => !ENV.has(e)).sort());
    return MAIN.get(f);
  }
  const subset = (a, b) => { const s = b instanceof Set ? b : new Set(b); return a.every(x => s.has(x)); };
  const show = ps => ps.map(p => B.names[p] ?? p).join(' + ');
  function cosine(a, b) {
    const num = pysum([...a].map(([k, v]) => v * (b.get(k) || 0)));
    return num / (Math.sqrt(pysum([...a.values()].map(v => v * v))) * Math.sqrt(pysum([...b.values()].map(v => v * v))));
  }
  function canon(f) {  // the proportionality class of an integer composition (reduced_formula of it), else null
    const c = comp(f), v = [...c.values()];
    if (!v.every(x => Math.abs(x - Math.round(x)) < 1e-8)) return null;
    const ints = [...c].map(([e, a]) => [e, Math.round(a)]);
    const gcd = (x, y) => (y ? gcd(y, x % y) : x);
    const g = ints.reduce((s, [, a]) => gcd(s, a), 0) || 1;
    return ints.map(([e, a]) => [e, a / g]).sort((x, y) => (x[0] < y[0] ? -1 : 1)).map(([e, a]) => e + a).join(' ');
  }
  function anionClass(p) {
    const els = [...comp(p).keys()].filter(e => ENV.has(e));
    if (!els.length) return 'element';
    const s = new Set(els), has = e => s.has(e);
    if (s.size === 1 && has('O')) return 'oxide';
    if (s.size === 2 && has('C') && has('O')) return 'carbonate';
    if (s.size === 2 && has('H') && has('O')) return 'hydroxide';
    if (has('N') && !has('C')) return amt(p, 'O') >= 3 * amt(p, 'N') ? 'nitrate' : 'ammonium';
    if (has('C')) return 'organic';
    return 'other';
  }

  // ---------- linear algebra ----------
  function lstsq(A, b) {  // minimum-norm least squares (numpy lstsq, rcond = eps * max(m, n)) by one-sided Jacobi SVD
    const m = A.length, n = A[0].length;
    const W = A.map(r => r.slice()), V = Array.from({length: n}, (_, i) => Array.from({length: n}, (_, j) => +(i === j)));
    for (let sweep = 0; sweep < 60; sweep++) {
      let off = 0;
      for (let p = 0; p < n - 1; p++) for (let q = p + 1; q < n; q++) {
        let al = 0, be = 0, ga = 0;
        for (let i = 0; i < m; i++) { al += W[i][p] ** 2; be += W[i][q] ** 2; ga += W[i][p] * W[i][q]; }
        if (Math.abs(ga) <= 1e-15 * Math.sqrt(al * be) || ga === 0) continue;
        off = Math.max(off, Math.abs(ga) / Math.sqrt(al * be));
        const z = (be - al) / (2 * ga), t = Math.sign(z || 1) / (Math.abs(z) + Math.sqrt(1 + z * z));
        const c = 1 / Math.sqrt(1 + t * t), s = c * t;
        for (let i = 0; i < m; i++) { const x = W[i][p], y = W[i][q]; W[i][p] = c * x - s * y; W[i][q] = s * x + c * y; }
        for (let i = 0; i < n; i++) { const x = V[i][p], y = V[i][q]; V[i][p] = c * x - s * y; V[i][q] = s * x + c * y; }
      }
      if (off < 1e-15) break;
    }
    const sig = Array.from({length: n}, (_, j) => Math.sqrt(W.reduce((s, r) => s + r[j] ** 2, 0)));
    const tol = Math.max(...sig) * 2.220446049250313e-16 * Math.max(m, n);
    const x = new Array(n).fill(0);
    for (let j = 0; j < n; j++) {
      if (!(sig[j] > tol)) continue;
      let wb = 0;
      for (let i = 0; i < m; i++) wb += W[i][j] * b[i];
      const f = wb / (sig[j] * sig[j]);
      for (let i = 0; i < n; i++) x[i] += V[i][j] * f;
    }
    return x;
  }
  const matvec = (A, x) => A.map(r => r.reduce((s, v, j) => s + v * x[j], 0));
  function solve(M, y) {  // Gaussian elimination with partial pivoting
    const n = y.length, A = M.map((r, i) => [...r, y[i]]);
    for (let c = 0; c < n; c++) {
      let p = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
      [A[c], A[p]] = [A[p], A[c]];
      if (A[c][c] === 0) throw new Error('singular matrix');
      for (let r = c + 1; r < n; r++) { const f = A[r][c] / A[c][c]; for (let k = c; k <= n; k++) A[r][k] -= f * A[c][k]; }
    }
    const x = new Array(n).fill(0);
    for (let r = n - 1; r >= 0; r--) { let s = A[r][n]; for (let k = r + 1; k < n; k++) s -= A[r][k] * x[k]; x[r] = s / A[r][r]; }
    return x;
  }
  function nnls(A, b) {  // Lawson-Hanson active set; returns [x, ||Ax - b||]
    const m = A.length, n = A[0].length, tol = 10 * Math.max(m, n) * 2.220446049250313e-16;
    const AtA = Array.from({length: n}, (_, i) => Array.from({length: n}, (_, j) => A.reduce((s, r) => s + r[i] * r[j], 0)));
    const Atb = Array.from({length: n}, (_, j) => A.reduce((s, r, i) => s + r[j] * b[i], 0));
    let x = new Array(n).fill(0), P = new Array(n).fill(false), w = Atb.slice(), it = 0;
    const sub = () => {
      const idx = P.map((v, i) => (v ? i : -1)).filter(i => i >= 0);
      const z = solve(idx.map(i => idx.map(j => AtA[i][j])), idx.map(i => Atb[i]));
      const s = new Array(n).fill(0);
      idx.forEach((i, k) => { s[i] = z[k]; });
      return s;
    };
    while (!P.every(Boolean) && w.some((v, i) => !P[i] && v > tol)) {
      let k = -1;
      for (let i = 0; i < n; i++) if (!P[i] && (k < 0 || w[i] > w[k])) k = i;
      P[k] = true;
      let s = sub();
      while (it < 3 * n && P.some((v, i) => v && s[i] < 0)) {
        it++;
        let alpha = Infinity;
        for (let i = 0; i < n; i++) if (P[i] && s[i] < 0) alpha = Math.min(alpha, x[i] / (x[i] - s[i]));
        x = x.map((v, i) => v * (1 - alpha) + alpha * s[i]);
        P = P.map((v, i) => v && x[i] > tol);
        s = sub();
      }
      x = s;
      const Ax = Array.from({length: n}, (_, i) => AtA[i].reduce((t, v, j) => t + v * x[j], 0));
      w = Atb.map((v, i) => v - Ax[i]);
    }
    const r = matvec(A, x).map((v, i) => v - b[i]);
    return [x, Math.sqrt(r.reduce((s, v) => s + v * v, 0))];
  }
  function feasible(A, b, lb, ub) {  // does A x = b, lb <= x <= ub have a solution (phase-1 simplex, Bland's rule)
    const m = A.length, cols = [], shift = new Array(A[0].length).fill(0);
    A[0].forEach((_, j) => {
      const lo = lb[j], hi = ub[j];
      if (lo !== null) { shift[j] = lo; cols.push([j, 1]); if (hi !== null) throw new Error('two-sided bound'); }
      else if (hi !== null) { shift[j] = hi; cols.push([j, -1]); }
      else { cols.push([j, 1]); cols.push([j, -1]); }
    });
    const nv = cols.length;
    const T = A.map((r, i) => {
      let rhs = b[i] - r.reduce((s, v, j) => s + v * shift[j], 0);
      let row = cols.map(([j, sg]) => r[j] * sg);
      if (rhs < 0) { row = row.map(v => -v); rhs = -rhs; }
      return [...row, ...Array.from({length: m}, (_, k) => +(k === i)), rhs];
    });
    const W = nv + m, basis = Array.from({length: m}, (_, i) => nv + i);
    const obj = new Array(W + 1).fill(0);  // reduced costs of the phase-1 objective (sum of artificials)
    for (let i = 0; i < m; i++) for (let j = 0; j <= W; j++) if (j < nv || j === W) obj[j] -= T[i][j];
    for (let it = 0; it < 5000; it++) {
      let e = -1;
      for (let j = 0; j < W; j++) if (obj[j] < -1e-12) { e = j; break; }
      if (e < 0) break;
      let r = -1, best = Infinity;
      for (let i = 0; i < m; i++) if (T[i][e] > 1e-12) {
        const q = T[i][W] / T[i][e];
        if (q < best - 1e-15 || (Math.abs(q - best) <= 1e-15 && basis[i] < basis[r])) { best = q; r = i; }
      }
      if (r < 0) break;
      const pv = T[r][e];
      for (let j = 0; j <= W; j++) T[r][j] /= pv;
      for (let i = 0; i < m; i++) if (i !== r && T[i][e] !== 0) { const f = T[i][e]; for (let j = 0; j <= W; j++) T[i][j] -= f * T[r][j]; }
      const f = obj[e];
      for (let j = 0; j <= W; j++) obj[j] -= f * T[r][j];
      basis[r] = e;
    }
    return -obj[W] <= 1e-7;
  }

  // ---------- scorecard.balance(target, set, air=True) is not None ----------
  function balances(target, precursors) {
    const t = comp(target), m = precursors.length;
    let els = new Set(t.keys());
    for (const p of precursors) for (const e of comp(p).keys()) els.add(e);
    let species = [...precursors, ...GAS_OF.filter(([e]) => els.has(e)).map(([, g]) => g)];
    for (const g of species.slice(m)) for (const e of comp(g).keys()) els.add(e);
    species = [...precursors, ...GAS_OF.filter(([e]) => els.has(e)).map(([, g]) => g)];
    const E = [...els].sort();
    const A = E.map(e => species.map(s => amt(s, e))), b = E.map(e => t.get(e) || 0);
    const x = lstsq(A, b), Ax = matvec(A, x);
    const bad = Ax.some((v, i) => Math.abs(v - b[i]) > 1e-6) || x.slice(0, m).some(v => v <= 1e-6) ||
      species.some((s, i) => i >= m && s !== 'O2' && x[i] > 1e-6);
    if (!bad) return true;
    const air = new Set(GAS_OF.filter(([e, g]) => t.has(e) && g !== 'O2').map(([, g]) => g));
    const tryAir = sp => {
      const up = sp.map((s, i) => (i >= m && air.has(s) ? i : -1)).filter(i => i >= 0);
      const A2 = E.map(e => [...sp.map(s => amt(s, e)), ...up.map(i => amt(sp[i], e))]);
      const lb = [...sp.map((s, i) => (i < m ? 1e-6 : null)), ...up.map(() => 0)];
      const ub = [...sp.map((s, i) => (i < m || s === 'O2' ? null : 0)), ...up.map(() => null)];
      return feasible(A2, b, lb, ub);
    };
    return tryAir(species) || tryAir([...species, ...VAPOURS.filter(v => t.has(v) && !precursors.includes(v))]);
  }

  // ---------- planner.Planner ----------
  class Planner {
    constructor(power, k) { this.power = power; this.k = k; this.nb = new Map(); }
    neighbors(target) {
      if (this.nb.has(target)) return this.nb.get(target);
      const v = fracVec(target), need = new Set(main(target)), scored = [];
      K.kb.forEach((t, i) => { if (K.main[i].some(e => need.has(e))) scored.push([cosine(v, K.vec[i]), t]); });
      const out = sortDesc(scored, x => x[0]).slice(0, this.k);
      this.nb.set(target, out);
      return out;
    }
    weighted(target) { return this.neighbors(target).map(([s, t]) => [s ** this.power, t]); }
    globalCounts(el, anions) {
      const c = new Map(anions.size && anions.has(el) ? [] : K.byEl.get(el) || []);
      if (anions.size) for (const [p, n] of K.multiByEl.get(el) || []) {
        const r = main(p).filter(e => !anions.has(e));
        if (r.length === 1 && r[0] === el) c.set(p, (c.get(p) || 0) + n);
      }
      return c;
    }
    elementPrior(target, anions = new Set()) {
      const allowed = new Set(main(target)), need = main(target).filter(e => !anions.has(e)), needS = new Set(need);
      const local = new Map(need.map(e => [e, new Map()]));
      for (const [sim, t] of this.weighted(target)) for (const r of t.recipes) for (const p of r.precursors) {
        const mp = main(p), key = mp.filter(e => !anions.has(e));
        if (key.length === 1 && needS.has(key[0]) && subset(mp, allowed)) {
          const L = local.get(key[0]);
          L.set(p, (L.get(p) || 0) + sim);
        }
      }
      const prior = new Map();
      for (const el of need) {
        const g = this.globalCounts(el, anions), L = local.get(el);
        const gt = pysum([...g.values()]) || 1, lt = pysum([...L.values()]) || 1;
        const cand = [...new Set([...g.keys(), ...L.keys()])];
        const sc = cand.map(p => [p, Math.log(0.7 * (L.get(p) || 0) / lt + 0.3 * (g.get(p) || 0) / gt + 1e-6)]);
        prior.set(el, sortDesc(sc, x => x[1]).slice(0, PER_ELEMENT));
      }
      return prior;
    }
    candidates(target) {
      const prior = [...this.elementPrior(target).values()];
      if (prior.some(v => !v.length)) return [];
      const combos = [];
      product(prior, choice => combos.push([[...new Set(choice.map(c => c[0]))].sort(), pysum(choice.map(c => c[1]))]));
      return dedupe(sortDesc(combos, x => x[1]), BEAM);
    }
    candidatesMulti(target) {
      const need = main(target), needS = new Set(need), prior = this.elementPrior(target);
      const local = new Map();
      let total = 0;
      for (const [sim, t] of this.weighted(target)) for (const r of t.recipes) {
        total += sim;
        for (const p of new Set(r.precursors)) {
          const mp = main(p);
          if (mp.length > 1 && subset(mp, needS)) local.set(p, (local.get(p) || 0) + sim);
        }
      }
      const n = Math.max(K.nRecipes, 1);
      const multi = sortDesc([...local].map(([p, w]) => [p, Math.log(0.7 * w / (total || 1) + 0.3 * (K.multi.get(p) || 0) / n + 1e-6)]),
        x => x[1]).slice(0, 5);
      const combos = [...prior.values()].every(v => v.length) ? this.candidates(target).slice() : [];
      for (const [mp, ms] of multi) {
        const ms_ = new Set(main(mp)), rest = need.filter(e => !ms_.has(e)).map(e => prior.get(e));
        if (rest.some(v => !v.length)) continue;
        product(rest, choice => combos.push([[...new Set([mp, ...choice.map(c => c[0])])].sort(), ms + pysum(choice.map(c => c[1]))]));
      }
      return dedupe(sortDesc(combos, x => x[1]), BEAM);
    }
    candidatesAnion(target) {
      const allowed = new Set(main(target)), anions = new Set(main(target).filter(e => ANIONS.has(e)));
      if (!anions.size || allowed.size === anions.size) return [];
      const prior = [...this.elementPrior(target, anions).values()];
      if (prior.some(v => !v.length)) return [];
      const sources = new Map();
      for (const [sim, t] of this.weighted(target)) for (const r of t.recipes) for (const p of r.precursors) {
        const mp = main(p);
        if (subset(mp, allowed)) for (const a of mp.filter(e => anions.has(e))) {
          if (!sources.has(a)) sources.set(a, new Map());
          const S = sources.get(a);
          S.set(p, (S.get(p) || 0) + sim);
        }
      }
      const out = new Map();
      product(prior, choice => {
        const s = new Set(choice.map(c => c[0]));
        let sc = pysum(choice.map(c => c[1]));
        const have = new Set([...s].flatMap(p => main(p)));
        for (const a of [...anions].filter(a => !have.has(a)).sort()) {
          const S = sources.get(a);
          if (!S || !S.size) return;
          let bp = null, bw = -Infinity;
          for (const [p, w] of S) if (w > bw) { bp = p; bw = w; }
          s.add(bp);
          sc += Math.log(bw / pysum([...S.values()]));
        }
        const key = [...s].sort(), kk = key.join('\u0000');
        out.set(kk, [key, Math.max(out.has(kk) ? out.get(kk)[1] : -1e9, sc)]);
      });
      return sortDesc([...out.values()], x => x[1]).slice(0, BEAM);
    }
    drivingForce(target, ps) {
      if (ps.some(p => !Object.hasOwn(B.E, p))) return null;
      const t = comp(target);
      let els = new Set(t.keys());
      for (const p of ps) for (const e of comp(p).keys()) els.add(e);
      const species = [...ps, ...GASES.filter(g => [...comp(g).keys()].some(e => els.has(e)))];
      for (const s of species) for (const e of comp(s).keys()) els.add(e);
      const E = [...els].sort(), A = E.map(e => species.map(s => amt(s, e))), b = E.map(e => t.get(e) || 0);
      const x = lstsq(A, b), Ax = matvec(A, x);
      if (Ax.some((v, i) => Math.abs(v - b[i]) > 1e-6) || x.slice(0, ps.length).some(v => v <= 1e-6)) return null;
      let de = 0;
      species.forEach((s, i) => { de += x[i] * B.E[s]; });
      return de / pysum([...t.values()]);
    }
    temperature(target) {
      const [ts, ws] = nbTemps(this, target);
      if (!ts.length) return npmedian(K.temps);
      const o = ts.map((_, i) => i).sort((i, j) => ts[i] - ts[j]);
      const cw = cumsum(o.map(i => ws[i]));
      return ts[o[searchsorted(cw, cw[cw.length - 1] / 2)]];
    }
  }
  function nbTemps(pl, target) {
    const ts = [], ws = [];
    for (const [sim, t] of pl.neighbors(target).slice(0, 10)) for (const r of t.recipes) if (r.temp_C) { ts.push(r.temp_C); ws.push(sim); }
    return [ts, ws];
  }
  function product(lists, fn) {
    const cur = [];
    const go = i => {
      if (i === lists.length) { fn(cur.slice()); return; }
      for (const x of lists[i]) { cur.push(x); go(i + 1); cur.pop(); }
    };
    go(0);
  }
  function dedupe(sorted, beam) {
    const seen = new Set(), out = [];
    for (const [s, sc] of sorted) {
      const k = s.join('\u0000');
      if (!seen.has(k)) { seen.add(k); out.push([s, sc]); }
      if (out.length >= beam) break;
    }
    return out;
  }

  // ---------- learned.Featurizer and plan.ranked_sets / plan_sets ----------
  function rows(pl, target) {
    const base = pl.candidatesMulti(target), seen = new Set(base.map(([s]) => s.join('\u0000')));
    let c = [...base.map(([s, sc]) => [s, sc, 0]),
      ...pl.candidatesAnion(target).filter(([s]) => !seen.has(s.join('\u0000'))).map(([s, sc]) => [s, sc, 1])];
    c = sortDesc(c.filter(x => x[0].length), x => x[1]);
    if (!c.length) return [[], []];
    const support = new Map(), cls = new Map();
    let total = 0, ncls = 0;
    for (const [sim, t] of pl.weighted(target)) for (const r of t.recipes) {
      total += sim;
      for (const p of new Set(r.precursors)) {
        support.set(p, (support.get(p) || 0) + sim);
        const a = anionClass(p);
        cls.set(a, (cls.get(a) || 0) + sim);
        ncls += sim;
      }
    }
    const nEl = main(target).length, N = K.nRecipes;
    const dfs = c.map(([s]) => pl.drivingForce(target, s)), known = dfs.filter(d => d !== null);
    const best = known.length ? Math.max(...known) : 0, fill = known.length ? npmedian(known) : 0;
    const GASC = new Set(['carbonate', 'nitrate', 'hydroxide', 'organic', 'ammonium']);
    const X = c.map(([s, sc, an], rank) => {
      const sup = s.map(p => Math.log((support.get(p) || 0) / (total || 1) + 1e-3));
      const pmi = [];
      for (let i = 0; i < s.length; i++) for (let j = i + 1; j < s.length; j++) {
        const nab = K.pair.get(s[i] + '\u0000' + s[j]) || 0, na = K.single.get(s[i]) || 0, nb = K.single.get(s[j]) || 0;
        pmi.push(Math.log((nab + 0.1) * N / ((na + 1) * (nb + 1))));
      }
      const classes = s.map(anionClass), hull = s.filter(p => Object.hasOwn(B.H, p)).map(p => B.H[p]), df = dfs[rank];
      return [sc, an, s.filter(p => main(p).length > 1).length, s.length - nEl, Math.min(...sup), npmean(sup),
        pmi.length ? npmean(pmi) : 0, pmi.length ? Math.min(...pmi) : 0,
        npmean(classes.map(k => Math.log((cls.get(k) || 0) / (ncls || 1) + 1e-3))),
        npmean(classes.map(k => +GASC.has(k))), npmean(classes.map(k => +(k === 'element'))),
        (df !== null ? df : fill) - best, +(df === null), hull.length ? Math.min(Math.max(...hull), 0.5) : 0,
        +(hull.length < s.length), Math.log1p(rank)];
    });
    return [c.map(x => x[0]), X];
  }
  function planSets(target, k) {
    const W = B.w_learned, [sets, X] = rows(K.P2, target);
    if (!sets.length) return [[], null];
    const z = X.map(r => r.reduce((s, v, j) => s + (v - W.mu[j]) / W.sd[j] * W.w[j], 0));
    const o = z.map((_, i) => i).sort((i, j) => z[j] - z[i] || i - j);
    const zmax = z[o[0]];
    let p = o.map(i => Math.exp(z[i] - zmax));
    const ps = npsum(p);
    p = p.map(v => v / ps);
    const keep = [];
    o.forEach((i, r) => { if (balances(target, sets[i])) keep.push(r); });
    let probs = null;
    if (keep.length) {
      const q = keep.map(r => p[r]), qs = npsum(q);
      probs = q.map(v => v / qs);
    }
    return [keep.slice(0, k).map(r => sets[o[r]]), probs];
  }

  // ---------- temperature.py, plan.firing_temperature, uncertainty ----------
  function firing(target) {
    const pl = K.P1, [ts, ws] = nbTemps(pl, target), med = pl.temperature(target);
    const els = main(target), tm = els.map(e => B.elements[e][1] ?? NaN), chi = els.map(e => B.elements[e][0] ?? NaN);
    const cands = pl.candidates(target), classes = cands.length ? cands[0][0].map(anionClass) : [];
    const frac = c => (classes.length ? npmean(classes.map(x => +(x === c))) : 0);
    const nanmean = a => { const v = a.filter(x => !Number.isNaN(x)); return v.length ? npsum(a.map(x => (Number.isNaN(x) ? 0 : x))) / v.length : NaN; };
    const nanmax = a => { const v = a.filter(x => !Number.isNaN(x)); return v.length ? Math.max(...v) : NaN; };
    const x = [med, ts.length ? npavg(ts, ws) : med, ts.length ? ts[0] : med, nanmean(tm), nanmax(tm), nanmean(chi),
      els.length, +els.some(e => ANIONS.has(e)), frac('carbonate'), frac('nitrate'), frac('element'), frac('organic')];
    const W = B.w_temperature;
    let f = 0;
    x.forEach((v, j) => { f += ((Number.isNaN(v) ? W.fill[j] : v) - W.mu[j]) / W.sd[j] * W.w[j]; });
    return f + W.w[x.length];
  }
  function spread(target) {
    const [ts, ws] = nbTemps(K.P1, target);
    if (ts.length < 2) return NaN;
    const m = npavg(ts, ws);
    return Math.sqrt(npavg(ts.map(t => (t - m) ** 2), ws));
  }

  // ---------- practical.py ----------
  function neighborConditions(target) {
    const atm = new Map(), box = new Map(), two = [], first = [], dwell = [];
    for (const [w, t] of K.P2.weighted(target)) for (const r of t.recipes) {
      const c = r.cond;
      if (!c) continue;
      const [a, ct, temps, dh] = c;
      if (a) atm.set(a, (atm.get(a) || 0) + w);
      if (ct) box.set(ct, (box.get(ct) || 0) + w);
      two.push([temps.length >= 2, w]);
      if (temps.length >= 2 && temps[0] < Math.max(...temps)) first.push([temps[0], w]);
      if (dh) dwell.push([dh, w]);
    }
    const norm = c => { const tot = pysum([...c.values()]) || 1; return sortDesc([...c], x => x[1]).map(([k, v]) => [k, v / tot]); };
    const out = {atmosphere: norm(atm), container: norm(box)};
    if (two.length) out.p_two_step = pysum(two.filter(x => x[0]).map(x => x[1])) / pysum(two.map(x => x[1]));
    if (first.length) out.first_step_C = wmedian(first.map(x => x[0]), first.map(x => x[1]));
    if (dwell.length) out.dwell_h = wmedian(dwell.map(x => x[0]), dwell.map(x => x[1]));
    return out;
  }
  function limitDenominator(v, maxd) {  // Fraction(v).limit_denominator(maxd).denominator, v a double
    const [m, e] = bits(Math.abs(v));
    let n = m, d = 1n;
    if (e >= 0) n <<= BigInt(e); else d <<= BigInt(-e);
    const g0 = (a, b) => (b ? g0(b, a % b) : a), g = g0(n, d) || 1n;
    n /= g; d /= g;
    const N = n, Dd = d, MD = BigInt(maxd);
    if (d <= MD) return Number(d);
    let p0 = 0n, q0 = 1n, p1 = 1n, q1 = 0n;
    for (;;) {
      const a = n / d, q2 = q0 + a * q1;
      if (q2 > MD) break;
      [p0, q0, p1, q1] = [p1, q1, p0 + a * p1, q2];
      [n, d] = [d, n - a * d];
    }
    const k = (MD - q0) / q1;
    return Number(2n * d * (q0 + k * q1) <= Dd ? q1 : q0 + k * q1);
  }
  const OXI = new Map();
  function oxi(formula) {  // practical.oxi: pymatgen oxi_state_guesses(max_sites=-20)[0], or null
    if (OXI.has(formula)) return OXI.get(formula);
    const r = (() => {
      const d = comp(formula);
      let den = 1;
      const gcd = (a, b) => (b ? gcd(b, a % b) : a);
      for (const v of d.values()) { const q = limitDenominator(v, 100); den = den * q / gcd(den, q); }
      const c = new Map();
      for (const [e, v] of d) { const n = pyroundInt(v * den); if (Math.abs(n) >= 1e-8) c.set(e, n); }
      let n = 0;
      for (const v of c.values()) n += v;
      if (n > 80) return null;
      if (c.size === 1) return new Map([[[...c.keys()][0], 0]]);
      let fac = [...c.values()].reduce(gcd, 0);
      const red = new Map([...c].map(([e, v]) => [e, v / fac]));
      if (isSpecial(red)) fac /= 2;
      const rc = new Map([...c].map(([e, v]) => [e, v / fac]));
      let na = 0;
      for (const v of rc.values()) na += v;
      if (na > 20) return null;
      if (![...rc.values()].every(v => v === Math.trunc(v))) return null;
      const els = [...rc.keys()], sums = [], best = [];
      for (const el of els) {
        const ox = B.elements[el][3].length ? B.elements[el][3] : B.elements[el][4], P = B.oxi_prob[el] || {};
        const S = [], sc = new Map();
        combos(ox, rc.get(el), cmb => {
          let s = 0;
          for (const o of cmb) s += o;
          if (!S.includes(s)) S.push(s);
          const score = pysum(cmb.map(o => P[String(o)] || 0));
          if (!sc.has(s) || score > (sc.get(s) || 0)) sc.set(s, score);
        });
        sums.push(S); best.push(sc);
      }
      const sols = [];
      product(sums, x => {
        if (x.reduce((a, b) => a + b, 0) === 0) {
          let score = 0;
          x.forEach((v, i) => { score += best[i].get(v); });
          sols.push([score, new Map(els.map((e, i) => [e, x[i] / rc.get(e)]))]);
        }
      });
      if (!sols.length) return null;
      return sortDesc(sols, s => s[0])[0][1];
    })();
    OXI.set(formula, r);
    return r;
  }
  function isSpecial(red) {
    const k = [...red].map(([e, v]) => e + ':' + v).sort().join(' ');
    return ['Li:1 O:1', 'Na:1 O:1', 'K:1 O:1', 'H:1 O:1', 'Cs:1 O:1', 'O:1 Rb:1', 'O:1', 'N:1', 'F:1', 'Cl:1', 'H:1'].includes(k);
  }
  function combos(ox, n, fn) {  // itertools.combinations_with_replacement(ox, n)
    const cur = [];
    const go = (start, left) => {
      if (!left) { fn(cur.slice()); return; }
      for (let i = start; i < ox.length; i++) { cur.push(ox[i]); go(i, left - 1); cur.pop(); }
    };
    go(0, n);
  }
  function pyroundInt(x) {  // round(x) for a float, half to even
    const f = Math.floor(x), d = x - f;
    if (d > 0.5) return f + 1;
    if (d < 0.5) return f;
    return f % 2 === 0 ? f : f + 1;
  }
  function atmosphere(target, nc) {
    const rate = (table, k, c) => {
      const row = table[k];
      if (!row) return 0.25;
      let tot = 0;
      for (const v of Object.values(row)) tot += v;
      return ((row[c] || 0) + 0.25) / (tot + 1.0);
    };
    const nb = new Map(nc.atmosphere), g = oxi(target) || new Map(), els = main(target);
    const x = CLASSES.map(c => nb.get(c) || 0);
    for (const c of CLASSES.slice(1)) x.push(Math.max(...els.map(e => rate(B.atm_el, e, c))));
    for (const c of CLASSES.slice(1)) x.push(Math.max(...els.map(e => {
      const k = g.has(e) ? `${e}|${pyroundInt(g.get(e))}` : null;
      return k && Object.hasOwn(B.atm_ox, k) ? rate(B.atm_ox, k, c) : rate(B.atm_el, e, c);
    })));
    const a = new Set(comp(target).keys()), any = l => l.some(e => a.has(e));
    x.push(+!g.size, +any(['S', 'Se', 'Te']), +any(['F', 'Cl', 'Br', 'I']), +(a.has('N') && !a.has('O')),
      +!any(['O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I', 'N']));
    const W = B.w_atmosphere, z = [...x.map((v, j) => (v - W.mu[j]) / W.sd[j]), 1.0];
    const s = CLASSES.map((_, k) => z.reduce((t, v, j) => t + v * W.W[j][k], 0)), sm = Math.max(...s);
    const e = s.map(v => Math.exp(v - sm)), es = npsum(e);
    return CLASSES.map((c, i) => [c, e[i] / es]);
  }
  function practicalLines(target, ps, fire, grams) {
    const lines = [];
    const els = [...new Set([...main(target), ...ps.flatMap(p => main(p))])].sort();
    const [x, res] = nnls(els.map(e => ps.map(p => amt(p, e))), els.map(e => amt(target, e)));
    const mol = grams / weight(target);
    const rws = ps.map((p, i) => [p, x[i], x[i] * mol * weight(p)]);
    const total = pysum(rws.map(r => r[2]));
    lines.push(`Batch for ${gfmt(grams)} g ${target} with set 1 (molar ratio, grams to weigh):`);
    for (const [p, c, g] of rws) lines.push(`  ${B.names[p] ?? p}: ${gfmt(c, 4)} mol per mol target, ${fixed(g, 3)} g`);
    if (total === 0) throw new Error('ZeroDivisionError: float division by zero');
    lines.push(`  total ${fixed(total, 3)} g; expected mass loss on firing ${fixed(total - grams, 3)} g ` +
      `(${pct((total - grams) / total)}) as gas or O2 exchange` + (res > 1e-6 ? `; WARNING element balance residual ${gfmt(res, 3)}` : ''));
    const gas = [];
    for (const [p, c] of rws) {
      const nco2 = amt(p, 'C') * c * mol;
      if (nco2 > 0) gas.push([`CO2 from ${p}`, nco2 * 44.009]);
      const nh2o = amt(p, 'H') / 2 * c * mol;
      if (nh2o > 0) gas.push([`H2O from ${p}`, nh2o * 18.015]);
      const nno2 = amt(p, 'N') * c * mol;
      if (nno2 > 0) gas.push([`NOx as NO2 from ${p}`, nno2 * 46.006]);
    }
    if (gas.length) {
      const gt = pysum(gas.map(g => g[1]));
      lines.push(`  off-gas about ${fixed(gt, 3)} g per ${gfmt(grams)} g ${target} (${fixed(gt / grams, 2)} kg per kg): ` +
        gas.map(([n, m]) => `${fixed(m, 3)} g ${n}`).join('; ') + `; atom economy ${pct(grams / total)}` +
        (gt > total - grams ? ' (includes oxygen taken up from the air)' : ''));
    } else {
      lines.push(`  off-gas: no carbonate/nitrate/hydroxide gas expected from set 1 (mass change is O2 exchange); atom economy ${pct(grams / total)}`);
    }
    const nc = neighborConditions(target);
    const p = sortDesc(atmosphere(target, nc), v => v[1]).filter(v => v[1] >= 0.05);
    const M = B.metrics.atmosphere;
    lines.push(`Atmosphere (learned; held-out accuracy ${pct(M[0])} vs ${pct(M[1])} for always-air): ` + p.map(([k, v]) => `${k} ${pct(v)}`).join(', '));
    if (nc.container.length) lines.push('Crucible in the most similar syntheses: ' + nc.container.slice(0, 3).map(([k, v]) => `${k} ${pct(v)}`).join(', '));
    if ('p_two_step' in nc) {
      let s = `Heating schedule: ${pct(nc.p_two_step)} of similar syntheses heat in two or more steps`;
      if (nc.p_two_step >= 0.5 && 'first_step_C' in nc) s += `; first (calcination) step near ${fixed(Math.min(nc.first_step_C, fire), 0)} C, regrind, then ${fixed(fire, 0)} C`;
      if ('dwell_h' in nc) {
        const f = Math.exp(B.uncertainty.dwell.q_log['0.9']);
        s += `; typical dwell at the top step ${fixed(nc.dwell_h, 0)} h (90% range ${fixed(nc.dwell_h / f, 1)}-${fixed(nc.dwell_h * f, 0)} h)`;
      }
      lines.push(s);
    }
    const warn = [];
    if (ps.some(q => ['Li', 'Na', 'K'].some(e => comp(q).has(e))) && fire >= 800)
      warn.push('alkali in set 1 (Li/Na/K) may attack an alumina crucible above ~800 C (alkali aluminates, Al contamination as in ' +
        'LLZO practice with LaAlO3); use Pt/ZrO2 or a mother-powder bed');
    const BASIC = new Set(['SrO', 'BaO', 'CaO', 'SrCO3', 'BaCO3', 'CaCO3', 'Sr(NO3)2', 'Ba(NO3)2']);
    const basic = [...new Set(ps.filter(q => BASIC.has(B.red[q] ?? q)))].sort();
    if (basic.length && fire >= 1300) warn.push('basic oxide from ' + basic.join('+') + ' may react with alumina at this temperature (SrAl2O4/BaAl2O4); prefer Pt/MgO');
    for (const q of ps) {
      const h = B.hazards[q];
      if (h && h[1].length) warn.push(`${q}: ${h[1].join('; ')} (PubChem CID ${h[0]}, GHS)`);
    }
    const have = new Set(ps), tc = canon(target);
    for (const q of B.pairs) {
      if (q[3].every(r => have.has(r)) && q[4] !== null && !q[4].some(f => canon(f) === tc && tc !== null))
        warn.push(`${q[0].join(' + ')} forms ${q[1]} at ${q[2]} C (A-Lab measured); an early intermediate can use up the driving force`);
    }
    if (warn.length) { lines.push('Warnings:'); for (const w of warn) lines.push('  ' + w); }
    return lines;
  }

  // ---------- impurity.py ----------
  class Space {
    constructor(targets) {
      this.targets = targets;
      const els = new Set(ENV);
      for (const t of targets) for (const e of comp(t).keys()) els.add(e);
      this.els = [...els].sort();
      this.ix = new Map(this.els.map((e, i) => [e, i]));
      this.V = targets.map(t => this.vec(t));
      this.M = targets.map(t => this.main(t));
    }
    vec(t) {
      const v = new Float64Array(this.els.length), n = nAtoms(t);
      for (const [e, a] of comp(t)) if (this.ix.has(e)) v[this.ix.get(e)] = a / n;
      const nr = Math.sqrt(v.reduce((s, x) => s + x * x, 0)) || 1;
      return v.map(x => x / nr);
    }
    main(t) { return [...comp(t).keys()].filter(e => !ENV.has(e) && this.ix.has(e)).map(e => this.ix.get(e)); }
    neighbours(q) {
      const v = this.vec(q), mq = new Set(this.main(q)), kk = Math.min(100, this.targets.length - 1);
      const S = this.V.map((u, i) => (this.M[i].some(j => mq.has(j)) ? u.reduce((s, x, j) => s + x * v[j], 0) : -1));
      const top = S.map((s, i) => [s, i]).sort((a, b) => b[0] - a[0] || a[1] - b[1]).slice(0, kk).filter(x => x[0] > 0);
      return [top.map(x => x[1]), top.map(x => x[0] ** 32)];
    }
  }
  function impurityLines(target, ps, fire) {
    const kb = B.impurity;
    if (!K.spI) { K.spI = new Space(Object.keys(kb.by_target).sort()); K.spA = new Space(Object.keys(kb.impure_rate_by_target).sort()); }
    const allowed = new Set(comp(target).keys());
    for (const p of ps) for (const e of comp(p).keys()) allowed.add(e);
    const cand = Object.keys(kb.global).filter(p => [...comp(p).keys()].every(e => allowed.has(e)) && p !== target);
    let [j, w] = K.spI.neighbours(target);
    const local = new Map();
    j.forEach((jj, i) => { for (const [p, n] of Object.entries(kb.by_target[K.spI.targets[jj]])) local.set(p, (local.get(p) || 0) + w[i] * n); });
    const g = cand.map(p => kb.global[p]), n = cand.map(p => local.get(p) || 0);
    const [wg, wn, wl] = kb.weights;
    let s = [], order = [];
    if (cand.length) {
      const gm = Math.max(...g) || 1, nm = Math.max(...n) || 1;
      s = cand.map((p, i) => wg * g[i] / gm + wn * n[i] / nm + wl * +ps.includes(p) + 1e-9 * g[i]);
      order = s.map((_, i) => i).sort((a, b) => s[b] - s[a] || a - b).slice(0, 5);
    }
    const out = [`Likely impurity phases (named as impurities in Lee et al. 2025 for the most similar targets; held out, a named impurity is in the top 3 for ${pct(B.metrics.impurity_hit3)} of reactions):`];
    if (!cand.length) out.push('  no impurity candidates in this chemistry');
    for (const i of order) if (s[i] > 1e-6) out.push(`  ${cand[i]}: named ${Math.trunc(g[i])} times in all, similarity-weighted ${gfmt(n[i], 2)}`);
    [j, w] = K.spA.neighbours(target);
    const rates = kb.impure_rate_by_target, m = kb.impure_model;
    const nb = j.length && npsum(w) > 0 ? npavg(j.map(jj => rates[K.spA.targets[jj]]), w) : m.base_rate;
    const tc = canon(target), h = tc !== null && K.Hc.has(tc) ? K.Hc.get(tc) : null;
    const x = [h !== null ? Math.min(h, 0.5) : m.fill_ehull, +(h === null), [...allowed].filter(e => !ENV.has(e)).length, ps.length, fire, 0, nb];
    let z = 0;
    x.forEach((v, i) => { z += (v - m.mu[i]) / m.sd[i] * m.w[i]; });
    z += m.w[x.length];
    const pr = 1 / (1 + Math.exp(-z));
    out.push(`  chance the paper reports an impurity: ${pct(pr)} (base rate ${pct(m.base_rate)}; held-out AUC ${fixed(B.metrics.impure_auc, 2)}; ` +
      'reporting bias: an impurity nobody looked for is not counted)');
    return out;
  }

  // ---------- the plan (lite.reference) ----------
  function sections(text) {
    const out = [];
    for (const line of text.split('\n')) {
      if (line.startsWith('  ') && out.length) out[out.length - 1].items.push(line.trim());
      else out.push({header: line, items: []});
    }
    return out;
  }
  function plan(target, grams = 10.0, k = 5) {
    try {
      if (!comp(target).size) throw new Error('empty formula');
    } catch (e) {
      return {target, error: `ValueError: ${e.message}`};
    }
    try { return planChecked(target, grams, k); } catch (e) { return {target, error: `${e.name}: ${e.message}`}; }
  }
  function planChecked(target, grams, k) {
    const U = B.uncertainty, T = U.temperature, P1 = K.P1;
    const [sets, probs] = planSets(target, k);
    const df = sets.map(s => P1.drivingForce(target, s));
    const rel = df.map(d => (d === null || df[0] === null ? null : df[0] - d));
    const lines = [`Target: ${target}`, 'Ranked precursor sets (reaction energy vs set 1, eV/atom, MP GGA/GGA+U; negative = more exothermic):'];
    const rowsOut = [];
    sets.forEach((s, i) => {
      lines.push(`  ${i + 1}. ${show(s)}   ` + (rel[i] === null ? 'n/a' : plus(rel[i], 3)));
      rowsOut.push({rank: i + 1, precursors: s, names: show(s), reaction_energy_vs_set1_eV_per_atom: rel[i] === null ? null : pyround(rel[i], 4)});
    });
    const fire = firing(target);
    let sp = spread(target);
    if (Number.isNaN(sp)) sp = T.fill_spread_C;
    const lv = Object.keys(T.q_scaled).sort();
    const rng = lv.map(l => `${Math.trunc(Number(l) * 100)}% range ${fixed(fire - T.q_scaled[l] * (sp + T.floor_C), 0)}-${fixed(fire + T.q_scaled[l] * (sp + T.floor_C), 0)} C`).join('; ');
    const tmp = {fire_C: fire, nn_similarity: null, conformal: {}};
    if (probs) {
      const n = searchsorted(cumsum(probs), U.sets.tau - 1e-12) + 1;
      lines.push(`Prediction set: the first ${n} ranked sets (of ${probs.length}) hold a reported set for ${pct(U.sets.level)} of held-out targets ` +
        `(conformal; test coverage ${pct(B.metrics.sets_coverage)}); set 1 has probability ${pct(probs[0])}`);
    }
    const nb = P1.neighbors(target), sim = nb.length ? nb[0][0] : 0.0;
    if (sim < U.abstain_below_nn_similarity) {
      const a = B.metrics.abstain_mae;
      lines.push(`Abstain: nearest known target has similarity ${fixed(sim, 3)}, below ${fixed(U.abstain_below_nn_similarity, 3)} ` +
        `(10th percentile of held-out targets); held-out temperature error on such targets ${fixed(a[0], 0)} C against ${fixed(a[1], 0)} C for the rest`);
    }
    lines.push(`Suggested maximum heating temperature: ${fixed(fire, 0)} C (${rng}; conformal, scaled by the spread of similar syntheses)`);
    if (sets.length) {
      let w = 0;  // practical.route_lines: every recipe of the solid-state knowledge base is solid-state
      for (const [s, t] of K.P2.weighted(target)) for (const r of t.recipes) w += s;
      if (w) lines.push(`Route of neighbour recipes (weighted): solid-state ${pct(w / w)}`);
      lines.push(...practicalLines(target, sets[0], fire, grams));
      lines.push(...impurityLines(target, sets[0], fire));
    }
    lines.push('Most similar known syntheses:');
    for (const [s, t] of nb.slice(0, 5)) for (const r of t.recipes.slice(0, 2)) {
      lines.push(`  ${t.target} (sim ${fixed(s, 2)}): ${show(r.precursors)}, ${r.temp_C ? fixed(r.temp_C, 0) + ' C' : 'temp n/a'}, doi:${r.doi}`);
    }
    const need = new Set(main(target));
    const tested = B.outcomes.filter(o => { const mo = main(o[0]); return mo.length && subset(mo, need) && mo.filter(e => need.has(e)).length >= 2; });
    if (tested.length) {
      lines.push('Robot-lab experiments in this chemical space (successes and failures):');
      const setKeys = new Set(sets.map(s => s.join('\u0000')));
      const match = tested.filter(o => o[1] && o[1].length && setKeys.has(o[1].join('\u0000')));
      const mk = new Set(match.map(o => JSON.stringify(o)));
      const rest = stableAsc(tested.filter(o => !mk.has(JSON.stringify(o))), o => -(o[3] || 0));
      for (const o of [...match.slice(0, 6), ...rest.slice(0, 3), ...rest.slice(3).slice(-3)]) {
        const y = o[3] === null ? '' : `, yield ${pct(o[3])}`;
        const pr = o[1] && o[1].length ? show(o[1]) : 'precursors not recorded';
        lines.push(`  ${o[0]}: ${pr}, ${o[2] ? fixed(o[2], 0) + ' C' : 'temp n/a'}, ${o[4]}${y} [${o[5]}, ${o[6]}]`);
      }
    }
    const used = new Set(sets.flat());
    const pairs = B.pairs.filter(q => q[3].every(r => used.has(r)));
    if (pairs.length) {
      lines.push('Measured pairwise reactions among these precursors (A-Lab SI Table 2):');
      for (const q of pairs.slice(0, 8)) lines.push(`  ${q[0].join(' + ')} -> ${q[1]} at ${q[2]} C`);
    }
    tmp.nn_similarity = nb.length ? pyround(nb[0][0], 4) : null;
    for (const l of lv) tmp.conformal[`${Math.trunc(Number(l) * 100)}%`] = [fire - T.q_scaled[l] * (sp + T.floor_C), fire + T.q_scaled[l] * (sp + T.floor_C)];
    if (probs) tmp.prediction_set = {n_sets: searchsorted(cumsum(probs), U.sets.tau - 1e-12) + 1, of: probs.length, level: U.sets.level, set1_probability: probs[0]};
    return {target, ranked_sets: rowsOut, fire_C: fire, temperature: tmp, sections: sections(lines.join('\n'))};
  }
  function stableAsc(arr, f) { return arr.map((x, i) => [x, i]).sort((a, b) => f(a[0]) - f(b[0]) || a[1] - b[1]).map(x => x[0]); }

  // ---------- knowledge ----------
  function init(bundle) {
    B = bundle;
    PARSE.clear(); MAIN.clear(); OXI.clear();
    K = {kb: B.kb.map(([target, rs]) => ({target, recipes: rs.map(([precursors, temp_C, doi, cond]) => ({precursors, temp_C, doi, cond}))}))};
    K.vec = K.kb.map(t => fracVec(t.target));
    K.main = K.kb.map(t => main(t.target));
    K.byEl = new Map(); K.multi = new Map(); K.multiByEl = new Map(); K.single = new Map(); K.pair = new Map();
    K.nRecipes = 0; K.temps = [];
    for (const t of K.kb) for (const r of t.recipes) {
      K.nRecipes++;
      if (r.temp_C) K.temps.push(r.temp_C);
      for (const p of r.precursors) {
        const m = main(p);
        if (m.length === 1) {
          if (!K.byEl.has(m[0])) K.byEl.set(m[0], new Map());
          const c = K.byEl.get(m[0]);
          c.set(p, (c.get(p) || 0) + 1);
        } else K.multi.set(p, (K.multi.get(p) || 0) + 1);
      }
      const ps = [...new Set(r.precursors)].sort();
      ps.forEach((a, i) => {
        K.single.set(a, (K.single.get(a) || 0) + 1);
        for (const b of ps.slice(i + 1)) K.pair.set(a + '\u0000' + b, (K.pair.get(a + '\u0000' + b) || 0) + 1);
      });
    }
    for (const [p, n] of K.multi) for (const e of main(p)) {
      if (!K.multiByEl.has(e)) K.multiByEl.set(e, []);
      K.multiByEl.get(e).push([p, n]);
    }
    K.Hc = new Map();
    for (const [f, v] of Object.entries(B.H)) { try { const c = canon(f); if (c !== null) K.Hc.set(c, v); } catch (e) { /* not a formula */ } }
    K.P1 = new Planner(1, 50);
    K.P2 = new Planner(32, 100);
  }
  async function load(url) {
    if (B) return;
    const r = await fetch(url);
    if (!r.ok) throw new Error(`${url}: HTTP ${r.status}`);
    const b = await r.arrayBuffer(), u = new Uint8Array(b);  // unzip unless the server already did (gzip starts 1f 8b)
    init(JSON.parse(u[0] === 31 && u[1] === 139
      ? await new Response(new Blob([b]).stream().pipeThrough(new DecompressionStream('gzip'))).text() : new TextDecoder().decode(b)));
  }
  const api = {init, load, plan, ready: () => !!B, stats: () => ({targets: K.kb.length, recipes: K.nRecipes}), _fmt: {fixed, gfmt, pct, plus}, _comp: f => [...comp(f)], _oxi: f => { const g = oxi(f); return g && [...g]; }};
  if (typeof module !== 'undefined' && module.exports) module.exports = api; else root.IGNITE_LITE = api;
})(typeof self !== 'undefined' ? self : this);
