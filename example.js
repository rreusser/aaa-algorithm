import AAA from "./aaa.js";
await AAA.ready();

const { VectorC, aaa } = AAA;
const { PI, sin, cos, sinh, cosh, exp } = Math;

const N = 1000;
const z = new VectorC(N);
const f = new VectorC(N);

for (let i = 0; i < N; i++) {
  const tr = -0.5 + i / (N - 1);
  const ti = 15.0 * PI * (i / (N - 1));
  const a = exp(tr) * cos(ti);
  const b = exp(tr) * sin(ti);

  z.set(i, [a, b]);
  f.set(i, [
    sin(PI * a) / (cosh(PI * b) + cos(PI * a)),
    sinh(PI * b) / (cosh(PI * b) + cos(PI * a)),
  ]);
}

const result = aaa(z, f, 1e-13, 100);

console.log(result.z.view);

z.delete();
f.delete();
result.delete();
