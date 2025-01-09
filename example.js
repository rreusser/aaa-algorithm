import AAA from "./aaa.js";
await AAA.ready();

const { VectorC, aaa } = AAA;
const { PI, sin, cos, sinh, cosh, exp } = Math;

const N = 1000;
const z = new VectorC(N);
const f = new VectorC(N);

for (let i = 0; i < N; i++) {
  // Z = exp(linspace(-.5, .5+15i*pi, 1000));
  const [tr, ti] = [-0.5 + i / (N - 1), 15.0 * PI * (i / (N - 1))];
  const [zr, zi] = [exp(tr) * cos(ti), exp(tr) * sin(ti)];
  z.set(i, [zr, zi]);

  // F = tan(pi*z/2);
  f.set(i, [
    sin(PI * zr) / (cosh(PI * zi) + cos(PI * zr)),
    sinh(PI * zi) / (cosh(PI * zi) + cos(PI * zr)),
  ]);
}

const result = aaa(z, f, 1e-13, 100);

for (let i = 0; i < result.z.length; i++) {
  const [zr, zi] = result.z.get(i);
  console.log(`z[${i}] = ${zr} + ${zi}i`);
}

z.delete();
f.delete();
result.delete();
