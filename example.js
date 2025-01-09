import AAA from "./aaa.js";
await AAA.ready();

const { VectorComplex, aaa } = AAA;
const { PI, sin, cos, sinh, cosh, exp } = Math;

const N = 1000;
const z = new VectorComplex(N);
const f = new VectorComplex(N);

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

// Compute the AAA approximation
const approx = aaa(z, f, 1e-13, 100);

// Evaluate the result at some points
const zEval = VectorComplex.fromArray([
  [-0.5, 0],
  [0, 0],
  [0.5, 0],
]);
const fEval = AAA.eval(zEval, approx);
for (let i = 0; i < zEval.length; i++) {
  const z = zEval.get(i);
  const f = fEval.get(i);
  console.log(`f(${z[0]} + ${z[1]}i) = ${f[0]} + ${f[1]}i`);
}

// Clean up
z.delete();
f.delete();
approx.delete();
fEval.delete();
zEval.delete();
