import AAA from "./aaa.js";
await AAA.ready();

const { VectorComplex, aaa } = AAA;
const { PI, sin, cos, sinh, cosh, exp } = Math;

// Z = exp(linspace(-.5, .5+15i*pi, 1000));
const N = 1000;
const z = [...Array(N).keys()].map((i) => {
  const [tr, ti] = [-0.5 + i / (N - 1), 15.0 * PI * (i / (N - 1))];
  return [exp(tr) * cos(ti), exp(tr) * sin(ti)];
});

// F = tan(pi*z/2);
const f = z.map(([zr, zi]) => [
  sin(PI * zr) / (cosh(PI * zi) + cos(PI * zr)),
  sinh(PI * zi) / (cosh(PI * zi) + cos(PI * zr)),
]);

// Compute the AAA approximation
const approx = aaa(z, f, 1e-13, 100);

// Evaluate the result
const [fr, fi] = approx.eval([0.5, 0]);
console.log(`tan(π/4) ~ ${fr} + ${fi}i`);

// Clean up
approx.delete();
