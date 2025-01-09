{
  const _ready = new Promise((resolve, reject) => {
    Module["onRuntimeInitialized"] = function () {
      const VectorComplex = Module["VectorComplex"];
      const VectorDouble = Module["VectorDouble"];

      function flatten(array) {
        if (Array.isArray(array)) return array.flat();
        return array;
      }

      VectorComplex["fromArray"] = function (array) {
        array = flatten(array);
        const v = new VectorComplex(array.length / 2);
        v.view.set(array);
        return v;
      };

      VectorDouble["fromArray"] = function (array) {
        array = flatten(array);
        const v = new VectorDouble(array.length);
        v.view.set(array);
        return v;
      };

      Module["aaa"] = function (z, f, tol = 1e-13, maxIter = 100) {
        if (z.length !== f.length) {
          throw new Error("z and f must have the same length");
        }
        const cleanup = [];

        if (!(z instanceof VectorComplex)) {
          z = VectorComplex.fromArray(z);
          cleanup.push(z);
        }
        if (!(f instanceof VectorComplex)) {
          f = VectorComplex.fromArray(f);
          cleanup.push(f);
        }

        const approx = Module["_aaa"].call(Module, z, f, tol, maxIter);

        while (cleanup.length) cleanup.pop().delete();

        return approx;
      };

      Module["eval"] = function (z, approx) {
        const cleanup = [];
        if (!(z instanceof VectorComplex)) {
          z = VectorComplex["fromArray"](z);
          cleanup.push(z);
        }
        const result = Module["_eval"].call(Module, z, approx);
        while (cleanup.length) cleanup.pop().delete();
        return result;
      };

      resolve();
    };
    Module["onAbort"] = function () {
      reject();
    };
  });
  Module["ready"] = () => _ready;
}
module.exports = Module;
