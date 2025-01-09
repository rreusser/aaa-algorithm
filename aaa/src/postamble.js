{
  const _ready = new Promise((resolve, reject) => {
    Module["onRuntimeInitialized"] = function () {
      const VectorC = Module["VectorC"];
      const VectorD = Module["VectorD"];

      VectorC["fromArray"] = function (array) {
        const v = new VectorC(array.length / 2);
        v.view.set(array);
        return v;
      };

      VectorD["fromArray"] = function (array) {
        const v = new VectorD(array.length);
        v.view.set(array);
        return v;
      };

      Module["aaa"] = function (z, f, tol = 1e-13, maxIter = 100) {
        if (z.length !== f.length) {
          throw new Error("z and f must have the same length");
        }
        const cleanup = [];

        if (!(z instanceof VectorC)) {
          z = VectorC.fromArray(z);
          cleanup.push(z);
        }
        if (!(f instanceof VectorC)) {
          f = VectorC.fromArray(f);
          cleanup.push(f);
        }

        const result = Module["_aaa"].call(Module, z, f, tol, maxIter);

        for (const obj of cleanup) obj.delete();

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
