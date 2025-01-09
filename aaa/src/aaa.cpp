#define EIGEN_STACK_ALLOCATION_LIMIT 1024
#define EIGEN_NO_DEBUG
#define EIGEN_NO_AUTOMATIC_RESIZING

//#define IO

#include <emscripten/bind.h>
#include <Eigen/Dense>
#include <complex>

#ifdef IO
#include <iostream>
#endif

using namespace emscripten;
using Eigen::Map;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::VectorXcd;
using Eigen::MatrixXcd;
using Eigen::VectorXi;
using Eigen::DiagonalWrapper;
using std::complex;

class VectorC {
  public:
    VectorC(int n): data(n) {
      data.setZero();
    }

    emscripten::val view() const {
      double* ptr = const_cast<double*>(reinterpret_cast<const double*>(data.data()));
      Map<VectorXd> vec(ptr, data.size() * 2);
      return emscripten::val(emscripten::typed_memory_view(vec.size(), vec.data()));
    }

    complex<double> get(int i) const {
      return static_cast<complex<double>>(data(i));
    }

    void set(int i, complex<double> val) {
      data(i) = static_cast<complex<double>>(val);
    }

    size_t length() const {
      return data.size();
    }

#ifdef IO
    std::string toString () {
      std::stringstream ss;
      ss << data;
      return ss.str();
    }
#endif

    VectorXcd data;
};

class VectorD {
  public:
    VectorD(int n): data(n) {
      data.setZero();
    }

    emscripten::val view() const {
      return emscripten::val(emscripten::typed_memory_view(data.size(), data.data()));
    }

    double get(int i) const {
      return data(i);
    }

    void set(int i, double val) {
      data(i) = val;
    }

    size_t length() const {
      return data.size();
    }

#ifdef IO
    std::string toString () {
      std::stringstream ss;
      ss << data;
      return ss.str();
    }
#endif

    Eigen::VectorXd data;
};

class AAAResult {
  public:
    AAAResult() : z(0), f(0), w(0), errvec(0) {}
    VectorC z;
    VectorC f;
    VectorC w;
    VectorD errvec;
};

std::unique_ptr<AAAResult> aaa (const VectorC& Z_, const VectorC& F_, double tol, uint32_t mmax) {
  const VectorXcd& Z = Z_.data;
  const VectorXcd& F = F_.data;
  int m = Z.size();
  std::unique_ptr<AAAResult> result = std::make_unique<AAAResult>();
  VectorXcd& w = result->w.data;
  VectorXcd& z = result->z.data;
  VectorXcd& f = result->f.data;
  VectorXd& errvec = result->errvec.data;
  
  VectorXi J = VectorXi::LinSpaced(m, 0, m - 1);
  MatrixXcd C(m, 0);
  VectorXcd R = Eigen::VectorXcd::Constant(m, F.mean());
  const DiagonalWrapper<const VectorXcd> SF = F.asDiagonal();

  for (int iter = 0; iter < mmax; ++iter) {
    // Select next support point
    Eigen::Index j;
    (F.array() - R.array()).cwiseAbs().maxCoeff(&j);

    // Update support points, data values
    z.conservativeResize(z.size() + 1);
    z(z.size() - 1) = Z(j);

    f.conservativeResize(f.size() + 1);
    f(f.size() - 1) = F(j);

    // Update index vector
    Eigen::Index idx = -1;
    for (Eigen::Index i = 0; i < J.size(); ++i) if (J(i) == j) { idx = i; break; }
    if (idx >= 0 && idx < J.size() - 1) J.segment(idx, J.size() - idx - 1) = J.tail(J.size() - idx - 1);
    J.conservativeResize(J.size() - 1);

    // Next column of Cauchy matrix
    C.conservativeResize(Eigen::NoChange, C.cols() + 1);
    C.col(C.cols() - 1) = (Z.array() - Z(j)).inverse();

    // Right scaling matrix
    Eigen::MatrixXcd Sf = f.asDiagonal().toDenseMatrix();

    // Loewner matrix
    Eigen::MatrixXcd A = SF * C - C * Sf;

    // SVD
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(A(J, Eigen::all), Eigen::ComputeFullV);
    w = svd.matrixV().col(iter);

    // Numerator and denominator
    Eigen::VectorXcd N = C * (w.array() * f.array()).matrix();
    Eigen::VectorXcd D = C * w;

    // Rational approximation
    R = F;
    for (Eigen::Index j = 0; j < J.size(); ++j) R(J(j)) = N(J(j)) / D(J(j));

    // Max error at sample points
    double err = (F - R).lpNorm<Eigen::Infinity>();
    errvec.conservativeResize(errvec.size() + 1);
    errvec(errvec.size() - 1) = err;

    // Stop if converged
#ifdef IO
    std::cout<<"Iteration "<<(iter+1)<<": max error = "<<err<<std::endl;
#endif
    if (err <= tol * F.lpNorm<Eigen::Infinity>()) {
#ifdef IO
      std::cout<<"Converged after "<<(iter + 1)<<" iterations"<<std::endl;
#endif
      break;
    }
  }

  return result;
}

/*
complex<double> feval (const complex<double> zz, std::unique_ptr<AAAResult>& result) {
  const VectorXcd& z = result.z.data;
  const VectorXcd& f = result.f.data;
  const VectorXcd& w = result.w.data;

  VectorXcd zv(1);
  zv(0) = zz;
  MatrixXcd CC = (zv.replicate(1, z.size()).array().colwise() - z.array()).inverse();
  VectorXcd r = (CC * (w.array() * f.array()).matrix()).array() / (CC * w).array();

  for (int i = 0; i < r.size(); ++i) {
    if (std::isnan(r(i).real()) || std::isnan(r(i).imag())) {
      for (int j = 0; j < z.size(); ++j) {
        if (zv(i) == z(j)) {
          r(i) = f(j);
          break;
        }
      }
    }
  }

  return r(0);
  //zv = zz(:);                               % vectorize zz if necessary
  //CC = 1./bsxfun(@minus,zv,z.');            % Cauchy matrix
  //r = (CC*(w.*f))./(CC*w);                  % AAA approx as vector
  //ii = find(isnan(r));                      % find values NaN = Inf/Inf if any
  //for j = 1:length(ii)
    //r(ii(j)) = f(find(zv(ii(j))==z));       % force interpolation there
  //end
  //r = reshape(r,size(zz));    
}
*/


/*class MatrixC {
  public:
    MatrixC(int m, int n): data(m, n) {
      data.setZero();
    }

    emscripten::val view() {
      double* ptr = reinterpret_cast<double*>(data.data());
      Map<MatrixXd> mat(ptr, data.rows(), data.cols() * 2);
      return emscripten::val(emscripten::typed_memory_view(mat.size(), mat.data()));
    }

    WireComplex get(int i, int j) {
      return static_cast<WireComplex>(data(i, j));
    }

    void set(int i, int j, WireComplex val) {
      data(i, j) = static_cast<complex<double>>(val);
    }

  private:
    MatrixXcd data;
};*/



EMSCRIPTEN_BINDINGS(aaa_wrapper) {
  /*class_<MatrixC>("MatrixC")
    .constructor<int, int>()
    .function("get", &MatrixC::get)
    .function("set", &MatrixC::set)
    .function("view", &MatrixC::view)
    ;*/
  class_<VectorC>("VectorC")
    .constructor<int>()
    .function("get", &VectorC::get)
    .function("set", &VectorC::set)
    .property("view", &VectorC::view)
    .property("length", &VectorC::length)
#ifdef IO
    .function("toString", &VectorC::toString)
#endif
    ;
  class_<VectorD>("VectorD")
    .constructor<int>()
    .function("get", &VectorD::get)
    .function("set", &VectorD::set)
    .property("view", &VectorD::view)
    .property("length", &VectorD::length)
#ifdef IO
    .function("toString", &VectorD::toString)
#endif
    ;
  class_<AAAResult>("AAAResult")
    .property("z", &AAAResult::z)
    .property("f", &AAAResult::f)
    .property("w", &AAAResult::w)
    .property("errvec", &AAAResult::errvec)
    ;
  value_array<complex<double>>("complex")
    .element(
      optional_override([](const complex<double>& z) { return z.real(); }),
      optional_override([](complex<double>& z, double value) { z.real(value); })
    )
    .element(
      optional_override([](const complex<double>& z) { return z.imag(); }),
      optional_override([](complex<double>& z, double value) { z.imag(value); })
    )
    ;
  
  function("_aaa", &aaa);
  //function("feval", &feval);
}

