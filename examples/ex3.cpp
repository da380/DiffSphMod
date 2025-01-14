

#include "GeoSphModel/GeoSphModel.hpp"
#include "NumericConcepts/Numeric.hpp"
#include <cstddef>
#include <iostream>

template <NumericConcepts::Real _Real> class SphericalViscoelastic;

namespace GeoSphModel::Internal {

template <NumericConcepts::Real _Real>
struct Traits<SphericalViscoelastic<_Real>> {
  using Int = std::ptrdiff_t;
  using Real = _Real;
};
} // namespace GeoSphModel::Internal

template <NumericConcepts::Real _Real>
class SphericalViscoelastic
    : public GeoSphModel::Viscoelastic<SphericalViscoelastic<_Real>> {

public:
  using Base = typename GeoSphModel::Viscoelastic<SphericalViscoelastic<_Real>>;
  using Int = typename Base::Int;
  using Real = typename Base::Real;
  using Complex = typename Base::Complex;
  using Vector = typename Base::Vector;
  using Matrix = typename Base::Matrix;

  // Define the scale-lengths.
  constexpr auto LengthScale() const { return static_cast<Real>(1); }
  constexpr auto DensityScale() const { return static_cast<Real>(1); }
  constexpr auto TimeScale() const { return static_cast<Real>(1); }
  constexpr auto TemperatureScale() const { return static_cast<Real>(1); }

  // Set the layering methods.
  Int NumberOfLayers() const { return 1; }
  Real ReferentialBoundaryRadii(Int i) const { return i == 0 ? 0 : 1; }

  // Set the radial mapping and its gradient.
  Real RadialMapping(Real r, Real theta, Real phi, Int i) const {
    return GeoSphModel::IdentityRadialMapping(r, theta, phi, i);
  }
  Vector RadialMappingGradient(Real r, Real theta, Real phi, Int i) const {
    return GeoSphModel::IdentityRadialMappingGradient(r, theta, phi, i);
  }

  // Set the referential density function.
  Real ReferentialDensity(Real r, Real theta, Real phi, Int i) const {
    return 1;
  }

  // Set the viscoelastic modulus
  Complex ReferentialViscoelasticModulus(Real r, Real theta, Real phi,
                                         Real omega, Int i, Int j, Int k, Int l,
                                         Int layer) const {
    using GeoSphModel::Delta;
    return _lambda * Delta(i, j) * Delta(k, l) +
           _mu * (Delta(i, k) * Delta(j, l) + Delta(i, l) * Delta(j, k));
  }

private:
  Real _lambda = 1;
  Real _mu = 1;
};

int main() {
  auto model = SphericalViscoelastic<double>();

  auto F = model.DeformationGradient(0.1, 0, 0, 0);
  auto J = model.Jacobian(0.1, 0, 0, 0);

  std::cout << F << std::endl;
  std::cout << J << std::endl;

  std::cout << model.ReferentialDensity(0.1, 0, 0, 0) << std::endl;
  std::cout << model.ReferentialViscoelasticModulus(0.1, 0, 0, 0, 1, 1, 1, 1, 0)
            << std::endl;
}