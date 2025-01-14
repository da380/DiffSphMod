

#include "GeoSphModel/GeoSphModel.hpp"
#include "NumericConcepts/Numeric.hpp"
#include <cstddef>
#include <iostream>

template <NumericConcepts::Real _Real> class SphericalDensity;

namespace GeoSphModel::Internal {

template <NumericConcepts::Real _Real> struct Traits<SphericalDensity<_Real>> {
  using Int = std::ptrdiff_t;
  using Real = _Real;
};
} // namespace GeoSphModel::Internal

template <NumericConcepts::Real _Real>
class SphericalDensity : public GeoSphModel::Density<SphericalDensity<_Real>> {

public:
  using Real = _Real;
  using Base = typename GeoSphModel::Density<SphericalDensity<Real>>;
  using Int = typename Base::Int;
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

private:
};

int main() {
  auto model = SphericalDensity<double>();

  auto F = model.DeformationGradient(0.1, 0, 0, 0);
  auto J = model.Jacobian(0.1, 0, 0, 0);

  std::cout << F << std::endl;
  std::cout << J << std::endl;

  std::cout << model.ReferentialDensity(0.1, 0, 0, 0) << std::endl;
}