

#include "DiffSphMod/DiffSphMod.hpp"
#include "NumericConcepts/Numeric.hpp"
#include <array>
#include <cstddef>
#include <iostream>

template <NumericConcepts::Real _Real> class SphericalGeometry;

namespace DiffSphMod::Internal {

template <NumericConcepts::Real _Real> struct Traits<SphericalGeometry<_Real>> {
  using Real = _Real;
};
} // namespace DiffSphMod::Internal

template <NumericConcepts::Real _Real>
class SphericalGeometry
    : public DiffSphMod::Geometry<SphericalGeometry<_Real>> {

public:
  using Base = typename DiffSphMod::Geometry<SphericalGeometry<_Real>>;
  using Real = typename Base::Real;
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
  Real ReferentialBoundaryRadius(Int i) const { return _radii[i]; }

  // Set the radial mapping and its gradient.
  Real RadialMapping(Real r, Real theta, Real phi, Int i) const {
    return 2 * r + 0.1 * r * r;
  }

  Vector RadialMappingGradient(Real r, Real theta, Real phi, Int i) const {
    return Vector(2 + 0.2 * r, 0, 0);
  }

private:
  std::array<Real, 2> _radii{0, 1};
};

int main() {
  auto model = SphericalGeometry<double>();

  auto r = 0.2;
  auto theta = 0.;
  auto phi = 0.;

  auto s = model.RadialMapping(r, theta, phi, 0);

  auto t = model.InverseRadialMapping(s, theta, phi, 0);

  std::cout << r << " " << t << std::endl;
}
