

#include "GeoSphModel/Geometry.hpp"
#include "NumericConcepts/Numeric.hpp"
#include <cstddef>
#include <iostream>

template <NumericConcepts::Real _Real> class Model;

namespace GeoSphModel::Internal {

template <NumericConcepts::Real _Real> struct Traits<Model<_Real>> {
  using Int = std::ptrdiff_t;
  using Real = _Real;
};
} // namespace GeoSphModel::Internal

template <NumericConcepts::Real _Real>
class Model : public GeoSphModel::Geometry<Model<_Real>> {

public:
  using Real = _Real;
  using Base = typename GeoSphModel::Geometry<Model<Real>>;
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

  // Set the mapping and its gradient.
  Real Mapping(Real r, Real theta, Real phi, Int i) const { return 0; }
  Vector MappingGradient(Real r, Real theta, Real phi, Int i) const {
    return Vector(0, 0, 0);
  }

private:
};

int main() {
  auto model = Model<double>();

  auto F = model.DeformationGradient(0.1, 0, 0, 0);
  auto J = model.Jacobian(0.1, 0, 0, 0);

  std::cout << F << std::endl;
  std::cout << J << std::endl;

  std::cout << J * F.inverse() * F.inverse().transpose() << std::endl;
}