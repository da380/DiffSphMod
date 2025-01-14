#pragma once

#include "Dimensions/Dimensions.hpp"
#include "NumericConcepts/Numeric.hpp"
#include "Utility.hpp"
#include <Eigen/Dense>
#include <cassert>
#include <cmath>
#include <functional>
#include <ranges>
#include <utility>

namespace GeoSphModel {

namespace Internal {

template <typename _Derived> struct Traits {};
} // namespace Internal

template <typename _Derived>
class Geometry : public Dimensions::Dimensions<Geometry<_Derived>> {

public:
  // Typedefs from traits
  using Int = typename Internal::Traits<_Derived>::Int;
  using Real = typename Internal::Traits<_Derived>::Real;
  using Vector = Eigen::Matrix<Real, 3, 1>;
  using Matrix = Eigen::Matrix<Real, 3, 3>;

  //  Methods for the Dimension class.
  constexpr auto LengthScale() const { return Derived().LengthScale(); }
  constexpr auto DensityScale() const { return Derived().DensityScale(); }
  constexpr auto TimeScale() const { return Derived().TimeScale(); }
  constexpr auto TemperatureScale() const {
    return Derived().TemperatureScale();
  }

  // Number of layers in the model.
  Int NumberOfLayers() const { return Derived().NumberOfLayers(); };

  // Return the radius of the ith referential layer boundary.
  Real ReferentialBoundaryRadius(Int i) const {
    return Derived().ReferentialBoundaryRadius();
  }

  // The radial , h, is such that in the ith layer the referential point
  // (r, \theta, \phi) is taken to (r + h(r,\theta,\phi,i), \theta, \phi).
  Real RadialMapping(Real r, Real theta, Real phi, Int i) const {
    return Derived().RadialMapping(r, theta, phi, i);
  }

  // The gradient of the mapping, h, defined relative to the spherical
  // polar unit vectors.
  Vector RadialMappingGradient(Real r, Real theta, Real phi, Int i) const {
    return Derived().RadialMappingGradient(r, theta, phi, i);
  }

  //-----------------------------------------------//
  //                Induced methods                //
  //-----------------------------------------------//

  // Returns the number of boundaries.
  auto NumberOfBoundaries() const { return NumberOfLayers() + 1; }

  // Returns a range over the layer Indices.
  auto LayerIndices() const {
    return std::ranges::views::iota(0, NumberOfLayers());
  }

  // Return a range of the referential layer radii
  auto ReferentialBoundaryRadii() const {
    return LayerIndices() | std::ranges::views::transform([this](auto i) {
             return ReferentialBoundaryRadius(i);
           });
  }

  // Return the boundary radii for the ith layer.
  auto ReferentialBoundaryRadiiForLayer(Int i) const {
    assert(i < NumberOfLayers());
    return std::pair{ReferentialBoundaryRadius(i),
                     ReferentialBoundaryRadius(i + 1)};
  }

  // Return the deformation gradient in the ith layer relative to the polar
  // basis in both the referential and spatial manifolds.
  Matrix DeformationGradient(Real r, Real theta, Real phi, Int i) const {
    Real ri = r > 0 ? static_cast<Real>(1) / r : 0;
    auto I = IdentityMatrix<Real>();
    auto h = RadialMapping(r, theta, phi, i);
    auto dh = RadialMappingGradient(r, theta, phi, i);
    auto x = PolarToCartesian(r, theta, phi);
    return (1 + h * ri) * I +
           x * (dh.transpose() - h * x.transpose() * ri * ri) * ri;
  }

  Matrix InverseDeformationGradient(Real r, Real theta, Real phi, Int i) const {
    Real ri = r > 0 ? static_cast<Real>(1) / r : 0;
    auto I = IdentityMatrix<Real>();
    auto h = RadialMapping(r, theta, phi, i);
    auto dh = RadialMappingGradient(r, theta, phi, i);
    auto x = PolarToCartesian(r, theta, phi);
    return (I - x * ri * (dh.transpose() - h * x.transpose() * ri * ri) /
                    (1 + dh(0))) /
           std::pow(1 + h * ri, 2);
  }

  // Return the Jacobian in the ith layer.
  Real Jacobian(Real r, Real theta, Real phi, Int i) const {
    Real ri = r > 0 ? static_cast<Real>(1) / r : 0;
    auto h = RadialMapping(r, theta, phi, i);
    auto dh = RadialMappingGradient(r, theta, phi, i);
    return std::pow(1 + h * ri, 2) * (1 + dh(0));
  }

private:
  constexpr auto &Derived() { return static_cast<_Derived &>(*this); }

  constexpr auto &Derived() const {
    return static_cast<const _Derived &>(*this);
  }
};

} // namespace GeoSphModel