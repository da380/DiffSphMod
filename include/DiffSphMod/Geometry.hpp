#pragma once

#include "Dimensions/Dimensions.hpp"
#include "Traits.hpp"
#include "Utility.hpp"
#include <Eigen/Core>
#include <boost/math/tools/roots.hpp>
#include <cstddef>
#include <iostream>
#include <limits>
#include <ranges>

namespace DiffSphMod {

template <typename _Derived>
class Geometry : public Dimensions::Dimensions<Geometry<_Derived>> {

public:
  // Typedefs from traits
  using Int = std::ptrdiff_t;
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
    return Derived().ReferentialBoundaryRadius(i);
  }

  // The radial , h, is such that in the ith layer the referential point
  // (r, \theta, \phi) is taken to (h(r,\theta,\phi,i), \theta, \phi).
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

  // Returns a range over the boundary Indices.
  auto BoundaryIndices() const {
    return std::ranges::views::iota(0, NumberOfBoundaries());
  }

  // Return a range of the referential layer radii
  auto ReferentialBoundaryRadii() const {
    return LayerIndices() | std::ranges::views::transform([this](auto i) {
             return ReferentialBoundaryRadius(i);
           });
  }

  // Return the boundary radii for the ith layer.
  auto ReferentialBoundaryRadiiForLayer(Int i) const {
    return std::pair<Real, Real>{ReferentialBoundaryRadius(i),
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
    // return h * ri * I + x * (dh.transpose() - h * x.transpose() * ri * ri) *
    // ri;
    return I;
  }

  Matrix InverseDeformationGradient(Real r, Real theta, Real phi, Int i) const {
    Real ri = r > 0 ? static_cast<Real>(1) / r : 0;
    auto I = IdentityMatrix<Real>();
    auto h = RadialMapping(r, theta, phi, i);
    auto dh = RadialMappingGradient(r, theta, phi, i);
    auto x = PolarToCartesian(r, theta, phi);
    // return (I - x * ri * (dh.transpose() - h * x.transpose() * ri * ri) /
    //                 (1 + dh(0))) /
    //        std::pow(1 + h * ri, 2);
    return I;
  }

  // Return the Jacobian in the ith layer.
  Real Jacobian(Real r, Real theta, Real phi, Int i) const {
    Real ri = r > 0 ? static_cast<Real>(1) / r : 0;
    auto h = RadialMapping(r, theta, phi, i);
    auto dh = RadialMappingGradient(r, theta, phi, i);
    return std::pow(h * ri, 2) * dh(0);
  }

  // Return the inverse radial mapping in the ith layer.
  Real InverseRadialMapping(Real s, Real theta, Real phi, Int i) const {

    // Set the radial bounds.
    auto [r0, r1] = ReferentialBoundaryRadiiForLayer(i);
    auto s0 = RadialMapping(r0, theta, phi, i);
    auto s1 = RadialMapping(r1, theta, phi, i);

    // Function for root finding.
    auto func = [this, theta, phi, i, s](auto r) {
      auto h = RadialMapping(r, theta, phi, i) - s;
      auto dh = RadialMappingGradient(r, theta, phi, i)(0);
      return std::pair{h, dh};
    };

    // Solve with Newton Raphson.
    return boost::math::tools::newton_raphson_iterate(
        func, 0.5 * (r0 + r1), r0, r1, std::numeric_limits<Real>::digits);
  }

private:
  constexpr auto &Derived() { return static_cast<_Derived &>(*this); }

  constexpr auto &Derived() const {
    return static_cast<const _Derived &>(*this);
  }
};

} // namespace DiffSphMod