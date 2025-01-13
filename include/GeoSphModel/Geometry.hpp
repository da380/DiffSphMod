#pragma once

#include "Dimensions/Dimensions.hpp"
#include "NumericConcepts/Numeric.hpp"
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
  constexpr auto MassScale() const { return static_cast<Real>(1); }

public:
  // Typedefs from traits
  using Int = typename Internal::Traits<_Derived>::Int;
  using Real = typename Internal::Traits<_Derived>::Real;
  using Vector = Eigen::Matrix<Real, 3, 1>;
  using Matrix = Eigen::Matrix<Real, 3, 3>;

  //------------------------------------------------//
  //      Defined methods in the derived class      //
  //------------------------------------------------//

  // Related to Dimensions class.
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
    assert(i < NumberOfBoundaries());
    return Derived().ReferentialBoundaryRadius();
  }

  Real Mapping(Real r, Real theta, Real phi, Int i) const {
    return Derived().Mapping(r, theta, phi, i);
  }

  Vector MappingGradient(Real r, Real theta, Real phi, Int i) const {
    return Derived().MappingGradient(r, theta, phi, i);
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
    auto I = IdentityMatrix();
    auto h = Mapping(r, theta, phi, i);
    auto dh = MappingGradient(r, theta, phi, i);
    auto x = PolarToCartesian(r, theta, phi);
    return (1 + h * ri) * I +
           x * (dh.transpose() - h * x.transpose() * ri * ri) * ri;
  }

  Matrix InverseDeformationGradient(Real r, Real theta, Real phi, Int i) const {
    Real ri = r > 0 ? static_cast<Real>(1) / r : 0;
    auto I = IdentityMatrix();
    auto h = Mapping(r, theta, phi, i);
    auto dh = MappingGradient(r, theta, phi, i);
    auto x = PolarToCartesian(r, theta, phi);
    return (I - x * ri * (dh.transpose() - h * x.transpose() * ri * ri) /
                    (1 + dh(0))) /
           std::pow(1 + h * ri, 2);
  }

  // Return the Jacobian in the ith layer.
  Real Jacobian(Real r, Real theta, Real phi, Int i) const {
    Real ri = r > 0 ? static_cast<Real>(1) / r : 0;
    auto h = Mapping(r, theta, phi, i);
    auto dh = MappingGradient(r, theta, phi, i);
    return std::pow(1 + h * ri, 2) * (1 + dh(0));
  }

private:
  constexpr auto &Derived() { return static_cast<_Derived &>(*this); }

  constexpr auto &Derived() const {
    return static_cast<const _Derived &>(*this);
  }

  // Identity matrix.
  Matrix IdentityMatrix() const {
    return Matrix{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
  }

  // Convert polar to cartesian coordinates.
  auto PolarToCartesian(Real r, Real theta, Real phi) const {
    return Vector(r * std::sin(theta) * std::cos(phi),
                  r * std::sin(theta) * std::sin(phi), r * std::cos(theta));
  }

  // Convert cartesian to polar coordinates.
  auto CartesianToPolar(const Vector &x) const {
    auto r = x.norm();
    auto theta = std::atan2(x(2), r);
    auto phi = std::atan2(x(0), x(1));
    return std::tuple{r, theta, phi};
  }
};
} // namespace GeoSphModel