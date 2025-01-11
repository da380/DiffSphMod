#pragma once

#include "Dimensions/Dimensions.hpp"
#include "NumericConcepts/Numeric.hpp"
#include <cassert>
#include <functional>
#include <ranges>
#include <utility>

namespace RadialModel {

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

  // Return radial mapping, \xi,  in the ith layer. This is such that
  // the point (r, \theta, \phi) in the reference model is mapped to the
  // point (\xi(r,ztheta, \phi), \theta, \phi) in the physical model.
  std::function<Real(Real, Real, Real)> RadialMapping(Int i) const {
    assert(i < NumberOfLayers());
    return Derived().RadialMapping(i);
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

  // Return a function that giving the spatial radius on the ith boundary.
  auto SpatialBoundaryRadiusFunction(Int i) const {
    assert(i < NumberOfBoundaries());
    return [this, i](auto theta, auto phi) {
      auto xi = RadialMapping(i);
      auto r = ReferentialBoundaryRadius(i);
      return xi(r, theta, phi);
    };
  }

  // Return functions giving the boundary radii for the ith layer.
  auto SpatialBoundaryRadiusFunctionsForLayer(Int i) const {
    assert(i < NumberOfLayers());
    return std::pair{SpatialBoundaryRadiusFunction(i),
                     SpatialBoundaryRadiusFunction(i + 1)};
  }

private:
  constexpr auto &Derived() { return static_cast<_Derived &>(*this); }

  constexpr auto &Derived() const {
    return static_cast<const _Derived &>(*this);
  }
};
} // namespace RadialModel