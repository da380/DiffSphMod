#pragma once

#include "IsotropicElastic.hpp"
#include "Utility.hpp"

namespace GeoSphModel {

template <typename _Derived> class SphericalIsotropicElastic;

namespace Internal {
template <typename _Derived>
struct Traits<SphericalIsotropicElastic<_Derived>> {
  using Int = typename Traits<_Derived>::Int;
  using Real = typename Traits<_Derived>::Real;
};
} // namespace Internal

template <typename _Derived>
class SphericalIsotropicElastic
    : public IsotropicElastic<SphericalIsotropicElastic<_Derived>> {
public:
  // Typedefs from traits
  using Int = typename Internal::Traits<_Derived>::Int;
  using Real = typename Internal::Traits<_Derived>::Real;
  using Complex = std::complex<Real>;
  using Vector = Eigen::Matrix<Real, 3, 1>;
  using Matrix = Eigen::Matrix<Real, 3, 3>;

  //-------------------------------------------------------------------//
  //                 Methods defined in derived class                  //
  //-------------------------------------------------------------------//

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
  Real BoundaryRadius(Int i) const { return Derived().BoundaryRadius(i); }

  // Return the  density in the ith layer.
  Real Density(Real r, Int i) const { return Derived().Density(r, i); }

  // Return the bulk modulus in the ith layer.
  Real BulkModulus(Real r, Int i) const { return Derived().BulkModulus(r, i); }

  // Return the shear modulus in the ith layer.
  Real ShearModulus(Real r, Int i) const {
    return Derived().ShearModulus(r, i);
  }

  //-------------------------------------------------------------------//
  //                          Induced methods                          //
  //-------------------------------------------------------------------//

  // Return the radius of the ith referential layer boundary.
  Real ReferentialBoundaryRadius(Int i) const { return BoundaryRadius(i); }

  // The mapping, h, is such that in the ith layer the referential point (r,
  // \theta, \phi) is taken to (r + h(r,\theta,\phi,i), \theta, \phi).
  Real RadialMapping(Real r, Real theta, Real phi, Int i) const {
    return IdentityRadialMapping(r, theta, phi, i);
  }

  // The gradient of the mapping, h, defined relative to the spherical
  // polar unit vectors.
  Vector RadialMappingGradient(Real r, Real theta, Real phi, Int i) const {
    return IdentityRadialMappingGradient(r, theta, phi, i);
  }

  // Return the referential density in the ith layer.
  Real ReferentialDensity(Real r, Real theta, Real phi, Int i) const {
    return Density(r, i);
  }

  // Return the bulk modulus in the ith layer.
  Real ReferentialBulkModulus(Real r, Real theta, Real phi, Int i) const {
    return BulkModulus(r, i);
  }

  // Return the shear modulus in the ith layer.
  Real ReferentialShearModulus(Real r, Real theta, Real phi, Int i) const {
    return ShearModulus(r, i);
  }

private:
  constexpr auto &Derived() { return static_cast<_Derived &>(*this); }

  constexpr auto &Derived() const {
    return static_cast<const _Derived &>(*this);
  }
};

} // namespace GeoSphModel