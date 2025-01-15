#pragma once

#include "IsotropicViscoelastic.hpp"
#include <complex>

namespace GeoSphModel {

template <typename _Derived> class MaxwellViscoelastic;

namespace Internal {
template <typename _Derived> struct Traits<MaxwellViscoelastic<_Derived>> {
  using Int = typename Traits<_Derived>::Int;
  using Real = typename Traits<_Derived>::Real;
};
} // namespace Internal

template <typename _Derived>
class MaxwellViscoelastic
    : public IsotropicViscoelastic<MaxwellViscoelastic<_Derived>> {
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
  Real ReferentialBoundaryRadius(Int i) const {
    return Derived().ReferentialBoundaryRadius(i);
  }

  // The mapping, h, is such that in the ith layer the referential point (r,
  // \theta, \phi) is taken to (r + h(r,\theta,\phi,i), \theta, \phi).
  Real RadialMapping(Real r, Real theta, Real phi, Int i) const {
    return Derived().RadialMapping(r, theta, phi, i);
  }

  // The gradient of the mapping, h, defined relative to the spherical
  // polar unit vectors.
  Vector RadialMappingGradient(Real r, Real theta, Real phi, Int i) const {
    return Derived().RadialMappingGradient(r, theta, phi, i);
  }

  // Return the referential density in the ith layer.
  Real ReferentialDensity(Real r, Real theta, Real phi, Int i) const {
    return Derived().ReferentialDensity(r, theta, phi, i);
  }

  // Return the bulk modulus in the ith layer.
  Real ReferentialBulkModulus(Real r, Real theta, Real phi, Int i) const {
    return Derived().ReferentialBulkModulus(r, theta, phi, i);
  }

  // Return the shear modulus in the ith layer.
  Real ReferentialShearModulus(Real r, Real theta, Real phi, Int i) const {
    return Derived().ReferentialShearModulus(r, theta, phi, i);
  }

  // Return the Maxwell time in the ith layer.
  Real ReferentialMaxwellTime(Real r, Real theta, Real phi, Int i) const {
    return Derived().ReferentialMaxwellTime(r, theta, phi, i);
  }

  //-------------------------------------------------------------------//
  //                          Induced methods                          //
  //-------------------------------------------------------------------//

  // Return the viscoelastic bulk modulus in the ith layer.
  Complex ReferentialViscoelasticBulkModulus(Real r, Real theta, Real phi,
                                             Real omega, Int i) const {
    return ReferentialBulkModulus(r, theta, phi, i);
  }

  // Return the viscoelastic shear modulus in the ith layer.
  Complex ReferentialViscoelasticShearModulus(Real r, Real theta, Real phi,
                                              Real omega, Int i) const {
    auto mu = ReferentialShearModulus(r, theta, phi, i);
    auto tau = ReferentialMaxwellTime(r, theta, phi, i);
    constexpr auto ii = std::complex<Real>{0, 1};
    return mu * (ii * omega * tau) / (ii * omega * tau + 1);
  }

  // Return the viscosity in the ith layer.
  Real ReferentialViscosity(Real r, Real theta, Real phi, Int i) const {
    return ReferentialShearModulus(r, theta, phi, i) *
           ReferentialMaxwellTime(r, theta, phi, i);
  }

private:
  constexpr auto &Derived() { return static_cast<_Derived &>(*this); }

  constexpr auto &Derived() const {
    return static_cast<const _Derived &>(*this);
  }
};

} // namespace GeoSphModel