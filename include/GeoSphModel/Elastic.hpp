#pragma once

#include "Viscoelastic.hpp"

namespace GeoSphModel {

template <typename _Derived> class Elastic;

namespace Internal {
template <typename _Derived> struct Traits<Elastic<_Derived>> {
  using Int = typename Traits<_Derived>::Int;
  using Real = typename Traits<_Derived>::Real;
};
} // namespace Internal

template <typename _Derived>
class Elastic : public Viscoelastic<Elastic<_Derived>> {
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

  // Return a component of the viscoelastic tensor in the ith layer.
  Real ReferentialSecondElasticTensor(Real r, Real theta, Real phi, Int i,
                                      Int j, Int k, Int l, Int layer) const {
    return Derived().ReferentialSecondElasticTensor(r, theta, phi, i, j, k, l,
                                                    layer);
  }

  //-------------------------------------------------------------------//
  //                          Induced methods                          //
  //-------------------------------------------------------------------//
  // Return a component of the viscoelastic tensor in the ith layer.
  Complex ReferentialSecondViscoelasticTensor(Real r, Real theta, Real phi,
                                              Real omega, Int i, Int j, Int k,
                                              Int l, Int layer) const {
    return ReferentialSecondElasticTensor(r, theta, phi, i, j, k, l, layer);
  }

private:
  constexpr auto &Derived() { return static_cast<_Derived &>(*this); }

  constexpr auto &Derived() const {
    return static_cast<const _Derived &>(*this);
  }
};

} // namespace GeoSphModel