#pragma once

#include "NumericConcepts/Numeric.hpp"
#include <Eigen/Dense>
#include <cmath>

namespace DiffSphMod {

// Identity matrix.
template <NumericConcepts::Real Real>
Eigen::Matrix<Real, 3, 3> IdentityMatrix() {
  return Eigen::Matrix<Real, 3, 3>{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
}

// Convert polar to cartesian coordinates.
template <NumericConcepts::Real Real>
auto PolarToCartesian(Real r, Real theta, Real phi) {
  return Eigen::Matrix<Real, 3, 1>(r * std::sin(theta) * std::cos(phi),
                                   r * std::sin(theta) * std::sin(phi),
                                   r * std::cos(theta));
}

// Convert cartesian to polar coordinates.
template <NumericConcepts::Real Real>
auto CartesianToPolar(const Eigen::Matrix<Real, 3, 1> &x) {
  auto r = x.norm();
  auto theta = std::atan2(x(2), r);
  auto phi = std::atan2(x(0), x(1));
  return std::tuple{r, theta, phi};
}

// Implementation of the identity radial mapping and its gradient for
// convenience.
template <NumericConcepts::Integral Int, NumericConcepts::Real Real>
Real IdentityRadialMapping(Real r, Real theta, Real phi, Int i) {
  return 0;
}
template <NumericConcepts::Integral Int, NumericConcepts::Real Real>
Eigen::Matrix<Real, 3, 1> IdentityRadialMappingGradient(Real r, Real theta,
                                                        Real phi, Int i) {
  return Eigen::Matrix<Real, 3, 1>(0, 0, 0);
}
} // namespace DiffSphMod