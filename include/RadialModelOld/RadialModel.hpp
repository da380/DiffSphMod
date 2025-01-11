#pragma once

#include <algorithm>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <format>
#include <fstream>
#include <functional>
#include <iostream>
#include <numbers>
#include <numeric>
#include <ranges>
#include <string>
#include <vector>

#include "Dimensions.hpp"
#include "NumericConcepts/Numeric.hpp"

namespace RadialModel {

template <NumericConcepts::Real Real>
class RadialModel : public Dimensions<Real> {
  using Int = std::ptrdiff_t;

public:
  using LengthScale = Dimension<Real>::LengthScale;
  // Constructor to be called from derived classes.
  RadialModel(const Dimensions<Real> &dimensions)
      : Dimensions<Real>(dimensions) {}

  // Return the number of layers.
  virtual Int NumberOfLayers() const = 0;

  // Return the bounding radii of the ith layer.
  virtual std::pair<Real, Real> LayerRadii(Int i) const = 0;

  // Return true if the ith layer is solid.
  virtual bool LayerIsSolid(Int i) const = 0;

  // Returns the number of boundaries.
  auto NumberOfBoundaries() const { return NumberOfLayers() + 1; }

  // Returns a range over the layer Indices.
  auto LayerIndices() const {
    return std::ranges::views::iota(0, NumberOfLayers());
  }

  void f() const { std::cout << LengthScale() << std::endl; }

  // Return the surface radius of the planet.
  auto SurfaceRadius() const { return LayerRadii(NumberOfLayers() - 1).second; }

  // Return the Jean length for a given degree.
  auto JeanLength(Int degree) const {
    return 2 * std::numbers::pi_v<Real> * SurfaceRadius() /
           static_cast<Real>(degree + 1);
  }

  // Write the model out in deck format.
  void WriteAsDeckModel(const std::string &fileName,
                        const std::array<std::string, 3> &header,
                        Real maximumKnotSpacing) const {
    // Open the file for output.
    auto modelFile = std::ofstream(fileName);

    // Write the header information.
    for (auto line : header) {
      modelFile << line << "\n";
    }

    // Loop over the layers
    for (auto i = 0; i < NumberOfLayers(); i++) {
      auto attribute = i + 1;
      auto [r0, r1] = LayerRadii(i);
      auto n = std::max<Int>(2, std::round((r1 - r0) / maximumKnotSpacing));
      auto dr = (r1 - r0) / (n - 1);
      for (auto j = 0; j < n; j++) {
        auto r = r0 + j * dr;
        /*
              modelFile << std::format(
                  "{:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f} "
                  "{:>11.3f} {:>11.3f} {:>11.3f} {:>11.3f}\n",
                  r * LengthScale(), Density()(r, attribute) * DensityScale(),
                  VerticalPVelocity()(r, attribute) * VelocityScale(),
                  VerticalSVelocity()(r, attribute) * VelocityScale(),
                  BulkQualityFactor()(r, attribute),
                  ShearQualityFactor()(r, attribute),
                  HorizontalPVelocity()(r, attribute) * VelocityScale(),
                  HorizontalSVelocity()(r, attribute) * VelocityScale(),
                  AnisotropicEtaParameter()(r, attribute));
      */
      }
    }
  }

  // Write the model out in deck format using default header lines.
  void WriteAsDeckModel(const std::string &fileName,
                        Real maximumKnotSpacing) const;

  // Material parameter functions for override.
  virtual std::function<Real(Real, Int)> Density() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusA() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusC() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusF() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusL() const = 0;
  virtual std::function<Real(Real, Int)> LoveModulusN() const = 0;
  virtual std::function<Real(Real, Int)> BulkQualityFactor() const = 0;
  virtual std::function<Real(Real, Int)> ShearQualityFactor() const = 0;

  // Derived material parameter functions.
  std::function<Real(Real, Int)> VerticalPVelocity() const;
  std::function<Real(Real, Int)> VerticalSVelocity() const;
  std::function<Real(Real, Int)> HorizontalPVelocity() const;
  std::function<Real(Real, Int)> HorizontalSVelocity() const;
  std::function<Real(Real, Int)> AnisotropicEtaParameter() const;
  std::function<Real(Real, Int)> BulkModulus() const;
  std::function<Real(Real, Int)> ShearModulus() const;
};

} // namespace RadialModel