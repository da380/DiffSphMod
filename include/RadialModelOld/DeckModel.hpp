#pragma once

#include "Dimensions.hpp"
#include "Interpolation/CubicSpline.hpp"
#include "NumericConcepts/Numeric.hpp"
#include "RadialModel.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace RadialModel {

template <NumericConcepts::Real Real>
class DeckModel : public RadialModel<Real> {
private:
  using Int = std::ptrdiff_t;

  // Store the header information.
  std::vector<std::string> _header;

  // Model values at the radial knots.
  std::vector<Real> _r;
  std::vector<Real> _density;
  std::vector<Real> _LoveModulusA;
  std::vector<Real> _LoveModulusC;
  std::vector<Real> _LoveModulusF;
  std::vector<Real> _LoveModulusL;
  std::vector<Real> _LoveModulusN;
  std::vector<Real> _bulkQualityFactor;
  std::vector<Real> _shearQualityFactor;

  // Layering information.
  std::vector<Int> _boundaryIndex;
  std::vector<Real> _boundaryRadius;
  std::vector<bool> _layerSolid;

  // Cubic spline interpolating functions within each layer.
  using Spline =
      Interpolation::CubicSpline<typename std::vector<Real>::iterator,
                                 typename std::vector<Real>::iterator>;

  std::vector<Spline> _densities;
  std::vector<Spline> _LoveModuliiA;
  std::vector<Spline> _LoveModuliiC;
  std::vector<Spline> _LoveModuliiF;
  std::vector<Spline> _LoveModuliiL;
  std::vector<Spline> _LoveModuliiN;
  std::vector<Spline> _bulkQualityFactors;
  std::vector<Spline> _shearQualityFactors;

public:
  DeckModel(const Dimensions<Real> &dimensions, const std::string &fileName);

  // Return the number of layers. Override of pure virtual function in
  // base class.
  Int NumberOfLayers() const override;

  // Return the bounding radii of the ith layer. Override of pure
  // virtual function in base class.
  std::pair<Real, Real> LayerRadii(Int i) const override;

  // Return true if ith layer is solid. Override of pure virtual function in
  // base class.
  bool LayerIsSolid(Int i) const override;

  // Return the number of knots in the model.
  Int NumberOfKnots() const;

  // Functions to return material parameters functions.
  std::function<Real(Real, Int)> Density() const override;
  std::function<Real(Real, Int)> LoveModulusA() const override;
  std::function<Real(Real, Int)> LoveModulusC() const override;
  std::function<Real(Real, Int)> LoveModulusF() const override;
  std::function<Real(Real, Int)> LoveModulusL() const override;
  std::function<Real(Real, Int)> LoveModulusN() const override;
  std::function<Real(Real, Int)> BulkQualityFactor() const override;
  std::function<Real(Real, Int)> ShearQualityFactor() const override;
};

} // namespace RadialModel