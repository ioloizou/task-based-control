#ifndef FEATURES_H
#define FEATURES_H

#include <visp/vpFeatureBuilder.h>
#include <visp/vpFeatureDepth.h>
#include <visp/vpFeatureTranslation.h>
#include <visp/vpFeatureThetaU.h>
#include <visp/vpFeaturePoint.h>
#include <visp/vpFeaturePointPolar.h>
#include <simulator.h>
#include <memory>
#include <vector>

enum class PointDescriptor {XY, Polar, Depth};
enum class TranslationDescriptor {NONE, cTo, cdTc};
enum class RotationDescriptor {NONE, cdRc, cRcd};

class FeatureStack
{
public:
  FeatureStack(Simulator &sim) :
    simulator(sim), cdMo(sim.cdMo),
    z_estim(sim.config_manager.read<double>("z_estim")),
    translation{{}, TranslationDescriptor::NONE},
    rotation{{}, RotationDescriptor::NONE}
  {
    setTranslation3D(config().read<std::string>("translation3D"));
    setRotation3D(config().read<std::string>("rotation3D"));
  }

  const log2plot::ConfigManager& config() const
  {
    return simulator.config_manager;
  }

  void addFeaturePoint(vpPoint P, PointDescriptor descriptor = PointDescriptor::XY);
  void setTranslation3D(std::string descriptor);
  void setRotation3D(std::string descriptor);

  void summary() const;

  void updateFeatures(const vpHomogeneousMatrix &cMo);

  vpColVector s() const  {return s_;}
  vpColVector sd();
  vpMatrix L() const  {return L_;}

protected:

  bool init_done = false;

  Simulator &simulator;

  const vpHomogeneousMatrix& cdMo;
  double z_estim = -1;

  vpColVector s_, sd_, e_;
  std::vector<double> eigvals_;
  uint dim_s = 0;
  vpMatrix L_, L_true_;

  std::vector<std::tuple<vpPoint, const PointDescriptor, double>> points3D;
  std::vector<vpFeaturePoint> pointsXY;
  std::vector<vpFeaturePointPolar> pointsPolar;
  std::vector<vpFeatureDepth> depths;

  template <class Feat, class Descriptor>
  struct Feature3D
  {
    Feat feature;
    Descriptor des;
  };

  Feature3D<vpFeatureTranslation, TranslationDescriptor> translation;
  Feature3D<vpFeatureThetaU, RotationDescriptor> rotation;

  uint update(uint row, vpBasicFeature &f, bool current, vpBasicFeature *f_true = nullptr)
  {
    if(current)
    {
      L_.insert(f.interaction(), row, 0);
      s_.insert(row, f.get_s());

      if(eigvals_.size())
      {
        if(f_true)
          L_true_.insert(f_true->interaction(), row, 0);
        else
          L_true_.insert(f.interaction(), row, 0);
      }
    }
    else
      sd_.insert(row, f.get_s());

    return row + f.dimension_s();
  }

  void computeFeatures(const vpHomogeneousMatrix &cMo, bool current = true);

  void initLog();
};

#endif // FEATURES_H
