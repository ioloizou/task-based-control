#ifndef VSSIMULATOR_H
#define VSSIMULATOR_H

#include <visp3/core/vpCameraParameters.h>
#include <visp3/core/vpHomogeneousMatrix.h>
#include <visp3/core/vpImage.h>
#include <visp3/core/vpMath.h>
#include <visp3/gui/vpDisplayOpenCV.h>
#include <visp3/io/vpImageIo.h>
#include <visp/vpRobotCamera.h>
#include <visp3/robot/vpWireFrameSimulator.h>
#include <log2plot/config_manager.h>
#include <log2plot/logger.h>
#include <chrono>

class FeatureStack;

class Simulator
{
public:
  Simulator();

  void setVelocity(const vpColVector &v);

  vpHomogeneousMatrix currentPose() const
  {
    return cMo;
  }
  vpHomogeneousMatrix desiredPose() const
  {
    return cdMo;
  }
  vpPoint cog() const
  {
    return center;
  }
  std::vector<vpPoint> observedPoints() const
  {
    return points;
  }

  void plot()
  {
    logger.plot();
    config_manager.saveConfig();
  }

  bool clicked() const
  {
    return  vpDisplay::getClick(Iint,false);
  }

  void waitForClick() const
  {
    std::cout << "Clic on the window to stop" << std::endl;
    vpDisplay::getClick(Iint,true);
  }

protected:
  friend class FeatureStack;
  vpImage<vpRGBa> Iint, Iext;
  vpDisplayOpenCV Dint, Dext;
  vpWireFrameSimulator sim;
  vpCameraParameters cam;
  vpRobotCamera robot;
  std::vector<vpImagePoint> history;

  std::vector<vpPoint> points;
  vpPoint center;
  vpColVector uv, vel;

  // time sampling
  const double dt_ms = 0.01 * 1000;  // ms
  double t0;


  // from configuration
  vpHomogeneousMatrix cMo, cdMo;
  vpPoseVector pose;
  log2plot::ConfigManager config_manager;
  log2plot::Logger logger;

  void initLog(const std::string &base_path, const std::string &exp_id, const std::string &legend);
  double computeV(const vpPoint &P) const
  {
    return Iint.getRows() - (cam.get_v0() + cam.get_py()*P.get_y());
  }
  double computeU(const vpPoint &P) const
  {
    return cam.get_u0() + cam.get_px()*P.get_x();
  }
};

#endif // VSSIMULATOR_H
