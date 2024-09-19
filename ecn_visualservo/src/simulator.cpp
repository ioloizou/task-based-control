#include <simulator.h>
#include <log2plot/logger.h>
#include <visp/vpMeterPixelConversion.h>

template <typename... Args> inline void UNUSED(Args&&...) {}

Simulator::Simulator() : config_manager(std::string(BASE_PATH) + "config.yaml")
{   
  log2plot::closePreviousPlots();

  // init images and display
  Iint.init(480, 640, 255);
  Iext.init(480, 640, 255);
  Dint.init(Iint, 100, 100, "Internal view");
  Dext.init(Iext, 100, 100, "External view");
  vpDisplay::setWindowPosition(Iint, 0, 0);
  vpDisplay::setWindowPosition(Iext, 700, 0);
  vpDisplay::display(Iint);
  vpDisplay::flush(Iint);
  vpDisplay::display(Iext);
  vpDisplay::flush(Iext);

  // config file
  // initial pose
  const auto start = config_manager.read<std::string>("startPos");
  const auto end = config_manager.read<std::string>("endPos");
  cMo = config_manager.read<vpHomogeneousMatrix>(start);
  cdMo = config_manager.read<vpHomogeneousMatrix>(end);
  robot.setSamplingTime(dt_ms/1000);
  robot.setPosition(cMo);
  robot.setMaxTranslationVelocity(10);
  robot.setMaxRotationVelocity(10);

  // simulator
  sim.initScene ( vpWireFrameSimulator::SQUARE_10CM, vpWireFrameSimulator::D_STANDARD );
  cam.initPersProjWithoutDistortion(1000, 1000, 320, 240);
  sim.setCurrentViewColor(vpColor::blue);
  sim.setDesiredViewColor(vpColor::red);
  sim.setCameraColor(vpColor::green);
  sim.setCameraPositionRelObj(cMo);
  sim.setDesiredCameraPosition(cdMo);
  sim.setExternalCameraPosition(vpHomogeneousMatrix(0.0, 0, 4.5, vpMath::rad(0), vpMath::rad(-30), 0));
  sim.setInternalCameraParameters(cam);
  sim.setExternalCameraParameters(cam);

  // 3D points
  const double ps = 0.05;
  points.emplace_back ( -ps,-ps,0 );
  points.emplace_back ( -ps,ps,0 );
  points.emplace_back ( ps,ps,0 );
  points.emplace_back ( ps,-ps,0 );

  // CoG
  double x(0), y(0), z(0);
  for(const auto &p: points)
  {
    x += p.get_oX();
    y += p.get_oY();
    z += p.get_oZ();
  }
  center = vpPoint(x/4, y/4, z/4);
}

void Simulator::setVelocity(const vpColVector &v)
{
  if(v.size() != 6)
  {
    std::cerr << "setVelocity: velocity vector is size "
              << v.size() << ", should be 6" << std::endl;
    return;
  }
  robot.setVelocity(vpRobotCamera::CAMERA_FRAME, v);
  vel = v;
  robot.getPosition(cMo);
  sim.setCameraPositionRelObj ( cMo ) ;
  sim.getInternalImage ( Iint );

  pose.buildFrom(cMo);

  uint idx(0);
  for(auto &P: points)
  {
    P.track(cMo);
    uv[2*idx] = computeU(P);
    uv[2*idx+1] = computeV(P);
    idx++;
    vpImagePoint ip;
    vpMeterPixelConversion::convertPoint(cam, P.get_x(), P.get_y(), ip);
    history.push_back(ip);

    vpImagePoint ipd;
    P.track(cdMo);
    vpMeterPixelConversion::convertPoint(cam, P.get_x(), P.get_y(), ipd);
    vpDisplay::displayLine(Iint, ip, ipd, vpColor::darkGray);
  }
  center.track(cMo);
  uv[2*idx] = computeU(center);
  uv[2*idx+1] = computeV(center);
  vpImagePoint ip;
  vpMeterPixelConversion::convertPoint(cam, center.get_x(), center.get_y(), ip);
  history.push_back(ip);

  for(const auto &past: history)
    vpDisplay::displayPoint(Iint, past, vpColor::darkRed, 2);


  vpDisplay::flush ( Iint );
  sim.getExternalImage(Iext);
  vpDisplay::flush ( Iext );
  vpTime::wait(t0, dt_ms);
  t0 = vpTime::measureTimeSecond();

  logger.update();
}

void Simulator::initLog(const std::string &base_path, const std::string &exp_id, const std::string &legend)
{
  std::string log_dir(BASE_PATH);
  const auto start = config_manager.read<std::string>("startPos");
  const auto end = config_manager.read<std::string>("endPos");

  log_dir += "results/" + start + "-" + end + "/" + base_path;
  config_manager.setDirName(log_dir);  
  config_manager.addNameElement(exp_id);

  log_dir = config_manager.fullName();

  logger.setSavePath(log_dir);

  auto rel = log_dir.find("ecn_visualservo");
  std::cout << "Saving to (...)/" << log_dir.substr(rel, log_dir.npos) << "_*" << std::endl;

  t0 = vpTime::measureTimeSecond();
  // 2D points XY (+center)
  // image limits
  const double w(Iint.getCols());
  const double h(Iint.getRows());

  const std::vector<std::vector<size_t>> rectangle{{0,1},{1,2},{2,3},{3,0}};

  uv.resize(10);
  logger.saveXY(uv, "_image", "[P_1, P_2, P_3, P_4, CoG]", "u", "v");
  logger.setLineType("[C0,C1,C2,C3,C4]");
  logger.showFixedShape(log2plot::Shape({{0,0},{0,h},{w,h},{w,0}},rectangle, "k-"));

  std::vector<std::vector<double>> des, cur;
  for(auto &P: points)
  {
    P.track(cdMo);
    des.push_back({computeU(P), computeV(P)});
    P.track(cMo);
    cur.push_back({computeU(P), computeV(P)});
  }

  logger.showFixedShape(log2plot::Shape(des, rectangle, "rs--"));
  logger.showFixedShape(log2plot::Shape(cur, rectangle, "bs--"));
  logger.setPlotArgs("--xLim -10 " + std::to_string(w+10) + " --yLim -10 " + std::to_string(h+10));

  // 3D pose
  logger.save3Dpose(pose, "_pose", "'"+legend+"'", true);
  logger.setLineType("b");
  vpPoseVector pose_d(cdMo.inverse());
  const auto cam{log2plot::Camera("b")};
  logger.showMovingShape(cam);
  logger.showFixedShape(cam.transform(pose_d, "r", "Desired pose"));
  logger.showFixedShape(log2plot::Box(-0.05, -0.05, 0, 0.05, 0.05, 0, "C0d", "Observed points"));


  // velocity
  vel.resize(6);
  logger.save(vel, "_v", "[v_x,v_y,v_z,\\omega_x,\\omega_y,\\omega_z]", "Velocity");
}
