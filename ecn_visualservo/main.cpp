#include <simulator.h>
#include <feature_stack.h>

using std::string;

int main()
{
  // the main simulator
  Simulator sim;

  // to group the features into a single vector
  FeatureStack stack(sim);

  // configuration file handle
  const auto config{stack.config()};

  // get considered 2D features from configuration
  // 3D features are loaded if not "none"
  const auto useXY{config.read<bool>("useXY")};
  const auto usePolar{config.read<bool>("usePolar")};
  const auto use2Half{config.read<bool>("use2Half")};

  // tuning
  const auto err_min{config.read<double>("errMin")};
  const auto lambda{config.read<double>("lambda")};
  const auto iter_max{config.read<uint>("iterMax")};

  if(useXY)
  {
    for(auto point: sim.observedPoints())
      stack.addFeaturePoint(point, PointDescriptor::XY);
  }

  if(usePolar)
  {
    for(auto point: sim.observedPoints())
      stack.addFeaturePoint(point, PointDescriptor::Polar);
  }

  // TODO add specific features that are used when doing 2 1/2 half VS
  if(use2Half)
  {
     stack.addFeaturePoint(sim.cog(), PointDescriptor::XY);
     stack.addFeaturePoint(sim.cog(), PointDescriptor::Depth);
     stack.setRotation3D("cdRc");
  }

  // 3D features are added anyway
  stack.setTranslation3D(config.read<string>("translation3D"));
  stack.setRotation3D(config.read<string>("rotation3D"));

  stack.summary();

  // loop variables
  uint iter(0);
  const auto sd{stack.sd()};
  vpColVector s(sd.size(), err_min);

  // main control loop
  while(iter++ < iter_max && (s-sd).frobeniusNorm() > err_min && !sim.clicked())
  {
    // update stack features from current simulation pose
    // this comment is useless, just read the code
    stack.updateFeatures(sim.currentPose());

    // TODO get the current features and their interaction matrix
    vpMatrix L = stack.L();
    vpColVector s = stack.s();
    vpColVector sd = stack.sd();



    // TODO compute velocity twist and send it to the simulation
    vpColVector v(6);

    v = -lambda*L.pseudoInverse()*(s - sd);

    sim.setVelocity(v);
  }

  // wait for a last clic before exiting
  sim.waitForClick();

  sim.plot();
}
