#include <feature_stack.h>
#include <Eigen/Eigenvalues>

void recap(const std::string s, int n, int dim)
{
  if(n)
    std::cout << s << ": " << n << " (dim " << dim*n << ")" << std::endl;
}

void FeatureStack::summary() const
{
  std::cout << "Total feature dimension: " << dim_s << std::endl;
  recap(" - Points XY", pointsXY.size(), 2);
  recap(" - Points Polar", pointsPolar.size(), 2);
  recap(" - Points depths", depths.size(), 1);
  if(translation.des != TranslationDescriptor::NONE)
    recap(" - 3D translations", 1, 3);
  if(rotation.des != RotationDescriptor::NONE)
    recap(" - 3D rotations", 1, 3);
}

vpColVector FeatureStack::sd()
{
  if(sd_.size() == 0)
  {
    sd_.resize(dim_s);
    computeFeatures(cdMo, false);
  }
  return sd_;
}

void FeatureStack::updateFeatures(const vpHomogeneousMatrix &cMo)
{
  if(!init_done)
    initLog();
  computeFeatures(cMo, true);
  e_ = s_-sd();
  if(eigvals_.size())
  {
    static Eigen::MatrixXd LLp(dim_s, dim_s);
    static vpMatrix LLp_visp;
    LLp_visp = L_true_*L_.pseudoInverse(1e-15);

    for(uint row=0; row < dim_s; ++row)
    {
      for(uint col = 0; col < dim_s; ++col)
      {
        LLp(row, col) = LLp_visp[row][col];
      }
    }

    const auto eigv{LLp.eigenvalues()};

    for(uint row = 0; row < dim_s; ++row)
      eigvals_[row] = eigv(row).real();
    std::sort(eigvals_.begin(), eigvals_.end());
  }
}

void FeatureStack::computeFeatures(const vpHomogeneousMatrix &cMo, bool current)
{
  uint row(0);
  if(translation.des != TranslationDescriptor::NONE)
  {
    if(translation.des == TranslationDescriptor::cTo)
      translation.feature.buildFrom(cMo);
    else
      translation.feature.buildFrom(cdMo*cMo.inverse());
    row = update(row, translation.feature, current);
  }

  if(rotation.des != RotationDescriptor::NONE)
  {
    if(rotation.des == RotationDescriptor::cRcd)
      rotation.feature.buildFrom(cMo*cdMo.inverse());
    else
      rotation.feature.buildFrom(cdMo*cMo.inverse());
    row = update(row, rotation.feature, current);
  }

  auto xy(pointsXY.begin());
  auto rt(pointsPolar.begin());
  auto d(depths.begin());
  for(auto &[P, descriptor, zd]: points3D)
  {
    P.track(cMo);
    const double z(P.get_Z());
    const auto P_true{P};
    // Z-estimation
    if(current)
    {
      if(z_estim > 0)
        P.set_Z(z_estim);
      else if(z_estim == 0)
        P.set_Z(zd);
    }
    if(descriptor == PointDescriptor::XY)
    {
      xy->buildFrom(P.get_x(), P.get_y(), P.get_Z());
      static vpFeaturePoint xy_true;
      vpFeatureBuilder::create(xy_true, P_true);
      row = update(row, *xy, current, &xy_true);
      xy++;
    }
    else if(descriptor == PointDescriptor::Polar)
    {
      vpFeatureBuilder::create(*rt, P);
      static vpFeaturePointPolar rt_true;
      vpFeatureBuilder::create(rt_true, P_true);
      row = update(row, *rt, current, &rt_true);
      rt++;
    }
    else
    {
      d->buildFrom(P.get_x(), P.get_y(), P.get_Z(), log(z/zd));
      static vpFeatureDepth d_true;
      d_true.buildFrom(P_true.get_x(), P_true.get_y(), P_true.get_Z(), log(P_true.get_Z()/zd));
      row = update(row, *d, current, &d_true);
      d++;
    }
  }
}

void FeatureStack::addFeaturePoint(vpPoint P, PointDescriptor descriptor)
{
  if(descriptor == PointDescriptor::XY)
  {
    pointsXY.push_back(vpFeaturePoint());
    dim_s += 2;
  }
  else if(descriptor == PointDescriptor::Polar)
  {
    pointsPolar.push_back(vpFeaturePointPolar());
    dim_s += 2;
  }
  else
  {
    depths.push_back(vpFeatureDepth());
    dim_s += 1;
  }
  P.track(cdMo);
  points3D.push_back({P, descriptor, P.get_Z()});
}

void FeatureStack::setTranslation3D(std::string descriptor)
{
  if(descriptor == "cTo")
  {
    translation.feature.setFeatureTranslationType(vpFeatureTranslation::cMo);    
    if(translation.des == TranslationDescriptor::NONE)
          dim_s += 3;
    translation.des = TranslationDescriptor::cTo;
  }
  else if(descriptor == "cdTc")
  {
    translation.feature.setFeatureTranslationType(vpFeatureTranslation::cdMc);
    if(translation.des == TranslationDescriptor::NONE)
          dim_s += 3;
    translation.des = TranslationDescriptor::cdTc;
  }
}

void FeatureStack::setRotation3D(std::string descriptor)
{
  if(descriptor == "cdRc")
  {
    rotation.feature.setFeatureThetaURotationType(vpFeatureThetaU::cdRc);
    if(rotation.des == RotationDescriptor::NONE)
          dim_s += 3;
    rotation.des = RotationDescriptor::cdRc;
  }
  else if(descriptor == "cRcd")
  {
    rotation.feature.setFeatureThetaURotationType(vpFeatureThetaU::cRcd);
    if(rotation.des == RotationDescriptor::NONE)
          dim_s += 3;
    rotation.des = RotationDescriptor::cRcd;
  }
}

void FeatureStack::initLog()
{
  init_done = true;

  // build sub-dir and legend
  std::stringstream ss;
  std::string base_path, legend3D;
  auto updatePath([&]()
  {
    const auto key(ss.str());
    auto legend_key(key);
    /* if(legend_key == "cdRc")
                          legend_key = "{}^{c*}\\mathbf{R}_c";
                        else if(legend_key == "cRcd")
                          legend_key = "{}^c\\mathbf{R}_{c*}";
                        else if(legend_key == "cTo")
                            legend_key = "{}^c\\mathbf{T}_o";
                        else if(legend_key == "cdTc")
                          legend_key = "{}^{c*}\\mathbf{T}_c";*/
    if(base_path == "")
    {
      base_path = key;
      legend3D = key;
    }
    else
    {
      base_path += "_" + key;
      legend3D += "+" + legend_key;
    }
    ss.str("");
  });

  std::stringstream legend;
  legend << "[";
  if(pointsXY.size())
  {
    ss << "XY" << pointsXY.size();
    updatePath();
    for(size_t i = 0; i < pointsXY.size(); ++i)
      legend << "'x_{" << i+1 << "}', 'y_{" << i+1 << "}',";
  }
  if(pointsPolar.size())
  {
    ss << "Polar" << pointsPolar.size();
    updatePath();
    for(size_t i = 0; i < pointsPolar.size(); ++i)
      legend << "'\\rho_{" << i+1 << "}', '\\theta_{" << i+1 << "}',";
  }
  if(depths.size())
  {
    ss << "Depth" << depths.size();
    updatePath();
    for(size_t i = 0; i < depths.size(); ++i)
      legend << "'\\log Z_{" << i+1 << "}/Z^*',";
  }
  if(translation.des != TranslationDescriptor::NONE)
  {
    ss << (translation.des == TranslationDescriptor::cTo ? "cTo" : "cdTc");
    updatePath();
    legend << "t_x, t_y, t_z,";
  }
  if(rotation.des != RotationDescriptor::NONE)
  {
    ss << (rotation.des == RotationDescriptor::cRcd ? "cRcd" : "cdRc");
    updatePath();
    legend << "\\theta u_x, \\theta u_y, \\theta u_z,";
  }


  // build ID for this experiment
  if(pointsXY.size() || pointsPolar.size() || depths.size())
  {
    if(z_estim == 0)
      ss << "Zd";
    else if(z_estim>0)
      ss << "Z" << z_estim;
    else
      ss << "Zc";
    ss << "_";
  }
  ss << "lambda" << simulator.config_manager.read<double>("lambda");

  simulator.initLog(base_path, ss.str(), legend3D);

  // plot feature error
  s_.resize(dim_s, false);
  e_.resize(dim_s, false);
  L_.resize(dim_s, 6, false);

  if(e_.size())
  {
    std::string legend_s(legend.str());
    legend_s.back() = ']';
    simulator.logger.save(e_, "_err", legend_s, "Feature error");

    if(dim_s > 6 || ((pointsXY.size() + pointsPolar.size() || depths.size()) && z_estim >= 0))
    {
      L_true_.resize(dim_s, 6, false);
      eigvals_.resize(dim_s);
      simulator.logger.save(eigvals_, "_eigvals", "\\lambda_", "<Eigen values of >L\\hat{L}^+");
    }
  }
}
