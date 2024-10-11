#include <visp/vpFeaturePoint.h>
#include <ecn_baxter_vs/baxter_arm.h>
#include <visp/vpSubMatrix.h>
#include <visp/vpSubColVector.h>
#include <ecn_common/visp_utils.h>

using namespace std;

int main(int argc, char** argv)
{
  BaxterArm arm(argc, argv, "right");    // defaults to right arm
  //  arm.detect(r, g, b, show, saturation, value); to pick detected color, otherwise default

  vpColVector q = arm.init();

  vpColVector q_4_7(3);
  q_4_7[0] = q[4];
  q_4_7[1] = q[5];
  q_4_7[2] = q[6];

  vpColVector qmin = arm.jointMin(), qmax = arm.jointMax();

  // define a simple 2D point feature and its desired value
  vpFeaturePoint p,pd;
  pd.set_xyZ(0,0,1);

  // the error
  vpColVector e(11);
  double x, y, area;
  // desired area
  const double area_d = arm.area_d();

  // loop variables
  vpColVector qdot(7);
  vpMatrix L(4, 6), Js(4,7);

  vpMatrix L_area(1,6);


  while(arm.ok())
  {
    std::cout << L << "-------------" << endl;

    auto bRc = arm.cameraPose().getRotationMatrix();
    double s_ro = bRc[2][0];

    // get point features
    x = arm.x();
    y = arm.y();
    area = arm.area();
    p.set_xyZ(x,y, 1);
    std::cout << "x: " << x << ", y: " << y << ", area: " << area << '\n';

    // update error vector e
    e[0] = p.get_x() - pd.get_x();
    e[1] = p.get_y() - pd.get_y();
    e[2] = area - area_d;
    e[3] = s_ro - 0;
    // update interaction matrix L
    vpMatrix L_k = p.interaction();

    ecn::putAt(L, L_k, 0, 0);

    // Make the interaction matrix for area
    L[2][2] = 1;

    // build H matrix (2nd section) using arm.rho()
    q = arm.jointPosition();

    // Perpendicular
    L[3][4] = -bRc[2][2];
    L[3][5] = bRc[2][1];

    // H matrix
    auto rho = arm.rho();

    auto q_plus = arm.jointMax();
    auto q_min = arm.jointMin();

    auto q_s_plus = q_plus - rho*(q_plus - q_min);
    auto q_s_min = q_min + rho*(q_plus - q_min);

    auto q_ref = 1/2*(q_plus + q_min);

    vpMatrix H(11,11);
    vpMatrix eye_4(4,4);

    eye_4.setIdentity();
    ecn::putAt(H, eye_4, 0, 0);

    for (int i=0; i<4; i++){
        auto h_max = ecn::weight(q[i], q_s_plus[i], q_plus[i]);
        auto h_min = ecn::weight(-q[i], -q_s_min[i], -q_min[i]);
        H[4+i][4+i] = h_max + h_min;
    }

    for (int i=0; i<3; i++){
        H[8+i][8+i] = 1;
    }


    // compute feature Jacobian from L and cameraJacobian
    const auto J = arm.cameraJacobian(q);

    Js = L*J;

    // Defining J extented with weights
    vpMatrix eye_n;
    eye_n.setIdentity(7);
    vpMatrix J_extented(11,7);

    ecn::putAt(J_extented, Js, 0, 0);
    ecn::putAt(J_extented, eye_n, 4, 0);

    // Defining error extented
    auto e_q = q - q_ref;

    ecn::putAt(e, e_q, 4);

    // But we need joint 3 to 7 to have init values as desired to not get rank deficient
    e[4] = q[4] - q_4_7[0];
    e[5] = q[5] - q_4_7[1];
    e[6] = q[6] - q_4_7[2];


    qdot = -arm.lambda()*(H*J_extented).pseudoInverse()*H*e;
    // send this command to the robot
    arm.setJointVelocity(qdot);

    // display current joint positions and VS error
    arm.plot(e);
  }
}
