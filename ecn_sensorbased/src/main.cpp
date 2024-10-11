#include <ecn_sensorbased/pioneer_cam.h>
#include <visp/vpFeaturePoint.h>
#include <visp/vpQuadProg.h>
#include <ecn_common/visp_utils.h>

using namespace std;


int main(int argc, char** argv)
{
    ros::init(argc, argv, "control_node");
    PioneerCam robot;

    // pose error gain
    const double lv = .5;
    // constraints gain
    const double lc = 2;
    geometry_msgs::Pose2D target;

    // robot command u = (v, omega, dot q_p, dot q_t)
    vpColVector u(4);

    // QP solver
    vpQuadProg qp;

    // Cost Function part of the QP
    int n=2;
    vpMatrix Q(n,4);
    Q[0][0] = 1;
    Q[1][1] = 1;
    vpColVector r(n);

    // Inequality part
    double radius=robot.radius();
    double base=robot.base();
    vpMatrix C(8,4);
    vpColVector d(8);

    C[0][0] = 1/radius;
    C[0][1] = base/radius;
    C[1][0] = 1/radius;
    C[1][1] = -base/radius;
    C[2][0] = -1/radius;
    C[2][1] = -base/radius;
    C[3][0] = -1/radius;
    C[3][1] = +base/radius;

    d[0] = robot.wmax();
    d[1] = robot.wmax();
    d[2] = robot.wmax();
    d[3] = robot.wmax();

    while(ros::ok())
    {
        cout << "-------------" << endl;

        if(robot.stateReceived())
        {
            // get robot and target positions to get position error
            target = robot.targetRelativePose();

            // linear velocity
            u[0] = lv*(target.x - .1);
            // angular velocity
            u[1] = lv*std::atan2(target.y, target.x);

            // Visual Servoing Constraint

            r[0] = u[0];
            r[1] = u[1];

            vpFeaturePoint sphere_cog = robot.imagePoint();
            vpColVector cam_bounds = robot.camLimits();
            vpMatrix cJc = robot.camJacobian();

            vpMatrix L = sphere_cog.interaction();

            // Size 2x4
            auto Js = L*cJc;

            ecn::putAt(C, Js, 4, 0);
            ecn::putAt(C, -Js, 6, 0);

            double a = 1;

            vpColVector sphere_cog_vec(2,1);
            sphere_cog_vec[0] = sphere_cog.get_x();
            sphere_cog_vec[1] = sphere_cog.get_y();

            vpColVector offset(2);
            offset[0]=0.5;
            offset[1]=0.5;

            vpColVector splus_side = a*(cam_bounds - sphere_cog_vec - offset);
            vpColVector sminus_side = -a*(cam_bounds - sphere_cog_vec - offset);
            ecn::putAt(d, splus_side, 4);
            ecn::putAt(d, sminus_side, 6);

            // Equality part
            vpMatrix A(1,4);
            A[0][0] = 1;
            A[0][1] = -r[0]/r[1];
            vpMatrix b(1,1);
            b[0][0] = 0;

            qp.solveQP( Q, r, A, b, C, d, u);
//            qp.solveQPi( Q, r, C, d, u);


            cout << "u: " << u.t() << endl;

            robot.sendVelocity(u);
        }
    }
}
