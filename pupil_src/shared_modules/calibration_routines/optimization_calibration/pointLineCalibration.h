


#include "common.h"
#include <vector>
#include <cstdio>
#include <limits>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Geometry>

using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::CauchyLoss;
using ceres::CostFunction;
using ceres::LossFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;

struct TransformationError {
    TransformationError(const Vector3 referencePoint,   const Vector3 gazePoint1, const Vector3 gazePoint2)
        : referencePoint(referencePoint), gazePoint1(gazePoint1), gazePoint2(gazePoint2) {}

    template <typename T>
    bool operator()(
        const T* const R,  // Rotation denoted by quaternion
        const T* const t,
        // followed with translation
        T* residuals) const
    {


        // Compute coordinates with current transformation matrix: x = RX + t.


        T g1[3] = {T(gazePoint1[0]), T(gazePoint1[1]), T(gazePoint1[2])};
        T g2[3] = {T(gazePoint2[0]), T(gazePoint2[1]), T(gazePoint2[2])};

        T p1[3];
        T p2[3];
        ceres::QuaternionRotatePoint(R, g1, p1 );
        ceres::QuaternionRotatePoint(R, g2, p2 );

        p1[0] += t[0];
        p1[1] += t[1];
        p1[2] += t[2];
        p2[0] += t[0];
        p2[1] += t[1];
        p2[2] += t[2];

        Eigen::Matrix<T, 3, 1> refP;
        refP << T(referencePoint[0]) , T(referencePoint[1]) , T(referencePoint[2]);

        Eigen::Matrix<T, 3, 1> ep1;
        ep1 << T(p1[0]) , T(p1[1]) , T(p1[2]);

        Eigen::Matrix<T, 3, 1> ep2;
        ep2 << T(p2[0]) , T(p2[1]) , T(p2[2]);

        // Equation 3: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
        // delta tells us if the point lies in direction "behind" or "infront" of ep1
        T delta = -(ep1 - refP ).dot(ep2 - ep1) / (ep2 - ep1).squaredNorm();

        // in our case the point should alway lay on the ray from p1 to p2
        if(  delta  >= 0.0 ){
            // now calculate the distance between the observed point and the nearest point on the line
            // Equation 10: http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
            // and divide by delta
            // by dividing by delta we actually optimize the Sine of the angle between these two lines
            residuals[0] = ((refP - ep1).cross(refP - ep2).squaredNorm() / (ep2 - ep1).squaredNorm()) / (delta*delta);
            return true;

        }
        return false;


    }

    const Vector3 referencePoint;
    const Vector3 gazePoint1;
    const Vector3 gazePoint2;
};


void pointLineCalibration(Vector3 spherePosition, const std::vector<Vector3>& refPoints, const std::vector<Vector3>& gazeDirections , double* orientation , double* translation )
{


    Problem problem;

    ceres::LocalParameterization* quaternion_parameterization = new ceres::QuaternionParameterization;

    // don't use Constructor 'Quaternion (const Scalar *data)' because the internal layout for coefficients is different from the one we use.
    // Memory Layout EIGEN: xyzw
    // Memory Layout CERES and the one we use: wxyz
    Eigen::Quaterniond q(orientation[0],orientation[1],orientation[2],orientation[3]);

    int i = 0;

    for (auto& p : refPoints) {
        auto g = gazeDirections.at(i);
        i++;

        g.normalize();

        // do a check to handle parameters we can't solve
        // First: the length of the line must be greater zero
        // Second: the angle between line direction and reference point direction must not be greater 90 degrees, considering the initial orientation
        auto v2 = q*g;
        std::cout << "v2: " << v2 << std::endl;
        std::cout << "dot: " << v2.dot(p) << std::endl;
        if( g.norm() >= std::numeric_limits<double>::epsilon() && v2.dot(p) > 0.0   ){
            CostFunction* cost = new AutoDiffCostFunction<TransformationError , 1, 4, 3 >(new TransformationError(p , Vector3::Zero() , g ));
            problem.AddResidualBlock(cost, nullptr, orientation,  translation);
        }else{
            std::cout << "no valid direction vector"  << std::endl;
        }
    }


    if( problem.NumResidualBlocks() == 0 ){

        std::cout << "nothing to solve"  << std::endl;
        return;
    }

    problem.SetParameterBlockConstant(translation);
    problem.SetParameterization(orientation, quaternion_parameterization);


    // Build and solve the problem.
    Solver::Options options;
    options.max_num_iterations = 1000;
    options.linear_solver_type = ceres::DENSE_QR;
    //options.parameter_tolerance = 1e-14;
    options.function_tolerance = 1e-10;
    options.gradient_tolerance = 1e-20;
    //options.minimizer_type = ceres::LINE_SEARCH;

    //options.use_nonmonotonic_steps = true;
    //options.minimizer_progress_to_stdout = true;
    //options.check_gradients = true;
    Solver::Summary summary;
    Solve(options, &problem, &summary);
    // // Recover r from m.

    // std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";
   for (int i = 0; i < 4; i++) {
        std::cout << "," << orientation[i] << std::endl;
    }
    for (int i = 0; i < 3; i++) {
        std::cout << "," << translation[i] << std::endl;
    }


    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> r;
    Eigen::Matrix<double, 3, 3> rc;
    ceres::QuaternionToRotation( orientation , r.data() );
    std::cout << "r: " << r << std::endl;
    std::cout << "det: " << r.determinant() << std::endl; // be sure we get a valid rotation matrix, det(R) == 1
    rc = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(r.data());


    // we need to take the sphere position into account
    // thus the actual translation is not right, because the local cordinate frame of the eye need to be translated in the opposite direction
    // of the sphere coordinates

    // since the actual translation is in world coordinates, the sphere translation needs to be calculated in world coordinates
    Eigen::Matrix4d eyeToWorld =  Eigen::Matrix4d::Identity();
    eyeToWorld.block<3,3>(0,0) = rc;
    eyeToWorld(0, 3) = translation[0];
    eyeToWorld(1, 3) = translation[1];
    eyeToWorld(2, 3) = translation[2];
    std::cout << "eyeToWorld: "  << eyeToWorld  << std::endl;

    Eigen::Vector4d sphereWorld = eyeToWorld * Eigen::Vector4d(spherePosition[0],spherePosition[1],spherePosition[2], 1.0 );
    std::cout << "sphere world: " << sphereWorld << std::endl;
    Vector3 sphereOffset =  sphereWorld.head<3>() - Vector3(translation);
    Vector3 actualtranslation =  Vector3(translation) - sphereOffset;
    translation[0] = actualtranslation[0];
    translation[1] = actualtranslation[1];
    translation[2] = actualtranslation[2];

    // Eigen::Matrix4d tt = Eigen::Map<Eigen::Matrix<double,4,4,Eigen::RowMajor> >(t.data());
    // std::cout << "transformation: "  << tt  << std::endl;

}
