


#include "common.h"
#include <vector>
#include <cstdio>
#include <limits>

#include <Eigen/Geometry>
#include <opencv2/calib3d.hpp>
#include <opencv2/core/eigen.hpp>
#include "utils.h"

// since for rayray distance both parameters s,t for eq r0=p0+s*d0 and r1=p1+t*d1 need to be positive
// ceres get's to know if the rays don't lie in the same direction with angle less than 90 degree
// Book: Geometric Tools for Computer Graphics, Side 413
template<typename Scalar, int Dim>
bool rayRayIntersectionInFront(const Eigen::ParametrizedLine<Scalar, Dim>& ray0, const Eigen::ParametrizedLine<Scalar, Dim>& ray1 )
{

    typedef typename Eigen::ParametrizedLine<Scalar, Dim>::VectorType Vector;

    Result<Scalar> result;
    result.valid = false;

    Vector diff = ray0.origin() - ray1.origin();
    Scalar a01 = - ray0.direction().dot(ray1.direction());
    Scalar b0 = diff.dot(ray0.direction());
    Scalar b1;
    Scalar s0, s1;

    if (ceres::abs(a01) < Scalar(1) )
    {
        // Rays are not parallel.
        b1 = -diff.dot(ray1.direction());
        s0 = a01 * b1 - b0;
        s1 = a01 * b0 - b1;

        if (s0 >= Scalar(0) )
        {
            if (s1 >= Scalar(0) )
            {
                return true;
            }
        }
    }
    // everything else is not valid
    return false;

}


bool fivePointAlgoCalibration(Vector3 spherePosition, const std::vector<cv::Point3d>& refDirections, const std::vector<cv::Point3d>& gazeDirections ,
    double* orientation , double* translation , bool fixTranslation = false ,
    Vector3 translationLowerBound = {15,5,5},Vector3 translationUpperBound = {15,5,5}
    )
{

    if( refDirections.size() < 5 || gazeDirections.size() < 5 )
        return false;

    cv::Mat refH , gazeH;
    cv::convertPointsFromHomogeneous(refDirections , refH);
    cv::convertPointsFromHomogeneous(gazeDirections , gazeH);

    auto essentialMatrix = cv::findEssentialMat( refH, gazeH );

    cv::Mat r1, r2, t;
    cv::decomposeEssentialMat( essentialMatrix , r1 , r2, t );

    std::cout << "r1: " << r1 << std::endl;
    std::cout << "r2: " << r2 << std::endl;
    std::cout << "t: " << t << std::endl;

    Eigen::Matrix<double, 3, 3> rotation1 ,rotation2;
    cv::cv2eigen(r1 , rotation1);
    cv::cv2eigen(r2 , rotation2);
    Eigen::Quaterniond q1 = Eigen::Quaterniond(rotation1);
    Eigen::Quaterniond q2 = Eigen::Quaterniond(rotation2);


    // find the solution which is physically plausible
    // by testing if the intersection point lie in front of the world camera
    Eigen::Matrix<double, 3, 1> origin = {0,0,0};
    Eigen::Matrix<double, 3, 1> translationEigen;
    cv::cv2eigen(t,translationEigen);

    bool valid = false;
    for(int i=0; i < refDirections.size(); i ++){

        auto refD = singleeyefitter::toEigen(refDirections.at(i));
        auto gazeD = singleeyefitter::toEigen(gazeDirections.at(i));

        refD.normalize();
        gazeD.normalize();
        auto gazeDTransformed = q1 * gazeD;

        Eigen::ParametrizedLine<double, 3> refLine = {origin, refD };
        Eigen::ParametrizedLine<double, 3> gazeLine = {translationEigen , gazeDTransformed};
        valid |= rayRayIntersectionInFront(refLine , gazeLine);
    }

    std::cout << "is valid: " <<  std::boolalpha << valid << std::endl;

    orientation[0] = q1.w();
    orientation[1] = q1.x();
    orientation[2] = q1.y();
    orientation[3] = q1.z();

    translation[0] = t.at<double>(0);
    translation[1] = t.at<double>(1);
    translation[2] = t.at<double>(2);

    //Ceres Matrices are RowMajor, where as Eigen is default ColumnMajor
    // Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation;
    // ceres::QuaternionToRotation( orientation , rotation.data() );
    // // ceres should always return a valid quaternion
    // // double det = r.determinant();
    // // std::cout << "det:: " << det << std::endl;
    // // if(  det == 1 ){
    // //     std::cout << "Error: No valid rotation matrix."   << std::endl;
    // //     return false;
    // // }


    // // we need to take the sphere position into account
    // // thus the actual translation is not right, because the local coordinate frame of the eye need to be translated in the opposite direction
    // // of the sphere coordinates

    // // since the actual translation is in world coordinates, the sphere translation needs to be calculated in world coordinates
    // Eigen::Matrix4d eyeToWorld =  Eigen::Matrix4d::Identity();
    // eyeToWorld.block<3,3>(0,0) = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(rotation.data());
    // eyeToWorld(0, 3) = translation[0];
    // eyeToWorld(1, 3) = translation[1];
    // eyeToWorld(2, 3) = translation[2];

    // Eigen::Vector4d sphereWorld = eyeToWorld * Eigen::Vector4d(spherePosition[0],spherePosition[1],spherePosition[2], 1.0 );
    // Vector3 sphereOffset =  sphereWorld.head<3>() - Vector3(translation);
    // Vector3 actualtranslation =  Vector3(translation) - sphereOffset;
    // // write the actual one back
    // translation[0] = actualtranslation[0];
    // translation[1] = actualtranslation[1];
    // translation[2] = actualtranslation[2];
    return true;

}

