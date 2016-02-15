#include "common.h"
#include "math/Intersect.h"
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

template<typename Scalar>
struct Result{
    Scalar distanceSquared;
    Scalar rayLength0;
    Scalar rayLength1;
    bool valid = false;
};

// since for rayray distance both parameters s,t for eq r0=p0+s*d0 and r1=p1+t*d1 need to be positive
// ceres get's to know if the rays don't lie in the same direction with angle less than 90 degree
// Book: Geometric Tools for Computer Graphics, Side 413
template<typename Scalar, int Dim>
Result<Scalar> ceresRayRayDistanceSquared(const Eigen::ParametrizedLine<Scalar, Dim>& ray0, const Eigen::ParametrizedLine<Scalar, Dim>& ray1 )
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
                // Minimum at two  points of rays.
                Scalar det = Scalar(1) - a01 * a01;
                s0 /= det;
                s1 /= det;

                Vector closestPoint0 = ray0.origin() + s0 * ray0.direction();
                Vector closestPoint1 = ray1.origin() + s1 * ray1.direction();
                diff = closestPoint0 - closestPoint1;

                result.distanceSquared =  diff.dot(diff);
                result.rayLength0 = s0;
                result.rayLength1 = s1;
                result.valid = true;
                return result;
            }
        }
    }
    // everything else is not valid
    return result;

}

template<typename Scalar>
struct ResultPointRay{
    Scalar distanceSquared;
    bool valid = false;
};

template<typename Scalar >
ResultPointRay<Scalar> ceresRayPointDistanceSquared(const Eigen::ParametrizedLine<Scalar, 3>& ray, const Eigen::Matrix<Scalar, 3, 1>& point )
{

    typedef typename Eigen::ParametrizedLine<Scalar, 3>::VectorType Vector;

    ResultPointRay<Scalar> result;
    result.valid = false;

    Vector diff = point - ray.origin();
    Scalar rayParameter = ray.direction().dot(diff);
    if (rayParameter > Scalar(0) )
    {
        Vector raypoint = ray.origin() + rayParameter * ray.direction();
        diff = point - raypoint;
        result.distanceSquared = diff.dot(diff);
        result.valid = true;
    }

    return result;

}

struct TransformationRayRayError {
    TransformationRayRayError(const Vector3 refDirection,   const Vector3 gazeDirection )
        : refDirection(refDirection), gazeDirection(gazeDirection) {}

    template <typename T>
    bool operator()(
        const T* const orientation,  // orientation denoted by quaternionParameterization
        const T* const translation,  // followed by translation
        T* residuals) const
    {

        // Compute coordinates with current transformation matrix: y = Rx + t.
        Eigen::Matrix<T, 3, 1> gazeP = {T(gazeDirection[0]), T(gazeDirection[1]), T(gazeDirection[2])};
        Eigen::Matrix<T, 3, 1> refP = {T(refDirection[0]) , T(refDirection[1]) , T(refDirection[2])};
        Eigen::Matrix<T, 3, 1> t = {T(translation[0]) , T(translation[1]) , T(translation[2])};

        Eigen::Matrix<T, 3, 1> gazeTransformed;

        //rotate
        ceres::QuaternionRotatePoint(orientation, gazeP.data(), gazeTransformed.data() );
        //ceres::QuaternionRotatePoint(orientation, l2.data(), p2.data() );

        //translate

        Eigen::Matrix<T, 3, 1> origin = {T(0),T(0),T(0)};
        Eigen::ParametrizedLine<T, 3> gazeLine = {t , gazeTransformed};
        Eigen::ParametrizedLine<T, 3> refLine = {origin, refP };

        Result<T> result = ceresRayRayDistanceSquared(gazeLine , refLine);

        if(  result.valid ){
            residuals[0] = result.distanceSquared / (result.rayLength0 * result.rayLength1);
            return true;
        }
        return false;

    }

    const Vector3 gazeDirection;
    const Vector3 refDirection;
};

struct TransformationBinocularRayRayError {
    TransformationBinocularRayRayError(const Vector3 refDirection,   const Vector3 gazeDirection0,  const Vector3 gazeDirection1 )
        : refDirection(refDirection), gazeDirection0(gazeDirection0),gazeDirection1(gazeDirection1)  {}

    template <typename T>
    bool operator()(
        const T* const orientation0,  // orientation denoted by quaternionParameterization
        const T* const translation0,  // followed by translation
        const T* const orientation1,  // orientation denoted by quaternionParameterization
        const T* const translation1,  // followed by translation
        T* residuals) const
    {

        typedef Eigen::Matrix<T, 3, 1> Vector;
        // Compute coordinates with current transformation matrix: y = Rx + t.
        Vector gaze0 = {T(gazeDirection0[0]), T(gazeDirection0[1]), T(gazeDirection0[2])};
        Vector gaze1 = {T(gazeDirection1[0]), T(gazeDirection1[1]), T(gazeDirection1[2])};
        Vector refP = {T(refDirection[0]) , T(refDirection[1]) , T(refDirection[2])};
        Vector t0 = {T(translation0[0]) , T(translation0[1]) , T(translation0[2])};
        Vector t1 = {T(translation1[0]) , T(translation1[1]) , T(translation1[2])};

        Vector gazeTransformed0;
        Vector gazeTransformed1;

        //rotate
        ceres::QuaternionRotatePoint(orientation0, gaze0.data(), gazeTransformed0.data() );
        ceres::QuaternionRotatePoint(orientation1, gaze1.data(), gazeTransformed1.data() );

        //translate
        typedef Eigen::ParametrizedLine<T, 3> Line;
        Vector origin = {T(0),T(0),T(0)};
        Line gazeLine0 = {t0 , gazeTransformed0};
        Line gazeLine1 = {t1 , gazeTransformed1};
        Line refLine = {origin, refP };

        std::vector<Line> lines = {refLine, gazeLine0 ,gazeLine1};

        Vector intersectionPoint = singleeyefitter::nearest_intersect(lines);


        ResultPointRay<T> result0 = ceresRayPointDistanceSquared(gazeLine0 , intersectionPoint);
        ResultPointRay<T> result1 = ceresRayPointDistanceSquared(gazeLine1 , intersectionPoint);

        if(  result0.valid && result1.valid ){
            residuals[0] = result0.distanceSquared;
            residuals[1] = result1.distanceSquared;
            return true;
        }
        return false;

    }

    const Vector3 gazeDirection0;
    const Vector3 gazeDirection1;
    const Vector3 refDirection;
};



bool lineLineCalibration(Vector3 spherePosition, const std::vector<Vector3>& refDirections, const std::vector<Vector3>& gazeDirections ,
    double* orientation , double* translation , bool fixTranslation = false ,
    Vector3 translationLowerBound = {15,5,5},Vector3 translationUpperBound = {15,5,5}
    )
{

    // don't use Constructor 'Quaternion (const Scalar *data)' because the internal layout for coefficients is different from the one we use.
    // Memory Layout EIGEN: xyzw
    // Memory Layout CERES and the one we use: wxyz
    Eigen::Quaterniond q(orientation[0],orientation[1],orientation[2],orientation[3]);

    Problem problem;
    double epsilon = std::numeric_limits<double>::epsilon();

    for(int i=0; i<refDirections.size(); i++) {

        auto gaze = gazeDirections.at(i);
        auto ref = refDirections.at(i);
        gaze.normalize(); //just to be sure
        ref.normalize(); //just to be sure
        i++;

        // do a check to handle parameters we can't solve
        // First: the length of the directions must not be zero
        // Second: the angle between gaze direction and reference direction must not be greater 90 degrees, considering the initial orientation

        bool valid = true;
        valid |= gaze.norm() >= epsilon;
        valid |= ref.norm() >= epsilon;
        valid |= (q*gaze).dot(ref) >= epsilon;

        if( valid ){
            CostFunction* cost = new AutoDiffCostFunction<TransformationRayRayError , 1, 4, 3 >(new TransformationRayRayError(ref, gaze ));
            // TODO use a loss function, to handle gaze point outliers
            problem.AddResidualBlock(cost, nullptr, orientation,  translation);
        }else{
            std::cout << "no valid direction vector"  << std::endl;
        }
    }

    if( problem.NumResidualBlocks() == 0 ){
        std::cout << "nothing to solve"  << std::endl;
        return false;
    }

    if (fixTranslation)
    {
        problem.SetParameterBlockConstant(translation);
    }else{

        Vector3 upperBound = Vector3(translation) + translationUpperBound;
        Vector3 lowerBound = Vector3(translation) - translationLowerBound;

        problem.SetParameterLowerBound(translation, 0 , lowerBound[0] );
        problem.SetParameterLowerBound(translation, 1 , lowerBound[1] );
        problem.SetParameterLowerBound(translation, 2 , lowerBound[2] );

        problem.SetParameterUpperBound(translation, 0 , upperBound[0] );
        problem.SetParameterUpperBound(translation, 1 , upperBound[1] );
        problem.SetParameterUpperBound(translation, 2 , upperBound[2] );
    }



    ceres::LocalParameterization* quaternionParameterization = new ceres::QuaternionParameterization; // owned by the problem
    problem.SetParameterization(orientation, quaternionParameterization);

    // Build and solve the problem.
    Solver::Options options;
    options.max_num_iterations = 1000;
    options.linear_solver_type = ceres::DENSE_QR;
    //options.parameter_tolerance = 1e-14;
    options.function_tolerance = 1e-10;
    options.gradient_tolerance = 1e-20;
    options.minimizer_progress_to_stdout = false;
    options.logging_type = ceres::SILENT;

    options.check_gradients = true;
    Solver::Summary summary;

    Solve(options, &problem, &summary);

    // std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";

    if( summary.termination_type != ceres::TerminationType::CONVERGENCE  ){
        std::cout << "Termination Error: " << ceres::TerminationTypeToString(summary.termination_type) << std::endl;
        return false;
    }

    //Ceres Matrices are RowMajor, where as Eigen is default ColumnMajor
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation;
    ceres::QuaternionToRotation( orientation , rotation.data() );
    // ceres should always return a valid quaternion
    // double det = r.determinant();
    // std::cout << "det:: " << det << std::endl;
    // if(  det == 1 ){
    //     std::cout << "Error: No valid rotation matrix."   << std::endl;
    //     return false;
    // }


    // we need to take the sphere position into account
    // thus the actual translation is not right, because the local coordinate frame of the eye need to be translated in the opposite direction
    // of the sphere coordinates

    // since the actual translation is in world coordinates, the sphere translation needs to be calculated in world coordinates
    Eigen::Matrix4d eyeToWorld =  Eigen::Matrix4d::Identity();
    eyeToWorld.block<3,3>(0,0) = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(rotation.data());
    eyeToWorld(0, 3) = translation[0];
    eyeToWorld(1, 3) = translation[1];
    eyeToWorld(2, 3) = translation[2];

    Eigen::Vector4d sphereWorld = eyeToWorld * Eigen::Vector4d(spherePosition[0],spherePosition[1],spherePosition[2], 1.0 );
    Vector3 sphereOffset =  sphereWorld.head<3>() - Vector3(translation);
    Vector3 actualtranslation =  Vector3(translation) - sphereOffset;
    // write the actual one back
    translation[0] = actualtranslation[0];
    translation[1] = actualtranslation[1];
    translation[2] = actualtranslation[2];
    return true;

}

bool lineLineCalibrationBinocular( Vector3 spherePosition0, Vector3 spherePosition1, const std::vector<Vector3>& refDirections, const std::vector<Vector3>& gazeDirections0 ,const std::vector<Vector3>& gazeDirections1,
    double* orientation0 ,double* orientation1 , double* translation0, double* translation1, bool fixTranslation = false ,
    Vector3 translationLowerBound = {15,5,5},Vector3 translationUpperBound = {15,5,5}
    )
{

    // don't use Constructor 'Quaternion (const Scalar *data)' because the internal layout for coefficients is different from the one we use.
    // Memory Layout EIGEN: xyzw
    // Memory Layout CERES and the one we use: wxyz
    Eigen::Quaterniond q0(orientation0[0],orientation0[1],orientation0[2],orientation0[3]);
    Eigen::Quaterniond q1(orientation1[0],orientation1[1],orientation1[2],orientation1[3]);

    Problem problem;
    double epsilon = std::numeric_limits<double>::epsilon();

    for(int i=0; i<refDirections.size(); i++) {

        auto gaze0 = gazeDirections0.at(i);
        auto gaze1 = gazeDirections1.at(i);
        auto ref = refDirections.at(i);
        gaze0.normalize(); //just to be sure
        gaze1.normalize(); //just to be sure
        ref.normalize(); //just to be sure
        i++;

        // do a check to handle parameters we can't solve
        // First: the length of the directions must not be zero
        // Second: the angle between gaze direction and reference direction must not be greater 90 degrees, considering the initial orientation

        bool valid = true;
        valid |= gaze0.norm() >= epsilon;
        valid |= gaze1.norm() >= epsilon;
        valid |= ref.norm() >= epsilon;
        valid |= (q0*gaze0).dot(ref) >= epsilon;
        valid |= (q1*gaze1).dot(ref) >= epsilon;

        if( valid ){
            CostFunction* cost = new AutoDiffCostFunction<TransformationBinocularRayRayError , 2, 4, 3, 4, 3 >(new TransformationBinocularRayRayError(ref, gaze0 , gaze1 ));
            // TODO use a loss function, to handle gaze point outliers
            problem.AddResidualBlock(cost, nullptr, orientation0,  translation0, orientation1,  translation1);
        }else{
            std::cout << "no valid direction vector"  << std::endl;
        }
    }

    if( problem.NumResidualBlocks() == 0 ){
        std::cout << "nothing to solve"  << std::endl;
        return false;
    }

    if (fixTranslation)
    {
        problem.SetParameterBlockConstant(translation0);
        problem.SetParameterBlockConstant(translation1);
    }else{

        Vector3 upperBound0 = Vector3(translation0) + translationUpperBound;
        Vector3 lowerBound0 = Vector3(translation0) - translationLowerBound;
        problem.SetParameterLowerBound(translation0, 0 , lowerBound0[0] );
        problem.SetParameterLowerBound(translation0, 1 , lowerBound0[1] );
        problem.SetParameterLowerBound(translation0, 2 , lowerBound0[2] );
        problem.SetParameterUpperBound(translation0, 0 , upperBound0[0] );
        problem.SetParameterUpperBound(translation0, 1 , upperBound0[1] );
        problem.SetParameterUpperBound(translation0, 2 , upperBound0[2] );

        Vector3 upperBound1 = Vector3(translation1) + translationUpperBound;
        Vector3 lowerBound1 = Vector3(translation1) - translationLowerBound;
        problem.SetParameterLowerBound(translation1, 0 , lowerBound1[0] );
        problem.SetParameterLowerBound(translation1, 1 , lowerBound1[1] );
        problem.SetParameterLowerBound(translation1, 2 , lowerBound1[2] );
        problem.SetParameterUpperBound(translation1, 0 , upperBound1[0] );
        problem.SetParameterUpperBound(translation1, 1 , upperBound1[1] );
        problem.SetParameterUpperBound(translation1, 2 , upperBound1[2] );
    }



    ceres::LocalParameterization* quaternionParameterization = new ceres::QuaternionParameterization; // owned by the problem
    problem.SetParameterization(orientation0, quaternionParameterization);
    problem.SetParameterization(orientation1, quaternionParameterization);

    // Build and solve the problem.
    Solver::Options options;
    options.max_num_iterations = 1000;
    options.linear_solver_type = ceres::DENSE_QR;
    //options.parameter_tolerance = 1e-14;
    options.function_tolerance = 1e-10;
    options.gradient_tolerance = 1e-20;
    options.minimizer_progress_to_stdout = false;
    options.logging_type = ceres::SILENT;

    options.check_gradients = true;
    Solver::Summary summary;

    Solve(options, &problem, &summary);

    // std::cout << summary.BriefReport() << "\n";
    std::cout << summary.FullReport() << "\n";

    if( summary.termination_type != ceres::TerminationType::CONVERGENCE  ){
        std::cout << "Termination Error: " << ceres::TerminationTypeToString(summary.termination_type) << std::endl;
        return false;
    }

    //Ceres Matrices are RowMajor, where as Eigen is default ColumnMajor
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation0;
    ceres::QuaternionToRotation( orientation0 , rotation0.data() );
    Eigen::Matrix<double, 3, 3, Eigen::RowMajor> rotation1;
    ceres::QuaternionToRotation( orientation1 , rotation1.data() );
    // ceres should always return a valid quaternion
    // double det = r.determinant();
    // std::cout << "det:: " << det << std::endl;
    // if(  det == 1 ){
    //     std::cout << "Error: No valid rotation matrix."   << std::endl;
    //     return false;
    // }


    // we need to take the sphere position into account
    // thus the actual translation is not right, because the local coordinate frame of the eye need to be translated in the opposite direction
    // of the sphere coordinates

    // since the actual translation is in world coordinates, the sphere translation needs to be calculated in world coordinates
    Eigen::Matrix4d eyeToWorld0 =  Eigen::Matrix4d::Identity();
    eyeToWorld0.block<3,3>(0,0) = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(rotation0.data());
    eyeToWorld0(0, 3) = translation0[0];
    eyeToWorld0(1, 3) = translation0[1];
    eyeToWorld0(2, 3) = translation0[2];

    Eigen::Vector4d sphereWorld0 = eyeToWorld0 * Eigen::Vector4d(spherePosition0[0],spherePosition0[1],spherePosition0[2], 1.0 );
    Vector3 sphereOffset0 =  sphereWorld0.head<3>() - Vector3(translation0);
    Vector3 actualtranslation0 =  Vector3(translation0) - sphereOffset0;
    // write the actual one back
    translation0[0] = actualtranslation0[0];
    translation0[1] = actualtranslation0[1];
    translation0[2] = actualtranslation0[2];

    Eigen::Matrix4d eyeToWorld1 =  Eigen::Matrix4d::Identity();
    eyeToWorld1.block<3,3>(0,0) = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(rotation1.data());
    eyeToWorld1(0, 3) = translation1[0];
    eyeToWorld1(1, 3) = translation1[1];
    eyeToWorld1(2, 3) = translation1[2];

    Eigen::Vector4d sphereWorld1 = eyeToWorld1 * Eigen::Vector4d(spherePosition1[0],spherePosition1[1],spherePosition1[2], 1.0 );
    Vector3 sphereOffset1 =  sphereWorld1.head<3>() - Vector3(translation1);
    Vector3 actualtranslation1 =  Vector3(translation1) - sphereOffset1;
    // write the actual one back
    translation1[0] = actualtranslation1[0];
    translation1[1] = actualtranslation1[1];
    translation1[2] = actualtranslation1[2];

    return true;

}

