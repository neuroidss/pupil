


#include "common.h"
#include <vector>
#include <cstdio>
#include <limits>

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <Eigen/Geometry>
#include "ceres/CeresParametrization.h"

using ceres::AutoDiffCostFunction;
using ceres::NumericDiffCostFunction;
using ceres::CauchyLoss;
using ceres::CostFunction;
using ceres::LossFunction;
using ceres::Problem;
using ceres::Solve;
using ceres::Solver;


struct CoplanarityError {
    CoplanarityError(const Vector3 refDirection,   const Vector3 gazeDirection )
        : refDirection(refDirection), gazeDirection(gazeDirection) {}

    template <typename T>
    bool operator()(
        const T* const orientation,  // orientation denoted by quaternion Parameterization
        const T* const translation,  // followed by translation
        T* residuals) const
    {

        Eigen::Matrix<T, 3, 1> gazeD = {T(gazeDirection[0]), T(gazeDirection[1]), T(gazeDirection[2])};
        Eigen::Matrix<T, 3, 1> refD = {T(refDirection[0]) , T(refDirection[1]) , T(refDirection[2])};

        //Ceres Matrices are RowMajor, where as Eigen is default ColumnMajor
        Eigen::Matrix<T, 3, 3, Eigen::RowMajor> rotationMatrix;
        ceres::QuaternionToRotation( orientation , rotationMatrix.data() );
         // cross-product matrix of the translation
        Eigen::Matrix<T, 3, 3 > translationMatrix;
        translationMatrix << T(0) , T(-translation[2]) , T(translation[1]),
                             T(translation[2]), T(0), T(-translation[0]),
                             T(-translation[1]), T(translation[0]), T(0);

        std::cout << "t: " << translation[0]<<translation[1]<<translation[2] << std::endl;
        Eigen::Matrix<T, 3, 3 > essentialMatrix = translationMatrix * rotationMatrix.transpose();

        //TODO add weighting factors to the residual , better approximation
        //coplanarity constraint  x1.T * E * x2 = 0
        auto res = refD.transpose() * essentialMatrix * gazeD;
        std::cout << "res: " << res[0] << std::endl;
        residuals[0] = res[0];
        return true;


    }

    const Vector3 gazeDirection;
    const Vector3 refDirection;
};



bool lineLineCalibration(Vector3 spherePosition, const std::vector<Vector3>& refDirections, const std::vector<Vector3>& gazeDirections ,
    double *const orientation , double *const translation , bool fixTranslation = false ,
    Vector3 translationLowerBound = {15,5,5},Vector3 translationUpperBound = {15,5,5}
    )
{

    // don't use Constructor 'Quaternion (const Scalar *data)' because the internal layout for coefficients is different from the one we use.
    // Memory Layout EIGEN: xyzw
    // Memory Layout CERES and the one we use: wxyz
    Eigen::Quaterniond q(orientation[0],orientation[1],orientation[2],orientation[3]); // don't mapp orientation
    Vector3 t =  Eigen::Map<Eigen::Matrix<double,3,1> >(translation);
    t.normalize();
    translation[0] = t[0];
    translation[1] = t[1];
    translation[2] = t[2];

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

            CostFunction* cost = new AutoDiffCostFunction<CoplanarityError , 1, 4, 3 >(new CoplanarityError(ref, gaze ));
            // TODO use a loss function, to handle gaze point outliers
            problem.AddResidualBlock(cost, nullptr, orientation,  translation );
        }else{
            std::cout << "no valid direction vector"  << std::endl;
        }
    }

    if( problem.NumResidualBlocks() == 0 ){
        std::cout << "nothing to solve"  << std::endl;
        return false;
    }

    // if (fixTranslation)
    // {
    //     problem.SetParameterBlockConstant(translation);
    // }
    problem.SetParameterBlockConstant(orientation);

    //else{

    //     Vector3 upperBound = Vector3(translation) + translationUpperBound;
    //     Vector3 lowerBound = Vector3(translation) - translationLowerBound;

    //     problem.SetParameterLowerBound(translation, 0 , lowerBound[0] );
    //     problem.SetParameterLowerBound(translation, 1 , lowerBound[1] );
    //     problem.SetParameterLowerBound(translation, 2 , lowerBound[2] );

    //     problem.SetParameterUpperBound(translation, 0 , upperBound[0] );
    //     problem.SetParameterUpperBound(translation, 1 , upperBound[1] );
    //     problem.SetParameterUpperBound(translation, 2 , upperBound[2] );
    // }



    ceres::LocalParameterization* quaternionParameterization = new ceres::QuaternionParameterization; // owned by the problem
    problem.SetParameterization(orientation, quaternionParameterization);

    ceres::LocalParameterization* normedTranslationParameterization = new pupillabs::Fixed3DNormParametrization(1.0); // owned by the problem
    problem.SetParameterization(translation, normedTranslationParameterization);

    // Build and solve the problem.
    Solver::Options options;
    options.max_num_iterations = 1000;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.parameter_tolerance = 1e-14;
    options.function_tolerance = 1e-10;
    options.gradient_tolerance = 1e-20;
    options.minimizer_progress_to_stdout = true;
    //options.logging_type = ceres::SILENT;

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
    double translationFactor = 30.0;

    // since the actual translation is in world coordinates, the sphere translation needs to be calculated in world coordinates
    Eigen::Matrix4d eyeToWorld =  Eigen::Matrix4d::Identity();
    eyeToWorld.block<3,3>(0,0) = Eigen::Map<Eigen::Matrix<double,3,3,Eigen::RowMajor> >(rotation.data());
    eyeToWorld(0, 3) = translation[0]*translationFactor;
    eyeToWorld(1, 3) = translation[1]*translationFactor;
    eyeToWorld(2, 3) = translation[2]*translationFactor;

    Eigen::Vector4d sphereWorld = eyeToWorld * Eigen::Vector4d(spherePosition[0],spherePosition[1],spherePosition[2], 1.0 );
    Vector3 sphereOffset =  sphereWorld.head<3>() - Vector3(translation);
    Vector3 actualtranslation =  Vector3(translation) - sphereOffset;
    // write the actual one back
    translation[0] = actualtranslation[0];
    translation[1] = actualtranslation[1];
    translation[2] = actualtranslation[2];
    return true;

}
