//
//  main.cpp
//  CubicSpline
//
//  Created by Henry Andrew Ewing Warhurst on 25/09/13.
//  Copyright (c) 2013 Henry Andrew Ewing Warhurst. All rights reserved.
//

#include <iostream>
#include <vector>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "Eigen/Dense"
#include "cubicspline.h"
#include "world.h"
#include "constants.h"
using namespace std;

int main()
{
    
    Eigen::MatrixXd xPoints(3, WORLD_WIDTH);
    Eigen::MatrixXd yPoints(3, WORLD_WIDTH);
    Eigen::MatrixXd yPointsSigma(3, WORLD_WIDTH);
    
    std::vector<double> pivotMeans(3);
    pivotMeans[0] = 8;
    pivotMeans[1] = 5;
    pivotMeans[2] = 2;
    
    // Randomly generate pivot points
    //static boost::mt19937 generator(static_cast<unsigned int>(std::time(0)));
    static boost::mt19937 generator(static_cast<unsigned int>(24));
    boost::normal_distribution<> normalDist(0.0, MEAN_SIGMA);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > numberGenerator(generator, normalDist);
    
    for (int i = 0; i < yPoints.rows(); i++) {
        for (int j = 0; j < yPoints.cols(); j++) {
            yPoints(i,j) = pivotMeans[i] + numberGenerator();
        }
    }
    
    //Set x points, same xPoints for all worlds
    for (int i = 0; i < xPoints.cols(); ++i) {
        xPoints(0, i) = i;
        xPoints(1, i) = i;
        xPoints(2, i) = i;
    }
    
    //Randomly generate the standard deviations of the y points
    boost::normal_distribution<> normalDist2(SIGMA_MEAN, SIGMA_SIGMA);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > numberGenerator2(generator, normalDist2);
    
    for (int i = 0; i < yPointsSigma.rows(); i++) {
        for (int j = 0; j < yPointsSigma.cols(); j++) {
            yPointsSigma(i,j) = numberGenerator2();
        }
    }
    
    vector<double> densityMeans(4);
    vector<double> densitySigmas(4);
    densityMeans[0] = 1.0;
    densityMeans[1] = 5.0;
    densityMeans[2] = 10.0;
    densityMeans[3] = 19.0;
    
    densitySigmas[0] = 0.5;
    densitySigmas[1] = 0.75;
    densitySigmas[2] = 1.0;
    densitySigmas[3] = 0.25;
    Random rand;
    rand.printAverageWorld(N_ITERATIONS, WORLD_HEIGHT, WORLD_WIDTH, densityMeans, densitySigmas, HEIGHT_RESOLUTION, WIDTH_RESOLUTION, xPoints, yPoints, yPointsSigma);
}
