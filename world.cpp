//
//  world.cpp
//  CubicSpline
//
//  Created by Henry Andrew Ewing Warhurst on 27/10/13.
//  Copyright (c) 2013 Henry Andrew Ewing Warhurst. All rights reserved.
//

#include <cmath>
#include <sstream>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include "world.h"


World::World(const Eigen::MatrixXd& xPoints, const Eigen::MatrixXd& yPoints, const std::vector<double>& densityMeans, const std::vector<double>& densitySigmas)
:   _xPoints(xPoints),
    _yPoints(yPoints),
    _nLayers((int) xPoints.col(0).size() + 1),
    _nControlPts((int) xPoints.row(0).size()),
    _densityMeans(densityMeans),
    _densitySigmas(densitySigmas)
{
    _splines.resize(_nLayers - 1);
    //Populate spline coefficient vectors
    getSplines();
    // Randomise density values using a normal distrobution
    randomiseDensities();
}

World::World()
{}

int World::getLayer(double x, double y) const
{
    // Find layer
    int rangeIdx = -1;
    int layer = 0;
    for (int i = 0; i < _nLayers - 1; i++) {
        for (int j = 0; j < _nControlPts - 1; j++) {
            if (_xPoints(i, j) <= x && _xPoints(i, j + 1) >= x) {
                rangeIdx = j;
                break;
            }
        }
        
        if (isBelowCurve(_splines[i][rangeIdx], x, y, _xPoints(i, rangeIdx))) {
            layer++;
        } else {
            return layer;
        }
    }
    return layer;
}

bool World::isBelowCurve(const CubicPolynomial& cubic, double x, double y, double x_i) const
{
    
    double expValue = cubic.getA()*pow((x - x_i),3) + cubic.getB()*pow((x - x_i), 2) + cubic.getC()*(x - x_i) + cubic.getD();
    
    //Check if below or above
    if (expValue > y)
        return true;
    else
        return false;
}

double World::layerToDensity(int layer) const
{
    //Check if layer is in the correct range
    assert(layer >= 0 && layer < _densityMeans.size());
    
    for (int i = 0; i < _densityMeans.size(); i++) {
        if (layer == i)
            return _densityMeans[i];
    }
    
    //If we reach here, it's not found for some reason. Return -1 to notify user.
    return -1.0;
}

void World::getSplines()
{
    CubicSpline cs;
    for (int i = 0; i < _nLayers - 1; i++) {
        
        std::vector<CubicPolynomial> curSpline;
        curSpline = cs.cubicSpline(_xPoints.row(i), _yPoints.row(i));
        _splines[i] = curSpline;
    }
}

void World::randomiseDensities()
{
    static boost::mt19937 generator(static_cast<unsigned int>(std::time(0)));
    
    
    for (int i = 0; i < _densityMeans.size(); i++) {
        // Generate random number
        boost::normal_distribution<> normalDist(0.0, _densitySigmas[i]);
        boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > numberGenerator(generator, normalDist);
        _densityMeans[i] += numberGenerator();
    }
}

void World::setXPoints(const Eigen::MatrixXd &xPoints)
{
    _xPoints = xPoints;
    //Recalculate parameters
    _nLayers = (int) xPoints.col(0).size() + 1;
    _nControlPts = (int) xPoints.row(0).size();
    //Regenerate the world
    getSplines();
}

void World::setYPoints(const Eigen::MatrixXd &yPoints)
{
    _yPoints = yPoints;
    //Regenerate the world
    getSplines();
}

Eigen::MatrixXd World::getXPoints() const
{
    return _xPoints;
}

Eigen::MatrixXd World::getYPoints() const
{
    return _yPoints;
}

void World::setNLayers(int nLayers)
{
    _nLayers = nLayers;
}

void World::setNControlPts(int nControlPts)
{
    _nControlPts = nControlPts;
}

Eigen::VectorXd Random::randYPoints(const Eigen::VectorXd& xPoints, const Eigen::VectorXd& yPoints, const Eigen::VectorXd& yPointsSigma)
{
    int nPoints = (int) xPoints.size();
    Eigen::VectorXd newYPoints(yPoints.rows(), yPoints.cols());
    
    static boost::mt19937 generator(static_cast<unsigned int>(std::time(0)));
    boost::normal_distribution<> normalDist(0.0, 1.0);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > numberGenerator(generator, normalDist);
    
    for (int i = 0; i < nPoints; ++i) {
        double temp = numberGenerator();
        newYPoints(i) = temp + yPoints(i);
    }
    return newYPoints;
}

World Random::randWorld(int width, const std::vector<double>& densityMeans, const std::vector<double>& densitySigmas, const Eigen::MatrixXd& xPoints, const Eigen::MatrixXd& yPoints, const Eigen::MatrixXd& yPointsSigma)
{
    Eigen::MatrixXd tempYPoints = yPoints;
    
    for (int i = 0; i < yPoints.rows(); i++) {
        tempYPoints.row(i) = randYPoints(xPoints.row(i), yPoints.row(i), yPointsSigma.row(i));
    }

    World curWorld(xPoints, tempYPoints, densityMeans, densitySigmas);
    return curWorld;
}

void Random::printAverageWorld(int nWorlds, int height, int width, const std::vector<double>& densityMeans, const std::vector<double>& densitySigmas, int heightRes, int widthRes, const Eigen::MatrixXd& xPoints, const Eigen::MatrixXd& yPoints, const Eigen::MatrixXd& yPointsSigma)
{
    Eigen::MatrixXd averageWorld = Eigen::MatrixXd::Zero(heightRes, widthRes);
    World curWorld;
    
    Random rand;
    std::string filename = "/Users/hwar/Documents/";
    std::ostringstream stream;
    
    // Find average layer for each cell
    for (int k = 0; k < nWorlds; k++) {
        curWorld = randWorld(width, densityMeans, densitySigmas, xPoints, yPoints, yPointsSigma);
        //printWorld(curWorld);
        
        for (int i = 0; i < widthRes; i++) {
            for (int j = 0; j < heightRes; j++) {
                float x = i * (width - 1) / (float) widthRes;
                float y = j * (height - 1) / (float) heightRes;
                
                double layer = (double) curWorld.getLayer(x, y);
                double curDensityMean = curWorld.layerToDensity(layer)/ (double) nWorlds;
                averageWorld(i, j) += curDensityMean;
            }
        }
    }
    
    rand.worldMatrixToCSV(averageWorld, "/Users/hwar/Documents/earth_final.csv");
}

void Random::printWorld(const World& world)
{
    for (int i = 10; i >= 0; --i) {
        for (int j = 0; j < world.getXPoints().row(0).size(); ++j) {
            double layer = (double) world.getLayer(j, i);
            double density = world.layerToDensity(layer);
            std::printf("%1.1f\t", density);
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void Random::worldMatrixToCSV(const Eigen::MatrixXd& world, const std::string& path)
{
    std::ofstream file(path);
    for (int i = (int) world.cols() - 1; i >= 0; i--) {
        for (int j = 0; j < world.rows(); j++) {
            file << world(j,i);
            // Put a comma except on the last column
            if (j != world.row(0).size())
                file << ",";
        }
        file << std::endl;
    }
}



