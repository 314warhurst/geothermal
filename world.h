//
//  world.h
//  CubicSpline
//
//  Created by Henry Andrew Ewing Warhurst on 27/10/13.
//  Copyright (c) 2013 Henry A.E. Warhurst. All rights reserved.
//  henrywarhurst@gmail.com

#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include "Eigen/Dense"
#include "cubicspline.h"

#define PRESCISION (1e-15);


class World {
public:
    
    // Constructs a virtual subterranean environment composed of a
    // number of layers. Each layer boundary is defined by a cubic
    // spline fitting though the provided x and y control points.
    //
    // Currently assumes each boundary has the same number of
    // control pts.
    //
    // Input Parameters:
    //
    // MatrixXd& xPoints            - Each row contains the x
    //                                control points for a
    //                                particular boundary.
    // MatrixXd& yPoints            - Each row contains the y
    //                                control points for a
    //                                particular boundary.
    // densityMeans                 - The density of each layer.
    //                                Layer zero's density maps
    //                                to layerDensities[0].
    // densitySigmas                - The standard deviation to
    //                                use when randomising the
    //                                layer densities.
    World(const Eigen::MatrixXd& xPoints, const Eigen::MatrixXd& yPoints, const std::vector<double>& densityMeans, const std::vector<double>& densitySigmas);
    
    // Default constructor
    World();
    
    // Returns the layer in which the given point resides.
    // Layer indexing starts from zero and from the top. I.e.
    // the closest layer to the surface will always be layer zero.
    //
    // Input Parameters:
    //
    // double x                     - x coordinate of ctrl pt.
    // double y                     - y coordinate of ctrl pt.
    int getLayer(double x, double y) const;
    
    // Returns the layer as a density. Density values are being
    // generated about the values provided in the constructor
    // using a normal distrubution with sd = sigma.
    //
    // Input Parameters:
    //
    // int layer                    - The layer for which to
    //                                look up the density.
    // Returns:
    //
    // double density               - The density at the
    //                                provided layer. Returns
    //                                -1 if the density is not
    //                                found.
    double layerToDensity(int layer) const;
    
    // Setters & getters, usage is fairly self explanitory.
    void setXPoints(const Eigen::MatrixXd& xPoints);
    void setYPoints(const Eigen::MatrixXd& yPoints);
    Eigen::MatrixXd getXPoints(void) const;
    Eigen::MatrixXd getYPoints(void) const;
    void setNLayers(int nLayers);
    void setNControlPts(int nControlPts);
    
private:
    Eigen::MatrixXd _xPoints;
    Eigen::MatrixXd _yPoints;
    int _nLayers;
    int _nControlPts;
    std::vector<double> _densityMeans;
    std::vector<double> _densitySigmas;
    std::vector<std::vector<CubicPolynomial> > _splines;
    
    // Randomises density values from the constructor using
    // a normal distrubution.
    void randomiseDensities();
    
    // Constructs spline boundaries of subterranean world.
    // This function is implicitly called in the
    // constructor.
    void getSplines();
    
    // Used implicitly for checking whether a given point
    // resides above or below a given spline.
    bool isBelowCurve(const CubicPolynomial& cubic, double x, double y, double x_i) const;
    
    
};


class Random {
public:
    
    // Generates a random world object. Currently uses Boost's
    // random number generator with a normal distrubution.
    //
    // Input Parameters:
    //
    // int width                   - Horizontal width of the
    //                               worlds to be generated.
    // densityMeans                - The densities of each layer.
    // densitySigmas               - The standard deviation to use
    //                               when normally distributing
    //                               density values.
    // xPoints                     - x values of ctrl pts.
    // yPoints                     - y values of ctrl pts.
    // yPointsSigma                - Standard deviations for each
    //                               y point.
    static World randWorld(int width, const std::vector<double>& densityMeans, const std::vector<double>& densitySigmas, const Eigen::MatrixXd& xPoints, const Eigen::MatrixXd& yPoints, const Eigen::MatrixXd& yPointsSigma);
    
    // Prints the average world of the worlds given.
    //
    // Input Parameters:
    //
    // int nWorlds                  - The number of worlds to
    //                                be averaged across.
    // int height                   - The upper extent of the
    //                                world to be printed.
    // int width                    - The rightmost extent of the
    //                                world to be printed.
    // meanDensities                - The mean density of each layer.
    // densitySigmas                - The standard deviation to use
    //                                when randomly generating
    //                                normally distributed densities.
    // int widthRes                 - The horizontal resolution.
    // int heightRes                - The vertical resolution.
    // xPoints                     - x values of ctrl pts.
    // yPoints                     - y values of ctrl pts.
    // yPointsSigma                - Standard deviations for each
    //                               y point.
    static void printAverageWorld(int nWorlds, int height, int width, const std::vector<double>& densityMeans, const std::vector<double>& densitySigmas, int widthRes, int heightRes, const Eigen::MatrixXd& xPoints, const Eigen::MatrixXd& ypoints, const Eigen::MatrixXd& yPointsSigma);
    
    // Prints out the given world to the command line.
    //
    // Input Parameters:
    //
    // const World& world           - The world for printing.
    static void printWorld(const World& world);
    
    static void worldMatrixToCSV(const Eigen::MatrixXd& world, const std::string& path);
private:
    Eigen::MatrixXd baseYPpoints;
    
    //Generates random y points for the given set of x points.
    //Normally distributed around the pivot.
    static Eigen::VectorXd randYPoints(const Eigen::VectorXd& xPoints, const Eigen::VectorXd& yPoints, const Eigen::VectorXd& yPointsSigma);
};

