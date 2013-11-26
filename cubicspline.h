//
//  cubicspline.h
//  CubicSpline
//
//  Created by Henry Andrew Ewing Warhurst on 13/10/13.
//  Copyright (c) 2013 Henry Andrew Ewing Warhurst. All rights reserved.
//  Email: henrywarhurst@gmail.com

#pragma once

#include <iostream>
#include <vector>
#include "Eigen/Dense"

class CubicPolynomial {
public:
    // Default constructor
    CubicPolynomial();
    
    // Constructs a cubic polynomial with the provided constants.
    // Please note that we are using the form:
    //
    //      a(x - x_i)^3 + b(x - x_i)^2 + c(x - x_i) + d
    //
    // Where x_i is the x value at the i-th pivot.
    //
    // Input Parameters:
    //
    // double a                     - x^3 coefficient
    // double b                     - x^2 coefficient
    // double c                     - x^1 coefficient
    // double d                     - x^0 coefficient
    CubicPolynomial(const double a, const double b, const double c, const double d);
    
    // Getters & setters
    double getA() const;
    double getB() const;
    double getC() const;
    double getD() const;
    void setA(double a);
    void setB(double b);
    void setC(double c);
    void setD(double d);
    
private:
    double _a;
    double _b;
    double _c;
    double _d;
};


class CubicSpline {
public:
    
    // Constructs a cubic spline running through the supplied
    // x and y values. Note: if the size of x and y do not match
    // this will crash.
    //
    // Input Parameters:
    //
    // VectorXd& xpoints                - X values of the knots.
    // VectorXd& ypoints                - Y values of the knots.
    //
    // Returns:
    //
    // vector<CubicPolynomial> spline   - Vector of spline
    //                                    coefficients.
    static std::vector<CubicPolynomial> cubicSpline(const Eigen::VectorXd& xpoints, const Eigen::VectorXd& ypoints);
    
    // Prints spline coefficients to a CSV.
    //
    // Input Parameters:
    //
    // vector<CubicPolynomial>& inputSplines - The spline coeffs
    //                                          to be printed.
    // string& path                          - The place to put
    //                                          the CSV.
    static void splinesToCSV(const std::vector<CubicPolynomial>& inputSplines, const std::string& path);
};

