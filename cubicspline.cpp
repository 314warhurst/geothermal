//
//  cubicspline.cpp
//  CubicSpline
//
//  Created by Henry Andrew Ewing Warhurst on 13/10/13.
//  Copyright (c) 2013 Henry Andrew Ewing Warhurst. All rights reserved.
//

#include <fstream>
#include <Eigen/Dense>
#include "cubicspline.h"

using namespace std;


vector<CubicPolynomial> CubicSpline::cubicSpline(const Eigen::VectorXd& xpoints, const Eigen::VectorXd& ypoints)
{
    //Crash if the sizes don't match up
    assert(xpoints.size() == ypoints.size());
    
    size_t n = xpoints.size();
    vector<CubicPolynomial> output(n - 1);
    vector<double> h(n - 1, 0);
    vector<double> mainDiag(n, 0);
    vector<double> upperDiag(n - 1, 0);
    vector<double> lowerDiag(n - 1, 0);
    Eigen::MatrixXd leftMatrix = Eigen::MatrixXd::Zero(n, n);
    Eigen::VectorXd ansVec = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd secondDerivs = Eigen::VectorXd::Zero(n);
    
    // Construct output
    for (int i = 0; i < n - 1; i++) {
        CubicPolynomial cur;
        output[i] = cur;
    }
    
    //Get h vals
    for (int i = 0; i < n - 1; i++) {
        h[i] = xpoints(i + 1) - xpoints(i);
    }
    
    //Main diagonal and right side column vector
    for (int i = 1; i < n - 1; i++) {
        mainDiag[i] = 2 * (h[i - 1] + h[i]);
        ansVec[i] = 6 * ((ypoints[i + 1] - ypoints[i])/h[i] - (ypoints[i] - ypoints[i - 1])/h[i - 1]);
    }
    mainDiag[0] = 1;
    mainDiag[n - 1] = 1;
    
    //Secondary Diagonals
    for (int i = 0; i < n - 1; i++) {
        upperDiag[i] = h[i];
        lowerDiag[i] = h[i];
    }
    upperDiag[0] = 0;
    lowerDiag[n - 2] = 0;
    
    //Fill left matrix
    for (int i = 0; i < n; i++) {
        leftMatrix(i, i) = mainDiag[i];
        if (i == 0) continue;
        leftMatrix(i - 1, i) = upperDiag[i - 1];
        leftMatrix(i, i - 1) = lowerDiag[i - 1];
    }
    
    //Solve linear system
    
    secondDerivs = leftMatrix.fullPivHouseholderQr().solve(ansVec);
    
    // Fill a, b, c, and d
    for (int i = 0; i < n - 1; i++) {
        output[i].setA((secondDerivs(i + 1) - secondDerivs(i))/(6.0 * h[i]));
        output[i].setB(secondDerivs(i)/2.0);
        output[i].setC((ypoints[i + 1] - ypoints[i])/h[i] - (secondDerivs(i + 1) * h[i])/6 - (secondDerivs(i) * h[i])/3.0);
        output[i].setD(ypoints[i]);
    }
    return output;
}

void CubicSpline::splinesToCSV(const vector<CubicPolynomial>& inputSplines, const string &path)
{
    
    ofstream file(path.c_str());
    if (file.is_open()) {
        for (int i = 0; i < inputSplines.size(); i++) {
            file << inputSplines[i].getA() << ",";
            file << inputSplines[i].getB() << ",";
            file << inputSplines[i].getC() << ",";
            file << inputSplines[i].getD() << endl;
        }
    }
    file.close();
}

CubicPolynomial::CubicPolynomial()
:   _a(0.0),
    _b(0.0),
    _c(0.0),
    _d(0.0)
{}

CubicPolynomial::CubicPolynomial(double a, double b, double c, double d)
:   _a(a),
    _b(b),
    _c(c),
    _d(d)
{}

double CubicPolynomial::getA() const
{
    return _a;
}

double CubicPolynomial::getB() const
{
    return _b;
}

double CubicPolynomial::getC() const
{
    return _c;
}

double CubicPolynomial::getD() const
{
    return _d;
}

void CubicPolynomial::setA(double a)
{
    _a = a;
}

void CubicPolynomial::setB(double b)
{
    _b = b;
}

void CubicPolynomial::setC(double c)
{
    _c = c;
}

void CubicPolynomial::setD(double d)
{
    _d = d;
}