//
//  constants.h
//  CubicSpline
//
//  Created by Henry Andrew Ewing Warhurst on 25/11/2013.
//  Copyright (c) 2013 Henry Andrew Ewing Warhurst. All rights reserved.
//

#pragma once

#include <string>
#include <iostream>

// SEE MAIN FOR THE SETTING OF DENSITY PARAMETERS

const int N_ITERATIONS = 10000;
const int WORLD_HEIGHT = 10;
const int WORLD_WIDTH = 10;
const int HEIGHT_RESOLUTION = 600;
const int WIDTH_RESOLUTION = 600;
// Sd to use when generating the means of the ctrl pts.
const double MEAN_SIGMA = 1.0;
// Mean to use when generating the standard deviations of
// the ctrl pts.
const double SIGMA_MEAN = 1.0;
// Sd to use when generating the standard deviations of the
// ctrl pts.
const double SIGMA_SIGMA = 0.1;
std::string CSV_PATH = "/Users/hwar/Documents/earth_final.csv";

