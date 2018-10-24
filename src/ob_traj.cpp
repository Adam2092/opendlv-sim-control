#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "data_structure.hpp"

vector<Obstacle> ob_traj(double t, bool isDynamic, Global_variables gl)
{
    vector<Obstacle> traj_ob;
    for (int i = 0; i < gl.no_ob; i++)
    {
        double v = isDynamic ? std::sqrt(t) : 0.0;
        Obstacle curr;
        curr.pos_x = gl.pos_ob_array_pre[i](0) + v * t;
        curr.pos_y = gl.pos_ob_array_pre[i](1);
        curr.vel_x = v;
        curr.vel_y = 0.0;
        curr.acc_x = 0.0;
        curr.acc_y = 0.0;
        curr.radius = gl.radius_pre[i];

        traj_ob.push_back(curr);
    }
    return traj_ob;
}