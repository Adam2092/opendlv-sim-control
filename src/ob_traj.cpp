#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "data_structure.hpp"

/*class Obstacle
{
public:
    double pos_x{0.0};
    double pos_y{0.0};
    double vel_x{0.0};
    double vel_y{0.0};
    double acc_x{0.0};
    double acc_y{0.0};
    double radius{0.0};
};*/

void ob_traj(double delta_t, bool isDynamic, Global_variables& gl) // update position of the obstacles
{
    for (int i = 0; i < gl.no_ob; i++)
    {
        double v = isDynamic ? std::sqrt(delta_t) * 2 : 0.0;

        gl.traj_ob[i].pos_x = gl.traj_ob[i].pos_x + v * delta_t;
        gl.traj_ob[i].pos_y = gl.traj_ob[i].pos_y;
        gl.traj_ob[i].vel_x = v;
        gl.traj_ob[i].vel_y = 0.0;
        gl.traj_ob[i].acc_x = 0.0;
        gl.traj_ob[i].acc_y = 0.0;
    }
}