#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>

// #include <Eigen/Dense>

#include "data_structure.hpp"

using namespace std;

/**
    vector<Eigen::Vector2d> pos_ob_array_pre;
    vector<double> radius_pre; 

class Obstacle
{
public:
    double pos_x{0.0};
    double pos_y{0.0};
    double vel_x{0.0};
    double vel_y{0.0};
    double acc_x{0.0};
    double acc_y{0.0};
    double radius{0.0};
};
*/

void generate_init_ob(Global_variables &gl) // pre_store variables not used
{
    int no_ob = gl.no_ob;
    bool flag_ok = false;
    srand((int)time(NULL));
    for (int i = 0; i < no_ob; i++)
    {
        Obstacle curr;
        do
        {
            curr.radius = 2 + 1.5 * (double)rand() / RAND_MAX;
            curr.pos_x = 80 + 100 * (double)rand() / RAND_MAX;
            curr.pos_y = -1.7 + 3.4 * (double)rand() / RAND_MAX;
        }
        while (curr.isConf(gl.traj_ob));
        gl.traj_ob.push_back(curr);
    }
}

