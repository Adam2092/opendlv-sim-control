#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "data_structure.hpp"

vector<Coefficient> constraint_obstacles_dynamics_complex(FB_state u, vector<Obstacle> traj_ob)
{
    vector<Coefficient> res;
    
    double dis_thresh = 600;
//    vector<Eigen::RowVector2d> pos_ob_array, vel_ob_array;
//    vector<double> radius_array;
    vector<Obstacle> ob_array;
    for (int i = 0; i<traj_ob.size(); i++)
    {
        Obstacle ob = traj_ob[i];
        if ((ob.pos_x >= (u.s - 3)) && (std::abs(u.s - ob.pos_x) <= dis_thresh))
            ob_array.push_back(ob);
    }
    if (ob_array.size() == 0)
        ob_array.push_back(traj_ob[traj_ob.size() - 1]);
    // line 77 in .m file so far
    
    vector<double> dis_2_vehicle;
    for (int i = 0; i < ob_array.size(); i++)
    {
        Eigen::Vector2d rel_pos;
        rel_pos << (u.s - ob_array[i].pos_x), (u.ey - ob_array[i].pos_y);
        dis_2_vehicle.push_back(rel_pos.norm());
    }
    // line 98 in .m so far
    
    //TODO: start from here
    
    return res;
}
