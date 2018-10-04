#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "data_structure.hpp"

std::vector<Coefficient> constraint_obstacles_dynamics_complex(FB_state u, std::vector<Obstacle> traj_ob, bool &dead)
{
    std::vector<Coefficient> res{};
    
    double dis_thresh = 600;
//    vector<Eigen::RowVector2d> pos_ob_array, vel_ob_array;
//    vector<double> radius_array;
    std::vector<Obstacle> ob_array;
    for (int i = 0; i<traj_ob.size(); i++)
    {
        Obstacle ob = traj_ob[i];
        if ((ob.pos_x >= (u.s - 3)) && (std::abs(u.s - ob.pos_x) <= dis_thresh))
            ob_array.push_back(ob);
    }
    if (ob_array.size() == 0)
        ob_array.push_back(traj_ob[traj_ob.size() - 1]);
    // line 77 in .m file so far
    
//    std::vector<double> dis_2_vehicle;
//    for (int i = 0; i < ob_array.size(); i++)
//    {
//        Eigen::Vector2d rel_pos;
//        rel_pos << (u.s - ob_array[i].pos_x), (u.ey - ob_array[i].pos_y);
//        dis_2_vehicle.push_back(rel_pos.norm());
//    }
//    // line 98 so far
    
    std::sort(ob_array.begin(), ob_array.end(), [u](Obstacle a, Obstacle b){
        // Lambda function for comparaison
        Eigen::Vector2d rp_a, rp_b;
        rp_a << (u.s - a.pos_x), (u.ey - a.pos_y);
        rp_b << (u.s - b.pos_x), (u.ey - b.pos_y);
        return (rp_a.norm() < rp_b.norm());
    }); // line 104 so far
    
    if (ob_array.size()>1)
    {
        if (std::abs(ob_array[0].pos_x - ob_array[1].pos_x) < 0.0001)
            ob_array[1].pos_x += (std::abs(ob_array[1].pos_x) >= (std::abs(ob_array[1].pos_y))) ? -0.01 : 0.01;
    } // line 117 in .m so far
    
    dead = false;
    
    // parameters and constants
    double ck = 1.0, ey_pos = 3.2, ey_neg = -3.2, a_m = 4.0;
    double a = 1.41, b = 1.576, mu = 0.5, Fzf = 21940.0/2, Fzr = 21940.0/2;
    double cf = 65000.0, cr = 65000.0, m = 2194.0, Iz = 4770.0;
    double psi_dot_com = 0.0, p = Iz / (m * b);
    // line 140 so far
    
    //TODO: start from here
    
    return res;
}
