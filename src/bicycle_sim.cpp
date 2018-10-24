#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "data_structure.hpp"

using namespace std;

void bicycle_sim(double time, Global_variables& gl) // "time" replaced "T" in .m file to avoid duplicated usage (T as template)
{
    if (time < 0.0) time = 12.0;
    FB_state y(19, 0, 0, 0, 0, 0, 0, 0);
    gl.dead = false;
    // t_ctrl and u_ctrl are now std::vector
    gl.u_global << 0.0, 0.0;
    gl.scale_tracking = 0;
    gl.scale_record = 0;
    gl.tra_com_pre << 0, 0, 0;
    gl.tra_com_dot_pre << 0, 0, 20.0;
    gl.tra_com_ddot_pre << 0, 0, 0;
    gl.brake_flag = false;
    gl.brake_flag_pre = false;

    //TODO: check global variables pos_ob_array_pre, radius_pre, traj_ob_seris

    //TODO: check the dependencies, parameters and return values of the following:
    //[t1,y1_nom, y1_actual, u1]=self_solverdynamics(@bicycle_nominal_ode, tspan, y0, options, current_hdl, @bicycle_actual_ode, current_hdl_actual);
}

Eigen::Vector2d virtual_control(std::vector<Eigen::Vector3d> tra, vector<Obstacle> traj_ob, Global_variables& gl)
{
    Eigen::Vector2d u;
    if (0 == gl.scale)
    {
        // "horizon = 1" unused
        //(FB_state u, std::vector<Eigen::Vector3d> tra, vector<Obstacle> traj_ob, Global_variables& gl)
        Output_safety correct = safety_certificate_complex(u, tra, traj_ob, &gl);
        gl.nosolution = !(correct.hasSolution);
        u = correct.x;
        gl.u_global = u;
    }
    else u = gl.u_global;

    gl.scale = (gl.scale < 20) ? gl.scale + 1 : 0;

    return u;
}

void bicycle_nominal_ode(double t, FB_state y, Global_variables &gl)
{
    if (0 == gl.scale_record)
        gl.trajd = traj_gen(t, y); // not yet written
    vector<Obstacle> traj_ob  = ob_traj(t); // not yet written


    // u = feval (ctrl_hdl, t, y, trajd, traj_ob);
}