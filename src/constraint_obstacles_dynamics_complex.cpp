#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "data_structure.hpp"

using namespace std;

vector<Coefficient> constraint_obstacles_dynamics_complex(FB_state u, Global_variables& gl)
{
    vector<Coefficient> res{};
    double xp_dot = u.xp_dot, yp_dot = u.yp_dot, psi_dot = u.psi_dot;
    double epsi = u.epsi, ey = u.ey, s = u.s; 
    vector<Obstacle> traj_ob = gl.traj_ob;
    
    double dis_thresh = 600;
//    vector<Eigen::RowVector2d> pos_ob_array, vel_ob_array;
//    vector<double> radius_array;
    vector<Obstacle> ob_array;
    for (int i = 0; i<traj_ob.size(); i++)
    {
        Obstacle ob = traj_ob[i];
        if ((ob.pos_x >= (u.s - 3)) && (abs(u.s - ob.pos_x) <= dis_thresh))
            ob_array.push_back(ob);
    }
    if (ob_array.size() == 0)
        ob_array.push_back(traj_ob[traj_ob.size() - 1]);
    
//    vector<double> dis_2_vehicle;
//    for (int i = 0; i < ob_array.size(); i++)
//    {
//        Eigen::Vector2d rel_pos;
//        rel_pos << (u.s - ob_array[i].pos_x), (u.ey - ob_array[i].pos_y);
//        dis_2_vehicle.push_back(rel_pos.norm());
//    }
    
    sort(ob_array.begin(), ob_array.end(), [u](Obstacle a, Obstacle b){
        // Lambda function for comparaison
        Eigen::Vector2d rp_a, rp_b;
        rp_a << (u.s - a.pos_x), (u.ey - a.pos_y);
        rp_b << (u.s - b.pos_x), (u.ey - b.pos_y);
        return (rp_a.norm() < rp_b.norm());
    }); // line 104 so far
    
    if (ob_array.size()>1)
    {
        if (abs(ob_array[0].pos_x - ob_array[1].pos_x) < 0.0001)
            ob_array[1].pos_x += (abs(ob_array[1].pos_x) >= (abs(ob_array[1].pos_y))) ? -0.01 : 0.01;
    }
    
    gl.dead = false;
    
    // parameters and constants
    double ck = 1.0, ey_pos = 3.2, ey_neg = -3.2, a_m = 4.0;
    double a = 1.41, b = 1.576, mu = 0.5, Fzf = 21940.0/2, Fzr = 21940.0/2;
    double cf = 65000.0, cr = 65000.0, m = 2194.0, Iz = 4770.0;
    double psi_dot_com = 0.0, p = Iz / (m * b);
    // line 140 so far
    
    double h_sid_pos = ey_pos - ey - pow((yp_dot * cos(epsi) + xp_dot * sin(epsi)), 2) / (2 * a_m);
    double L_f_h_sid_pos = ((yp_dot * cos(epsi) + xp_dot * sin(epsi)) * ( 2 * cf * yp_dot * cos(epsi) 
        + 2 * cr * yp_dot * cos(epsi) - a_m * m * xp_dot + 2 * a * cf * psi_dot * cos(epsi) 
        - 2 * b * cr * psi_dot * cos(epsi) + m * psi_dot_com * pow(xp_dot, 2) * cos(epsi) 
        - m * psi_dot_com * xp_dot * yp_dot * sin(epsi))) / (a_m * m * xp_dot);
    Eigen::Vector2d L_g_h_sid_pos;
    L_g_h_sid_pos << -( 2 * cf * cos(epsi) * (yp_dot*cos(epsi) + xp_dot * sin(epsi))) / (a_m * m),
        -(sin(epsi) * (yp_dot * cos(epsi) + xp_dot * sin(epsi))) / a_m;

    double h_sid_neg = ey - ey_neg - pow((yp_dot * cos(epsi) + xp_dot * sin(epsi)), 2) / (2 * a_m);
    double L_f_h_sid_neg = ((yp_dot * cos(epsi) + xp_dot * sin(epsi)) * (2 * cf * yp_dot * cos(epsi) 
        + 2 * cr * yp_dot * cos(epsi) + a_m * m * xp_dot + 2 * a * cf * psi_dot * cos(epsi) 
        - 2 * b * cr * psi_dot * cos(epsi) + m * psi_dot_com * pow(xp_dot, 2) * cos(epsi) 
        - m * psi_dot_com * xp_dot * yp_dot * sin(epsi))) / (a_m * m * xp_dot);
    Eigen::Vector2d L_g_h_sid_neg;
    L_g_h_sid_neg << -( 2 * cf * cos(epsi) * (yp_dot * cos(epsi) + xp_dot * sin(epsi))) / (a_m * m), 
        -(sin(epsi) * (yp_dot * cos(epsi) + xp_dot * sin(epsi))) / a_m;
    // line 154 so far

    double temp_d = yp_dot * cos(epsi) + xp_dot * sin(epsi);
    Eigen::Vector2d A_n_side_pos, A_n_side_neg;
    double b_n_side_pos, b_n_side_neg;

    if ((h_sid_pos > 0) && (temp_d < 0))
    {
        A_n_side_pos << 0, 0;
        b_n_side_pos - 1;
    }
    else
    {
        A_n_side_pos = -L_g_h_sid_pos;
        b_n_side_pos = L_f_h_sid_pos + 0.5 * h_sid_pos;
    }
    
    if ((h_sid_neg > 0) && (temp_d > 0))
    {
        A_n_side_neg << 0, 0;
        b_n_side_neg - 1;
    }
    else
    {
        A_n_side_neg = -L_g_h_sid_neg;
        b_n_side_neg = L_f_h_sid_neg + 0.5 * h_sid_neg;
    }
    if (abs(ey) >= 3.7) gl.dead = true;
    // line 171 so far

    for (int i = 0; i < ob_array.size(); i++)
    {
        double pos_ob_x = ob_array[i].pos_x;
        double pos_ob_y = ob_array[i].pos_y;
        double vel_ob_x = ob_array[i].vel_x;
        double vel_ob_y = ob_array[i].vel_y;
        double acc_ob_x = ob_array[i].acc_x;
        double acc_ob_y = ob_array[i].acc_y;

        Eigen::Vector2d v_vehicle, rel_pos, rel_vel;
        v_vehicle << xp_dot * cos(epsi) - yp_dot * sin(epsi), yp_dot * cos(epsi) + xp_dot * sin(epsi);
        rel_pos << s - pos_ob_x, ey - pos_ob_y;
        rel_vel << v_vehicle(0) - vel_ob_x, v_vehicle(1) - vel_ob_y;
        // line 201 so far

        double Ds = ob_array[i].radius + 0.5; 
        double cos_rel_ang = -(rel_pos.adjoint() * rel_vel)(0);
        cos_rel_ang /= rel_pos.norm();
        cos_rel_ang /= rel_vel.norm();

        double h_ang = 3.0;
        Eigen::Vector2d A_n_angle_fix;
        A_n_angle_fix << 0, 0;
        double b_n_angle_fix = 1.0;

        if (cos_rel_ang > -0.99)
        {
            cos_rel_ang = (cos_rel_ang >= 1)? 1: cos_rel_ang;
            double rel_ang = acos(cos_rel_ang);
            double ratio = Ds / rel_pos.norm();
            if (ratio >= 1)
            {
                ratio = 0.9999999;
                gl.dead = true;
            }
            h_ang = rel_ang - asin(ratio);
            // line 238 so far

            pos_ob_y = (abs(pos_ob_y) < 0.0001) ? pos_ob_y / abs(pos_ob_y) * 0.0001 : pos_ob_y;
            epsi = (abs(epsi) < 0.00001) ? epsi / abs(epsi) * 0.00001 : epsi;
            // line 249 so far

            double L_f_h_ang_part1 = -((psi_dot - psi_dot_com) 
                    * (pow(xp_dot, 2) + pow(yp_dot, 2) - vel_ob_x * xp_dot * cos(epsi) - vel_ob_y * yp_dot * cos(epsi) 
                        - vel_ob_y * xp_dot * sin(epsi) + vel_ob_x * yp_dot * sin(epsi)) 
                    * (ey * vel_ob_x + pos_ob_x * vel_ob_y - pos_ob_y * vel_ob_x - s * vel_ob_y 
                        - ey * xp_dot * cos(epsi) + pos_ob_y * xp_dot * cos(epsi) 
                        - pos_ob_x * yp_dot * cos(epsi) + s * yp_dot * cos(epsi) 
                        + ey * yp_dot * sin(epsi) - pos_ob_x * xp_dot * sin(epsi) 
                        - pos_ob_y * yp_dot * sin(epsi) + s * xp_dot * sin(epsi))) 
                / (sqrt(1 - pow(((pos_ob_x - s) * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi)) 
                    + (ey - pos_ob_y) * (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi))), 2) 
                    / ((pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                        * (pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                        + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2))))
                    * sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                        * pow(pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                        + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2), 1.5)) 
                - ((ey - pos_ob_y) * (xp_dot * cos(epsi) - yp_dot * sin(epsi)) 
                    * (ey * vel_ob_x + pos_ob_x * vel_ob_y - pos_ob_y * vel_ob_x - s * vel_ob_y 
                        - ey * xp_dot * cos(epsi) + pos_ob_y * xp_dot * cos(epsi) 
                        - pos_ob_x * yp_dot * cos(epsi) + s * yp_dot * cos(epsi) 
                        + ey * yp_dot * sin(epsi) - pos_ob_x * xp_dot * sin(epsi) 
                        - pos_ob_y * yp_dot * sin(epsi) + s * xp_dot * sin(epsi)))
                    / (sqrt(1 - pow((pos_ob_x - s) * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi)) 
                        + (ey - pos_ob_y) * (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi)), 2) 
                        / ((pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                            * (pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                                + (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2)))) 
                * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5) 
                * sqrt(pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                    + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2))) 
                - ((pos_ob_x - s) * (yp_dot * cos(epsi) + xp_dot * sin(epsi)) 
                    * (ey * vel_ob_x + pos_ob_x * vel_ob_y - pos_ob_y * vel_ob_x - s * vel_ob_y 
                        - ey * xp_dot * cos(epsi) + pos_ob_y * xp_dot * cos(epsi) 
                        - pos_ob_x * yp_dot * cos(epsi) + s * yp_dot * cos(epsi) 
                        + ey * yp_dot * sin(epsi) - pos_ob_x * xp_dot * sin(epsi) 
                        - pos_ob_y * yp_dot * sin(epsi) + s * xp_dot * sin(epsi)))
                / (sqrt(1 - pow((pos_ob_x - s) * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi)) 
                    + (ey - pos_ob_y) * (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi)), 2)
                / ((pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                    * (pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                        + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2))))
                * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5) 
                * sqrt(pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                    + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2))) 
                - ((vel_ob_x * cos(epsi) - xp_dot + vel_ob_y * sin(epsi)) 
                    * (m * psi_dot * pow(xp_dot, 2) + 2 * cf * yp_dot + 2 * cr * yp_dot 
                        + 2 * a * cf * psi_dot - 2 * b * cr * psi_dot) 
                    * (ey * vel_ob_x + pos_ob_x * vel_ob_y - pos_ob_y * vel_ob_x - s * vel_ob_y 
                        - ey * xp_dot * cos(epsi) + pos_ob_y * xp_dot * cos(epsi) 
                        - pos_ob_x * yp_dot * cos(epsi) + s * yp_dot * cos(epsi) 
                        + ey * yp_dot * sin(epsi) - pos_ob_x * xp_dot * sin(epsi) 
                        - pos_ob_y * yp_dot * sin(epsi) + s * xp_dot * sin(epsi)))
                    / (m * xp_dot * sqrt(1 - pow((pos_ob_x - s) * (vel_ob_x - xp_dot * cos(epsi) 
                            + yp_dot * sin(epsi)) + (ey - pos_ob_y) * (yp_dot * cos(epsi) 
                            - vel_ob_y + xp_dot * sin(epsi)), 2) 
                        /((pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                            * (pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                                + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2))))
                    * sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                    * pow(pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                        + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2), 1.5));
            // line 269 so far
            double L_t_h_ang_part1 = -((ey * vel_ob_x + pos_ob_x * vel_ob_y - pos_ob_y * vel_ob_x - s * vel_ob_y 
                    - ey * xp_dot * cos(epsi) + pos_ob_y * xp_dot * cos(epsi) - pos_ob_x * yp_dot * cos(epsi) 
                    + s * yp_dot * cos(epsi) + ey * yp_dot * sin(epsi) - pos_ob_x * xp_dot * sin(epsi) 
                    - pos_ob_y * yp_dot * sin(epsi) + s * xp_dot * sin(epsi)) 
                * (pos_ob_y * pow(vel_ob_x, 3) - pos_ob_x * pow(vel_ob_y, 3) - ey * pow(vel_ob_x, 3) 
                    + s * pow(vel_ob_y, 3) - acc_ob_x * pow(pos_ob_x, 2) * vel_ob_y 
                    + acc_ob_y * pow(pos_ob_x, 2) * vel_ob_x - acc_ob_x * pow(pos_ob_y, 2) * vel_ob_y 
                    + acc_ob_y * pow(pos_ob_y,2 ) * vel_ob_x - acc_ob_x * pow(s, 2) * vel_ob_y 
                    + acc_ob_y * pow(s, 2) * vel_ob_x - ey * vel_ob_x * pow(vel_ob_y, 2) 
                    - ey * vel_ob_x * pow(xp_dot, 2) - ey * vel_ob_x * pow(yp_dot, 2) 
                    - pos_ob_x * pow(vel_ob_x, 2) * vel_ob_y + pos_ob_y * vel_ob_x * pow(vel_ob_y, 2) 
                    - pos_ob_x * vel_ob_y * pow(xp_dot, 2) + pos_ob_y * vel_ob_x * pow(xp_dot, 2) 
                    + s * pow(vel_ob_x, 2) * vel_ob_y - pos_ob_x * vel_ob_y * pow(yp_dot, 2) 
                    + pos_ob_y * vel_ob_x * pow(yp_dot, 2) + s * vel_ob_y * pow(xp_dot, 2) 
                    + s * vel_ob_y * pow(yp_dot, 2) - acc_ob_x * pow(ey, 2) * vel_ob_y 
                    + acc_ob_y * pow(ey, 2) * vel_ob_x + 2 * acc_ob_x * ey * pos_ob_y * vel_ob_y 
                    - 2 * acc_ob_y * ey * pos_ob_y * vel_ob_x + 2 * acc_ob_x * pos_ob_x * s * vel_ob_y 
                    - 2 * acc_ob_y * pos_ob_x * s * vel_ob_x - acc_ob_y * pow(ey, 2) * xp_dot * cos(epsi) 
                    + acc_ob_x * pow(ey, 2) * yp_dot * cos(epsi) - acc_ob_y * pow(pos_ob_x, 2) * xp_dot * cos(epsi) 
                    - acc_ob_y * pow(pos_ob_y, 2) * xp_dot * cos(epsi) + acc_ob_x * pow(pos_ob_x, 2) * yp_dot * cos(epsi) 
                    + acc_ob_x * pow(pos_ob_y, 2) * yp_dot * cos(epsi) - acc_ob_y * pow(s, 2) * xp_dot * cos(epsi) 
                    + acc_ob_x * pow(s, 2) * yp_dot * cos(epsi) + acc_ob_x * pow(ey, 2) * xp_dot * sin(epsi) 
                    + acc_ob_y * pow(ey, 2) * yp_dot * sin(epsi) + 2 * ey * pow(vel_ob_x, 2) * xp_dot * cos(epsi) 
                    + acc_ob_x * pow(pos_ob_x, 2) * xp_dot * sin(epsi) + acc_ob_x * pow(pos_ob_y, 2) * xp_dot * sin(epsi) 
                    + acc_ob_y * pow(pos_ob_x, 2) * yp_dot * sin(epsi) + acc_ob_y * pow(pos_ob_y, 2) * yp_dot * sin(epsi) 
                    + acc_ob_x * pow(s, 2) * xp_dot * sin(epsi) + acc_ob_y * pow(s, 2) * yp_dot * sin(epsi) 
                    - 2 * pos_ob_y * pow(vel_ob_x, 2) * xp_dot * cos(epsi) 
                    + 2 * pos_ob_x * pow(vel_ob_y, 2) * yp_dot * cos(epsi) 
                    - 2 * s * pow(vel_ob_y, 2) * yp_dot * cos(epsi) - 2 * ey * pow(vel_ob_x, 2) * yp_dot * sin(epsi) 
                    + 2 * pos_ob_x * pow(vel_ob_y, 2) * xp_dot * sin(epsi) 
                    + 2 * pos_ob_y * pow(vel_ob_x, 2) * yp_dot * sin(epsi) 
                    - 2 * s * pow(vel_ob_y, 2) * xp_dot * sin(epsi) - 2 * s * vel_ob_x * vel_ob_y * xp_dot * cos(epsi) 
                    + 2 * ey * vel_ob_x * vel_ob_y * xp_dot * sin(epsi) 
                    - 2 * pos_ob_y * vel_ob_x * vel_ob_y * xp_dot * sin(epsi) 
                    - 2 * pos_ob_x * vel_ob_x * vel_ob_y * yp_dot * sin(epsi) 
                    + 2 * s * vel_ob_x * vel_ob_y * yp_dot * sin(epsi) + 2 * acc_ob_y * ey * pos_ob_y * xp_dot * cos(epsi) 
                    - 2 * acc_ob_x * ey * pos_ob_y * yp_dot * cos(epsi) + 2 * acc_ob_y * pos_ob_x * s * xp_dot * cos(epsi) 
                    - 2 * acc_ob_x * pos_ob_x * s * yp_dot * cos(epsi) - 2 * acc_ob_x * ey * pos_ob_y * xp_dot * sin(epsi) 
                    - 2 * acc_ob_y * ey * pos_ob_y * yp_dot * sin(epsi) + 2 * ey * vel_ob_x * vel_ob_y * yp_dot * cos(epsi) 
                    - 2 * acc_ob_x * pos_ob_x * s * xp_dot * sin(epsi) - 2 * acc_ob_y * pos_ob_x * s * yp_dot * sin(epsi) 
                    + 2 * pos_ob_x * vel_ob_x * vel_ob_y * xp_dot * cos(epsi) 
                    - 2 * pos_ob_y * vel_ob_x * vel_ob_y * yp_dot * cos(epsi))) 
                / (sqrt(1 - pow((pos_ob_x - s) * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi)) 
                        + (ey - pos_ob_y) * (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi)), 2) 
                    / ((pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                        * (pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                            + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2)))) 
                * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5) 
                * pow(pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                    + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2), 1.5));
            // line 270 so far
            double L_g_h_ang_part1 = (2 * cf * (vel_ob_x * cos(epsi) - xp_dot + vel_ob_y * sin(epsi)) 
                    * (ey * vel_ob_x + pos_ob_x * vel_ob_y - pos_ob_y * vel_ob_x - s * vel_ob_y - ey * xp_dot * cos(epsi) 
                        + pos_ob_y * xp_dot * cos(epsi) - pos_ob_x * yp_dot*cos(epsi) + s * yp_dot * cos(epsi) 
                        + ey * yp_dot * sin(epsi) - pos_ob_x * xp_dot * sin(epsi) - pos_ob_y * yp_dot * sin(epsi) 
                        + s * xp_dot * sin(epsi))) 
                / (m * sqrt(1 - pow((pos_ob_x - s) * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi)) 
                            + (ey - pos_ob_y) * (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi)), 2) 
                        / ((pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                            * (pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                                + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2))))
                    * sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                    * pow(pow(vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi), 2) 
                        + pow(yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi), 2), 1.5));
            // line 271 so far

            double L_f_h_ang_part2 = ((ey - pos_ob_y) * (yp_dot * cos(epsi) + xp_dot * sin(epsi)))
                    / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                - ((pos_ob_x - s) * (xp_dot * cos(epsi) - yp_dot * sin(epsi))) 
                    / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2));
            double L_g_h_ang_part2 = 0.0;
            double L_t_h_ang_part2 = -(ey * vel_ob_y - pos_ob_x * vel_ob_x - pos_ob_y * vel_ob_y + s * vel_ob_x) 
                / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2));
            //line 279 so far

            temp_d = - Ds / pow(rel_pos.norm(), 2);
            double asin_dot = 1 / sqrt(1 - ratio * ratio);
            double L_f_h_ang = L_f_h_ang_part1 - asin_dot * temp_d * L_f_h_ang_part2;
            double L_g_h_ang = L_g_h_ang_part1 - asin_dot * temp_d * L_g_h_ang_part2;
            double L_t_h_ang = L_t_h_ang_part1 - asin_dot * temp_d * L_t_h_ang_part2;
            // line 288 so far

            A_n_angle_fix(0) = - L_g_h_ang;
            b_n_angle_fix = b_n_angle_fix = L_f_h_ang + L_t_h_ang + 3 * h_ang;
            // line 297 so far
        } // line 303 so far

        bool alert = false;
        double dis_maxacc_sqr = 2 * a_m * (rel_pos.norm() - Ds);
        if (dis_maxacc_sqr <= 0)
        {
            alert = true;
            gl.dead = true;
            dis_maxacc_sqr = 0.000001;
        }
        double h_vel = sqrt(dis_maxacc_sqr) + (rel_pos.adjoint() * rel_vel)(0) / rel_vel.norm();
        // line 318 so far

        double L_f_h_norm = ((ey - pos_ob_y) * (yp_dot * cos(epsi) + xp_dot * sin(epsi))) 
                / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
            - ((pos_ob_x - s) * (xp_dot * cos(epsi) - yp_dot * sin(epsi))) 
                / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2));
        Eigen::Vector2d L_g_h_norm;
        L_g_h_norm << 0, 0;
        double L_t_h_norm = -(ey * vel_ob_y - pos_ob_x * vel_ob_x - pos_ob_y * vel_ob_y + s * vel_ob_x) 
            / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)); 
        // line 345 so far

        double L_f_h_vel = (((ey - pos_ob_y) * (xp_dot * cos(epsi) - yp_dot * sin(epsi))) 
                        /sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                    + ((pos_ob_x - s) * (yp_dot * cos(epsi) + xp_dot * sin(epsi))) 
                        /sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2))) 
                * (psi_dot - psi_dot_com) 
            - (yp_dot * cos(epsi) + xp_dot * sin(epsi)) 
                * (((ey - pos_ob_y) * (2 * ey - 2 * pos_ob_y) * (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi))) 
                    / (2 * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5)) 
                        - (yp_dot*cos(epsi) - vel_ob_y + xp_dot * sin(epsi)) 
                    / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y,2 )) 
                        + ((pos_ob_x - s) * (2 * ey - 2 * pos_ob_y) 
                            * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi))) 
                    / ( 2 * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5)))
            + (xp_dot * cos(epsi) - yp_dot * sin(epsi)) 
                * (((2 * pos_ob_x - 2 * s) * (ey - pos_ob_y) * (yp_dot*cos(epsi) - vel_ob_y + xp_dot * sin(epsi))) 
                    / (2 * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5)) 
                        - (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi)) 
                    / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                        + ((2 * pos_ob_x - 2 * s) * (pos_ob_x - s) 
                            * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi)))
                    / (2 * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5))) 
            - ((cos(epsi) * (ey - pos_ob_y)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                    + (sin(epsi) * (pos_ob_x - s)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2))) 
                * (psi_dot * xp_dot + (psi_dot * (2 * a * cf - 2 * b * cr)) / (m * xp_dot) 
                    + (yp_dot * (2 * cf + 2 * cr)) / (m * xp_dot)) 
            - psi_dot * yp_dot * ((cos(epsi) * (pos_ob_x - s)) 
                    / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                - (sin(epsi) * (ey - pos_ob_y)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)));
        // line 352 so far
        Eigen::Vector2d L_g_h_vel;
        L_g_h_vel << (2 * cf * ((cos(epsi) * (ey - pos_ob_y)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                + (sin(epsi) * (pos_ob_x - s)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2))))/m, 
            (sin(epsi) * (ey - pos_ob_y)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y,2 )) 
                - (cos(epsi) * (pos_ob_x - s)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2));
        // line 353 so far
        double L_t_h_vel = vel_ob_y * (((ey - pos_ob_y) * (2 * ey - 2 * pos_ob_y) 
                        * (yp_dot * cos(epsi) - vel_ob_y + xp_dot*sin(epsi))) 
                    / (2 * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5)) 
                - (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi)) 
                    / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                + ((pos_ob_x - s) * (2 * ey - 2 * pos_ob_y) * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi))) 
                    / (2 * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5))) 
            - vel_ob_x * (((2 * pos_ob_x - 2 * s) * (ey - pos_ob_y) 
                        * (yp_dot * cos(epsi) - vel_ob_y + xp_dot * sin(epsi))) 
                    /(2 * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5)) 
                - (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi)) 
                    / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
                + ((2 * pos_ob_x - 2 * s) * (pos_ob_x - s) * (vel_ob_x - xp_dot * cos(epsi) + yp_dot * sin(epsi))) 
                    /(2 * pow(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2), 1.5))) 
            - (acc_ob_y * (ey - pos_ob_y)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2)) 
            + (acc_ob_x * (pos_ob_x - s)) / sqrt(pow(pos_ob_x - s, 2) + pow(ey - pos_ob_y, 2));
        // line 354 so far

        temp_d = a_m / sqrt(dis_maxacc_sqr);
        L_f_h_vel += temp_d * L_f_h_norm;
        L_g_h_vel += temp_d * L_g_h_norm;
        L_t_h_vel += temp_d * L_t_h_norm;

        Eigen::Vector2d A_n_vel = -L_g_h_vel;
        double b_n_vel = L_f_h_vel + L_t_h_vel + 3 * h_vel;
        // line 368 so far
        
        Coefficient c;
        c.norm_relpos = rel_pos.norm();
        c.h_angle_moving = -2;
        Eigen::Vector2d temp_v;
        temp_v << 0, 0;
        c.A_n_angle_moving = temp_v;
        c.b_n_angle_moving = 1;
        c.h_angle_fix = h_ang;
        c.A_n_angle_fix = A_n_angle_fix;
        c.b_n_angle_fix = b_n_angle_fix;
        c.h_dis = h_vel; // needs double-check
        c.A_n_dis = A_n_vel; // needs double-check
        c.b_n_dis = b_n_vel; // needs double-check
        c.alert = alert;
        c.h_sid_pos = h_sid_pos;
        c.A_n_side_pos = A_n_side_pos;
        c.b_n_side_pos = b_n_side_pos;
        c.h_sid_neg = h_sid_neg;
        c.A_n_side_neg = A_n_side_neg;
        c.b_n_side_neg = b_n_side_neg;
        c.radius = ob_array[i].radius;

        res.push_back(c);

    } // for each ob in ob_array
    return res;
}
