/** This header file contains necessary data structures
*   for code transplant from Matlab to C++
*
*   Yue Kang and Yushu Yu, 
*   Revere Lab, Chalmers / GU, 2018
*/

struct FB_state // feedback states, line 7-12 in .m file
{
    double xp_dot;  // longitudinal speed
    double yp_dot;  // lateral speed
    double psi_dot; //
    double epsi;    //
    double lp;      // lateral position
    double s;       // longitudinal position
};

struct Obstacle // line 18-26
{
    double pos_x;
    double pos_y;
    double vel_x;
    double vel_y;
    double acc_x;
    double acc_y;
    double radius;
}

struct Coefficient // the return value "out", line  463-480
{
    double norm_relpos;
    double h_angle_moving;
    double A_n_angle_moving; 
    double B_n_angle_moving;    
    double h_angle_fix; 
    double A_n_angle_fix;
    double B_n_angle_fix;
    double h_dis;
    double A_n_dis;
    double B_n_dis;
    double alert;
    double h_sid_pos;
    double A_n_side_pos;
    double b_n_side_pos;
    double h_sid_pos;
    double A_n_side_neg;
    double b_n_side_neg;
    double radius;
}
