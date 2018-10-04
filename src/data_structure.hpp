/** This header file contains necessary data structures
*   for code transplant from Matlab to C++
*
*   Yue Kang and Yushu Yu, 
*   Revere Lab, Chalmers / GU, 2018
*/

class FB_state // feedback states, line 7-12 in .m file
{
public:
    double xp_dot{0.0};  // longitudinal speed
    double yp_dot{0.0};  // lateral speed
    double psi_dot{0.0}; //
    double epsi{0.0};    //
    double ey{0.0};      // lateral position
    double s{0.0};       // longitudinal position
};

class Obstacle // line 18-26
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

class Coefficient // the return value "out", line  463-480
{
public:
    double norm_relpos{0.0};
    double h_angle_moving{0.0};
    double A_n_angle_moving{0.0}; 
    double B_n_angle_moving{0.0};    
    double h_angle_fix{0.0}; 
    double A_n_angle_fix{0.0};
    double B_n_angle_fix{0.0};
    double h_dis{0.0};
    double A_n_dis{0.0};
    double B_n_dis{0.0};
    double alert{0.0};
    double h_sid_pos{0.0};
    double A_n_side_pos{0.0};
    double b_n_side_pos{0.0};
    double h_sid_neg{0.0};
    double A_n_side_neg{0.0};
    double b_n_side_neg{0.0};
    double radius{0.0};
};
