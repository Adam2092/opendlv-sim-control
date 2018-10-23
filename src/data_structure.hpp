/** This header file contains necessary data structures
*   for code transplant from Matlab to C++
*
*   Yue Kang and Yushu Yu, 
*   Revere Lab, Chalmers / GU, 2018
*/

#include <iostream>
#include <vector>
#include <Eigen/Dense>

class FB_state
{
public:
    double xp_dot{0.0};  // longitudinal speed
    double yp_dot{0.0};  // lateral speed
    double psi_dot{0.0}; //
    double epsi{0.0};    //
    double ey{0.0};      // lateral position
    double s{0.0};       // longitudinal position
    double steer{0.0};
    double acc{0.0};
    FB_state(double, double, double, double, double, double, double, double);
};

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

class Coefficient // the return value "out", line 463-480 in constraint_obs~.m
{
public:
    double norm_relpos{0.0};
    double h_angle_moving{0.0};
    Eigen::Vector2d A_n_angle_moving; 
    double b_n_angle_moving{0.0};    
    double h_angle_fix{0.0}; 
    Eigen::Vector2d A_n_angle_fix;
    double b_n_angle_fix{0.0};
    double h_dis{0.0};
    Eigen::Vector2d A_n_dis;
    double b_n_dis{0.0};
    bool alert{false};
    double h_sid_pos{0.0};
    Eigen::Vector2d A_n_side_pos;
    double b_n_side_pos{0.0};
    double h_sid_neg{0.0};
    Eigen::Vector2d A_n_side_neg;
    double b_n_side_neg{0.0};
    double radius{0.0};
};

class Output_safety // the return value output_safety.out, line 355-359 in safety_cert~.m
{
public:
    Eigen::Vector2d x;
    Coefficient coef;
    double value_min{100000000.0};
    bool hasSolution{false}; // line 298 in safety_cert~.m, INVERSED BOOLEAN VALUE of "nosolution"
};

class Global_variables
{
public:
    int scale;
    Eigen::Vector2d u_global;
    int scale_tracking;
    Eigen::Vector2d u_tracking_global;
    int scale_record;

    Eigen::Vector3d tra_com_pre, tra_com_dot_pre, tra_com_ddot_pre;

    bool brake_flag, brake_flag_pre;
    bool nosolution;

    // int no_ob;
    /* vector pos_ob_array_pre, radius_pre; */ // Dynamic size, not to be global

    vector<bool> beta_2;
    double dt;

    vector<double> t_ctrl;
    vector<Eigen::VectorXd> u_ctrl;

    FB_state trajd;

    // Unlisted global variables
    bool dead;
    FB_state state_brakeini;
    
}

std::vector<Coefficient> constraint_obstacles_dynamics_complex(FB_state, std::vector<Obstacle>, Global_variables);

Output_safety safety_certificate_complex(FB_state, std::vector<Eigen::Vector3d>, std::vector<Obstacle>, Global_variables);

Eigen::Vector2d virtual_control(FB_state, std::vector<Eigen::Vector3d>, vector<Obstacle>, Global_variables)

void bicycle_model(double);