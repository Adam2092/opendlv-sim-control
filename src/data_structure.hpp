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
    FB_state(double a, double b, double c, double d, double e, double f, double g, double h):
        xp_dot(a), yp_dot(b), psi_dot(c), epsi(d), ey(e), s(f), steer(g), acc(h){}
    FB_state operator +(const FB_state& other)
    {
        return FB_state(
            this->xp_dot + other.xp_dot,
            this->yp_dot + other.yp_dot,
            this->psi_dot + other.psi_dot,
            this->epsi + other.epsi,
            this->ey + other.ey,
            this->s + other.s,
            this->steer + other.steer,
            this->acc + other.acc
        );
    }
    FB_state operator +(const Eigen::VectorXd& other)
    {
        return FB_state( 
            this->xp_dot + other(0),
            this->yp_dot + other(1),
            this->psi_dot + other(2)
            this->epsi + other(3),
            this->ey + other(4),
            this->s + other(5),
            this->steer + other(6),
            this->acc + other(7)
        );
    }
    FB_state& operator +=(const FB_state& other)
    {
        this->xp_dot += other.xp_dot;
        this->yp_dot += other.yp_dot;
        this->psi_dot += other.psi_dot;
        this->epsi += other.epsi;
        this->ey += other.ey;
        this->s += other.s;
        this->steer += other.steer;
        this->acc += other.acc;
        return *this;
    }
    FB_state& operator +=(const Eigen::VectorXd& other)
    {
        this->xp_dot += other(0);
        this->yp_dot += other(1);
        this->psi_dot += other(2);
        this->epsi += other(3);
        this->ey += other(4);
        this->s += other(5);
        this->steer += other(6);
        this->acc += other(7);
        return *this;
    }
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

    bool isConf(Obstacle that)
    {
        Eigen::Vector2d temp;
        temp << this->pos_x - that.pos_x, this->pos_y - that.pos_y;
        return (temp.norm() <= (this->radius + that.radius)) ? true : false;
    }

    bool isConf(std::vector<Obstacle> list)
    {
        if (0 == list.size()) return false;
        for (int i = 0; i < list.size(); i++)
        {
            if (this->isConf(list[i])) return true;
        }
    return false;
    }
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

    std::vector<Eigen::Vector3d> trajd; // this contains all the 3 variables in the following line
    //Eigen::Vector3d tra_com_pre, tra_com_dot_pre, tra_com_ddot_pre;

    bool brake_flag, brake_flag_pre;
    bool nosolution;

    // the following three are initialised by constraint_obstacles.m
    std::vector<Obstacle> traj_ob;
    // int no_ob;
    // vector<Eigen::Vector2d> pos_ob_array_pre;
    // vector<double> radius_pre;

    std::vector<bool> beta_2;
    double dt;

    // // For data sample and visualisation
    // vector<double> t_ctrl;
    // vector<Eigen::VectorXd> u_ctrl;

    // Unlisted global variables
    bool dead;
    FB_state state_brakeini;
    
}

std::vector<Coefficient> constraint_obstacles_dynamics_complex(FB_state, Global_variables &);

Output_safety safety_certificate_complex(FB_state, Global_variables &);

Eigen::Vector2d virtual_control(Global_variables &)

void bicycle_model(double);

void ob_traj(double, bool, &Global_variables);