#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <cstdlib>

#include <Eigen/Dense>

#include "data_structure.hpp"

using namespace std;

void bicycle_sim(double time, Global_variables& gl) // "time" replaced "T" in .m file to avoid duplicated usage (T as template)
{
    if (time < 0.0) time = 12.0; // TODO: this should be replaced by delta_t, determined by freq
    FB_state y(19, 0, 0, 0, 0, 0, 0, 0);
    gl.dead = false;
    // t_ctrl and u_ctrl are now std::vector
    gl.u_global << 0.0, 0.0;
    // gl.scale_tracking = 0;
    // gl.scale_record = 0;
    gl.tra_com_pre << 0, 0, 0;
    gl.tra_com_dot_pre << 0, 0, 20.0;
    gl.tra_com_ddot_pre << 0, 0, 0;
    gl.brake_flag = false;
    gl.brake_flag_pre = false;

    /*TODO: check the dependencies, parameters and return values of the following:
        [t1,y1_nom, y1_actual, u1]=self_solverdynamics(@bicycle_nominal_ode,  tspan, y0, options, current_hdl, 
            @bicycle_actual_ode,    current_hdl_actual);
    with the following realisation:
        [t1,y1_nom, y1_actual, u1]=self_solverdynamics(ode_nom,               tspan, y0, options, current_hdl_nom, 
            ode_actual,             current_hdl_actual)
    */
    while(1)
    {
        /*
            [dyk1, u_nomk1] = feval(ode_nom, t, y_nom,  current_hdl_nom);
            [dyk2, u_nomk2] = feval(ode_nom, t+0.5*dt, y_nom+0.5*dt*dyk1,  current_hdl_nom);
            [dyk3, u_nomk3] = feval(ode_nom, t+0.5*dt, y_nom+0.5*dt*dyk2,   current_hdl_nom);
            [dyk4, u_nomk4] = feval(ode_nom, t+dt, y_nom+dt*dyk3,  current_hdl_nom);
            y_nom=y_nom+ dt/6*(dyk1+2*dyk2+2*dyk3+dyk4);
            u_nom = u_nomk4;
        */
        /*
            Eigen::VectorXd bicycle_nominal_ode(double delta_t, FB_state y, Eigen::Vector2d &u, Global_variables &gl)
            [dy, u] = bicycle_nominal_ode(t, y, ctrl_hdl)
            [u] = virtual_Control(t, y, trajd, traj_ob)
            // ctrl_hdl == current_hdl = virtual_Control
        */
        FB_state y_nom = y, y_actual = y;

        // 4th-order Runge-Kutta method
        Eigen::VectorXd dyk1(8), dyk2(8), dyk3(8), dyk4(8);
        Eigen::Vector2d u_nom;
        dyk1 = bicycle_nominal_ode(delta_t,         y_nom,                          u_nom, gl);
        dyk2 = bicycle_nominal_ode(delta_t * 1.5,   y_nom + 0.5 * delta_t * dyk1,   u_nom, gl);
        dyk3 = bicycle_nominal_ode(delta_t * 1.5,   y_nom + 0.5 * delta_t * dyk2,   u_nom, gl);
        dyk4 = bicycle_nominal_ode(delta_t * 2,     y_nom + delta_t * dyk3,         u_nom, gl);

        y_nom += delta_t * (dyk1 + 2 * dyk2 + 2 * dyk3 + dyk4) / 6;
        if (y_nom.xp_dot < 1e-2)
        {
            y_nom.xp_dot = 0.0;
            u_nom << 0.0, 0.0;
        }

        // 1st-order solver

    }
}

Eigen::Vector2d virtual_control(Global_variables& gl)
{
    Eigen::Vector2d u;
    if (0 == gl.scale)
    {
        // "horizon = 1" unused
        Output_safety correct = safety_certificate_complex(u, gl);
        gl.nosolution = !(correct.hasSolution);
        u = correct.x;
        gl.u_global = u;
    }
    else u = gl.u_global;

    gl.scale = (gl.scale < 20) ? gl.scale + 1 : 0;

    return u;
}

Eigen::VectorXd bicycle_nominal_ode(double delta_t, FB_state y, Eigen::Vector2d &u, Global_variables &gl)
{
    // if (0 == gl.scale_record)
    //    gl.trajd = traj_gen(t, y); // not yet written

    // TODO: the following "false" should be determined by terminal command
    ob_traj(delta_t, false, gl);
    // u = feval (ctrl_hdl, t, y, trajd, traj_ob);
    u = virtual_control(gl);

    // if (gl.nosolution) // commented out

    Eigen::VectorXd dy(8);

    if (y.xp_dot <= 1e-2)
    {
        dy = Eigen::VectorXd::Zero(8);
    }
    else
    {
        // double delta_f = u(0), a_x = u(1);
        double xp_dot = (y.xp_dot <= 1e-2) ? 1e-2 : y.xp_dot;
        double yp_dot = y.yp_dot, psi_dot = y.psi_dot;
        double epsi = y.epsi, ey = y.ey, s = y.s; 
        double steer = y.steer, acc = y.acc;

        // constants
        double a = 1.41, b = 1.576, mu = 0.5, Fzf = 21940.0/2, Fzr = 21940.0/2;
        double cf = 65000.0, cr = 65000.0, m = 2194.0, Iz = 4770.0;
        double psi_dot_com = 0.0, p = Iz / (m * b);

        Eigen::Matrix2d ka;
        ka << 10.0, 0.0, 0.0, 10.0;

        Eigen::VectorXd f_x(6);
        f_x <<  yp_dot * psi_dot,
                -2 * (cf + cr) / (m * xp_dot) * yp_dot - 2 * (a * cf - b * cr) / m / xp_dot * psi_dot - xp_dot * psi_dot,
                -2 * (a * cf - b * cr) / Iz / xp_dot * yp_dot - 2 * (a * a * cf + b * b * cr) / Iz / xp_dot * psi_dot,
                psi_dot - psi_dot_com,
                yp_dot * cos(epsi) + xp_dot * sin(epsi),
                xp_dot * cos(epsi) - yp_dot * sin(epsi);

        Eigen::Matrix<double, 6, 2> g_x;
        g_x <<  0, 1,
                2 * cf / m, 0,
                2 * a * cf / Iz, 0,
                Eigen::MatrixXd::Zero(3, 2);
        Eigen::Vector2d tempV2d;
        tempV2d << steer, acc;
        dy.head(6) = f_x + g_x * tempV2d;
        dy.tail(2) = -ka * tempV2d + ka * u;
    }

    // if (0 == gl.scale_record)
    // {
    //     gl.t_ctrl.push_back(t);
    //     Eigen::VectorXd tempVd;
    //     tempVd << u, trajd;
    // } // scale_record, t_ctrl and u_ctrl are not used, therefore not implemented here

    return dy;
}

Eigen::Vector2d tracking_control(double delta_t, FB_state y, FB_state y_nom, Eigen::Vector2d u_nom)
{
    // TODO: finish this function
    Eigen::Vector2d u;
    u << 0.0, 0.0;
    return u;
}

Eigen::VectorXd bicycle_actual_ode(double delta_t, FB_state y, FB_state y_nom, Eigen::Vector2d &u, Eigen::Vector2d u_nom, Global_variables &gl)
{
    Eigen::VectorXd dy(8);

    u = tracking_control(delta_t, y, y_nom, u_nom);
    // distur = 2*[1*rand(1,1)-0.5; 1*rand(1,1)-0.5; 1*rand(1,1)-0.5; 0.2*rand(1,1)-0.1; 2*rand(1,1)-1; 2*rand(1,1)-1];

    if (y.xp_dot <= 1e-2)
    {
        dy = Eigen::VectorXd::Zero(8);
    }
    else
    {
        // double delta_f = u(0), a_x = u(1);
        double xp_dot = (y.xp_dot <= 1e-2) ? 1e-2 : y.xp_dot;
        double yp_dot = y.yp_dot, psi_dot = y.psi_dot;
        double epsi = y.epsi, ey = y.ey, s = y.s; 
        double steer = y.steer, acc = y.acc;

        // constants
        double a = 1.41, b = 1.576, mu = 0.5, Fzf = 21940.0/2, Fzr = 21940.0/2;
        double cf = 65000.0, cr = 65000.0, m = 2194.0, Iz = 4770.0;
        double psi_dot_com = 0.0, p = Iz / (m * b);

        Eigen::Matrix2d ka;
        ka << 10.0, 0.0, 0.0, 10.0;

        Eigen::VectorXd f_x(6);
        f_x <<  yp_dot * psi_dot,
                -2 * (cf + cr) / (m * xp_dot) * yp_dot - 2 * (a * cf - b * cr) / m / xp_dot * psi_dot - xp_dot * psi_dot,
                -2 * (a * cf - b * cr) / Iz / xp_dot * yp_dot - 2 * (a * a * cf + b * b * cr) / Iz / xp_dot * psi_dot,
                psi_dot - psi_dot_com,
                yp_dot * cos(epsi) + xp_dot * sin(epsi),
                xp_dot * cos(epsi) - yp_dot * sin(epsi);

        Eigen::Matrix<double, 6, 2> g_x;
        g_x <<  0, 1,
                2 * cf / m, 0,
                2 * a * cf / Iz, 0,
                Eigen::MatrixXd::Zero(3, 2);
        Eigen::Vector2d tempV2d;
        tempV2d << steer, acc;
        dy.head(6) = f_x + g_x * tempV2d;
        dy.tail(2) = -ka * tempV2d + ka * u;

        srand((int)time(NULL));
        Eigen::VectorXd distur(6);
        distur << 2 * (double)rand() / RAND_MAX - 1.0,
            2 * (double)rand() / RAND_MAX - 1.0,
            2 * (double)rand() / RAND_MAX - 1.0,
            0.4 * (double)rand() / RAND_MAX - 0.2,
            4 * (double)rand() / RAND_MAX - 2.0,
            4 * (double)rand() / RAND_MAX - 2.0;
        dy.head(6) += distur;
    }
    return dy;
}