#include <iostream>
#include <vector>
#include <cmath>

#include <Eigen/Dense>

#include "data_structure.hpp"
#include "qpOASES/qpOASES.hpp"

using namespace std;

Output_safety safety_certificate_complex(FB_state u, std::vector<Eigen::Vector3d> tra, vector<Obstacle> traj_ob, vector<bool> &beta_2)
{
    Output_safety out;

    double xp_dot = u.xp_dot, yp_dot = u.yp_dot, psi_dot = u.psi_dot;
    double epsi = u.epsi, ey = u.ey, s = u.s; 

    Eigen::Vector3d tra_com = tra[0], tra_com_dot = tra[1], tra_com_ddot = tra[2];

    // constants
    double a = 1.41, b = 1.576, mu = 0.5, Fzf = 21940.0/2, Fzr = 21940.0/2;
    double cf = 65000.0, cr = 65000.0, m = 2194.0, Iz = 4770.0;
    double psi_dot_com = 0.0, p = Iz / (m * b);
    // line 31 so far

    Eigen::Matrix2d k1, k2;
    k1 << 9, 0, 0, 3;
    k2 << 6 * sqrt(2), 0, 0, 2 * sqrt(6);

    Eigen::Vector2d L_f_out, L_f_f_output;
    L_f_out << yp_dot * cos(epsi) + xp_dot * sin(epsi), xp_dot * cos(epsi) - yp_dot * sin(epsi);
    double tempD = (psi_dot * xp_dot + (psi_dot * 2 * (a * cf - b * cr)) / (m * xp_dot) 
        + (yp_dot * 2 * (cf + cr)) / (m * xp_dot));
    double err_psi_dot = psi_dot - psi_dot_com;
    L_f_f_output << err_psi_dot * (xp_dot * cos(epsi) - yp_dot * sin(epsi)) 
                        - cos(epsi) * tempD + psi_dot * yp_dot * sin(epsi),
                    sin(epsi) * tempD - err_psi_dot * (yp_dot * cos(epsi) + xp_dot * sin(epsi)) 
                        + psi_dot * yp_dot * cos(epsi);

    Eigen::Matrix2d L_g_f_output;
    L_g_f_output << 2 * cf * cos(epsi) / m, sin(epsi), -2 * cf * sin(epsi) / m, cos(epsi);

    Eigen::Vector2d u_nom_lin, u_nom;
    //u_nom_lin=tra_com_ddot(2:3)-k1*([ ey; s]-tra_com(2:3))-k2*(L_f_output-tra_com_dot(2:3));
    Eigen::Vector2d tempV2d, tempTail0, tempTail1, tempTail2;
    tempV2d << ey, s;
    tempTail0 << tra_com(1), tra_com(2);
    tempTail1 << tra_com_dot(1), tra_com_dot(2);
    tempTail2 << tra_com_ddot(1), tra_com_ddot(2);
    // u_nom_lin = tra_com_ddot.tail(2) - k1 * (tempV2d - tra_com.tail(2)) - k2 * (L_g_f_output - tra_com_dot.tail(2));
    u_nom_lin = tempTail2 - k1 * (tempV2d - tempTail0) - k2 * (tempV2d - tempTail1);
    // line 67 so far

    Eigen::Vector2d alpha;
    alpha << 1.0, 4.0;
    tempD = (yp_dot + a * psi_dot) / xp_dot;
    double delta_min = (tempD - 0.5 > -1.)? tempD - 0.5 : -1.0;
    double delta_max = (tempD + 0.5 < 1.0)? tempD + 0.5 : 1.0;

    bool flag_bound = false, dead = false, alert = false;
    // line 83 so far

    vector<Coefficient> results_2 = constraint_obstacles_dynamics_complex(u, traj_ob, &dead);
    int no_ob_active = results_2.size();
    int nu_combine = 1;

    double shreshold_movingangle = pow(10, -20);
    // Eigen::Vector2d slack_mult;
    Eigen::MatrixXi slack_mult(2, no_ob_active);

    for (int i = 0; i < no_ob_active; ++i)
    {
        double Ds = results_2[i].radius + 0.5;
        double theta_d_big = asin(Ds / results_2[i].norm_relpos) - asin((Ds - 0.1) / results_2[i].norm_relpos);
        double theta_d_small = theta_d_big / 1000;

        // TODO: double-check if beta_2 is initialized as "all false"
        if ((!beta_2[i]) && (results_2[i].h_angle_fix > -theta_d_small)) beta_2[i] = true;
        else if (beta_2[i] && (results_2[i].h_angle_fix <= -theta_d_big)) beta_2[i] = false;

        if (results_2[i].alert) alert = true;
        // line 112 so far


        if (beta_2[i] && (results_2[i].h_angle_moving <= shreshold_movingangle))
        {
            theta_d_small = theta_d_big / 2;
            slack_mult(0, i) = (results_2[i].h_angle_fix >= -theta_d_big) ? 1 : 0;
            slack_mult(1, i) = (results_2[i].h_dis >= 0) ? 1 : 0;
            if ((0 == slack_mult(0, i)) && (0 == slack_mult(1, i)))
            {
                slack_mult(1, i) = 1; // at least one should be 1
                alert = true;
            }
            //line 144 so far, "i" used to be "no_ob_active" in the original .m file

            ////////////////////////////////////
            // line 146 - 153 need discussion //
            ////////////////////////////////////

            nu_combine *= (slack_mult(0, i) + slack_mult(1, i));
        }
    }

    Eigen::MatrixXi order(nu_combine, no_ob_active);
    order.setOnes();
    int tempI = 0;
    for (int i = 0; i < no_ob_active; ++i)
    {
        if ((0 == slack_mult(0, i)) && (1 == slack_mult(1, i)))
        {
            order.col(no_ob_active) *= 2;
        }
        else if ((1 == slack_mult(0, i)) && (1 == slack_mult(1, i)))
        {
            tempI++; // the I-th case with 2 possibilities
            int length = nu_combine / pow(2, tempI);
            for (int j = 0; j < pow(2, tempI - 1); j++)
            {
                // fill in the sub-matrices by 2s
                order.block((2 * j + 1) * length, no_ob_active, length, 1) *= 2;
            }
        }
    }
    // line 167 so far

    double value_min = 100000000;
    // x_min = [0;0];
    double x_min[2] = {0.0, 0.0};

    for (int i = 0; i < nu_combine; i++) // i <-> i_combine in .m file
    {
        vector<Eigen::Vector2d> A_n_and{}, A_n_or{};
        vector<double> b_n_and{}, b_n_or{};

        A_n_and.push_back(results_2[0].A_n_side_pos);
        A_n_and.push_back(results_2[0].A_n_side_neg);
        b_n_and.push_back(results_2[0].b_n_side_pos);
        b_n_and.push_back(results_2[0].b_n_side_neg);

        for (int j = 0; j < no_ob_active; j++) // j <-> aa in .m file
        {
            if (!beta_2[j])
            {
                A_n_and.push_back(results_2[j].A_n_angle_fix);
                A_n_and.push_back(results_2[j].A_n_dis);
                b_n_and.push_back(results_2[j].b_n_angle_fix);
                b_n_and.push_back(results_2[j].b_n_dis);
            }
            else if (results_2[j].h_angle_moving > shreshold_movingangle)
            {
                A_n_and.push_back(results_2[j].A_n_angle_fix);
                A_n_and.push_back(results_2[j].A_n_angle_moving);
                b_n_and.push_back(results_2[j].b_n_angle_fix);
                b_n_and.push_back(results_2[j].b_n_angle_moving);
            } // line 192 so far
            else
            {
                if (order(i, j) == 1)
                {
                    A_n_or.push_back(results_2[j].A_n_angle_fix);
                    b_n_or.push_back(results_2[j].b_n_angle_fix);
                }
                else if (order(i, j) == 2)
                {
                    A_n_or.push_back(results_2[j].A_n_dis);
                    b_n_or.push_back(results_2[j].b_n_dis);
                }
            }
        } // for (j) 
        // line 211 so far
        {
            using namespace qpOASES;
            SQProblem qp(2, 1);
            real_t H[4] = {1.0, 0.0, 0.0, 1.0};
            real_t f2[2] = {-2 * u_nom(0), -2 * u_nom(1)};
            real_t rtA_n_and[A_n_and.size() * 2], rtb_n_and[b_n_and.size()];
            real_t *rtNullprt = NULL;
            real_t rtOut[2];
            for (int i = 0; i < A_n_and.size(); i++)
            {
                rtA_n_and[2 * i] = A_n_and[i][0];
                rtA_n_and[2 * i + 1] = A_n_and[i][1];
            }
            for (int i = 0; i < b_n_and.size(); i++)
                rtb_n_and[i] = b_n_and[i];
            real_t rtlb[2]{delta_min, -alpha(1)}, rtub[2]{delta_max, alpha(1)};
            int nWSR = 10;
            if (A_n_and.size() > 0)
            {
                if (!flag_bound) 
                    qp.init(H, f2, rtA_n_and, rtlb, rtub, rtNullprt, rtb_n_and, nWSR, 0);
                // (Currently) NEVER goes to "else"
                qp.getPrimalSolution(rtOut);
                // line 224 so far
                if (!qp.isSolved())
                {
                    A_n_and.clear();
                    b_n_and.clear();
                    A_n_and.push_back(results_2[0].A_n_side_pos);
                    A_n_and.push_back(results_2[0].A_n_side_neg);
                    b_n_and.push_back(results_2[0].b_n_side_pos);
                    b_n_and.push_back(results_2[0].b_n_side_neg);
                    for (int j = 0; j < no_ob_active; j++) // j <-> aa in .m file
                    {
                        if (!beta_2[j])
                        {
                            if (0 == j)
                            {
                                A_n_and.push_back(results_2[j].A_n_angle_fix);
                                A_n_and.push_back(results_2[j].A_n_dis);
                                b_n_and.push_back(results_2[j].b_n_angle_fix);
                                b_n_and.push_back(results_2[j].b_n_dis);
                            }
                            else
                            {
                                A_n_and.push_back(results_2[j].A_n_dis);
                                b_n_and.push_back(results_2[j].b_n_dis);
                            }
                        }
                        else if (results_2[j].h_angle_moving > shreshold_movingangle)
                        {
                            A_n_and.push_back(results_2[j].A_n_angle_fix);
                            A_n_and.push_back(results_2[j].A_n_angle_moving);
                            b_n_and.push_back(results_2[j].b_n_angle_fix);
                            b_n_and.push_back(results_2[j].b_n_angle_moving);
                        } 
                    } // for (j) 
                } // if (!isSolved)
            }// if (size > 0)
            // line 251 so far

            A_n_and.insert(A_n_and.end(), A_n_or.begin(), A_n_or.end());
            b_n_and.insert(b_n_and.end(), b_n_or.begin(), b_n_or.end());

            real_t rtA_new[A_n_and.size() * 2], rtb_new[b_n_and.size()];
            for (int i = 0; i < A_n_and.size(); i++)
            {
                rtA_new[2 * i] = A_n_and[i][0];
                rtA_new[2 * i + 1] = A_n_and[i][1];
            }
            for (int i = 0; i < b_n_and.size(); i++)
                rtb_new[i] = b_n_and[i];
            qp.init(H, f2, rtA_new, rtlb, rtub, rtNullprt, rtb_new, nWSR, 0);
            // line 269 so far
            qp.getPrimalSolution(rtOut);
            double FVAL = (qp.isSolved()) ? (double)(qp.getObjVal()) : 100000000.0;

            if (FVAL < value_min)
            {
                value_min = FVAL;
                x_min[0] = rtOut[0];
                x_min[1] = rtOut[1];
                // Removed unused matrices A_min and b_min
            }
        } // using namespace qpOASES
    } // for (i = 0 to nu_combine)
    // line 288 so far

    if((value_min < 100000000.0) && (!alert) && (!dead))
    {
        out.x[0] = x_min[0];
        out.x[1] = x_min[1];
        out.hasSolution = true;
    }
    else
    {
        out.x[0] = 0.0;
        out.x[1] = alpha(1);
        out.hasSolution = false;
    }
    out.value_min = value_min;
    out.coef = results_2[0];
    return out;
}