function output_safety = safety_certificate_complex(u, traj_ob)
%calculate the safet_certificate of the system

%u: the input and the state of the system and the obstacles 

%feedback states:
xp_dot = u(1);  %lateral speed
yp_dot = u(2);  %longitudinal speed
psi_dot = u(3); 
epsi = u(4);
ey= u(5);  %lateral position
s = u(6);  %logitudinal position 

%reference trajectory
tra_com = [u(7);u(8);u(9)];  %epsi, ey, s, notice the variables are in this order, velocity and acc are the same 
tra_com_dot = [u(10);u(11);u(12)];
tra_com_ddot = [u(13);u(14);u(15)]; 
time = u(16); 

%constants: 
a = 1.41; 
b = 1.576; 
mu =0.5; 
Fzf = 21940/2; 
Fzr = 21940/2; 
cf = 65000; 
cr = 65000; 
m = 2194; 
Iz = 4770; 
psi_dot_com = 0;
p =Iz/(m*b);

%nominal input, if the output is the position x, y and the yaw angle,
%contains 3 variables: 
% feedback linearization control to obtain the nominal input: 
% k1= 9;
% k2=2*1.5*sqrt(k1);    %very important, k2 and k1 must be tuned together. 
% %2018-08-15, tunning parameters: 
% k1= diag([50,20,9]);
% k2=2*diag([4, 1.414, 1.414])*sqrt(k1);    %very important, k2 and k1 must be tuned together. 
% % k1= diag([40,2,9]);
% % k2=2*diag([4, 1.414, 1.414])*sqrt(k1);    %very important, k2 and k1 must be tuned together. 
%  
% L_f_output =[ psi_dot - psi_dot_com
%  yp_dot*cos(epsi) + xp_dot*sin(epsi)
%  xp_dot*cos(epsi) - yp_dot*sin(epsi)];
% L_f_f_output = [- (yp_dot*(2*a*cf - 2*b*cr))/(Iz*xp_dot) - (psi_dot*(2*cf*a^2 + 2*cr*b^2))/(Iz*xp_dot)
%  (psi_dot - psi_dot_com)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) + psi_dot*yp_dot*sin(epsi)
%  sin(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - (psi_dot - psi_dot_com)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + psi_dot*yp_dot*cos(epsi)];
% L_g_f_output = [[         (2*a*cf)/Iz,         0]
% [  (2*cf*cos(epsi))/m, sin(epsi)]
% [ -(2*cf*sin(epsi))/m, cos(epsi)]]; 
% u_nom_lin=tra_com_ddot-k1*([epsi; ey; s]-tra_com)-k2*(L_f_output-tra_com_dot);
% u_nom = pinv(L_g_f_output)*(u_nom_lin-L_f_f_output);   %feedback linearization


%if the output is only the x and y position, 20180815:
k1= diag([9,3]);
k2=2*diag([1.414, 1.414])*sqrt(k1);    %very important, k2 and k1 must be tuned together. 
L_f_output = [ yp_dot*cos(epsi) + xp_dot*sin(epsi)
 xp_dot*cos(epsi) - yp_dot*sin(epsi) ];
L_f_f_output = [ (psi_dot - psi_dot_com)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) + psi_dot*yp_dot*sin(epsi)
 sin(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - (psi_dot - psi_dot_com)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + psi_dot*yp_dot*cos(epsi)];
L_g_f_output = [[(2*cf*cos(epsi))/m, sin(epsi)]
[ -(2*cf*sin(epsi))/m, cos(epsi)]];
u_nom_lin=tra_com_ddot(2:3)-k1*([ ey; s]-tra_com(2:3))-k2*(L_f_output-tra_com_dot(2:3));
u_nom = inv(L_g_f_output)*(u_nom_lin-L_f_f_output);   %feedback linearization

%bound for control
alpha=[1.0; 4];

min_ = -0.5;
max_= 0.5;

delta_min = max([-1; -max_ + (yp_dot + a*psi_dot)/xp_dot ]);
delta_max = min([1; -min_ + (yp_dot + a*psi_dot)/xp_dot ]);

% delta_min = -1;
% delta_max = 1;


%flag to determin if use the bounded input:
flag_bound = 0; 


%% for dynamic obstacles:  
global results_2;   %output of   constraint_obstacles_dynamics_complex
% results_2 = constraint_obstacles_dynamics([p_x; p_y; v; psi], time); 
results_2 = constraint_obstacles_dynamics_complex([xp_dot; yp_dot; psi_dot; epsi; ey; s], time, traj_ob ); 

no_ob_active = size(results_2,2); 
global beta_2;   %initial value is 0, 100-by-1
% global  radius_pre; 
alert = 0;
for i_ob = 1:no_ob_active
    %justify when jump: 
    Ds = 1.1; 
    Ds = results_2(i_ob).radius+0.5;
    theta_d_big = asin((Ds)/results_2(i_ob).norm_relpos) - asin( (Ds -0.1) /results_2(i_ob).norm_relpos);
%     theta_d_big =0.1;
    theta_d_small = theta_d_big/1000; 
    if (beta_2(i_ob) == 0 ) && (results_2(i_ob).h_angle_fix > -theta_d_small)
        beta_2(i_ob) = 1;
    elseif (beta_2(i_ob) == 1 ) && (results_2(i_ob).h_angle_fix <=  -theta_d_big)
        beta_2(i_ob) = 0;
    end
    
    if (results_2(i_ob).alert==1)
        alert=1;
    end
    
end
  
shreshold_movingangle = 1e-20; 

%the variable determine which CBF is active now 
slack_mult = zeros(2,no_ob_active);
%number of the possible conditions: 
nu_combine =  1; 
order = [];
for aa = 1:no_ob_active
    if(beta_2(aa)==1) && (results_2(aa).h_angle_moving<=shreshold_movingangle)
        Ds = 1.1; 
        Ds = results_2(aa).radius+0.5;
        theta_d_big = asin((Ds)/results_2(aa).norm_relpos) - asin( (Ds-0.1) /results_2(aa).norm_relpos);
%         theta_d_big = 0.1;
        theta_d_small = theta_d_big/2;
        %if does not point to the obstacle:
    
        if  (results_2(aa).h_angle_fix>= -theta_d_big)   %pointing constraint
            slack_mult(1, no_ob_active) = 1;   %active
        else
            slack_mult(1, no_ob_active) = 0;
        end

        if (  results_2(aa).h_dis >= 0)   %distance constraint
            slack_mult(2, no_ob_active) = 1;   %active
        else
            slack_mult(2, no_ob_active) = 0 ;
        end 
        if (  slack_mult(1, no_ob_active) == 0) && (  slack_mult(2, no_ob_active) == 0)
            slack_mult(2, no_ob_active) =1; %at least one should be 1 
            alert = 1;
        end
    
        nu_combine = nu_combine*sum(slack_mult(:,no_ob_active));   %number of the possible conditions
    
        row_order = size(order,1);
        if (row_order == 0)
            row_order = 1;
        end
        
        if(slack_mult(1, no_ob_active) == 1) && (slack_mult(2, no_ob_active) ==1)
            order = [order, ones(row_order, 1); order, 2*ones(row_order, 1)];  %record the place 
        elseif(slack_mult(1, no_ob_active) == 1) && (slack_mult(2, no_ob_active) == 0)
            order = [order, ones(row_order, 1)];  %record the place 
        elseif(slack_mult(1, no_ob_active) == 0) && (slack_mult(2, no_ob_active) == 1)
            order = [order, 2*ones(row_order, 1)];  %record the place 
        end
    else
        row_order = size(order,1);
        if (row_order == 0)
            row_order = 1;
        end
        order = [order, zeros(row_order, 1)];
    end
end

%the minmal value and the corresponding solution, notice there may be no
%solution:
value_min = 100000000;
x_min = [0;0];

for i_combine = 1:nu_combine
    A_n_and = [];
    b_n_and = [];
    A_n_and = [results_2(aa).A_n_side_pos; results_2(aa).A_n_side_neg];
    b_n_and = [results_2(aa).b_n_side_pos; results_2(aa).b_n_side_neg];
    A_n_or = [];
    b_n_or = [];
    
    for aa = 1:no_ob_active
        if (beta_2(aa) == 0 ) 
% %             if pointing to the obstacle, both conditions should be satisfied 
            A_n_and = [A_n_and;  results_2(aa).A_n_angle_fix; results_2(aa).A_n_dis];
            b_n_and = [b_n_and;   results_2(aa).B_n_angle_fix;  results_2(aa).B_n_dis ]; 
        elseif (results_2(aa).h_angle_moving > shreshold_movingangle)
            %if the vehicle and the obstacle have been in the opposite
            %direction
            A_n_and = [A_n_and;  results_2(aa).A_n_angle_fix; results_2(aa).A_n_angle_moving];
            b_n_and = [b_n_and;   results_2(aa).B_n_angle_fix;  results_2(aa).B_n_angle_moving ];
            
        else            
            if(order(i_combine, aa) ==1)  %pointing constraint
%                  A_n_or = [A_n_or;  results_2(aa).A_n_angle_fix; ]; 
%                  b_n_or = [b_n_or;   results_2(aa).B_n_angle_fix;  ];
%                  
%                  A_n_or = [A_n_or;  results_2(aa).A_n_angle_fix;  results_2(aa).A_n_dis;]; 
%                  b_n_or = [b_n_or;   results_2(aa).B_n_angle_fix; results_2(aa).B_n_dis; ];
                 A_n_or = [A_n_or;  results_2(aa).A_n_angle_fix  ]; 
                 b_n_or = [b_n_or;   results_2(aa).B_n_angle_fix ]; 
                 
            elseif (order(i_combine, aa) == 2)    %distance constraint
                 A_n_or = [A_n_or;  results_2(aa).A_n_dis  ]; 
                 b_n_or = [b_n_or;   results_2(aa).B_n_dis ];    
%                  A_n_or = [A_n_or;  results_2(aa).A_n_dis; results_2(aa).A_n_angle_fix ]; 
%                  b_n_or = [b_n_or;   results_2(aa).B_n_dis; results_2(aa).B_n_angle_fix]; 
            end 
        end
 
    end
   
     %solve QP at the end, see if the angle constraints for multiple
     %obstacles solvable 
     H= diag([1;1]);
     f2 = -2* u_nom;  %the optimal goal is for the entire control
     optoption_1 = optimset('Display', 'off', 'TolFun', 1e-6);
     if (size(A_n_and,1)>0)
        if(flag_bound ==0)
%              [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n_and, b_n_and, [], [], -alpha, alpha, [], optoption_1);
            [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n_and, b_n_and, [], [], [delta_min; -alpha(2)], [delta_max; alpha(2)], [], optoption_1);
        else
            [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n_and, b_n_and, [], [], [], [], [], optoption_1);
        end        
%         delta_just = width_control(A_n_and, b_n_and);   %calucate the width of the feasible control  
        if (EXITFLAG<0)  
            %qp  has no solution 
            A_n_and = [];
            b_n_and = [];
            A_n_and = [results_2(aa).A_n_side_pos; results_2(aa).A_n_side_neg];
            b_n_and = [results_2(aa).b_n_side_pos; results_2(aa).b_n_side_neg];
            for aa = 1:no_ob_active
                if (beta_2(aa) == 0 ) && (aa == 1)
        % %             if angle CBF for multiple obstacles does not solvable,
        % then only consider the angle constraints for the main obstacle 
                    A_n_and = [A_n_and;  results_2(aa).A_n_angle_fix;  results_2(aa).A_n_dis];
                    b_n_and = [b_n_and;   results_2(aa).B_n_angle_fix;  results_2(aa).B_n_dis ];
                elseif (beta_2(aa) == 0 ) && (aa ~= 1) 
                    A_n_and = [A_n_and;  results_2(aa).A_n_dis; ];
                    b_n_and = [b_n_and;   results_2(aa).B_n_dis; ];
                elseif (results_2(aa).h_angle_moving > shreshold_movingangle) 
                    %if the vehicle and the obstacle have been in the opposite
                    %direction
                    A_n_and = [A_n_and;  results_2(aa).A_n_angle_fix; results_2(aa).A_n_angle_moving];
                    b_n_and = [b_n_and;   results_2(aa).B_n_angle_fix;  results_2(aa).B_n_angle_moving ];
                end
            end

        end
     end
     
    global A_n b_n;    
    A_n = [A_n_and; A_n_or];
    b_n = [b_n_and; b_n_or];
    
     %solve QP at the end

%     delta_just2 = width_control(A_n, b_n);   %calucate the width of the feasible control 
    
    %see if solvable 
%     if(delta_just2.max< 0 )
        %if the QP is solvable, then solve it 
        if(flag_bound ==0)
%             [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n, b_n, [], [], -alpha-u_nom, alpha-u_nom, [], optoption_1);
%             [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n, b_n, [], [], -alpha, alpha, [], optoption_1);
            [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n, b_n, [], [], [delta_min; -alpha(2)], [delta_max; alpha(2)], [], optoption_1);
        else
            [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n, b_n, [], [], [], [], [], optoption_1);
        end
        if (EXITFLAG < 0)  
            %qp  has no solution 
            FVAL = 100000000; %no solution, set a very big value 
        end

        if (FVAL < value_min)
            value_min = FVAL; %update
            x_min = x; 
            A_min = A_n; 
            b_min = b_n;
        end   
 
%     end
    
%     if(delta_just2.max < max_delta_lb) 
%         max_delta_lb = delta_just2.max; 
%         x_min = x; 
%     end
end

global dead; 
% the output: 
if (value_min~=100000000) && (alert==0) &&(dead == 0)
% if (max_delta_lb<0)
    %select the minimal solution 
     out= x_min; 
%      A_min = A_n; 
%      b_min = b_n;
    output_safety.nosolution = 0;
else
    %no solution exsits, brake  
    out = [0; -alpha(2)];
    output_safety.nosolution = 1;
end






%if test
testflag = 0;
if (testflag ==1)
%     global results_2; 
%     results_2 = constraint_obstacles_dynamics([p_x; p_y; v; psi], 0); 
%     result_dynamic = constraint_obstacles_dynamics([p_x; p_y; v; psi], time); 
%     no_ob_dynamics = size(result_dynamic,2);
% 
%     %define empty matrix 
%     A_n = [];
%     b_n =[];    
%     for aa = 1:no_ob_dynamics
%         A_n = [A_n; result_dynamic(aa).A_n_dis];
%         b_n = [b_n; result_dynamic(aa).B_n_dis];  
%     end
% 
%     f2 = -2* u_nom;  %the optimal goal is for the entire control: norm(x-bar x), notice 
% 
%     H =eye(2,2);
% flag_bound = 0; 
% alpha= [4;4];
%     if(flag_bound ==0)
%         [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n, b_n, [], [], -alpha, alpha, []);
%     else
%         [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n, b_n, [], [], [], [], []);
%     end
% 
%     value_min = 100;
%     if (EXITFLAG < 0)  
%     %qp  has no solution 
%         out = [-alpha(1); 0];   %no solution, brake  
%     else 
%         out = x; 
%     end
end 

%  A_n = [  results_2(1).A_n_angle_fix; results_2(1).A_n_angle_moving];
%  b_n = [ results_2(1).B_n_angle_fix;  results_2(1).B_n_angle_moving ]; 
%  [x, FVAL, EXITFLAG] = quadprog(H, f2, A_n, b_n, [], [], [], [], [], optoption_1);

% if (time>6.8)
%     out = u_nom;
% end

% global output_safety; 
output_safety.out = [out; ...
    results_2(1).A_n_angle_fix(1); results_2(1).A_n_angle_fix(2); results_2(1).B_n_angle_fix; ...
    results_2(1).A_n_angle_moving(1); results_2(1).A_n_angle_moving(2); results_2(1).B_n_angle_moving; ...
    results_2(1).A_n_dis(1); results_2(1).A_n_dis(2); results_2(1).B_n_dis; ...
    results_2(1).h_angle_fix; results_2(1).h_angle_moving;  results_2(1).h_dis; beta_2(1); value_min];
% output_safety = [out; h_move_angle; A_n(1,1); A_n(1,2); b_n(1); zeros(5,1)];
 

