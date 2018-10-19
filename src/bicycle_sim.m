%% CLF-CBF-QP geometric control of a single 3D-moving quadrotor
%Yu Yushu
% ------------------------------------------------
% The corresponding data of system trajectory
% is stored in mat file 
% ------------------------------------------------
% T: the time duration of simulation 
%
% 

function bicycle_sim(T, count)
clc;
% close all;
if nargin < 1
    T = 12; 
end


%% Simulation based on ODE function
% control function handles
ctrl_hdl1 = @virtual_Control;
current_hdl = ctrl_hdl1;
ctrl_hdl_str = func2str(current_hdl);

%actual tracking control: 
ctrl_hdl1 = @tracking_Control;
current_hdl_actual = ctrl_hdl1;
ctrl_hdl_str_actual = func2str(current_hdl_actual);

% initial condition of different trials
% -------------------------------------------------
 %20181011: state order: 8-by-1, including the steering angle and acc. 
 %x dot, y dot, psi dot, e_psi, psi, ey, s, bar steer, bar acc. 

% y0=[0.5;0; -1; -1.0; 2; -4];
y0=[19; 0; 0; 0; 0; 0; 0; 0];

global dead;
dead = 0 ; %once there is vilation of constraints, trigger dead,
global t_ctrl; 
global u_ctrl; 
t_ctrl = [0];
u_ctrl = [0; 0; zeros(9, 1)]; % the control value and also the reference value 


global scale  u_global scale_tracking u_tracking_global scale_record trajd; 
scale = 0;
u_global = 0;
scale_tracking = 0;
scale_record= 0;
global tra_com_pre tra_com_dot_pre tra_com_ddot_pre; %record the previous trajectory point 
tra_com_pre = [0;0;0];
tra_com_dot_pre = [0;0;20]; %notice the initial value
tra_com_ddot_pre = [0;0;0];

global brake_flag brake_flag_pre; %initial time of brake 
brake_flag = 0;
brake_flag_pre = 0;

%generate the initial position of the obstacles, used in simulation, note
%the initial position of the obstacles are global variables  
global pos_ob_array_pre radius_pre;
% generate_init_ob(count);
global traj_ob_seris; %the trajectory of the obstacles 
traj_ob_seris = pos_ob_array_pre;

% option of ode function 
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

% simulation process
disp('The simulation process has started.');
disp(strcat('Controller: ', ctrl_hdl_str));
disp(strcat('Controller: ', ctrl_hdl_str_actual));
disp('---------------------------------------------');
tspan = [0 T];
% [t1, y1] = ode15s(@quad_3d_ode, tspan, y0, options, current_hdl);
% [t1, y1] = ode45(@quad_3d_ode, tspan, y0, options, current_hdl);
% [t1, y1] = self_solverdynamics(@bicycle_nominal_ode, tspan, y0, options, current_hdl);

%2018-09-24, modifiy, adding the real plant, the previous is the nominal
%value: 
[t1,y1_nom, y1_actual, u1]=self_solverdynamics(@bicycle_nominal_ode, tspan, y0, options, current_hdl, @bicycle_actual_ode, current_hdl_actual);


% save all the data into .mat file  
save('sim_data.mat');
clear u_global flag_mode scale; 
global name;  

save(name);
disp('Data successfully stored!');


end


%% nominal  control
function [u] = virtual_Control(t, y, trajd, traj_ob)
global scale  u_global;
global x_nom u_nom; 
global flag_mode; 

%trajectory of obstacles: traj_ob

if (scale==0)
%run mpc 
    %time horizon: 
    horizon = 1; 
    correct =  safety_certificate_complex([y; trajd; t], traj_ob);
    global nosolution; 
    nosolution = correct.nosolution;
    u = correct.out(1:2);
    u_global =u; 

    % check simulation time for stability property
    debug = 1;
    if debug == 1
        disp(['The current time is ', num2str(t)]);

    end
% -----------------------------------------------
else 
    u = u_global;   
end
 
scale = scale+1;
if (scale == 20)
    scale = 0;
end
end



%% Ode Function of this vehicle
function [dy, u] = bicycle_nominal_ode(t, y, ctrl_hdl)

dy =zeros(8,1);

global scale_record trajd;  %used in 4th-order Range-Kutta method, because traj_gen() can only run once in each step 

if(scale_record==0)    
    % get the current reference and control input 
    trajd = traj_gen(t, y);
end
%get the current information of the obstacles: 
traj_ob = ob_traj(t); 
u = feval(ctrl_hdl, t, y, trajd, traj_ob);

global nosolution; 
 
no_iter = 0; 
if(no_iter==0) && (nosolution==1)
    %repeat again if there is no solution: 
%     trajd = traj_gen(t, y);    
%     scale = 0; %back 
%     u = feval(ctrl_hdl, t, y, trajd, traj_ob);  %be careful, may induce error in time elapse 
end


%%2018-08-04, bicycle model: 
%input signal 
delta_f = u(1);   %steering angle 
a_x = u(2);    %acc 

if (y(1)<= 1e-2)
    y(1) = 1e-2;
end 

%states:
xp_dot = y(1);  %lateral speed
yp_dot = y(2);  %longitudinal speed
psi_dot = y(3); 
epsi = y(4);
ey= y(5);  %lateral position
s = y(6);  %logitudinal position 

%20181011:
steer = y(7);   %steering angle 
acc = y(8);  %longitudinal acc 

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

ka = diag([10, 10]); %parameters in the actuator dynamics, added 20181011

%state equation: 
f_x = [ yp_dot*psi_dot;... 
    -2*(cf+cr)/(m*xp_dot)*yp_dot-2*(a*cf-b*cr)/m/xp_dot*psi_dot-xp_dot*psi_dot; ...
     -2*(a*cf-b*cr)/Iz/xp_dot*yp_dot-2*(a*a*cf+b*b*cr)/Iz/xp_dot*psi_dot;...
     psi_dot - psi_dot_com;...
     yp_dot*cos(epsi) + xp_dot*sin(epsi); ...
     xp_dot*cos(epsi)-yp_dot*sin(epsi)];
 
g_x = [0, 1; ...
    2*cf/m, 0; ...
    2*a*cf/Iz, 0;...
    0, 0;...
    0, 0;...
    0, 0];

%20181011, consider actuator dynamics 
% dy = f_x + g_x*u;
dy = [f_x + g_x*[steer; acc];  -ka*[steer; acc]] + [zeros(6,2); ka]*u; 

if (y(1)<= 1e-2)
    dy = zeros(8,1); 
end 

%record control:  
global t_ctrl; 
global u_ctrl;
global traj_ob_seris;

if (scale_record == 0)
    %%record data,
    t_ctrl = [t_ctrl; t];
    u_ctrl = [u_ctrl, [u;trajd] ];
    traj_ob_seris = [traj_ob_seris; traj_ob.pos];
end
scale_record = scale_record+1;
if (scale_record==4)
    scale_record=0;
end
end

%% actual tracking control: 
% u = feval(current_hdl_actual, t, y, y_nom, u_nom );
function [u] = tracking_Control(t, y, y_nom, u_nom )
%actual tracking controller
%t: time
%y: state of the real system, n-by-1, here 6-by-1
%y_nom: nominal state: n-by-1, here 6-by-1
%u_nom: nominal control: m-by-1, here 2-by-1

global scale_tracking u_tracking_global;

if (scale_tracking==0)
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

    %feedback states:
    xp_dot = y(1);  %lateral speed
    yp_dot = y(2);  %longitudinal speed
    psi_dot = y(3); 
    epsi = y(4);
    ey= y(5);  %lateral position
    s = y(6);  %logitudinal position 

    if (xp_dot>1e-1)  %longitudinal velocity should be bigger than 0, otherwise feedback linearization will be fail
        L_f_output = [ yp_dot*cos(epsi) + xp_dot*sin(epsi)
         xp_dot*cos(epsi) - yp_dot*sin(epsi) ];
        L_f_f_output = [ (psi_dot - psi_dot_com)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) + psi_dot*yp_dot*sin(epsi)
         sin(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - (psi_dot - psi_dot_com)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + psi_dot*yp_dot*cos(epsi)];
        L_g_f_output = [[(2*cf*cos(epsi))/m, sin(epsi)]
        [ -(2*cf*sin(epsi))/m, cos(epsi)]];
        p = [ y(5); y(6)];  %output
        p_dot = L_f_output;  %derivative of output

        %nominal state:
        xp_dot = y_nom(1);  %lateral speed
        yp_dot = y_nom(2);  %longitudinal speed
        psi_dot = y_nom(3); 
        epsi = y_nom(4);
        ey= y_nom(5);  %lateral position
        s = y_nom(6);  %logitudinal position 

        if (xp_dot>1e-1)  %longitudinal velocity should be bigger than 0, otherwise feedback linearization will be fail
            L_f_output_nom = [ yp_dot*cos(epsi) + xp_dot*sin(epsi)
             xp_dot*cos(epsi) - yp_dot*sin(epsi) ];
            L_f_f_output_nom  = [ (psi_dot - psi_dot_com)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) + psi_dot*yp_dot*sin(epsi)
             sin(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - (psi_dot - psi_dot_com)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + psi_dot*yp_dot*cos(epsi)];
            L_g_f_output_nom  = [[(2*cf*cos(epsi))/m, sin(epsi)]
            [ -(2*cf*sin(epsi))/m, cos(epsi)]];
            p_nom = [y_nom(5); y_nom(6)];  % output, nominal 
            p_nom_dot = L_f_output_nom;  %derivative of output, nominal 

            %tracking error for the nomnial position: 
            p_err = p-p_nom;
            p_err_dot = p_dot - p_nom_dot; 

            
%             y_nom = [20; 0; 0; 0; 0; 20*t]; %test only
            xi_err = y-y_nom;  %state error


            %if the output is only the x and y position, 20180815:
            k1= diag([1,3]);
            k2=2*diag([1.414, 1.414])*sqrt(k1);    %very important, k2 and k1 must be tuned together. 
            
%             k1= diag([0.04,0.3]);
            k2=2*diag([1.414, 1.414])*sqrt(k1);    

            %actual control: 
            u= L_g_f_output\( -k1*p_err -k2*p_err_dot - L_f_f_output + L_g_f_output_nom*u_nom + L_f_f_output_nom);
            u= ( -k1*p_err -k2*p_err_dot  + u_nom );

            %calculated from LQR: 
            K =[ 0.0000    0.6257    0.6908    4.9809    1.0000    0.0000
                 1.7321   -0.0000    0.0000    0.0000    0.0000    1.0000];
%             u = -diag([1.5 0.9])*K*xi_err + u_nom;

    %         u= L_g_f_output\(  -k1*p_err   -k2*p_err_dot   + L_g_f_output_nom*u_nom );
        else 
            u = u_nom;
        end
    else
        u = u_nom;
    end
    u_tracking_global = u; 
else
    u = u_tracking_global; 
end


% scale_tracking = scale_tracking+1;
if (scale_tracking == 10)
    scale_tracking = 0;
end

end


%%  [dy, u] = feval(ode_actual, t, y_actual, y_nom, u_nom, current_hdl_actual);
function [dy, u] = bicycle_actual_ode(t, y, y_nom, u_nom, current_hdl_actual)


dy =zeros(8,1);

%actual tracking control: 
u = feval(current_hdl_actual, t, y, y_nom, u_nom);

%disturbance:
distur = 2*[1*rand(1,1)-0.5; 1*rand(1,1)-0.5; 1*rand(1,1)-0.5; 0.2*rand(1,1)-0.1; 2*rand(1,1)-1; 2*rand(1,1)-1];
 
 

% 2018-09-24, bicycle model: 
%input signal 
% delta_f = u(1);   %steering angle 
% a_x = u(2);    %acc 

if (y(1)<= 1e-2)
    y(1) = 1e-2;
end 

%states:
xp_dot = y(1);  %lateral speed
yp_dot = y(2);  %longitudinal speed
psi_dot = y(3); 
epsi = y(4);
ey= y(5);  %lateral position
s = y(6);  %logitudinal position 

%20181011:
steer = y(7);   %steering angle 
acc = y(8);  %longitudinal acc 

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

ka = diag([10, 10]); %parameters in the actuator dynamics, added 20181011

%state equation: 
f_x = [ yp_dot*psi_dot;... 
    -2*(cf+cr)/(m*xp_dot)*yp_dot-2*(a*cf-b*cr)/m/xp_dot*psi_dot-xp_dot*psi_dot; ...
     -2*(a*cf-b*cr)/Iz/xp_dot*yp_dot-2*(a*a*cf+b*b*cr)/Iz/xp_dot*psi_dot;...
     psi_dot - psi_dot_com;...
     yp_dot*cos(epsi) + xp_dot*sin(epsi); ...
     xp_dot*cos(epsi)-yp_dot*sin(epsi)];
 
g_x = [0, 1; ...
    2*cf/m, 0; ...
    2*a*cf/Iz, 0;...
    0, 0;...
    0, 0;...
    0, 0];

 

%20181011, consider actuator dynamics 
% dy = f_x + g_x*u  + distur;
dy = [f_x + g_x*[steer; acc];  -ka*[steer; acc]] + [zeros(6,2); ka]*u + [distur; 0; 0]; 

if (y(1)<= 1e-2)
    dy = zeros(8,1); 
end 

end

% %% 
% function out = traj_gen(t,y)
% 
% % r=7.5;
% % pos = [r*cos(1/r*t); r*sin(1/r*t); 0];
% % vel = [-sin(1/r*t); cos(1/r*t); 0];
% % acc = [-1/r*cos(1/r*t); -1/r*sin(1/r*t); 0];
% % dacc =[1/r^2*sin(1/r*t); -1/r^2*cos(1/r*t); 0]; 
% % d2acc = [1/r^3*cos(1/r*t); 1/r^3*sin(1/r*t); 0];
% 
% %reference trajectory
% % tra_com = [u(7);u(8);u(9)];  %epsi, ey, s, notice the variables are in this order, velocity and acc are the same 
% % tra_com_dot = [u(10);u(11);u(12)];
% % tra_com_ddot = [u(13);u(14);u(15)]; 
% % time = u(16); 
% 
% 
% global tra_com_pre tra_com_dot_pre tra_com_ddot_pre; 
% xp_dot = y(1);  %lateral speed
% yp_dot = y(2);  %longitudinal speed
% psi_dot = y(3); 
% epsi = y(4);
% ey= y(5);  %lateral position
% s = y(6);  %logitudinal position 
% psi_dot_com = 0;
% L_f_output =[ psi_dot - psi_dot_com
%  yp_dot*cos(epsi) + xp_dot*sin(epsi)
%  xp_dot*cos(epsi) - yp_dot*sin(epsi)];
% err = [([epsi; ey; s]-tra_com_pre); 0*(L_f_output-tra_com_dot_pre)];
% err_dot = [0*([epsi; ey; s]-tra_com_pre); 1*(L_f_output-tra_com_dot_pre)];
% norm_err = norm(err);
% global nosolution;
% if (nosolution == 0)
%     ks = 0.0;
% else
%     ks = 0.0;
%     nosolution = 0;
% end
% virtual_time = exp(-ks*norm_err*norm_err);
% virtual_time_dot = -ks*virtual_time*2*err'*err_dot;
% 
% v=20*virtual_time;
% v_dot = 20*virtual_time_dot;
% p = [v*t; 0]; 
% psi = 0;
% 
% tra_com = [0; 0 ; v*t];  %epsi, ey, s, notice the variables are in this order, velocity and acc are the same 
% tra_com_dot = [0; 0; v];
% tra_com_ddot = [0; 0; 0]; 
% global dt;
% deltat = dt;
% tra_com = tra_com_pre+[0; 0 ; v]*deltat;  %epsi, ey, s, notice the variables are in this order, velocity and acc are the same 
% tra_com_dot = [0; 0; v];
% % tra_com_ddot = (tra_com_dot - tra_com_dot_pre)/deltat; 
% tra_com_ddot = [0; 0; v_dot];
% 
% %update: 
% tra_com_pre = tra_com;
% tra_com_dot_pre = tra_com_ddot;
% tra_com_ddot_pre = tra_com_ddot;
%  
% 
% out=[tra_com; tra_com_dot; tra_com_ddot];
% end


 
