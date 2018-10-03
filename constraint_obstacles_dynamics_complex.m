function out = constraint_obstacles_dynamics_complex(u, t, traj_ob)
%calculate the coefficients from the multiple static obstacles 

%carefully check this file

%feedback states:
xp_dot = u(1);  %longitudinal speed
yp_dot = u(2);  %lateral speed
psi_dot = u(3); 
epsi = u(4);
ey= u(5);  %lateral position
s = u(6);  %logitudinal position 
 
%static obstacles: 
% global no_ob;
global dead;  %once there is vilation of constraints, trigger dead,
% global pos_ob_array_pre radius_pre;
no_ob_input = traj_ob.no;
pos_ob_array_pre_input = zeros(2,no_ob_input);
vel_ob_array_pre = zeros(2,no_ob_input);
acc_ob_array = zeros(2,no_ob_input);

pos_ob_array_pre_input = traj_ob.pos;
vel_ob_array_pre = traj_ob.vel;
acc_ob_array = traj_ob.acc;
radius_pre = traj_ob.rad; 
% radius = ones(no_ob,2);

%%note that if two of the obstacles have the same x-coordinates, then the
%%solver will let the vehicle go into the gap between the two obstacles,
%%and if the gap is two small or even negative, then the QP problem will
%%become non solvable. 
%in this case, we can modify the center of the obstacle a little, so that
%the solver will not let the vehicle go to the gap between the two
%obstacles, but go to another side which is not between the
%two obstacles. 
% pos_ob_array_pre(:,1) = [200+0*t; 0.0];
% pos_ob_array_pre(:,2) = [500; -0.0];
% pos_ob_array_pre(:,3) = [800; -0.1];
% pos_ob_array_pre(:,4) = [800;  1.5];
% pos_ob_array_pre(:,5) = [1000;  0.1];
% pos_ob_array_pre(:,6) = [1200; -0.1];
% pos_ob_array_pre(:,7) = [1400;  0.1];
% pos_ob_array_pre(:,8) = [1400; -0.5];
% pos_ob_array_pre(:,9) = [1700; -0.1];
% pos_ob_array_pre(:,10) = [1900; -0.1];

% pos_ob = [46.4496317612391,49.4606128179940,50.1520441588954,56.2263831741385,58.8805838768333;0.414589645301348,2.00584160428506,1.55964509106729,1.77552917183758,2.57119129022772];
% % pos_ob(1,:) = sort(pos_ob(1,:));
%     pos_ob_array_pre = pos_ob;
% vel_ob_array_pre(:,1) = [ 0; 0];
% vel_ob_array_pre(:,2) = [0; 0];


%the size of the output depends on the number of the obstacles and the
%number of the constraints:
 
%select the active obstacles, which are in the front of the vehicle, and the
%distance is less than a shreshold. 
dis_shresh = 600;
pos_ob_array = [];
vel_ob_array =[];
radius_array = [];
for i = 1:no_ob_input 
    if (pos_ob_array_pre_input(1,i)>=(s-3)) && (abs(s-pos_ob_array_pre_input(1,i))<=dis_shresh)  %be careful, should be very careful
        pos_ob_array = [pos_ob_array, pos_ob_array_pre_input(:,i)];
        vel_ob_array = [vel_ob_array, vel_ob_array_pre(:,i)];
        radius_array = [radius_array; radius_pre(i)];
    end
end
no_ob_active = size(pos_ob_array,2); 
if(no_ob_active == 0)
    pos_ob_array = pos_ob_array_pre_input(:,end);
    vel_ob_array = vel_ob_array_pre(:,end); 
    radius_array = radius_pre(end);%very important, cause a lot of errors 
    no_ob_active=1;
end

dis_2_vehicle = zeros(no_ob_active,1);
for i_ob =1:no_ob_active
%% states of the vehicle and each obstacle:
    pos_ob = pos_ob_array(:,i_ob);  
    vel_ob = vel_ob_array(:,i_ob);  
    acc_ob = acc_ob_array(:,i_ob);  
    
    pos_ob_x = pos_ob(1);
    pos_ob_y = pos_ob(2);
    vel_ob_x = vel_ob(1);
    vel_ob_y = vel_ob(2);
    acc_ob_x = acc_ob(1);
    acc_ob_y = acc_ob(2);

    v_vehicle = [xp_dot*cos(epsi)-yp_dot*sin(epsi); yp_dot*cos(epsi) + xp_dot*sin(epsi)];
    rel_pos = [s; ey] - pos_ob;
    rel_vel = v_vehicle - [vel_ob_x; vel_ob_y]; %assume epsi is small, 
    dis_2_vehicle(i_ob) = -norm(rel_pos)^2/(rel_pos'*rel_vel);  %distance from the vehicle
    dis_2_vehicle(i_ob) = norm(rel_pos) ;  %distance from the vehicle
end 

%order the obstacles, from nearest to farest, considering the velocity: 
[dis, i_dis] =sort(dis_2_vehicle);
pos_ob_array = pos_ob_array(:, i_dis);
vel_ob_array = vel_ob_array(:, i_dis); 
radius_array = radius_array( i_dis);


if (no_ob_active >=2)
    %if the x-coordinates of the first two obstacles are the same, then
    %modified the second obstacles a little, see the comments previously 
    if(pos_ob_array(1,1) == pos_ob_array(1,2))
       if(abs(pos_ob_array(2,1)) >= abs(pos_ob_array(2,2)))
            pos_ob_array(1,2) = pos_ob_array(1,2) - 0.01;
       else
           pos_ob_array(1,2) = pos_ob_array(1,2) + 0.01;
       end
    end
end

dead = 0; %unlock 

%On 20180815, I found a mistake in the dynamics
%%  for the road side constraints: 
ck = 1;  %for relax cbf 
ey_pos = 3.2;
ey_neg = -3.2; 

a_m = 4;  %maximum acc 

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

h_sid_pos = ey_pos - ey - (yp_dot*cos(epsi) + xp_dot*sin(epsi))^2/(2*a_m);
% L_f_h_sid_pos = ((yp_dot*cos(epsi) + xp_dot*sin(epsi))*(2*cf*yp_dot*cos(epsi) + 2*cr*yp_dot*cos(epsi) - a_m*m*xp_dot + 2*a*cf*psi_dot*cos(epsi) - 2*b*cr*psi_dot*cos(epsi) + m*psi_dot_com*xp_dot^2*cos(epsi) - m*psi_dot_com*xp_dot*yp_dot*sin(epsi)))/(a_m*m*xp_dot); 
% L_g_h_sid_pos = [ -(2*cf*cos(epsi)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/(a_m*m), -(sin(epsi)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/a_m];

L_f_h_sid_pos = ((yp_dot*cos(epsi) + xp_dot*sin(epsi))*(2*cf*yp_dot*cos(epsi) + 2*cr*yp_dot*cos(epsi) - a_m*m*xp_dot + 2*a*cf*psi_dot*cos(epsi) - 2*b*cr*psi_dot*cos(epsi) + m*psi_dot_com*xp_dot^2*cos(epsi) - m*psi_dot_com*xp_dot*yp_dot*sin(epsi)))/(a_m*m*xp_dot);
L_g_h_sid_pos = [ -(2*cf*cos(epsi)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/(a_m*m), -(sin(epsi)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/a_m];

h_sid_neg = ey - ey_neg - (yp_dot*cos(epsi) + xp_dot*sin(epsi))^2/(2*a_m);
% L_f_h_sid_neg = ((yp_dot*cos(epsi) + xp_dot*sin(epsi))*(2*cf*yp_dot*cos(epsi) + 2*cr*yp_dot*cos(epsi) + a_m*m*xp_dot + 2*a*cf*psi_dot*cos(epsi) - 2*b*cr*psi_dot*cos(epsi) + m*psi_dot_com*xp_dot^2*cos(epsi) - m*psi_dot_com*xp_dot*yp_dot*sin(epsi)))/(a_m*m*xp_dot);
% L_g_h_sid_neg = [ -(2*cf*cos(epsi)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/(a_m*m), -(sin(epsi)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/a_m];

L_f_h_sid_neg = ((yp_dot*cos(epsi) + xp_dot*sin(epsi))*(2*cf*yp_dot*cos(epsi) + 2*cr*yp_dot*cos(epsi) + a_m*m*xp_dot + 2*a*cf*psi_dot*cos(epsi) - 2*b*cr*psi_dot*cos(epsi) + m*psi_dot_com*xp_dot^2*cos(epsi) - m*psi_dot_com*xp_dot*yp_dot*sin(epsi)))/(a_m*m*xp_dot);
L_g_h_sid_neg = [ -(2*cf*cos(epsi)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/(a_m*m), -(sin(epsi)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/a_m];

A_n_side_pos =  -L_g_h_sid_pos;
b_n_side_pos = L_f_h_sid_pos + 0.5*h_sid_pos ;
A_n_side_neg=  [-L_g_h_sid_neg];
b_n_side_neg = L_f_h_sid_neg + 0.5*h_sid_neg ;

if(h_sid_pos>0) && ((yp_dot*cos(epsi) + xp_dot*sin(epsi))<0) %if y-velocity is negative, do not need to constraint 
    A_n_side_pos =  [0, 0];
    b_n_side_pos = 1;
end
if(h_sid_neg>0) && ((yp_dot*cos(epsi) + xp_dot*sin(epsi))>0) %if y-velocity is positive, do not need to constraint 
    A_n_side_neg =  [0, 0];
    b_n_side_neg = 1;
end
if (ey >= 3.7) || (ey  <= -3.7)
    dead = 1;  %kill all 
end

% if (h_sid_neg <=0) || (h_sid_pos  <= 0)
%     dead = 1;  %kill all 
% end

% tuning: 
% A_n_side_pos = [-0,  0];
% b_n_side_pos = 1;
% A_n_side_neg=  [0, 0];
% b_n_side_neg = 1;
% dead = 0;   



for i_ob =1:no_ob_active
%% states of the vehicle and each obstacle:
    pos_ob = pos_ob_array(:,i_ob);  
    vel_ob = vel_ob_array(:,i_ob);  
    acc_ob = acc_ob_array(:,i_ob);  
    
    pos_ob_x = pos_ob(1);
    pos_ob_y = pos_ob(2);
    vel_ob_x = vel_ob(1);
    vel_ob_y = vel_ob(2);
    acc_ob_x = acc_ob(1);
    acc_ob_y = acc_ob(2);
 
    v_vehicle = [xp_dot*cos(epsi)-yp_dot*sin(epsi); yp_dot*cos(epsi) + xp_dot*sin(epsi)];
    rel_pos = [s; ey] - pos_ob;
    rel_vel = v_vehicle - [vel_ob_x; vel_ob_y]; %assume epsi is small, 
    
%     Ds = 1.2;  %the radius of obstacle 
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
        
 %% August 3th, angle constraint for obstacles, 
        %constraint 1, angle constraint, test:  
    %notice there may be virtrual number, nan, of inf, you should avoid these
    %conditions: 
%     Ds = 1.1;
    Ds = radius_array(i_ob)+0.5; %0.2
    cos_rel_ang = (-rel_pos.'*rel_vel) /norm(rel_pos)/norm(rel_vel); 
    
    if (cos_rel_ang>=-0.99) 
        if (cos_rel_ang>=1)
            cos_rel_ang = 1;
        end
        rel_ang = acos(cos_rel_ang); 
        % h_ang = rel_ang-asin(Ds/norm(rel_pos));

        ratio = Ds/norm(rel_pos);
        %notice if do not deal carefully, there maybe virtual number appearing. 
        if ratio>=1
            ratio = 0.9999999;
            dead = 1;  %kill all 
        end    
        h_ang = rel_ang - asin(ratio); 

        %very important, if this is zero, may make the qp infeasible due to the
        %zeros coefficient matrix, but the right matrix is negative 
        %sometimes, it may be NaN of inf, so this constraint actually do not work
        %noticed on July, 18th, 2018
        if(abs(pos_ob_y)<1e-4)    
            pos_ob_y = sign(pos_ob_y)*1e-4;
        end
        if (abs(epsi) <= 1e-5)
            epsi = sign(epsi)*1e-5;
        end
        
%         if (pos_ob_y == 0)
%             pos_ob_y =   1e-4;
%         end
%         if (epsi == 0 )
%             epsi =   1e-5;
%         end

        %very important, due to calculation errors, there maybe sometimes the
        %results are virtual, should be treated carefully. 

        %notice, the longitudinal velocity is not be considered as state here 
        % L_f_h_ang = (yp_dot*cos(epsi) + xp_dot*sin(epsi))*((Ds*(ey - pos_ob_y))/((1 - Ds^2/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - ((pos_ob_x - s)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2))) - cos(epsi)*(xp_dot - yp_dot)*((Ds*(pos_ob_x - s))/((1 - Ds^2/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) + ((ey - pos_ob_y)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2))) + ((psi_dot - psi_dot_com)*(xp_dot*yp_dot - xp_dot^2 + vel_ob_x*xp_dot*cos(epsi) + vel_ob_y*xp_dot*sin(epsi) - vel_ob_x*yp_dot*sin(epsi) - vel_ob_y*yp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)) - (cos(epsi)*(vel_ob_x + vel_ob_y - xp_dot*cos(epsi) - xp_dot*sin(epsi))*(m*psi_dot*xp_dot^2 + 2*cf*yp_dot + 2*cr*yp_dot + 2*a*cf*psi_dot - 2*b*cr*psi_dot)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*xp_dot*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)); 
        % L_g_h_ang = (2*cf*cos(epsi)*(vel_ob_x + vel_ob_y - xp_dot*cos(epsi) - xp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)); 

%         L_f_h_ang_part1 = ((psi_dot - psi_dot_com)*(xp_dot*yp_dot - xp_dot^2 + vel_ob_x*xp_dot*cos(epsi) + vel_ob_y*xp_dot*sin(epsi) - vel_ob_x*yp_dot*sin(epsi) - vel_ob_y*yp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)) - ((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)) - (cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)) - (cos(epsi)*(vel_ob_x + vel_ob_y - xp_dot*cos(epsi) - xp_dot*sin(epsi))*(m*psi_dot*xp_dot^2 + 2*cf*yp_dot + 2*cr*yp_dot + 2*a*cf*psi_dot - 2*b*cr*psi_dot)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*xp_dot*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)); 
%         L_t_h_ang_part1 = -((ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi))*(pos_ob_y*vel_ob_x^3 - pos_ob_x*vel_ob_y^3 - ey*vel_ob_x^3 + s*vel_ob_y^3 - acc_ob_x*pos_ob_x^2*vel_ob_y + acc_ob_y*pos_ob_x^2*vel_ob_x - acc_ob_x*pos_ob_y^2*vel_ob_y + acc_ob_y*pos_ob_y^2*vel_ob_x - acc_ob_x*s^2*vel_ob_y + acc_ob_y*s^2*vel_ob_x - ey*vel_ob_x*vel_ob_y^2 - ey*vel_ob_x*xp_dot^2 - pos_ob_x*vel_ob_x^2*vel_ob_y + pos_ob_y*vel_ob_x*vel_ob_y^2 - pos_ob_x*vel_ob_y*xp_dot^2 + pos_ob_y*vel_ob_x*xp_dot^2 + s*vel_ob_x^2*vel_ob_y + s*vel_ob_y*xp_dot^2 - acc_ob_x*ey^2*vel_ob_y + acc_ob_y*ey^2*vel_ob_x + 2*acc_ob_x*ey*pos_ob_y*vel_ob_y - 2*acc_ob_y*ey*pos_ob_y*vel_ob_x + 2*acc_ob_x*pos_ob_x*s*vel_ob_y - 2*acc_ob_y*pos_ob_x*s*vel_ob_x - 2*ey*vel_ob_x*yp_dot^2*cos(epsi)^2 - 2*pos_ob_x*vel_ob_y*yp_dot^2*cos(epsi)^2 + 2*pos_ob_y*vel_ob_x*yp_dot^2*cos(epsi)^2 + 2*s*vel_ob_y*yp_dot^2*cos(epsi)^2 - acc_ob_y*ey^2*xp_dot*cos(epsi) + acc_ob_x*ey^2*yp_dot*cos(epsi) + acc_ob_y*ey^2*yp_dot*cos(epsi) - acc_ob_y*pos_ob_x^2*xp_dot*cos(epsi) - acc_ob_y*pos_ob_y^2*xp_dot*cos(epsi) + acc_ob_x*pos_ob_x^2*yp_dot*cos(epsi) + acc_ob_x*pos_ob_y^2*yp_dot*cos(epsi) + acc_ob_y*pos_ob_x^2*yp_dot*cos(epsi) + acc_ob_y*pos_ob_y^2*yp_dot*cos(epsi) - acc_ob_y*s^2*xp_dot*cos(epsi) + acc_ob_x*s^2*yp_dot*cos(epsi) + acc_ob_y*s^2*yp_dot*cos(epsi) + acc_ob_x*ey^2*xp_dot*sin(epsi) + 2*ey*vel_ob_x^2*xp_dot*cos(epsi) - 2*ey*vel_ob_x^2*yp_dot*cos(epsi) + acc_ob_x*pos_ob_x^2*xp_dot*sin(epsi) + acc_ob_x*pos_ob_y^2*xp_dot*sin(epsi) + acc_ob_x*s^2*xp_dot*sin(epsi) - 2*pos_ob_y*vel_ob_x^2*xp_dot*cos(epsi) + 2*pos_ob_x*vel_ob_y^2*yp_dot*cos(epsi) + 2*pos_ob_y*vel_ob_x^2*yp_dot*cos(epsi) - 2*s*vel_ob_y^2*yp_dot*cos(epsi) + 2*pos_ob_x*vel_ob_y^2*xp_dot*sin(epsi) - 2*s*vel_ob_y^2*xp_dot*sin(epsi) - 2*s*vel_ob_x*vel_ob_y*xp_dot*cos(epsi) + 2*s*vel_ob_x*vel_ob_y*yp_dot*cos(epsi) + 2*ey*vel_ob_x*vel_ob_y*xp_dot*sin(epsi) - 2*pos_ob_y*vel_ob_x*vel_ob_y*xp_dot*sin(epsi) + 2*ey*vel_ob_x*xp_dot*yp_dot*cos(epsi)^2 + 2*pos_ob_x*vel_ob_y*xp_dot*yp_dot*cos(epsi)^2 - 2*pos_ob_y*vel_ob_x*xp_dot*yp_dot*cos(epsi)^2 - 2*s*vel_ob_y*xp_dot*yp_dot*cos(epsi)^2 - ey*vel_ob_x*xp_dot*yp_dot*sin(2*epsi) - pos_ob_x*vel_ob_y*xp_dot*yp_dot*sin(2*epsi) + pos_ob_y*vel_ob_x*xp_dot*yp_dot*sin(2*epsi) + s*vel_ob_y*xp_dot*yp_dot*sin(2*epsi) + 2*acc_ob_y*ey*pos_ob_y*xp_dot*cos(epsi) - 2*acc_ob_x*ey*pos_ob_y*yp_dot*cos(epsi) - 2*acc_ob_y*ey*pos_ob_y*yp_dot*cos(epsi) + 2*acc_ob_y*pos_ob_x*s*xp_dot*cos(epsi) - 2*acc_ob_x*pos_ob_x*s*yp_dot*cos(epsi) - 2*acc_ob_y*pos_ob_x*s*yp_dot*cos(epsi) - 2*acc_ob_x*ey*pos_ob_y*xp_dot*sin(epsi) + 2*ey*vel_ob_x*vel_ob_y*yp_dot*cos(epsi) - 2*acc_ob_x*pos_ob_x*s*xp_dot*sin(epsi) + 2*pos_ob_x*vel_ob_x*vel_ob_y*xp_dot*cos(epsi) - 2*pos_ob_x*vel_ob_x*vel_ob_y*yp_dot*cos(epsi) - 2*pos_ob_y*vel_ob_x*vel_ob_y*yp_dot*cos(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2));
%         L_g_h_ang_part1 = (2*cf*cos(epsi)*(vel_ob_x + vel_ob_y - xp_dot*cos(epsi) - xp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)); 
       
        L_f_h_ang_part1 = - ((psi_dot - psi_dot_com)*(xp_dot^2 + yp_dot^2 - vel_ob_x*xp_dot*cos(epsi) - vel_ob_y*yp_dot*cos(epsi) - vel_ob_y*xp_dot*sin(epsi) + vel_ob_x*yp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/((1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)) - ((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/((1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)) - ((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/((1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)) - ((vel_ob_x*cos(epsi) - xp_dot + vel_ob_y*sin(epsi))*(m*psi_dot*xp_dot^2 + 2*cf*yp_dot + 2*cr*yp_dot + 2*a*cf*psi_dot - 2*b*cr*psi_dot)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*xp_dot*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2));
        L_t_h_ang_part1 = -((ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi))*(pos_ob_y*vel_ob_x^3 - pos_ob_x*vel_ob_y^3 - ey*vel_ob_x^3 + s*vel_ob_y^3 - acc_ob_x*pos_ob_x^2*vel_ob_y + acc_ob_y*pos_ob_x^2*vel_ob_x - acc_ob_x*pos_ob_y^2*vel_ob_y + acc_ob_y*pos_ob_y^2*vel_ob_x - acc_ob_x*s^2*vel_ob_y + acc_ob_y*s^2*vel_ob_x - ey*vel_ob_x*vel_ob_y^2 - ey*vel_ob_x*xp_dot^2 - ey*vel_ob_x*yp_dot^2 - pos_ob_x*vel_ob_x^2*vel_ob_y + pos_ob_y*vel_ob_x*vel_ob_y^2 - pos_ob_x*vel_ob_y*xp_dot^2 + pos_ob_y*vel_ob_x*xp_dot^2 + s*vel_ob_x^2*vel_ob_y - pos_ob_x*vel_ob_y*yp_dot^2 + pos_ob_y*vel_ob_x*yp_dot^2 + s*vel_ob_y*xp_dot^2 + s*vel_ob_y*yp_dot^2 - acc_ob_x*ey^2*vel_ob_y + acc_ob_y*ey^2*vel_ob_x + 2*acc_ob_x*ey*pos_ob_y*vel_ob_y - 2*acc_ob_y*ey*pos_ob_y*vel_ob_x + 2*acc_ob_x*pos_ob_x*s*vel_ob_y - 2*acc_ob_y*pos_ob_x*s*vel_ob_x - acc_ob_y*ey^2*xp_dot*cos(epsi) + acc_ob_x*ey^2*yp_dot*cos(epsi) - acc_ob_y*pos_ob_x^2*xp_dot*cos(epsi) - acc_ob_y*pos_ob_y^2*xp_dot*cos(epsi) + acc_ob_x*pos_ob_x^2*yp_dot*cos(epsi) + acc_ob_x*pos_ob_y^2*yp_dot*cos(epsi) - acc_ob_y*s^2*xp_dot*cos(epsi) + acc_ob_x*s^2*yp_dot*cos(epsi) + acc_ob_x*ey^2*xp_dot*sin(epsi) + acc_ob_y*ey^2*yp_dot*sin(epsi) + 2*ey*vel_ob_x^2*xp_dot*cos(epsi) + acc_ob_x*pos_ob_x^2*xp_dot*sin(epsi) + acc_ob_x*pos_ob_y^2*xp_dot*sin(epsi) + acc_ob_y*pos_ob_x^2*yp_dot*sin(epsi) + acc_ob_y*pos_ob_y^2*yp_dot*sin(epsi) + acc_ob_x*s^2*xp_dot*sin(epsi) + acc_ob_y*s^2*yp_dot*sin(epsi) - 2*pos_ob_y*vel_ob_x^2*xp_dot*cos(epsi) + 2*pos_ob_x*vel_ob_y^2*yp_dot*cos(epsi) - 2*s*vel_ob_y^2*yp_dot*cos(epsi) - 2*ey*vel_ob_x^2*yp_dot*sin(epsi) + 2*pos_ob_x*vel_ob_y^2*xp_dot*sin(epsi) + 2*pos_ob_y*vel_ob_x^2*yp_dot*sin(epsi) - 2*s*vel_ob_y^2*xp_dot*sin(epsi) - 2*s*vel_ob_x*vel_ob_y*xp_dot*cos(epsi) + 2*ey*vel_ob_x*vel_ob_y*xp_dot*sin(epsi) - 2*pos_ob_y*vel_ob_x*vel_ob_y*xp_dot*sin(epsi) - 2*pos_ob_x*vel_ob_x*vel_ob_y*yp_dot*sin(epsi) + 2*s*vel_ob_x*vel_ob_y*yp_dot*sin(epsi) + 2*acc_ob_y*ey*pos_ob_y*xp_dot*cos(epsi) - 2*acc_ob_x*ey*pos_ob_y*yp_dot*cos(epsi) + 2*acc_ob_y*pos_ob_x*s*xp_dot*cos(epsi) - 2*acc_ob_x*pos_ob_x*s*yp_dot*cos(epsi) - 2*acc_ob_x*ey*pos_ob_y*xp_dot*sin(epsi) - 2*acc_ob_y*ey*pos_ob_y*yp_dot*sin(epsi) + 2*ey*vel_ob_x*vel_ob_y*yp_dot*cos(epsi) - 2*acc_ob_x*pos_ob_x*s*xp_dot*sin(epsi) - 2*acc_ob_y*pos_ob_x*s*yp_dot*sin(epsi) + 2*pos_ob_x*vel_ob_x*vel_ob_y*xp_dot*cos(epsi) - 2*pos_ob_y*vel_ob_x*vel_ob_y*yp_dot*cos(epsi)))/((1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2));
        L_g_h_ang_part1 = (2*cf*(vel_ob_x*cos(epsi) - xp_dot + vel_ob_y*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) + ey*yp_dot*sin(epsi) - pos_ob_x*xp_dot*sin(epsi) - pos_ob_y*yp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2));
 
%         L_f_h_ang_part2 = ((ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - ((xp_dot*cos(epsi) - yp_dot*cos(epsi))*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
%         L_g_h_ang_part2 = 0;
%         L_t_h_ang_part2 = -(ey*vel_ob_y - pos_ob_x*vel_ob_x - pos_ob_y*vel_ob_y + s*vel_ob_x)/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);

        L_f_h_ang_part2 = ((ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - ((pos_ob_x - s)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
        L_g_h_ang_part2 = 0;
        L_t_h_ang_part2 = -(ey*vel_ob_y - pos_ob_x*vel_ob_x - pos_ob_y*vel_ob_y + s*vel_ob_x)/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
 
        
        ratio_pf = -Ds/ (norm(rel_pos))^2 *L_f_h_ang_part2; 
        ratio_pg = -Ds/ (norm(rel_pos))^2 *L_g_h_ang_part2; 
        ratio_pt = -Ds/ (norm(rel_pos))^2 *L_t_h_ang_part2; 
        asin_dot = 1/(1 - ratio^2)^(1/2); 
        L_f_h_ang=  L_f_h_ang_part1   -asin_dot*ratio_pf;
        L_g_h_ang=  L_g_h_ang_part1   -asin_dot*ratio_pg;
        L_t_h_ang=  L_t_h_ang_part1   -asin_dot*ratio_pt;

        A_n_angle_fix = [-L_g_h_ang, 0];
%         b_n_angle_fix = L_f_h_ang + L_t_h_ang + 20*h_ang^7;
        b_n_angle_fix = L_f_h_ang + L_t_h_ang + 1*h_ang^5;
        
        %20180815 tuning parameters and k-function
        b_n_angle_fix = L_f_h_ang + L_t_h_ang + 1*h_ang^5;
        b_n_angle_fix = L_f_h_ang + L_t_h_ang + 1*h_ang;
        b_n_angle_fix = L_f_h_ang + L_t_h_ang + 3*h_ang;
%         b_n_angle_fix = L_f_h_ang + L_t_h_ang + 3*h_ang^5;
    elseif (cos_rel_ang<=-0.99)
        h_ang = 3; 
        A_n_angle_fix = [0, 0];
        b_n_angle_fix = 1; 
    end
 
%%  August 3rd, constraint 2, velocity constraint, test:  
    %notice there may be virtrual number, nan, of inf, you should avoid these
    %conditions: 
%     Ds = 1;    
    Ds = radius_array(i_ob)+0.5;
%     am = 1; 
    alert = 0 ;
    dis_maxacc_sqr = 2*a_m*(norm(rel_pos)-Ds); 
    if (dis_maxacc_sqr<=0)
        alert = 1;
        dead = 1;  %kill all 
        dis_maxacc_sqr = 1e-6;
    end
    h_vel  = sqrt(dis_maxacc_sqr) + rel_pos'/norm(rel_pos)*rel_vel; 
    
	if (h_vel<= 0.1) && (h_ang<=-0.01)
%         alert = 1; 
    end

    %very important, if this is zero, may make the qp infeasible due to the
    %zeros 
    %coefficient matrix, but the right matrix is negative 
    %sometimes, it may be NaN of inf, so this constraint actually do not work
    %noticed on July, 18th, 2018
%     if(abs(pos_ob_y)<1e-4)    
%         pos_ob_y = sign(pos_ob_y)*1e-4;
%     end
%     if (abs(epsi) <= 1e-5)
%         epsi = sign(epsi)*1e-5;
%     end

    % L_f_h_ang = (yp_dot*cos(epsi) + xp_dot*sin(epsi))*((Ds*(ey - pos_ob_y))/((1 - Ds^2/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - ((pos_ob_x - s)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2))) - cos(epsi)*(xp_dot - yp_dot)*((Ds*(pos_ob_x - s))/((1 - Ds^2/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) + ((ey - pos_ob_y)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2))) + ((psi_dot - psi_dot_com)*(xp_dot*yp_dot - xp_dot^2 + vel_ob_x*xp_dot*cos(epsi) + vel_ob_y*xp_dot*sin(epsi) - vel_ob_x*yp_dot*sin(epsi) - vel_ob_y*yp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)) - (cos(epsi)*(vel_ob_x + vel_ob_y - xp_dot*cos(epsi) - xp_dot*sin(epsi))*(m*psi_dot*xp_dot^2 + 2*cf*yp_dot + 2*cr*yp_dot + 2*a*cf*psi_dot - 2*b*cr*psi_dot)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*xp_dot*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)); 
    % L_g_h_ang = (2*cf*cos(epsi)*(vel_ob_x + vel_ob_y - xp_dot*cos(epsi) - xp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)); 

%     L_f_h_norm = ((ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - ((xp_dot*cos(epsi) - yp_dot*cos(epsi))*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
%     L_g_h_norm = [ 0, 0];
%     L_t_h_norm = -(ey*vel_ob_y - pos_ob_x*vel_ob_x - pos_ob_y*vel_ob_y + s*vel_ob_x)/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);

    L_f_h_norm = ((ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - ((pos_ob_x - s)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
    L_g_h_norm = [ 0, 0];
    L_t_h_norm = -(ey*vel_ob_y - pos_ob_x*vel_ob_x - pos_ob_y*vel_ob_y + s*vel_ob_x)/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2); 
    
    
%     L_f_h_vel = (((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((pos_ob_x - s)*(xp_dot*sin(epsi) - yp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2))*(psi_dot - psi_dot_com) + (xp_dot*cos(epsi) - yp_dot*cos(epsi))*(((2*pos_ob_x - 2*s)*(ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - (vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((2*pos_ob_x - 2*s)*(pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2))) - (yp_dot*cos(epsi) + xp_dot*sin(epsi))*(((ey - pos_ob_y)*(2*ey - 2*pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((pos_ob_x - s)*(2*ey - 2*pos_ob_y)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2))) - psi_dot*yp_dot*((cos(epsi)*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - (sin(epsi)*(ey - pos_ob_y))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)) - (cos(epsi)*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot))*(ey + pos_ob_x - pos_ob_y - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
%     L_g_h_vel = [ (2*cf*cos(epsi)*(ey + pos_ob_x - pos_ob_y - s))/(m*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)), (sin(epsi)*(ey - pos_ob_y))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - (cos(epsi)*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)];
%     L_t_h_vel = vel_ob_y*(((ey - pos_ob_y)*(2*ey - 2*pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((pos_ob_x - s)*(2*ey - 2*pos_ob_y)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2))) - vel_ob_x*(((2*pos_ob_x - 2*s)*(ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - (vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((2*pos_ob_x - 2*s)*(pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2))) - (acc_ob_y*(ey - pos_ob_y))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + (acc_ob_x*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);

    L_f_h_vel = (((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2))*(psi_dot - psi_dot_com) - (yp_dot*cos(epsi) + xp_dot*sin(epsi))*(((ey - pos_ob_y)*(2*ey - 2*pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((pos_ob_x - s)*(2*ey - 2*pos_ob_y)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2))) + (xp_dot*cos(epsi) - yp_dot*sin(epsi))*(((2*pos_ob_x - 2*s)*(ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - (vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((2*pos_ob_x - 2*s)*(pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2))) - ((cos(epsi)*(ey - pos_ob_y))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + (sin(epsi)*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2))*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - psi_dot*yp_dot*((cos(epsi)*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - (sin(epsi)*(ey - pos_ob_y))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2));
    L_g_h_vel = [ (2*cf*((cos(epsi)*(ey - pos_ob_y))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + (sin(epsi)*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)))/m, (sin(epsi)*(ey - pos_ob_y))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - (cos(epsi)*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)];
    L_t_h_vel = vel_ob_y*(((ey - pos_ob_y)*(2*ey - 2*pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((pos_ob_x - s)*(2*ey - 2*pos_ob_y)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2))) - vel_ob_x*(((2*pos_ob_x - 2*s)*(ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - (vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + ((2*pos_ob_x - 2*s)*(pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*sin(epsi)))/(2*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2))) - (acc_ob_y*(ey - pos_ob_y))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) + (acc_ob_x*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
    
    temp_h_vel = 1/2*1/sqrt(dis_maxacc_sqr)*2*a_m; 

    L_f_h_vel= temp_h_vel* L_f_h_norm+L_f_h_vel;
    L_g_h_vel= temp_h_vel* L_g_h_norm+L_g_h_vel;
    L_t_h_vel= temp_h_vel* L_t_h_norm+L_t_h_vel;

    A_n_vel = -L_g_h_vel;
    b_n_vel = L_f_h_vel + L_t_h_vel + 2*h_vel^5; %seems it is safe it the gain is large? 
    
    %20180815 tuning parameters and k-function
    b_n_vel = L_f_h_vel + L_t_h_vel + 2*h_vel; %seems it is safe it the gain is large? 
    b_n_vel = L_f_h_vel + L_t_h_vel + 2*h_vel^5; %seems it is safe it the gain is large? 
    b_n_vel = L_f_h_vel + L_t_h_vel + 3*h_vel; 
    
     %%  August 10th, constraint 2, velocity constraint, 2nd type, only constraint the acc, tested not good, always induces to infeasible solution
    %notice there may be virtrual number, nan, of inf, you should avoid these
    %conditions: 
%     Ds = 1; 
%     Ds = radius_array(i_ob)+0.1;
%     alert = 0;
%     dis_maxacc_sqr = 2*a_m*(norm(rel_pos)-Ds); 
%     if (dis_maxacc_sqr<= 1e-4)
%         dis_maxacc_sqr = 1e-4;
%         alert = 1; 
%     end
%     h_vel  = dis_maxacc_sqr - xp_dot^2; 
%     
%     if (h_vel<= 0.1)
% %         dis_maxacc_sqr = 1e-4;
% %         alert = 1; 
%     end
%     
%  
%     %very important, if this is zero, may make the qp infeasible due to the
%     %coefficient matrix be zeros, but the right matrix is negative 
%     %sometimes, it may be NaN of inf, so this constraint actually do not work
%     %noticed on July, 18th, 2018
%     if(abs(pos_ob_y)<1e-4)    
%         pos_ob_y = sign(pos_ob_y)*1e-4;
%     end
%     if (abs(epsi) <= 1e-5)
%         epsi = sign(epsi)*1e-5;
%     end
% 
%     % L_f_h_ang = (yp_dot*cos(epsi) + xp_dot*sin(epsi))*((Ds*(ey - pos_ob_y))/((1 - Ds^2/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) - ((pos_ob_x - s)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2))) - cos(epsi)*(xp_dot - yp_dot)*((Ds*(pos_ob_x - s))/((1 - Ds^2/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2))^(1/2)*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)) + ((ey - pos_ob_y)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(3/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(1/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2))) + ((psi_dot - psi_dot_com)*(xp_dot*yp_dot - xp_dot^2 + vel_ob_x*xp_dot*cos(epsi) + vel_ob_y*xp_dot*sin(epsi) - vel_ob_x*yp_dot*sin(epsi) - vel_ob_y*yp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)) - (cos(epsi)*(vel_ob_x + vel_ob_y - xp_dot*cos(epsi) - xp_dot*sin(epsi))*(m*psi_dot*xp_dot^2 + 2*cf*yp_dot + 2*cr*yp_dot + 2*a*cf*psi_dot - 2*b*cr*psi_dot)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*xp_dot*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)); 
%     % L_g_h_ang = (2*cf*cos(epsi)*(vel_ob_x + vel_ob_y - xp_dot*cos(epsi) - xp_dot*sin(epsi))*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y - ey*xp_dot*cos(epsi) + ey*yp_dot*cos(epsi) + pos_ob_y*xp_dot*cos(epsi) - pos_ob_x*yp_dot*cos(epsi) - pos_ob_y*yp_dot*cos(epsi) + s*yp_dot*cos(epsi) - pos_ob_x*xp_dot*sin(epsi) + s*xp_dot*sin(epsi)))/(m*((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)^(3/2)*(1 - ((pos_ob_x - s)*(vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi)) + (ey - pos_ob_y)*(yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi)))^2/(((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)*((vel_ob_x - xp_dot*cos(epsi) + yp_dot*cos(epsi))^2 + (yp_dot*cos(epsi) - vel_ob_y + xp_dot*sin(epsi))^2)))^(1/2)); 
% 
%     L_f_h_norm = ((ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2) - ((xp_dot*cos(epsi) - yp_dot*cos(epsi))*(pos_ob_x - s))/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
%     L_g_h_norm = [ 0, 0];
%     L_t_h_norm = -(ey*vel_ob_y - pos_ob_x*vel_ob_x - pos_ob_y*vel_ob_y + s*vel_ob_x)/((pos_ob_x - s)^2 + (ey - pos_ob_y)^2)^(1/2);
% 
%     L_f_h_vel = -2*xp_dot*yp_dot*psi_dot;
%     L_g_h_vel = [0, -2*xp_dot];
%     L_t_h_vel =  0; 
%     
%     temp_h_vel =  2*a_m; 
% 
%     L_f_h_vel= temp_h_vel* L_f_h_norm+L_f_h_vel;
%     L_g_h_vel= temp_h_vel* L_g_h_norm+L_g_h_vel;
%     L_t_h_vel= temp_h_vel* L_t_h_norm+L_t_h_vel;
%   
%     A_n_vel = -L_g_h_vel;
%     b_n_vel = L_f_h_vel + L_t_h_vel + 100*h_vel^7;
    

%% August, 3rd, constraint 3, velocity direction constraint for moving obstacles, test:  

    %notice there may be virtrual number, nan, of inf, you should avoid these
    %conditions: 
   
    if (rel_pos(1)> 1) && (rel_pos'* v_vehicle* rel_pos'*vel_ob < 0) 
        h_move_angle = -2; 
        A_n_anglemoving = [0, 0];
        B_n_anglemoving = 1; 
    else
        rel_pos_vert =  [-rel_pos(2); rel_pos(1)];  %normal to rel_pos
        % h_move_angle = -rel_pos_vert'* v_vehicle* rel_pos_vert'*vel_ob/norm(rel_pos_vert)/norm(rel_pos_vert) ;
        h_move_angle = -rel_pos_vert'* v_vehicle* rel_pos_vert'*vel_ob  ;
        %very important, if this is zero, may make the qp infeasible due to the
        %coefficient matrix be zeros, but the right matrix is negative 
        %sometimes, it may be NaN of inf, so this constraint actually do not work
        %noticed on July, 18th, 2018
        % if(abs(pos_ob_y)<1e-4)    
        %     pos_ob_y = sign(pos_ob_y)*1e-4;
        % end
        % if (abs(epsi) <= 1e-5)
        %     epsi = sign(epsi)*1e-5;
        % end

%         L_f_h_move_angle = cos(epsi)*(xp_dot - yp_dot)*(vel_ob_y*((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot)) + vel_ob_x*(ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + vel_ob_y*(pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) - (yp_dot*cos(epsi) + xp_dot*sin(epsi))*(vel_ob_x*((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot)) + vel_ob_x*cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot) + vel_ob_y*cos(epsi)*(pos_ob_x - s)*(xp_dot - yp_dot)) - (vel_ob_y*((pos_ob_x - s)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - sin(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot))*(pos_ob_x - s) + vel_ob_x*((pos_ob_x - s)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) - sin(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot))*(ey - pos_ob_y))*(psi_dot - psi_dot_com) - psi_dot*yp_dot*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y)*(ey*cos(epsi) - pos_ob_y*cos(epsi) + pos_ob_x*sin(epsi) - s*sin(epsi)) - (cos(epsi)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y)*(ey - pos_ob_x - pos_ob_y + s)*(m*psi_dot*xp_dot^2 + 2*cf*yp_dot + 2*cr*yp_dot + 2*a*cf*psi_dot - 2*b*cr*psi_dot))/(m*xp_dot);
%         L_g_h_move_angle = [ (2*cf*cos(epsi)*(ey*vel_ob_x + pos_ob_x*vel_ob_y - pos_ob_y*vel_ob_x - s*vel_ob_y)*(ey - pos_ob_x - pos_ob_y + s))/m, - vel_ob_x*(ey - pos_ob_y)*(cos(epsi)*(ey - pos_ob_y) + sin(epsi)*(pos_ob_x - s)) - vel_ob_y*(pos_ob_x - s)*(cos(epsi)*(ey - pos_ob_y) + sin(epsi)*(pos_ob_x - s))];
%         L_t_h_move_angle = vel_ob_y*(vel_ob_x*((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot)) + vel_ob_x*cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot) + vel_ob_y*cos(epsi)*(pos_ob_x - s)*(xp_dot - yp_dot)) - vel_ob_x*(vel_ob_y*((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot)) + vel_ob_x*(ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + vel_ob_y*(pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) - acc_ob_x*((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot))*(ey - pos_ob_y) - acc_ob_y*((pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + cos(epsi)*(ey - pos_ob_y)*(xp_dot - yp_dot))*(pos_ob_x - s);
        
        L_f_h_move_angle = (xp_dot*cos(epsi) - yp_dot*sin(epsi))*(vel_ob_y*((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) + (pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) + vel_ob_x*(ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + vel_ob_y*(pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) - (yp_dot*cos(epsi) + xp_dot*sin(epsi))*(vel_ob_x*((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) + (pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) + vel_ob_x*(ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) + vel_ob_y*(pos_ob_x - s)*(xp_dot*cos(epsi) - yp_dot*sin(epsi))) + (psi_dot - psi_dot_com)*(vel_ob_x*(ey - pos_ob_y)*((ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) - (pos_ob_x - s)*(xp_dot*cos(epsi) - yp_dot*sin(epsi))) + vel_ob_y*(pos_ob_x - s)*((ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) - (pos_ob_x - s)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)))) + (vel_ob_x*(ey - pos_ob_y)*(cos(epsi)*(pos_ob_x - s) - sin(epsi)*(ey - pos_ob_y)) + vel_ob_y*(pos_ob_x - s)*(cos(epsi)*(pos_ob_x - s) - sin(epsi)*(ey - pos_ob_y)))*(psi_dot*xp_dot + (psi_dot*(2*a*cf - 2*b*cr))/(m*xp_dot) + (yp_dot*(2*cf + 2*cr))/(m*xp_dot)) - psi_dot*yp_dot*(vel_ob_x*(ey - pos_ob_y)*(cos(epsi)*(ey - pos_ob_y) + sin(epsi)*(pos_ob_x - s)) + vel_ob_y*(pos_ob_x - s)*(cos(epsi)*(ey - pos_ob_y) + sin(epsi)*(pos_ob_x - s)));
        L_g_h_move_angle = [ -(2*cf*(vel_ob_x*(ey - pos_ob_y)*(cos(epsi)*(pos_ob_x - s) - sin(epsi)*(ey - pos_ob_y)) + vel_ob_y*(pos_ob_x - s)*(cos(epsi)*(pos_ob_x - s) - sin(epsi)*(ey - pos_ob_y))))/m, - vel_ob_x*(ey - pos_ob_y)*(cos(epsi)*(ey - pos_ob_y) + sin(epsi)*(pos_ob_x - s)) - vel_ob_y*(pos_ob_x - s)*(cos(epsi)*(ey - pos_ob_y) + sin(epsi)*(pos_ob_x - s))];
        L_t_h_move_angle = vel_ob_y*(vel_ob_x*((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) + (pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) + vel_ob_x*(ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) + vel_ob_y*(pos_ob_x - s)*(xp_dot*cos(epsi) - yp_dot*sin(epsi))) - vel_ob_x*(vel_ob_y*((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) + (pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) + vel_ob_x*(ey - pos_ob_y)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)) + vel_ob_y*(pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) - acc_ob_x*(ey - pos_ob_y)*((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) + (pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi))) - acc_ob_y*(pos_ob_x - s)*((ey - pos_ob_y)*(xp_dot*cos(epsi) - yp_dot*sin(epsi)) + (pos_ob_x - s)*(yp_dot*cos(epsi) + xp_dot*sin(epsi)));  
 
        A_n_anglemoving =  -L_g_h_move_angle; 
        B_n_anglemoving = L_f_h_move_angle + L_t_h_move_angle +  3*h_move_angle;
        
        h_move_angle = -2; 
        A_n_anglemoving = [0, 0];
        B_n_anglemoving = 1; 
    end 
 
        
%% output, notice the define of the output variables: 
    out(i_ob).norm_relpos = norm(rel_pos);
    out(i_ob).h_angle_moving= h_move_angle; 
    out(i_ob).A_n_angle_moving  = A_n_anglemoving; 
    out(i_ob).B_n_angle_moving = B_n_anglemoving;    
    out(i_ob).h_angle_fix =  h_ang; 
    out(i_ob).A_n_angle_fix=  A_n_angle_fix; 
    out(i_ob).B_n_angle_fix = b_n_angle_fix;
    out(i_ob).h_dis = h_vel; 
    out(i_ob).A_n_dis = A_n_vel; 
    out(i_ob).B_n_dis = b_n_vel; 
    out(i_ob).alert = alert;
    out(i_ob).h_sid_pos = h_sid_pos;
    out(i_ob).A_n_side_pos = A_n_side_pos;
    out(i_ob).b_n_side_pos = b_n_side_pos;
    out(i_ob).h_sid_pos = h_sid_neg;
    out(i_ob).A_n_side_neg = A_n_side_neg;
    out(i_ob).b_n_side_neg = b_n_side_neg;
    out(i_ob).radius = radius_array(i_ob);
        
%     %test:
%     out(i_ob).h_angle_moving= NaN;  %cannot be 0?
%     out(i_ob).A_n_angle_moving  = A_n_anglemoving; 
%     out(i_ob).B_n_angle_moving = B_n_anglemoving;    
%     out(i_ob).h_angle_fix =  h_ang; 
%     out(i_ob).A_n_angle_fix=  A_n_angle_fix; 
%     out(i_ob).B_n_angle_fix = b_n_angle_fix;
%     out(i_ob).h_dis = h_vel; 
%     out(i_ob).A_n_dis = A_n_vel; 
%     out(i_ob).B_n_dis = b_n_vel; 
end




