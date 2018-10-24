%% simulator for obstacles: 
function traj = ob_traj(t)
 %calculate the trajectory of obstacles: 

global pos_ob_array_pre radius_pre;  
global no_ob;

for i_h=1:no_ob
    v= sqrt(t);
    v= 0; 
 

    %x,y notice the variables are in this order, velocity and acc are the same 
    %different from the orders in the trajectory generator for the vehicles
    traj.pos(:,i_h)  = pos_ob_array_pre(:,i_h)  + [ v*t ; 0];  
    traj.vel(:,i_h)  = [ v; 0];
    traj.acc(:,i_h)  = [ 0; 0];     
end
traj.rad = radius_pre;
traj.no = no_ob; 
end