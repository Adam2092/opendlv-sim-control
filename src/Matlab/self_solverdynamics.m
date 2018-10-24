function [t1,y1_nom, y1_actual, u1]=self_solverdynamics(ode_nom, tspan, y0, options, current_hdl_nom, ode_actual, current_hdl_actual)

 
% [t1, y1] = ode45(@quad_3d_ode, tspan, y0, options, current_hdl);
%ode_nom: ode of nominal system
%ode_actual: ode of actual system 

global dt;
% dt=0.01;

n_state=length(y0);   %dimension of the states

t1=tspan(1):dt:tspan(end);
n_time=length(t1);

%nominal state: 
y1_nom=zeros(n_time, n_state);
y1_nom(1,:) = y0';

%actual state: 
y1_actual=zeros(n_time, n_state);
y1_actual(1,:) = y0';

%state of each step  
y_nom=y0;
y_actual = y0;

%control: 
u1 = zeros(n_time, 2);
u1(1,:) = [0; 0];

%loop, in each step, update the nominal state and the actual state: 
% for i=2:n_time    
for i=1:n_time-1
    t=t1(i);   

%4th-order solver: 
    [dyk1, u_nomk1] = feval(ode_nom, t, y_nom,  current_hdl_nom);
    [dyk2, u_nomk2] = feval(ode_nom, t+0.5*dt, y_nom+0.5*dt*dyk1,  current_hdl_nom);
    [dyk3, u_nomk3] = feval(ode_nom, t+0.5*dt, y_nom+0.5*dt*dyk2,   current_hdl_nom);
    [dyk4, u_nomk4] = feval(ode_nom, t+dt, y_nom+dt*dyk3,  current_hdl_nom);
    y_nom=y_nom+ dt/6*(dyk1+2*dyk2+2*dyk3+dyk4);
    u_nom = u_nomk4;
    
    if(y_nom(1)<1e-2)
        y_nom(1)=0;
        u_nom = zeros(2,1);
    end
 
% 1st order solver: 
%     [dy, u_nom] = feval(ode_nom, t+dt, y_nom,  current_hdl_nom);
%     y_nom=y_nom+dy*dt; 
    
    y1_nom(i+1,:)=y_nom';
    
    %update actual state and record it: 
    [dy, u] = feval(ode_actual, t, y_actual, y_nom, u_nom, current_hdl_actual);
    y_actual=y_actual+dy*dt; 
    
    if(y1_actual(1)<1e-2)
        %stop:
        y1_actual(1)=0;
        u=zeros(2,1);
    end
    
    y1_actual(i+1,:)=y_actual';
    u1(i+1,:) = u; 
    
end