global scale  u_global scale_tracking u_tracking_global scale_record; 
% scale 1-by-1, int
% u_global  2-by-1, double
% scale_tracking  1-by-1, int
% u_tracking_global 2-by-1,double
% scale_record  1-by-1, int

global tra_com_pre tra_com_dot_pre tra_com_ddot_pre;  
% tra_com_pre  3-by-1, double
% tra_com_dot_pre  3-by-1, double
% tra_com_ddot_pre  3-by-1, double

global brake_flag brake_flag_pre;  
% brake_flag   1-by-1, int
% brake_flag_pre  1-by-1, int

global nosolution;
% nosolution 1-by-1, int

global no_ob; 
% no_ob: 1-by-1, int
global pos_ob_array_pre radius_pre;
% pos_ob_array_pre: 2-by-no_ob, double
% radius_pre: no_ob-by-1, double

global beta_2; 
global dt;  
% dt 1-by-1, double 
% beta_2  100-by-1, int

%suppose number of time steps: n
global t_ctrl; 
%t_ctrl:  n-by-1, double
global u_ctrl;
%u_ctrl:  11-by-n, double
global traj_ob_seris; 
% traj_ob_seris: (2*n)-by-no_ob, double


global trajd;
%trajd: 9-by-1, double 
 
