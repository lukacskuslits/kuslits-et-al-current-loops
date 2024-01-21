function set_initial_params(nr_of_loops)
if exist(['correct_pars_',num2str(nr_of_loops),'.mat'],'file')==2
    correct_pars = load(['correct_pars_',num2str(nr_of_loops),'.mat']);
    correct_pars = correct_pars.correct_pars;
    save('correct_pars.mat','correct_pars')
else
%% Loop parameters
% Our method
mu = 4*pi*1e-7;
n = nr_of_loops;  
rmin = 3.25; %3.25; %original upper and lower vals for R
rmax = 10; %10; %
Imin = 2.2e8;
Imax = 1e9;
depth_min = 2.6; %2.6; ICB
depth_max = 3.2; %3.48; CMB
dtB_max = 7e-13; %4e-11 with the lowest std
%Highest possible rate of change in current:
%dI/dt = 2*Rmax*(dB/dt)max/mu
dtI_max = 2*(((depth_max-depth_min)*1e6)^2+(3.25*1e5)^2)^(3/2)*(mu*(3.25*1e5)^2)^(-1)*dtB_max;
dtI_max = dtI_max/1100;
disp('dIt')
disp(dtI_max)
dtI_min = 0;
v_max = 0;
trials = 1000;
for i=1:trials
disp(i)
limits_static = [1;90;1;180;rmin;rmax;Imin;Imax;-1;1;0;0;0;0;depth_min;depth_max;dtI_min;0;-v_max;v_max];
try
    correct_pars = loop_gen_constr_rad_total(limits_static,n,3,3,16,16,[2,2]); %read_results
catch ME
   disp(ME)
   continue
end
break
end

%deg_res = 2;
%cords=pp_static(1:2,:)/deg_res;
correct_pars(5,:)=correct_pars(1,:)*pi()/180;
correct_pars(6,:)=correct_pars(2,:)*pi()/180;
correct_pars(1,:)=correct_pars(1,:)*pi()/180;
correct_pars(2,:)=correct_pars(2,:)*pi()/180;
%correct_pars = compute_alldredge_params(pp_static);
save('correct_pars.mat','correct_pars')
save(['correct_pars_',num2str(nr_of_loops),'.mat'],'correct_pars')
end