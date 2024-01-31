function RPE=calculate_parameter_error_loop(est_loops,true_loops,deg_res)
%%Calculates the relative parameter error for a given GA estimation of the loop parameters 
%%using the normalized parameter values
%Input params:
%-------------
%est_loops: array containing the estimated loop parameters
%true_loops: array containing the true loop parameters
%deg_res: angular resolution of images (maps) in degrees
%Output values:
%--------------
%RPE: array containing normalized relative parameter estimation errors

nearest_true_loops = get_nearest_true_loops(est_loops, true_loops, deg_res);
abs_esttrue_differences = calculate_abs_esttrue_differences(est_loops, nearest_true_loops, deg_res);
[Imax, dtI_max] = get_scaling_values();
%Tolerance to avoid very large relative errors in case of small current or dtI values (these values can cross zero)
tolerance = [0;0;0;Imax/100;0;0;0;dtI_max/100]; %for small true current or dtI values
RPE = abs_esttrue_differences./bsxfun(@plus,tolerance, abs(nearest_true_loops));


function nearest_true_loops = get_nearest_true_loops(est_loops, true_loops, deg_res)
nr_est_loops = length(est_loops(1,:));
nr_true_loops = length(true_loops(1,:));

fis_e=deg_res*pi*est_loops(1,:)/180;
lambdas_e=deg_res*pi*est_loops(2,:)/180;
fis_t=deg_res*pi*true_loops(1,:)/180;
lambdas_t=deg_res*pi*true_loops(2,:)/180;

distance_earth_centre_toloop_e = est_loops(7,:);
distance_earth_centre_toloop_t = true_loops(7,:);

xe=distance_earth_centre_toloop_e.*sin(fis_e).*cos(lambdas_e);
ye=distance_earth_centre_toloop_e.*sin(fis_e).*sin(lambdas_e);
ze=distance_earth_centre_toloop_e.*cos(fis_e);
xt=distance_earth_centre_toloop_t.*sin(fis_t).*cos(lambdas_t);
yt=distance_earth_centre_toloop_t.*sin(fis_t).*sin(lambdas_t);
zt=distance_earth_centre_toloop_t.*cos(fis_t);

nearest_true_loops=zeros(size(est_loops));

for est_loop=1:nr_est_loops
    minimum_distance = min((xe(est_loop)-xt).^2+(ye(est_loop)-yt).^2+(ze(est_loop)-zt).^2);
    for true_loop=1:nr_true_loops
        distance = (xe(est_loop)-xt(true_loop))^2+(ye(est_loop)-yt(true_loop))^2+(ze(est_loop)-zt(true_loop))^2;
        if distance == minimum_distance          
            nearest_true_loops(:,est_loop) = true_loops(:,true_loop);
        end
    end
end

function [Imax, dtI_max] = get_scaling_values()

mu = 4*pi*1e-7;
rmin = 3.25e5;
Imax = 1e9;
depth_max = 3.48e6;
depth_min = 2.7e6;
dtB_max = 7.33e-13; 
dtI_max = 2*((depth_max-depth_min)^2+(rmin)^2)^(3/2)*(mu*(rmin)^2)^(-1)*dtB_max; %1.8182 [A/s]
attenuate=1100;
dtI_max = dtI_max/attenuate; %

function abs_esttrue_differences = calculate_abs_esttrue_differences(est_loops, nearest_true_loops, deg_res)
se = size(est_loops);
abs_esttrue_differences=zeros(size(est_loops));


for ii = 1:se(1)
    for jj = 1:se(2)
     abs_esttrue_differences(ii,jj)=abs(nearest_true_loops(ii,jj)-est_loops(ii,jj));
     %Correction should be made if the closest source in space is on the opposite edge of the map in longitude.
     if abs_esttrue_differences(ii,jj)>=180/deg_res && (ii==2)
        abs_esttrue_differences(ii,jj) = abs(360/deg_res-abs_esttrue_differences(ii,jj));
     end
    end
end
