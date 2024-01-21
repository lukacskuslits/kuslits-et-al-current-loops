function loop=loop_gen_constr_rad_total(limits,nloops,dist_lat,dist_lon,dist_Ilat,dist_Ilon,deg_res)
%% generates parameters of radially aligned loops and their randomly
%% for nloop current loops
%Input params:
% limits: vector containing the lower and upper limit for each
    % corresponding loop parameter - the ranges were derived from REF.m
% nloops: number of loops to generate the loop parameters for
% dist_lat: nearest distance possible between neighboring loops in latitude (degrees)
% dist_lon: nearest distance possible between neighboring loops in longitude (degrees)
% dist_Ilat: KELL?
% dist_Ilon: KELL?
% deg_res: lateral resolution for generating the geographic positions of the loops (in degrees)

d1=limits(1);
d2=limits(2);
a1=limits(3);
a2=limits(4);
R1=limits(5);
R2=limits(6);
I1=limits(7);
I2=limits(8);
Isign1=limits(9);
Isign2=limits(10);
t1=limits(11);
t2=limits(12);
l1=limits(13);
l2=limits(14);
z1=limits(15);
z2=limits(16);
tI1=limits(17);
tI2=limits(18);
v1 = limits(19);
v2 = limits(20);


loop=zeros((length(limits)-2)/2,nloops);
% sr_number(1) = n;
% sr_number = double(sr_number);
ii=1;
retry = 0;
%outl = 0;
while ii<nloops+1
    loop(1,ii)=randi([d1,d2]);
    loop(2,ii)=randi([a1,a2]);
    loop(3,ii)=randi([R1*1e5,R2*1e5]);%;%abs(round(normrnd((R2+R1)/2,(R2-R1)/500)*1e5));
    %---------------------------------------------------
    Isigns = [Isign1,Isign2];
    loop(4,ii)=randi([I1,I2])*Isigns(randi(2));%randi([I1,I2])*Isigns(randi(2));%normrnd((I2+I1)/50,(I2-I1)/500)*Isign;
    %---------------------------------------------------
    loop(5,ii)=deg_res(1)*loop(1,ii); %;
    loop(6,ii)=deg_res(2)*loop(2,ii); %;


    loop(7,ii)=randi([z1*1e6,z2*1e6]);%;%normrnd((z2+z1)/2.1,(z2-z1)/10)*1e6;
    %---------------------------------------------------
    dtI_signs = [-1,1];
    loop(8,ii)=randi([tI1,round(tI2*1e4)])*1e-4*dtI_signs(randi(2));%*dtI_signs(randi(2));%normrnd((tI2+tI1)/2,(tI2-tI1)/100)*dtI_signs(randi(2));%(round(normrnd(3/2.3,1/6)));normrnd(tImean,tIstd);
    %---------------------------------------------------

    loop(9,ii)=round(randi([v1,v2]))/(3.15e7*1e2);
    loop(10,ii)=round(randi([v1,v2]))/(3.15e7*1e2);  
    loop(1,ii)=deg_res(1)*loop(1,ii); %;
    loop(2,ii)=deg_res(2)*loop(2,ii); %;
    Roc = 3.48e6;
    mu = 4*pi*1e-7;
    I = loop(4,ii);
    r = loop(3,ii);
    z = Roc -loop(7,ii);
    dtI = loop(8,ii);
    Br_static = mu*I*r^2*(r^2 + z^2)^(-3/2);
    Br_induced = dtI*ind_approx_single(0,z/8e5,r/8e5);
    Br_above =  Br_static + Br_induced;
    if (Br_above>0.0032/log(nloops)) || (Br_above<-0.0031/log(nloops))
        %disp('Loop field too strong, discarded')
        ii = ii-1;
        retry = retry + 1;
    elseif (Br_above/Br_static)<1/exp(1)
        %disp('Too large induction')
        ii = ii-1;
        retry = retry + 1;
    elseif (r - sqrt(Roc^2 - loop(7,ii)^2)>-5e4)
        %disp('Too close to surf')
        ii = ii-1;
        retry = retry + 1;
    elseif ii>1
        diff_arr= abs(loop(1,ii)-loop(1,1:ii-1))+abs(loop(2,ii)-loop(2,1:ii-1));
        dist_min=min(diff_arr);
        nearest_neighbor = diff_arr==dist_min;
        compare_signs = loop(4,ii)*prod(loop(4,nearest_neighbor));
        dmin = abs(loop(1,ii) - loop(1,nearest_neighbor));
        amin = abs(loop(2,ii) - loop(2,nearest_neighbor));
        if sum(dmin<dist_Ilat)>0 && sum(amin<dist_Ilon)>0
            if compare_signs<0
               ii=ii-1;
               retry = retry+1;
               %disp('Too close to each other')
            elseif sum(dmin<dist_lat)>0 && sum(amin<dist_lon)>0
               ii=ii-1;
               retry = retry+1;
               %disp('Too close to each other')
           end
        end
        clear diff_arr
        clear nearest_neighbor
    end
    ii = ii+1;
    if retry>100
        error('Allj le!')
        break
    end
end
end




