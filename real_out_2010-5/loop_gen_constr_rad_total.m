function loop=loop_gen_constr_rad_total(limits,n,dist_lat, dist_lon, dist_Ilat,dist_Ilon,radial,cone_angle)

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


loop=zeros((length(limits)-2)/2,n);
ii=1;
retry = 0;
%outl = 0;
while ii<n+1
    loop(1,ii)=randi([d1,d2]);
    loop(2,ii)=randi([a1,a2]);
    %loop(3,ii)=randi([R1*1e5,R2*1e5]);
    %loop(3,ii) = round(val_from_dist(R1,R2,1e5));
    
    %---------------------------------------------------
    %loop3 generated with params: 
    %2.2 500
    %2.3, 500
    %2, 500
    %2, 700
    loop(3,ii)=randi([R1*1e5,R2*1e5]);%;%abs(round(normrnd((R2+R1)/2,(R2-R1)/500)*1e5));
    %---------------------------------------------------
    
    Isigns = [Isign1,Isign2];
    %loop(4,ii) = round(val_from_dist(I1,I2,1e8)*Isigns(randi(2)));
    %loop(4,ii)=randi([I1,I2])*Isigns(randi(2));
%     if loop(1,ii)<=22
%         Isign = Isigns(round(normrnd(3/1.7,1/6))); %(randi(2));
%     else
%         Isign = Isigns(round(normrnd(3/2.3,1/6)));
%     end
    %---------------------------------------------------
    %loop4 generated with params: 
    %45 500
    %120 500
    %100 500
    %30, 700
    loop(4,ii)=randi([I1,I2])*Isigns(randi(2));%randi([I1,I2])*Isigns(randi(2));%normrnd((I2+I1)/50,(I2-I1)/500)*Isign;
    %---------------------------------------------------
    
    if radial==0
       while 1==1
       loop(5,ii)=randi([t1,t2]);
       loop(6,ii)=randi([l1,l2]);
       dev_angle = deviation_angle(loop(:,ii));
       %disp(dev_angle) 
       if (abs(dev_angle) <= cone_angle(1)) || (abs(180-dev_angle) <= cone_angle(1))
       if length(cone_angle) > 1
          if (abs(dev_angle) >= cone_angle(2)) || (abs(180-dev_angle) >= cone_angle(2))
              break
          end
       else 
           break
       end 
       end
       end
    else
    loop(5,ii)=4*loop(1,ii); %;
    loop(6,ii)=4*loop(2,ii); %;
    end
      
    %loop(7,ii)=randi([z1*1e6,z2*1e6]);
    %loop(7,ii)=val_from_dist(z1,z2,1e6);
    
    %---------------------------------------------------
    %loop7 generated with params: 
    %2.1 70
    %2.2 100
    %2.1, 500
    %2.1, 500
    loop(7,ii)=randi([z1*1e6,z2*1e6]);%;%normrnd((z2+z1)/2.1,(z2-z1)/10)*1e6;
    %---------------------------------------------------
    
    dtI_signs = [-1,1];
    %loop(8,ii)=val_from_dist(0,tI2,1e-6)*dtI_signs(randi(2));
    %loop(8,ii)=randi([tI1,tI2])/1e5;
    
    %---------------------------------------------------
    %loop8 generated with params:
    %15 500
    %25 500
    %50, 500
    %50, 700
%     disp(tI1)
%     disp(round(tI2*1e4))
    loop(8,ii)=randi([tI1,round(tI2*1e4)])*1e-4*dtI_signs(randi(2));%*dtI_signs(randi(2));%normrnd((tI2+tI1)/2,(tI2-tI1)/100)*dtI_signs(randi(2));%(round(normrnd(3/2.3,1/6)));normrnd(tImean,tIstd);
    %---------------------------------------------------
    
    %loop(9,ii)=round(normrnd(v1,v2/100))/(3.15e7*1e1);
    %loop(10,ii)=round(normrnd(v1,v2/100))/(3.15e7*1e1);
    loop(9,ii)=round(randi([v1,v2]))/(3.15e7*1e2);
    loop(10,ii)=round(randi([v1,v2]))/(3.15e7*1e2);  
    if radial~=0
       loop(1,ii)=4*loop(1,ii); %;
       loop(2,ii)=4*loop(2,ii); %;
    end
    Roc = 3.48e6;
    mu = 4*pi*1e-7;
    I = loop(4,ii);
    r = loop(3,ii);
    z = Roc -loop(7,ii);
    dtI = loop(8,ii);
%     disp('dtI')
%     disp(dtI)
    Br_static = mu*I*r^2*(r^2 + z^2)^(-3/2);
    Br_induced = dtI*ind_approx_single(0,z/8e5,r/8e5);
    Br_above =  Br_static + Br_induced;
%     disp('Induced field')
%     disp(dtI*ind_approx_single(0,z/8e5,r/8e5))
    if (Br_above>0.0032/log(n)) || (Br_above<-0.0031/log(n))
        %disp('Loop field too strong, discarded')
        ii = ii-1;
        retry = retry + 1;
        %outl = outl + 1;
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
               %disp('Tul kozeli')
            elseif sum(dmin<dist_lat)>0 && sum(amin<dist_lon)>0
               ii=ii-1;
               retry = retry+1;
               %disp('Tul kozeli')
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


