function loop=loop_gen_constr_rad_just_coords(limits,n, dist_lat,dist_lon)

d1=limits(1);
d2=limits(2);
a1=limits(3);
a2=limits(4);


loop=zeros((length(limits)-2)/2,n);
ii=1;
retry = 0;
%outl = 0;
while ii<n+1
    loop(1,ii)=randi([d1,d2]);
    loop(2,ii)=randi([a1,a2]);

    if ii>1
        diff_arr= abs(loop(1,ii)-loop(1,1:ii-1))+abs(loop(2,ii)-loop(2,1:ii-1));
        dist_min=min(diff_arr);
        nearest_neighbor = diff_arr==dist_min;      
        dmin = abs(loop(1,ii) - loop(1,nearest_neighbor));
        amin = abs(loop(2,ii) - loop(2,nearest_neighbor));
   
        if sum(dmin<dist_lat)>0 && sum(amin<dist_lon)>0
               ii=ii-1;
               retry = retry+1;
               disp('Tul kozeli')
        end
    end
    ii = ii+1;
    if retry>100
        error('Allj le!')
        break
    end
end
end


