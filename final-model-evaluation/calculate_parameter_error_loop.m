function [err,prefc]=calculate_parameter_error_loop(p,pref,m,n)
mu = 4*pi*1e-7;
dmax = 90;
amax = 180;
% rmin = 3.25e5;
rmax = 9e5;
% Imin = 2.2e8;
Imax = 1e9;
depth_max = (3.48-0.15)*1e6; %3.48e6;
depth_min = 2.7e6;
dtB_max = 7.33e-13; %4e-11 with the lowest std
dtI_max = 2*(((depth_max-depth_min)*1e6)^2+(3.25*1e5)^2)^(3/2)*(mu*(3.25*1e5)^2)^(-1)*dtB_max; %1.8182 [A/s]
attenuate=1100;
dtI_max = dtI_max/attenuate; %Bind(dIt_max) ~ 5.59*Brmax


p(1,:)=exp(p(1,:));
p(2,:)=exp(p(2,:));
pref(1,:)=exp(pref(1,:));
pref(2,:)=exp(pref(2,:));

pd(1,:)=2*pi*p(1,:)/180;
pd(2,:)=2*pi*p(2,:)/180;
pdr(1,:)=2*pi*pref(1,:)/180;
pdr(2,:)=2*pi*pref(2,:)/180;

xp=p(7,:).*sin(pd(1,:)).*cos(pd(2,:));
yp=p(7,:).*sin(pd(1,:)).*sin(pd(2,:));
zp=p(7,:).*cos(pd(1,:));
xr=pref(7,:).*sin(pdr(1,:)).*cos(pdr(2,:));
yr=pref(7,:).*sin(pdr(1,:)).*sin(pdr(2,:));
zr=pref(7,:).*cos(pdr(1,:));

err=zeros(m,n);
prefc=zeros(m,n);
nn = zeros(1,n);
% 
% xc=zeros(1,n);
% yc=zeros(1,n);
% zc=zeros(1,n);
% disp(nn)
disp('REFCOORDS')

for pp=1:length(xp)
    disp(pp)
    xp(pp) = round(xp(pp));
    yp(pp) = round(yp(pp));
    zp(pp) = round(zp(pp));
    %disp(min((xp(pp)-xr).*(xp(pp)-xr)+(yp(pp)-yr).*(yp(pp)-yr)+(zp(pp)-zr).*(zp(pp)-zr)))
    for rr=1:length(xr)
        xr(rr) = round(xr(rr));
        yr(rr) = round(yr(rr));
        zr(rr) = round(zr(rr));
        aa = (xp(pp)-xr(rr))*(xp(pp)-xr(rr))+(yp(pp)-yr(rr))*(yp(pp)-yr(rr))+(zp(pp)-zr(rr))*(zp(pp)-zr(rr));
        bb = min((xp(pp)-xr).*(xp(pp)-xr)+(yp(pp)-yr).*(yp(pp)-yr)+(zp(pp)-zr).*(zp(pp)-zr));
        if abs(aa-bb)<20          
            prefc(:,pp) = pref(:,rr);
        end
    end
end

disp(prefc)
disp(p)

delta=zeros(m,n);
MAXI = [dmax;amax;rmax;Imax;dmax;amax;depth_max;dtI_max];
for ii = 1:n
    for jj = 1:m
     delta(jj,ii)=abs(prefc(jj,ii)-p(jj,ii));
     if delta(jj,ii)>=MAXI(jj)/2 && (jj~=4 || jj~=7)
         disp(jj)
        delta(jj,ii) = MAXI(jj)-delta(jj,ii);
     end
    end
end

for ii = 1:n
    disp(delta(:,ii))
    disp(prefc(:,ii))
    err(:,ii) = delta(:,ii)./(2*abs(prefc(:,ii)));
end
err(1,:)=err(1,:).*prefc(1,:)/90;
err(2,:)=err(2,:).*prefc(2,:)/180;
err(5,:)=err(5,:).*prefc(5,:)/90;
err(6,:)=err(6,:).*prefc(6,:)/180;
