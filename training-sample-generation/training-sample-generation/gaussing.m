function gsurf=gaussing(surface,sigma,v,h)
%% Generates a 2d gaussian blurring where sources can be found in 
%% the ground truth data in the training set
% Input parameters
%-----------------
% surface: the original unblurred (sparse) matrix of source locations
% sigma: standard deviation of the applied gaussian filter
% v: meridional kernel size
% h: longitudinal kernel size
% Output parameters
%-----------------
% gsurf: preprocessed (blurred) 2d surface (map) of ground truth source
% parameters

ssize=size(surface);
%disp(ssize(1))
%disp(ssize(2))
gsurf=surface;
for ii=1:ssize(1)
    for jj=1:ssize(2)
        if abs(surface(ii,jj))>0
            
            if ii-v(1)>=1
            bl=ii-v(1)+1;
            else
            bl=1;
            end
            if ii+v(2)<=ssize(1)
            bu=ii+v(2);
            else
            bu=ssize(1);
            bl=ssize(1)-v(1)+1;
            end
            if jj-h(1)>=1
            dl=jj-h(1)+1;
            else
            dl=1;
            end
            if jj+h(2)<=ssize(2)
            du=jj+h(2);
            else
            du=ssize(2);
            dl=ssize(2)-h(1)-1;
            end
            %disp(bl), disp(bu), disp(dl), disp(du)
            mat=surface(bl:bu,dl:du);
            if surface(ii,jj)>0
               gsurf(bl:bu,dl:du)=...
               max(gsurf(bl:bu,dl:du),gauss2d(surface(ii,jj),mat,sigma,[ii-bl+1, jj-dl+1]));
            else
               gsurf(bl:bu,dl:du)=...
               -max(abs(gsurf(bl:bu,dl:du)),gauss2d(abs(surface(ii,jj)),mat,sigma,[ii-bl+1, jj-dl+1]));
            end
        end
    end
end
end