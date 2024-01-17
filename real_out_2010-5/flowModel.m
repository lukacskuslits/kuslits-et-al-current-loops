clear all, close all

load('output_2010.mat')
est_pos = squeeze(vals(1,1,:,:));
%est_depth = squeeze(vals(1,2,:,:));
%est_depth(est_depth>0.7)=0;
%est_post = est_pos.*(1-est_depth);
%est_pos(est_pos<0.9)=0;

est_pos = mirror(est_pos,'row','surfa');
    
%scatter3(gridpoints(2,:), gridpoints(1,:), Z, 100, c, 'filled');

load('FlowModel2010_lat_lon_angle_length_4deg.dat')

flow=FlowModel2010_lat_lon_angle_length_4deg;

flow(:,1)=(flow(:,1)-min(flow(:,1)))/4;
flow(:,2)=(flow(:,2)-min(flow(:,2)))/4;

L=length(flow(:,1));
flowq=zeros(L,4);
MF=max(flow(:,2));

%flow(:,2)=flipdim(flowq(:,2),1);
flowt(:,2)=flow(:,2);
flow(:,2)=sort(flow(:,2),'ascend');
% for jj=1:L
% flowq(jj,2)=MF-flow(jj,2)+1;
% end

flowq(:,1)=flow(:,1);

flowq(:,2)=flow(:,2);

%flowq(:,1)=flipdim(flowq(:,1),1);
%flowq(:,1:2)=flipdim(flowq(:,1:2),2);

flowq(:,1)=shiftvec(flowq(:,1),23)';

flowq(:,3)=flow(:,4).*sin(deg2rad(flow(:,3)));
flowq(:,4)=-flow(:,4).*cos(deg2rad(flow(:,3)));


figure(5), 
imagesc(est_pos)
hold on
quiver(flowq(:,1),flowq(:,2),10*flowq(:,3),10*flowq(:,4), 'black', 'AutoScaleFactor',2)


load('input_2010.mat')
true_field = squeeze(vals(1,:,:));
%figure(4), hold on, contour(true_field, 'green')
% 
% figure(5), imagesc(true_field)
% figure(5), hold on, quiver(flowq(:,1),flowq(:,2),flowq(:,3),flowq(:,4), 'black')

flowMagMap = zeros(45,90);
for ii=1:length(flowq(:,1))
        flowMagMap(flowq(ii,2)+1,flowq(ii,1)+1) = flowq(ii,3)^2+flowq(ii,4)^2;
end

norm_logmap = log(flowMagMap);
norm_logmap(1,:) = norm_logmap(2,:);
norm_logmap(:,1) = norm_logmap(:,2);
norm_logmap(45,:) = norm_logmap(44,:);
norm_logmap = abs(norm_logmap);
norm_logmap = (norm_logmap-min(min(norm_logmap)))/(max(max(norm_logmap))-min(min(norm_logmap)));
norm_logmap = gradient(norm_logmap);
figure(5), imagesc(norm_logmap)
norm_logmap = mirror(norm_logmap,'row','surfa');
norm_logmap = abs(norm_logmap);
norm_logmap = (norm_logmap-min(min(norm_logmap)))/(max(max(norm_logmap))-min(min(norm_logmap)));
norm_logmap(norm_logmap>0.75)=0;
figure(6), imagesc(norm_logmap)
EstCoords = est_pos;
EstCoords = coordDerive(mirror(EstCoords,'row','surfa'),0.2,0.001,6);

%% Vajon tényleg nem véletlenszerűen rakosgatja le a forrásokat?
compare_grad = 0;
for ii=1:length(EstCoords(1,:))
    compare_grad = compare_grad+norm_logmap(EstCoords(1,ii),EstCoords(2,ii));
end
disp('compare grad')
disp(compare_grad)

random_grads = zeros(1000,1);
for jj = 1:1000
limits = [1,45,1,90];
RandomCoords = loop_gen_constr_rad_just_coords(limits,length(EstCoords(1,:)),2,2);

random_grad = 0;
for ii=1:length(EstCoords(1,:))
    random_grad = random_grad+abs(norm_logmap(RandomCoords(1,ii),RandomCoords(2,ii)));
end
random_grads(jj) = random_grad;
disp('random_grad')
disp(random_grad)
end

figure(1), hist(random_grads)
xlabel('Anyagáramlási gradiensek összege a véletlen mintában')
ylabel('Minta darabszáma')
%'$(\sum_i{abs(\nabla{log(abs(u_h))}})$'
% BR = load('input_2010.mat');
% BR = squeeze(BR.vals(1,:,:));
% 
% 
% [X, Y] = meshgrid(45,90);
% X = X';
% Y = Y';
% X1(:,:,1:3) = X;
% Y1(:,:,1:3) = Y;
% Z1(:,:,1:3) = 3.3e6*ones(45,90,1);
% Z1(:,:,1) = Z1(:,:,2)-1e5;
% Z1(:,:,3) = Z1(:,:,2)+1e5;
% U = zeros(45,90,3);
% V = zeros(45,90,3);
% W(:,:,1) = BR-0.1;
% W(:,:,2) = BR;
% W(:,:,3) = BR+0.1;




