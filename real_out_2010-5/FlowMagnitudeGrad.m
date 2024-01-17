clear all, close all
%[pd_o,B,CoordShallow]=eval_net_meas('net_10_150_shallow_aa.mat','B_R_core20102.txt',[10,10],2,45,90);

%load('CoordShallowMeas.mat')
load('FlowModel2010_lat_lon_angle_length_4deg.dat')
flow=FlowModel2010_lat_lon_angle_length_4deg;
flow(:,1)=(flow(:,1)-min(flow(:,1)))/4;
flow(:,2)=(flow(:,2)-min(flow(:,2)))/4;
flow(:,2)=sort(flow(:,2),'ascend');
flowq=[flow(:,1:2)+1,flow(:,4)];
flowq(:,1)=shiftvec(flowq(:,1),23)';

sf=size(flowq);


F=zeros(max(flowq(:,2)),max(flowq(:,1)));
for ii=1:sf(1)

    F(flowq(ii,2),flowq(ii,1))=flowq(ii,3);
    
end   


figure(6), quiver(F)


GF=gradient(F);

figure(7), imagesc(GF)
hold on, load('output_2010.mat')
figure(6), contour(squeeze(vals(1,1,:,:)), 'black')
