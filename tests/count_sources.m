clear all, close all
source_nr = 110;
max_sample = 50;
indexes = load(['relevant_index_',num2str(source_nr),'.txt'], '-ascii');
CrossCorrs = zeros(5,max_sample-1);
ii = 1;
for index = indexes
    disp('Index:')
    disp(index)
    label = load(['label_', num2str(index), '.mat']);
    output = load(['output_', num2str(index), '.mat']);
    label_data = label.vals;
    output_data = output.vals;

    outp_pos = squeeze(output_data(1,1,:,:));
    outp_rad = squeeze(output_data(1,2,:,:));
    outp_I = squeeze(output_data(1,3,:,:));
    outp_depth = squeeze(output_data(1,4,:,:));
    outp_dtI = squeeze(output_data(1,5,:,:));
    label_pos = squeeze(label_data(1,:,:));
    label_rad = squeeze(label_data(2,:,:));
    label_I = squeeze(label_data(3,:,:));
    label_depth = squeeze(label_data(4,:,:));
    label_dtI = squeeze(label_data(5,:,:));

    EstCoords=coordDerive(label_pos,.2,.001);
    sourceCount = length(EstCoords);
    CrossCorrs(1,ii) = corr2(outp_pos, label_pos);
    CrossCorrs(2,ii) = corr2(outp_rad, label_rad);
    CrossCorrs(3,ii) = corr2(outp_I, label_I);
    CrossCorrs(4,ii) = corr2(outp_depth, label_depth);
    CrossCorrs(5,ii) = corr2(outp_dtI, label_dtI);
    ii = ii + 1;
    disp(ii)
end

disp('Relevant indices:')
disp(indexes)

disp('Cross-correlations:')
disp(CrossCorrs)

% Load a test example
output = load('output_161.mat');
label = load('label_161.mat');
outp_pos1 = squeeze(output.vals(1,1,:,:));
figure(1), imagesc(outp_pos1)
label_pos1 = squeeze(label.vals(1,:,:));
figure(2), imagesc(label_pos1)

