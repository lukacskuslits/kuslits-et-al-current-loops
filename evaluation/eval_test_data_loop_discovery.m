%% Count how many sources were 'detected' using Deep Learning
clear all, close all
source_nr = 110;
indexes = load(['relevant_index_',num2str(source_nr),'.txt'], '-ascii');

discovered_loops_all = zeros(1,length(indexes));
ii = 1;
for index = indexes
    disp('Index:')
    disp(index)
    label = load(['label_', num2str(index), '.mat']);
    output = load(['output_', num2str(index), '.mat']);
    label_data = label.vals;
    output_data = output.vals;

    outp_pos = squeeze(output_data(1,1,:,:));
    label_pos = squeeze(label_data(1,:,:));
    EstCoords=coordDerive(outp_pos,.2,.001);
    TrueCoords=coordDerive(label_pos,.2,.001);
    nr_detected_loops = 0;
    for col = 1:length(EstCoords(1,:))
        est_coord = EstCoords(:,col);
        est_coord_arr = bsxfun(@times,est_coord,ones(size(TrueCoords)));
        diff_arr = est_coord_arr-TrueCoords;
        diff_arr = sum(abs(diff_arr));
        nr_detected_loops = nr_detected_loops+min(1,sum(diff_arr<2));
    end
    discovered_loops_all(ii) = nr_detected_loops;
    ii = ii+1;
end

average_nr_of_identified_loops = mean(discovered_loops_all);
average_discovery_rate = mean(discovered_loops_all/source_nr);
set_size = length(indexes);

%% In the manuscript
% Nr.| Dl   | Dr    | Set
% 30 | 25   | 83.3% | 20
% 50 | 35.9 | 71.8% | 50
% 70 | 47   | 67.1% | 50
% 90 | 57.3 | 64%   | 50
% 110| 64.9 | 60%   | 50
% 130

%% Count how many sources were 'detected' using Alldredge's method by 
%% selecting local maxima on the radial field input maps
clear all, close all
source_nr = 110;
indexes = load(['relevant_index_',num2str(source_nr),'.txt'], '-ascii');

discovered_loops_all = zeros(1,length(indexes));
ii = 1;
for index = indexes
    disp('Index:')
    disp(index)
    label = load(['label_', num2str(index), '.mat']);
    input = load(['input_', num2str(index), '.mat']);
    label_data = label.vals;
    input_data = input.vals;

    input_field = squeeze(input_data(1,:,:));
    label_pos = squeeze(label_data(1,:,:));
    nr_detected_loops = 0;
    try
    EstCoords=coordDerive(input_field,.2,.001);
    TrueCoords=coordDerive(label_pos,.2,.001);
    for col = 1:length(EstCoords(1,:))
        est_coord = EstCoords(:,col);
        est_coord_arr = bsxfun(@times,est_coord,ones(size(TrueCoords)));
        diff_arr = est_coord_arr-TrueCoords;
        diff_arr = sum(abs(diff_arr));
        nr_detected_loops = nr_detected_loops+min(1,sum(diff_arr<2));
    end
    catch ME
      disp(ME)
    end
    discovered_loops_all(ii) = nr_detected_loops;
    ii = ii+1;
end

average_nr_of_identified_loops = mean(discovered_loops_all);
average_discovery_rate = mean(discovered_loops_all/source_nr);
set_size = length(indexes);

%% In the manuscript
% Nr.| Dl   | Dr    | Set
% 30 | 2.05 | 6.8%  | 20
% 50 | 4.1  | 8.2%  | 50
% 70 | 5.14 | 7.3%  | 50
% 90 | 5.28 | 5.9%  | 50
% 110| 5.92 | 5.4%  | 50
% 130
