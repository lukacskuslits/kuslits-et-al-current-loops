% %% Count how many sources were 'detected' using Deep Learning
clear all, close all
source_nr = 20;
indexes = load(['relevant_index_',num2str(source_nr),'.txt'], '-ascii');
indexes = indexes(indexes>0);
% Thershold value to define the successful identifitacion of a loop
threshold = 2; 

discovered_loops_all = zeros(1,length(indexes));
ii = 1;
for index = indexes
    disp('Index:')
    disp(index)
    label = load(['label_', num2str(index), '.mat']);
    output = load(['output_', num2str(index), '.mat']);
    if source_nr<30
     label_data = label.label;
    else
     label_data = label.vals;
    end

    output_data = output.vals;

    outp_pos = squeeze(output_data(1,1,:,:));
    label_pos = squeeze(label_data(1,:,:));
    EstCoords=coordDerive(outp_pos,.2,.001);
    TrueCoords=coordDerive(label_pos,.2,.001);
    nr_detected_loops = 0;
    for col = 1:length(TrueCoords(1,:))
        true_coord = TrueCoords(:,col);
        true_coord_arr = bsxfun(@times,true_coord,ones(size(EstCoords)));
        diff_arr = true_coord_arr-EstCoords;
        diff_arr = sum(abs(diff_arr));
        nr_detected_loops = nr_detected_loops+min(1,sum(diff_arr<threshold));
    end
    discovered_loops_all(ii) = nr_detected_loops;
    ii = ii+1;
end

average_nr_of_identified_loops = mean(discovered_loops_all);
average_discovery_rate = mean(discovered_loops_all/source_nr);
set_size = length(indexes);

%% In the manuscript
% Nr.| Dl   | Dr    | Set
% 5  | 1    | 20%   | 2
% 10 | 0    | 0%    | 2
% 15 | 0    | 0%    | 2
% 20 | 1    | 5%    | 2
% 30 | 24.05| 80.2% | 20
% 50 | 35.9 | 71.8% | 50
% 70 | 46.2 | 66%   | 50
% 90 | 56.5 | 62.8% | 50
% 110| 64.16| 58.3% | 50
% 130| 72.73| 56%   | 15
% 150|------------->| 0

%(With a threshold value of 10:
% 5  | 3    | 60%   | 2
% 10 | 7    | 70%   | 2
% 15 | 9    | 60%   | 2
% 20 | 13   | 65%   | 2
% 30 | 27.65| 92.2% | 20
% 50 | 44.98| 90%   | 50
% 70 | 63.6 | 90.9% | 50
% 90 | 83.4 | 92.7% | 50
% 110| 101.6| 92.4% | 50
% 130| 120.6| 93%   | 15
% 150|------------->| 0)

%% Count how many sources were 'detected' using Alldredge's method by 
%% selecting local maxima on the radial field input maps
clear all, close all
source_nr = 130; %110;
indexes = load(['relevant_index_',num2str(source_nr),'.txt'], '-ascii');
indexes = indexes(indexes>0);

% Thershold value to define the successful identifitacion of a loop
threshold = 10; 

discovered_loops_all = zeros(1,length(indexes));
ii = 1;
for index = indexes
    disp('Index:')
    disp(index)
    label = load(['label_', num2str(index), '.mat']);
    input = load(['input_', num2str(index), '.mat']);
    %variable names were assigned differently for models having a number of loops
    %less then those in the Deep Learning training set
    if source_nr<30
     label_data = label.label;
     input_data = input.input;
    else
     label_data = label.vals;
     input_data = input.vals;
    end

    input_field = squeeze(input_data(1,:,:));
    
    %Get ormalized absolute input field data for the gradient method
    input_field = (abs(input_field)-min(min(abs(input_field))))/ (max(max(abs(input_field)))-min(min(abs(input_field))));
    label_pos = squeeze(label_data(1,:,:));
    nr_detected_loops = 0;
    try
    EstCoords=coordDerive(abs(input_field),.01,.001);
    TrueCoords=coordDerive(label_pos,.01,.001);
    for col = 1:length(TrueCoords(1,:))
        true_coord = TrueCoords(:,col);
        true_coord_arr = bsxfun(@times,true_coord,ones(size(EstCoords)));
        diff_arr = true_coord_arr-EstCoords;
        diff_arr = sum(abs(diff_arr));
        nr_detected_loops = nr_detected_loops+min(1,sum(diff_arr<threshold));
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

zed = 10*ones(size(EstCoords(1,:)));
figure(1), imagesc(input_field) 
figure(1), hold on, scatter3(EstCoords(2,:), EstCoords(1,:), zed, 100, 'black', 'filled')

figure(2), imagesc(label_pos)
%% In the manuscript
%Threshold is 2 pixels (2/90 ~ within 2% error)
% Nr.| Dl   | Dr    | Set
% 5  | 2    | 40%   | 2
% 10 | 2.5  | 25%   | 2
% 15 | 5.5  | 37%   | 2
% 20 | 3    | 15%   | 2
% 30 | 2.05 | 6.8%  | 20
% 50 | 4.1  | 8.2%  | 50
% 70 | 5.14 | 7.3%  | 50
% 90 | 5.28 | 5.9%  | 50
% 110| 5.92 | 5.4%  | 50
% 130| 6.47 | 5%    | 15
% 150|------------->| 0

%Threshold is 10 pixels (10/90 ~ within 10% error)
% Nr.| Dl   | Dr    | Set
% 5  | 5    | 100%  | 2
% 10 | 9    | 90%   | 2
% 15 | 12.5 | 83%   | 2
% 20 | 16.5 | 82.5% | 2
% 30 | 18.2 | 60.7% | 20
% 50 | 27.8 | 55.5% | 50
% 70 | 33.9 | 56.3% | 50
% 90 | 51.2 | 56.8% | 50
% 110| 62.6 | 56.9% | 50
% 130| 77.4 | 59.6% | 15
% 150|------------->| 0