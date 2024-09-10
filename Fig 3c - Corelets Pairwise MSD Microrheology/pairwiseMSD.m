% 8/3/2023 JZ
%
% Adapted from Natalia's R script for RPA tracking
%
% input: filename of spots statistics file from trackmate
%
% outputs: cell array containing MSD info for each pair tracked, with each
% table saved to current directory as .csv file

function output = pairwiseMSD(filename)
% Read files:
frametosecond = 0.06; 

spots = sortrows(readtable(filename),[3,9]);

% Set a temporal track cutoff
lengthTrack = size(unique(spots.TRACK_ID),1); %returns number of tracks
[C,~,ic] = unique(spots.TRACK_ID);
track_label = C;
track_duration = accumarray(ic,1);
min_duration = 200; %in frames
real_tracks = [];
for i = 1 : lengthTrack
    if track_duration(i) >= min_duration
        real_tracks(end+1) = track_label(i);
    end
end
track_count = length(real_tracks);
duration_cutoff = min_duration;
max_tau = min_duration;
max_duration = max(track_duration);

%% For each of the tracks IDed as real, find the coexistence matrix
% matrix(first_track, second_track) = number of frames first and second track are co-existing
matrix = zeros(track_count);

% Tally coexistence periods of different long-lived tracks
for frame = 1:max_duration
    current_frame = spots(spots.FRAME == frame, :); % extracting ALL rows at current frame!
    for first_track = 1:track_count
        for second_track = first_track:track_count
            if any(real_tracks(first_track) == current_frame.TRACK_ID) && ...
                    any(real_tracks(second_track) == current_frame.TRACK_ID)
                matrix(first_track, second_track) = matrix(first_track, second_track) + 1;
                % matrix(second_track, first_track) = matrix(second_track, first_track) + 1;
            end
        end
    end
end

% Returns row and column indices that will let us locate decently co-existing pairs
first_entries = [];
second_entries = [];
for i = 1:size(matrix, 1)
    for j = 1:size(matrix, 2)
        if matrix(i, j) >= duration_cutoff && i ~= j
            first_entries = [first_entries, i];
            second_entries = [second_entries, j];
        end
    end
end

%% MSD Part
Taus = {};
Pair_IDs = {};
MSDs = {};
Mean_Dists = {};
Delta_Dists = {};

if isempty(first_entries) 
    output = {};
else
    % For a given pair, compute starting distance
    for index = 1:length(first_entries)
        first_spot = first_entries(index);
        second_spot = second_entries(index);
        
        % Pull out relevant position data
        first_xs = spots(spots.TRACK_ID == real_tracks(first_spot) & spots.FRAME > 0, :).POSITION_X;
        first_ys = spots(spots.TRACK_ID == real_tracks(first_spot) & spots.FRAME > 0, :).POSITION_Y;
        first_times = spots(spots.TRACK_ID == real_tracks(first_spot) & spots.FRAME > 0, :).FRAME;
        second_xs = spots(spots.TRACK_ID == real_tracks(second_spot) & spots.FRAME > 0, :).POSITION_X;
        second_ys = spots(spots.TRACK_ID == real_tracks(second_spot) & spots.FRAME > 0, :).POSITION_Y;
        second_times = spots(spots.TRACK_ID == real_tracks(second_spot) & spots.FRAME > 0, :).FRAME;
        
        % Identify the relevant overlapping times:
        start_time = max(min(first_times), min(second_times));
        end_time = min(max(first_times), max(second_times));
        first_start = find(first_times == start_time, 1);
        first_end = find(first_times == end_time, 1);
        second_start = find(second_times == start_time, 1);
        second_end = find(second_times == end_time, 1);
        relevant_first_xs = first_xs(first_start:first_end);
        relevant_first_ys = first_ys(first_start:first_end);
        relevant_second_xs = second_xs(second_start:second_end);
        relevant_second_ys = second_ys(second_start:second_end);

        % Generate a squared displacement vector; for each overlapping time point, find how far apart the dots are:
        distances = zeros(1, end_time - start_time);
        for timepoint = 1:(end_time - start_time)
            distances(timepoint) = sqrt((relevant_first_xs(timepoint) - relevant_second_xs(timepoint))^2 + ...
                (relevant_first_ys(timepoint) - relevant_second_ys(timepoint))^2);
        end
        
        % MSD it!
        msd_t = zeros(1, max_tau);
        for dt = 1:max_tau
            displacement = distances((1 + dt):end) - distances(1:(end - dt));
            sqrdispl = displacement.^2;
            msd_t(dt) = mean(sqrdispl);
        end
        
        % Tack onto the output frame
        Taus{index} = (1:max_tau)*frametosecond; 
        Pair_IDs{index} = repelem(index, max_tau);
        MSDs{index} = msd_t;
        mean_distance = mean(distances);
        Mean_Dists{index} = repelem(mean_distance, max_tau);
        delt_dist = abs(distances(1) - distances(end));
        Delta_Dists{index} = repelem(delt_dist, max_tau);
        
        % Generate output for each pair
        output{index} = table([Taus{index}]', [MSDs{index}]', [Pair_IDs{index}]', [Mean_Dists{index}]', 'VariableNames', {'Lagtime', 'MSD', 'Pair ID', 'Mean Distance'});
        %writetable(output{index}, sprintf('Corelets_PairwiseMSD_%s_Pair%d.csv', filename, index));
    end
end

end

