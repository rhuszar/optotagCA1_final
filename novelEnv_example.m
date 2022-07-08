



basepaths_famnov = [{'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e13\e13_26m1\e13_26m1_211019'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e13\e13_26m1\e13_26m1_211116'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e13\e13_26m1\e13_26m1_211119'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e13\e13_26m1\e13_26m1_211124'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211210'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211211'} ; ...    % 2 novel
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211212'} ; ...    % 2 novel
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211213'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220117'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220118'} ; ...  % 2 novel
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220119'} ; ...  % 2 novel
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220120'}];
                                                                                                         % #PYR       #tagPYR    proportion
basepaths_2nov =    [{'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211211'} ; ...   % 164.0000   18.0000    0.1098
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211212'} ; ...    % 165.0000   14.0000    0.0848
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220118'} ; ...  % 318.0000   12.0000    0.0377
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220119'}];      % 293.0000   13.0000    0.0444
                
%%
for kp = 1:length(basepaths_famnov)
    getThetaCycles( basepaths_famnov{kp} )
end

%% Calculate number of tagged cells in sessions with multiple novel sessions

vals = [];
for kp = 1:length(basepaths_2nov)
   basepath =  basepaths_2nov{kp};
   basename = bz_BasenameFromBasepath(basepath);
   load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
   
    opto_id_log = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id_log = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id_log & pyr_id_log;
    
    vals = [vals ; sum(pyr_id_log) sum(pyr_opto_id) sum(pyr_opto_id)./sum(pyr_id_log)];
end

%% Session of interest selected !! 
% Manu linear maze conv factor 
% 116.3406

basepath = alterPath( 'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211211', true);
basename = bz_BasenameFromBasepath(basepath);

load(fullfile(basepath, [basename '.Behavior.mat']))
load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
load(fullfile(basepath, [basename '.Tracking.Behavior.mat']))
load(fullfile(basepath, [basename '.pulses.events.mat']))
mapsLin = load(fullfile(basepath, [basename '.firingMapsAvg_multSess.cellinfo.mat']));
mapsOF = load(fullfile(basepath, [basename '.firingMapsAvg2D.cellinfo.mat']));

opto_ind = (cell_metrics.optoTag.p_salt<0.001) & cell_metrics.optoTag.h_reliable;
pyr_ind = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_ind = find( opto_ind & pyr_ind );

% Obtain indices for good trials
% Familiar


%%
pos = [];
recN = [];
raw_pos = {};
speed_threshold = 3;


% Run optotagged cells
cell_set = pyr_opto_ind;
% Run untagged cells
%cell_set = setdiff( find(pyr_ind), pyr_opto_ind );

clear celltune
% Number of cells for a session of interest
celltune(length(cell_set)) = struct('pos', [], 'recN', [], 'uid', []);

% Recording id for the different mazes - 1 is always familiar
recN_lin = unique( behavior.masks.recording )';
recN_maze = setdiff([1 2 3], recN_lin);
x_all_lin = behavior.position.x;
y_all_lin = behavior.position.y;
% Go over all the linear mazes
for recordingN = recN_lin
    
    x_subset = x_all_lin(behavior.masks.recording == recordingN);
    y_subset = y_all_lin(behavior.masks.recording == recordingN);
    raw_pos{recordingN} = [x_subset(1:20000) y_subset(1:20000)];
    
    recordingN_trials = ~isnan( behavior.trials.trial_ints{recordingN}(:,1) );
    visitedArm = behavior.trials.visitedArm(behavior.trials.recordings == recordingN);
    leftArm = visitedArm == 0 & recordingN_trials;
    rightArm = visitedArm == 1 & recordingN_trials;

    % Indicies of left arm trials
    leftArm_inds = arrayfun(@(x) behavior.trials.position_trcat(recordingN).timestamps >= behavior.trials.trial_ints{recordingN}(x,1) & ...
                                 behavior.trials.position_trcat(recordingN).timestamps <= behavior.trials.trial_ints{recordingN}(x,2), find(leftArm)', 'UniformOutput', false);
    % Indices of right arm trials
    rightArm_inds = arrayfun(@(x) behavior.trials.position_trcat(recordingN).timestamps >= behavior.trials.trial_ints{recordingN}(x,1) & ...
                                 behavior.trials.position_trcat(recordingN).timestamps <= behavior.trials.trial_ints{recordingN}(x,2), find(rightArm)', 'UniformOutput', false);
    leftArm = find(leftArm);
    rightArm = find(rightArm);

    for kp = 1:length(cell_set)

        fprintf('%d / %d \n', kp, length(cell_set))
        spk_t = spikes.times{cell_set(kp)};

        % LEFT
        % Get rate maps in individual trials
        leftArm_logical = logical( sum( cell2mat( leftArm_inds ),2));
        % Assembly expression times in given trial
        spk_t_in = spk_t( InIntervals(spk_t, behavior.trials.trial_ints{recordingN}(leftArm,:)) );
        spk_t_in = spk_t_in( interp1( behavior.trials.position_trcat(recordingN).timestamps( leftArm_logical ), ...
                                           behavior.trials.position_trcat(recordingN).v( leftArm_logical ), spk_t_in) > speed_threshold );
        % Associated positions
        xl = interp1( behavior.timestamps( behavior.masks.recording == recordingN ), behavior.position.x( behavior.masks.recording == recordingN ), spk_t_in );
        yl = interp1( behavior.timestamps( behavior.masks.recording == recordingN ), behavior.position.y( behavior.masks.recording == recordingN ), spk_t_in );

        % RIGHT
        % Get rate maps in individual trials
        rightArm_logical = logical( sum( cell2mat( rightArm_inds ), 2) );
        % assembly expression times in given trial
        spk_t_right = spk_t( InIntervals(spk_t, behavior.trials.trial_ints{recordingN}(rightArm,:)) );
        spk_t_right = spk_t_right( interp1( behavior.trials.position_trcat(recordingN).timestamps( rightArm_logical ), ...
                                               behavior.trials.position_trcat(recordingN).v( rightArm_logical ), spk_t_right) > speed_threshold );
        % associated positions
        % Associated positions
        xr = interp1( behavior.timestamps( behavior.masks.recording == recordingN ), behavior.position.x( behavior.masks.recording == recordingN ), spk_t_right );
        yr = interp1( behavior.timestamps( behavior.masks.recording == recordingN ), behavior.position.y( behavior.masks.recording == recordingN ), spk_t_right );

    %     figure
    %     h = plot(behavior.position.x(1:20000), behavior.position.y(1:20000, 1), '-b'); hold on
    %     h.Color(4) = 0.1;

        celltune(kp).pos = [celltune(kp).pos ; [[ xl ; xr ] [ yl ; yr ]]];
        % celltune(kp).recN = [celltune(kp).recN ; repmat(recordingN, size([ xl ; xr ]))];
        celltune(kp).recN = [celltune(kp).recN ; repmat(recordingN*10+1, size( xl )) ; repmat(recordingN*10+2, size( xr )) ];
         celltune(kp).uid = cell_set(kp);
    %     hold on
    %     plot([ xl ; xr ], [ yl ; yr ], '.r')
    %     uiwait

    end
end

x = tracking.position.x( tracking.events.subSessionsMask == recN_maze );
y = tracking.position.y( tracking.events.subSessionsMask == recN_maze );
tvec = tracking.timestamps( tracking.events.subSessionsMask == recN_maze );

raw_pos{recN_maze} = [x y ];

% Get the velocity
[~,~,~,vx,vy,~,~] = trajectory_kalman_filter(x, y, tvec,2);
v = [vx vy]';
v = arrayfun(@(x) norm(v(:,x)), 1:length(v)); 

xedges = 0:1:ceil( max(x) );
yedges = 0:1:ceil( max(y) );
camsr = 25;
for kp = 1:length(cell_set)

    fprintf('%d / %d \n', kp, length(cell_set))
    spk_t = spikes.times{cell_set(kp)};

    spkt_maze = histcounts(spk_t, [ tvec ; tvec(end)+( 1/camsr )]);
    % Get rid of spikes during immobility
    spkt_maze( v < speed_threshold ) = 0;

    spk_pos = arrayfun(@(ip) [repmat(x(ip), spkt_maze(ip), 1) repmat(y(ip), spkt_maze(ip), 1)], find( spkt_maze ), 'UniformOutput', false );
    spk_pos = cell2mat( spk_pos' );
    
    celltune(kp).pos = [celltune(kp).pos ; spk_pos];
    celltune(kp).recN = [celltune(kp).recN ; repmat(recN_maze, size(spk_pos, 1), 1)];

end

%%

celltune_all = celltune;


%%
xlims = [min(behavior.position.x) max(behavior.position.x)];
ylims = [min(behavior.position.y) max(behavior.position.y)];
close all
outpath = fullfile(basepath, [basename '_tagRatemaps']);
if ~exist(outpath, 'dir')
    mkdir( outpath )
end

ePulses = pulses.timestamps( pulses.duration < 0.006, :);
ePulse_on = ePulses(:,1);
% Some parameters for plotting the rasters
interval = [0 round( median( diff( ePulses, [], 2) ), 3)];
lm = [1 size(ePulses,1)];


for kp = 1:length( celltune )
    close all
%     
%     
%     figure
%     [ ccgs, t ] = CCG([ePulse_on ; spikes.times{pyr_opto_ind(kp)}], [ones(size(ePulse_on)) ; 2*ones(size(spikes.times{pyr_opto_ind(kp)}))], 'duration', 0.05, 'binsize', .0008);
%     spk_t = spikes.times{pyr_opto_ind(kp)}( ~InIntervals(  spikes.times{pyr_opto_ind(kp)}, pulses.intsPeriods )  );
%     [ acgs, t1 ] = CCG(spk_t, ones(size(spk_t)), 'duration', 0.1, 'binsize', .0008);
%     set(gcf,'Position', [ 657         119        1139         389])
%     set(gcf,'Position', [656   445   951   472])
%     
%     subplot(1,3,1)
%     plot(t1, acgs./ length(spk_t) ./ .0008, '-b', 'linewidth', 1.5)
%     xlabel('Time (s)', 'fontsize', 14)
%     ylabel('Firing rate (Hz)', 'fontsize', 14)
%     title('Autocorrelation', 'fontsize', 14)
%     
%     subplot(1,3,2)
%     plot(cell_metrics.optoTag.tvec_raster, cell_metrics.optoTag.pulse_raster{ pyr_opto_ind(kp) }', '.k'); hold on
%     ylim(lm)
%     patch([interval flip(interval) ],[lm(1) lm(1) lm(2) lm(2)],'b', 'FaceAlpha', 0.2, 'linestyle', 'none')
%     xlim([-.01 0.02])
%     box off
% 
%     subplot(1,3,3)
%     plot(t, ccgs(:,1,2) ./ length(ePulse_on) ./ .0008, '-b', 'linewidth', 1.5)
%     xlabel('Time (s)', 'fontsize', 14)
%     ylabel('Firingrate', 'fontsize', 14)
%     title('Raster', 'fontsize', 14)
%     xlim([-.01 0.02])
    
    
    h = {};
    figure
    % set(gcf,'Position',[444         353        1348         712])
    set(gcf,'Position',[444         686        1348         379])
    % Familiar
    % subplot(2,6,[1 2])
    subplot(1,3,1)
    recordingN = recN_lin(1);
    h{1} = plot(raw_pos{recordingN}(:,1), raw_pos{recordingN}(:,2), '-b'); hold on
    h{1}.Color(4) = 0.1;
    plot( celltune(kp).pos(celltune(kp).recN==recordingN,1), celltune(kp).pos(celltune(kp).recN==recordingN,2), '.r', 'MarkerSize', 3 );
    xlim(xlims)
    ylim(ylims)
    xticks([0 50])
    %axis off
%     subplot(2,6,7)
%     imagesc( mapsLin.firingMaps(recordingN).rateMaps_trial{pyr_opto_ind(kp),1} )
%     subplot(2,6,8)
%     imagesc( mapsLin.firingMaps(recordingN).rateMaps_trial{pyr_opto_ind(kp),2} )
    
    
    % Novel 1
    recordingN = recN_lin(2);
    % subplot(2,6,[3 4])
    subplot(1,3,2)
    h{2} = plot(raw_pos{recordingN}(:,1), raw_pos{recordingN}(:,2), '-b'); hold on
    h{2}.Color(4) = 0.1;
    plot( celltune(kp).pos(celltune(kp).recN==recordingN,1), celltune(kp).pos(celltune(kp).recN==recordingN,2), '.r', 'MarkerSize', 3 );
    xlim(xlims)
    ylim(ylims)
    axis off
    
    % Novel 2
        recordingN = recN_maze;
    % subplot(2,6,[5 6])
    subplot(1,3,3)
    h{3} = plot(raw_pos{recordingN}(:,1), raw_pos{recordingN}(:,2), '-b'); hold on
    h{3}.Color(4) = 0.1;
    plot( celltune(kp).pos(celltune(kp).recN==recordingN,1), celltune(kp).pos(celltune(kp).recN==recordingN,2), '.r', 'MarkerSize', 3 );
    xlim(xlims)
    ylim(ylims)
    axis off
    outfile = fullfile(outpath, sprintf('pyrtag_%d.png', pyr_opto_ind(kp)));
%     saveas(gcf,outfile)
    
    uiwait
end

%% Plot all pairs of tagged cells in given session
%  UPDATE 3/31 - do this for both novel environment sessions 
%  This won't work in sessions with a single novel exploration
binsize = .0008;
duration = 0.05;
pairs = nchoosek(1:length(celltune), 2);
max_spkt = max( cellfun( @max, spikes.times) );
% Recording duration minus pulse intervals
recDur = max_spkt - sum( diff( pulses.intsPeriods') );
recN_linNov = min( recN_lin );
msk = find( tracking.events.subSessionsMask == recordingN );

% ints = [tracking.timestamps(msk(1)) tracking.timestamps(msk(end))];
% recDur = diff(ints');
for kp = 112:length(pairs)    % 112 is the one we actually chose
    
    
    close all
    
    spk1 = spikes.times{pyr_opto_ind(pairs(kp,1))}; spk1 = spk1( ~InIntervals(spk1, pulses.intsPeriods ) );
    spk2 = spikes.times{pyr_opto_ind(pairs(kp,2))}; spk2 = spk2( ~InIntervals(spk2, pulses.intsPeriods ) );
    [ ccg, t ] = CCG([spk1 ; spk2], [ones(size(spk1)) ; 2.*ones(size(spk2))], 'binsize', binsize, 'duration', duration);

    fprintf('%d\n', kp)
    figure(1)
    set(gcf,'Position', [319   699   625   399])
    subplot(2,2,1)
    plot(t, ccg(:,1,1) ./ length(spk1) ./ binsize, 'r' )
    title(sprintf('%.2f', length(spk1) ./ recDur))
    subplot(2,2,4)
    plot(t, ccg(:,2,2) ./ length(spk2) ./ binsize, 'g' )
    title(sprintf('%.2f', length(spk2) ./ recDur))
    subplot(2,2,2)
    plot(t, smooth( ccg(:,1,2) ./ length(spk1) ./ binsize ), 'k' )
        
    % Maze explorations
    
    figure(2)
    set(gcf,'Position',[952 630 1365 638])
    % Familir left
    subplot(2,3,1)
    h1 = plot(raw_pos{1}(:,2), raw_pos{1}(:,1), '-b'); hold on
    h1.Color(4) = 0.1;
    % Plot the left
    plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==11,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==11,1), '.r', 'MarkerSize', 6 );
    hold on
    plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==11,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==11,1), '.g', 'MarkerSize', 6 );
    % Familiar right
    subplot(2,3,4)
    h2 = plot(raw_pos{1}(:,2), raw_pos{1}(:,1), '-b'); hold on
    h2.Color(4) = 0.1;
    plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==12,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==12,1), '.r', 'MarkerSize', 6 );
    hold on
    plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==12,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==12,1), '.g', 'MarkerSize', 6 );

    % First novel
    subplot(2,3,2)
    h2 = plot(raw_pos{2}(:,2), raw_pos{2}(:,1), '-b'); hold on
    h2.Color(4) = 0.1;
    if ismember( 21, unique( celltune(pairs(kp,1)).recN ) )   % Is this a novel linear track ? 
        % Plot left
        plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==21,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==21,1), '.r', 'MarkerSize', 6 );
        hold on
        plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==21,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==21,1), '.g', 'MarkerSize', 6 );
        % Plot right
        subplot(2,3,5)
        hx = plot(raw_pos{1}(:,2), raw_pos{1}(:,1), '-b'); hold on
        hx.Color(4) = 0.1;
        plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==22,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==22,1), '.r', 'MarkerSize', 6 );
        hold on
        plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==22,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==22,1), '.g', 'MarkerSize', 6 );
    else
        plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==2,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==2,1), '.r', 'MarkerSize', 6 );
        hold on
        plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==2,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==2,1), '.g', 'MarkerSize', 6 );
    end

    % Second novel

    if ismember( 31, unique( celltune(pairs(kp,1)).recN ) )   % Is this a novel linear track ?
        subplot(2,3,3)
        h3 = plot(raw_pos{3}(:,2), raw_pos{3}(:,1), '-b'); hold on
        h3.Color(4) = 0.1;
        % Plot left
        plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==31,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==31,1), '.r', 'MarkerSize', 6 );
        hold on
        plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==31,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==31,1), '.g', 'MarkerSize', 6 );
        % Plot right
        subplot(2,3,6)
        hy = plot(raw_pos{3}(:,2), raw_pos{3}(:,1), '-b'); hold on
        hy.Color(4) = 0.1;
        plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==32,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==32,1), '.r', 'MarkerSize', 6 );
        hold on
        plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==32,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==32,1), '.g', 'MarkerSize', 6 );
    elseif ismember( 3, unique( celltune(pairs(kp,1)).recN ) )
        subplot(2,3,3)
        h3 = plot(raw_pos{3}(:,2), raw_pos{3}(:,1), '-b'); hold on
        h3.Color(4) = 0.1;
        plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==3,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==3,1), '.r', 'MarkerSize', 6 );
        hold on
        plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==3,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==3,1), '.g', 'MarkerSize', 6 );
    end
    uiwait
    
end

%% For a fixed pair, loop over untagged cells to find an example of something that overlaps less well to illustrate the example


pair_id = 112;

max_spkt = max( cellfun( @max, spikes.times) );
% Recording duration minus pulse intervals
recDur = max_spkt - sum( diff( pulses.intsPeriods') );
recordingN = recN_maze;

% Reextract spike trains
spk_red = spikes.times{pyr_opto_ind(pairs(pair_id,1))}; spk_red = spk_red( ~InIntervals(spk_red, pulses.intsPeriods ) );
spk_green = spikes.times{pyr_opto_ind(pairs(pair_id,2))}; spk_green = spk_green( ~InIntervals(spk_green, pulses.intsPeriods ) );

% Loop over all remaining cells
for jp = 70:length( celltune_all )
    
    fprintf('UID: %d\n', celltune_all(jp).uid)

    % Display temporal cofiring patterns 
    spk_t = spikes.times{ celltune_all(jp).uid }; spk_t = spk_t( ~InIntervals(spk_t, pulses.intsPeriods ) );
    [ ccg, t ] = CCG([spk_red ; spk_green ; spk_t], [ones(size(spk_red)) ; 2.*ones(size(spk_green)) ; 3.*ones(size(spk_t))], 'binsize', binsize, 'duration', duration);
    
    figure(1)
    set(gcf,'Position', [183   567   700   433])
    % Red
    subplot(3,3,1)
    plot(t, ccg(:,1,1) ./ length(spk_red) ./ binsize, 'r' )
    title(sprintf('%.2f', length(spk_red) ./ recDur))
    xlim([ min( t ) max( t )])
    % Green
    subplot(3,3,5)
    plot(t, ccg(:,2,2) ./ length(spk_green) ./ binsize, 'g' )
    title(sprintf('%.2f', length(spk_green) ./ recDur))
    xlim([ min( t ) max( t )])
    % Blue
    subplot(3,3,9)
    plot(t, ccg(:,3,3) ./ length(spk_t) ./ binsize, 'b' )
    title(sprintf('%.2f', length(spk_t) ./ recDur))
    xlim([ min( t ) max( t )])
    
    h_cc = {};
    % Show the cross correlograms
    h_cc{1} = subplot(3,3,2);
    plot(t, smooth( ccg(:,1,2) ./ length(spk_red) ./ binsize ), 'k' );
    xlim([ min( t ) max( t )])
    h_cc{2} = subplot(3,3,3);
    plot(t, smooth( ccg(:,1,3) ./ length(spk_red) ./ binsize ), 'k' );
    xlim([ min( t ) max( t )])
    h_cc{3} = subplot(3,3,6);
    plot(t, smooth( ccg(:,2,3) ./ length(spk_red) ./ binsize ), 'k' );
    xlim([ min( t ) max( t )])
    ylim(h_cc{2}, h_cc{1}.YLim);
    ylim(h_cc{3}, h_cc{1}.YLim);
    
    % Show the raw spikes plotted on top of the position
    figure(2)
    
    % Familir left
    subplot(2,3,1)
    % Plot the left
    hold on
    v1 = plot( celltune_all(jp).pos( celltune_all(jp).recN==11,2 ), celltune_all(jp).pos( celltune_all(jp).recN==11,1), '.k', 'MarkerSize', 6 );
    % Plot the  right
    subplot(2,3,4)
    hold on
    v2 = plot( celltune_all(jp).pos( celltune_all(jp).recN==12,2 ), celltune_all(jp).pos( celltune_all(jp).recN==12,1), '.k', 'MarkerSize', 6 );

    % First novel
    subplot(2,3,2)
    if ismember( 21, unique( celltune_all(jp).recN ) )   % Is this a novel linear track ? 
        % Plot left
        hold on
        v3 = plot( celltune_all(jp).pos( celltune_all(jp).recN==21,2), celltune_all(jp).pos( celltune_all(jp).recN==21,1), '.k', 'MarkerSize', 6 );
        % Plot right
        subplot(2,3,5)
        hold on
        v4 = plot( celltune_all(jp).pos( celltune_all(jp).recN==22,2), celltune_all(jp).pos( celltune_all(jp).recN==22,1), '.k', 'MarkerSize', 6 );
    else
        v3 = plot( celltune_all(jp).pos( celltune_all(jp).recN==2,2), celltune_all(jp).pos( celltune_all(jp).recN==2,1), '.k', 'MarkerSize', 6 );
    end

    % Second novel
    if ismember( 31, unique( celltune_all(jp).recN ) )   % Is this a novel linear track ?
        subplot(2,3,3)
        % Plot left
        hold on
        v5 =  plot( celltune_all(jp).pos( celltune_all(jp).recN==31,2), celltune_all(jp).pos( celltune_all(jp).recN==31,1), '.k', 'MarkerSize', 6 );
        % Plot right
        subplot(2,3,6)
        hold on
        v6 = plot( celltune_all(jp).pos( celltune_all(jp).recN==32,2), celltune_all(jp).pos( celltune_all(jp).recN==32,1), '.k', 'MarkerSize', 6 );
    elseif ismember( 3, unique( celltune_all(jp).recN ) )
        subplot(2,3,3)
        hold on
        v5 = plot( celltune_all(jp).pos( celltune_all(jp).recN==3,2), celltune_all(jp).pos( celltune_all(jp).recN==3,1), '.k', 'MarkerSize', 6 );
    end
    
    figure(1)
    shg
    uiwait
    if exist('v1', 'var'); delete(v1); end
    if exist('v2', 'var'); delete(v2); end
    if exist('v3', 'var'); delete(v3); end
    if exist('v4', 'var'); delete(v4); end
    if exist('v5', 'var'); delete(v5); end
    if exist('v6', 'var'); delete(v6); end
        
%     f = {};
%     f{1} = plot(raw_pos{recordingN}(:,2), raw_pos{recordingN}(:,1), '-b'); hold on
%     f{1}.Color(4) = 0.1;
%     plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==recordingN,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==recordingN,1), '.r', 'MarkerSize', 6 );
%     hold on
%     plot( celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==recordingN,2), celltune(pairs(kp,2)).pos(celltune(pairs(kp,2)).recN==recordingN,1), '.g', 'MarkerSize', 6 );
%     
%     spkt_maze = histcounts(spk_t, [ tvec ; tvec(end)+( 1/camsr )]);
%     % Get rid of spikes during immobility
%     spkt_maze( v < speed_threshold ) = 0;
% 
%     spk_pos = arrayfun(@(ip) [repmat(x(ip), spkt_maze(ip), 1) repmat(y(ip), spkt_maze(ip), 1)], find( spkt_maze ), 'UniformOutput', false );
%     spk_pos = cell2mat( spk_pos' );
%     subplot(1,2,2)
%     f{2} = plot(raw_pos{recordingN}(:,2), raw_pos{recordingN}(:,1), '-b'); hold on
%     f{2}.Color(4) = 0.1;
%     plot( celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==recordingN,2), celltune(pairs(kp,1)).pos(celltune(pairs(kp,1)).recN==recordingN,1), '.r', 'MarkerSize', 6 );
%     hold on
%     plot( spk_pos(:,2), spk_pos(:,1), '.k', 'MarkerSize', 6 );
%     uiwait
    
end
    
%% manu's linear maze example lot

n = 184;

x1 = 1:110;
x2 = 110:-1:1;

t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1, firingMaps(2).rateMaps{n}{2})         % outbound
xlim([1 max(x1)])
ax1.XTick = x1(10:10:end);

ax2 = axes(t);
plot(ax2, flip( firingMaps(2).rateMaps{n}{1} ), 'r' )   % inbound
xlim([1 max(x1)])
ax2.XTick = [ 1 10:10:100 ];
ax2.XTickLabel = arrayfun(@(x ) num2str(x), 110:-10:10 , 'UniformOutput', false);
ax2.XAxisLocation = 'top';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';

%% Thomas' linear maze example lot

%%

x1 = 0:0.1:40;
y1 = 4.*cos(x1)./(x1+2);
x2 = 1:0.2:20;
y2 = x2.^2./x2.^3;

t = tiledlayout(1,1);
ax1 = axes(t);
plot(ax1,x1,y1,'-r')

ax2 = axes(t);
plot(ax2,x2,y2,'-k')
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
ax1.Box = 'off';
ax2.Box = 'off';


