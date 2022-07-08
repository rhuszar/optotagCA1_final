



nrand = 100;
fwhm = 0.025;    % fwhm of gaussian for smoothing
dt = 0.001;
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
LIM = 120*60;   % 2 hours

basepath = 'C:\Users\Roman\Google Drive\buzz_workspace\optotag\DATA\e13\e13_16f1\e13_16f1_210328';
basename = bz_BasenameFromBasepath(basepath);


load( fullfile(basepath, [ basename '.spikes.cellinfo.mat']) )         
load( fullfile(basepath, [ basename '.cell_metrics.cellinfo.mat']) )
load( fullfile(basepath, [ basename '.ripples.events.mat']) )
load( fullfile(basepath, [ basename '.pulses.events.mat']) )
load( fullfile(basepath, [basename '.session.mat']) );



% Get the last spike in the recording
reclen = session.timeSeries.adc.nSamples / session.timeSeries.adc.sr;
try
% Behavior
    load( fullfile(basepath, [basename '.Behavior.mat']) )
    if exist( fullfile(basepath, [basename '.Tracking.Behavior.mat']) ) 
        % If multiple behavioral sessions, consider homecage data surrounding
        % the first session only
        load( fullfile(basepath, [basename '.Tracking.Behavior.mat'])  )
        beh_ts = tracking.events.subSessions(1,:);
        discard = tracking.events.subSessions(2,1);
    else
        load( fullfile(basepath, [basename '.MergePoints.events.mat']) )
        [~, ff] = InIntervals( behavior.timestamps(1), MergePoints.timestamps  );
        beh_ts = MergePoints.timestamps(ff,:);
        discard = reclen;
    end
    disp('Behavior session!')
catch
% No behavior
    beh_ts = [];
    if LIM > reclen
        discard = reclen;
    else
        discard = LIM;
    end
    disp('No behavior session!')
end

%%


nonpulse_period = SubtractIntervals( [0 discard], [ pulses.intsPeriods ; beh_ts]);

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);


% Indices into pyrInds which are tagged
pyrInds_opto = arrayfun(@(x) find( x == pyrInds ), pyrInds_opto );
pyrInds_nopto = setdiff( 1:length(pyrInds), pyrInds_opto);

% Intervals of interest for state-dependent expression strenth
% ( Remove any overlapping stim intervals )
ripple_int = ripples.timestamps;
ripple_int = ripple_int( InIntervals(ripple_int(:,1), nonpulse_period) & InIntervals(ripple_int(:,2), nonpulse_period), : );
ripple_int_pad = [ ripple_int(:,1)-0.2 ripple_int(:,2)+0.05 ];


%%

fprintf('\nGet smoothed matrices in the homecage, outside pulses\n')
spkCnts_smooth = struct('data', [], 'timestamps', [], 'dt', [], 'binsize', []);
ripple_indicator = struct('indicator', []);
for perN = 1:size(ripple_int_pad,1)
    spkCnts_smooth(perN) = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  dt, 'binsize', fwhm, 'win', ripple_int_pad(perN,:),'bintype','gaussian');
    ripple_indicator(perN).indicator = perN.*ones(size( spkCnts_smooth(perN).data,1 ), 1);
end
spkCnts_smooth_cat = cell2mat( {spkCnts_smooth.data}' );
ripple_indicator_cat = cell2mat( {ripple_indicator.indicator}' );
timestamps_cat = cell2mat( { spkCnts_smooth.timestamps }' );

% Define an assembly consisting of all tagged neurons
assembly_tag = zeros(length(pyrInds),1); assembly_tag(pyrInds_opto) = 1; 
assembly_tag = assembly_tag ./ norm(assembly_tag);

% Spontaneous expression around SPW-Rs
activities_spont_tag = assembly_activity(assembly_tag, spkCnts_smooth_cat' );

%% Find a ripple we like

lfp = ripples.detectorinfo.detectionparms.lfp;
ts = 0:(1/1250):( length(lfp)*(1/1250) - (1/1250) );
for kp = 1:length(spkCnts_smooth)
    
   tv = [ spkCnts_smooth(kp).timestamps(1) spkCnts_smooth(kp).timestamps(end) ];
   indicator = InIntervals(ts, tv);
   
   assembly_rand_all = {};
    activities_rand_all = {};
    for jp = 1:nrand
        % Randomly generate assembly template
        assembly_rand = zeros(length(pyrInds),1); assembly_rand( pyrInds_nopto( randi( length(pyrInds_nopto), length(pyrInds_opto), 1 ) ) ) = 1; 
        assembly_rand_all{jp} = assembly_rand ./ norm(assembly_rand);
        activities_rand_all{jp} = assembly_activity(assembly_rand, spkCnts_smooth_cat(ripple_indicator_cat==kp,:)' );
    end

%    subplot(2,1,1)
%    plot(ts(indicator), lfp(indicator))
%    subplot(2,1,2)
%    plot(timestamps_cat( ripple_indicator_cat == kp ), activities_spont_tag( ripple_indicator_cat == kp) )
%    uiwait
 %  close all
    set(gcf,'Position', [1000         645         560         693])
    subplot(6,1,[1 2])
    plot(ts(indicator), lfp(indicator))
    subplot(6,1,3)
    plot(timestamps_cat( ripple_indicator_cat == kp ), activities_spont_tag( ripple_indicator_cat == kp) )
    %ylim([-1 10.5])
    subplot(6,1,4)
    plot(timestamps_cat( ripple_indicator_cat == kp ), activities_rand_all{1} )
    %ylim([-1 10.5])
    subplot(6,1,5)
    plot(timestamps_cat( ripple_indicator_cat == kp ), activities_rand_all{2} )
    %ylim([-1 10.5])
    subplot(6,1,6)
    plot(timestamps_cat( ripple_indicator_cat == kp ), activities_rand_all{3} )
    %ylim([-1 10.5])
    shg
    uiwait

end

%%
aa = zeros(1, length(assembly_rand_all{3}));
aa(1:length(pyrInds_opto)) = assembly_rand_all{3}(pyrInds_opto);
aa(length(pyrInds_opto)+1:end) = assembly_rand_all{3}(pyrInds_nopto);
set(gcf,'Position',[160         810        2144         173])
stem(aa, 'filled')


%% 26
rip_index = kp;

nrand = 3;
assembly_rand_all = {};
activities_rand_all = {};
for jp = 1:nrand
% Randomly generate assembly template
assembly_rand = zeros(length(pyrInds),1); assembly_rand( pyrInds_nopto( randi( length(pyrInds_nopto), length(pyrInds_opto), 1 ) ) ) = 1; 
assembly_rand_all{jp} = assembly_rand ./ norm(assembly_rand);
activities_rand_all{jp} = assembly_activity(assembly_rand, spkCnts_smooth_cat(ripple_indicator_cat==rip_index,:)' );
end

tv = [ spkCnts_smooth(rip_index).timestamps(1) spkCnts_smooth(rip_index).timestamps(end) ];
indicator = InIntervals(ts, tv);

close all
set(gcf,'Position', [1000         645         560         693])
subplot(6,1,[1 2])
plot(ts(indicator), lfp(indicator))
subplot(6,1,3)
plot(timestamps_cat( ripple_indicator_cat == rip_index ), activities_spont_tag( ripple_indicator_cat == rip_index) )
ylim([-1 10.5])
subplot(6,1,4)
plot(timestamps_cat( ripple_indicator_cat == rip_index ), activities_rand_all{1} )
ylim([-1 10.5])
subplot(6,1,5)
plot(timestamps_cat( ripple_indicator_cat == rip_index ), activities_rand_all{2} )
ylim([-1 10.5])
subplot(6,1,6)
plot(timestamps_cat( ripple_indicator_cat == rip_index ), activities_rand_all{3} )
ylim([-1 10.5])
shg


