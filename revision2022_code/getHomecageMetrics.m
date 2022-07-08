


function getHomecageMetrics(basepath)

% SWR-related firing rate and participation probability
% Spike transmission probability as a function of ISI
% Spike transmission probability baseline
% Proportion of PYRs tagged
% Fraction of REM shifters
% Burst index
% Theta depth modulation 

basename = bz_BasenameFromBasepath(basepath); 
disp(basename)
load( fullfile(basepath, [basename '.pulses.events.mat']) )
load( fullfile(basepath, [basename '.spikes.cellinfo.mat']) )
load( fullfile(basepath, [ basename '.cell_metrics.cellinfo.mat']) )
load( fullfile(basepath, [basename '.ripples.events.mat']) )
load( fullfile(basepath, [basename '.mono_res.cellinfo.mat']) )
load( fullfile(basepath, [basename '.deepSuperficialfromRipple.channelinfo.mat']) )
load( fullfile(basepath, [basename '.ThetaModulation.cellinfo.mat']) );
load( fullfile(basepath, [basename '.session.mat']) );

% Recording length, in seconds
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
    discard = reclen;
    disp('No behavior session!')
end

% Get non pulse period as interval and as duration of time
nonpulse_period = SubtractIntervals( [0 discard], [ pulses.intsPeriods ; beh_ts]);
nonpulse_period_t = sum( diff(nonpulse_period') );

% Get ripples outside pulses, occurring in the homecage
ripple_ts = ripples.timestamps;
ripple_ts( ripple_ts(:,1) > discard,: ) = [];
rip_in_pulse = InIntervals( ripple_ts(:,1), [ pulses.intsPeriods ; beh_ts ] ) | InIntervals( ripple_ts(:,2), [ pulses.intsPeriods ; beh_ts ] );
ripple_ts(rip_in_pulse,:) = [];

% Tagged pyramidal cells
opto_ind = (cell_metrics.optoTag.p_salt<0.001) & cell_metrics.optoTag.h_reliable;
pyr_ind = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_ind_tagged = find( opto_ind & pyr_ind );
n_int = sum( cellfun(@(x) ~isempty(regexp(x, 'Inter', 'once')), cell_metrics.putativeCellType) );
% Store proportion tagged
prop_tagged = length(pyr_ind_tagged) / sum(pyr_ind);

% Preallocate variables we wish to store
fr = nan(length(pyr_ind_tagged), 1);
burst =  nan(length(pyr_ind_tagged), 1);
div_index = nan(length(pyr_ind_tagged), 1);
ripple_pp =  nan(length(pyr_ind_tagged), 1);
ripple_fr =  nan(length(pyr_ind_tagged), 1);
radial_depth = deepSuperficialfromRipple.channelDistance_estimated( spikes.maxWaveformCh1(pyr_ind_tagged) );
theta_phase_wake = ThetaMod.th_phaselock_wake.phasestats.m( pyr_ind_tagged )';
theta_mDepth_wake = ThetaMod.th_phaselock_wake.phasestats.r( pyr_ind_tagged )';
theta_sig_wake = ThetaMod.th_phaselock_wake.phasestats.p( pyr_ind_tagged )';

% Take into account sessions that don't have REM sleep
if isempty( ThetaMod.th_phaselock_rem.phasedistros )
    theta_phase_rem = nan(length( pyr_ind_tagged ), 1);
    theta_mDepth_rem = nan(length( pyr_ind_tagged ), 1);
    theta_sig_rem = nan(length( pyr_ind_tagged ), 1);
else
    theta_phase_rem = ThetaMod.th_phaselock_rem.phasestats.m( pyr_ind_tagged )';
    theta_mDepth_rem = ThetaMod.th_phaselock_rem.phasestats.r( pyr_ind_tagged )';
    theta_sig_rem = ThetaMod.th_phaselock_rem.phasestats.p( pyr_ind_tagged )';
end

for kp = 1:length( pyr_ind_tagged )
    spk = spikes.times{pyr_ind_tagged( kp) };
    % Exclude spikes outside of interval of consideration
    spk = spk(InIntervals(spk, nonpulse_period));
    % Firing rate
    fr(kp) = length( spk ) / nonpulse_period_t;
    % Burst index
    [ acg, t] = CCG(spk, ones(size(spk)), 'duration', 0.6, 'binsize', 0.001);
    burst(kp) = sum( acg(t >= 0.003 & t <= 0.005) ) / sum( acg(t >= 0.2 & t <= 0.3) );
    % PYR-INT divergence index - fraction of all interneurons you project to
    div_index(kp) = sum( mono_res.pyr2int(:,1) == pyr_ind_tagged(kp) ) ./ n_int;

    [status, intt ] = InIntervals( spk, ripple_ts);                               
    rip_p = unique( intt ); rip_p(rip_p == 0) = [];
    ripple_pp(kp) = length( rip_p ) / size(ripple_ts,1);  
    ripple_fr(kp) = sum( status ) / sum( diff(ripple_ts,[],2) );

end

save(fullfile( basepath, 'HomecageMetrics.mat'), 'fr', 'burst', 'div_index', 'ripple_pp', 'ripple_fr', 'radial_depth', ...
                'theta_phase_wake', 'theta_mDepth_wake', 'theta_sig_wake', 'theta_phase_rem', 'theta_mDepth_rem', 'theta_sig_rem', 'prop_tagged');
    
 
    
end