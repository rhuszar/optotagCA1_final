function runAssembly_Vrev(index)
% Discover assemblies in the PRE & POST (outside stim) - track expression throughout
% pre and post ; ignore everything on the track
% If there is no behavior, detect / track assemblies throughout the entire
% recording, or in the first 4h
%
% Crosscorrelate the spike times of each neuron with times of assembly
% expression
% 
% Threshold depends on baseline expression of each assembly.


animal_exclude = [{'e16_3m2'} {'e15_13f1'}];

% Add buzcode to path
addpath(genpath('/gpfs/home/rh2618/buzcode'))
addpath(genpath('/gpfs/scratch/rh2618/code'))
% Path to sessions folder on the HPC
datpath = '/gpfs/scratch/rh2618/rh_data';

fwhm = 0.025;    % fwhm of gaussian for smoothing
dt = 0.001;
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
LIM = 120*60;   % 2 hours

index = str2double(index);
% load('basepaths_mac.mat')
% load('basepaths_mac_rev1.mat')
load('basepaths_rerun.mat')
basepath = basepaths_all{ index };

% Basepepath on HPC
path_parts = strsplit(basepath,'/');
basepath_hpc = fullfile(datpath,path_parts{end-2},path_parts{end-1},path_parts{end});
basename = bz_BasenameFromBasepath(basepath_hpc);
disp(basepath_hpc)

if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), animal_exclude) )
    disp('This animal is not being considered, returning...')
    return
end

% Define
st = dbstack;
funcname = st.name;
outpath = fullfile(datpath, [funcname '_OUT']);
if ~exist(outpath, 'dir')
    mkdir(outpath);
end

outfile = fullfile(outpath,[basename '.mat']);
outfile_err = fullfile(outpath,[basename '_ERROR.mat']);
if exist(outfile)
    disp('Already fit')
    return
end

%% Load data

load( fullfile(basepath_hpc, [ basename '.spikes.cellinfo.mat']) )         
load( fullfile(basepath_hpc, [ basename '.cell_metrics.cellinfo.mat']) )
load( fullfile(basepath_hpc, [ basename '.ripples.events.mat']) )
load( fullfile(basepath_hpc, [ basename '.pulses.events.mat']) )
load( fullfile(basepath_hpc, [basename '.SleepState.states.mat']) )
load( fullfile(basepath_hpc, [basename '.session.mat']) );

%% Calculate intervals for spont. homecage depending on whether we recorded multiple sessions
% Across all relevant sessions, we calculate the ripple / theta
% cofiring in the homecage 

% Get the last spike in the recording
reclen = session.timeSeries.adc.nSamples / session.timeSeries.adc.sr;
try
% Behavior
    load( fullfile(basepath_hpc, [basename '.Behavior.mat']) )
    if exist( fullfile(basepath_hpc, [basename '.Tracking.Behavior.mat']) ) 
        % If multiple behavioral sessions, consider homecage data surrounding
        % the first session only
        load( fullfile(basepath_hpc, [basename '.Tracking.Behavior.mat'])  )
        beh_ts = tracking.events.subSessions(1,:);
        discard = tracking.events.subSessions(2,1);
    else
        load( fullfile(basepath_hpc, [basename '.MergePoints.events.mat']) )
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

%%  Obtain the various intervals / indicators into cell types

nonpulse_period = SubtractIntervals( [0 discard], [ pulses.intsPeriods ; beh_ts]);
discard_period = ConsolidateIntervals( [ pulses.intsPeriods ; beh_ts ; [ discard reclen ] ]);

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);

% Intervals of interest for state-dependent expression strenth
% ( Remove any overlapping stim intervals )
ripple_int = ripples.timestamps;
ripple_int = ripple_int( InIntervals(ripple_int(:,1), nonpulse_period) & InIntervals(ripple_int(:,2), nonpulse_period), : );
% Pre NREM
nrem_int = SubtractIntervals( SleepState.ints.NREMstate, discard_period );
% Pre REM
rem_int = SubtractIntervals( SleepState.ints.REMstate, discard_period );

%% Error handling - this code may crash for various reasons (insufficient memory requested, etc.) ; store the error message in this case
try
%% Obtain binned / smoothed spike matrices and run ICA
fprintf('\nGet smoothed matrices in the homecage, outside pulses\n')
spkCnts_smooth = struct('data', [], 'timestamps', [], 'dt', [], 'binsize', []);
for perN = 1:size(nonpulse_period,1)
    spkCnts_smooth(perN) = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  dt, 'binsize', fwhm, 'win', nonpulse_period(perN,:),'bintype','gaussian');
end
spkCnts_smooth_cat = cell2mat( {spkCnts_smooth.data}' );
timestamps_cat = cell2mat( { spkCnts_smooth.timestamps }' );

% Detect assembly expression
fprintf('\nExtract ICs:\n')
tic
assemblyTemplates = assembly_patterns( spkCnts_smooth_cat',opts );
elapsed_t = toc;
fprintf('Duration: %.3f seconds\n', elapsed_t)
% Normalize the vectors
[~, b] = max(abs(assemblyTemplates));          % sign of largest absolute weight connection positive
lininds = sub2ind( size(assemblyTemplates), b, 1:size(assemblyTemplates,2));
sgn = sign(assemblyTemplates(lininds));
assemblyTemplates = assemblyTemplates .* repmat(sgn, [size(assemblyTemplates,1), 1]);
assemblyTemplates = assemblyTemplates ./ sqrt(sum(assemblyTemplates.^2));

% Spontaneous expression strength in homecage
fprintf('\nExpression strength in homecage:\n')
tic
activities_spont = assembly_activity(assemblyTemplates, spkCnts_smooth_cat' );
elapsed_t = toc;
fprintf('Duration: %.3f seconds\n', elapsed_t)

%% Identify peaks in the signal as significant assembly expression 

% Compute assembly activations on the track and in homecage
assembly_act_spont = cell( size(activities_spont,1), 1 );
ccg_pyrAssembly = cell( size(activities_spont,1), 1 );
pks_spont = cell( size(assemblyTemplates,2), 1 );

Rthresh = nan( size(assemblyTemplates,2), 1 );
for kp = 1:size(assemblyTemplates,2)
    
    % Adapt threshold between assemblies
    Rthresh(kp) = 2 * std( activities_spont(kp,:) ) ;
    % Detect assembly peaks
    [pks_spont{kp},locs] = findpeaks( activities_spont(kp,:) );       
    assembly_act_spont{kp} = timestamps_cat( locs( pks_spont{kp} > Rthresh(kp) ) );
    
    % Compute the normalized CCG between assembly expressions on the
    % track, and single unit spikes on the track
    ccg_pyrAssembly{kp} = [];
    for jj = pyrInds
        curr_spk = spikes.times{jj}( InIntervals( spikes.times{jj}, nonpulse_period ) );
        ccg = CCG([ curr_spk ; assembly_act_spont{kp} ], [ ones(size( curr_spk )) ; 2*ones(size( assembly_act_spont{kp} )) ], 'binsize', .002, 'duration', .1 );
        ccg_pyrAssembly{kp} = [ ccg_pyrAssembly{kp} ; ccg(:,1,2)' ./ length(curr_spk) ];
    end
    
end

% Average expression strength within REM / NREM / RIPPLE while in the homecage
% For behavior sessions, consider both pre/post in the case of ripples
try
    rem_spont_strength_mu = nanmean( activities_spont(:, InIntervals( timestamps_cat, rem_int ) ), 2);
    rem_spont_strength_sem = nanstd( activities_spont(:, InIntervals( timestamps_cat, rem_int ) ), [], 2) ./ sqrt( sum( InIntervals( timestamps_cat, rem_int ) ) );
catch
    rem_spont_strength_mu = [];
    rem_spont_strength_sem = [];
end
try
    nrem_spont_strength_mu = mean( activities_spont(:, InIntervals( timestamps_cat, nrem_int ) ), 2);
    nrem_spont_strength_sem = nanstd( activities_spont(:, InIntervals( timestamps_cat, nrem_int ) ),[], 2) ./ sqrt( sum( InIntervals( timestamps_cat, nrem_int ) ) );
catch
    nrem_spont_strength_mu = [];
    nrem_spont_strength_sem = [];
end
try
    ripple_spont_strength_mu = mean( activities_spont(:, InIntervals( timestamps_cat, ripple_int ) ), 2);
    ripple_spont_strength_sem = nanstd( activities_spont(:, InIntervals( timestamps_cat, ripple_int ) ),[], 2) ./ sqrt( sum( InIntervals( timestamps_cat, ripple_int ) ) );
catch
    ripple_spont_strength_mu = [];
    ripple_spont_strength_sem = [];
end

catch err_msg
    save( outfile_err, 'basepath', 'err_msg' )
end

%%

save( outfile, 'basepath', 'basename', 'assemblyTemplates', 'assembly_act_spont','pyrInds', 'pyrInds_opto', 'rem_spont_strength_mu', 'rem_spont_strength_sem', ...
      'nrem_spont_strength_mu', 'nrem_spont_strength_sem', 'ripple_spont_strength_mu', 'ripple_spont_strength_sem', 'ripple_int', ...
      'nrem_int', 'rem_int', 'pks_spont', 'Rthresh', 'nonpulse_period', 'discard_period', 'beh_ts', 'discard', 'LIM', 'ccg_pyrAssembly');

fprintf('Done!\n')


end