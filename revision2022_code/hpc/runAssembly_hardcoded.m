function runAssembly_hardcoded(index)
% For a given session, hardcode an assembly comprising all N tagged cells,
% and track its expression around SPW-Rs ; as control, make
% assemblies comprising of random subsets of N unresponsive neurons, and
% track their expression around SPW-Rs. 
% Ripple rate of our assembly should be expressed as a z-score 


animal_exclude = [{'e16_3m2'} {'e15_13f1'}];

% Add buzcode to path
addpath(genpath('/gpfs/home/rh2618/buzcode'))
addpath(genpath('/gpfs/scratch/rh2618/code'))
% Path to sessions folder on the HPC
datpath = '/gpfs/scratch/rh2618/rh_data';

nrand = 100;
fwhm = 0.025;    % fwhm of gaussian for smoothing
dt = 0.001;
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;
LIM = 120*60;   % 2 hours

index = str2double(index);
load('basepaths_mac_rev1.mat')
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

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);

% Only include sessions with >3 tagged neurons
if length( pyrInds_opto ) < 3
    fprintf('\nInsufficient number of tagged neuorons\n')
    return
end
    

% Indices into pyrInds which are tagged
pyrInds_opto = arrayfun(@(x) find( x == pyrInds ), pyrInds_opto );
pyrInds_nopto = setdiff( 1:length(pyrInds), pyrInds_opto);

% Intervals of interest for state-dependent expression strenth
% ( Remove any overlapping stim intervals )
ripple_int = ripples.timestamps;
ripple_int = ripple_int( InIntervals(ripple_int(:,1), nonpulse_period) & InIntervals(ripple_int(:,2), nonpulse_period), : );
ripple_int_pad = [ ripple_int(:,1)-0.2 ripple_int(:,2)+0.05 ];

%% Error handling - this code may crash for various reasons (insufficient memory requested, etc.) ; store the error message in this case
try
%% Obtain binned / smoothed spike matrices and run ICA
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

% Find peaks, and get rid of edge peaks
[pks_tag,locs_tag] = findpeaks( activities_spont_tag ); 
nonedge_locs = 3 == sum( ripple_indicator_cat( [locs_tag'-1 locs_tag' locs_tag'+1 ] ) == ...
                                repmat(ripple_indicator_cat( locs_tag' ), 1, 3), 2);
locs_tag(~nonedge_locs) = [];
pks_tag(~nonedge_locs) = [];
    

% Define assemblies from random subsets of untagged cells - size of subset
% equals number of tagged neurons ; for each randomly defined assembly,
% track expression, and identify peaks
locs_rand_all = cell(1,nrand);
pks_rand_all = cell(1,nrand);
mus = nan(1,nrand);
sds = nan(1,nrand);
for jp = 1:nrand
    fprintf('%d/%d\n', jp, nrand)
    tic
    % Randomly generate assembly template
    assembly_rand = zeros(length(pyrInds),1); assembly_rand( pyrInds_nopto( randi( length(pyrInds_nopto), length(pyrInds_opto), 1 ) ) ) = 1; 
    assembly_rand = assembly_rand ./ norm(assembly_rand);
    % Assembly activity
    activities_rand = assembly_activity(assembly_rand, spkCnts_smooth_cat' );
    % Find activation peaks and exclude edges
    [pks_rand,locs_rand] = findpeaks( activities_rand ); 
    nonedge_locs = 3 == sum( ripple_indicator_cat( [locs_rand'-1 locs_rand' locs_rand'+1 ] ) == ...
                                repmat(ripple_indicator_cat( locs_rand' ), 1, 3), 2);
    locs_rand(~nonedge_locs) = [];
    pks_rand(~nonedge_locs) = [];
    
    % Store the means and sd of activation strength
    locs_rand_all{jp} = locs_rand;
    pks_rand_all{jp} = pks_rand;
    mus(jp) = mean( activities_rand );
    sds(jp) = std(activities_rand);
    toc
end

%% Identify significant peaks as assembly expression 

Rthresh = nanmean( mus ) + 2*nanmean(sds);
assembly_tag_act = timestamps_cat( locs_tag( pks_tag > Rthresh ) );
ripple_rate_tag = sum( InIntervals( assembly_tag_act, ripple_int ) ) ./ sum( diff(ripple_int') );
ripple_rate_null = nan(nrand,1);
for kp = 1:nrand
    % Detect significant peaks
    assembly_act_rand = timestamps_cat( locs_rand_all{kp}( pks_rand_all{kp} > Rthresh ) );
    ripple_rate_null(kp) = sum( InIntervals( assembly_act_rand, ripple_int ) ) ./ sum( diff(ripple_int') );
end


catch err_msg
    save( outfile_err, 'basepath', 'err_msg' )
end

%%

save( outfile, 'basepath', 'basename', 'ripple_rate_null', 'ripple_rate_tag','pyrInds', 'pyrInds_opto','mus', 'sds', ...
      'pyrInds_nopto', 'ripple_int', 'ripple_int_pad', 'Rthresh', 'nonpulse_period', 'discard', 'beh_ts', 'LIM', 'reclen');

fprintf('Done!\n')


end