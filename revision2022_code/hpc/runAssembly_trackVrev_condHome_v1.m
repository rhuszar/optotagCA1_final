function runAssembly_trackVrev_condHome_v1(index)
% Conditional assembly detection in the homecage ; specifically only in
% ripples - detect assemblies in a spike matrix made of concatenated
% ripples (these are binned) ; then detect smoothed assembly expression
% around ripples

% The hope is that the assemblies we extract will more faithfully reflect
% ripple coactivations

% pyrInds(pyrOnlyInds_opto)    ... UIDs of optotagged PYRs
% pyrInds_opto                 ... UIDs of optotagged PYRs
% pyrOnlyInds_opto             ... IDs of optotagged PYRs in an array of
%                                  only PYRs
% assemblies_cond( 1 ).pyrUID_cond   ... UIDs of PYRs after one has been
%                                        taken out
% assemblies_cond( 1 ).indOpto       ... indices into pyrUID_cond - optotagged units 
% assemblies_cond( 1 ).pyrUID_opto   ... UID of optotagged cell being
%                                        conditioned upon
% Update 3/4/22 - make the code robust to multiple session recordings ;
%                 only considering recordings on the maze

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
load('basepaths_mac_rev1.mat')
%load('basepaths_rerun.mat')
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

outfile = fullfile(outpath,[basename '.mat']);
outfile_err = fullfile(outpath,[basename '_ERROR.mat']);
if exist(outfile)
    disp('Already fit')
    return
end

% If some error occurs, we want to catch it and know about it
try

%% Load data    

load( fullfile(basepath_hpc, [ basename '.spikes.cellinfo.mat']) )         
load( fullfile(basepath_hpc, [ basename '.cell_metrics.cellinfo.mat']) )
load( fullfile(basepath_hpc, [ basename '.pulses.events.mat']) )
load( fullfile(basepath_hpc, [ basename '.ripples.events.mat']) )
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
% discard_period = ConsolidateIntervals( [ pulses.intsPeriods ; beh_ts ; [ discard reclen ] ]);

% Intervals of interest for state-dependent expression strenth
% ( Remove any overlapping stim intervals )
ripple_int = ripples.timestamps;
ripple_int = ripple_int( InIntervals(ripple_int(:,1), nonpulse_period) & InIntervals(ripple_int(:,2), nonpulse_period), : );
ripple_dur = diff( ripple_int' );
ripple_int_pad = [ ripple_int(:,1)-0.2 ripple_int(:,2)+0.05 ];

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);

pyrOnlyInds_opto = arrayfun(@(x) find( x==pyrInds ), pyrInds_opto);


%%

% Smooth / bin data for assembly detection
fprintf('Get smoothed matrices\n')

fprintf('\nGet smoothed matrices in the homecage, outside pulses\n')
spkCounts_coarse = struct('data', [], 'timestamps', [], 'dt', [], 'binsize', []);
spkCounts_smooth = struct('data', [], 'timestamps', [], 'dt', [], 'binsize', []);
ripple_indicator = struct('indicator', []);
for perN = 1:size(ripple_int,1)
    if ripple_dur(perN) < fwhm
        spkCounts_coarse(perN) = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  ripple_dur(perN), 'binsize', ripple_dur(perN), 'win', ripple_int(perN,:));
    else
        spkCounts_coarse(perN) = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  fwhm, 'binsize', fwhm, 'win', ripple_int(perN,:));
    end
    spkCounts_smooth(perN) = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  dt, 'binsize', fwhm, 'win', ripple_int_pad(perN,:),'bintype','gaussian');
    ripple_indicator(perN).indicator = perN.*ones(size( spkCounts_smooth(perN).data,1 ), 1);
end
spkCnts_coarse_cat = cell2mat( {spkCounts_coarse.data}' );
%timestamps_coarse_cat = cell2mat( { spkCounts_coarse.timestamps }' );
spkCnts_smooth_cat = cell2mat( {spkCounts_smooth.data}' );
timestamps_smooth_cat = cell2mat( { spkCounts_smooth.timestamps }' );
ripple_indicator_cat = cell2mat( {ripple_indicator.indicator}' );


assemblies_cond( length(pyrOnlyInds_opto) ) = struct();

%%
% Detect / track assemblies conditioned on track related firing of tagged neurons
index = 1;
indicator_val = 0;
for kp = pyrOnlyInds_opto
    indicator_val = indicator_val+1;
    fprintf('Processing cell %d/%d\n', indicator_val, length(pyrOnlyInds_opto))
    cond_dat = spkCnts_coarse_cat( spkCnts_coarse_cat(:,kp) > 0,: );
    cond_dat(:,kp) = [];
    
    cond_dat1 = spkCnts_smooth_cat; cond_dat1(:,kp) = [];
    
    % We need at least the same number of observations as number of neurons
    % Low rate cells on the track will be discounted..
    [nobs, nsig]  = size(cond_dat);
    if nobs < nsig
        continue
    end
    
    % Keep track of UIDs, which neurons are tagged, which neuron's spikes we conditioned on
    assemblies_cond(index).pyrUID_cond = pyrInds; assemblies_cond(index).pyrUID_cond(kp) = [];
    assemblies_cond(index).indOpto = cell2mat( arrayfun(@(x) find( x == assemblies_cond(index).pyrUID_cond ), pyrInds_opto, 'UniformOutput', false) );
    assemblies_cond(index).pyrUID_opto = pyrInds( kp );
    
    % Extract assemblies
    assemblies = assembly_patterns( cond_dat', opts );
    [~, b] = max(abs(assemblies));          % sign of largest absolute weight connection positive
    lininds = sub2ind( size(assemblies), b, 1:size(assemblies,2));
    sgn = sign(assemblies(lininds));
    assemblies = assemblies .* repmat(sgn, [size(assemblies,1), 1]);
    assemblies = assemblies ./ sqrt(sum(assemblies.^2));
    
    assemblies_cond(index).assemblies = assemblies;
    

    % Expression strength at fine resolution
    activities_smooth = assembly_activity(assemblies, cond_dat1' );
    assemblies_cond(index).activities = assembly_activity(assemblies, cond_dat' );

    % Compute assembly activations on the track and in homecage
    assembly_act = cell( size(activities_smooth,1), 1 );
    pks = cell( size(assemblies,2), 1 );

    Rthresh = nan( size(assemblies,2), 1 );
    % Go over assemblies
    for jp = 1:size(assemblies,2)

        % Adapt threshold between assemblies
        Rthresh(jp) = 2 * std( activities_smooth(jp,:) ) ;
        [pks{jp},locs] = findpeaks( activities_smooth(jp,:) );
        
        % Find peaks within intervals of interest, rather than at the edges between
        nonedge_locs = 3 == sum( ripple_indicator_cat( [locs'-1 locs' locs'+1 ] ) == ...
                                repmat(ripple_indicator_cat( locs' ), 1, 3), 2);
        locs(~nonedge_locs) = [];
        pks{jp}(~nonedge_locs) = [];
        assembly_act{jp} = timestamps_smooth_cat( locs( pks{jp} > Rthresh(jp) ) );

    end
    assemblies_cond(index).assembly_act = assembly_act;
    index = index+1;
end

%%
disp('Attempting to save')
save( outfile, 'basepath', 'assemblies_cond','pyrInds', 'pyrInds_opto', 'pyrOnlyInds_opto', '-v7.3');
catch err_msg

disp('ERROR OCCURRED')
save( outfile, 'basepath', 'err_msg');

end

end
