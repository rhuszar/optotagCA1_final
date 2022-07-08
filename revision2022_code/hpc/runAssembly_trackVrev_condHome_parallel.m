function runAssembly_trackVrev_condHome_parallel(index, cellID)
% Conditional assembly detection in the homecage

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
%load('basepaths_mac_rev1.mat')
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
% outpath = fullfile(datpath, [funcname '_OUT']);
outpath = fullfile(datpath, 'runAssembly_trackVrev_condHome_OUT', basename);
if ~exist(outpath, 'dir')
    mkdir(outpath);
end


outfile = fullfile(outpath, [basename sprintf('_%d.mat', cellID)] );
outfile_err = fullfile(outpath,[basename sprintf('_%d_ERROR.mat', cellID)]);

% If some error occurs, we want to catch it and know about it
try

%% Load data    

load( fullfile(basepath_hpc, [ basename '.spikes.cellinfo.mat']) )         
load( fullfile(basepath_hpc, [ basename '.cell_metrics.cellinfo.mat']) )
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
% discard_period = ConsolidateIntervals( [ pulses.intsPeriods ; beh_ts ; [ discard reclen ] ]);

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);

fprintf('Index %d/%d of %s', cellID, length(pyrInds_opto), basename)
if cellID > length(pyrInds_opto)
    disp('Already fit everything...')
    return
end
pyrOnlyInds_opto = arrayfun(@(x) find( x==pyrInds ), pyrInds_opto);


%%

% Smooth / bin data for assembly detection
fprintf('Get smoothed matrices\n')

fprintf('\nGet smoothed matrices in the homecage, outside pulses\n')
spkCounts_coarse = struct('data', [], 'timestamps', [], 'dt', [], 'binsize', []);
spkCounts_smooth = struct('data', [], 'timestamps', [], 'dt', [], 'binsize', []);
for perN = 1:size(nonpulse_period,1)
    spkCounts_coarse(perN) = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  fwhm, 'binsize', fwhm, 'win', nonpulse_period(perN,:),'bintype','gaussian');
    spkCounts_smooth(perN) = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  dt, 'binsize', fwhm, 'win', nonpulse_period(perN,:),'bintype','gaussian');
end
spkCnts_coarse_cat = cell2mat( {spkCounts_coarse.data}' );
%timestamps_coarse_cat = cell2mat( { spkCounts_coarse.timestamps }' );
spkCnts_smooth_cat = cell2mat( {spkCounts_smooth.data}' );
timestamps_smooth_cat = cell2mat( { spkCounts_smooth.timestamps }' );


assemblies_cond = struct('pyrUID_cond', [], 'indOpto', [], 'pyrUID_opto', [], 'assemblies', [], 'activities', [], 'assembly_act_track', []);

%%
% Detect / track assemblies conditioned on track related firing of tagged neurons

kp = pyrOnlyInds_opto(cellID);

cond_dat = spkCnts_coarse_cat( spkCnts_coarse_cat(:,kp) > 0,: );
cond_dat(:,kp) = [];

cond_dat1 = spkCnts_smooth_cat; cond_dat1(:,kp) = [];

% We need at least the same number of observations as number of neurons
% Low rate cells on the track will be discounted..
[nobs, nsig]  = size(cond_dat);
if nobs < nsig
    return
end

% Keep track of UIDs, which neurons are tagged, which neuron's spikes we conditioned on
assemblies_cond.pyrUID_cond = pyrInds; assemblies_cond.pyrUID_cond(kp) = [];
assemblies_cond.indOpto = cell2mat( arrayfun(@(x) find( x == assemblies_cond.pyrUID_cond ), pyrInds_opto, 'UniformOutput', false) );
assemblies_cond.pyrUID_opto = pyrInds( kp );

% Extract assemblies
assemblies = assembly_patterns( cond_dat', opts );
[~, b] = max(abs(assemblies));          % sign of largest absolute weight connection positive
lininds = sub2ind( size(assemblies), b, 1:size(assemblies,2));
sgn = sign(assemblies(lininds));
assemblies = assemblies .* repmat(sgn, [size(assemblies,1), 1]);
assemblies = assemblies ./ sqrt(sum(assemblies.^2));

assemblies_cond.assemblies = assemblies;


% Expression strength at fine resolution
activities_track = assembly_activity(assemblies, cond_dat1' );
assemblies_cond.activities = assembly_activity(assemblies, cond_dat' );

% Compute assembly activations on the track and in homecage
assembly_act_track = cell( size(activities_track,1), 1 );
pks_track = cell( size(assemblies,2), 1 );

Rthresh = nan( size(assemblies,2), 1 );
% Go over assemblies
for jp = 1:size(assemblies,2)

    % Adapt threshold between assemblies
    Rthresh(jp) = 2 * std( activities_track(jp,:) ) ;
    [pks_track{jp},locs] = findpeaks( activities_track(jp,:) );       
    assembly_act_track{jp} = timestamps_smooth_cat( locs( pks_track{jp} > Rthresh(jp) ) );

end
assemblies_cond.assembly_act_track = assembly_act_track;

%%

save( outfile, 'basepath', 'assemblies_cond','pyrInds', 'pyrInds_opto', 'pyrOnlyInds_opto','cellID', '-v7.3');


catch err_msg

disp('ERROR OCCURRED')
save( outfile, 'basepath', 'err_msg');

end

end