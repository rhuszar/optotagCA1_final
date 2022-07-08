function runAssembly_trackVrev_condTrack(index)
% Conditional assembly detection on the track

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

index = str2double(index);
load('basepaths_mac_rev1.mat')
basepath = basepaths_beh{ index };
% basepath = basepaths_all{ index };

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
disp(funcname)
if ~exist(outpath, 'dir')
    mkdir(outpath);
end

outfile = fullfile(outpath,[basename '.mat']);

% If some error occurs, we want to catch it and know about it
try

load( fullfile(basepath_hpc, [ basename '.spikes.cellinfo.mat']) )         
load( fullfile(basepath_hpc, [ basename '.cell_metrics.cellinfo.mat']) )
load( fullfile(basepath_hpc, [basename '.Behavior.mat']) )
load( fullfile(basepath_hpc, [basename '.thetaCycles.events.mat']) )

outfile = fullfile(outpath,[basename '.mat']);
if exist(outfile)
    disp('Already fit')
    return
end

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);

pyrOnlyInds_opto = arrayfun(@(x) find( x==pyrInds ), pyrInds_opto);

% Behavioral timestamps - robust to multiple session recordings
beh_tvec = behavior.timestamps( behavior.masks.recording == 1 );
beh_int = [ beh_tvec(1) beh_tvec(end) ];

theta_cycle_ints = Theta.cycles( InIntervals( Theta.cycles(:,1), beh_int) & InIntervals( Theta.cycles(:,2), beh_int), : );

% Smooth / bin data for assembly detection
fprintf('Get smoothed matrices\n')
spkCounts_smooth = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  dt, 'binsize', fwhm, 'win', beh_int,'bintype','gaussian');

spkCounts_coarse = struct('data', [], 'timestamps', [], 'dt', [], 'binsize', []);
for perN = 1:size(theta_cycle_ints,1)
    spkCounts_coarse(perN) = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  fwhm, 'binsize', fwhm, 'win', theta_cycle_ints(perN,:));
end
spkCounts_coarse_cat = cell2mat( {spkCounts_coarse.data}' );

assemblies_cond( length(pyrOnlyInds_opto) ) = struct();

%%
% Detect / track assemblies conditioned on track related firing of tagged neurons
index = 1;
for kp = pyrOnlyInds_opto
    
    cond_dat = spkCounts_coarse_cat( spkCounts_coarse_cat(:,kp) > 0,: );
    cond_dat(:,kp) = [];
    
    cond_dat1 = spkCounts_smooth.data; cond_dat1(:,kp) = [];
    
    % We need at least the same number of observations as number of neurons
    % Low rate cells on the track will be discounted..
    [nobs, nsig]  = size(cond_dat);
    if nobs < nsig
        continue
    end
    
    % Keep track of UIDs, which neurons are tagged, which neuron's spikes we conditioned on
    assemblies_cond(index).pyrUID_cond = pyrInds; assemblies_cond(index).pyrUID_cond(kp) = [];
    assemblies_cond(index).indOpto = cell2mat( arrayfun(@(x) find( x== assemblies_cond(index).pyrUID_cond ), pyrInds_opto, 'UniformOutput', false) );
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
    activities_track = assembly_activity(assemblies, cond_dat1' );
    assemblies_cond(index).activities = assembly_activity(assemblies, cond_dat' );

    % Compute assembly activations on the track and in homecage
    assembly_act_track = cell( size(activities_track,1), 1 );
    pks_track = cell( size(assemblies,2), 1 );

    Rthresh = nan( size(assemblies,2), 1 );
    % Go over assemblies
    for jp = 1:size(assemblies,2)

        % Adapt threshold between assemblies
        Rthresh(jp) = 2 * std( activities_track(jp,:) ) ;
        [pks_track{jp},locs] = findpeaks( activities_track(jp,:) );       
        assembly_act_track{jp} = spkCounts_smooth.timestamps( locs( pks_track{jp} > Rthresh(jp) ) );

    end
    assemblies_cond(index).assembly_act_track = assembly_act_track;
    index = index+1;
end

%%
save( outfile, 'basepath', 'assemblies_cond','pyrInds', 'pyrInds_opto', 'pyrOnlyInds_opto');
catch err_msg
    
save( outfile, 'basepath', 'err_msg');

end

end






