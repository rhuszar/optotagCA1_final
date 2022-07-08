function estimate_k(index)
% For a given session, calculate rate of k coactive cells in SPW-Rs ;
% normalize by the number of bins... 
% Take sessions with at least 2 tagged cells


animal_exclude = [{'e16_3m2'} {'e15_13f1'}];

% Add buzcode to path
addpath(genpath('/gpfs/home/rh2618/buzcode'))
addpath(genpath('/gpfs/scratch/rh2618/code'))
% Path to sessions folder on the HPC
datpath = '/gpfs/scratch/rh2618/rh_data';

nrand = 1000;

index = str2double(index);
load('basepaths_mac_rev1.mat')
basepath = basepaths_all{ index };

% Basepepath on HPC
path_parts = strsplit(basepath,'/');
basepath_hpc = fullfile(datpath,path_parts{end-2},path_parts{end-1},path_parts{end});
basename = bz_BasenameFromBasepath(basepath_hpc);
disp(basepath_hpc)

% if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), animal_exclude) )
%     disp('This animal is not being considered, returning...')
%     return
% end

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

% rip_ints = cell2mat( ripples.cycle_ints_fixed');
% rip_ints( isnan( rip_ints(:,1) ), : ) = [];
gr = ~InIntervals(ripples.timestamps(:,1), pulses.intsPeriods) & ...
                ~InIntervals(ripples.timestamps(:,2), pulses.intsPeriods);
rip_ints = ripples.timestamps(gr,:);
%%  Obtain the various intervals / indicators into cell types

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);

% Only include sessions with >2 tagged neurons
if length( pyrInds_opto ) < 2
    fprintf('\nInsufficient number of tagged neuorons\n')
    return
end
kvals = 2:length(pyrInds_opto);
    
% Indices into pyrInds which are tagged
pyrInds_opto = arrayfun(@(x) find( x == pyrInds ), pyrInds_opto );
pyrInds_nopto = setdiff( 1:length(pyrInds), pyrInds_opto);

%% Error handling - this code may crash for various reasons (insufficient memory requested, etc.) ; store the error message in this case
try
    
%% Obtain spike matrix binned into ripple cycles

% Obtain the spike matrix - for every ripples cycle, see if the cell spiked
% or not
spkCnts = zeros( length(pyrInds), length(rip_ints) );
for kp = 1:length(pyrInds)
    [ ~, int_ind, ~ ] = InIntervals( spikes.times{pyrInds(kp)}, rip_ints );
    int_ind = unique( int_ind ); int_ind(1) = [];
    spkCnts(kp,int_ind) = 1;

end

% For every possible number of coactive neurons, get its probability of occurrence size
k_tag = nan(size(kvals));
for kp = 1:length(kvals)
    k_tag(kp) = sum( sum( spkCnts(pyrInds_opto,:) ) == kvals(kp) ) ./ length(rip_ints);
end

% Repeat the same in random subsets of non active neurons
k_rand = nan(nrand, length(kvals));
for jp = 1:nrand
    fprintf('%d/%d\n', jp, nrand)
    randInds = pyrInds_nopto( randi( length(pyrInds_nopto), length(pyrInds_opto), 1 ) );
    for kp = 1:length(kvals)
        k_rand(jp, kp) = sum( sum( spkCnts(randInds,:) ) == kvals(kp) ) ./ length(rip_ints);
    end
end


catch err_msg
    save( outfile_err, 'basepath', 'err_msg' )
end

%%

save( outfile, 'basepath', 'basename', 'pyrInds', 'pyrInds_opto', ...
      'pyrInds_nopto', 'k_tag', 'k_rand');

fprintf('Done!\n')


end