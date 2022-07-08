function estimate_k_v2(index)
% For a given session, calculate probability of k coactive cells in SPW-R cycles ;
% normalize by the number of bins... 
% Take sessions with at least 2 tagged cells


animal_exclude = [{'e16_3m2'} {'e15_13f1'}];

% Add buzcode to path
addpath(genpath('/gpfs/home/rh2618/buzcode'))
addpath(genpath('/gpfs/scratch/rh2618/code'))
% Path to sessions folder on the HPC
datpath = '/gpfs/scratch/rh2618/rh_data';

nrand = 100;
nit = 500;

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

rip_cycle_ints = cell2mat( ripples.cycle_ints_fixed');
rip_cycle_ints( isnan( rip_cycle_ints(:,1) ), : ) = [];

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
disp('Binning the data')
spkCnts = zeros( length(pyrInds), length(rip_cycle_ints) );
for kp = 1:length(pyrInds)
    [ ~, int_ind, ~ ] = InIntervals( spikes.times{pyrInds(kp)}, rip_cycle_ints );
    int_ind = unique( int_ind ); int_ind(1) = [];
    spkCnts(kp,int_ind) = 1;
end



% For every possible number of coactive neurons, get its probability of occurrence size
k_tag = nan(size(kvals));
for kp = 1:length(kvals)
    k_tag(kp) = sum( sum( spkCnts(pyrInds_opto,:) ) == kvals(kp) ) ./ length(rip_cycle_ints);
end
lambda = mean(spkCnts(pyrInds_opto,:),2);
% Repeat the same in random subsets of non active neurons
k_poiss = nan(nrand, length(kvals));
for jp = 1:nrand
    % Simulate spike counts based on average rates
    spkCnts_null = zeros( length(lambda), length(rip_cycle_ints) );
    for kp = 1:length(lambda)
        spkCnts_null(kp, :) = poissrnd(lambda(kp), 1,  length(rip_cycle_ints));
    end
    spkCnts_null( spkCnts_null > 1 ) = 1;
    for kp = 1:length(kvals)
        k_poiss(jp, kp) = sum( sum( spkCnts_null ) == kvals(kp) ) ./ length(rip_cycle_ints);
    end
end 
% Find nonzero k probabilities for this subset
k_nonzero = find( k_tag > 0 );
prob_zval_tag = [];
for jp = 1:length(k_nonzero)
    % Control distribution - log transform
    d =  k_poiss(:,k_nonzero(jp) ) ;
    mu = mean( d ); sd = std( d );
    prob_zval_tag(jp) = ( k_tag(k_nonzero(jp)) - mu ) ./ sd;
end
prob_zval_tag( isnan( prob_zval_tag ) ) = [];



prob_zval_null = {};
for ip = 1:nit
    fprintf('%d/%d\n', ip, nit)
    randInds = pyrInds_nopto( randi( length(pyrInds_nopto), length(pyrInds_opto), 1 ) );
    k_rand = zeros(1, length(kvals) );
    for kp = 1:length(kvals)
        k_rand(kp) = sum( sum( spkCnts(randInds,:) ) == kvals(kp) ) ./ length(rip_cycle_ints);
    end
    % Get the average rate
    lambda = mean(spkCnts(randInds,:),2);

    % Repeat the same in random subsets of non active neurons
    k_poiss = nan(nrand, length(kvals));
    for jp = 1:nrand
        % Simulate spike counts based on average rates
        spkCnts_null = zeros( length(lambda), length(rip_cycle_ints) );
        for kp = 1:length(lambda)
            spkCnts_null(kp, :) = poissrnd(lambda(kp), 1,  length(rip_cycle_ints));
        end
        spkCnts_null( spkCnts_null > 1 ) = 1;
        for kp = 1:length(kvals)
            k_poiss(jp, kp) = sum( sum( spkCnts_null ) == kvals(kp) ) ./ length(rip_cycle_ints);
        end
    end
    
    % Find nonzero k probabilities for this subset
    k_nonzero = find( k_rand > 0 );
    prob_zval = [];
    for jp = 1:length(k_nonzero)
        % Control distribution - log transform
        d =  k_poiss(:,k_nonzero(jp) ) ;
        mu = mean( d ); sd = std( d );
        prob_zval(jp) = ( k_rand(k_nonzero(jp)) - mu ) ./ sd;
    end
    prob_zval( isnan( prob_zval ) ) = [];
    prob_zval_null{ip} = prob_zval;
end

catch err_msg
    save( outfile_err, 'basepath', 'err_msg' )
end

%%

save( outfile, 'basepath', 'basename', 'pyrInds', 'pyrInds_opto', ...
      'pyrInds_nopto', 'prob_zval_null', 'prob_zval_tag');

fprintf('Done!\n')


end