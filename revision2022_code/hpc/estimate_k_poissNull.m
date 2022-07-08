
% For a given session, calculate rate of k coactive cells in SPW-Rs ;
% normalize by the number of bins... 
% Take sessions with at least 2 tagged cells

animal_exclude = [{'e16_3m2'} {'e15_13f1'}];
outpath = 'C:\Users\Roman\Documents\DATA\estimate_k_poissNull_OUT';
nrand = 500;


for bp = 1:length(basepaths_all)
    
    
basepath = alterPath( basepaths_all{bp}, true );
basename = bz_BasenameFromBasepath(basepath);
disp(basename)


if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), animal_exclude) )
    disp('This animal is not being considered, returning...')
    continue
end


outfile = fullfile(outpath,[basename '.mat']);
if exist(outfile)
    disp('Already computed')
    continue
end

% Load data

load( fullfile(basepath, [ basename '.spikes.cellinfo.mat']) )         
load( fullfile(basepath, [ basename '.cell_metrics.cellinfo.mat']) )
load( fullfile(basepath, [ basename '.ripples.events.mat']) )

rip_cycle_ints = cell2mat( ripples.cycle_ints_fixed');
rip_cycle_ints( isnan( rip_cycle_ints(:,1) ), : ) = [];

%  Obtain the various intervals / indicators into cell types

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);

% Only include sessions with >2 tagged neurons
if length( pyrInds_opto ) < 2
    fprintf('\nInsufficient number of tagged neuorons\n')
    continue
end
kvals = 2:length(pyrInds_opto);

% Obtain spike matrix binned into ripple cycles

% Obtain the spike matrix - for every ripples cycle, see if the cell spiked
% or not
spkCnts = zeros( length(pyrInds_opto), length(rip_cycle_ints) );
for kp = 1:length(pyrInds_opto)
    [ ~, int_ind, ~ ] = InIntervals( spikes.times{pyrInds_opto(kp)}, rip_cycle_ints );
    int_ind = unique( int_ind ); int_ind(1) = [];
    spkCnts(kp,int_ind) = 1;

end
lambda = mean(spkCnts,2);

% For every possible number of coactive neurons, get its probability of occurrence size
k_tag = nan(size(kvals));
for kp = 1:length(kvals)
    k_tag(kp) = sum( sum( spkCnts ) == kvals(kp) ) ./ length(rip_cycle_ints);
end

% Repeat the same in random subsets of non active neurons
k_rand = nan(nrand, length(kvals));
for jp = 1:nrand
    fprintf('%d/%d\n', jp, nrand)
    % Simulate spike counts based on average rates
    spkCnts_null = zeros( length(pyrInds_opto), length(rip_cycle_ints) );
    for kp = 1:length(pyrInds_opto)
        spkCnts_null(kp, :) = poissrnd(lambda(kp), 1,  length(rip_cycle_ints));
    end
    spkCnts_null( spkCnts_null > 1 ) = 1;
    for kp = 1:length(kvals)
        k_rand(jp, kp) = sum( sum( spkCnts_null ) == kvals(kp) ) ./ length(rip_cycle_ints);
    end
end

save( outfile, 'basepath', 'basename', 'pyrInds', 'pyrInds_opto', ...
      'pyrInds_nopto', 'k_tag', 'k_rand');

end


%%



fprintf('Done!\n')
