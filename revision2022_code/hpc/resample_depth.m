function resample_depth(index)


index = str2double(index);
% Add buzcode to path
addpath(genpath('/gpfs/home/rh2618/buzcode'))
addpath(genpath('/gpfs/scratch/rh2618/code'))
% Path to sessions folder on the HPC
datpath = '/gpfs/scratch/rh2618/rh_data';


% Define
st = dbstack;
funcname = st.name;
outpath = fullfile(datpath, [funcname '_OUT']);
if ~exist(outpath, 'dir')
    mkdir(outpath);
end

outfile = fullfile(outpath,sprintf('resample_%d.mat', index));
if exist(outfile)
    disp('Already fit !')
    return
end

% Load the reference distribution and target distributions
load('refdist_depth.mat')
load('targetdist_depth.mat')   
load('targetdist_var.mat')   

fprintf('%d\n', randval(index))
rng( randval(index) )

pv = [];
while isempty( pv )
    
    % Equalize the depth distributions across birthdates
    indices_e13 = subsample1D(refdist_depth, e13_depth, 'logTransform',false);
    indices_e14 = subsample1D(refdist_depth, e14_depth, 'logTransform',false);
    indices_e15 = subsample1D(refdist_depth, e15_depth, 'logTransform',false);
    indices_e16 = subsample1D(refdist_depth, e16_depth, 'logTransform',false);

    a = [e13_depth(indices_e13) ; e14_depth(indices_e14) ; e15_depth(indices_e15) ; e16_depth(indices_e16)];
    b = [ones( sum( length(e13_depth) ), 1) ; 2.*ones(length(e14_depth), 1) ; 3.*ones(length(e15_depth), 1) ; 4.*ones(length(e16_depth), 1)];
    try
        stats1 = boot_anova1(a,b,'classical', false);
    catch
        continue
    end
    % Check that distribution is indeed equalized
    if stats1.p <= 0.4 || stats1.p >= 0.6
        continue
    end


    a = [e13_var(indices_e13) ; e14_var(indices_e14) ; e15_var(indices_e15) ; e16_var(indices_e16) ];
    b = [ones( length(e13_depth), 1) ; 2.*ones(length(e14_depth), 1) ; 3.*ones(length(e15_depth), 1) ; 4.*ones(length(e16_depth), 1)];
    try
        stats = boot_anova1(a,b,'classical', false,'alpha',0.1);
    catch
        continue
    end
    pv = [pv ; stats.p ];

end

means = [ nanmean(e13_var(indices_e13)) ; nanmean(e14_var(indices_e14)) ; nanmean( e15_var(indices_e15)) ; nanmean( e16_var(indices_e16)) ];

disp('Saving !')
save(outfile, 'pv', 'means')
