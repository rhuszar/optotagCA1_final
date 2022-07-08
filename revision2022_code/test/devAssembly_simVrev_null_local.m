

% cd('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/DATA/devNull')
% load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/analyze/final_code/revision2022_code/devAssembly_paramsVrev.mat')
% load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/analyze/final_code/revision2022_code/assembly_ratesVrev.mat')
path_to_data = 'C:\Users\Roman\Documents\DATA\devNull';
cd(path_to_data)
load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\final_code\revision2022_code\devAssembly_paramsVrev.mat")
load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\final_code\revision2022_code\assembly_ratesVrev.mat")

nsamples = 10000;
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.method = 'ICA';
opts.Patterns.number_of_iterations = 500;

min_nAssemblies = 10;
nas = [];
for index = 1:420
    mkdir(sprintf('index_%d', index))
for nit = 1:100

outfile = fullfile(path_to_data, sprintf('index_%d', index),sprintf('devSim_%d_num%d.mat', index, nit));
% Time span of correlation
l = l_mat(index);
% Magnitude of excess correlation
varMAX = var_mat(index);


% Birthdate timestamps
dt = .005;
timestamps_birthdate = 12.5:dt:17.5;
peak_neurogenesis = 14.5;
sigma_neurogenesis = 1;
NN = 350;
%%%%%%%%%%%%%%%%%%

% Parameters for controling the shape of the neurogenesis curve
% TODO - see if we can estimate this from data)

% Covariance function
kSQ = @(r) varMAX*exp(-0.5 / (l^2) * r^2 );

assembly_rates_sim = [];
assembly_bdays_sim = [];

% Gamma timescale correlation
tstep = 0.025;
timestamps = [ 0:tstep:(nsamples*tstep-tstep) ]';

rng('shuffle')
flag = false;
slow_index = [];
% Run as many simulations as we need to get ~200 assemblies
while length(assembly_rates_sim) < 200
    
    fprintf('number of assemblies: %d\n', length(assembly_rates_sim))
    % Select the rate of neurogenesis curve
    if strcmp(neurogenesis_shape, 'gaussian')
        n_birthdates = sort( normrnd(peak_neurogenesis, sigma_neurogenesis, NN, 1) );
    elseif strcmp(neurogenesis_shape, 'uniform')
        n_birthdates = sort( min(timestamps_birthdate) + (max(timestamps_birthdate)-min(timestamps_birthdate))*rand(NN,1) );
    else
        error(sprintf('Unrecognized neurogenesis shape ''%s''\n', neurogenesis_shape))
    end

    % Generate the covariance matrix
    Covv = nan(NN,NN);
    if varMAX == 0
        Covv = eye(NN);
    else
        for kp = 1:NN
            for jp = 1:NN
                Covv(kp,jp) = kSQ( n_birthdates(kp)-n_birthdates(jp) );
            end
        end
    end

    Covv = nearestSPD(Covv);
    try
        L = chol(Covv, 'lower');
    catch
        disp('crash')
        continue
    end

    % Samples from standard normals
    s_wn = randn(NN,nsamples);
    % Samples with desired covariance
    s_corr = L * s_wn;

    % Set the target rates
    target_lambda = nan( NN, 1) ;
    target_lambda( n_birthdates < 14 ) = datasample(rates{1}, sum(n_birthdates < 14));
    target_lambda( n_birthdates >= 14 & n_birthdates < 15 ) = datasample(rates{2}, sum( n_birthdates >= 14 & n_birthdates < 15 ));
    target_lambda( n_birthdates >= 15 & n_birthdates < 16 ) = datasample(rates{3}, sum( n_birthdates >= 15 & n_birthdates < 16 ));
    target_lambda( n_birthdates >= 16 ) = datasample(rates{4}, sum( n_birthdates >= 16 ));
    
    u = target_lambda ./ exp( varMAX/2 );
    v = repmat(u,1,nsamples) .* exp(s_corr);

    % Similate data
    sim_spk = poissrnd( v.*tstep );

    assemblyTemplates = assembly_patterns( sim_spk,opts );
    %if size(assemblyTemplates,2) < min_nAssemblies
        flag = true;
        nas = [nas ; size(assemblyTemplates,2)];
        disp('This index is too slow, requires cluster')
        break
    %end
    [~, b] = max(abs(assemblyTemplates));         
    lininds = sub2ind( size(assemblyTemplates), b, 1:size(assemblyTemplates,2));
    sgn = sign(assemblyTemplates(lininds));
    assemblyTemplates = assemblyTemplates .* repmat(sgn, [size(assemblyTemplates,1), 1]);
    assemblyTemplates = assemblyTemplates ./ sqrt(sum(assemblyTemplates.^2));

    activities = assembly_activity(assemblyTemplates, sim_spk );
    assembly_act = cell( size(activities,1), 1 );
    % Potentially could change
    thresh = std( activities(:) );

    for kp = 1:size(assemblyTemplates,2)
        [pks,locs] = findpeaks(activities(kp,:));     % Detect assembly peaks on the track
        assembly_act{kp} = timestamps( locs( pks > 2*thresh) );
    end

    ar_sim = cellfun(@length, assembly_act) ./ timestamps(end);
    assembly_bdays = nan(size(assemblyTemplates,2), 1);
    for kp = 1:size(assemblyTemplates,2)
        assembly_bdays(kp) = mean( n_birthdates( assemblyTemplates(:,kp) > 2*std(assemblyTemplates(:,kp)) ) );
    end

    assembly_rates_sim = [assembly_rates_sim ; ar_sim];
    assembly_bdays_sim = [assembly_bdays_sim ; assembly_bdays];

end

if flag == true
    break
end
% Comparison with data

assembly_rates_simz = zscore( assembly_rates_sim ) + mean(assembly_rates);

meanfunc = [];                    % empty: don't use a mean function
covfunc = @covSEiso;              % Squared Exponental covariance function
likfunc = @likGauss;              % Gaussian likelihood
hyp = struct('mean', [], 'cov', [0 0], 'lik', -1);
hyp2 = minimize(hyp, @gp, -100, @infGaussLik, meanfunc, covfunc, likfunc, assembly_bdays_sim, assembly_rates_simz);
% Obtain the predictive distribution to visualize our posterior
[ymu, ys2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, assembly_bdays_sim, assembly_rates_simz, timestamps_birthdate');

tvec = 13.5:16.5;
assembly_rates_mu_orig = arrayfun(@(x) mean( assembly_rates( birthdates == x ) ), 13:16);
% Obtain the log predictive probabilities of the data - means
[~, ~, ~, ~, lp1] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, assembly_bdays_sim, assembly_rates_simz, tvec', assembly_rates_mu_orig');
% Obtain the log predictive probabilities of the data - all data
[~, ~, ~, ~, lp2] = gp(hyp2, @infGaussLik, meanfunc, covfunc, likfunc, assembly_bdays_sim, assembly_rates_simz, birthdates+.5, assembly_rates);

parameters = struct();
parameters.dt = dt;
parameters.peak_neurogenesis = peak_neurogenesis;
parameters.sigma_neurogenesis = sigma_neurogenesis;
parameters.tvec = tvec;


% f = [ymu+2*sqrt(ys2); flip(ymu-2*sqrt(ys2))];
%   fill([timestamps_birthdate'; flip(timestamps_birthdate')], f, [7 7 7]/8, 'linestyle', 'none')
%   
%   hold on
% plot(assembly_bdays_sim, assembly_rates_simz, '.b')
% plot(timestamps_birthdate, ymu, 'k')

% Save the output
save(outfile ,'lp1', 'lp2', 'assembly_rates_sim','assembly_rates_simz', 'assembly_bdays_sim', 'varMAX', 'l', 'index', 'ymu','ys2', 'timestamps_birthdate', 'parameters')
end
end