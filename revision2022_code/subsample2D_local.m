
%{
To resample the firing rates we can just take the 2D histogram from one 
population and use it as a PMF for sampling the 2D histogram of the other population

ref_dist    - 2D reference distribution according to which we sample (tag)
target_dist - target distribution we sample (untagged)
%}


%%




ref_dist = [ theta_tag_rates{2}( ~( isnan( pfcorr_tag_all(:,2) ) ), : ) ; theta_tag_rates{3}( ~( isnan( pfcorr_tag_all(:,3) ) ),: ) ] ;
target_dist = [ theta_ntag_rates{2}( ~( isnan( pfcorr_ntag_all(:,2) ) ), : ) ; theta_ntag_rates{3}( ~( isnan( pfcorr_ntag_all(:,3) ) ),: ) ] ;
outpath = fullfile('/Users/romanhuszar/Documents/rhData/resampleRate/theta_behnov');

% Number of resamplings we do
ns = 500;
%%



for kp = 1:ns

fprintf('%d/%d\n', kp, ns)
outfile = fullfile(outpath, sprintf('subsample_%d.mat', kp));
if exist(outfile)
    continue
end
rng('shuffle')

sz = 100;

target_dist_log = log10( target_dist );
target_dist_log(isinf( target_dist_log(:,1) ), 1) = -2;
target_dist_log(isinf( target_dist_log(:,2) ), 2) = -2;

ref_dist_log = log10( ref_dist );
ref_dist_log(isinf( ref_dist_log(:,1) ), 1) = -2;
ref_dist_log(isinf( ref_dist_log(:,2) ), 2) = -2;

nsamp = 2*length(ref_dist_log);

% Bin the target distribution
[~,Xedges,Yedges,binX_target,binY_target] = histcounts2( target_dist_log(:,1), target_dist_log(:,2),sz, 'Normalization', 'probability' );
% Get bin probabilities in the reference distribution
[p_ref,~,~,binX_ref,binY_ref] = histcounts2( ref_dist_log(:,1), ref_dist_log(:,2), Xedges, Yedges, 'Normalization', 'probability' );


% Transform to linear indices
%tag_lin = sub2ind([sz sz],binX_tag,binY_tag);
ntag_lin = sub2ind([sz sz],binX_target,binY_target);

% Sample
p = p_ref(:);
lin = [ 1:length(p) ]';
[ p_sort, ip ] = sort(p);
lin_sort = lin(ip);

c = cumsum([0,p_sort(:).']);
c = c/c(end);
[~,i] = histc(rand(1,nsamp),c);

% Sampled bins
sampled_bins = lin_sort(i);

indices = [];
for ip = 1:length(sampled_bins)
    v = find( sampled_bins(ip) == ntag_lin );
    if isempty(v)
        continue
    end
    indices = [indices ; v( randi(length(v)) ) ];
end

indices(length(ref_dist_log)+1:end) = [];

save(outfile,'indices')

end

%%

cd('/Users/romanhuszar/Documents/rhData/resampleRate/rip_homecage')
% ntag_resampledVar = ntag_lRat;
% tag_refVar = tagtag_lRat;


% tag_refVar = [ theta_tag_rates{2}( ~( isnan( pfcorr_tag_all(:,2) ) ), : ) ; theta_tag_rates{3}( ~( isnan( pfcorr_tag_all(:,3) ) ),: ) ] ;
% ntag_resampledVar = [ theta_ntag_rates{2}( ~( isnan( pfcorr_ntag_all(:,2) ) ), : ) ; theta_ntag_rates{3}( ~( isnan( pfcorr_ntag_all(:,3) ) ),: ) ] ;
tag_refVar = [ lRat_tag( ~( isnan( pfcorr_tag_all(:,2) ) ), : ) ; lRat_tag( ~( isnan( pfcorr_tag_all(:,3) ) ),: ) ] ;
ntag_resampledVar = [ lRat_ntag( ~( isnan( pfcorr_ntag_all(:,2) ) ), : ) ; lRat_ntag( ~( isnan( pfcorr_ntag_all(:,3) ) ),: ) ] ;
%% Plot the 2D distribution
loglog(ntag_resampledVar(:,1), ntag_resampledVar(:,2), '.r')
hold on
loglog(tag_refVar(:,1), tag_refVar(:,2), '.b')
axis square

%% Plot the 2D distribution

loglog(ntag_resampledVar(indices,1), ntag_resampledVar(indices,2), '.r')
hold on
loglog(tag_refVar(:,1), tag_refVar(:,2), '.b')
axis square




%% Ripple correlation, resampled according to some variable

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( ripcorr_ntag(indices)), [ 1:length(ripcorr_ntag(indices)) ] ./ length(ripcorr_ntag(indices)), '-r');
hold on
plot(sort( ripcorr_tagtag), [ 1:length(ripcorr_tagtag) ] ./ length(ripcorr_tagtag), '-b')
xlim( [-0.075 0.15] )
xlabel('Cofiring (\rho) in SPW-Rs', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)

figure
set(gcf,'Position', [1077 630 560 420])
boxplot([ripcorr_ntag(indices) ; ripcorr_tagtag], [ones(size(ripcorr_ntag(indices))) ; 2*ones(size(ripcorr_tagtag))],'notch','on', 'whisker',inf,...
        'labels',{'different birthdate','same birthdate'})
    ylim([-0.01 0.07])
    
ranksum(ripcorr_ntag(indices), ripcorr_tagtag)

%% Nov env place field corr, resampled according to some variable

tag_dist = [ pfcorr_tag_all( ~( isnan( pfcorr_tag_all(:,2) ) ), 2 ) ; pfcorr_tag_all( ~( isnan( pfcorr_tag_all(:,3) ) ),3 ) ] ;
ntag_dist = [ pfcorr_ntag_all( ~( isnan( pfcorr_ntag_all(:,2) ) ), 2 ) ; pfcorr_ntag_all( ~( isnan( pfcorr_ntag_all(:,3) ) ),3 ) ] ;

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( ntag_dist(indices)), [ 1:length(ntag_dist(indices)) ] ./ length(ntag_dist(indices)), '-r');
hold on
plot(sort( tag_dist), [ 1:length(tag_dist) ] ./ length(tag_dist), '-b')
xlim( [-0.4 0.8] )
xlabel('Spatial ratemap correlation', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)


figure
set(gcf,'Position', [1077 630 560 420])
boxplot([ntag_dist(indices) ; tag_dist], [ones(size(ntag_dist(indices))) ; 2*ones(size(tag_dist))],'notch','on', 'whisker',inf,...
        'labels',{'different birthdate','same birthdate'})
    ylim([-0.2 0.4])

ranksum( ntag_dist(indices), tag_dist)



%% P-value distribution
fils = dir;
fils(1:2) = [];
fils = {fils.name};
pv = [];

% ntag_dist = pfcorr_ntag_all;
% tag_dist = pfcorr_tagtag_all;
tag_dist = [ pfcorr_tag_all( ~( isnan( pfcorr_tag_all(:,2) ) ), 2 ) ; pfcorr_tag_all( ~( isnan( pfcorr_tag_all(:,3) ) ),3 ) ] ;
ntag_dist = [ pfcorr_ntag_all( ~( isnan( pfcorr_ntag_all(:,2) ) ), 2 ) ; pfcorr_ntag_all( ~( isnan( pfcorr_ntag_all(:,3) ) ),3 ) ] ;
for kp = 1:length(fils)
    load( fils{kp} )
    pv(kp) = ranksum(ntag_dist(indices), tag_dist);
end
%%
[~,edges] = histcounts(log10(pv), 20);
histogram(pv,10.^edges,'linestyle', 'none', 'Normalization', 'probability')
set(gca, 'xscale','log')
xline(.05)


%%

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( pfcorr_ntag_all(indices)), [ 1:length(pfcorr_ntag_all(indices)) ] ./ length(pfcorr_ntag_all(indices)), '-r');
hold on
plot(sort( pfcorr_tagtag_all), [ 1:length(pfcorr_tagtag_all) ] ./ length(pfcorr_tagtag_all), '-b')
xlim( [-0.4 0.6] )


