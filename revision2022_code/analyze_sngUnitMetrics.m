

%% Rerun getHomecageMetrics, if needed for something

for kp = 1:length(basepaths_all)
    basepath = alterPath( basepaths_all{kp}, true );
    fprintf( '%d/%d\n', kp, length(basepaths_all) )
    basename = bz_BasenameFromBasepath(basepath);
    getHomecageMetrics(basepath)
    
end


%%
metrics_cat = struct('fr',[],'burst',[], 'div_index',[], 'ripple_pp', [], 'ripple_fr', [],...
                     'radial_depth', [], 'theta_phase_wake', [], 'theta_mDepth_wake', [], 'theta_sig_wake', [],...
                     'theta_phase_rem', [], 'theta_mDepth_rem', [], 'theta_sig_rem', [], 'isec',[], 'spatialCov', [],'basename', []);
%load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac_rev1.mat")
%load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/basepaths_mac_rev1.mat')

basename_prop = [];
prop_tagged = [];
% Load the homecage metrics
for kp = 1:length(basepaths_all)
   basepath = alterPath(basepaths_all{kp}, true);
   basename = bz_BasenameFromBasepath(basepath);
   fprintf('%d/%d\n', kp, length(basepaths_all))
   v = load( fullfile(basepath, 'HomecageMetrics.mat') );
   metrics_cat.fr = [ metrics_cat.fr ; v.fr ];
   metrics_cat.burst = [ metrics_cat.burst ; v.burst ];
   metrics_cat.div_index = [ metrics_cat.div_index ; v.div_index ];
   metrics_cat.ripple_pp = [ metrics_cat.ripple_pp ; v.ripple_pp ];
   metrics_cat.ripple_fr = [ metrics_cat.ripple_fr ; v.ripple_fr ];
   metrics_cat.radial_depth = [ metrics_cat.radial_depth ; v.radial_depth(:) ];
   metrics_cat.theta_phase_wake = [ metrics_cat.theta_phase_wake ; v.theta_phase_wake ];
   metrics_cat.theta_mDepth_wake = [ metrics_cat.theta_mDepth_wake ; v.theta_mDepth_wake ];
   metrics_cat.theta_sig_wake = [ metrics_cat.theta_sig_wake ; v.theta_sig_wake ];
   metrics_cat.theta_phase_rem = [ metrics_cat.theta_phase_rem ; v.theta_phase_rem ];
   metrics_cat.theta_mDepth_rem = [ metrics_cat.theta_mDepth_rem ; v.theta_mDepth_rem ];
   metrics_cat.theta_sig_rem = [ metrics_cat.theta_sig_rem ; v.theta_sig_rem ];
   metrics_cat.basename = [ metrics_cat.basename ; repmat( {basename}, length(v.fr), 1 ) ];
   if exist( fullfile(basepath, 'BehaviorMetrics.mat')  )
       v1 = load( fullfile(basepath, 'BehaviorMetrics.mat') );
       metrics_cat.isec = [ metrics_cat.isec ; v1.isec ];
       metrics_cat.spatialCov = [ metrics_cat.spatialCov ; v1.spatialCov ];
   else
       metrics_cat.isec = [ metrics_cat.isec ; nan( length(v.fr), 1) ];
       metrics_cat.spatialCov = [ metrics_cat.spatialCov ; nan( length(v.fr), 1) ];
   end
   if length( metrics_cat.isec ) ~= length( metrics_cat.fr )
       error('err')
   end
   basename_prop = [basename_prop ; {basename}];
   prop_tagged = [prop_tagged v.prop_tagged];
end

%%
e13_tag = cellfun(@(x) ~isempty(regexp(x, 'e13', 'once')), basename_prop );
e14_tag = cellfun(@(x) ~isempty(regexp(x, 'e14', 'once')), basename_prop );
e15_tag = cellfun(@(x) ~isempty(regexp(x, 'e15', 'once')), basename_prop ) & ~cellfun(@(x) ~isempty(regexp(x, 'e15_13f1', 'once')), basename_prop );
e16_tag = cellfun(@(x) ~isempty(regexp(x, 'e16', 'once')), basename_prop ) & ~cellfun(@(x) ~isempty(regexp(x, 'e16_3m2', 'once')), basename_prop );

s = 0.1;

set(gcf,'Position', [945   642   366   540])
plot(repmat( 13.5, 1, sum(e13_tag) )+s.*randn(1,sum(e13_tag) ), prop_tagged(e13_tag),'.b' )
hold on
plot(repmat( 14.5, 1, sum(e14_tag) )+s.*randn(1,sum(e14_tag) ), prop_tagged(e14_tag),'.r' )
plot(repmat( 15.5, 1, sum(e15_tag) )+s.*randn(1,sum(e15_tag) ), prop_tagged(e15_tag),'.m' )
plot(repmat( 16.5, 1, sum(e16_tag) )+s.*randn(1,sum(e16_tag) ), prop_tagged(e16_tag),'.y' )

x = 13.5:1:16.5;
y = [nanmean(prop_tagged(e13_tag)) nanmean(prop_tagged(e14_tag)) nanmean(prop_tagged(e15_tag)) nanmean(prop_tagged(e16_tag))];
err = [nanstd(prop_tagged(e13_tag)) nanstd(prop_tagged(e14_tag)) nanstd(prop_tagged(e15_tag)) nanstd(prop_tagged(e16_tag))];
errorbar(x,y,err,'-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red')

% print the means and standard deviations
fprintf('E13.5: mean %.2f, std %.2f\n', nanmean( prop_tagged( e13_tag ) )*100, std( prop_tagged( e13_tag ) )*100 )
fprintf('E14.5: mean %.2f, std %.2f\n', nanmean( prop_tagged( e14_tag ) )*100, std( prop_tagged( e14_tag ) )*100 )
fprintf('E15.5: mean %.2f, std %.2f\n', nanmean( prop_tagged( e15_tag ) )*100, std( prop_tagged( e15_tag ) )*100 )
fprintf('E16.5: mean %.2f, std %.2f\n', nanmean( prop_tagged( e16_tag ) )*100, std( prop_tagged( e16_tag ) )*100 )
xticks([13.5 14.5 15.5 16.5])
ylabel('Proportion of PYRs tagged')
xlabel('Birthdate')
a = [ prop_tagged(e13_tag) prop_tagged(e14_tag) prop_tagged(e15_tag) prop_tagged(e16_tag) ];
b = [ ones( size( prop_tagged(e13_tag) )) 2.*ones( size(prop_tagged(e14_tag) )) 3.*ones( size(prop_tagged(e15_tag) )) 4.*ones( size(prop_tagged(e16_tag) )) ];
stats = boot_anova1(a, b)
%% Label the populations

e13_tag = cellfun(@(x) ~isempty(regexp(x, 'e13', 'once')), metrics_cat.basename );
e14_tag = cellfun(@(x) ~isempty(regexp(x, 'e14', 'once')), metrics_cat.basename );
e15_tag = cellfun(@(x) ~isempty(regexp(x, 'e15', 'once')), metrics_cat.basename ) & ~cellfun(@(x) ~isempty(regexp(x, 'e15_13f1', 'once')), metrics_cat.basename );
e16_tag = cellfun(@(x) ~isempty(regexp(x, 'e16', 'once')), metrics_cat.basename ) & ~cellfun(@(x) ~isempty(regexp(x, 'e16_3m2', 'once')), metrics_cat.basename );

%% Isec

beh_ind = ~isnan(metrics_cat.isec);
a = [ metrics_cat.isec(beh_ind & e13_tag) ; metrics_cat.isec(beh_ind & e14_tag) ; metrics_cat.isec(beh_ind & e15_tag) ; metrics_cat.isec(beh_ind & e16_tag) ];
b = [ ones(sum( beh_ind & e13_tag ), 1) ; 2.*ones(sum( beh_ind & e14_tag ), 1) ; 3.*ones(sum( beh_ind & e15_tag ), 1)  ; 4.*ones(sum( beh_ind & e16_tag ), 1) ];
stats = boot_anova1(a, b, 'classical', false)
stats.multcomp

%% Spatial coverage

beh_ind = ~isnan(metrics_cat.spatialCov);
a = [ metrics_cat.spatialCov(beh_ind & e13_tag) ; metrics_cat.spatialCov(beh_ind & e14_tag) ; metrics_cat.spatialCov(beh_ind & e15_tag) ; metrics_cat.spatialCov(beh_ind & e16_tag) ];
b = [ ones(sum( beh_ind & e13_tag ), 1) ; 2.*ones(sum( beh_ind & e14_tag ), 1) ; 3.*ones(sum( beh_ind & e15_tag ), 1)  ; 4.*ones(sum( beh_ind & e16_tag ), 1) ];
stats = boot_anova1(a, b, 'classical', false)
stats.multcomp

%% PYR->INT divergence
a = [metrics_cat.div_index(e13_tag) ; metrics_cat.div_index(e14_tag) ; metrics_cat.div_index(e15_tag) ; metrics_cat.div_index(e16_tag)];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];
boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Firing rate (Hz)')
set(gca,'yscale', 'log')
% Run the stats
boot_anova1(a,b,'classical', false)
ylabel('Divergence index')
shg

%% Average homecage firing rate

a = [metrics_cat.fr(e13_tag) ; metrics_cat.fr(e14_tag) ; metrics_cat.fr(e15_tag) ; metrics_cat.fr(e16_tag)];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Firing rate (Hz)')
set(gca,'yscale', 'log')
ylim([0.75    2.8])
% Run the stats
stats = boot_anova1(a,b,'classical', false)

%% Homecage burst index

a = [metrics_cat.burst(e13_tag) ; metrics_cat.burst(e14_tag) ; metrics_cat.burst(e15_tag) ; metrics_cat.burst(e16_tag)];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Burst index')
set(gca,'yscale', 'log')
ylim([0.07    0.301])
% Run the stats
stats = boot_anova1(a,b,'classical', false)

%% Homecage ripple firing rate

a = [metrics_cat.ripple_fr(e13_tag) ; metrics_cat.ripple_fr(e14_tag) ; metrics_cat.ripple_fr(e15_tag) ; metrics_cat.ripple_fr(e16_tag)];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Firing rate (Hz) in SPW-Rs')
set(gca,'yscale', 'log')
% Run the stats
stats = boot_anova1(a,b,'classical', false)
ylim([2.6   14])

%% Homecage ripple participation probability

a = [metrics_cat.ripple_pp(e13_tag) ; metrics_cat.ripple_pp(e14_tag) ; metrics_cat.ripple_pp(e15_tag) ; metrics_cat.ripple_pp(e16_tag)];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('SPW-R participation probability')
%set(gca,'yscale', 'log')
% Run the stats
ylim([0   0.4])
stats = boot_anova1(a,b,'classical', false)



%% Split the cells into deep / superficial

to_include = ~cellfun(@(x) ~isempty(regexp(x, 'e15_13f1', 'once')), metrics_cat.basename ) & ~cellfun(@(x) ~isempty(regexp(x, 'e16_3m2', 'once')), metrics_cat.basename );
histogram(metrics_cat.radial_depth( to_include ), 'Normalization', 'probability', 'linestyle', 'none')
xline(nanmedian(metrics_cat.radial_depth(to_include)  ))
xlim([-250   250])
xlabel('Estimated radial depth (um)')

% sup = metrics_cat.radial_depth > 0;
sup = metrics_cat.radial_depth > nanmedian(metrics_cat.radial_depth(to_include) );
deep = ~sup;

%% Firing rate, deep versus superficial

a = [metrics_cat.fr(deep & to_include) ; metrics_cat.fr(sup & to_include) ];
b = [ones( sum( deep & to_include ), 1) ; 2.*ones(sum( sup & to_include ), 1) ];

boxplot(a,b,'labels', {'deep', 'sup'}, 'notch', 'on', 'whisker', 1000)
set(gca,'yscale', 'log')
ylabel('Firing rate (Hz)')
ylim([0.65 2.5])
ranksum(metrics_cat.fr(deep & to_include),  metrics_cat.fr(sup & to_include ))

%metrics_cat.fr(deep & e13_tag)
%%
clc
fprintf( 'Deep: %.3f\n', median( metrics_cat.fr(deep & e13_tag) ) ) 
fprintf( 'Sup: %.3f\n', median( metrics_cat.fr(sup & e13_tag) ) )
%% Burst index, deep versus superficial

a = [metrics_cat.burst(deep & to_include) ; metrics_cat.burst(sup & to_include) ];
b = [ones( sum( deep & to_include ), 1) ; 2.*ones(sum( sup & to_include ), 1) ];

boxplot(a,b,'labels', {'deep', 'sup'}, 'notch', 'on', 'whisker', 1000)
set(gca,'yscale', 'log')
ylim([0.075    0.25])
ylabel('Burst index')
ranksum(metrics_cat.burst(deep & to_include),  metrics_cat.burst(sup & to_include))

%%
clc
fprintf( 'Deep: %.3f\n', median( metrics_cat.burst(deep & e16_tag) ) ) 
fprintf( 'Sup: %.3f\n', median( metrics_cat.burst(sup & e16_tag) ) )
%% Ripple rate, deep versus superficial
a = [metrics_cat.ripple_fr(deep & to_include) ; metrics_cat.ripple_fr(sup & to_include) ];
b = [ones( sum( deep & to_include ), 1) ; 2.*ones(sum( sup & to_include ), 1) ];

boxplot(a,b,'labels', {'deep', 'sup'}, 'notch', 'on', 'whisker', 1000)
ylabel('Firing rate (Hz) in SPW-Rs')
set(gca,'yscale', 'log')
ylim([3.5   16.5])
ranksum(metrics_cat.ripple_fr(deep & to_include),metrics_cat.ripple_fr(sup & to_include))
%%
clc
fprintf( 'Deep: %.3f\n', median( metrics_cat.ripple_fr(deep & e16_tag) ) ) 
fprintf( 'Sup: %.3f\n', median( metrics_cat.ripple_fr(sup & e16_tag) ) )

%% Ripple participation probability, deep versus superficial
a = [metrics_cat.ripple_pp(deep & to_include) ; metrics_cat.ripple_pp(sup & to_include) ];
b = [ones( sum( deep & to_include ), 1) ; 2.*ones(sum( sup & to_include ), 1) ];

boxplot(a,b,'labels', {'deep', 'sup'}, 'notch', 'on', 'whisker', 1000)
ylim([ 0.1    0.4])
ylabel('SPW-R participation probability')
ranksum(metrics_cat.ripple_pp(deep & to_include) ,metrics_cat.ripple_pp(sup & to_include) )
%%
clc
fprintf( 'Deep: %.3f\n', median( metrics_cat.ripple_pp(deep & e16_tag) ) ) 
fprintf( 'Sup: %.3f\n', median( metrics_cat.ripple_pp(sup & e16_tag) ) )

%%

clc
fprintf( 'Deep: %.3f\n', nanmedian( metrics_cat.spatialCov(deep & e16_tag) ) ) 
fprintf( 'Sup: %.3f\n', nanmedian( metrics_cat.spatialCov(sup & e16_tag) ) )

%% Estimated depths of cells at various birthdates

a = [metrics_cat.radial_depth(e13_tag) ; metrics_cat.radial_depth(e14_tag) ; metrics_cat.radial_depth(e15_tag) ; metrics_cat.radial_depth(e16_tag)];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Estimated radial depth (um)')
ylim([-60 20])
set(gca, 'YDir','reverse')

% set(gca,'yscale', 'log')
% Run the stats
stats = boot_anova1(a,b,'classical', false)

%%

alpha = 0.05;
% indicator of significant phase locking in WAKE / REM
siglocked_indicator = [ metrics_cat.theta_sig_wake<alpha   metrics_cat.theta_sig_rem<alpha ];

% theta phase preference in degrees - only keep significantly locked
thetaphase_pref = [ rad2deg( metrics_cat.theta_phase_wake ) rad2deg( metrics_cat.theta_phase_rem )];

% Cells locked to trough in wake
wake_trough = thetaphase_pref(:,1) >= 120 & thetaphase_pref(:,1) <= 300;
% Cells locked to peak in REM
rem_peak = ( thetaphase_pref(:,2) >= 0 & thetaphase_pref(:,2) < 120 ) |  ( thetaphase_pref(:,2) > 300 & thetaphase_pref(:,2) < 360 );

a = sum( wake_trough( e16_tag & deep  & all( siglocked_indicator, 2)) & rem_peak( e16_tag & deep  & all( siglocked_indicator, 2) ) ) ./ sum(e16_tag & deep & all(siglocked_indicator,2)   )
b = sum( wake_trough( e16_tag & sup & all( siglocked_indicator, 2)) & rem_peak( e16_tag & sup & all( siglocked_indicator, 2) ) ) ./ sum(e16_tag & sup & all(siglocked_indicator,2)  )

%% Save the reference depths - NOTE, this is for homecage data
cd('C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\final_code\revision2022_code\hpc')
e13_depth = metrics_cat.radial_depth(e13_tag & ~isnan( metrics_cat.spatialCov ));
e14_depth = metrics_cat.radial_depth(e14_tag & ~isnan( metrics_cat.spatialCov ));
e15_depth = metrics_cat.radial_depth(e15_tag & ~isnan( metrics_cat.spatialCov ));
e16_depth = metrics_cat.radial_depth(e16_tag & ~isnan( metrics_cat.spatialCov ));

refdist_depth = metrics_cat.radial_depth( to_include );

rng('shuffle')
randval = randi(20000, 1,20000); randval = unique(randval);
randval = randval( randperm(length(randval)) );
randval = randval(1:500);

save('refdist_depth.mat', 'refdist_depth', 'randval')
save('targetdist_depth.mat', 'e13_depth', 'e14_depth', 'e15_depth', 'e16_depth') 

%% Prepare variable for the hpc

e13_var = metrics_cat.spatialCov(e13_tag & ~isnan( metrics_cat.spatialCov ) );
e14_var = metrics_cat.spatialCov(e14_tag & ~isnan( metrics_cat.spatialCov ) );
e15_var = metrics_cat.spatialCov(e15_tag & ~isnan( metrics_cat.spatialCov ));
e16_var = metrics_cat.spatialCov(e16_tag & ~isnan( metrics_cat.spatialCov ));
save('targetdist_var.mat', 'e13_var', 'e14_var', 'e15_var', 'e16_var') 

%%
close all
cd('C:\Users\Roman\Documents\DATA\resample_depth_OUT')
fils = dir; fils = {fils.name};
fils(1:2) = [];
pv_all = [];
means_all = [];
h = {};

set(gcf,'Position', [863   416   391   420])
hold on
for kp = 1:length(fils)
    fprintf('%d/%d\n', kp, length(fils))
    load( fils{kp} )
    pv_all(kp) = pv;
    means_all = [means_all means];
    h{kp} = plot(means,'k');
    h{kp}.Color(4) = 0.03;
end
plot(mean( means_all' ), 'b', 'linewidth', 3)

figure
%pv_all( pv_all == 0 ) = 0.00000000000000000001;
[~,edges] = histcounts(log10(pv_all));
 histogram(pv_all,10.^edges, 'linestyle', 'none',  'normalization', 'probability')
 set(gca, 'xscale','log')
 xline(0.05)
 
 %%
cd('C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\figs\revision2022_figs\Supp_depth\data')

save('spatialCoverageResample.mat','means_all', 'pv_all')
%% Resample according to the generate distribution
e13_depth = metrics_cat.radial_depth(e13_tag);
e14_depth = metrics_cat.radial_depth(e14_tag);
e15_depth = metrics_cat.radial_depth(e15_tag);
e16_depth = metrics_cat.radial_depth(e16_tag);
indices_e13 = subsample1D(metrics_cat.radial_depth( to_include ), e13_depth, 'logTransform',false);
indices_e14 = subsample1D(metrics_cat.radial_depth( to_include ), e14_depth, 'logTransform',false);
indices_e15 = subsample1D(metrics_cat.radial_depth( to_include ), e15_depth, 'logTransform',false);
indices_e16 = subsample1D(metrics_cat.radial_depth( to_include ), e16_depth, 'logTransform',false);

a = [e13_depth(indices_e13) ; e14_depth(indices_e14) ; e15_depth(indices_e15) ; e16_depth(indices_e16)];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Estimated radial depth (um)')
ylim([-60 20])
set(gca, 'YDir','reverse')
shg

stats = boot_anova1(a,b,'classical', false,'alpha',0.1)

%% Resampled bursting
e13_burst = metrics_cat.burst(e13_tag);
e14_burst = metrics_cat.burst(e14_tag);
e15_burst = metrics_cat.burst(e15_tag);
e16_burst = metrics_cat.burst(e16_tag);
a = [e13_burst(indices_e13) ; e14_burst(indices_e14) ; e15_burst(indices_e15) ; e16_burst(indices_e16) ];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Burst index')
set(gca,'yscale', 'log')
ylim([0.07    0.35])

stats = boot_anova1(a,b,'classical', false,'alpha',0.1)

%% Resampled sharp wave ripple participation

e13_burst = metrics_cat.ripple_pp(e13_tag);
e14_burst = metrics_cat.ripple_pp(e14_tag);
e15_burst = metrics_cat.ripple_pp(e15_tag);
e16_burst = metrics_cat.ripple_pp(e16_tag);
a = [e13_burst(indices_e13) ; e14_burst(indices_e14) ; e15_burst(indices_e15) ; e16_burst(indices_e16) ];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Burst index')
%set(gca,'yscale', 'log')
ylim([0.07    0.5])

stats = boot_anova1(a,b,'classical', false,'alpha',0.1)


%%
% Resample
pv_all = [];
while length(pv_all) < 500

    fprintf('%d\n', length(pv_all))
indices_e13 = subsample1D(metrics_cat.radial_depth( to_include ), e13_depth, 'logTransform',false);
indices_e14 = subsample1D(metrics_cat.radial_depth( to_include ), e14_depth, 'logTransform',false);
indices_e15 = subsample1D(metrics_cat.radial_depth( to_include ), e15_depth, 'logTransform',false);
indices_e16 = subsample1D(metrics_cat.radial_depth( to_include ), e16_depth, 'logTransform',false);

a = [e13_depth(indices_e13) ; e14_depth(indices_e14) ; e15_depth(indices_e15) ; e16_depth(indices_e16)];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];
try
stats1 = boot_anova1(a,b,'classical', false);
catch
    continue
end

if stats1.p <= 0.3
    continue
end

e13_burst = metrics_cat.burst(e13_tag);
e14_burst = metrics_cat.burst(e14_tag);
e15_burst = metrics_cat.burst(e15_tag);
e16_burst = metrics_cat.burst(e16_tag);

a = [e13_burst(indices_e13) ; e14_burst(indices_e14) ; e15_burst(indices_e15) ; e16_burst(indices_e16) ];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];
try
stats = boot_anova1(a,b,'classical', false,'alpha',0.1);
catch
    continue
end
pv_all = [pv_all ; stats.p ];

end
    

%%
alpha = 0.01;
depthMod_wake = metrics_cat.theta_mDepth_wake  ( metrics_cat.theta_sig_wake < alpha );
depthMod_rem = metrics_cat.theta_mDepth_rem ( metrics_cat.theta_sig_rem < alpha );


%% Theta depth modulation - mean resultant vector
%  Across brain states

alpha = 0.01;
e13_depthMod_wake = metrics_cat.theta_mDepth_wake ( e13_tag & metrics_cat.theta_sig_wake < alpha );
e14_depthMod_wake = metrics_cat.theta_mDepth_wake  ( e14_tag & metrics_cat.theta_sig_wake < alpha );
e15_depthMod_wake = metrics_cat.theta_mDepth_wake  ( e15_tag & metrics_cat.theta_sig_wake < alpha );
e16_depthMod_wake = metrics_cat.theta_mDepth_wake  ( e16_tag & metrics_cat.theta_sig_wake < alpha );


e13_depthMod_rem = metrics_cat.theta_mDepth_rem ( e13_tag & metrics_cat.theta_sig_rem < alpha );
e14_depthMod_rem = metrics_cat.theta_mDepth_rem ( e14_tag & metrics_cat.theta_sig_rem < alpha );
e15_depthMod_rem = metrics_cat.theta_mDepth_rem ( e15_tag & metrics_cat.theta_sig_rem < alpha );
e16_depthMod_rem = metrics_cat.theta_mDepth_rem ( e16_tag & metrics_cat.theta_sig_rem < alpha );


figure
set(gcf,'Position', [1000         518         352         820])
subplot(2,1,1)
a = [e13_depthMod_wake ; e14_depthMod_wake ; e15_depthMod_wake ; e16_depthMod_wake];
b = [ones(size(e13_depthMod_wake)) ; 2*ones(size(e14_depthMod_wake)) ; 3*ones(size(e15_depthMod_wake)) ; 4*ones(size(e16_depthMod_wake))];
boxplot(a, b,'notch','on', 'whisker',inf)
ylim([0.05 0.17])
xticks(1:4)
xticklabels({'E13.5', 'E14.5', 'E15.5', 'E16.6'})
ylabel('Theta depth modulation', 'fontsize', 18)
xlabel('Birthdate', 'fontsize', 18)
title('WAKE','fontsize', 18)
stats1 = boot_anova1(a,b, 'classical', false);

c = [e13_depthMod_rem ; e14_depthMod_rem ; e15_depthMod_rem ; e16_depthMod_rem];
d = [ones(size(e13_depthMod_rem)) ; 2*ones(size(e14_depthMod_rem)) ; 3*ones(size(e15_depthMod_rem)) ; 4*ones(size(e16_depthMod_rem))];
subplot(2,1,2)
boxplot( c, d, 'notch','on', 'whisker',inf)
xticks(1:4)
ylim([0.05 0.25])
xticklabels({'E13.5', 'E14.5', 'E15.5', 'E16.6'})
ylabel('Theta depth modulation', 'fontsize', 18)
xlabel('Birthdate', 'fontsize', 18)
title('REM','fontsize', 18)

stats2 = boot_anova1(c,d,'classical', false);


%% Theta rem shift fraction

alpha = 0.05;
% indicator of significant phase locking in WAKE / REM
siglocked_indicator = [ metrics_cat.theta_sig_wake<alpha   metrics_cat.theta_sig_rem<alpha ];
siglocked_indicator( isnan( metrics_cat.theta_sig_rem ) | isnan( metrics_cat.theta_sig_wake ), : ) = [];

% theta phase preference in degrees - only keep significantly locked
thetaphase_pref = [ rad2deg( metrics_cat.theta_phase_wake ) rad2deg( metrics_cat.theta_phase_rem )];
thetaphase_pref( isnan( metrics_cat.theta_sig_rem ) | isnan( metrics_cat.theta_sig_wake ), : ) = [];
basenames = metrics_cat.basename; basenames( isnan( metrics_cat.theta_sig_rem ) | isnan( metrics_cat.theta_sig_wake ) ) = [];

% Extract significantly phase locked cells in both brain states
thetaphase_pref_siglock = thetaphase_pref( all(siglocked_indicator,2),: );
% Cells locked to trough in wake
wake_trough = thetaphase_pref_siglock(:,1) >= 120 & thetaphase_pref_siglock(:,1) <= 300;
% Cells locked to peak in REM
rem_peak = ( thetaphase_pref_siglock(:,2) >= 0 & thetaphase_pref_siglock(:,2) < 120 ) |  ( thetaphase_pref_siglock(:,2) > 300 & thetaphase_pref_siglock(:,2) < 360 );

% keep again only relevant basenames
basenames = basenames( all(siglocked_indicator,2) );
% Get indicators for different birthdates
e13_tag = cellfun(@(x) ~isempty(regexp(x, 'e13', 'once')), basenames );% & ~cellfun(@(x) ~isempty(regexp(x, 'e13_26m1', 'once')), basenames );
e14_tag = cellfun(@(x) ~isempty(regexp(x, 'e14', 'once')), basenames );
e15_tag = cellfun(@(x) ~isempty(regexp(x, 'e15', 'once')), basenames ) & ~cellfun(@(x) ~isempty(regexp(x, 'e15_13f1', 'once')), basenames );
e16_tag = cellfun(@(x) ~isempty(regexp(x, 'e16', 'once')), basenames ) & ~cellfun(@(x) ~isempty(regexp(x, 'e16_3m2', 'once')), basenames );

remshift_e13 = sum( wake_trough(e13_tag) & rem_peak(e13_tag) ) ./ sum(e13_tag);
remshift_e14 = sum( wake_trough(e14_tag) & rem_peak(e14_tag) ) ./ sum(e14_tag);
remshift_e15 = sum( wake_trough(e15_tag) & rem_peak(e15_tag) ) ./ sum(e15_tag);
remshift_e16 = sum( wake_trough(e16_tag) & rem_peak(e16_tag) ) ./ sum(e16_tag);

%% Reorganize how we save the cells
gp = [ e13_tag + 2.*e14_tag + 3.*e15_tag + 4.*e16_tag ];
y = wake_trough & rem_peak;
y(gp == 0) = [];
gp(gp == 0) = [];
% load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/analyze/figs/shortFormat/Figure1/remshift_indicator.mat')

%%
nresample = 1000;
% Resample the proportions in a loop
% For each resampling, store the fraction, compute the leas
yres = nan(nresample, 4);
slp_res = nan(nresample, 1);
frac = @(x) sum(x) / length(x);

% Bootstrapped estimate of the 
for kp = 1:nresample
    
    ysur = [];
    for jp = 1:4
        gpt = find( gp==jp );
        ysur = [ ysur {y( gpt( randi( length( gpt ), length( gpt ), 1 ) ) ) }];
    end
    coeff = polyfit([1:4], cellfun(frac, ysur ), 1);
    slp_res(kp) = coeff(1);
    yres(kp,:) = cellfun(frac, ysur );
    
end


data = [mean( y(gp==1) ) mean( y(gp==2) ) mean( y(gp==3) ) mean( y(gp==4) )];

coeff = polyfit([1:4], data, 1);
%%

reg = polyval(coeff, 1:4);
upper = polyval([prctile(slp_res, 97.5) coeff(2)], 1:4);
lower = polyval([prctile(slp_res, 2.5) coeff(2)], 1:4);


f = [upper flip(lower)];
figure
set(gcf,'Position', [1038         742         465         291])
fill([1:4 flip(1:4)], f, 'b', 'FaceAlpha', 0.2, 'linestyle', 'none' )
hold on
plot(1:4, reg, '--b', 'Linewidth', 2)

xticks([1 2 3 4])
xticklabels({'E13.5', 'E14.5', 'E15.5', 'E16.5'})
yticks([0 0.1 0.2 0.3 0.4])
yticklabels({'0', '10', '20', '30', '40'})
ylabel('Percentage of REM-shifters')
% Plot bootstrapped error bar
ci_pos = arrayfun(@(x) prctile( yres(:,x), 97.5 ), 1:4 );
ci_neg = arrayfun(@(x) prctile( yres(:,x), 2.5 ), 1:4 );
errorbar(1:4,data,abs( ci_neg-data ),abs( ci_pos-data ),'k+')
xlim([0.8 4.2])

%% Resample the proportions

% Total number of samples
ns = length(y);
% Group sizes
gp_ids_uniq = unique( gp ) ;
gp_sizes = arrayfun(@(x) sum( x == gp ), gp_ids_uniq);
t_surr = nan(1, nresample);
for kp = 1:nresample
    y_boot =  arrayfun(@(x) y(randi(ns, 1, gp_sizes(x))), 1:length(gp_sizes), 'UniformOutput', false);
    lm_boot = fitlm( 1:4, cellfun(@mean, y_boot ) );
    t_surr(kp) = lm_boot.Coefficients{'x1', 'tStat'};
end

%%

lm = fitlm( 1:4, data );
histogram(t_surr); 
xline( lm.Coefficients{'x1', 'tStat'} )
sum( t_surr < lm.Coefficients{'x1', 'tStat'} ) / 1000 
