
% Here we do whatever constitutes a pairwise analysis

%% Get cofiring in the homecage 

%load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\final_code\revision2022_code\corrHome_v7.mat")
load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/analyze/final_code/revision2022_code/corrHome_v7.mat')

%% Ripple cofiring (fig 3)

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( ripcorr_ntag), [ 1:length(ripcorr_ntag) ] ./ length(ripcorr_ntag), '-r');
hold on
plot(sort( ripcorr_tagtag), [ 1:length(ripcorr_tagtag) ] ./ length(ripcorr_tagtag), '-b')
xlim( [-0.075 0.15] )
xlabel('Cofiring (\rho) in SPW-Rs', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)

figure
set(gcf,'Position', [1077 630 560 420])
boxplot([ripcorr_ntag ; ripcorr_tagtag], [ones(size(ripcorr_ntag)) ; 2*ones(size(ripcorr_tagtag))],'notch','on', 'whisker',inf,...
        'labels',{'different birthdate','same birthdate'})
    ylim([-0.01 0.07])
    
    
%% Ripple cofiring- individual animals

uniq_animals = unique(animal_tag);
figure(1)
set(gcf,'Position', [ 161          87         810        1219 ])
figure(2)
set(gcf,'Position', [1598         106         810        1219])
for kp = 1:length( uniq_animals )
    

    tt_fam = ripcorr_tagtag(cellfun(@(x) strcmp(x, uniq_animals{kp}), animal_tag),1);
    nt_fam = ripcorr_ntag(cellfun(@(x) strcmp(x, uniq_animals{kp}), animal_ntag),1);
    
    figure(1)
    subplot(4,4, kp)
    plot(sort(tt_fam), [ 1:length(tt_fam) ] ./length(tt_fam), 'b' )
    hold on
    plot(sort(nt_fam), [ 1:length(nt_fam) ] ./length(nt_fam), 'r' )
    xlim( [-0.075 0.15] )
    title(sprintf('%s', uniq_animals{kp}),  'Interpreter', 'none')
    
    figure(2)
    subplot(4,4, kp)
    boxplot([nt_fam ; tt_fam], [ones(size(nt_fam)) ; 2.*ones(size(tt_fam))], 'notch', 'on', 'whisker', inf)
    ylim([-0.02 0.09])
    title(sprintf('%e', ranksum(tt_fam, nt_fam) ))
    
end
    
    
%% Convergence onto interneurons, SBD v. DBD (fig 2B)

figure
set(gcf,'Position', [1077 630 560 420])
boxplot([converge_ntag ; converge_tagtag], [ones(size(converge_ntag)) ; 2*ones(size(converge_tagtag))],'notch','on', 'whisker',inf,...
        'labels',{'different birthdate','same birthdate'})
ylim([-0.025    0.55])
ylabel('Convergence onto INTs')
    
%% Convergence onto interneurons, split by birthdate (fig 2C)

e13_tag = cellfun(@(x) ~isempty( regexp( x, 'e13', 'once' )), basenames_tagtag );
e14_tag = cellfun(@(x) ~isempty( regexp( x, 'e14', 'once' )), basenames_tagtag );
e15_tag = cellfun(@(x) ~isempty( regexp( x, 'e15', 'once' )), basenames_tagtag );
e16_tag = cellfun(@(x) ~isempty( regexp( x, 'e16', 'once' )), basenames_tagtag );

e13_ntag = cellfun(@(x) ~isempty( regexp( x, 'e13', 'once' )), basenames_ntag );
e14_ntag = cellfun(@(x) ~isempty( regexp( x, 'e14', 'once' )), basenames_ntag );
e15_ntag = cellfun(@(x) ~isempty( regexp( x, 'e15', 'once' )), basenames_ntag );
e16_ntag = cellfun(@(x) ~isempty( regexp( x, 'e16', 'once' )), basenames_ntag );

% Correct convergences by median of DBD population
a = [converge_tagtag(e13_tag) - nanmedian( converge_ntag(e13_ntag) ); converge_tagtag(e14_tag) - nanmedian( converge_ntag(e14_ntag) ) ; ...
    converge_tagtag(e15_tag) - nanmedian( converge_ntag(e15_ntag) ) ; converge_tagtag(e16_tag) - nanmedian( converge_ntag(e16_ntag) )];
b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylim([-0.2 0.4])
ylabel('Convergence onto INTs (corrected)')

% Express convergence as a z-score with respect to the DBD convergence distribution at each birthdate population
% a = [ ( log10( converge_tagtag(e13_tag) ) - nanmean( log10( converge_ntag(e13_ntag) ) ) ) ./ nanstd( log10( converge_ntag(e13_ntag) ) ); ...
%       ( log10( converge_tagtag(e14_tag) ) - nanmean( log10( converge_ntag(e14_tag) ) ) ) ./ nanstd( log10( converge_ntag(e14_tag) ) ) ; ...
%       ( log10( converge_tagtag(e15_tag) ) - nanmean( log10( converge_ntag(e15_tag) ) ) ) ./ nanstd( log10( converge_ntag(e15_tag) ) ) ; ...
%       ( log10( converge_tagtag(e16_tag) ) - nanmean( log10( converge_ntag(e16_tag) ) ) ) ./ nanstd( log10( converge_ntag(e16_tag) ) ) ];
% b = [ones( sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1)];

%boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)

stats = boot_anova1(a,b, 'classical', false);



%% Theta cycle cofiring (fig 6C)

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( thetacorr_ntag), [ 1:length(thetacorr_ntag) ] ./ length(thetacorr_ntag), '-r');
hold on
plot(sort( thetacorr_tagtag), [ 1:length(thetacorr_tagtag) ] ./ length(thetacorr_tagtag), '-b')
xlim( [-0.075 0.15] )
xlabel('Theta cycle cofiring (\rho)', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)

figure
set(gcf,'Position', [1077 630 560 420])
boxplot([thetacorr_ntag ; thetacorr_tagtag], [ones(size(thetacorr_ntag)) ; 2*ones(size(thetacorr_tagtag))],'notch','on', 'whisker',inf,...
        'labels',{'different birthdate','same birthdate'})
    ylim([-0.015 0.055])
    

%% Theta cycle cofiring - individual animals

uniq_animals = unique(animal_tag);
figure(1)
set(gcf,'Position', [ 161          87         810        1219 ])
figure(2)
set(gcf,'Position', [1598         106         810        1219])
for kp = 1:length( uniq_animals )
    

    tt_fam = thetacorr_tagtag(cellfun(@(x) strcmp(x, uniq_animals{kp}), animal_tag),1);
    nt_fam = thetacorr_ntag(cellfun(@(x) strcmp(x, uniq_animals{kp}), animal_ntag),1);
    
    figure(1)
    subplot(4,4, kp)
    plot(sort(tt_fam), [ 1:length(tt_fam) ] ./length(tt_fam), 'b' )
    hold on
    plot(sort(nt_fam), [ 1:length(nt_fam) ] ./length(nt_fam), 'r' )
    xlim( [-0.075 0.15] )
    title(sprintf('%s', uniq_animals{kp}),  'Interpreter', 'none')
    
    figure(2)
    subplot(4,4, kp)
    boxplot([nt_fam ; tt_fam], [ones(size(nt_fam)) ; 2.*ones(size(tt_fam))], 'notch', 'on', 'whisker', inf)
    ylim([-0.02 0.09])
    title(sprintf('%e', ranksum(tt_fam, nt_fam) ))
    
end
    
%% PYR->INT convergence without labels
% Show the distribution, excluding cell pairs that don't converge on any
% interneuron
% Divide distribution into percentiles

datatype = 'place'; % 'theta', 'place'

if strcmp(datatype, 'ripple') || strcmp(datatype, 'theta')
    load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\final_code\revision2022_code\corrHome_v7.mat")
elseif strcmp(datatype, 'place')
    load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\final_code\revision2022_code\corrBeh_v7.mat")
end

npercentiles = 5;

converge_nolab_nonzero_ind = converge_nolab > 0;
converge_nolab_nonzero = converge_nolab( converge_nolab_nonzero_ind );

perc_bins = prctile( converge_nolab_nonzero, 100.* [ 0:1 ./ npercentiles:1 ] );


pd = fitdist(converge_nolab_nonzero,'kernel','Width',0.05);
x = -0.05:0.025:1.1;
y = pdf(pd,x);

figure
set(gcf,'Position', [1000         918         560         281])
histogram(converge_nolab_nonzero,20, 'Normalization', 'pdf', 'DisplayStyle', 'stairs' );
hold on

plot(x,y  ,'Color','k','LineStyle','--')
xlim([min(x) max(x)])

ylim_vals = ylim;
xlim_vals = xlim;


boundaries = perc_bins(2:end-1); % we use edges defined by the plot for visualization
b1 = rectangle('Position', [ xlim_vals(1), ylim_vals(1), boundaries(1) - xlim_vals(1), ylim_vals(2) ], 'FaceColor',[0.9290 0.6940 0.1250], 'linestyle', 'none');
b1.FaceColor(4) = 0.2;
b2 = rectangle('Position', [ boundaries(1), ylim_vals(1), boundaries(2) - boundaries(1), ylim_vals(2) ], 'FaceColor',[0.9290 0.6940 0.1250], 'linestyle', 'none');
b2.FaceColor(4) = 0.3;
b3 = rectangle('Position', [ boundaries(2), ylim_vals(1), boundaries(3) - boundaries(2), ylim_vals(2) ], 'FaceColor',[0.9290 0.6940 0.1250], 'linestyle', 'none');
b3.FaceColor(4) = 0.4;
b4 = rectangle('Position', [ boundaries(3), ylim_vals(1), boundaries(4) - boundaries(3), ylim_vals(2) ], 'FaceColor',[0.9290 0.6940 0.1250], 'linestyle', 'none');
b4.FaceColor(4) = 0.5;
b5 = rectangle('Position', [ boundaries(4), ylim_vals(1), xlim_vals(2) - boundaries(4), ylim_vals(2) ], 'FaceColor',[0.9290 0.6940 0.1250], 'linestyle', 'none');
b5.FaceColor(4) = 0.6;
% Now divide this distribution into percentiles
xlabel('Convergence onto INTs', 'fontsize', 14)
ylabel('Probability density', 'fontsize', 14)

figure
bar( 100.*[ mean(converge_nolab == 0) mean(converge_nolab > 0)] )
ylabel('Percentage of PYR pairs')

%% We could do the same statistical analysis we did for the depth in the pyramidal layer


% Choose variable of interest - ripple, theta, or place related
% correlations
datatype = 'place'; % 'theta', 'place'

if strcmp(datatype, 'ripple')
    y = ripcorr_nolab(converge_nolab_nonzero_ind);
    ytag = ripcorr_tagtag;
elseif strcmp(datatype, 'theta')
    y = thetacorr_nolab(converge_nolab_nonzero_ind);
    ytag = thetacorr_tagtag;
elseif strcmp(datatype, 'place')
    y = pfcorr_nolab_all(converge_nolab_nonzero_ind);
    ytag = pfcorr_tagtag_all;
end
[~, ~, gp] = histcounts(converge_nolab_nonzero, perc_bins);

nresample = 1000;
% Resample the proportions in a loop
% For each resampling, store the fraction, compute the leas
yres = nan(nresample, npercentiles);
slp_res = nan(nresample, 1);
yres_tag = nan(nresample, 1);
% Bootstrapped estimate of the 
for kp = 1:nresample
    
    ysur = [];
    for jp = 1:npercentiles
        gpt = find( gp==jp );
        ysur = [ ysur {y( gpt( randi( length( gpt ), length( gpt ), 1 ) ) ) }];
    end
    coeff = polyfit([1:npercentiles], cellfun(@nanmean, ysur ), 1);
    slp_res(kp) = coeff(1);
    yres(kp,:) = cellfun(@nanmean, ysur );
    yres_tag(kp) = nanmean( ytag(randi( length(ytag), length(ytag), 1 )) );
end


data = arrayfun(@(x) nanmean(y(gp == x)), 1:npercentiles);

coeff = polyfit([1:npercentiles], data, 1);


reg = polyval(coeff, 1:npercentiles);
upper = polyval([prctile(slp_res, 97.5) coeff(2)], 1:npercentiles);
lower = polyval([prctile(slp_res, 2.5) coeff(2)], 1:npercentiles);


f = [upper flip(lower)];
figure
set(gcf,'Position', [913   715   345   336])
fill([1:npercentiles flip(1:npercentiles)], f, 'b', 'FaceAlpha', 0.2, 'linestyle', 'none' )
hold on
plot(1:npercentiles, reg, '--b', 'Linewidth', 2)

if strcmp(datatype, 'ripple')
    ylabel('Ripple cofiring', 'fontsize', 14)
elseif strcmp(datatype, 'theta')
    ylabel('Theta cofiring', 'fontsize', 14)
elseif strcmp(datatype, 'place')
    ylabel('Spatial ratemap correlation', 'fontsize', 14)
end

% Plot bootstrapped error bar
ci_pos = arrayfun(@(x) prctile( yres(:,x), 97.5 ), 1:npercentiles );
ci_neg = arrayfun(@(x) prctile( yres(:,x), 2.5 ), 1:npercentiles );
ci_pos_tag = prctile( yres_tag, 97.5 );
ci_neg_tag = prctile( yres_tag, 2.5 );
errorbar(1:npercentiles,data,abs( ci_neg-data ),abs( ci_pos-data ),'k+')
errorbar(npercentiles+1,nanmean(ytag),abs( ci_neg_tag-nanmean(ytag) ),abs( ci_pos_tag-nanmean(ytag) ),'k+')
xlim([0.05 npercentiles+1.5])

%%
ns = length(y);
% Group sizes
gp_ids_uniq = unique( gp ) ;
gp_sizes = arrayfun(@(x) sum( x == gp ), gp_ids_uniq);
t_surr = nan(1, nresample);
for kp = 1:nresample
    y_boot =  arrayfun(@(x) y(randi(ns, 1, gp_sizes(x))), 1:length(gp_sizes), 'UniformOutput', false);
    lm_boot = fitlm( 1:npercentiles, cellfun(@nanmean, y_boot ) );
    t_surr(kp) = lm_boot.Coefficients{'x1', 'tStat'};
end

figure
lm = fitlm( 1:npercentiles, data );
histogram(t_surr); 
xline( lm.Coefficients{'x1', 'tStat'} )
pval = sum( abs( t_surr ) > lm.Coefficients{'x1', 'tStat'} ) / nresample;

    
%% Load up behavior specific stuff

load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/analyze/final_code/revision2022_code/corrBeh_v7.mat')

%% Population vector decorrelation across left/right trial types
 
cc_all = [];

animal_exclude = [{'e16_3m2'} {'e15_13f1'}];
for kp = 1:length(basepaths_beh)
    fprintf('%d/%d\n',kp,length(basepaths_beh))
    basepath = alterPath( basepaths_beh{kp}, true );
    basename = bz_BasenameFromBasepath(basepath);
    
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), animal_exclude) )
        continue
    end
    
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    try
        load(fullfile(basepath, [basename '.firingMapsAvg_v2.cellinfo.mat']))
    catch
        load(fullfile(basepath, [basename '.firingMapsAvg_multSess.cellinfo.mat']))
    end
    pyr_id_log = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyrRatemaps = firingMaps(1).rateMaps(pyr_id_log);
    
    % pyrRatemaps = firingMaps.rateMaps( setdiff( firingMaps.UID, firingMaps.pyr_ind ) );
    % Combine the ratemaps
    ratemaps_left = [];
    ratemaps_right = [];
    for kp = 1:length(pyrRatemaps)
        ratemaps_left = [ ratemaps_left ; pyrRatemaps{kp}{1}  ];
        ratemaps_right = [ ratemaps_right ; pyrRatemaps{kp}{2}  ];
    end
    cc = nan(1,205);
    for kp = 1:205
        t = corr( [ ratemaps_left(:,kp) ratemaps_right(:,kp)] );
        cc(kp) = t(1,2);
    end
    cc_all = [cc_all ; cc];
        
end

%%

subplot(3,1,1)
imagesc( cc_all )
colormap(cool)
xlabel('Linearized position (cm)', 'fontsize', 16)
ylabel('Session #', 'fontsize', 16)
axis xy
colorbar
title('Population vector decorrelation')

subplot(3,1,2)


xv = 1:205;
fnt = [nanmean( cc_all )+std( cc_all )  flip( nanmean( cc_all )-std( cc_all ) )];
fill( [xv  flip(xv)], fnt, 'k', 'facealpha', 0.1, 'linestyle', 'none' )
hold on
%plot(tt_mu, '-b', 'linewidth', 1)
plot( nanmean( cc_all ), 'k', 'linewidth', 2 )
xlim([1 205])
colorbar

xline([74 111 185 222])

subplot(3,1,3)
cc_stem = cc_all(:,1:74); cc_stem = cc_stem(:);
cc_turn1 = cc_all(:,75:111); cc_turn1 = cc_turn1(:);
cc_arm = cc_all(:,112:185); cc_arm = cc_arm(:);
cc_turn2 = cc_all(:,186:205); cc_turn2 = cc_turn2(:);
boxplot([cc_stem ; cc_turn1 ; cc_arm ; cc_turn2], [ones(size(cc_stem)) ; 2*ones(size(cc_turn1)) ; 3*ones(size(cc_arm)) ; 4*ones(size(cc_turn2))],...
    'notch','on','whisker',3,'labels',{'stem','turn1', 'arm', 'turn2'})

a = [cc_stem ; cc_turn1 ; cc_arm ; cc_turn2];
b = [ones(size(cc_stem)) ; 2*ones(size(cc_turn1)) ; 3*ones(size(cc_arm)) ; 4*ones(size(cc_turn2))];

stats = boot_anova1(a,b,'classical', false)

%% Spatial ratemap correlation

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( pfcorr_ntag_all), [ 1:length(pfcorr_ntag_all) ] ./ length(pfcorr_ntag_all), '-r');
hold on
plot(sort( pfcorr_tagtag_all), [ 1:length(pfcorr_tagtag_all) ] ./ length(pfcorr_tagtag_all), '-b')
xlim( [-0.4 0.6] )
xlabel('Spatial ratemap correlation', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)


figure
set(gcf,'Position', [730   345   435   420])
boxplot([pfcorr_ntag_all ; pfcorr_tagtag_all], [ones(size(pfcorr_ntag_all)) ; 2*ones(size(pfcorr_tagtag_all))],'notch','on', 'whisker',inf,...
        'labels',{'different birthdate','same birthdate'})
ylim([-0.1725    0.2500])

%% Spatial ratemap correlation as function of shank distance

tagtag_mu = [ nanmean( pfcorr_tagtag_all( pfcorr_tagtag_shankdist == 0 ) );
      nanmean( pfcorr_tagtag_all( pfcorr_tagtag_shankdist == 1  ) ) ;
      nanmean( pfcorr_tagtag_all(pfcorr_tagtag_shankdist == 2 ) ) ];
tagtag_sem = [ sem( pfcorr_tagtag_all( pfcorr_tagtag_shankdist == 0  )) ; ...
               sem( pfcorr_tagtag_all( pfcorr_tagtag_shankdist == 1  ) ) ;
               sem( pfcorr_tagtag_all( pfcorr_tagtag_shankdist == 2  ) ) ];

ntag_mu = [ nanmean( pfcorr_ntag_all( pfcorr_ntag_shankdist == 0 ) ) ;
            nanmean( pfcorr_ntag_all( pfcorr_ntag_shankdist == 1  ) );
             nanmean( pfcorr_ntag_all( pfcorr_ntag_shankdist == 2 ) ) ];
ntag_sem = [ sem( pfcorr_ntag_all( pfcorr_ntag_shankdist == 0) ) ; ...
             sem( pfcorr_ntag_all( pfcorr_ntag_shankdist == 1) ) ; 
             sem( pfcorr_ntag_all( pfcorr_ntag_shankdist == 2) ) ];

         
 errorbar(1:3,tagtag_mu,tagtag_sem, '-s', 'MarkerSize',10,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','CapSize',18')
hold on
errorbar(1:3,ntag_mu,ntag_sem,'-s', 'MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',18)
xticks([1 2 3])
xticklabels({'0', '1', '2'})

set(gcf,'Position',[440 378 404 420])
legend('same birthdate', 'different birthdate')
xlabel('Shank distance', 'fontsize', 20)
ylabel('Spatial ratemap correlation', 'fontsize', 20)


rm = isnan( pfcorr_tagtag_right ) | isnan( pfcorr_tagtag_left ); 
pf_tagtag_right = pfcorr_tagtag_right(~rm);
pf_tagtag_left = pfcorr_tagtag_left(~rm);

rm = isnan( pfcorr_ntag_right ) | isnan( pfcorr_ntag_left );
pf_ntag_right = pfcorr_ntag_right(~rm);
pf_ntag_left = pfcorr_ntag_left(~rm);

%% Spatial ratemap correlation - individual animals

uniq_animals = unique(animal_tag);
figure
set(gcf,'Position', [59         631        2441         321])
for kp = 1:length( uniq_animals )
    
    subplot(1,length( uniq_animals ), kp)
    tt_fam = pfcorr_tagtag_all(cellfun(@(x) strcmp(x, uniq_animals{kp}), animal_tag),1);
    nt_fam = pfcorr_ntag_all(cellfun(@(x) strcmp(x, uniq_animals{kp}), animal_ntag),1);
    
    plot(sort(tt_fam), [ 1:length(tt_fam) ] ./length(tt_fam), 'b' )
    hold on
    plot(sort(nt_fam), [ 1:length(nt_fam) ] ./length(nt_fam), 'r' )
    xlim([-0.4 0.6])
    title(sprintf('%s, n=%d SBD', uniq_animals{kp}, length(tt_fam)),  'Interpreter', 'none')
    legend('same birthdate', 'different birthdate', 'location', 'best')
    
end

%%

pfcorr_by_shankdist_tag = [];
pfcorr_by_shankdist_ntag = [];
uniq_bp = unique( basenames_tagtag );
for kp = 1:length(uniq_bp)
    indic_tag = cellfun(@(x) strcmp(x, uniq_bp{kp}), basenames_tagtag);
    indic_ntag = cellfun(@(x) strcmp(x, uniq_bp{kp}), basenames_ntag);
    cc_tag = pfcorr_tagtag_all(indic_tag);
    cc_ntag = pfcorr_ntag_all(indic_ntag);
    shankd_tag = pfcorr_tagtag_shankdist(indic_tag);
    shankd_ntag = pfcorr_ntag_shankdist(indic_ntag);
    pfcorr_by_shankdist_tag = [pfcorr_by_shankdist_tag ; [ nanmean( cc_tag(shankd_tag==0) ) nanmean( cc_tag(shankd_tag==1) ) nanmean( cc_tag(shankd_tag==2) ) ] ];
    pfcorr_by_shankdist_ntag = [pfcorr_by_shankdist_ntag ; [ nanmean( cc_ntag(shankd_ntag==0) ) nanmean( cc_ntag(shankd_ntag==1) ) nanmean( cc_ntag(shankd_ntag==2) ) ] ];
end
%%
figure
set(gcf,'Position', [455   289   282   326])
plot( pfcorr_by_shankdist_tag', 'b' )
hold on
plot( nanmean( pfcorr_by_shankdist_tag ), 'k', 'LineWidth', 2 )
plot( pfcorr_by_shankdist_ntag', 'r' )
plot( nanmean( pfcorr_by_shankdist_ntag ), 'm', 'LineWidth', 2 )
ylim([-0.3 0.35])
%% Ratemap remap of SBD v. DBD across left and right trial types

minvals = 100;
step = 33;
edges = prctile([ pfcorr_tagtag_left ; pfcorr_tagtag_right ] ,0:step:100);
%edges = prctile([ pfcorr_ntag_left ; pfcorr_ntag_right ] ,0:step:100);
% Left
[cnts_tt_left, ~, ib_tt_left] = histcounts( pfcorr_tagtag_left, edges );
[cnts_nt_left, ~, ib_nt_left] = histcounts( pfcorr_ntag_left, edges );
% Right
[cnts_tt_right, ~, ib_tt_right] = histcounts( pfcorr_tagtag_right, edges );
[cnts_nt_right, ~, ib_nt_right] = histcounts( pfcorr_ntag_right, edges );

good_bins = find( cnts_tt_left > minvals );

tt_mu = []; tt_sem  = [];
nt_mu = []; nt_sem = [];
tt = [];
nt = [];
for kp = good_bins
   tt = [tt {[ pfcorr_tagtag_right( ib_tt_left == kp ) ; pfcorr_tagtag_left( ib_tt_right == kp ) ]}];
   nt = [nt {[ pfcorr_ntag_right( ib_nt_left == kp ) ; pfcorr_ntag_left( ib_nt_right == kp ) ]}]; 
   tt_mu = [ tt_mu ; nanmean( [ pfcorr_tagtag_right( ib_tt_left == kp ) ; pfcorr_tagtag_left( ib_tt_right == kp ) ] ) ]; tt_sem = [tt_sem ; nansem([ pfcorr_tagtag_right( ib_tt_left == kp ) ; pfcorr_tagtag_left( ib_tt_right == kp ) ]  )];
   nt_mu = [ nt_mu ; nanmean( [ pfcorr_ntag_right( ib_nt_left == kp ) ; pfcorr_ntag_left( ib_nt_right == kp ) ] ) ]; nt_sem = [nt_sem ; nansem( [ pfcorr_ntag_right( ib_nt_left == kp ) ; pfcorr_ntag_left( ib_nt_right == kp ) ] )];
end

xv = [1:length(good_bins)]';
ftt = [tt_mu+tt_sem ; flip( tt_mu-tt_sem )];
fnt = [nt_mu+nt_sem ; flip( nt_mu-nt_sem )];

figure
set(gcf,'Position',[440   377   306   420])
fill( [xv ; flip(xv)], ftt, 'b', 'facealpha', 0.1, 'linestyle', 'none' )
hold on
fill( [xv ; flip(xv)], fnt, 'r', 'facealpha', 0.1, 'linestyle', 'none' )
plot(tt_mu, '-b', 'linewidth', 1)
plot(nt_mu, '-r', 'linewidth', 1)
% errorbar(xv,tt_mu,tt_sem, '-s', 'MarkerSize',8,...
% 'MarkerEdgeColor','blue','MarkerFaceColor','blue','CapSize',18')
% hold on
% errorbar(xv,nt_mu,nt_sem, '-s', 'MarkerSize',8,...
% 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',18')
% axis square
% xticks(1:length(good_bins))

xlabel('Ratemap corr. on Trial A (percentiles)', 'fontsize', 14)
ylabel('Ratemap corr. on Trial B', 'fontsize', 14)


%%

%vals = [ pfcorr_tag_all(:) ; pfcorr_ntag_all(:) ];
vals = [ pfcorr_tagtag_left ; pfcorr_tagtag_right ];
nanvals = isnan(vals);
pd = fitdist(vals(~nanvals),'kernel','Width',0.05);
x = min( vals(~nanvals) ):0.01:max( vals(~nanvals) );
y = pdf(pd,x);

figure
set(gcf,'Position', [1000         918         560         281])
histogram(vals(~nanvals),25, 'Normalization', 'pdf','linestyle','none' );
hold on

plot(x,y  ,'Color','k','LineStyle','--')
xlim([min(x) max(x)])

ylim_vals = ylim;
xlim_vals = xlim;

boundaries = edges(2:end-1); % we use edges defined by the plot for visualization
b1 = rectangle('Position', [ xlim_vals(1), ylim_vals(1), boundaries(1) - xlim_vals(1), ylim_vals(2) ], 'FaceColor','b', 'linestyle', 'none');
b1.FaceColor(4) = 0.4;
b2 = rectangle('Position', [ boundaries(1), ylim_vals(1), boundaries(2) - boundaries(1), ylim_vals(2) ], 'FaceColor','r', 'linestyle', 'none');
b2.FaceColor(4) = 0.6;
b3 = rectangle('Position', [ boundaries(2), ylim_vals(1), xlim_vals(2) - boundaries(2), ylim_vals(2) ], 'FaceColor','g', 'linestyle', 'none');
b3.FaceColor(4) = 0.8;
ylim(ylim_vals)
% Now divide this distribution into percentiles
xlabel('Ratemap corr. in FAM', 'fontsize', 14)
ylabel('Probability density', 'fontsize', 14)

%%

y= [pfcorr_ntag_all( pfcorr_ntag_shankdist == 0) ; ...
    pfcorr_ntag_all( pfcorr_ntag_shankdist == 1) ; ...
    pfcorr_ntag_all( pfcorr_ntag_shankdist == 2) ; ...
    pfcorr_tagtag_all( pfcorr_tagtag_shankdist == 0) ; ...
    pfcorr_tagtag_all( pfcorr_tagtag_shankdist == 1) ; ...
    pfcorr_tagtag_all( pfcorr_tagtag_shankdist == 2)];
g1 = [repmat({'ntag'}, sum(pfcorr_ntag_shankdist == 0)+sum(pfcorr_ntag_shankdist == 1)+sum(pfcorr_ntag_shankdist == 2),1); ...
      repmat({'tag'},sum(pfcorr_tagtag_shankdist == 0)+sum(pfcorr_tagtag_shankdist == 1)+sum(pfcorr_tagtag_shankdist == 2),1)];
g2 = [zeros( sum(pfcorr_ntag_shankdist == 0), 1 ) ; ones( sum(pfcorr_ntag_shankdist == 1), 1 ) ; 2.*ones( sum(pfcorr_ntag_shankdist == 2), 1 ) ;...
      zeros( sum(pfcorr_tagtag_shankdist == 0), 1 ) ; ones( sum(pfcorr_tagtag_shankdist == 1), 1 ) ; 2.*ones( sum(pfcorr_tagtag_shankdist == 2), 1 ) ];

[p,tbl,stats] = anovan(y,{g1,g2},'model','interaction',...
    'varnames',{'g1','g2'});
[c,~,~,gnames] = multcompare(stats,'Dimension',[1 2]);

[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))]

%%
y = [ cell2mat( tt') ; cell2mat(nt') ];
gp1 = [repmat({'tag'}, size(cell2mat( tt'))) ; repmat({'ntag'}, size(cell2mat( nt'))) ];
gp2 = [ ones(size( tt{1} ) ) ; 2.*ones(size( tt{2} )) ; 3.*ones(size( tt{3} )) ; ...];
        ones(size( nt{1} ) ) ; 2.*ones(size( nt{2} )) ; 3.*ones(size( nt{3} ))  ];

[p,tbl,stats] = anovan(y,{gp1,gp2},'model','interaction',...
    'varnames',{'g1','g2'});
[c,~,~,gnames] = multcompare(stats,'Dimension',[1 2]);

[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))]


    