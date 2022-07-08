

%load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac_rev1.mat")
load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/basepaths_mac_rev1.mat')

metrics_cat = struct('int_ratemod',[],'spktrans_gamma',[], 'spktrans_mean',[], 'spktrans_splitISI', [], 'uid_pair', [], 'basename', []);
for kp = 1:length(basepaths_all)
    fprintf('%d/%d\n', kp, length(basepaths_all))
    basepath = alterPath(basepaths_all{kp}, true);
    basename = bz_BasenameFromBasepath(basepath);
    v = load(fullfile(basepath, 'monoConnections_pretag.mat') );
    metrics_cat.int_ratemod = [ metrics_cat.int_ratemod ; v.int_ratemod ];
    metrics_cat.spktrans_gamma = [ metrics_cat.spktrans_gamma ; v.spktrans_gamma ];
    metrics_cat.spktrans_mean = [ metrics_cat.spktrans_mean ; v.spktrans_mean ];
    metrics_cat.spktrans_splitISI = [ metrics_cat.spktrans_splitISI ; v.spktrans_splitISI ];
    metrics_cat.uid_pair = [ metrics_cat.uid_pair ; v.uid_pair ];
    metrics_cat.basename = [ metrics_cat.basename ; repmat({basename}, length(v.spktrans_mean), 1)];
end


%%
brown = [133 87 34]/254;
e13_tag = cellfun(@(x) ~isempty(regexp(x, 'e13', 'once')), metrics_cat.basename );
e14_tag = cellfun(@(x) ~isempty(regexp(x, 'e14', 'once')), metrics_cat.basename );
e15_tag = cellfun(@(x) ~isempty(regexp(x, 'e15', 'once')), metrics_cat.basename ) & ~cellfun(@(x) ~isempty(regexp(x, 'e15_13f1', 'once')), metrics_cat.basename );
e16_tag = cellfun(@(x) ~isempty(regexp(x, 'e16', 'once')), metrics_cat.basename ) & ~cellfun(@(x) ~isempty(regexp(x, 'e16_3m2', 'once')), metrics_cat.basename );

mod = logical( metrics_cat.int_ratemod(:,1) );
mod = metrics_cat.int_ratemod(:,2) > 0 & metrics_cat.int_ratemod(:,3) < 0.05;
% mod = true( size( metrics_cat.int_ratemod(:,1) ) );

isi = [ logspace(log10(5),log10(1000),15)/1000 ];


figure
hold on

spktrans_e13 = metrics_cat.spktrans_splitISI(e13_tag ,:);
mu_e13 = nanmean(spktrans_e13 );
se = nanstd(spktrans_e13) ./ sqrt( sum( ~isnan( spktrans_e13 ) ) );
fe13 = [mu_e13+se flip( mu_e13-se )];

spktrans_e14 = metrics_cat.spktrans_splitISI(e14_tag ,:);
mu_e14 = nanmean(spktrans_e14);
se = nanstd(spktrans_e14) ./ sqrt( sum( ~isnan( spktrans_e14 ) ) );
fe14 = [mu_e14+se flip( mu_e14-se )];

spktrans_e15 = metrics_cat.spktrans_splitISI(e15_tag ,:);
mu_e15 = nanmean(spktrans_e15);
se = nanstd(spktrans_e15) ./ sqrt( sum( ~isnan( spktrans_e15 ) ) );
fe15 = [mu_e15+se flip( mu_e15-se )];

spktrans_e16 =metrics_cat.spktrans_splitISI(e16_tag ,:);
mu_e16 = nanmean(spktrans_e16);
se = nanstd(spktrans_e16) ./ sqrt( sum( ~isnan( spktrans_e16 ) ) );
fe16 = [mu_e16+se flip( mu_e16-se )];



plot(isi, mu_e13, 'Color', brown, 'linewidth', 2)
plot(isi, mu_e14, 'r', 'linewidth', 2)
plot(isi, mu_e15, 'b', 'linewidth', 2)
plot(isi, mu_e16, 'k', 'linewidth', 2)
fill( [isi flip(isi)], fe13, brown, 'facealpha', 0.1, 'linestyle', 'none' )
fill( [isi flip(isi)], fe14, 'r', 'facealpha', 0.1, 'linestyle', 'none' )
fill( [isi flip(isi)], fe15, 'b', 'facealpha', 0.1, 'linestyle', 'none' )
fill( [isi flip(isi)], fe16, 'k', 'facealpha', 0.1, 'linestyle', 'none' )

legend('E13.5', 'E14.5', 'E15.5', 'E16.5', 'location', 'best')


set(gca,'xscale', 'log')
set(gca,'yscale', 'log')

xlabel('presynaptic ISI (s)', 'fontsize', 20)
ylabel('Spike transmission probability', 'fontsize', 20)

%%
set(gcf,'Position', [1000         918         382         420])
a = [ metrics_cat.spktrans_gamma(e13_tag) ; metrics_cat.spktrans_gamma(e14_tag) ; metrics_cat.spktrans_gamma(e15_tag) ; metrics_cat.spktrans_gamma(e16_tag) ];
b = [ ones(sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1) ];
boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Spk. transmission probability')
set(gca,'yscale', 'log')
ylim([0.008    0.05])

stats = boot_anova1(a,b, 'classical', false)

%%

a = [ metrics_cat.spktrans_mean(e13_tag ) ; metrics_cat.spktrans_mean(e14_tag) ; metrics_cat.spktrans_mean(e15_tag) ; metrics_cat.spktrans_mean(e16_tag) ];
b = [ ones(sum( e13_tag ), 1) ; 2.*ones(sum( e14_tag ), 1) ; 3.*ones(sum( e15_tag ), 1) ; 4.*ones(sum( e16_tag ), 1) ];
boxplot(a,b,'labels', {'E13.5', 'E14.5', 'E15.5', 'E16.5'}, 'notch', 'on', 'whisker', 1000)
ylabel('Spike transmission probability')
set(gca,'yscale', 'log')
ylim([0.008    0.05])

stats = boot_anova1(a,b, 'classical', false)
