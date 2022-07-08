

load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/analyze/final_code/revision2022/corrHome_v7.mat')
%% Switch out 


idx = converge_nolab > 0;
cnt = converge_nolab(idx);

% Switch between "ripcorr_nolab" & 

if ripflag
    var_ntag = ripcorr_nolab(idx,:);
    var_tagtag = ripcorr_tagtag;
else
    var_ntag = thetacorr_nolab(idx,:);
    var_tagtag = thetacorr_tagtag;
end
binN = 12;
perc_step = 20;
perc_all = 0:perc_step:100;
prc_vals = prctile( cnt, perc_all );

[~, ~, ip] = histcounts( cnt, prc_vals );

figure
set(gcf,'Position',[890   780   341   420])
err = arrayfun(@(x) sem( var_ntag( ip == x) ), 1:length( perc_all ) - 1 );
mu = arrayfun(@(x) nanmean( var_ntag( ip == x) ), 1:length( perc_all ) - 1 );
f = [mu+err flip(mu-err)];
fill([1:length( perc_all ) - 1 flip(1:length( perc_all ) - 1)], f , 'b', 'facealpha', 0.2, 'linestyle', 'none')
hold on
plot(1:length( perc_all ) - 1, mu , '-b' )
hold on
%plot( 5, tht(kp), 'kx')
errorbar(length( perc_all ) - .5, nanmean( var_tagtag ), sem(var_tagtag)) 
plot( length( perc_all ) - .5, nanmean( var_tagtag ), '.r')
xticks( 1:length( perc_all ) - 1)
xlim([ 1 length( perc_all )+1 ])
xline(length( perc_all )-1)
xlabel('Convergence index percentile')
if ripflag
    ylabel('Co-firing in SPW-Rs')
else
    ylabel('Theta cycle co-firing in SPW-Rs')
end

shg

figure
histogram( cnt, binN, 'linestyle', 'none' )
hold on
xline( prc_vals(2:end-1) )

%% Alternative visualization using box plots

figure
set(gcf,'Position', [440   377   312   420])
boxplot([ var_ntag ; var_tagtag], [ ip ; (max(ip)+1).*ones(size(var_tagtag))], 'whisker',1000, 'notch', 'on')
ylim([-.05 0.1])

xticks( 1:length( perc_all ) - 1)
xlabel('Convergence index percentile')
ylabel('Cofiring in SPW-Rs')

