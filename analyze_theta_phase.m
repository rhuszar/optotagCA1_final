
%% Average the LFP

load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/e16/e16_3m1/e16_3m1_210111/e16_3m1_210111.thetaCycles.events.mat')
load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/e16/e16_3m1/e16_3m1_210111/e16_3m1_210111.ripples.events.mat')

passband = [6 14];
data = ripples.detectorinfo.detectionparms.lfp;
timestamps  = [ 0:(1/1250):length( data ) * (1/1250) ]'; timestamps(end) = [];
lfp = struct();
lfp.data = data;
lfp.timestamps = timestamps;
lfp.channels = 1;

lfp_th = bz_Filter( lfp, 'passband', passband );

% Pick the onsets of detected a subset of theta cycles
% Find nearest time in lfp, and go in time / phase
cyc = Theta.cycles( Theta.run,1 );
cyc_subset = cyc(1:5:end);
%%
dat = [];
for ip = 1:length(cyc_subset)
    [~, ind] = min( abs( timestamps - cyc_subset(ip) ) );
    dat = [ dat ; data((ind-200):( ind+380-1 ))' ];
end

%%
basepaths_all = [ basepaths_e13_26m1 ; basepaths_e15_13f1 ; basepaths_e16_3m2 ];
thetaPhase_chan = [123 122 116 116 27 27 66 67 31 31 2 2];
%%
%load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\thetaPhase_chan.mat")
% failed = [5     8     9    10    11    12    15    16    17    29    30    31    32    33    37    38    41    45    46    48    51    52    54    56    59    60    61    63    64    66    69    75    79    82    87    88    90    92];
failed = [];
% for kp = failed
for kp = 1:length(basepaths_all) 
%for kp = 29:length(basepaths_all)
   fprintf('%d/%d\n', kp,length(basepaths_all) )
   basepath_remote = alterPath( basepaths_all{kp} , false );
   basepath_local = alterPath( basepaths_all{kp} , true );
   basename = bz_BasenameFromBasepath( basepath_remote );
   
   % if exist(fullfile(basepath_remote, [basename '.lfp'])) && exist( fullfile(basepath_local, 'SleepScoreLFP.mat') )
   try
       lfp = bz_GetLFP(thetaPhase_chan(kp),'basepath', basepath_remote);
       save(fullfile(basepath_local, 'SleepScoreLFP_v1.mat'), 'lfp')
       delete(fullfile(basepath_local, 'SleepScoreLFP.mat'))
   %else
   catch
       failed = [failed ; kp];
   end
end

%% Get theta modulation

for kp = 1:length(basepaths_all)
    fprintf('%d/%d\n', kp,length(basepaths_all) )
    basepath_local = alterPath( basepaths_all{kp} , true );
    getThetaMod(basepath_local)
end

%% Get autocorrelogram frequency
for kp = 1:length(basepaths_all)
    fprintf('%d/%d\n', kp,length(basepaths_all) )
    basepath_local = alterPath( basepaths_all{kp} , true );
    getUnitThetaFreq(basepath)
end
     
%%
isLocal = true;
alpha = .01;
failed = [];
phase_wake = []; phase_rem = [];
mDepth_wake = []; mDepth_rem = [];
sig_wake = []; sig_rem = [];
phaseRate_rem = []; phaseRate_wake = [];
tag_indicator = [];
session_id = []; 
acdata_freq = []; acdata_power = [];

for kp = 1:length( basepaths_all )
    
    fprintf('%d/%d\n', kp, length( basepaths_all ))
    basepath = alterPath(basepaths_all{kp}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
%     try
        load( fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']) );
       %  load( fullfile(basepath, [basename '.ThetaModulation_hc.cellinfo.mat']) );
       load( fullfile(basepath, [basename '.ThetaModulation.cellinfo.mat']) );
       % load( fullfile(basepath, [basename '.ACData.cellinfo.mat']) );
%     catch
%         failed = [failed ; kp];
%         continue
%     end

    tmp = strsplit( basename, '_' );
    
    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id & pyr_id;
    
    % Gather all the pyramidal cells, keep track of which are tagged / theta modulated
    tag = arrayfun(@(x) find( x == find( pyr_id ) ), find(pyr_opto_id)) ;
    t = false( sum(pyr_id), 1 ); t(tag) = true;
    % Store the phases / depth / unit+lfp frequency
    phase_wake = [phase_wake ; ThetaMod.th_phaselock_wake.phasestats.m( pyr_id )' ];
    mDepth_wake = [mDepth_wake ; ThetaMod.th_phaselock_wake.phasestats.r( pyr_id )' ];
    
    ave_fr = cell2mat( cellfun(@(x) nanmean(x(11:end,:)), ThetaMod.th_powerphasefr_wake.ratemap(pyr_id), 'UniformOutput', false)' );
    phaseRate_wake = [phaseRate_wake ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
    
    sig_wake = [sig_wake ; ThetaMod.th_phaselock_wake.phasestats.p( pyr_id )' ];
    tag_indicator = [tag_indicator ; t];
    session_id= [session_id ; repmat(tmp(1), sum( pyr_id ), 1)];
    % Autocorrelogram features
%     acdata_freq = [acdata_freq ; [ACData.unitfreq' ACData.lfpfreq']];
%     acdata_power = [acdata_power ; ACData.unitpower'];
    % Get theta modulated PYRs in REM
    if ~isempty( ThetaMod.th_phaselock_rem.phasestats )
        % Store the phases / depth
        phase_rem = [phase_rem ; ThetaMod.th_phaselock_rem.phasestats.m( pyr_id )' ];
        mDepth_rem = [mDepth_rem ; ThetaMod.th_phaselock_rem.phasestats.r( pyr_id )' ];
        sig_rem = [sig_rem ; ThetaMod.th_phaselock_rem.phasestats.p( pyr_id )' ];
        
        ave_fr = cell2mat( cellfun(@(x) nanmean(x(11:end,:)), ThetaMod.th_powerphasefr_rem.ratemap(pyr_id), 'UniformOutput', false)' );
        phaseRate_rem = [phaseRate_rem ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        
    else
         phase_rem = [phase_rem ; nan(sum(pyr_id),1) ];
         mDepth_rem = [mDepth_rem ; nan(sum(pyr_id),1)];
         sig_rem = [sig_rem ; nan(sum(pyr_id),1) ];
         phaseRate_rem = [ phaseRate_rem ; nan(sum(pyr_id),40) ];
    end
end
tag_indicator = logical( tag_indicator );
%% Theta phase distribution of all cells ; split by birthdate
alpha = 0.01;
brown = [133 87 34]/254;
% Identify birthdated neurons
e13_tag = cellfun(@(x) strcmp(x, 'e13'), session_id ) & tag_indicator; 
e14_tag = cellfun(@(x) strcmp(x, 'e14'), session_id ) & tag_indicator; 
e15_tag = cellfun(@(x) strcmp(x, 'e15'), session_id ) & tag_indicator; 
e16_tag = cellfun(@(x) strcmp(x, 'e16'), session_id ) & tag_indicator;

% Untagged PYR theta phases
wake_untagged = [ rad2deg( phase_wake(~tag_indicator & sig_wake < alpha) ) ; 360+rad2deg( phase_wake(~tag_indicator & sig_wake < alpha) ) ];
rem_untagged = [ rad2deg( phase_rem(~tag_indicator & sig_rem < alpha) ) ; 360+rad2deg( phase_rem(~tag_indicator & sig_rem < alpha) ) ];

% Tagged PYR theta phases
e13_wake_tagged = [ rad2deg( phase_wake(e13_tag & sig_wake < alpha) ) ; 360+rad2deg( phase_wake(e13_tag & sig_wake < alpha) ) ]; 
e13_rem_tagged = [ rad2deg( phase_rem(e13_tag & sig_rem < alpha) ) ; 360+rad2deg( phase_rem(e13_tag & sig_rem < alpha) ) ];
e14_wake_tagged = [ rad2deg( phase_wake(e14_tag & sig_wake < alpha) ) ; 360+rad2deg( phase_wake(e14_tag & sig_wake < alpha) ) ]; 
e14_rem_tagged = [ rad2deg( phase_rem(e14_tag & sig_rem < alpha) ) ; 360+rad2deg( phase_rem(e14_tag & sig_rem < alpha) ) ];
e15_wake_tagged = [ rad2deg( phase_wake(e15_tag & sig_wake < alpha) ) ; 360+rad2deg( phase_wake(e15_tag & sig_wake < alpha) ) ]; 
e15_rem_tagged = [ rad2deg( phase_rem(e15_tag & sig_rem < alpha) ) ; 360+rad2deg( phase_rem(e15_tag & sig_rem < alpha) ) ];
e16_wake_tagged = [ rad2deg( phase_wake(e16_tag & sig_wake < alpha) ) ; 360+rad2deg( phase_wake(e16_tag & sig_wake < alpha) ) ]; 
e16_rem_tagged = [ rad2deg( phase_rem(e16_tag & sig_rem < alpha) ) ; 360+rad2deg( phase_rem(e16_tag & sig_rem < alpha) ) ];
nb = 18;

% Bin the untagged neurons in each state 
[prob_wake, edges_wake] = histcounts(wake_untagged,nb,'Normalization', 'probability');
[prob_rem, edges_rem] = histcounts(rem_untagged,nb,'Normalization', 'probability');

% Bin the tagged neurons in wake
[e13_prob_wake, e13_edges_wake] = histcounts(e13_wake_tagged,nb,'Normalization', 'probability'); 
e13_edges_wake = [0 e13_edges_wake e13_edges_wake(end)+diff(e13_edges_wake(1:2))]; e13_prob_wake = [0 e13_prob_wake 0];
[e14_prob_wake, e14_edges_wake] = histcounts(e14_wake_tagged,nb,'Normalization', 'probability'); 
e14_edges_wake = [0 e14_edges_wake e14_edges_wake(end)+diff(e14_edges_wake(1:2))]; e14_prob_wake = [0 e14_prob_wake 0];
[e15_prob_wake, e15_edges_wake] = histcounts(e15_wake_tagged,nb,'Normalization', 'probability'); 
e15_edges_wake = [0 e15_edges_wake e15_edges_wake(end)+diff(e15_edges_wake(1:2))]; e15_prob_wake = [0 e15_prob_wake 0];
[e16_prob_wake, e16_edges_wake] = histcounts(e16_wake_tagged,nb,'Normalization', 'probability'); 
e16_edges_wake = [0 e16_edges_wake e16_edges_wake(end)+diff(e16_edges_wake(1:2))]; e16_prob_wake = [0 e16_prob_wake 0];

% Bin the tagged neurons in REM
[e13_prob_rem, e13_edges_rem] = histcounts(e13_rem_tagged,nb,'Normalization', 'probability'); 
e13_edges_rem = [0 e13_edges_rem e13_edges_rem(end)+diff(e13_edges_rem(1:2))]; e13_prob_rem = [0 e13_prob_rem 0];
[e14_prob_rem, e14_edges_rem] = histcounts(e14_rem_tagged,nb,'Normalization', 'probability'); 
e14_edges_rem = [0 e14_edges_rem e14_edges_rem(end)+diff(e14_edges_rem(1:2))]; e14_prob_rem = [0 e14_prob_rem 0];
[e15_prob_rem, e15_edges_rem] = histcounts(e15_rem_tagged,nb,'Normalization', 'probability'); 
e15_edges_rem = [0 e15_edges_rem e15_edges_rem(end)+diff(e15_edges_rem(1:2))]; e15_prob_rem = [0 e15_prob_rem 0];
[e16_prob_rem, e16_edges_rem] = histcounts(e16_rem_tagged,nb,'Normalization', 'probability'); 
e16_edges_rem = [0 e16_edges_rem e16_edges_rem(end)+diff(e16_edges_rem(1:2))]; e16_prob_rem = [0 e16_prob_rem 0];

% WAKE
subplot(3,1,1)
plot(edges_wake(1:end-1)+diff(edges_wake)./2, prob_wake, '--'); hold on
plot(e13_edges_wake(1:end-1)+diff(e13_edges_wake)./2, e13_prob_wake, 'Color', brown);
plot(e14_edges_wake(1:end-1)+diff(e14_edges_wake)./2, e14_prob_wake, '-r');
plot(e15_edges_wake(1:end-1)+diff(e15_edges_wake)./2, e15_prob_wake, 'b');
plot(e16_edges_wake(1:end-1)+diff(e16_edges_wake)./2, e16_prob_wake, 'k');
xticks([0:90:720])
xlim([0 720])

% REM
subplot(3,1,2)
plot(edges_rem(1:end-1)+diff(edges_rem)./2, prob_rem, '--'); hold on
plot(e13_edges_rem(1:end-1)+diff(e13_edges_rem)./2, e13_prob_rem, 'Color', brown);
plot(e14_edges_rem(1:end-1)+diff(e14_edges_rem)./2, e14_prob_rem, '-r');
plot(e15_edges_rem(1:end-1)+diff(e15_edges_rem)./2, e15_prob_rem, 'b');
plot(e16_edges_rem(1:end-1)+diff(e16_edges_rem)./2, e16_prob_rem, 'k');
xticks([0:90:720])
xlim([0 720])
xlabel('Preferred theta phase', 'fontsize', 18)
ylabel('Proportion of cells', 'fontsize', 18)

load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\figs\shortFormat\Figure4\thetaPhase_example.mat")
% load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/analyze/figs/Figure4/thetaPhase_example.mat')
subplot(3,1,3)
plot( theta_phase, lfp )
xlim([0 720])

%% Theta phase distribution of all cells ; tagged versus untagged
alpha = 0.01;
brown = [133 87 34]/254;
% Identify birthdated neurons
e13_tag = cellfun(@(x) strcmp(x, 'e13'), session_id ) & tag_indicator; 
e14_tag = cellfun(@(x) strcmp(x, 'e14'), session_id ) & tag_indicator; 
e15_tag = cellfun(@(x) strcmp(x, 'e15'), session_id ) & tag_indicator; 
e16_tag = cellfun(@(x) strcmp(x, 'e16'), session_id ) & tag_indicator;

% Untagged PYR theta phases
wake_untagged = [ rad2deg( phase_wake(~tag_indicator & sig_wake < alpha) ) ; 360+rad2deg( phase_wake(~tag_indicator & sig_wake < alpha) ) ];
rem_untagged = [ rad2deg( phase_rem(~tag_indicator & sig_rem < alpha) ) ; 360+rad2deg( phase_rem(~tag_indicator & sig_rem < alpha) ) ];

% Tagged PYR theta phases
wake_tagged = [ rad2deg( phase_wake(tag_indicator & sig_wake < alpha) ) ; 360+rad2deg( phase_wake(tag_indicator & sig_wake < alpha) ) ];
rem_tagged = [ rad2deg( phase_rem(tag_indicator & sig_rem < alpha) ) ; 360+rad2deg( phase_rem(tag_indicator & sig_rem < alpha) ) ];

nb = 18;
% Bin the untagged neurons in each state 
[prob_wake_nt, edges_wake_nt] = histcounts(wake_untagged,nb,'Normalization', 'probability');
[prob_rem_nt, edges_rem_nt] = histcounts(rem_untagged,nb,'Normalization', 'probability');
% Bin the tagged neurons in each state 
[prob_wake_tt, edges_wake_tt] = histcounts(wake_tagged,nb,'Normalization', 'probability');
[prob_rem_tt, edges_rem_tt] = histcounts(rem_tagged,nb,'Normalization', 'probability');



% WAKE
subplot(3,1,1)
plot(edges_wake_nt(1:end-1)+diff(edges_wake_nt)./2, prob_wake_nt, '-r'); hold on
plot(edges_wake_tt(1:end-1)+diff(edges_wake_tt)./2, prob_wake_tt, '-b')

xticks([0:90:720])
xlim([0 720])

% REM
subplot(3,1,2)
plot(edges_rem_nt(1:end-1)+diff(edges_rem_nt)./2, prob_rem_nt, '-r'); hold on
plot(edges_rem_tt(1:end-1)+diff(edges_rem_tt)./2, prob_rem_tt, '-b'); 
xticks([0:90:720])
xlim([0 720])
xlabel('Preferred theta phase', 'fontsize', 18)
ylabel('Proportion of cells', 'fontsize', 18)

load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\figs\shortFormat\Figure4\thetaPhase_example.mat")
subplot(3,1,3)
plot( theta_phase, lfp* 0.195 )
xlim([0 720])

%% Theta depth modulation - mean resultant vector
%  Across brain states
e13_depthMod_wake = mDepth_wake ( e13_tag & sig_wake < alpha );
e14_depthMod_wake = mDepth_wake ( e14_tag & sig_wake < alpha );
e15_depthMod_wake = mDepth_wake ( e15_tag & sig_wake < alpha );
e16_depthMod_wake = mDepth_wake ( e16_tag & sig_wake < alpha );


e13_depthMod_rem = mDepth_rem ( e13_tag & sig_rem < alpha );
e14_depthMod_rem = mDepth_rem ( e14_tag & sig_rem < alpha );
e15_depthMod_rem = mDepth_rem ( e15_tag & sig_rem < alpha );
e16_depthMod_rem = mDepth_rem ( e16_tag & sig_rem < alpha );


figure
subplot(2,1,1)
boxplot([e13_depthMod_wake ; e14_depthMod_wake ; e15_depthMod_wake ; e16_depthMod_wake], ...
    [ones(size(e13_depthMod_wake)) ; 2*ones(size(e14_depthMod_wake)) ; 3*ones(size(e15_depthMod_wake)) ; 4*ones(size(e16_depthMod_wake))],'notch','on', 'whisker',inf)
ylim([0.05 0.17])
xticks(1:4)
xticklabels({'E13.5', 'E14.5', 'E15.5', 'E16.6'})
ylabel('Theta depth modulation', 'fontsize', 18)
xlabel('Birthdate', 'fontsize', 18)
title('WAKE','fontsize', 18)

subplot(2,1,2)
boxplot([e13_depthMod_rem ; e14_depthMod_rem ; e15_depthMod_rem ; e16_depthMod_rem], ...
    [ones(size(e13_depthMod_rem)) ; 2*ones(size(e14_depthMod_rem)) ; 3*ones(size(e15_depthMod_rem)) ; 4*ones(size(e16_depthMod_rem))],'notch','on', 'whisker',inf)
xticks(1:4)
ylim([0.05 0.25])
xticklabels({'E13.5', 'E14.5', 'E15.5', 'E16.6'})
ylabel('Theta depth modulation', 'fontsize', 18)
xlabel('Birthdate', 'fontsize', 18)
title('REM','fontsize', 18)

%% Unit oscillation frequency 

e13_thetaFreq = acdata_freq( e13_tag & sig_wake < alpha,: );
e14_thetaFreq = acdata_freq( e14_tag & sig_wake < alpha,: );
e15_thetaFreq = acdata_freq( e15_tag & sig_wake < alpha,: );
e16_thetaFreq = acdata_freq( e16_tag & sig_wake < alpha,: );


%% Distribution of unit oscillaton frequencies
figure
set(gcf,'Position', [2002 6 584 703])
subplot(2,1,1)
boxplot([e13_thetaFreq(:,1) ; e14_thetaFreq(:,1) ; e15_thetaFreq(:,1) ; e16_thetaFreq(:,1)], ...
    [ones(size(e13_thetaFreq(:,1) )) ; 2*ones(size(e14_thetaFreq(:,1))) ; 3*ones(size(e15_thetaFreq(:,1))) ; 4*ones(size(e16_thetaFreq(:,1)))],'notch','on', 'whisker',inf)
xticks(1:4)
ylim([7.5 10])
xticklabels({'E13.5', 'E14.5', 'E15.5', 'E16.6'})
ylabel('Unit oscillation frequency', 'fontsize', 20)

subplot(2,1,2)
boxplot([-diff(e13_thetaFreq')' ; -diff(e14_thetaFreq')' ; -diff(e15_thetaFreq')' ; -diff(e16_thetaFreq')'], ...
    [ones(size(e13_thetaFreq(:,1) )) ; 2*ones(size(e14_thetaFreq(:,1))) ; 3*ones(size(e15_thetaFreq(:,1))) ; 4*ones(size(e16_thetaFreq(:,1)))],'notch','on', 'whisker',inf)
ylabel('Difference between unit and LFP theta frequency', 'fontsize', 18)
ylim([-1 2])
xticks(1:4)
xticklabels({'E13.5', 'E14.5', 'E15.5', 'E16.6'})

%%

e13_tag = cellfun(@(x) strcmp(x, 'e13'), session_id ) & tag_indicator; 
e14_tag = cellfun(@(x) strcmp(x, 'e14'), session_id ) & tag_indicator; 
e15_tag = cellfun(@(x) strcmp(x, 'e15'), session_id ) & tag_indicator; 
e16_tag = cellfun(@(x) strcmp(x, 'e16'), session_id ) & tag_indicator;

th_wakephase_e13 = [ sig_wake(e13_tag) phase_wake(e13_tag)]; th_remphase_e13 = [ sig_rem(e13_tag) phase_rem(e13_tag)]; 
th_wakephase_e14 = [ sig_wake(e14_tag) phase_wake(e14_tag)]; th_remphase_e14 = [ sig_rem(e14_tag) phase_rem(e14_tag)]; 
th_wakephase_e15 = [ sig_wake(e15_tag) phase_wake(e15_tag)]; th_remphase_e15 = [ sig_rem(e15_tag) phase_rem(e15_tag)]; 
th_wakephase_e16 = [ sig_wake(e16_tag) phase_wake(e16_tag)]; th_remphase_e16 = [ sig_rem(e16_tag) phase_rem(e16_tag)];

th_wake_e13 = phaseRate_wake(e13_tag,:);
th_wake_e14 = phaseRate_wake(e14_tag,:);
th_wake_e15 = phaseRate_wake(e15_tag,:);
th_wake_e16 = phaseRate_wake(e16_tag,:);

th_rem_e13 = phaseRate_rem(e13_tag,:);
th_rem_e14 = phaseRate_rem(e14_tag,:);
th_rem_e15 = phaseRate_rem(e15_tag,:);
th_rem_e16 = phaseRate_rem(e16_tag,:);

%%

histogram( rad2deg( th_wakephase_e13( th_wakephase_e13(:,1) < .05, 2) ), 'Normalization', 'probability' )
hold on
histogram( rad2deg( th_wakephase_e14( th_wakephase_e14(:,1) < .05, 2) ), 'Normalization', 'probability' )
histogram( rad2deg( th_wakephase_e15( th_wakephase_e15(:,1) < .05, 2) ), 'Normalization', 'probability' )
histogram( rad2deg( th_wakephase_e16( th_wakephase_e16(:,1) < .05, 2) ), 'Normalization', 'probability' )

%% Mizuseki REM shift plots
%  Basically this means it seems that all of my tagged neurons have a
%  fraction of REM-shifters. 

sig_e13 = th_remphase_e13(:,1) < .05 & th_wakephase_e13(:,1) <.05;
sig_e13( isnan(th_remphase_e13(:,1)) | isnan(th_wakephase_e13(:,1)) ) = false;
remphase_e13_rad = th_remphase_e13(sig_e13,2); 
wakephase_e13_rad = th_wakephase_e13(sig_e13,2);
remphase_e13 = [ rad2deg(remphase_e13_rad) ; 360+rad2deg(remphase_e13_rad) ; rad2deg(remphase_e13_rad) ; 360+rad2deg(remphase_e13_rad) ];
wakephase_e13 = [ rad2deg(wakephase_e13_rad) ; 360+rad2deg(wakephase_e13_rad) ; 360+rad2deg(wakephase_e13_rad); rad2deg(wakephase_e13_rad) ];

% E13
figure
set(gcf,'Position', [1486 357 560 420])
subplot(2,2,3)
plot(remphase_e13, wakephase_e13,'.', 'Color', brown, 'MarkerSize', 10)
xlabel('REM theta phase')
ylabel('WAKE theta phase')
axis square
xticks([0 180 360 540 720])
yticks([0 180 360 540 720])
xlim([0 720])
ylim([0 720])

edges = 0:40:720;
subplot(2,2,1)

[vals, edges] = histcounts( [ rad2deg(remphase_e13_rad) ; 360+rad2deg(remphase_e13_rad) ], edges, 'Normalization', 'probability' );
plot(edges(1:end-1)+diff( edges )/2, vals )
xticks([0 180 360 540 720])
xlim([0 720])
ylim([0 0.15])
axis square

subplot(2,2,4)
[vals, edges] = histcounts( [ rad2deg(wakephase_e13_rad) ; 360+rad2deg(wakephase_e13_rad) ], edges, 'Normalization', 'probability' );
plot(vals, edges(1:end-1)+diff( edges )/2 )
yticks([0 180 360 540 720])
ylim([0 720])
axis xy
axis square

% E14

sig_e14 = th_remphase_e14(:,1) < .05 & th_wakephase_e14(:,1) <.05;
remphase_e14_rad = th_remphase_e14(sig_e14,2); 
wakephase_e14_rad = th_wakephase_e14(sig_e14,2);
remphase_e14 = [ rad2deg(remphase_e14_rad) ; 360+rad2deg(remphase_e14_rad) ; rad2deg(remphase_e14_rad) ; 360+rad2deg(remphase_e14_rad) ];
wakephase_e14 = [ rad2deg(wakephase_e14_rad) ; 360+rad2deg(wakephase_e14_rad) ; 360+rad2deg(wakephase_e14_rad); rad2deg(wakephase_e14_rad) ];

figure
%plot(remphase_e15, wakephase_e15, '.')
set(gcf,'Position', [1554 -158 560 420])
subplot(2,2,3)
plot(remphase_e14, wakephase_e14,'.r', 'MarkerSize', 10)
xlabel('REM theta phase')
ylabel('WAKE theta phase')
axis square
xticks([0 180 360 540 720])
yticks([0 180 360 540 720])
xlim([0 720])
ylim([0 720])

edges = 0:40:720;
subplot(2,2,1)

[vals, edges] = histcounts( [ rad2deg(remphase_e14_rad) ; 360+rad2deg(remphase_e14_rad) ], edges, 'Normalization', 'probability' );
plot(edges(1:end-1)+diff( edges )/2, vals )
xlim([0 720])
axis square

subplot(2,2,4)
[vals, edges] = histcounts( [ rad2deg(wakephase_e14_rad) ; 360+rad2deg(wakephase_e14_rad) ], edges, 'Normalization', 'probability' );
plot(vals, edges(1:end-1)+diff( edges )/2 )
ylim([0 720])
axis xy
axis square

% E15

sig_e15 = th_remphase_e15(:,1) < .05 & th_wakephase_e15(:,1) <.05;
remphase_e15_rad = th_remphase_e15(sig_e15,2); 
wakephase_e15_rad = th_wakephase_e15(sig_e15,2);
remphase_e15 = [ rad2deg(remphase_e15_rad) ; 360+rad2deg(remphase_e15_rad) ; rad2deg(remphase_e15_rad) ; 360+rad2deg(remphase_e15_rad) ];
wakephase_e15 = [ rad2deg(wakephase_e15_rad) ; 360+rad2deg(wakephase_e15_rad) ; 360+rad2deg(wakephase_e15_rad); rad2deg(wakephase_e15_rad) ];


figure
subplot(2,2,3)
set(gcf,'Position', [2067 373 560 420])
plot(remphase_e15, wakephase_e15, '.b')
xlabel('REM theta phase')
ylabel('WAKE theta phase')
axis square
xticks([0 180 360 540 720])
yticks([0 180 360 540 720])
xlim([0 720])
ylim([0 720])

edges = 0:40:720;
subplot(2,2,1)

[vals, edges] = histcounts( [ rad2deg(remphase_e15_rad) ; 360+rad2deg(remphase_e15_rad) ], edges, 'Normalization', 'probability' );
plot(edges(1:end-1)+diff( edges )/2, vals )
xlim([0 720])
axis square

subplot(2,2,4)
[vals, edges] = histcounts( [ rad2deg(wakephase_e15_rad) ; 360+rad2deg(wakephase_e15_rad) ], edges, 'Normalization', 'probability' );
plot(vals, edges(1:end-1)+diff( edges )/2 )
ylim([0 720])
axis xy
axis square

% E16

sig_e16 = th_remphase_e16(:,1) < .05 & th_wakephase_e16(:,1) <.05;
remphase_e16_rad = th_remphase_e16(sig_e16,2); 
wakephase_e16_rad = th_wakephase_e16(sig_e16,2);
remphase_e16 = [ rad2deg(remphase_e16_rad) ; 360+rad2deg(remphase_e16_rad) ; rad2deg(remphase_e16_rad) ; 360+rad2deg(remphase_e16_rad) ];
wakephase_e16 = [ rad2deg(wakephase_e16_rad) ; 360+rad2deg(wakephase_e16_rad) ; 360+rad2deg(wakephase_e16_rad); rad2deg(wakephase_e16_rad) ];


figure
subplot(2,2,3)
set(gcf,'Position', [2067 373 560 420])
plot(remphase_e16, wakephase_e16, '.k')
xlabel('REM theta phase')
ylabel('WAKE theta phase')
axis square
xticks([0 180 360 540 720])
yticks([0 180 360 540 720])
xlim([0 720])
ylim([0 720])

edges = 0:40:720;
subplot(2,2,1)

[vals, edges] = histcounts( [ rad2deg(remphase_e16_rad) ; 360+rad2deg(remphase_e16_rad) ], edges, 'Normalization', 'probability' );
plot(edges(1:end-1)+diff( edges )/2, vals )
xlim([0 720])
axis square

subplot(2,2,4)
[vals, edges] = histcounts( [ rad2deg(wakephase_e16_rad) ; 360+rad2deg(wakephase_e16_rad) ], edges, 'Normalization', 'probability' );
plot(vals, edges(1:end-1)+diff( edges )/2 )
ylim([0 720])
axis xy
axis square


%%
sig = 0.05;
% Pull out the theta phase locked cells significance values
sig_e13 = [ th_wakephase_e13(:,1) < sig th_remphase_e13(:,1) < sig ];
sig_e13( isnan(th_remphase_e13(:,1)) | isnan(th_wakephase_e13(:,1)), : ) = [];
phase_e13 = [ rad2deg( th_wakephase_e13(:,2) ) rad2deg( th_remphase_e13(:,2) )];
phase_e13( isnan(th_remphase_e13(:,1)) | isnan(th_wakephase_e13(:,1)), : ) = [];
phaserate_e13_wake = th_wake_e13; phaserate_e13_wake( isnan(th_remphase_e13(:,1)) | isnan(th_wakephase_e13(:,1)),: ) = [];
phaserate_e13_rem = th_rem_e13; phaserate_e13_rem( isnan(th_remphase_e13(:,1)) | isnan(th_wakephase_e13(:,1)),: ) = [];

sig_e14 = [ th_wakephase_e14(:,1) < sig th_remphase_e14(:,1) < sig ];
sig_e14( isnan(th_remphase_e14(:,1)) | isnan(th_wakephase_e14(:,1)), : ) = [];
phase_e14 = [ rad2deg( th_wakephase_e14(:,2) ) rad2deg( th_remphase_e14(:,2) )];
phase_e14( isnan(th_remphase_e14(:,1)) | isnan(th_wakephase_e14(:,1)), : ) = [];
phaserate_e14_wake = th_wake_e14; phaserate_e14_wake( isnan(th_remphase_e14(:,1)) | isnan(th_wakephase_e14(:,1)),: ) = [];
phaserate_e14_rem = th_rem_e14; phaserate_e14_rem( isnan(th_remphase_e14(:,1)) | isnan(th_wakephase_e14(:,1)),: ) = [];

sig_e15 = [ th_wakephase_e15(:,1) < sig th_remphase_e15(:,1) < sig];
sig_e15( isnan(th_remphase_e15(:,1)) | isnan(th_wakephase_e15(:,1)), : ) = [];
phase_e15 = [ rad2deg( th_wakephase_e15(:,2) ) rad2deg( th_remphase_e15(:,2) )];
phase_e15( isnan(th_remphase_e15(:,1)) | isnan(th_wakephase_e15(:,1)), : ) = [];
phaserate_e15_wake = th_wake_e15; phaserate_e15_wake( isnan(th_remphase_e15(:,1)) | isnan(th_wakephase_e15(:,1)),: ) = [];
phaserate_e15_rem = th_rem_e15; phaserate_e15_rem( isnan(th_remphase_e15(:,1)) | isnan(th_wakephase_e15(:,1)),: ) = [];

sig_e16 = [ th_wakephase_e16(:,1) < sig th_remphase_e16(:,1) < sig ];
sig_e16( isnan(th_remphase_e16(:,1)) | isnan(th_wakephase_e16(:,1)), : ) = [];
phase_e16 = [ rad2deg( th_wakephase_e16(:,2) ) rad2deg( th_remphase_e16(:,2) )];
phase_e16( isnan(th_remphase_e16(:,1)) | isnan(th_wakephase_e16(:,1)), : ) = [];
phaserate_e16_wake = th_wake_e16; phaserate_e16_wake( isnan(th_remphase_e16(:,1)) | isnan(th_wakephase_e16(:,1)),: ) = [];
phaserate_e16_rem = th_rem_e16; phaserate_e16_rem( isnan(th_remphase_e16(:,1)) | isnan(th_wakephase_e16(:,1)),: ) = [];


prop_wake = [sum(sig_e13(:,1)) ./ length( sig_e13 ) sum(sig_e14(:,1)) ./ length( sig_e14 ) ...
             sum(sig_e15(:,1)) ./ length( sig_e15 ) sum(sig_e16(:,1)) ./ length( sig_e16 )];
         
prop_rem = [sum(sig_e13(:,2)) ./ length( sig_e13 ) sum(sig_e14(:,2)) ./ length( sig_e14 ) ...
             sum(sig_e15(:,2)) ./ length( sig_e15 ) sum(sig_e16(:,2)) ./ length( sig_e16 )];
         
prop_both = [sum(all(sig_e13,2)) ./ length( sig_e13 ) sum(all(sig_e14,2)) ./ length( sig_e14 ) ...
             sum(all(sig_e15,2)) ./ length( sig_e15 ) sum(all(sig_e16,2)) ./ length( sig_e16 ) ];

%%
phasebins = [-2.9845   -2.6704   -2.3562   -2.0420   -1.7279   -1.4137   -1.0996   -0.7854   -0.4712   -0.1571    0.1571    0.4712    0.7854   1.0996    1.4137    1.7279    2.0420    2.3562    2.6704    2.9845];

phasebins_2cycles = [ rad2deg( phasebins ) 360+rad2deg(phasebins) ];
imagesc( phasebins_2cycles, 1:size( phaserate_e16_wake,1 ), phaserate_e16_wake)


%% REM shifting neurons
% Consider group of neurons that are phaselocked in both REM and WAKE -
% the subset that locks to the trough in WAKE and shifts to the peak in REM
% is the REM shifting population

figure
set(gcf,'Position', [55 264 2144 434])

% E13
phase_e13_sig = phase_e13( all(sig_e13,2),: );
rate_e13_wake = phaserate_e13_wake(all(sig_e13,2),:) ;
rate_e13_rem = phaserate_e13_rem(all(sig_e13,2),:) ;
wake_trough_e13 = phase_e13_sig(:,1) >= 120 & phase_e13_sig(:,1) <= 300;
rem_peak_e13 = ( phase_e13_sig(:,2) >= 0 & phase_e13_sig(:,2) < 120 ) |  ( phase_e13_sig(:,2) > 300 & phase_e13_sig(:,2) < 360 );
remshift_e13 = sum( wake_trough_e13 & rem_peak_e13 ) ./ length(phase_e13_sig);%sum(wake_trough_e13);

subplot(2,8,1)
imagesc(phasebins_2cycles, 1:sum(wake_trough_e13 & rem_peak_e13),imgaussfilt( rate_e13_wake(wake_trough_e13 & rem_peak_e13,:) ) )
colormap(redblue)
xline([120 300], 1)
title('E13')
ylabel('WAKE')
subplot(2,8,9)
imagesc(phasebins_2cycles, 1:sum(wake_trough_e13 & rem_peak_e13), imgaussfilt( rate_e13_rem(wake_trough_e13 & rem_peak_e13,:) ) )
colormap(redblue)
xline([120 300], 1)
ylabel('REM')

% E14
phase_e14_sig = phase_e14( all(sig_e14,2),: );
rate_e14_wake = phaserate_e14_wake(all(sig_e14,2),:) ;
rate_e14_rem = phaserate_e14_rem(all(sig_e14,2),:) ;
wake_trough_e14 = phase_e14_sig(:,1) >= 120 & phase_e14_sig(:,1) <= 300;
rem_peak_e14 = ( phase_e14_sig(:,2) >= 0 & phase_e14_sig(:,2) < 120 ) |  ( phase_e14_sig(:,2) > 300 & phase_e14_sig(:,2) < 360 );
remshift_e14 = sum( wake_trough_e14 & rem_peak_e14 ) ./ length(phase_e14_sig);%sum(wake_trough_e14);

subplot(2,8,2)
imagesc(phasebins_2cycles, 1:sum(wake_trough_e14 & rem_peak_e14),imgaussfilt( rate_e14_wake(wake_trough_e14 & rem_peak_e14,:) ) )
colormap(redblue)
xline([120 300], 1)
title('E14')
subplot(2,8,10)
imagesc(phasebins_2cycles, 1:sum(wake_trough_e14 & rem_peak_e14), imgaussfilt( rate_e14_rem(wake_trough_e14 & rem_peak_e14,:) ) )
colormap(redblue)
xline([120 300], 1)

% E15
phase_e15_sig = phase_e15( all(sig_e15,2),: );
rate_e15_wake = phaserate_e15_wake(all(sig_e15,2),:) ;
rate_e15_rem = phaserate_e15_rem(all(sig_e15,2),:) ;
wake_trough_e15 = phase_e15_sig(:,1) >= 120 & phase_e15_sig(:,1) <= 300;
rem_peak_e15 = ( phase_e15_sig(:,2) >= 0 & phase_e15_sig(:,2) < 120 ) |  ( phase_e15_sig(:,2) > 300 & phase_e15_sig(:,2) < 360 );
remshift_e15 = sum( wake_trough_e15 & rem_peak_e15 ) ./ length(phase_e15_sig);%sum(wake_trough_e15);

subplot(2,8,3)
imagesc(phasebins_2cycles, 1:sum(wake_trough_e15 & rem_peak_e15),imgaussfilt( rate_e15_wake(wake_trough_e15 & rem_peak_e15,:)  ) )
colormap(redblue)
xline([120 300], 1)
title('E15')
subplot(2,8,11)
imagesc(phasebins_2cycles, 1:sum(wake_trough_e15 & rem_peak_e15), imgaussfilt( rate_e15_rem(wake_trough_e15 & rem_peak_e15,:) ) )
colormap(redblue)
xline([120 300], 1)

% E16
phase_e16_sig = phase_e16( all(sig_e16,2),: );
rate_e16_wake = phaserate_e16_wake(all(sig_e16,2),:) ;
rate_e16_rem = phaserate_e16_rem(all(sig_e16,2),:) ;
wake_trough_e16 = phase_e16_sig(:,1) >= 120 & phase_e16_sig(:,1) <= 300;
rem_peak_e16 = ( phase_e16_sig(:,2) >= 0 & phase_e16_sig(:,2) < 120 ) |  ( phase_e16_sig(:,2) > 300 & phase_e16_sig(:,2) < 360 );
remshift_e16 = sum( wake_trough_e16 & rem_peak_e16 ) ./ length(phase_e16_sig);%sum(wake_trough_e16);

subplot(2,8,4)
imagesc(phasebins_2cycles, 1:sum(wake_trough_e16 & rem_peak_e16),imgaussfilt( rate_e16_wake(wake_trough_e16 & rem_peak_e16,:)  ) )
colormap(redblue)
xline([120 300], 1)
title('E16')
subplot(2,8,12)
imagesc(phasebins_2cycles, 1:sum(wake_trough_e16 & rem_peak_e16), imgaussfilt( rate_e16_rem(wake_trough_e16 & rem_peak_e16,:) ) )
colormap(redblue)
xline([120 300], 1)


subplot(2,8,[5 6 13 14])
bar([remshift_e13 remshift_e14 remshift_e15 remshift_e16])
ylabel('Fraction of REM-shifters (trough in WAKE, peak in REM)', 'fontsize', 14)
xticks([1 2 3 4])
xticklabels([{'E13', 'E14', 'E15', 'E16'}])
title('Consider subgroup of PYRs that are phase locked in REM and WAKE')


%% Get a bootstrapped estimate of the error

% E13
ns = 100;
remshift_e13_dist = []; remshift_e14_dist = []; remshift_e15_dist = []; remshift_e16_dist = [];
for kp = 1:ns
phase_e13_sig = phase_e13( all(sig_e13,2),: );
phase_e13_sig = phase_e13_sig( randi(length( phase_e13_sig ), 1,length( phase_e13_sig )), : );
wake_trough_e13 = phase_e13_sig(:,1) >= 120 & phase_e13_sig(:,1) <= 300;
remshift_e13_dist(kp) = sum( wake_trough_e13 & rem_peak_e13 ) ./ length(phase_e13_sig);

phase_e14_sig = phase_e14( all(sig_e14,2),: );
phase_e14_sig = phase_e14_sig( randi(length( phase_e14_sig ), 1,length( phase_e14_sig )), : );
wake_trough_e14 = phase_e14_sig(:,1) >= 120 & phase_e14_sig(:,1) <= 300;
remshift_e14_dist(kp) = sum( wake_trough_e14 & rem_peak_e14 ) ./ length(phase_e14_sig);

phase_e15_sig = phase_e15( all(sig_e15,2),: );
phase_e15_sig = phase_e15_sig( randi(length( phase_e15_sig ), 1,length( phase_e15_sig )), : );
wake_trough_e15 = phase_e15_sig(:,1) >= 120 & phase_e15_sig(:,1) <= 300;
remshift_e15_dist(kp) = sum( wake_trough_e15 & rem_peak_e15 ) ./ length(phase_e15_sig);

phase_e16_sig = phase_e16( all(sig_e16,2),: );
phase_e16_sig = phase_e16_sig( randi(length( phase_e16_sig ), 1,length( phase_e16_sig )), : );
wake_trough_e16 = phase_e16_sig(:,1) >= 120 & phase_e16_sig(:,1) <= 300;
remshift_e16_dist(kp) = sum( wake_trough_e16 & rem_peak_e16 ) ./ length(phase_e16_sig);
end

data = [remshift_e13 remshift_e14 remshift_e15 remshift_e16]';
err = [std(remshift_e13_dist') std(remshift_e14_dist') std(remshift_e15_dist') std(remshift_e16_dist')];

bar(1:4,data)
hold on
er = errorbar(1:4,data,err, err);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';  

hold off
%%

load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/analyze/cellphys_v7.mat')

nbins = 40;     % total number of phase bins

th_rem_e13 = []; th_wake_e13 = [];
th_remphase_e13 = []; th_wakephase_e13 = [];
th_rem_e14 = []; th_wake_e14 = [];
th_remphase_e14 = []; th_wakephase_e14 = [];
th_rem_e15 = []; th_wake_e15 = [];
th_remphase_e15 = []; th_wakephase_e15 = [];
th_rem_e16 = []; th_wake_e16 = [];
th_remphase_e16 = []; th_wakephase_e16 = [];

% E13
th_powerphasefr_rem_sessions = pyr_opto_phys(1).th_powerphasefr_rem;
th_powerphasefr_wake_sessions = pyr_opto_phys(1).th_powerphasefr_wake;
th_lock_rem_sessions = pyr_opto_phys(1).th_phaselock_rem;
th_lock_wake_sessions = pyr_opto_phys(1).th_phaselock_wake;
% Loop over sessions, each of has a struct of theta phase data for each cells 
for kp = 1:length(th_powerphasefr_rem_sessions)
    
    % average REM theta ratemaps
    if ~isempty( th_powerphasefr_rem_sessions(kp).ratemap )
        ave_fr = cell2mat( cellfun(@nanmean, th_powerphasefr_rem_sessions(kp).ratemap, 'UniformOutput', false)' );
        th_rem_e13 = [th_rem_e13 ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        th_remphase_e13 = [th_remphase_e13 ; [ th_lock_rem_sessions(kp).phasestats.p' th_lock_rem_sessions(kp).phasestats.m']];
    else
        th_rem_e13 = [ th_rem_e13 ; nan( length( th_powerphasefr_rem_sessions(kp).UID), nbins ) ];
        th_remphase_e13 = [th_remphase_e13  ; nan( length( th_powerphasefr_rem_sessions(kp).UID), 2 ) ];
    end
    % average WAKE theta ratemaps
    if ~isempty( th_powerphasefr_wake_sessions(kp).ratemap )
        ave_fr = cell2mat( cellfun(@nanmean, th_powerphasefr_wake_sessions(kp).ratemap, 'UniformOutput', false)' );
        th_wake_e13 = [th_wake_e13 ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        th_wakephase_e13 = [th_wakephase_e13 ; [ th_lock_wake_sessions(kp).phasestats.p' th_lock_wake_sessions(kp).phasestats.m']];
    else
        th_wake_e13 = [ th_wake_e13 ; nan( length( th_powerphasefr_wake_sessions(kp).UID), nbins ) ];
        th_wakephase_e13 = [th_wakephase_e13  ; nan( length( th_powerphasefr_wake_sessions(kp).UID), 2 ) ];
    end   
    
end


% 14
th_powerphasefr_rem_sessions = pyr_opto_phys(2).th_powerphasefr_rem;
th_powerphasefr_wake_sessions = pyr_opto_phys(2).th_powerphasefr_wake;
th_lock_rem_sessions = pyr_opto_phys(2).th_phaselock_rem;
th_lock_wake_sessions = pyr_opto_phys(2).th_phaselock_wake;
% Loop over sessions, each of has a struct of theta phase data for each cells 
for kp = 1:length(th_powerphasefr_rem_sessions)
    
    % average REM theta ratemaps
    if ~isempty( th_powerphasefr_rem_sessions(kp).ratemap )
        ave_fr = cell2mat( cellfun(@nanmean, th_powerphasefr_rem_sessions(kp).ratemap, 'UniformOutput', false)' );
        th_rem_e14 = [th_rem_e14 ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        th_remphase_e14 = [th_remphase_e14 ; [ th_lock_rem_sessions(kp).phasestats.p' th_lock_rem_sessions(kp).phasestats.m']];
    else
        th_rem_e14 = [ th_rem_e14 ; nan( length( th_powerphasefr_rem_sessions(kp).UID), nbins ) ];
        th_remphase_e14 = [th_remphase_e14  ; nan( length( th_powerphasefr_rem_sessions(kp).UID), 2 ) ];
    end
    % average WAKE theta ratemaps
    if ~isempty( th_powerphasefr_wake_sessions(kp).ratemap )
        ave_fr = cell2mat( cellfun(@nanmean, th_powerphasefr_wake_sessions(kp).ratemap, 'UniformOutput', false)' );
        th_wake_e14 = [th_wake_e14 ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        th_wakephase_e14 = [th_wakephase_e14 ; [ th_lock_wake_sessions(kp).phasestats.p' th_lock_wake_sessions(kp).phasestats.m']];
    else
        th_wake_e14 = [ th_wake_e14 ; nan( length( th_powerphasefr_wake_sessions(kp).UID), nbins ) ];
        th_wakephase_e14 = [th_wakephase_e14  ; nan( length( th_powerphasefr_wake_sessions(kp).UID), 2 ) ];
    end   
    
end

% 15
th_powerphasefr_rem_sessions = pyr_opto_phys(3).th_powerphasefr_rem;
th_powerphasefr_wake_sessions = pyr_opto_phys(3).th_powerphasefr_wake;
th_lock_rem_sessions = pyr_opto_phys(3).th_phaselock_rem;
th_lock_wake_sessions = pyr_opto_phys(3).th_phaselock_wake;
% Loop over sessions, each of has a struct of theta phase data for each cells 
for kp = 1:length(th_powerphasefr_rem_sessions)
    
    % average REM theta ratemaps
    if ~isempty( th_powerphasefr_rem_sessions(kp).ratemap )
        ave_fr = cell2mat( cellfun(@nanmean, th_powerphasefr_rem_sessions(kp).ratemap, 'UniformOutput', false)' );
        th_rem_e15 = [th_rem_e15 ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        th_remphase_e15 = [th_remphase_e15 ; [ th_lock_rem_sessions(kp).phasestats.p' th_lock_rem_sessions(kp).phasestats.m']];
    else
        th_rem_e15 = [ th_rem_e15 ; nan( length( th_powerphasefr_rem_sessions(kp).UID), nbins ) ];
        th_remphase_e15 = [th_remphase_e15  ; nan( length( th_powerphasefr_rem_sessions(kp).UID), 2 ) ];
    end
    % average WAKE theta ratemaps
    if ~isempty( th_powerphasefr_wake_sessions(kp).ratemap )
        ave_fr = cell2mat( cellfun(@nanmean, th_powerphasefr_wake_sessions(kp).ratemap, 'UniformOutput', false)' );
        th_wake_e15 = [th_wake_e15 ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        th_wakephase_e15 = [th_wakephase_e15 ; [ th_lock_wake_sessions(kp).phasestats.p' th_lock_wake_sessions(kp).phasestats.m']];
    else
        th_wake_e15 = [ th_wake_e15 ; nan( length( th_powerphasefr_wake_sessions(kp).UID), nbins ) ];
        th_wakephase_e15 = [th_wakephase_e15  ; nan( length( th_powerphasefr_wake_sessions(kp).UID), 2 ) ];
    end   
    
end

% 16
th_powerphasefr_rem_sessions = pyr_opto_phys(4).th_powerphasefr_rem;
th_powerphasefr_wake_sessions = pyr_opto_phys(4).th_powerphasefr_wake;
th_lock_rem_sessions = pyr_opto_phys(4).th_phaselock_rem;
th_lock_wake_sessions = pyr_opto_phys(4).th_phaselock_wake;
% Loop over sessions, each of has a struct of theta phase data for each cells 
for kp = 1:length(th_powerphasefr_rem_sessions)
    
    % average REM theta ratemaps
    if ~isempty( th_powerphasefr_rem_sessions(kp).ratemap )
        ave_fr = cell2mat( cellfun(@nanmean, th_powerphasefr_rem_sessions(kp).ratemap, 'UniformOutput', false)' );
        th_rem_e16 = [th_rem_e16 ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        th_remphase_e16 = [th_remphase_e16 ; [ th_lock_rem_sessions(kp).phasestats.p' th_lock_rem_sessions(kp).phasestats.m']];
    else
        th_rem_e16 = [ th_rem_e16 ; nan( length( th_powerphasefr_rem_sessions(kp).UID), nbins ) ];
        th_remphase_e16 = [th_remphase_e16  ; nan( length( th_powerphasefr_rem_sessions(kp).UID), 2 ) ];
    end
    % average WAKE theta ratemaps
    if ~isempty( th_powerphasefr_wake_sessions(kp).ratemap )
        ave_fr = cell2mat( cellfun(@nanmean, th_powerphasefr_wake_sessions(kp).ratemap, 'UniformOutput', false)' );
        th_wake_e16 = [th_wake_e16 ; [ave_fr ave_fr] ./ max(ave_fr,[],2) ];
        th_wakephase_e16 = [th_wakephase_e16 ; [ th_lock_wake_sessions(kp).phasestats.p' th_lock_wake_sessions(kp).phasestats.m']];
    else
        th_wake_e16 = [ th_wake_e16 ; nan( length( th_powerphasefr_wake_sessions(kp).UID), nbins ) ];
        th_wakephase_e16 = [th_wakephase_e16  ; nan( length( th_powerphasefr_wake_sessions(kp).UID), 2 ) ];
    end   
    
end
