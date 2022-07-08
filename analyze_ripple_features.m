
%% Collect ripple features based on wavelet decomposition
%  All this is stored in the .ripple_spectral.mat (computed by get_ripple_features.m), 
%  so it's just a matter  of pulling it out

clear

isLocal = true;

nfreqs = 100;
frange = [60 250];
path_to_basepaths = '/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat';
load(path_to_basepaths)
 
cellpops = [{'basepaths_e13'}, {'basepaths_e14'}, {'basepaths_e15'}, {'basepaths_e16'} ];

ripple_specdist( 1:length(cellpops) ) = struct('powers_pulse_any', [], 'freqs_pulse_any', [], 'powers_pulse_any_uncorr', [], 'freqs_pulse_any_uncorr', [], 'powers_pulse_max1', [],...
                                               'freqs_pulse_max1', [],'powers_pulse_max1_uncorr', [],'freqs_pulse_max1_uncorr', [], 'powers_pulse_max2', [], 'freqs_pulse_max2', [],...
                                               'powers_pulse_max2_uncorr', [], 'freqs_pulse_max2_uncorr', [], 'durs_any', [], 'durs_max1', [], 'durs_max2', []);

bp_counter = 1;
for cellpop = cellpops
    basepaths = eval(cellpop{1});
    fprintf('Processing %s ...\n', cellpop{1})
    for ii = 1:length(basepaths)
        
        fprintf('%d/%d \n', ii, length(basepaths))
        % Get basepath
        basepath = alterPath( basepaths{ii}, isLocal );
        basename = bz_BasenameFromBasepath(basepath);
        load( fullfile(basepath, [basename '.ripple_spectral.mat']) )
        
        % PS within ripples to which tagged neurons contribute (ANY)
        ps = ripple_spectral.ripple_ps(:,ripple_spectral.ripp_any);
        powers_pulse_any = [];
        freqs_pulse_any = [];
        powers_pulse_any_uncorr = [];
        freqs_pulse_any_uncorr = [];
        for kp = 1:size(ps,2)
            p = ps(:,kp);
            plin = interp1([ ripple_spectral.freqs(1) ripple_spectral.freqs(end) ],[p(1) p(end)],ripple_spectral.freqs); 
            if any(isnan(p))
                continue
            end
            if any( p - plin' > 0 )
                % Linearly correct for 1/f
                [~, arg] = max( p - plin');
                powers_pulse_any = [ powers_pulse_any p(arg) ];
                freqs_pulse_any = [ freqs_pulse_any ripple_spectral.freqs(arg) ];
                % Uncorrected
                [~, arg] = max( p );
                powers_pulse_any_uncorr = [ powers_pulse_any_uncorr p(arg) ];
                freqs_pulse_any_uncorr = [ freqs_pulse_any_uncorr ripple_spectral.freqs(arg) ];
            end
        end
        
        % PS within ripples to which tagged neurons contribute (MAX1)
        ps = ripple_spectral.ripple_ps(:,ripple_spectral.ripp_max1);
        powers_pulse_max1 = [];
        freqs_pulse_max1 = [];
        powers_pulse_max1_uncorr = [];
        freqs_pulse_max1_uncorr = [];
        for kp = 1:size(ps,2)
            p = ps(:,kp);
            plin = interp1([ ripple_spectral.freqs(1) ripple_spectral.freqs(end) ],[p(1) p(end)],ripple_spectral.freqs); 
            if any(isnan(p))
                continue
            end
            if any( p - plin' > 0 )
                % Linearly correct for 1/f
                [~, arg] = max( p - plin');
                powers_pulse_max1 = [ powers_pulse_max1 p(arg) ];
                freqs_pulse_max1 = [ freqs_pulse_max1 ripple_spectral.freqs(arg) ];
                % Uncorrected
                [~, arg] = max( p );
                powers_pulse_max1_uncorr = [ powers_pulse_max1_uncorr p(arg) ];
                freqs_pulse_max1_uncorr = [ freqs_pulse_max1_uncorr ripple_spectral.freqs(arg) ];
            end
        end
        
        % PS within ripples to which tagged neurons contribute (MAX2)
        ps = ripple_spectral.ripple_ps(:,ripple_spectral.ripp_max2);
        powers_pulse_max2 = [];
        freqs_pulse_max2 = [];
        powers_pulse_max2_uncorr = [];
        freqs_pulse_max2_uncorr = [];
        for kp = 1:size(ps,2)
            p = ps(:,kp);
            plin = interp1([ ripple_spectral.freqs(1) ripple_spectral.freqs(end) ],[p(1) p(end)],ripple_spectral.freqs); 
            if any(isnan(p))
                continue
            end
            if any( p - plin' > 0 )
                % Linearly correct for 1/f
                [~, arg] = max( p - plin');
                powers_pulse_max2 = [ powers_pulse_max2 p(arg) ];
                freqs_pulse_max2 = [ freqs_pulse_max2 ripple_spectral.freqs(arg) ];
                % Uncorrected
                [~, arg] = max( p );
                powers_pulse_max2_uncorr = [ powers_pulse_max2_uncorr p(arg) ];
                freqs_pulse_max2_uncorr = [ freqs_pulse_max2_uncorr ripple_spectral.freqs(arg) ];
            end
        end
        
        
        ripple_specdist(bp_counter).powers_pulse_any = [ ripple_specdist(bp_counter).powers_pulse_any powers_pulse_any ];
        ripple_specdist(bp_counter).freqs_pulse_any = [ ripple_specdist(bp_counter).freqs_pulse_any freqs_pulse_any ];
        ripple_specdist(bp_counter).powers_pulse_any_uncorr = [ ripple_specdist(bp_counter).powers_pulse_any_uncorr powers_pulse_any_uncorr ];
        ripple_specdist(bp_counter).freqs_pulse_any_uncorr = [ ripple_specdist(bp_counter).freqs_pulse_any_uncorr freqs_pulse_any_uncorr ];
        
        ripple_specdist(bp_counter).powers_pulse_max1 = [ ripple_specdist(bp_counter).powers_pulse_max1 powers_pulse_max1 ];
        ripple_specdist(bp_counter).freqs_pulse_max1 = [ ripple_specdist(bp_counter).freqs_pulse_max1 freqs_pulse_max1 ];
        ripple_specdist(bp_counter).powers_pulse_max1_uncorr = [ ripple_specdist(bp_counter).powers_pulse_max1_uncorr powers_pulse_max1_uncorr ];
        ripple_specdist(bp_counter).freqs_pulse_max1_uncorr = [ ripple_specdist(bp_counter).freqs_pulse_max1_uncorr freqs_pulse_max1_uncorr ];
        
        ripple_specdist(bp_counter).powers_pulse_max2 = [ ripple_specdist(bp_counter).powers_pulse_max2 powers_pulse_max2 ];
        ripple_specdist(bp_counter).freqs_pulse_max2 = [ ripple_specdist(bp_counter).freqs_pulse_max2 freqs_pulse_max2 ];
        ripple_specdist(bp_counter).powers_pulse_max2_uncorr = [ ripple_specdist(bp_counter).powers_pulse_max2_uncorr powers_pulse_max2_uncorr ];
        ripple_specdist(bp_counter).freqs_pulse_max2_uncorr = [ ripple_specdist(bp_counter).freqs_pulse_max2_uncorr freqs_pulse_max2_uncorr ];
        
        ripple_specdist(bp_counter).durs_any = [ ripple_specdist(bp_counter).durs_any ; ripple_spectral.ripple_durs( ripple_spectral.ripp_any ) ];
        ripple_specdist(bp_counter).durs_max1 = [ ripple_specdist(bp_counter).durs_max1 ;  ripple_spectral.ripple_durs( ripple_spectral.ripp_max1 ) ];
        ripple_specdist(bp_counter).durs_max2 = [ ripple_specdist(bp_counter).durs_max2 ; ripple_spectral.ripple_durs( ripple_spectral.ripp_max2 ) ];
    end
    bp_counter = bp_counter+1;
end

%% Collect ripple features based on hilbert transform
%  All this is stored in the ripples.events.mat file, so it's just a matter
%  of pulling it out

clear

isLocal = true;
nfreqs = 100;
frange = [60 250];
path_to_basepaths = '/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat';
load(path_to_basepaths)
 
cellpops = [{'basepaths_e13'}, {'basepaths_e14'}, {'basepaths_e15'}, {'basepaths_e16'} ];

ripple_specdist( 1:length(cellpops) ) = struct('powers_pulse_any', [], 'freqs_pulse_any', [], 'powers_pulse_max1', [],...
                                               'freqs_pulse_max1', [], 'durs_any', [], 'durs_max1', []);
                                               
bp_counter = 1;
for cellpop = cellpops
    basepaths = eval(cellpop{1});
    fprintf('Processing %s ...\n', cellpop{1})
    for ii = 1:length(basepaths)
        
        fprintf('%d/%d \n', ii, length(basepaths))
        % Get basepath
        basepath = alterPath( basepaths{ii}, isLocal );
        basename = bz_BasenameFromBasepath(basepath);
        load( fullfile(basepath, [basename '.ripple_spectral.mat']) )
        
        % Get Kilosort folder
        load( fullfile(basepath, [basename '.ripples.events.mat']) )                                % ripples
        load( fullfile(basepath, [basename '.pulses.events.mat']) )
        load( fullfile(basepath, [ basename '.cell_metrics.cellinfo.mat']) )
        load( fullfile(basepath, [ basename '.spikes.cellinfo.mat']) )
        
        % Optotagged pyramidal cells
        opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
        pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
        pyr_opto_id = opto_id & pyr_id;
        
        pulses_l_ints = pulses.timestamps( pulses.duration > 0.005, : ); 
        ripple_indicator = ~( InIntervals(ripples.timestamps(:,1), pulses_l_ints) |...   % Spontaneous ripples (exclude opto)
                                          InIntervals(ripples.timestamps(:,2), pulses_l_ints) );
                                      
        % Get indicators for each tagged cell
        log_arr = false( size(ripples.timestamps,1), sum(pyr_opto_id) ); counter = 1;
        for kp = find( pyr_opto_id )
            [status, intt ] = InIntervals( spikes.times{kp}, ripples.timestamps);                               
            rip_p = unique( intt ); rip_p(rip_p == 0) = [];
            log_arr(rip_p,counter) = true;
            counter = counter+1;
        end
        % Ripples in which at least one tagged cell spiked 
        ripp_any = logical( sum(log_arr, 2) > 0 ); ripp_any = ripp_any & ripple_indicator;
        % Ripples in which at one tagged cell spiked
        ripp_max1 = logical( sum(log_arr, 2) > 0 & sum(log_arr, 2) < 2 ); ripp_max1 = ripp_max1 & ripple_indicator;
        
        ripple_specdist(bp_counter).powers_pulse_any = [ripple_specdist(bp_counter).powers_pulse_any ; ripples.data.peakAmplitude(ripp_any) ];
        ripple_specdist(bp_counter).freqs_pulse_any = [ripple_specdist(bp_counter).freqs_pulse_any ; ripples.data.peakFrequency(ripp_any) ];
        ripple_specdist(bp_counter).powers_pulse_max1 = [ripple_specdist(bp_counter).powers_pulse_max1 ; ripples.data.peakAmplitude(ripp_max1) ];
        ripple_specdist(bp_counter).freqs_pulse_max1 = [ripple_specdist(bp_counter).freqs_pulse_max1 ; ripples.data.peakFrequency(ripp_max1) ];
        ripple_specdist(bp_counter).durs_any = [ripple_specdist(bp_counter).durs_any ; ripples.data.duration(ripp_any) ];
        ripple_specdist(bp_counter).durs_max1 = [ripple_specdist(bp_counter).durs_max1 ; ripples.data.duration(ripp_max1) ];
        
    end
    bp_counter = bp_counter+1;
end


%% Pull out (average) ripple-triggered wavelet spectrograms 
%  in example sessions across birthdates

base_e13 = alterPath( '/Volumes/Data/e13/e13_1m1/e13_1m1_200224', true );
base_e14 = alterPath( '/Volumes/RomanExternal/e14/e14_1m1/e14_1m1_200201', true );
base_e15 = alterPath( '/Volumes/Data1/e15/e15_5m1/e15_5m1_200313_1', true );
base_e16 = alterPath( '/Volumes/Data2/e16/e16_1f1/e16_1f1_200918', true );

% E13

figure
set(gcf,'Position', [175         527        1227         271])

subplot(1,4,1)
basename = bz_BasenameFromBasepath(base_e13);
load(fullfile(base_e13, [ basename  '.ripple_spectral.mat']))

imagesc(ripple_spectral.wavespec_tvec, log2(ripple_spectral.freqs), log10(ripple_spectral.wavespec_mean_max1'))
LogScale('y',2)
axis xy
caxis([1.7753    3.6])
xlim([-0.05 0.1])

ylabel('Frequency', 'fontsize', 14)
xlabel('Time (s)', 'fontsize', 14)
colorbar

% E14

subplot(1,4,2)
basename = bz_BasenameFromBasepath(base_e14);
load(fullfile(base_e14, [ basename  '.ripple_spectral.mat']))

imagesc(ripple_spectral.wavespec_tvec, log2(ripple_spectral.freqs), log10(ripple_spectral.wavespec_mean_max1'))
LogScale('y',2)
axis xy
caxis([1.7753    3.62])
xlim([-0.05 0.1])
xlabel('Time (s)', 'fontsize', 14)
colorbar

% E14

subplot(1,4,3)
basename = bz_BasenameFromBasepath(base_e15);
load(fullfile(base_e15, [ basename  '.ripple_spectral.mat']))

imagesc(ripple_spectral.wavespec_tvec, log2(ripple_spectral.freqs), log10(ripple_spectral.wavespec_mean_max1'))
LogScale('y',2)
axis xy
caxis([1.7753    3.62])
xlim([-0.05 0.1])
xlabel('Time (s)', 'fontsize', 14)
colorbar

subplot(1,4,4)
basename = bz_BasenameFromBasepath(base_e16);
load(fullfile(base_e16, [ basename  '.ripple_spectral.mat']))

imagesc(ripple_spectral.wavespec_tvec, log2(ripple_spectral.freqs), log10(ripple_spectral.wavespec_mean_max1'))
LogScale('y',2)
axis xy
caxis([1.7753    3.62])
xlim([-0.05 0.1])
xlabel('Time (s)', 'fontsize', 14)
colorbar

%% Pull out the individual-event / cross-event mean power spectra associated 
%  with the example sessions above

brown = [133 87 34]/254;
nlines = 50;

% E13

basename = bz_BasenameFromBasepath(base_e13);
load(fullfile(base_e13, [ basename  '.ripple_spectral.mat']))

ps = ripple_spectral.ripple_ps(:,ripple_spectral.ripp_max1);
ps_inds = randi(sum(ripple_spectral.ripp_max1), 1, nlines);

figure; hold on
gg = cell(1,nlines);
for kp = 1:nlines
    gg{kp} = plot( log2(ripple_spectral.freqs), ps(:,ps_inds(kp)), 'Color', brown);
    gg{kp}.Color(4) = 0.2;
end
plot( log2(ripple_spectral.freqs), mean( ps, 2), 'Color', brown, 'linewidth', 3)

% E14

basename = bz_BasenameFromBasepath(base_e14);
load(fullfile(base_e14, [ basename  '.ripple_spectral.mat']))

ps = ripple_spectral.ripple_ps(:,ripple_spectral.ripp_max1);
ps_inds = randi(sum(ripple_spectral.ripp_max1), 1, nlines);

gg = cell(1,nlines);
for kp = 1:nlines
    gg{kp} = plot( log2(ripple_spectral.freqs), ps(:,ps_inds(kp)), 'r' );
    gg{kp}.Color(4) = 0.2;
end
plot( log2(ripple_spectral.freqs), mean( ps, 2), '-r', 'linewidth', 3)

% E15

basename = bz_BasenameFromBasepath(base_e15);
load(fullfile(base_e15, [ basename  '.ripple_spectral.mat']))

ps = ripple_spectral.ripple_ps(:,ripple_spectral.ripp_max1);
ps_inds = randi(sum(ripple_spectral.ripp_max1), 1, nlines);

gg = cell(1,nlines);
for kp = 1:nlines
    gg{kp} = plot( log2(ripple_spectral.freqs), ps(:,ps_inds(kp)), 'b' );
    gg{kp}.Color(4) = 0.2;
end
plot( log2(ripple_spectral.freqs), mean( ps, 2), '-b', 'linewidth', 3)
xlim([min( log2(ripple_spectral.freqs) ) max( log2(ripple_spectral.freqs) )])
LogScale('x', 2)

% E16

basename = bz_BasenameFromBasepath(base_e16);
load(fullfile(base_e16, [ basename  '.ripple_spectral.mat']))

ps = ripple_spectral.ripple_ps(:,ripple_spectral.ripp_max1);
ps_inds = randi(sum(ripple_spectral.ripp_max1), 1, nlines);

gg = cell(1,nlines);
for kp = 1:nlines
    gg{kp} = plot( log2(ripple_spectral.freqs), ps(:,ps_inds(kp)), 'k' );
    gg{kp}.Color(4) = 0.2;
end
plot( log2(ripple_spectral.freqs), mean( ps, 2), '-k', 'linewidth', 3)
xlim([min( log2(ripple_spectral.freqs) ) max( log2(ripple_spectral.freqs) )])
LogScale('x', 2)

xlabel('Frequency','fontsize',20)
ylabel('log(amplitude)','fontsize',20)
title('Ripple power (example sessions)', 'fontsize', 20)

%% AMPLITUDES and FREQUENCIES (MAX1, corrected)
%  peak amplitude within ripple events

%%% AMPLITUDES

figure
subplot(2,2,1)
hold on
e13_amp =  log10( ripple_specdist(1).powers_pulse_max1 );
e14_amp = log10(  ripple_specdist(2).powers_pulse_max1 );
e15_amp =  log10( ripple_specdist(3).powers_pulse_max1 );
e16_amp =  log10( ripple_specdist(4).powers_pulse_max1 );

[cnts_e13, edges_e13] = histcounts( e13_amp, 14, 'Normalization', 'probability' );
[cnts_e14, edges_e14] = histcounts( e14_amp,14, 'Normalization', 'probability' );
[cnts_e15, edges_e15] = histcounts( e15_amp,20, 'Normalization', 'probability' );
[cnts_e16, edges_e16] = histcounts( e16_amp,17, 'Normalization', 'probability' );


plot(edges_e13(1:end-1)+median( diff(edges_e13) )/2, cnts_e13, 'Color',brown, 'linewidth',2)
plot(edges_e14(1:end-1)+median( diff(edges_e14) )/2, cnts_e14, 'Color','r', 'linewidth',2)
plot(edges_e15(1:end-1)+median( diff(edges_e15) )/2, cnts_e15, 'Color','b', 'linewidth',2)
plot(edges_e16(1:end-1)+median( diff(edges_e16) )/2, cnts_e16, 'Color','k', 'linewidth',2)

xlabel('log(amplitude)', 'fontsize', 20)
ylabel('Proportion', 'fontsize', 20)
title('Peak amplitude distributions', 'fontsize', 20)

subplot(2,2,3)
boxplot([e13_amp ; e14_amp; e15_amp ; e16_amp], [ones(size(e13_amp)) ; 2*ones(size(e14_amp)) ; 3*ones(size(e15_amp)) ; 4*ones(size(e16_amp))],'notch','on', 'whisker',inf,...
         'labels',{'E13','E14', 'E15', 'E16'})
ylim([2.75 3.25])
ylabel('log(amplitude)', 'fontsize', 20)
%%% FREQUENCIES

subplot(2,2,2)
hold on

e13_freq = ripple_specdist(1).freqs_pulse_max1;
e14_freq = ripple_specdist(2).freqs_pulse_max1;
e15_freq = ripple_specdist(3).freqs_pulse_max1;
e16_freq = ripple_specdist(4).freqs_pulse_max1;

[cnts_e13, edges_e13] = histcounts( e13_freq,15, 'Normalization', 'probability' );
[cnts_e14, edges_e14] = histcounts( e14_freq,15, 'Normalization', 'probability' );
[cnts_e15, edges_e15] = histcounts( e15_freq,20, 'Normalization', 'probability' );
[cnts_e16, edges_e16] = histcounts( e16_freq,22, 'Normalization', 'probability' ); 

plot(edges_e13(1:end-1)+median( diff(edges_e13) )/2, cnts_e13, 'Color',brown, 'linewidth',2)
plot(edges_e14(1:end-1)+median( diff(edges_e14) )/2, cnts_e14, 'Color','r', 'linewidth',2)
plot(edges_e15(1:end-1)+median( diff(edges_e15) )/2, cnts_e15, 'Color','b', 'linewidth',2)
plot(edges_e16(1:end-1)+median( diff(edges_e16) )/2, cnts_e16, 'Color','k', 'linewidth',2)
xlim([130 190])
xlabel('Frequency (Hz)', 'fontsize', 20)
ylabel('Proportion', 'fontsize', 20)
title('Frequency at peak amplitude', 'fontsize', 20)


subplot(2,2,4)
boxplot([e13_freq ; e14_freq; e15_freq ; e16_freq], [ones(size(e13_freq)) ; 2*ones(size(e14_freq)) ; 3*ones(size(e15_freq)) ; 4*ones(size(e16_freq))],'notch','on', 'whisker',inf,...
         'labels',{'E13','E14', 'E15', 'E16'})
ylim([145 170])
ylabel('Frequency (Hz)', 'fontsize', 20)

%% Durations 
%  We only look at ripples that had exactly one optotagged neuron
%  contributing
brown = [133 87 34]/254;

e13_dur = ripple_specdist(1).durs_max1 ;
e14_dur = ripple_specdist(2).durs_max1 ;
e15_dur = ripple_specdist(3).durs_max1 ;
e16_dur = ripple_specdist(4).durs_max1 ;

[cnts_e13, edges_e13] = histcounts( log10( e13_dur ),12, 'Normalization', 'probability' );
[cnts_e14, edges_e14] = histcounts( log10( e14_dur ),12, 'Normalization', 'probability' );
[cnts_e15, edges_e15] = histcounts( log10( e15_dur ),13, 'Normalization', 'probability' );
[cnts_e16, edges_e16] = histcounts( log10( e16_dur ),12, 'Normalization', 'probability' );

figure; 
set(gcf,'Position', [172 335 1254 378])

subplot(1,3,1)
hold on
plot(edges_e13(1:end-1)+median( diff(edges_e13) )/2, cnts_e13, 'Color',brown, 'linewidth',2)
plot(edges_e14(1:end-1)+median( diff(edges_e14) )/2, cnts_e14, 'Color','r', 'linewidth',2)
plot(edges_e15(1:end-1)+median( diff(edges_e15) )/2, cnts_e15, 'Color','b', 'linewidth',2)
plot(edges_e16(1:end-1)+median( diff(edges_e16) )/2, cnts_e16, 'Color','k', 'linewidth',2)
LogScale('x', 10)


xlabel('Ripple duration', 'fontsize', 20)
ylabel('Proportion', 'fontsize', 20)


subplot(1,3,2)
hold on
plot(sort(  log10( e13_dur +0.006*rand(size(e13_dur)) ) ), [ 1:length( e13_dur ) ] ./ length( e13_dur ) , 'Color', brown)
plot(sort( log10( e14_dur +0.006*rand(size(e14_dur)) ) ), [ 1:length( e14_dur ) ] ./ length( e14_dur ) , 'Color', 'r')
plot(sort( log10( e15_dur +0.006*rand(size(e15_dur))) ), [ 1:length( e15_dur ) ] ./ length( e15_dur ) , 'Color', 'b')
plot(sort( log10( e16_dur +0.006*rand(size(e16_dur)) )), [ 1:length(e16_dur) ] ./ length(e16_dur) , 'Color', 'k')
xlim( [-1.79 -0.84 ])
LogScale('x', 10)

xlabel('Ripple duration (s)', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)


subplot(1,3,3)
hold on

boxplot([log10( e13_dur ) ; log10(e14_dur ); log10(e15_dur) ; log10(e16_dur)], [ones(size(e13_dur)) ; 2*ones(size(e14_dur)) ; 3*ones(size(e15_dur)) ; 4*ones(size(e16_dur))],'notch','on', 'whisker',inf,...
         'labels',{'E13','E14', 'E15', 'E16'})
LogScale('y', 10)

ylabel('Ripple duration (s)', 'fontsize', 20)
xlabel('Birthdate', 'fontsize', 20)

%% Resample frequency according to E13 duration distribution

% How many p-values do we look at
NREP = 100;
pval = [];

tagDate = 3;   % 1 == E13 ; 2 == E14 ; 3 == E15

% Ripple durations / frequencies pooled
ripdurs_pooled = [ ripple_specdist(1).durs_max1 ; ripple_specdist(2).durs_max1 ; ripple_specdist(3).durs_max1 ];
ripfreqs_pooled = [ ripple_specdist(1).freqs_pulse_max1 ; ripple_specdist(2).freqs_pulse_max1 ; ripple_specdist(3).freqs_pulse_max1 ];

% Get the CDF for ripple durations in which E13 neurons participate
tag_durs_sort = sort( ripple_specdist(tagDate).durs_max1 );
probs = [ 1:length(e13_durs_sort) ] / length(tag_durs_sort);

% for kp = 1:NREP
    
fprintf('%d / %d \n', kp, NREP)

rng('shuffle')
indices = [];
for ii = 1:1000

    % Sample values from the pooled distribution according to the CDF of the
    % e13-participant distribution
    sampled_prob = rand;
    ind1 = find(probs >= sampled_prob); ind2 = find(probs < sampled_prob);
    if ~isempty(ind1) && ~isempty(ind2)
        inds = find( ripdurs_pooled >= tag_durs_sort(ind2(end)) & ripdurs_pooled <= tag_durs_sort(ind1(1)) );
        if isempty(inds)
            continue
        else
            indices = [ indices ; inds(randi(length(inds),1)) ];
        end
    end

end

% % If distribution of e13 participant rippled durations no different from
% % resampled pooled distribiton, test the frequencies of the e13 participant
% % against the frequencies of resampled pooled.
% if ranksum(ripple_specdist(tagDate).durs_max1, ripdurs_pooled(indices)) > .05
%     [~, pval(kp)] = ttest2( ripple_specdist(tagDate).freqs_pulse_max1, ripfreqs_pooled(indices));%, 'tail', 'left' );
% end
% 
% end
% 
% %%
% histogram(pval, 'normalization', 'probability','linestyle', 'none')
% xline(.01)
% xlim([0 .015])
% xlabel('P-value', 'fontsize', 16)
% ylabel('Proportion', 'fontsize', 16)
% 
% 
% %% This is for the plotting

close all

[cnts_tag_dur, edges_tag_dur] = histcounts( log10( ripple_specdist(tagDate).durs_max1 ),13, 'Normalization', 'probability' );
[cnts_pool_dur, edges_pool_dur] = histcounts( log10( ripdurs_pooled(indices) ),12, 'Normalization', 'probability' );

[cnts_tag, edges_tag] = histcounts( ripple_specdist(tagDate).freqs_pulse_max1,13, 'Normalization', 'probability' );
[cnts_pool, edges_pool] = histcounts( ripfreqs_pooled(indices),10, 'Normalization', 'probability' );

figure
set(gcf,'Position',[440   378   890   420])


subplot(1,2,1)
hold on
area(edges_tag_dur(1:end-1)+median( diff(edges_tag_dur) )/2, cnts_tag_dur, 'facecolor','k','linestyle', 'none', 'facealpha', 0.3)
area(edges_pool_dur(1:end-1)+median( diff(edges_pool_dur) )/2, cnts_pool_dur, 'facecolor','b','linestyle', 'none', 'facealpha', 0.3)
LogScale('x', 10)
xlabel('Ripple duration','fontsize', 16)
ylabel('Proportion','fontsize', 16)
% Adjust as we change time point
legend('E15', 'pooled, resampled')

subplot(1,2,2)
hold on

area(edges_tag(1:end-1)+median( diff(edges_tag) )/2, cnts_tag, 'facecolor','k','linestyle', 'none', 'facealpha', 0.3)
area(edges_pool(1:end-1)+median( diff(edges_pool) )/2, cnts_pool, 'facecolor','b','linestyle', 'none', 'facealpha', 0.3)
% xline([ mean(ripple_specdist(1).freqs_pulse_max1), mean(ripfreqs_pooled(indices)) ])
xlim([130 190])
xlabel('Ripple frequency','fontsize', 16)
ylabel('Proportion','fontsize', 16)
legend('E15', 'pooled, resampled')
fprintf('%f \n %f \n', ranksum(ripple_specdist(tagDate).durs_max1, ripdurs_pooled(indices)), ranksum(ripple_specdist(tagDate).freqs_pulse_max1, ripfreqs_pooled(indices)))


%% AMPLITUDES (MAX1, uncorrected)
%  peak amplitude within ripple events

[cnts_e13, edges_e13] = histcounts( ripple_specdist(1).powers_pulse_max1_uncorr,12, 'Normalization', 'probability' );
[cnts_e14, edges_e14] = histcounts( ripple_specdist(2).powers_pulse_max1_uncorr,12, 'Normalization', 'probability' );
[cnts_e15, edges_e15] = histcounts( ripple_specdist(3).powers_pulse_max1_uncorr,30, 'Normalization', 'probability' );
[edges_e16, cnts_e16] = histcounts( ripple_specdist(4).powers_pulse_max1_uncorr,30, 'Normalization', 'probability' );

figure; hold on
area(edges_e15(1:end-1)+median( diff(edges_e15) )/2, cnts_e15, 'facecolor','b','linestyle', 'none', 'facealpha', 0.3)
area(edges_e14(1:end-1)+median( diff(edges_e14) )/2, cnts_e14, 'facecolor','r','linestyle', 'none', 'facealpha', 0.3)
area(edges_e13(1:end-1)+median( diff(edges_e13) )/2, cnts_e13, 'facecolor',browb,'linestyle', 'none', 'facealpha', 0.3)

xlim([2.8 4])
xlabel('log(amplitude)', 'fontsize', 20)
ylabel('Proportion', 'fontsize', 20)
title('Peak amplitude distributions', 'fontsize', 20)



%% FREQUENCIES at peak amplitude (MAX1, uncorrected)

[cnts_tag, edges_tag] = histcounts( ripple_specdist(1).freqs_pulse_max1_uncorr,20, 'Normalization', 'probability' );
[cnts_e14, edges_e14] = histcounts( ripple_specdist(2).freqs_pulse_max1_uncorr,20, 'Normalization', 'probability' );
[cnts_e15, edges_e15] = histcounts( ripple_specdist(3).freqs_pulse_max1_uncorr,22, 'Normalization', 'probability' );
[cnts_e16, edges_e16] = histcounts( ripple_specdist(4).freqs_pulse_max1_uncorr,22, 'Normalization', 'probability' );

figure; hold on
area(edges_tag(1:end-1)+median( diff(edges_tag) )/2, cnts_tag, 'facecolor',brown,'linestyle', 'none', 'facealpha', 0.3)
area(edges_e14(1:end-1)+median( diff(edges_e14) )/2, cnts_e14, 'facecolor','r','linestyle', 'none', 'facealpha', 0.3)
area(edges_e15(1:end-1)+median( diff(edges_e15) )/2, cnts_e15, 'facecolor','b','linestyle', 'none', 'facealpha', 0.3)
area(edges_e16(1:end-1)+median( diff(edges_e16) )/2, cnts_e16, 'facecolor','k','linestyle', 'none', 'facealpha', 0.3)

xlim([100 200])
xlabel('Frequency', 'fontsize', 20)
ylabel('Proportion', 'fontsize', 20)
title('Frequency at peak amplitude', 'fontsize', 20)





