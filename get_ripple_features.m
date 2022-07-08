%% Extract features (amplitude/frequency/wavelet) of E13/14/15/16 participating ripples
%  This is a time consuming routine based on wavelets - computes the
%  ripple triggered wavelet spectrogram for each event, and stores this.

clear

isLocal = true;
nfreqs = 100;
frange = [60 250];
path_to_basepaths = '/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat';
load(path_to_basepaths)
 
cellpops = [{'basepaths_e13'}, {'basepaths_e14'}, {'basepaths_e15'}, {'basepaths_e16'} ];
 

for cellpop = cellpops
    basepaths = eval(cellpop{1});
    fprintf('Processing %s ...\n', cellpop{1})
    for ii = 1:length(basepaths)
        
        fprintf('%d/%d \n', ii, length(basepaths))
        % Get basepath
        basepath = alterPath( basepaths{ii}, isLocal );
        basename = bz_BasenameFromBasepath(basepath);

        % Get Kilosort folder
        load( fullfile(basepath, [basename '.ripples.events.mat']) )                               
        load( fullfile(basepath, [basename '.pulses.events.mat']) )
        load( fullfile(basepath, [ basename '.cell_metrics.cellinfo.mat']) )
        load( fullfile(basepath, [ basename '.spikes.cellinfo.mat']) )
        
        try
            load( fullfile(basepath, [ basename '.Behavior.mat']) )
            beh_ts = [behavior.timestamps(1) behavior.timestamps(end)];
        catch
            beh_ts = [];
        end
        % Optotagged pyramidal cells
        opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
        pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
        pyr_opto_id = opto_id & pyr_id;
        
        ripple_ints = ripples.timestamps( ~( InIntervals(ripples.timestamps(:,1), pulses.intsPeriods) |...   % Spontaneous ripples (exclude opto)
                                          InIntervals(ripples.timestamps(:,2), pulses.intsPeriods) ), : );
        ripple_ints = SubtractIntervals(ripple_ints, beh_ts);
        ripple_durs = diff(ripple_ints,[],2);
       
%         if isempty(ripple_ints)
%             continue
%         end
        
        % Get indicators for each tagged cell
        log_arr = false( size(ripple_ints,1), sum(pyr_opto_id) ); counter = 1;
        for kp = find( pyr_opto_id )
            [status, intt ] = InIntervals( spikes.times{kp}, ripple_ints);                               
            rip_p = unique( intt ); rip_p(rip_p == 0) = [];
            log_arr(rip_p,counter) = true;
            counter = counter+1;
        end
        % Ripples in which at least one tagged cell spiked
        ripp_any = logical( sum(log_arr, 2) > 0 );
        % Ripples in which exactly one tagged cell spiked
        ripp_max1 = logical( sum(log_arr, 2) > 0 & sum(log_arr, 2) < 2 );
        % Ripples in which at least but no more than two tagged cells
        % spiked
        ripp_max2 = logical( sum(log_arr, 2) > 0 & sum(log_arr, 2) < 3 );
        
        sr = ripples.detectorinfo.detectionparms.frequency;
        dt = 1/sr;

        % Intervals for computing wavelets
        p_vis = [ ripple_ints(:,1) - 0.1 ripple_ints(:,1) + 0.250 ];     % For visualization
        
        % RIPPLE WAVELET ! 

        % Wavelet within pulses
        wavespec = bz_WaveSpec(ripples.detectorinfo.detectionparms.lfp,'samplingRate', sr, 'frange', frange, 'intervals', ripple_ints, 'nfreqs', nfreqs);
        % length of longest 
        [~, lli] = max( diff(ripple_ints,[],2) );
        ni = sum( InIntervals(wavespec.timestamps, ripple_ints(lli,:)) );

        % Pull out individual wavelet spectrograms
        wavespec_ints = nan( ni, nfreqs,  length(ripple_ints) );
        for kp = 1:length(ripple_ints)
            tt = InIntervals( wavespec.timestamps, ripple_ints(kp,:) );
            wavespec_ints(1:sum(tt),:,kp) = wavespec.data( tt , :);
        end

        wavespec_ints_amp = abs(wavespec_ints);
        % Power spectrum in each ripple
        ps_rip =  nanmean(wavespec_ints_amp,1);
        ps_rip =  log10( squeeze(ps_rip) );

        % VISUALIZATION WAVELET AVERAGE ! 

        % Wavelet within pulses
        wavespec = bz_WaveSpec(ripples.detectorinfo.detectionparms.lfp,'samplingRate', sr, 'frange', frange, 'intervals', p_vis, 'nfreqs', nfreqs);
        % length of longest 
        [~, lli] = max( diff(p_vis,[],2) );
        ni = sum( InIntervals(wavespec.timestamps, p_vis(lli,:)) );

        % Pull out individual wavelet spectrograms
        wavespec_ints = nan( ni, nfreqs,  length(p_vis) );
        for kp = 1:length(p_vis)
            tt = InIntervals( wavespec.timestamps, p_vis(kp,:) );
            wavespec_ints(1:sum(tt),:,kp) = wavespec.data( tt, :);
        end
        wavespec_ints_amp = abs(wavespec_ints);
        wavespec_tvec = 0:dt:(ni*dt)-dt; wavespec_tvec = wavespec_tvec-0.1;
        
        % Pull out wavelep spectrograms from ripples in which tagged cells
        % participate
        tmp = wavespec_ints_amp(:,:,ripp_any);
        wavespec_mean_any = nanmean(tmp,3);
        tmp = wavespec_ints_amp(:,:,ripp_max1);
        wavespec_mean_max1 = nanmean(tmp,3);
        tmp = wavespec_ints_amp(:,:,ripp_max2);
        wavespec_mean_max2 = nanmean(tmp,3);

        ripple_spectral = struct();

        % to save
        ripple_spectral.ripple_ints = ripple_ints;
        ripple_spectral.ripple_durs = ripple_durs;
        ripple_spectral.freqs = wavespec.freqs;
        ripple_spectral.ripple_ps = ps_rip;
        ripple_spectral.wavespec_mean_any = wavespec_mean_any;
        ripple_spectral.wavespec_mean_max1 = wavespec_mean_max1;
        ripple_spectral.wavespec_mean_max2 = wavespec_mean_max2;
        ripple_spectral.ripp_any = ripp_any;
        ripple_spectral.ripp_max1 = ripp_max1;
        ripple_spectral.ripp_max2 = ripp_max2;
        ripple_spectral.wavespec_tvec = wavespec_tvec;

        save( fullfile(basepath, [basename '.ripple_spectral.mat']), 'ripple_spectral')

 
    end
 
end