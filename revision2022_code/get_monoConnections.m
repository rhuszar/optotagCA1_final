%% Update the cell metrics
%  See if interneurons are rate modulated around the time of short light pulses

load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac_rev1.mat")
postPulse_win = 0.02;  
isLocal = true;
for ip = 79%:length(basepaths_all)
    
    basepath = alterPath( basepaths_all{ip}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    
    fprintf('%d/%d: %s\n', ip, length(basepaths_all), basename)
    
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    % Check if we've extracted interneuron modulation already
    if isfield( cell_metrics.optoTag, 'int_ratemod' ) && ...
            isfield( cell_metrics.optoTag, 'int_ratemod_pval' )
        disp('int modulation extracted')
        continue
    end
    
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    % load(fullfile(basepath, [basename '.mono_res.cellinfo.mat']))
    load(fullfile(basepath, [basename '.pulses.events.mat']))
    
    pulses_sh_ints = pulses.timestamps( diff( pulses.timestamps' ) <0.005,: );
    post_pulse_sh_ints = [ pulses_sh_ints(:,2) pulses_sh_ints(:,2) + postPulse_win ];
    pre_pulses_sh_ints = post_pulse_sh_ints - diff(post_pulse_sh_ints,[],2) - .01;

    
    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    int_id = cellfun(@(x) ~isempty(regexp(x, 'Interneuron', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id & pyr_id;
    
    cell_metrics.optoTag.int_ratemod = nan(1,length(int_id));
    cell_metrics.optoTag.int_ratemod_pval = nan(1,length(int_id)); 
    
    for kp = find( int_id )
        spk = spikes.times{kp};
        [~, interval] = InIntervals( spk, post_pulse_sh_ints );
        stim_rate = histoc(interval(interval>0), 1:size(pre_pulses_sh_ints,1)) ./ diff(pre_pulses_sh_ints,[],2);
        [~, interval] = InIntervals( spk, pre_pulses_sh_ints );
        pre_stim_rate = histoc(interval(interval>0), 1:size(pre_pulses_sh_ints,1)) ./ diff(pre_pulses_sh_ints,[],2);
        % compute sign test and store
        [p,~] = signtest(stim_rate, pre_stim_rate);
        cell_metrics.optoTag.int_ratemod(kp) =  ( nanmean(stim_rate)  - nanmean( pre_stim_rate ) ) / nanmean( pre_stim_rate );
        cell_metrics.optoTag.int_ratemod_pval(kp) =  p;
    
    end
    
    save(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']),'cell_metrics', '-v7.3')

end

%%

for ip = 1:length(basepaths_all)
    basepath = alterPath( basepaths_all{ip}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    if exist( fullfile(basepath, 'monoConnections_pretag_v2.mat') )
        movefile(fullfile(basepath, 'monoConnections_pretag_v2.mat'), fullfile(basepath, 'monoConnections_pretag.mat'), 'f')
    end
end


%% Calculate drive from tagged PYR onto interneurons - both average and 
%  split by ISI ; furthermore, store information about rate changes to each
%  interneuron following pulses
%  Store this information with each session


isLocal = true;
binsize = .0008;
duration = 0.2;

for ip = 1:length(basepaths_all)
    
    uid_pair = [];
    int_ratemod = [];
    spktrans_mean = [];
    spktrans_gamma = [];
    spktrans_splitISI = [];
    
    
    fprintf('%d/%d\n', ip, length(basepaths_all))
    basepath = alterPath( basepaths_all{ip}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    
    outfile = fullfile(basepath, 'monoConnections_pretag_v2.mat');
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.mono_res.cellinfo.mat']))
    load(fullfile(basepath, [basename '.pulses.events.mat']))
    load( fullfile(basepath, [basename '.session.mat']) );
    reclen = session.timeSeries.adc.nSamples / session.timeSeries.adc.sr;
    try
    % Behavior
        load( fullfile(basepath, [basename '.Behavior.mat']) )
        if exist( fullfile(basepath, [basename '.Tracking.Behavior.mat']) ) 
            % If multiple behavioral sessions, consider homecage data surrounding
            % the first session only
            load( fullfile(basepath, [basename '.Tracking.Behavior.mat'])  )
            beh_ts = tracking.events.subSessions(1,:);
            discard = tracking.events.subSessions(2,1);
        else
            load( fullfile(basepath, [basename '.MergePoints.events.mat']) )
            [~, ff] = InIntervals( behavior.timestamps(1), MergePoints.timestamps  );
            beh_ts = MergePoints.timestamps(ff,:);
            discard = reclen;
        end
        disp('Behavior session!')
    catch
    % No behavior
        beh_ts = [];
        discard = reclen;
        disp('No behavior session!')
    end
    
    % Get non pulse period as interval
    nonpulse_period = SubtractIntervals( [0 discard], [ pulses.intsPeriods ; beh_ts]);
    
    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    int_id = cellfun(@(x) ~isempty(regexp(x, 'Interneuron', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = find( opto_id & pyr_id );
    int_opto_id = opto_id & int_id;
   
    
    presyn_tagged = arrayfun(@(x) x == mono_res.pyr2int(:,1), pyr_opto_id, 'UniformOutput', false );
    % Loop over connections where presynaptic neuron is tagged
    pretag_id = find( sum( cell2mat( presyn_tagged ), 2) )';
    for jj = pretag_id
        
        pyr2int = mono_res.pyr2int(jj,:);
        % Presynaptic spikes outside pulses
        pre = spikes.times{pyr2int(1)}( InIntervals( spikes.times{pyr2int(1)}, nonpulse_period ) );
        post = spikes.times{pyr2int(2)};
        ses = ShortTermCCG(pre, post, binsize, duration, 'time', [ logspace(log10(5),log10(1000),15)/1000 inf]);
        
        % Also get the average spike transmission probability
        [ccg,~] = CCG([pre ; post], [ones(size(pre)) ; 2*ones(size(post))], 'binsize', binsize, 'duration', duration);
        [ccgSpkTrans,~,~,~] = GetTransProb(ccg(:,1,2),length(pre),binsize);
        % Average spike transmission when pre is firing at gamma
        % gammaSpk = find( diff( pre ) >= 0.0156 &  diff( pre ) <= 0.0484 )+1;
        gammaSpk = find( diff( pre ) >= 0.01 &  diff( pre ) <= 0.04 )+1;
        [ccg_gamma,~] = CCG([pre(gammaSpk) ; post], [ones(length(gammaSpk),1) ; 2*ones(size(post))], 'binsize', binsize, 'duration', duration);
        [ccgSpkTrans_gamma,~,~,~] = GetTransProb(ccg_gamma(:,1,2),length(gammaSpk),binsize);
        
        % Store everything
        spktrans_splitISI = [ spktrans_splitISI ; ses.trans ];
        spktrans_mean = [spktrans_mean ; ccgSpkTrans]; 
        spktrans_gamma = [spktrans_gamma ; ccgSpkTrans_gamma]; 
        uid_pair = [ uid_pair ; pyr2int ];
        % Rate modulation of the interneuron around pulse
        int_ratemod = [ int_ratemod ; [ double( int_opto_id( pyr2int(2) ) )  cell_metrics.optoTag.int_ratemod( pyr2int(2) ) cell_metrics.optoTag.int_ratemod_pval( pyr2int(2) ) ] ];
        
    end
    if ~isempty( pretag_id )
        save(outfile, 'spktrans_splitISI', 'spktrans_mean', 'spktrans_gamma', 'uid_pair', 'int_ratemod')
    end
end
