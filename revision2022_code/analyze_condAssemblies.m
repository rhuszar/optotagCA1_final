%%
%{

Conditional assembly expression analyses
- in each session, discover assemblies conditioned on spiking of opto
tagged neurons
- for each assembly extracted, track it throughout the entire track running
period

From this setup, it isn't immediately obvious where in space the assemblies
will be active, to what extent this activity overlaps with the spiking of
the target neuron, etc.


% Conditional assembly detection on the track

% pyrInds(pyrOnlyInds_opto)    ... UIDs of optotagged PYRs
% pyrInds_opto                 ... UIDs of optotagged PYRs
% pyrOnlyInds_opto             ... IDs of optotagged PYRs in an array of
%                                  only PYRs
% assemblies_cond( 1 ).pyrUID_cond   ... UIDs of PYRs after one has been
%                                        taken out
% assemblies_cond( 1 ).indOpto       ... indices into pyrUID_cond - optotagged units 
% assemblies_cond( 1 ).pyrUID_opto   ... UID of optotagged cell being
%                                        conditioned upon

%}
%% Extract place maps in conditional assemblies

%path = '/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV5_OUT_V1';
path = 'C:\Users\Roman\Google Drive\buzz_workspace\optotag\DATA\runAssembly_trackV5_OUT_V1';
%%

fils = dir('C:\Users\Roman\Documents\DATA\runAssembly_trackVrev_condTrack_OUT');
fils = {fils.name}; fils(1:2) = [];
isLocal = true;
exception_animal = {'e13_26m1'};
for kp = 1:length(basepaths_beh)
    
    fprintf('\n%d/%d\n', kp, length(basepaths_beh))
    basepath_loc = alterPath( basepaths_beh{kp}, isLocal );
    basename = bz_BasenameFromBasepath( basepath_loc );
    
    
    disp(basename)
    fils_id = find( cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils) );
    if ~any( fils_id )
        continue
    end
    
    % Load assembly output
    load( fils{fils_id} );
    
    % If firing maps have been extracted, skip
    if isfield(assemblies_cond, 'firingMaps')
        fprintf('Conditional maps extracted...\n');
        continue
    end
    
    for jp = 1:length( assemblies_cond )
        % Whenever we have fewer observations than number of neurons, ICA
        % can't be run so these entries will be empty
        if isempty( assemblies_cond(jp).assembly_act_track )
            continue
        end
        if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exception_animal) )
            firingMaps = getFiringMaps_rh_multSess('basepath', basepath_loc, 'spike_times', assemblies_cond(jp).assembly_act_track, 'forceReload', true);
        else
            firingMaps = getFiringMaps_rh_v2('basepath', basepath_loc, 'spike_times', assemblies_cond(jp).assembly_act_track, 'forceReload', true);
        end
        assemblies_cond(jp).firingMaps = firingMaps;
%         placeFieldStats_trial = [];
%         for ip = 1:length( assemblies_cond(jp).assembly_act_track )
%             ratemap = [];
%             ratemap(1,:,:) = firingMaps.rateMaps_trial{ip,1};
%             [ fields_left ] = bz_getPlaceFields1D('ratemap', ratemap,'percentThreshold',.5);
%             ratemap = [];
%             ratemap(1,:,:) = firingMaps.rateMaps_trial{ip,2};
%             [ fields_right ] = bz_getPlaceFields1D('ratemap', ratemap,'percentThreshold',.5);
%             placeFieldStats_trial = [placeFieldStats_trial ; [fields_left fields_right]];
%         end
%         assemblies_cond(jp).placeFieldStats_trial = placeFieldStats_trial;
    end
    
    save( fils{fils_id}, 'assemblies_cond', '-append' );

end

%% Visualize the raster to show conditional assembly detection

datpath = '/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA';

fwhm = 0.025;
MAX = 101;


basepath_loc = basepaths_beh{10};
basename = bz_BasenameFromBasepath(basepath_loc);
pathparts = strsplit(basepath_loc,filesep);
basepath_loc = fullfile(datpath, pathparts{end-2}, pathparts{end-1}, pathparts{end});

v = load( [ basename '.mat'] );

load(fullfile(basepath_loc, [basename '.spikes.cellinfo.mat']))
load(fullfile(basepath_loc, [basename '.cell_metrics.cellinfo.mat']))
load(fullfile(basepath_loc, [basename '.Behavior.mat']))

interval = [behavior.timestamps(1) behavior.timestamps(end)];

% Cell type ids
opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
pyr_opto_id = opto_id & pyr_id;
pyrInds = find(pyr_id);
pyrInds_opto = find(pyr_opto_id);
pyrOnlyInds_opto = arrayfun(@(x) find( x==pyrInds ), pyrInds_opto);

spkCounts_track = bz_SpktToSpkmat(spikes.times(pyrInds), 'dt',  fwhm, 'binsize', fwhm, 'win', interval);

%% Visualize raster

kp = 1;

vals = spkCounts_track.data(1:MAX, pyrOnlyInds_opto(kp));
tvec = spkCounts_track.timestamps(1:MAX);


subplot(2,1,1)
imagesc(tvec,1:length(pyrInds), spkCounts_track.data(1:MAX,:)' )
colormap( flip( bone ))
xlabel('Time (s)', 'fontsize', 20)
ylabel('PYR', 'fontsize', 20)

subplot(2,1,2)
imagesc(tvec,1:length(pyrInds), spkCounts_track.data(1:MAX,:)' )
colormap( flip( bone ))
hold on
xline(tvec( logical( vals )))
xlabel('Time (s)', 'fontsize', 20)
ylabel('PYR', 'fontsize', 20)

%% Visualize assemblies as a stem plot


assemblyStemPlot(v.assemblies_cond(kp).assemblies, v.assemblies_cond(kp).indOpto)



%% Ratemap correlation between optotagged neurons, and assemblies detected conditioned
%  on the spiking of these neurons
%  Code for figure 7C (original submission, long format)

%  TREAT LEFT AND RIGHT RATEMAPS SEPARATELY
%  USEFUL FOR PLOTTING

% cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV5_OUT')
cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV5_OUT_V1')
load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
fils = dir;
fils = {fils.name}; fils(1:2) = [];
nsd = 2;

% Indices of what's after the maze stem


corr_tag_leftALL = [];
corr_tag_rightALL = [];
corr_nontag_leftALL = [];
corr_nontag_rightALL = [];
ntag_leftALL = [];
ntag_rightALL = [];
for ip = 1:length(basepaths_beh)
    
    fprintf('%d/%d\n', ip, length(basepaths_beh))
    basepath_loc = alterPath( basepaths_beh{ip}, true );
    basename = bz_BasenameFromBasepath( basepath_loc );
    
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    
    disp(basename)
    fils_id = find( cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils) );
    if ~any( fils_id )
        continue
    end
    
    % Load assembly output
    v = load( fils{fils_id} );

    load(fullfile(basepath_loc, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath_loc, [basename '.firingMapsAvg_v2.cellinfo.mat']))

    for kp = 1:length( v.assemblies_cond )

        % UID of optotagged cell being conditioned upon 
        % for assembly extraction
        optoUID = v.assemblies_cond(kp).pyrUID_opto;

        % Find assemblies with at least one tagged neuron as a high contributor
        HC_indicator = zscore( v.assemblies_cond(kp).assemblies ) > repmat( nsd, size(v.assemblies_cond(kp).assemblies));
        %taggedAssembly_id = logical( sum( HC_indicator( v.assemblies_cond(kp).indOpto, : ),1 ))';
        taggedAssembly_id = sum( HC_indicator( v.assemblies_cond(kp).indOpto, : ),1 );
        
%         % Skip if no tagged neurons are present
        if ~any(taggedAssembly_id)
            continue
        end
        subplot(1,2,1)
        plot( firingMaps.rateMaps{optoUID}{1}, '-k', 'linewidth', 1.5 )
        hold on
        tag_maps_left = [];
        % Overlay the opto tagged assemblies
        for jp = find( taggedAssembly_id )
            plot( v.assemblies_cond(kp).firingMaps.rateMaps{jp}{1}, '-b' )
            tag_maps_left = [ tag_maps_left v.assemblies_cond(kp).firingMaps.rateMaps{jp}{1}' ];
        end
        corr_tag_left = corr( [ firingMaps.rateMaps{optoUID}{1}' tag_maps_left ] ) ;
        corr_tag_leftALL = [corr_tag_leftALL ; {corr_tag_left(1,2:end)'}];
        ntag_leftALL = [ ntag_leftALL ; {taggedAssembly_id( find( taggedAssembly_id ) )'}];
        
        nontag_maps_left = [];
        for jp = find( ~taggedAssembly_id )
            plot( v.assemblies_cond(kp).firingMaps.rateMaps{jp}{1}, '-r' )
            nontag_maps_left = [ nontag_maps_left v.assemblies_cond(kp).firingMaps.rateMaps{jp}{1}' ];
        end
        corr_nontag_left = corr( [ firingMaps.rateMaps{optoUID}{1}' nontag_maps_left ] ) ;
        corr_nontag_leftALL = [corr_nontag_leftALL ; {corr_nontag_left(1,2:end)'}];
        xlim([0 205])
        
        %
        subplot(1,2,2)
        plot( firingMaps.rateMaps{optoUID}{2}, '-k', 'linewidth', 1.5 )
        hold on
        tag_maps_right = [];
        % Overlay the opto tagged assemblies
        for jp = find( taggedAssembly_id )
            plot( v.assemblies_cond(kp).firingMaps.rateMaps{jp}{2}, '-b' )
            tag_maps_right = [ tag_maps_right v.assemblies_cond(kp).firingMaps.rateMaps{jp}{2}' ];
        end
        corr_tag_right = corr( [ firingMaps.rateMaps{optoUID}{2}' tag_maps_right ] ) ;
        corr_tag_rightALL = [corr_tag_rightALL ; {corr_tag_right(1,2:end)'}];
        ntag_rightALL = [ ntag_rightALL ; {taggedAssembly_id( find( taggedAssembly_id ) )'}];
        
        nontag_maps_right = [];
        for jp = find( ~taggedAssembly_id )
            plot( v.assemblies_cond(kp).firingMaps.rateMaps{jp}{2}, '-r' )
            nontag_maps_right = [ nontag_maps_right v.assemblies_cond(kp).firingMaps.rateMaps{jp}{2}' ];
        end
        corr_nontag_right = corr( [ firingMaps.rateMaps{optoUID}{2}' nontag_maps_right ] ) ;
        corr_nontag_rightALL = [corr_nontag_rightALL ; {corr_nontag_right(1,2:end)'}];
        xlim([0 205])
        
        uiwait
    end


end

%%

tag_corr = [cell2mat(corr_tag_rightALL) ; cell2mat(corr_tag_leftALL)];
nontag_corr = [cell2mat(corr_nontag_rightALL) ; cell2mat(corr_nontag_leftALL)];


plot(sort( tag_corr ), [ 1:length(tag_corr) ] ./ length(tag_corr) ); hold on
plot(sort( nontag_corr ), [ 1:length(nontag_corr) ] ./ length(nontag_corr) )

%% Ratemap correlation between optotagged neurons, and assemblies detected conditioned
%  Code for figure 7D (original submission, long format)

%  on the spiking of these neurons
%  COMBINE LEFT AND RIGHT RATEMAPS

% cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV5_OUT')
cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV5_OUT_V1')
% load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/basepaths_mac_rev1.mat')
%%
fils = dir;
fils = {fils.name}; fils(1:2) = [];
nsd = 2;

post_stem_inds = 75:205;

corr_tag_ALL = [];
corr_nontag_ALL = [];
ntag_ALL = [];
exclude_animal = {'e15_13f1', 'e16_3m2'};
for ip = 1:length(basepaths_beh)
    
    fprintf('%d/%d\n', ip, length(basepaths_beh))
    basepath_loc = alterPath( basepaths_beh{ip}, true );
    basename = bz_BasenameFromBasepath( basepath_loc );
    
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    
    disp(basename)
    fils_id = find( cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils) );
    if ~any( fils_id )
        continue
    end
    
    
    % Load assembly output
    v = load( fils{fils_id} );

    load(fullfile(basepath_loc, [basename '.spikes.cellinfo.mat']))
    % Check if we're dealing with a session that has multiple beh sessions
    if exist(fullfile(basepath_loc, [basename '.Tracking.Behavior.mat']))
        load(fullfile(basepath_loc, [basename '.firingMapsAvg_multSess.cellinfo.mat']))
    else
        load(fullfile(basepath_loc, [basename '.firingMapsAvg_v2.cellinfo.mat']))
    end

    for kp = 1:length( v.assemblies_cond )

        if isempty( v.assemblies_cond(kp).assembly_act_track )
            continue
        end
        % UID of optotagged cell being conditioned upon 
        % for assembly extraction
        optoUID = v.assemblies_cond(kp).pyrUID_opto;

        % Find assemblies with at least one tagged neuron as a high contributor
        HC_indicator = zscore( v.assemblies_cond(kp).assemblies ) > repmat( nsd, size(v.assemblies_cond(kp).assemblies));
        %taggedAssembly_id = logical( sum( HC_indicator( v.assemblies_cond(kp).indOpto, : ),1 ))';
        taggedAssembly_id = sum( HC_indicator( v.assemblies_cond(kp).indOpto, : ),1 );
        
%         % Skip if no tagged neurons are present
%         if ~any(taggedAssembly_id)
%             continue
%         end
%         subplot(1,2,1)
%         plot( firingMaps.rateMaps{optoUID}{1}, '-k', 'linewidth', 1.5 )
%         hold on
        tag_maps = [];
        % Overlay the opto tagged assemblies
        for jp = find( taggedAssembly_id )
            % plot( v.assemblies_cond(kp).firingMaps.rateMaps{jp}{1}, '-b' )
            tag_maps = [ tag_maps [ v.assemblies_cond(kp).firingMaps(1).rateMaps{jp}{1}(post_stem_inds)' ; ...
                                    v.assemblies_cond(kp).firingMaps(1).rateMaps{jp}{2}(post_stem_inds)' ] ];
        end
        corr_tag = corr( [ [ firingMaps(1).rateMaps{optoUID}{1}(post_stem_inds)' ; firingMaps(1).rateMaps{optoUID}{2}(post_stem_inds)' ] tag_maps ] ) ;
        corr_tag_ALL = [corr_tag_ALL ; {corr_tag(1,2:end)'}];
        ntag_ALL = [ ntag_ALL ; {taggedAssembly_id( find( taggedAssembly_id ) )'}];
        
        nontag_maps = [];
        for jp = find( ~taggedAssembly_id )
            % plot( v.assemblies_cond(kp).firingMaps.rateMaps{jp}{1}, '-r' )
            nontag_maps = [ nontag_maps [ v.assemblies_cond(kp).firingMaps(1).rateMaps{jp}{1}(post_stem_inds)' ; ...
                                          v.assemblies_cond(kp).firingMaps(1).rateMaps{jp}{2}(post_stem_inds)' ]  ];
        end
        corr_nontag = corr( [ [ firingMaps(1).rateMaps{optoUID}{1}(post_stem_inds)' ; firingMaps(1).rateMaps{optoUID}{2}(post_stem_inds)' ] nontag_maps ] ) ;
        corr_nontag_ALL = [corr_nontag_ALL ; {corr_nontag(1,2:end)'}];

   
        
        %uiwait
    end
end

%%

corr_tagtag = cell2mat( corr_tag_ALL );
corr_nontag = cell2mat( corr_nontag_ALL );

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( corr_tagtag), [ 1:length(corr_tagtag) ] ./ length(corr_tagtag), '-b');
hold on
plot(sort( corr_nontag), [ 1:length(corr_nontag) ] ./ length(corr_nontag), '-r')
xlim( [-0.4 0.8] )
xlabel('Rate-assembly map correlation','fontsize', 14)
ylabel('Cumulative proportion','fontsize', 14)


figure
boxplot([corr_nontag ; corr_tagtag ], [ones(size(corr_nontag)) ; 2*ones(size(corr_tagtag))],'notch','on', 'whisker',inf)
ylim([-0.15 0.3])

%% Assembly expression conditioned on spikes of optotagged neurons

tag_corr = cell2mat(corr_tag_ALL);
nontag_corr = cell2mat(corr_nontag_ALL);

plot(sort( tag_corr ), [ 1:length(tag_corr) ] ./ length(tag_corr) ); hold on
plot(sort( nontag_corr ), [ 1:length(nontag_corr) ] ./ length(nontag_corr) )

%% Assembly expression strength around the spikes of held out neurons
%  Code for figure 7B (original submission, long format)

%  Split assemblies by whether they contain coborn neurons or not
cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV5_OUT_V1')
%load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/basepaths_mac_rev1.mat')
%%
fils = dir;
fils = {fils.name}; fils(1:2) = [];
nsd = 2;

binsize = .002;
duration = 1;

expr_tag_ALL = [];
cc_tag_ALL = [];
expr_nontag_ALL = [];
cc_nontag_ALL = [];
exclude_animal = {'e15_13f1', 'e16_3m2'};
for ip = 1:length(basepaths_beh)
    
    fprintf('%d/%d\n', ip, length(basepaths_beh))
    basepath_loc = alterPath( basepaths_beh{ip}, true );
    basename = bz_BasenameFromBasepath( basepath_loc );
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    disp(basename)
    fils_id = find( cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils) );
    if ~any( fils_id )
        continue
    end
    
    % Load assembly output
    v = load( fils{fils_id} );

    load(fullfile(basepath_loc, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath_loc, [basename '.Behavior.mat']))
    beh_tvec = behavior.timestamps( behavior.masks.recording == 1 );
    beh_int = [beh_tvec(1) beh_tvec(end)];

    for kp = 1:length( v.assemblies_cond )
        
        if isempty( v.assemblies_cond(kp).assembly_act_track )
            continue
        end
        % UID of optotagged cell being conditioned upon 
        % for assembly extraction
        optoUID = v.assemblies_cond(kp).pyrUID_opto;
        % Spikes of optotagged cell while on the track
        spk = spikes.times{optoUID};
        spk = spk( InIntervals(spk, beh_int) );

        % Find assemblies with at least one tagged neuron as a high contributor
        HC_indicator = zscore( v.assemblies_cond(kp).assemblies ) > repmat( nsd, size(v.assemblies_cond(kp).assemblies));
        %taggedAssembly_id = logical( sum( HC_indicator( v.assemblies_cond(kp).indOpto, : ),1 ))';
        taggedAssembly_id = sum( HC_indicator( v.assemblies_cond(kp).indOpto, : ),1 );
        
        for jp = find( taggedAssembly_id )
            assembly_spk = v.assemblies_cond(kp).assembly_act_track{jp};
            ccg = CCG([spk ; assembly_spk], [ones(size(spk)) ; 2.*ones(size(assembly_spk))], 'binsize', binsize, 'duration', duration);
            cc_tag_ALL = [ cc_tag_ALL ; ccg(:,1,2)' ./ length(spk) ];
        end
        expr_tag_ALL = [ expr_tag_ALL ; mean( v.assemblies_cond(kp).activities(logical( taggedAssembly_id ),:), 2)];
        
        for jp = find( ~taggedAssembly_id )
            assembly_spk = v.assemblies_cond(kp).assembly_act_track{jp};
            [ ccg, tvec ] = CCG([spk ; assembly_spk], [ones(size(spk)) ; 2.*ones(size(assembly_spk))], 'binsize', binsize, 'duration', duration);
            cc_nontag_ALL = [ cc_nontag_ALL ; ccg(:,1,2)' ./ length(spk) ];
        end
        expr_nontag_ALL = [ expr_nontag_ALL ; mean( v.assemblies_cond(kp).activities( ~taggedAssembly_id,:), 2)];
    end
end

%%
plot( tvec, mean( cc_tag_ALL ), '-b' )
hold on
plot( tvec, smooth( mean( cc_tag_ALL ), 5 ), '-b' )
plot( tvec,mean( cc_nontag_ALL ), '-r' )
plot( tvec,smooth( mean( cc_nontag_ALL ), 5), '-r' )
ylim([0.001    0.0028])
xlabel('Lag (s)', 'fontsize', 14)
ylabel('Prob. of assembly expression', 'fontsize', 14)

%%

[counts,edges] = histcounts(log10(expr_tag_ALL), 14,'normalization', 'probability');
plot(10.^edges(1:end-1), counts, '-b')
set(gca, 'xscale','log')

hold on

[counts,edges] = histcounts(log10(expr_nontag_ALL), 15,'normalization', 'probability');
plot(10.^edges(1:end-1), counts, '-r')
set(gca, 'xscale','log')

xlabel('Assembly expression (projected z-score)', 'fontsize', 14)
ylabel('Proportion', 'fontsize', 14)

%% Compute theta cofiring with tagged / untagged assembly members, and assembly nonmembers
%  Code for figure 7E (original submission, long format)

cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV5_OUT_V1')
%load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
load('/Volumes/GoogleDrive/My Drive/buzz_workspace/optotag/basepaths_mac_rev1.mat')
%%
fils = dir;
fils = {fils.name}; fils(1:2) = [];
nsd = 2;

thetacorr_tagtag = [];
theta_tagtag_rates = [];
id_tagtag = [];

thetacorr_ntagMem = [];
theta_ntagMem_rates = [];
id_ntagMem = [];

thetacorr_ntagNMem = [];
theta_ntagNMem_rates = [];
id_ntagNMem = [];

randLabel = 2000*rand;
exclude_animal = {'e15_13f1', 'e16_3m2'};
for ip = 1:length(basepaths_beh)
    
    fprintf('%d/%d\n', ip, length(basepaths_beh))
    basepath_loc = alterPath( basepaths_beh{ip}, true );
    basename = bz_BasenameFromBasepath( basepath_loc );
    tmp = strsplit( basename, '_' );
    
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    
    disp(basename)
    fils_id = find( cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils) );
    if ~any( fils_id )
        continue
    end
    
    % Load assembly output
    v = load( fils{fils_id} );

    load(fullfile(basepath_loc, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath_loc, [basename '.thetaCycles.events.mat']))
    load(fullfile(basepath_loc, [basename '.SleepState.states.mat']))
    
    if exist( fullfile(basepath_loc, [basename '.Tracking.Behavior.mat']) ) 
        % If multiple behavioral sessions, consider homecage data surrounding
        % the first session only
        load( fullfile(basepath_loc, [basename '.Tracking.Behavior.mat'])  )
        beh_ts = tracking.events.subSessions(1,:);
    else
        load( fullfile(basepath_loc, [basename '.MergePoints.events.mat']) )
        [~, ff] = InIntervals( behavior.timestamps(1), MergePoints.timestamps  );
        beh_ts = MergePoints.timestamps(ff,:);
    end

    theta_ints_run = Theta.cycles( InIntervals( Theta.cycles(:,1), beh_ts ) & InIntervals( Theta.cycles(:,2), beh_ts ), : ); 
    theta_ints_run = theta_ints_run( InIntervals( theta_ints_run(:,1), SleepState.ints.WAKEstate ) & InIntervals( theta_ints_run(:,2), SleepState.ints.WAKEstate ), : );
    theta_len = sum( diff( theta_ints_run' ) );
    
    % This is the heavy lifting
    sameBD_members = []; diffBD_members = []; nonmembers = [];
    for kp = 1:length( v.assemblies_cond )
        
        if isempty( v.assemblies_cond(kp).assembly_act_track )
            continue
        end
        % tagged cell that was conditioned on - index into PYRs only 
        optoUID = v.pyrOnlyInds_opto(kp);

        % Find assemblies with at least one tagged neuron as a high contributor
        HC_indicator = zscore( v.assemblies_cond(kp).assemblies ) > repmat( nsd, size(v.assemblies_cond(kp).assemblies));
        % assembly members
        members = find( any( HC_indicator, 2) );
        nmembers = setdiff( 1:size(HC_indicator,1), members');
        
        t_uid =v.assemblies_cond(kp).pyrUID_cond(intersect( members, v.assemblies_cond(kp).indOpto' ) );
        t_PYRuid = arrayfun(@(x) find( x==v.pyrInds ), t_uid);
        sameBD_members = [sameBD_members ; [repmat(optoUID, length(t_PYRuid),1 ) t_PYRuid'] ];
        
        nt_uid = v.assemblies_cond(kp).pyrUID_cond(setdiff( members, v.assemblies_cond(kp).indOpto' ) );
        nt_PYRuid = arrayfun(@(x) find( x==v.pyrInds ), nt_uid);
        diffBD_members = [diffBD_members ; [repmat(optoUID, length(nt_PYRuid),1 ) nt_PYRuid'] ];
        
        nm_uid = v.assemblies_cond(kp).pyrUID_cond( nmembers );
        nm_PYRuid = arrayfun(@(x) find( x==v.pyrInds ), nm_uid);
        nonmembers = [nonmembers ; [repmat(optoUID, length(nm_PYRuid),1 ) nm_PYRuid'] ];
    end
    
    if ~isempty(sameBD_members)
        sw_ind = sameBD_members(:,2) <  sameBD_members(:,1);
        sameBD_members(sw_ind,:) = [sameBD_members(sw_ind,2) sameBD_members(sw_ind,1)];
        ind = sub2ind([length( v.pyrInds ) length( v.pyrInds )], sameBD_members(:,1), sameBD_members(:,2) );
        [~, ia] = unique(ind);
        sameBD_members = sameBD_members(ia,:);
    end
    if ~isempty(diffBD_members)
        sw_ind = diffBD_members(:,2) <  diffBD_members(:,1);
        diffBD_members(sw_ind,:) = [diffBD_members(sw_ind,2) diffBD_members(sw_ind,1)];
        ind = sub2ind([length( v.pyrInds ) length( v.pyrInds )], diffBD_members(:,1), diffBD_members(:,2) );
        [~, ia] = unique(ind);
        diffBD_members = diffBD_members(ia,:);
    end
    if ~isempty(nonmembers)
        sw_ind = nonmembers(:,2) <  nonmembers(:,1);
        nonmembers(sw_ind,:) = [nonmembers(sw_ind,2) nonmembers(sw_ind,1)];
        ind = sub2ind([length( v.pyrInds ) length( v.pyrInds )], nonmembers(:,1), nonmembers(:,2) );
        [~, ia] = unique(ind);
        nonmembers = nonmembers(ia,:);
    end
    
    spkcount_theta = zeros(length( theta_ints_run ), length(v.pyrInds));
%     rateTheta_pyr = nan(length(v.pyrInds),1);
    % Number of spikes per ripple / theta cycle 
    for jp = 1:length(v.pyrInds)
        [theta_spk,int_indices_theta,~] = InIntervals( spikes.times{ v.pyrInds(jp) }, theta_ints_run );
        % Compute ripple rates
%         rateTheta_pyr(jp) = sum(theta_spk) / theta_len;
        % Spike count in each ripple / theta cycle
        int_indices_theta = int_indices_theta( int_indices_theta>0 );
        spkcount_theta(unique(int_indices_theta), jp) = [ cellfun(@(x) sum(x == int_indices_theta), num2cell( unique(int_indices_theta') )  ) ];
    end
    
    corrmat_theta = corr( spkcount_theta );
%     rateHor_theta = repmat( rateTheta_pyr, 1, length(rateTheta_pyr) );
%     rateVert_theta = repmat( rateTheta_pyr', length(rateTheta_pyr), 1 );

    % Gather cofiring between assembly members with the same birthdate
    if ~isempty(sameBD_members) && size(sameBD_members,1) > 1
        toremove = sameBD_members;
        
        thetacorr_tagtag = [ thetacorr_tagtag ; corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) ];
%         theta_tagtag_rates = [theta_tagtag_rates ; [rateHor_theta(sub2ind( size(rateHor_theta), toremove(:,1), toremove(:,2) ))...
%                                       rateVert_theta(sub2ind( size(rateVert_theta), toremove(:,1), toremove(:,2) ))  ] ];
        id_tagtag = [ id_tagtag ; repmat(tmp(1), size( toremove, 1), 1)];

        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
%         rateHor_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) = randLabel;
%         rateHor_theta( sub2ind( size(corrmat_theta), toremove(:,2), toremove(:,1) ) ) = randLabel;
%         rateVert_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) = randLabel;
%         rateVert_theta( sub2ind( size(corrmat_theta), toremove(:,2), toremove(:,1) ) ) = randLabel;
        corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) = randLabel;
        corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,2),  toremove(:,1) ) ) = randLabel;
    end
    % Gather cofiring between assembly members with different birthdates
    if ~isempty(diffBD_members) && size(diffBD_members,1) > 1
        toremove = diffBD_members;
 
        thetacorr_ntagMem = [ thetacorr_ntagMem ; corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) ];
%         theta_ntagMem_rates = [theta_ntagMem_rates ; [rateHor_theta(sub2ind( size(rateHor_theta), toremove(:,1), toremove(:,2) ))...
%                                       rateVert_theta(sub2ind( size(rateVert_theta), toremove(:,1), toremove(:,2) ))  ] ];
        id_ntagMem = [ id_ntagMem ; repmat(tmp(1), size( toremove, 1), 1)];

        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
%         rateHor_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) = randLabel;
%         rateHor_theta( sub2ind( size(corrmat_theta), toremove(:,2), toremove(:,1) ) ) = randLabel;
%         rateVert_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) = randLabel;
%         rateVert_theta( sub2ind( size(corrmat_theta), toremove(:,2), toremove(:,1) ) ) = randLabel;
        corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) = randLabel;
        corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,2),  toremove(:,1) ) ) = randLabel;
    end
    % Gather cofiring between held out neuron and assembly nonmebers
    if ~isempty(nonmembers) && size(nonmembers,1) > 1
        toremove = nonmembers;
 
        thetacorr_ntagNMem = [ thetacorr_ntagNMem ; corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) ];
%         ripple_ntagMem_rates = [ripple_ntagMem_rates ; [rateHor_ripple(sub2ind( size(rateHor_ripple), toremove(:,1), toremove(:,2) ))...
%                                       rateVert_ripple(sub2ind( size(rateVert_ripple), toremove(:,1), toremove(:,2) ))  ] ];
        id_ntagNMem = [ id_ntagNMem ; repmat(tmp(1), size( toremove, 1), 1)];

        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
%         rateHor_ripple( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
%         rateHor_ripple( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
%         rateVert_ripple( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
%         rateVert_ripple( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) = randLabel;
        corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,2),  toremove(:,1) ) ) = randLabel;
    end
    
%     corrmat_theta(1:length(v.pyrInds)+1:end) = randLabel;
%     rateHor_theta(1:length(v.pyrInds)+1:end) = randLabel;
%     rateVert_theta(1:length(v.pyrInds)+1:end) = randLabel;
%     
%     rateHor_theta = rateHor_theta(v.pyrOnlyInds_opto,:);       rateHor_theta = rateHor_theta(:);
%     rateVert_theta = rateVert_theta(v.pyrOnlyInds_opto,:);     rateVert_theta = rateVert_theta(:);
%     corrvec_theta = corrmat_theta(v.pyrOnlyInds_opto,:);       corrvec_theta = corrvec_theta(:);
%     
%     corrvec_theta( corrvec_theta == randLabel ) = [];
%     rateHor_theta(rateHor_theta == randLabel) = [];
%     rateVert_theta(rateVert_theta == randLabel) = [];
%     
%     thetacorr_ntagNMem = [thetacorr_ntagNMem ; corrvec_theta];
%     theta_ntagNMem_rates = [theta_ntagNMem_rates ; [rateHor_theta rateVert_theta]];
%     id_ntagNMem = [id_ntagNMem ; repmat(tmp(1), length(corrvec_theta), 1)];
end

%%

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( thetacorr_tagtag), [ 1:length(thetacorr_tagtag) ] ./ length(thetacorr_tagtag), '-b');
hold on
plot(sort( thetacorr_ntagMem), [ 1:length(thetacorr_ntagMem) ] ./ length(thetacorr_ntagMem), '-r')
plot(sort( thetacorr_ntagNMem), [ 1:length(thetacorr_ntagNMem) ] ./ length(thetacorr_ntagNMem), '-k')
xlim( [-0.1 0.16] )
xlabel('Theta cycle co-firing (\rho)', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)
legend('assembly members, same birthdate', 'assembly members, different birthdate', 'assembly non-members', 'location', 'best')

%%
figure
boxplot([thetacorr_ntagNMem ; thetacorr_ntagMem ; thetacorr_tagtag], [ones(size(thetacorr_ntagNMem)) ; 2*ones(size(thetacorr_ntagMem)) ; 3*ones(size(thetacorr_tagtag))],'notch','on', 'whisker',inf)
ylim([-0.02 0.06])

%%

