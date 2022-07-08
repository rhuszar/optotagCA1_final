

saveMat = true;
% cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV3_OUT')
%cd("C:\Users\Roman\Google Drive\buzz_workspace\optotag\DATA\runAssembly_trackV3_OUT")
cd('C:\Users\Roman\Documents\DATA\runAssembly_trackV3_OUT_rev')
% Interested in everything
%load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac.mat")
basepaths = basepaths_all;
% basepaths = basepaths_beh;

% Interested in behavior
% load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_beh.mat')
% basepaths = basepaths_beh;

fils = dir;
fils = {fils.name}; fils(1:2) = [];
                   
nsd = 2;

% Distance estimates
ntag_dist = [];
ntag_shankdist = [];
% Correlations
ripcorr_ntag = [];
% Ripple rates
ripple_ntag_rates = [];
id_ntag = [];
animal_ntag = [];
basenames_ntag = [];

% Distance estimates
tagtag_dist = [];
tagtag_shankdist = [];
% Correlations
ripcorr_tagtag = [];
% Ripple rates
ripple_tagtag_rates = [];
id_tag = [];
animal_tag = [];
basenames_tag = [];

randLabel = 2000*rand;

% Each file is a single session 
% Single out assemblies with at least one high contributing opto tagged
% neuron
exclude_animal = {'e15_13f1', 'e16_3m2'};
for ii = 1:length( basepaths )
    
    basepath = alterPath( basepaths{ii}, true );
    basename = bz_BasenameFromBasepath( basepath );
    
    disp(basename)
    fils_id = find( cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils) );
    if ~any( fils_id )
        continue
    end
    
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    
    % Load assembly output
    v = load( fils{fils_id} );
    % What optotag session are we in?
    tmp = strsplit( fils{fils_id}, '_');    
    
    pyrInds_opto = arrayfun(@(x) find( x == v.pyrInds) ,  v.pyrInds_opto);
    % Number of "high contributing" neurons (weights > 2SD)
    HC_indicator = zscore(v.assemblyTemplates_pre) > repmat( nsd, size(v.assemblyTemplates_pre));
    % Find assemblies with at least one tagged neuron as a high contributor
    taggedAssembly_id = find( sum( HC_indicator( pyrInds_opto, : ),1 ))';
    
    if isempty(taggedAssembly_id)
        continue
    else
        %disp('We have something')
    end
    
    
    load(fullfile(basepath, [basename '.ripples.events.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.pulses.events.mat']))
    load( fullfile(basepath, [basename '.session.mat']) );

    load(fullfile(basepath, 'chanMap.mat'))
    v.chanMap = [xcoords ; ycoords];
    v.kcoords = kcoords;
    
%     if exist( fullfile(basepath, [basename '.Behavior.mat']) )
%         load(fullfile(basepath, [basename '.Behavior.mat']))
%         beh_ints = [behavior.timestamps(1) behavior.timestamps(end)];
%     else
%         beh_ints = [];
%     end
%     ripple_ints = SubtractIntervals( ripples.timestamps,  [ pulses.intsPeriods ; beh_ints ] );
%     ripple_time = sum( diff( ripple_ints' ) );

    reclen = session.timeSeries.adc.nSamples / session.timeSeries.adc.sr;
    % Across all relevant sessions, we calculate the ripple / theta
    % cofiring in the homecage 
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
    ripple_ints = ripples.timestamps;
    ripple_ints( ripple_ints(:,1) > discard,: ) = [];
    rip_in_pulse = InIntervals( ripple_ints(:,1), [ pulses.intsPeriods ; beh_ts ] ) | InIntervals( ripple_ints(:,2), [ pulses.intsPeriods ; beh_ts ] );
    ripple_ints(rip_in_pulse,:) = [];
    ripple_time = sum( diff( ripple_ints' ) );
    
    % Obtain all pairs of high contributing neurons
    hc_pairs = [];
    hc_pyrs = [];
    for jp = 1:length(taggedAssembly_id)
        hcp = find( HC_indicator(:,taggedAssembly_id(jp)) );
        hc_pyrs = [hc_pyrs ; hcp];
        if length(hcp) > 1
            hc_pairs = [ hc_pairs ; nchoosek(  hcp', 2)];
        end
    end
    % Make sure we only keep unique pairs - express these as indices into an upper triangular matrix
    if ~isempty(hc_pairs)
        sw_ind = hc_pairs(:,2) <  hc_pairs(:,1);
        hc_pairs(sw_ind,:) = [hc_pairs(sw_ind,2) hc_pairs(sw_ind,1)];
        ind = sub2ind([length( v.pyrInds ) length( v.pyrInds )], hc_pairs(:,1), hc_pairs(:,2) );
        [~, ia] = unique(ind);
        hc_pairs = hc_pairs(ia,:);
    end
    hc_pyrs = unique(hc_pyrs);
    
    dist_mat = nan(length(v.pyrInds), length(v.pyrInds));
    shankdist_mat = nan(length(v.pyrInds), length(v.pyrInds));
    for kp = 1:length(v.pyrInds)
        for jp = 1:length(v.pyrInds)
            a = [ v.chanMap(1, spikes.maxWaveformCh1(v.pyrInds(kp)) ) ;
            v.chanMap(2, spikes.maxWaveformCh1(v.pyrInds(kp)) )  ];
            b = [ v.chanMap(1, spikes.maxWaveformCh1(v.pyrInds(jp)) ) ;
            v.chanMap(2, spikes.maxWaveformCh1(v.pyrInds(jp)) )  ];
            dist_mat(kp,jp) = norm(a - b);
            shankdist_mat(kp, jp) = abs( v.kcoords( spikes.maxWaveformCh1(v.pyrInds(kp)) ) - v.kcoords( spikes.maxWaveformCh1(v.pyrInds(jp)) ) );
        end
    end
    
    spkcount_ripple = zeros(length( ripple_ints ), length(v.pyrInds));
    rateRipple_pyr = nan(length(v.pyrInds),1);
    % Number of spikes per ripple / theta cycle 
    for jp = 1:length(v.pyrInds)
        [rip_spk,int_indices_ripple,~] = InIntervals( spikes.times{ v.pyrInds(jp) }, ripple_ints );
        % Compute ripple rates
        rateRipple_pyr(jp) = sum(rip_spk) / ripple_time;
        % Spike count in each ripple / theta cycle
        int_indices_ripple = int_indices_ripple( int_indices_ripple>0 );
        spkcount_ripple(unique(int_indices_ripple), jp) = [ cellfun(@(x) sum(x == int_indices_ripple), num2cell( unique(int_indices_ripple') )  ) ];
    end
    
    corrmat_ripple = corr( spkcount_ripple );
    rateHor_rip = repmat( rateRipple_pyr, 1, length(rateRipple_pyr) );
    rateVert_rip = repmat( rateRipple_pyr', length(rateRipple_pyr), 1 );
    
    if ~isempty(hc_pairs) && size(hc_pairs,1) > 1
        toremove = hc_pairs;
        
        ripcorr_tagtag = [ ripcorr_tagtag ; corrmat_ripple( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) ];
        ripple_tagtag_rates = [ripple_tagtag_rates ; [rateHor_rip(sub2ind( size(rateHor_rip), toremove(:,1), toremove(:,2) ))...
                                      rateVert_rip(sub2ind( size(rateVert_rip), toremove(:,1), toremove(:,2) ))  ] ];
        id_tag = [ id_tag ; repmat(tmp(1), size( toremove, 1), 1)];
        animal_tag = [ animal_tag ; repmat({[ tmp{1} '_' tmp{2} ]}, size( toremove, 1), 1)];
        tagtag_dist = [ tagtag_dist ; dist_mat( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) ];
        tagtag_shankdist = [ tagtag_shankdist ; shankdist_mat( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) ];
        basenames_tag = [basenames_tag ; repmat( {basename}, size( toremove, 1), 1) ];
        
        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
        
        rateHor_rip( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        rateHor_rip( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        rateVert_rip( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        rateVert_rip( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        
        dist_mat( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        dist_mat( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        shankdist_mat( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        shankdist_mat( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        corrmat_ripple( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        corrmat_ripple( sub2ind( size(corrmat_ripple), toremove(:,2),  toremove(:,1) ) ) = randLabel;
    end
    
    dist_mat(1:length(v.pyrInds)+1:end) = randLabel;
    shankdist_mat(1:length(v.pyrInds)+1:end) = randLabel;
    corrmat_ripple(1:length(v.pyrInds)+1:end) = randLabel;
    rateHor_rip(1:length(v.pyrInds)+1:end) = randLabel;
    rateVert_rip(1:length(v.pyrInds)+1:end) = randLabel;
    
    rateHor_rip = rateHor_rip(hc_pyrs,:);           rateHor_rip = rateHor_rip(:);
    rateVert_rip = rateVert_rip(hc_pyrs,:);         rateVert_rip = rateVert_rip(:);
    corrvec_ripple = corrmat_ripple(hc_pyrs,:);     corrvec_ripple = corrvec_ripple(:);
    distvec = dist_mat(hc_pyrs,:);                  distvec = distvec(:);
    shankdistvec = shankdist_mat(hc_pyrs,:);         shankdistvec = shankdistvec(:);
    
    corrvec_ripple( corrvec_ripple == randLabel ) = [];
    distvec( distvec == randLabel ) = [];
    shankdistvec( shankdistvec == randLabel ) = [];
    rateHor_rip(rateHor_rip == randLabel) = [];
    rateVert_rip(rateVert_rip == randLabel) = [];
    
    ripple_ntag_rates = [ripple_ntag_rates ; [rateHor_rip rateVert_rip]];
    ntag_dist = [ ntag_dist ; distvec ];
   ntag_shankdist = [ ntag_shankdist ; shankdistvec ];
    ripcorr_ntag = [ ripcorr_ntag ; corrvec_ripple ];
    id_ntag = [ id_ntag ; repmat(tmp(1), length(corrvec_ripple), 1)];
    animal_ntag = [animal_ntag ; repmat({[ tmp{1} '_' tmp{2} ]}, length(corrvec_ripple), 1)];
    basenames_ntag = [basenames_ntag ; repmat( {basename}, length(corrvec_ripple), 1) ];
    
    
end

%% Ripple cofiring of assembly members / nonmembers

figure
set(gcf,'Position',[440 378 396 420])
plot(sort( ripcorr_ntag), [ 1:length(ripcorr_ntag) ] ./ length(ripcorr_ntag), '-r');
hold on
plot(sort( ripcorr_tagtag), [ 1:length(ripcorr_tagtag) ] ./ length(ripcorr_tagtag), '-b')
xlim( [-0.1 0.3] )
xlabel('Ripple cofiring (\rho)', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)
legend('Assembly non-members', 'Assembly members', 'Location', 'Best')

figure
boxplot([ripcorr_ntag ; ripcorr_tagtag], [ones(size(ripcorr_ntag)) ; 2*ones(size(ripcorr_tagtag))],'notch','on', 'whisker',inf,...
        'labels',{'Assembly non-members','Assembly members'})
ylim([-.02 0.11])
%% Ripple cofiring of assembly members / nonmembers, split by birthdate of tagged assembly members

e13_tag = cellfun(@(x) strcmp(x, 'e13'), id_tag); e13_ntag = cellfun(@(x) strcmp(x, 'e13'), id_ntag);
e14_tag = cellfun(@(x) strcmp(x, 'e14'), id_tag); e14_ntag = cellfun(@(x) strcmp(x, 'e14'), id_ntag);
e15_tag = cellfun(@(x) strcmp(x, 'e15'), id_tag); e15_ntag = cellfun(@(x) strcmp(x, 'e15'), id_ntag);
e16_tag = cellfun(@(x) strcmp(x, 'e16'), id_tag); e16_ntag = cellfun(@(x) strcmp(x, 'e16'), id_ntag);

tagtag_mu = [nanmean(ripcorr_tagtag(e13_tag)) ; nanmean(ripcorr_tagtag(e14_tag)) ; nanmean(ripcorr_tagtag(e15_tag)) ; nanmean(ripcorr_tagtag(e16_tag)) ; ];
tagtag_sem = [sem(ripcorr_tagtag(e13_tag)) ; sem(ripcorr_tagtag(e14_tag)) ; sem(ripcorr_tagtag(e15_tag)) ; sem(ripcorr_tagtag(e16_tag)); ];

ntag_mu = [nanmean(ripcorr_ntag(e13_ntag)) ; nanmean(ripcorr_ntag(e14_ntag)) ; nanmean(ripcorr_ntag(e15_ntag)) ; nanmean(ripcorr_ntag(e16_ntag)) ; ];
ntag_sem = [sem(ripcorr_ntag(e13_ntag)) ; sem(ripcorr_ntag(e14_ntag)) ; sem(ripcorr_ntag(e15_ntag)) ; sem(ripcorr_ntag(e16_ntag)); ];

errorbar(1:4,tagtag_mu,tagtag_sem, '-s', 'MarkerSize',10,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','CapSize',18')
hold on
errorbar(1:4,ntag_mu,ntag_sem,'-s', 'MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',18)

set(gcf,'Position',[440 378 404 420])
legend('Assembly members', 'Assembly non-members', 'Location', 'Best')

xticks([1 2 3 4])
xticklabels({'E13.5', 'E14.5', 'E15.5', 'E16.5'})
ylabel('Ripple cofiring (\rho)', 'fontsize', 20)
xlabel('Birthdate','fontsize', 20)

