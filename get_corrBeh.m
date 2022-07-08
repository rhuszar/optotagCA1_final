%% Correlations between spatial firing maps of optotagged neurons!


% same birthdate versus different birthdate
load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
%load('C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac.mat')
% Do we wish to load data from local hard drive? 
isLocal = true;

% Do we consider the stem of the maze, or not
consider_stem = false;
% Indices of what's after the maze stem
post_stem_inds = 75:205;

ntag_wfoverlap = [];
pfcorr_ntag_all = [];
pfcorr_ntag_left = [];
pfcorr_ntag_right = [];
theta_ntag_rates = [];
pfcorr_ntag_dist = [];
pfcorr_ntag_shankdist = [];
thetacorr_ntag = [];            % pre / behavior / post
ntag_lRat = [];
id_ntag = [];
animal_ntag = [];

tagtag_wfoverlap = [];
pfcorr_tagtag_all = [];
pfcorr_tagtag_left = [];
pfcorr_tagtag_right = [];
theta_tagtag_rates = [];
pfcorr_tagtag_dist = [];
pfcorr_tagtag_shankdist = [];
thetacorr_tagtag = [];
tagtag_lRat = [];
id_tag = [];
animal_tag = [];


ttcorr_ntag = [];
id_ntag_cells = [];
ttcorr_tag = [];
id_tag_cells = [];

rateWake_tag = [];
rateTheta_tag = [];
rateWake_ntag = [];
rateTheta_ntag = [];

randLabel = 2000*rand;


for ip = 1:length(basepaths_beh)
    
    fprintf('%d/%d\n', ip, length(basepaths_beh))
    basepath = alterPath( basepaths_beh{ip}, isLocal );
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    mouse = [ tmp{1} '_' tmp{2} ];
    
    load(fullfile(basepath, [basename '.firingMapsAvg_v2.cellinfo.mat']))
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.thetaCycles.events.mat']))
    load(fullfile(basepath, [basename '.SleepState.states.mat']))
    cluQual = load(fullfile(basepath, 'PYR_clusterQuality.mat'));
    load( fullfile(basepath, 'wf1.mat') )
    load(fullfile(basepath, 'chanMap.mat'))
    v.chanMap = [xcoords ; ycoords];
    v.kcoords = kcoords;

    opto_id_log = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id_log = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id_log & pyr_id_log;
    pyr_id = find(pyr_id_log);
    
    % Get the L ratio
    [~, loc] = ismember(pyr_id, cluQual.pyr_UID);
    L_Rat = cluQual.L_Rat(loc); L_Rat = L_Rat(:);
    
    % Prepare theta intervals
    
    % Theta intervals in the PRE homecage
    theta_ints_pre = Theta.cycles(Theta.pre,:); 
    theta_ints_pre = theta_ints_pre( InIntervals( theta_ints_pre(:,1), SleepState.ints.WAKEstate ) & InIntervals( theta_ints_pre(:,2), SleepState.ints.WAKEstate ), : );
    theta_pre_len = sum( diff( theta_ints_pre' ) );
    % Theta intervals in the RUN period
    theta_ints_run = Theta.cycles(Theta.run,:); 
    theta_ints_run = theta_ints_run( InIntervals( theta_ints_run(:,1), SleepState.ints.WAKEstate ) & InIntervals( theta_ints_run(:,2), SleepState.ints.WAKEstate ), : );
    theta_run_len = sum( diff( theta_ints_run' ) );
    % Theta intervals in the POST homecage
    theta_ints_post = Theta.cycles(Theta.post,:);
    theta_ints_post = theta_ints_post( InIntervals( theta_ints_post(:,1), SleepState.ints.WAKEstate ) & InIntervals( theta_ints_post(:,2), SleepState.ints.WAKEstate ), : );
    theta_post_len = sum( diff( theta_ints_post' ) );
    
    % There are indices into pyr_id - i.e., if 5th index in pyr_id is optotagged, you want to store it 
    pyr_opto_id = arrayfun(@(x) find(x==pyr_id),  find(pyr_opto_id) );
    pyr_nopto_id = setdiff(1:length(pyr_id), pyr_opto_id);
    
    rateWAKE_pyr = cell_metrics.firingRate_WAKEstate(pyr_id)';
    % Below is code that can be used to extend physiology_process_v3.m
    
    dist_mat = nan(length(pyr_id), length(pyr_id));
    shankdist_mat = nan(length(pyr_id), length(pyr_id));
    for kp = 1:length(pyr_id)
        for jp = 1:length(pyr_id)
            a = [ v.chanMap(1, spikes.maxWaveformCh1(pyr_id(kp)) ) ;
            v.chanMap(2, spikes.maxWaveformCh1(pyr_id(kp)) )  ];
            b = [ v.chanMap(1, spikes.maxWaveformCh1(pyr_id(jp)) ) ;
            v.chanMap(2, spikes.maxWaveformCh1(pyr_id(jp)) )  ];
            dist_mat(kp,jp) = norm(a - b);
            shankdist_mat(kp, jp) = abs( v.kcoords( spikes.maxWaveformCh1(pyr_id(kp)) ) - v.kcoords( spikes.maxWaveformCh1(pyr_id(jp)) ) );
        end
    end
    
    % Quantify waveform overlap
    wf_overlap = nan(length(pyr_id), length(pyr_id));
    for kp = 1:length(pyr_id)
        for jp = 1:length(pyr_id)
            % Is kp neuron opto tagged
            if ismember( pyr_id(kp), find( opto_id_log & pyr_id_log ) ) && shankdist_mat(kp, jp)==0
                ind1 = find( wf.pyr_opto_id == pyr_id(kp) );
                if any( wf.pyr_opto_id == pyr_id(jp) )
                    ind2 = find( wf.pyr_opto_id == pyr_id(jp) );
                    a = cell2mat( wf.opto_spont_wf{ ind1 }( wf.opto_goodchan{ ind1 }) );
                    b = cell2mat( wf.opto_spont_wf{ ind2 }( wf.opto_goodchan{ ind2 }) );
                    wf_overlap(kp,jp) = norm( a-b ) / sum(wf.opto_goodchan{ ind1 });
                elseif any( wf.other_pyrid == pyr_id(jp) )
                    ind2 = find( wf.other_pyrid == pyr_id(jp) );
                    a = cell2mat( wf.opto_spont_wf{ ind1 }( wf.opto_goodchan{ ind1 }) );
                    b = cell2mat( wf.other_spont_wf{ ind2 }( wf.other_goodchan{ ind2 }) );
                    wf_overlap(kp,jp) = norm( a-b ) / sum(wf.other_goodchan{ ind2 });
                end
            end
        end
    end
    
    % Spike counts in theta cycles across experimental epochs
    spkcount_theta_pre = zeros(length( theta_ints_pre ), length(pyr_id));
    spkcount_theta_run = zeros(length( theta_ints_run ), length(pyr_id));
    spkcount_theta_post = zeros(length( theta_ints_post ), length(pyr_id));
    % Store ratemaps
    if consider_stem
        ratemaps_all = nan( 2*205, length(pyr_id));
        ratemaps_left = nan( 205, length(pyr_id));
        ratemaps_right = nan( 205, length(pyr_id));
    else
        ratemaps_all = nan( 2*length(post_stem_inds), length(pyr_id));
        ratemaps_left = nan( length(post_stem_inds), length(pyr_id));
        ratemaps_right = nan( length(post_stem_inds), length(pyr_id));
    end
    tt_corr_all =  nan( length(pyr_id), 1);
    
    thetaCycle_rates = nan(length(pyr_id), 3);
    % Number of spikes per theta cycle
    for jp = 1:length(pyr_id)
        % Spike counts in theta cycles
        [theta_pre_spk,int_indices_theta_pre,~] = InIntervals( spikes.times{ pyr_id(jp) }, theta_ints_pre );
        [theta_run_spk,int_indices_theta_run,~] = InIntervals( spikes.times{ pyr_id(jp) }, theta_ints_run );
        [theta_post_spk,int_indices_theta_post,~] = InIntervals( spikes.times{ pyr_id(jp) }, theta_ints_post );
        
        int_indices_theta_pre = int_indices_theta_pre( int_indices_theta_pre>0 );
        int_indices_theta_run = int_indices_theta_run( int_indices_theta_run>0 );
        int_indices_theta_post = int_indices_theta_post( int_indices_theta_post>0 );
        spkcount_theta_pre(unique(int_indices_theta_pre), jp) = [ cellfun(@(x) sum(x == int_indices_theta_pre), num2cell( unique(int_indices_theta_pre') )  ) ];
        spkcount_theta_run(unique(int_indices_theta_run), jp) = [ cellfun(@(x) sum(x == int_indices_theta_run), num2cell( unique(int_indices_theta_run') )  ) ];
        spkcount_theta_post(unique(int_indices_theta_post), jp) = [ cellfun(@(x) sum(x == int_indices_theta_post), num2cell( unique(int_indices_theta_post') )  ) ];
        % Ratemap correlations
        if consider_stem
            ratemaps_all(:,jp) = cell2mat( firingMaps.rateMaps{pyr_id(jp)} );
            ratemaps_left(:,jp) = firingMaps.rateMaps{pyr_id(jp)}{1};
            ratemaps_right(:,jp) = firingMaps.rateMaps{pyr_id(jp)}{2};
        else
            ratemaps_left(:,jp) = firingMaps.rateMaps{pyr_id(jp)}{1}(post_stem_inds);
            ratemaps_right(:,jp) = firingMaps.rateMaps{pyr_id(jp)}{2}(post_stem_inds);
            ratemaps_all(:,jp) = [ratemaps_left(:,jp) ; ratemaps_right(:,jp)];
        end
        cc = corr( [ratemaps_left(:,jp) ratemaps_right(:,jp) ] );
        tt_corr_all(jp) = cc(1,2);
        % Theta cycle rates
        thetaCycle_rates(jp,:) = [sum(theta_pre_spk) ./ theta_pre_len  sum(theta_run_spk) ./ theta_run_len sum(theta_post_spk) ./ theta_post_len];
    end
    
    % Store correlations between left/right trail type ratemaps
    % NOTE - this is one cc per unit, not per pair of units 
    ttcorr_tag = [ttcorr_tag ; tt_corr_all(pyr_opto_id)];
    ttcorr_ntag = [ttcorr_ntag ; tt_corr_all(pyr_nopto_id ) ];

    rateWake_tag = [rateWake_tag ; rateWAKE_pyr(pyr_opto_id)];
    rateWake_ntag = [rateWake_ntag ; rateWAKE_pyr(pyr_nopto_id)];
    rateTheta_tag = [rateTheta_tag ; thetaCycle_rates(pyr_opto_id,:)];
    rateTheta_ntag = [rateTheta_ntag ; thetaCycle_rates(pyr_nopto_id,:)];
    
    id_tag_cells = [ id_tag_cells ; repmat(tmp(1), length(pyr_opto_id), 1)];
    id_ntag_cells = [ id_ntag_cells ; repmat(tmp(1), length(pyr_nopto_id), 1)];
    
    % Ratemap correlation matrices
    corrmat_all = corr( ratemaps_all );
    corrmat_left = corr( ratemaps_left );
    corrmat_right = corr( ratemaps_right );
    
    rateHor_theta = repmat( thetaCycle_rates(:,2), 1, length(thetaCycle_rates(:,2)) );
    rateVert_theta = repmat( thetaCycle_rates(:,2)', length(thetaCycle_rates(:,2)), 1 );
    lHor = repmat( L_Rat, 1, length(L_Rat) );
    lVert = repmat( L_Rat', length(L_Rat), 1 );
    % Theta cycle correlation matrices
    corrmat_theta_pre = corr( spkcount_theta_pre );
    corrmat_theta_run = corr( spkcount_theta_run );
    corrmat_theta_post = corr( spkcount_theta_post );
    
    % Label duplicate pairs (e.g., if neuron 5 and 10 are tagged,  
    % pairs (5,10) and (10,5) will have the same correlation)
    label = 2000*rand;
    if length(pyr_opto_id)>1
        toremove = nchoosek( pyr_opto_id, 2 );
        
        pfcorr_tagtag_all = [ pfcorr_tagtag_all ; corrmat_all( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) ];
        pfcorr_tagtag_left = [ pfcorr_tagtag_left ; corrmat_left( sub2ind( size(corrmat_left), toremove(:,1), toremove(:,2) ) ) ];
        pfcorr_tagtag_right = [ pfcorr_tagtag_right ; corrmat_right( sub2ind( size(corrmat_right), toremove(:,1), toremove(:,2) ) ) ];
        pfcorr_tagtag_dist = [ pfcorr_tagtag_dist ; dist_mat( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) ];
        pfcorr_tagtag_shankdist = [ pfcorr_tagtag_shankdist ; shankdist_mat( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) ];
        % Here is where I store the rates of the tag tag
        theta_tagtag_rates = [theta_tagtag_rates ; [rateHor_theta(sub2ind( size(rateHor_theta), toremove(:,1), toremove(:,2) ))...
                                                      rateVert_theta(sub2ind( size(rateVert_theta), toremove(:,1), toremove(:,2) ))  ] ];
        thetacorr_tagtag = [thetacorr_tagtag ; [corrmat_theta_pre( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) ...
                                                  corrmat_theta_run( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) ...
                                                  corrmat_theta_post( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) )]];
        tagtag_lRat = [tagtag_lRat ; [lHor(sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ))...
                                              lVert(sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ))  ]];
        tagtag_wfoverlap = [tagtag_wfoverlap ; wf_overlap( sub2ind( size( corrmat_all ), toremove(:,1), toremove(:,2) ) ) ];                                  
        id_tag = [ id_tag ; repmat(tmp(1), size( toremove, 1), 1)];
        animal_tag = [animal_tag ; repmat({mouse}, size( toremove, 1), 1)];
        
        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
        lHor( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        lHor( sub2ind( size(corrmat_all), toremove(:,2), toremove(:,1) ) ) = label;
        lVert( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        lVert( sub2ind( size(corrmat_all), toremove(:,2), toremove(:,1) ) ) = label;
        
        wf_overlap( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        wf_overlap( sub2ind( size(corrmat_all), toremove(:,2), toremove(:,1) ) ) = label;
        
        rateHor_theta( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        rateHor_theta( sub2ind( size(corrmat_all), toremove(:,2), toremove(:,1) ) ) = label;
        rateVert_theta( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        rateVert_theta( sub2ind( size(corrmat_all), toremove(:,2), toremove(:,1) ) ) = label;
        dist_mat( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        dist_mat( sub2ind( size(corrmat_all), toremove(:,2), toremove(:,1) ) ) = label;
        shankdist_mat( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        shankdist_mat( sub2ind( size(corrmat_all), toremove(:,2), toremove(:,1) ) ) = label;
        corrmat_left( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        corrmat_left( sub2ind( size(corrmat_all), toremove(:,2),  toremove(:,1) ) ) = label;
        corrmat_right( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        corrmat_right( sub2ind( size(corrmat_all), toremove(:,2),  toremove(:,1) ) ) = label;
        corrmat_all( sub2ind( size(corrmat_all), toremove(:,1), toremove(:,2) ) ) = label;
        corrmat_all( sub2ind( size(corrmat_all), toremove(:,2),  toremove(:,1) ) ) = label;
        corrmat_theta_pre( sub2ind( size(corrmat_theta_pre), toremove(:,2),  toremove(:,1) ) ) = label;
        corrmat_theta_run( sub2ind( size(corrmat_theta_run), toremove(:,2),  toremove(:,1) ) ) = label;
        corrmat_theta_post( sub2ind( size(corrmat_theta_post), toremove(:,2),  toremove(:,1) ) ) = label;

    end
    % Label the diagonal 
    wf_overlap(1:length(pyr_id)+1:end) = label;
    dist_mat(1:length(pyr_id)+1:end) = label;
    shankdist_mat(1:length(pyr_id)+1:end) = label;
    corrmat_all(1:length(pyr_id)+1:end) = label;
    corrmat_left(1:length(pyr_id)+1:end) = label;
    corrmat_right(1:length(pyr_id)+1:end) = label;
    corrmat_theta_pre(1:length(pyr_id)+1:end) = label;
    corrmat_theta_run(1:length(pyr_id)+1:end) = label;
    corrmat_theta_post(1:length(pyr_id)+1:end) = label; 
    rateHor_theta(1:length(pyr_id)+1:end) = label;
    rateVert_theta(1:length(pyr_id)+1:end) = label;
    lHor(1:length(pyr_id)+1:end) = label;
    lVert(1:length(pyr_id)+1:end) = label;
    
    wf_overlap_vec = wf_overlap(pyr_opto_id,:);               wf_overlap_vec = wf_overlap_vec(:);
    corrvec_all = corrmat_all(pyr_opto_id,:);                 corrvec_all = corrvec_all(:);
    corrvec_left = corrmat_left(pyr_opto_id,:);               corrvec_left = corrvec_left(:); 
    corrvec_right = corrmat_right(pyr_opto_id,:);             corrvec_right = corrvec_right(:); 
    corrvec_theta_pre = corrmat_theta_pre(pyr_opto_id, :);    corrvec_theta_pre = corrvec_theta_pre(:);
    corrvec_theta_run = corrmat_theta_run(pyr_opto_id, :);    corrvec_theta_run = corrvec_theta_run(:);
    corrvec_theta_post = corrmat_theta_post(pyr_opto_id, :);  corrvec_theta_post = corrvec_theta_post(:);
    rateHor_theta = rateHor_theta(pyr_opto_id,:);               rateHor_theta = rateHor_theta(:);
    rateVert_theta = rateVert_theta(pyr_opto_id,:);             rateVert_theta = rateVert_theta(:);
    distvec = dist_mat(pyr_opto_id,:);                        distvec = distvec(:);
    shankdistvec = shankdist_mat(pyr_opto_id,:);              shankdistvec = shankdistvec(:);
    lHor = lHor(pyr_opto_id,:);                               lHor = lHor(:);
    lVert = lVert(pyr_opto_id,:);                             lVert = lVert(:);
    
    corrvec_all( corrvec_all == label ) = [];
    corrvec_left( corrvec_left == label ) = [];
    corrvec_right( corrvec_right == label ) = [];
    corrvec_theta_pre( corrvec_theta_pre == label ) = [];
    corrvec_theta_run( corrvec_theta_run == label ) = [];
    corrvec_theta_post( corrvec_theta_post == label ) = [];
    rateHor_theta(rateHor_theta == label) = [];
    rateVert_theta(rateVert_theta == label) = [];
    distvec( distvec == label ) = [];
    shankdistvec( shankdistvec == label ) = [];
    wf_overlap_vec(wf_overlap_vec == label) = [];
    lHor(lHor == label) = [];
    lVert(lVert == label) = [];
    
    ntag_wfoverlap = [ntag_wfoverlap ; wf_overlap_vec];
    pfcorr_ntag_dist = [ pfcorr_ntag_dist ; distvec ];
    pfcorr_ntag_shankdist = [ pfcorr_ntag_shankdist ; shankdistvec ];
    pfcorr_ntag_all = [ pfcorr_ntag_all ; corrvec_all ];
    pfcorr_ntag_left = [ pfcorr_ntag_left ; corrvec_left ];
    pfcorr_ntag_right = [ pfcorr_ntag_right ; corrvec_right ];
    theta_ntag_rates = [theta_ntag_rates ; [rateHor_theta rateVert_theta]];
    thetacorr_ntag = [thetacorr_ntag ; [corrvec_theta_pre corrvec_theta_run corrvec_theta_post]];
    id_ntag = [ id_ntag ; repmat(tmp(1), length(corrvec_all), 1)];
    ntag_lRat = [ntag_lRat ; [lHor lVert]];
    animal_ntag = [animal_ntag ; repmat({mouse}, length(corrvec_all), 1)];
    

end

save('pfcorr_v6.mat', 'pfcorr_ntag_all', 'pfcorr_ntag_left', 'pfcorr_ntag_right', 'pfcorr_ntag_dist', 'pfcorr_ntag_shankdist', ...
             'theta_ntag_rates', 'thetacorr_ntag', 'id_ntag', 'pfcorr_tagtag_all', 'pfcorr_tagtag_left', 'pfcorr_tagtag_right', ...
          'pfcorr_tagtag_dist', 'pfcorr_tagtag_shankdist', 'theta_tagtag_rates', 'thetacorr_tagtag', 'id_tag', ...
       'ttcorr_tag', 'ttcorr_ntag', 'rateWake_tag', 'rateWake_ntag', 'rateTheta_tag', 'rateTheta_ntag',...
     'id_tag_cells', 'id_ntag_cells', 'ntag_lRat', 'tagtag_lRat', 'animal_ntag', 'animal_tag', 'ntag_wfoverlap', 'tagtag_wfoverlap')

 
 %% Distance between place fields - 

% same birthdate versus different birthdate
% load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
load('C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac.mat')
% Do we wish to load data from local hard drive? 
isLocal = true;

% Do we consider the stem of the maze, or not
consider_stem = false;
% Indices of what's after the maze stem
post_stem_inds = 75:205;

pfDist_ntag_all = [];
pfDist_ntag_dist = [];
pfDist_ntag_shankdist = [];
pfDist_ntag_left = [];
pfDist_ntag_right = [];
theta_ntag_rates = [];
ntag_lRat = [];
id_ntag = [];
animal_ntag = [];

pfDist_tagtag_all = [];
pfDist_tagtag_left = [];
pfDist_tagtag_right = [];
pfDist_tagtag_dist = [];
pfDist_tagtag_shankdist = [];
theta_tagtag_rates = [];
tagtag_lRat = [];
id_tag = [];
animal_tag = [];

randLabel = 2000*rand;


for ip = 1:length(basepaths_beh)
    
    fprintf('%d/%d\n', ip, length(basepaths_beh))
    basepath = alterPath( basepaths_beh{ip}, isLocal );
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    mouse = [ tmp{1} '_' tmp{2} ];
    
%     pathparts = strsplit(basepath,filesep);
%     basepath = fullfile(datpath, pathparts{end-2}, pathparts{end-1}, pathparts{end});
    
    load(fullfile(basepath, [basename '.firingMapsAvg_v2.cellinfo.mat']))
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath, [basename '.SleepState.states.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.thetaCycles.events.mat']))
    cluQual = load(fullfile(basepath, 'PYR_clusterQuality.mat'));
    load(fullfile(basepath, 'chanMap.mat'))
    v.chanMap = [xcoords ; ycoords];
    v.kcoords = kcoords;

    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id & pyr_id;
    pyr_id = find(pyr_id);
    
    % Get the L ratio
    [~, loc] = ismember(pyr_id, cluQual.pyr_UID);
    L_Rat = cluQual.L_Rat(loc); L_Rat = L_Rat(:);
    
    % Prepare theta intervals

    % Theta intervals in the RUN period
    theta_ints_run = Theta.cycles(Theta.run,:); 
    theta_ints_run = theta_ints_run( InIntervals( theta_ints_run(:,1), SleepState.ints.WAKEstate ) & InIntervals( theta_ints_run(:,2), SleepState.ints.WAKEstate ), : );
    theta_run_len = sum( diff( theta_ints_run' ) );

    % There are indices into pyr_id - i.e., if 5th index in pyr_id is optotagged, you want to store it 
    pyr_opto_id = arrayfun(@(x) find(x==pyr_id),  find(pyr_opto_id) );
    pyr_nopto_id = setdiff(1:length(pyr_id), pyr_opto_id);
    
    % Below is code that can be used to extend physiology_process_v3.m
    
    dist_mat = nan(length(pyr_id), length(pyr_id));
    shankdist_mat = nan(length(pyr_id), length(pyr_id));
    for kp = 1:length(pyr_id)
        for jp = 1:length(pyr_id)
            a = [ v.chanMap(1, spikes.maxWaveformCh1(pyr_id(kp)) ) ;
            v.chanMap(2, spikes.maxWaveformCh1(pyr_id(kp)) )  ];
            b = [ v.chanMap(1, spikes.maxWaveformCh1(pyr_id(jp)) ) ;
            v.chanMap(2, spikes.maxWaveformCh1(pyr_id(jp)) )  ];
            dist_mat(kp,jp) = norm(a - b);
            shankdist_mat(kp, jp) = abs( v.kcoords( spikes.maxWaveformCh1(pyr_id(kp)) ) - v.kcoords( spikes.maxWaveformCh1(pyr_id(jp)) ) );
        end
    end
    
    fields_left = {}; fields_right = {};
    % Extract the place fields
    for kp = 1:length(pyr_id)
        [ fields_left(kp) ] = returnPlaceFields( firingMaps.rateMaps{pyr_id(kp)}{1}(post_stem_inds), firingMaps.pos_info{pyr_id(kp),1}(post_stem_inds));
        [ fields_right(kp) ] = returnPlaceFields(firingMaps.rateMaps{pyr_id(kp)}{2}(post_stem_inds), firingMaps.pos_info{pyr_id(kp),2}(post_stem_inds));
    end
    
    % Extract place field distances
    pf_dist_left = nan(length(fields_left), length(fields_left));
    pf_dist_right = nan(length(fields_right), length(fields_left));
    pf_dist_all = nan(length(fields_right), length(fields_left));
    for kp = 1:length(fields_left)
        for jp = 1:length(fields_left)
            pf_dist_left(kp,jp) = pf_dist(fields_left{kp}, fields_left{jp});
            pf_dist_right(kp,jp) = pf_dist(fields_right{kp}, fields_right{jp});
            pf_dist_all(kp,jp) = min([ pf_dist_left(kp,jp), pf_dist_right(kp,jp) ]); 
        end
    end
    
    % Spike counts in theta cycles across experimental epochs
    spkcount_theta_run = zeros(length( theta_ints_run ), length(pyr_id));

    thetaCycle_rates = nan(length(pyr_id), 1);
    % Theta cycle related firing rate
    for jp = 1:length(pyr_id)
        % Spike counts in theta cycles
        [theta_run_spk,int_indices_theta_run,~] = InIntervals( spikes.times{ pyr_id(jp) }, theta_ints_run );
        % Theta cycle rates
        thetaCycle_rates(jp) = sum(theta_run_spk) ./ theta_run_len ;
    end
    
    rateHor_theta = repmat( thetaCycle_rates, 1, length(thetaCycle_rates) );
    rateVert_theta = repmat( thetaCycle_rates', length(thetaCycle_rates), 1 );
    lHor = repmat( L_Rat, 1, length(L_Rat) );
    lVert = repmat( L_Rat', length(L_Rat), 1 );

    % Label duplicate pairs (e.g., if neuron 5 and 10 are tagged,  
    % pairs (5,10) and (10,5) will have the same correlation)
    label = 2000*rand;
    if length(pyr_opto_id)>1
        toremove = nchoosek( pyr_opto_id, 2 );

        pfDist_tagtag_left = [ pfDist_tagtag_left ; pf_dist_left( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) ];
        pfDist_tagtag_right = [ pfDist_tagtag_right ; pf_dist_right( sub2ind( size(pf_dist_right), toremove(:,1), toremove(:,2) ) ) ];
        pfDist_tagtag_all = [ pfDist_tagtag_all ; pf_dist_all( sub2ind( size(pf_dist_right), toremove(:,1), toremove(:,2) ) ) ];
        pfDist_tagtag_dist = [ pfDist_tagtag_dist ; dist_mat( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) ];
        pfDist_tagtag_shankdist = [ pfDist_tagtag_shankdist ; shankdist_mat( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) ];
        % Here is where I store the rates of the tag tag
        theta_tagtag_rates = [theta_tagtag_rates ; [rateHor_theta(sub2ind( size(rateHor_theta), toremove(:,1), toremove(:,2) ))...
                                                      rateVert_theta(sub2ind( size(rateVert_theta), toremove(:,1), toremove(:,2) ))  ] ];
        tagtag_lRat = [tagtag_lRat ; [lHor(sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ))...
                                              lVert(sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ))  ]];
                                          
        id_tag = [ id_tag ; repmat(tmp(1), size( toremove, 1), 1)];
        animal_tag = [animal_tag ; repmat({mouse}, size( toremove, 1), 1)];
        
        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
        lHor( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        lHor( sub2ind( size(pf_dist_left), toremove(:,2), toremove(:,1) ) ) = label;
        lVert( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        lVert( sub2ind( size(pf_dist_left), toremove(:,2), toremove(:,1) ) ) = label;
        
        rateHor_theta( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        rateHor_theta( sub2ind( size(pf_dist_left), toremove(:,2), toremove(:,1) ) ) = label;
        rateVert_theta( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        rateVert_theta( sub2ind( size(pf_dist_left), toremove(:,2), toremove(:,1) ) ) = label;
        dist_mat( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        dist_mat( sub2ind( size(pf_dist_left), toremove(:,2), toremove(:,1) ) ) = label;
        shankdist_mat( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        shankdist_mat( sub2ind( size(pf_dist_left), toremove(:,2), toremove(:,1) ) ) = label;
        pf_dist_left( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        pf_dist_left( sub2ind( size(pf_dist_left), toremove(:,2),  toremove(:,1) ) ) = label;
        pf_dist_right( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        pf_dist_right( sub2ind( size(pf_dist_left), toremove(:,2),  toremove(:,1) ) ) = label;
        pf_dist_all( sub2ind( size(pf_dist_left), toremove(:,1), toremove(:,2) ) ) = label;
        pf_dist_all( sub2ind( size(pf_dist_left), toremove(:,2),  toremove(:,1) ) ) = label;

    end
    % Label the diagonal 
    dist_mat(1:length(pyr_id)+1:end) = label;
    shankdist_mat(1:length(pyr_id)+1:end) = label;
    pf_dist_left(1:length(pyr_id)+1:end) = label;
    pf_dist_right(1:length(pyr_id)+1:end) = label;
    pf_dist_all(1:length(pyr_id)+1:end) = label;
    rateHor_theta(1:length(pyr_id)+1:end) = label;
    rateVert_theta(1:length(pyr_id)+1:end) = label;
    lHor(1:length(pyr_id)+1:end) = label;
    lVert(1:length(pyr_id)+1:end) = label;
    
    pfdistvec_left = pf_dist_left(pyr_opto_id,:);               pfdistvec_left = pfdistvec_left(:); 
    pfdistvec_right = pf_dist_right(pyr_opto_id,:);             pfdistvec_right = pfdistvec_right(:);
    pfdistvec_all = pf_dist_all(pyr_opto_id,:);                 pfdistvec_all = pfdistvec_all(:); 
    rateHor_theta = rateHor_theta(pyr_opto_id,:);               rateHor_theta = rateHor_theta(:);
    rateVert_theta = rateVert_theta(pyr_opto_id,:);             rateVert_theta = rateVert_theta(:);
    distvec = dist_mat(pyr_opto_id,:);                          distvec = distvec(:);
    shankdistvec = shankdist_mat(pyr_opto_id,:);                shankdistvec = shankdistvec(:);
    lHor = lHor(pyr_opto_id,:);                                 lHor = lHor(:);
    lVert = lVert(pyr_opto_id,:);                               lVert = lVert(:);
    
    pfdistvec_left( pfdistvec_left == label ) = [];
    pfdistvec_right( pfdistvec_right == label ) = [];
    pfdistvec_all( pfdistvec_all == label ) = [];
    rateHor_theta(rateHor_theta == label) = [];
    rateVert_theta(rateVert_theta == label) = [];
    distvec( distvec == label ) = [];
    shankdistvec( shankdistvec == label ) = [];
    lHor(lHor == randLabel) = [];
    lVert(lVert == randLabel) = [];
    
    pfDist_ntag_dist = [ pfDist_ntag_dist ; distvec ];
    pfDist_ntag_shankdist = [ pfDist_ntag_shankdist ; shankdistvec ];
    pfDist_ntag_left = [ pfDist_ntag_left ; pfdistvec_left ];
    pfDist_ntag_right = [ pfDist_ntag_right ; pfdistvec_right ];
    pfDist_ntag_all = [ pfDist_ntag_all ; pfdistvec_all ];
    theta_ntag_rates = [theta_ntag_rates ; [rateHor_theta rateVert_theta]];
    id_ntag = [ id_ntag ; repmat(tmp(1), length(pfdistvec_right), 1)];
    ntag_lRat = [ntag_lRat ; [lHor lVert]];
    animal_ntag = [animal_ntag ; repmat({mouse}, length(pfdistvec_right), 1)];

end

save('pfdist_v1.mat', 'pfDist_ntag_dist', 'pfDist_ntag_shankdist', 'pfDist_ntag_left', 'pfDist_ntag_right', 'pfDist_ntag_all', ...
             'theta_ntag_rates', 'id_ntag', 'ntag_lRat', 'animal_ntag', ...
          'pfDist_tagtag_left', 'pfDist_tagtag_right', 'pfDist_tagtag_all', 'pfDist_tagtag_dist', 'pfDist_tagtag_shankdist', ...
       'theta_tagtag_rates', 'tagtag_lRat', 'id_tag', 'animal_tag')

%%


