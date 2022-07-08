%%

% same birthdate versus different birthdate

load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
% load('C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac.mat')

% Do we wish to load data from local hard drive? 
isLocal = true;

% Distance estimates
ntag_dist = [];
ntag_shankdist = [];
ntag_lRat = [];
% Correlations
ripcorr_ntag = [];
thetacorr_ntag = [];
% NREM / WAKE / ripple rates
theta_ntag_rates = [];
ripple_ntag_rates = [];
id_ntag = [];
animal_ntag = [];
ntag_wfoverlap = [];
converge_ntag = [];

% Distance estimates
tagtag_dist = [];
tagtag_shankdist = [];
tagtag_lRat = [];
% Correlations
ripcorr_tagtag = [];
thetacorr_tagtag = [];
% NREM / WAKE / ripple rates
theta_tagtag_rates = [];
ripple_tagtag_rates = [];
id_tag = [];
animal_tag = [];
tagtag_wfoverlap = [];
converge_tagtag = [];

id_session = [];
ripcorr_tagtag_session = [];
ripcorr_ntag_session = [];
propotion_tagged = [];


% Store unlabelled as nolab
converge_nolab = [];
nolab_wfoverlap = [];
ripple_nolab_rates = [];
theta_nolab_rates = [];
nolab_lRat = [];
nolab_dist = [];
nolab_shankdist = [];
ripcorr_nolab = [];
thetacorr_nolab = [];
id_nolab = [];
animal_nolab = [];

ripcorr_nolab_session = [];

randLabel = 2000*rand;

for ip = 1:length(basepaths_all)
    flag = false;
    fprintf('%d/%d\n', ip, length(basepaths_all))
    basepath = alterPath( basepaths_all{ip}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    mouse = [ tmp{1} '_' tmp{2} ];
    
    load(fullfile(basepath, [basename '.ripples.events.mat']))
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.pulses.events.mat']))
    load(fullfile(basepath, [basename '.thetaCycles.events.mat']))
    load(fullfile(basepath, [basename '.SleepState.states.mat']))
    load(fullfile(basepath, [basename '.mono_res.cellinfo.mat']))
    cluQual = load(fullfile(basepath, 'PYR_clusterQuality.mat'));
    load( fullfile(basepath, 'wf1.mat') )
    if exist( fullfile(basepath, [basename '.Behavior.mat']) )
        load(fullfile(basepath, [basename '.Behavior.mat']))
        beh_ints = [behavior.timestamps(1) behavior.timestamps(end)];
    else
        beh_ints = [];
    end
    
    load(fullfile(basepath, 'chanMap.mat'))
    v.chanMap = [xcoords ; ycoords];
    v.kcoords = kcoords;
    
    % Prepare ripple intervals
    ripple_ints = SubtractIntervals( ripples.timestamps,  [ pulses.intsPeriods ; beh_ints ] );
    ripple_time = sum( diff( ripple_ints' ) );
    % Prepare theta intervals
    theta_ints = SubtractIntervals( Theta.cycles, [ pulses.intsPeriods ; beh_ints ] );
    theta_ints = theta_ints( InIntervals( theta_ints(:,1), SleepState.ints.WAKEstate ) & InIntervals( theta_ints(:,2), SleepState.ints.WAKEstate ), : );
    theta_time = sum( diff( theta_ints' ) );
    
    opto_id_log = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id_log = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id_log & pyr_id_log;
    
    pyr_id = find(pyr_id_log);
    pyr_opto_id = arrayfun(@(x) find(x==pyr_id),  find(pyr_opto_id) );
    pyr_nopto_id = setdiff(1:length(pyr_id), pyr_opto_id);
    
    % Get the L ratio
    [~, loc] = ismember(pyr_id, cluQual.pyr_UID);
    L_Rat = cluQual.L_Rat(loc); L_Rat = L_Rat(:);
    
    % Obtain shank / channel distances
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
    
    convergence_mat = nan(length(pyr_id),length(pyr_id));
    % Number of spikes per ripple / theta cycle 
    for jp = 1:length(pyr_id)
        for kp = 1:length(pyr_id)
            n1 = unique( mono_res.pyr2int( mono_res.pyr2int(:,1) == pyr_id(jp), 2) );
            n2 = unique( mono_res.pyr2int( mono_res.pyr2int(:,1) == pyr_id(kp), 2) );
            convergence_mat(jp,kp) = length( intersect(n1, n2) ) / length( unique( [n1 ; n2] ) );
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
    
    spkcount_ripple = zeros(length( ripple_ints ), length(pyr_id));
    spkcount_theta = zeros(length( theta_ints ), length(pyr_id));
    rateRipple_pyr = nan(length(pyr_id),1);
    rateTheta_pyr = nan(length(pyr_id),1);
    % Number of spikes per ripple / theta cycle 
    for jp = 1:length(pyr_id)
        [rip_spk,int_indices_ripple,~] = InIntervals( spikes.times{ pyr_id(jp) }, ripple_ints );
        [theta_spk,int_indices_theta,~] = InIntervals( spikes.times{ pyr_id(jp) }, theta_ints );
        % Compute ripple / theta rates
        rateRipple_pyr(jp) = sum(rip_spk) / ripple_time;
        rateTheta_pyr(jp) = sum(theta_spk) / theta_time;
        % Spike count in each ripple / theta cycle
        int_indices_ripple = int_indices_ripple( int_indices_ripple>0 );
        int_indices_theta = int_indices_theta( int_indices_theta>0 );
        spkcount_ripple(unique(int_indices_ripple), jp) = [ cellfun(@(x) sum(x == int_indices_ripple), num2cell( unique(int_indices_ripple') )  ) ];
        spkcount_theta(unique(int_indices_theta), jp) = [ cellfun(@(x) sum(x == int_indices_theta), num2cell( unique(int_indices_theta') )  ) ];
    end
    corrmat_ripple = corr( spkcount_ripple );
    corrmat_theta = corr( spkcount_theta );
    rateHor_theta = repmat( rateTheta_pyr, 1, length(rateTheta_pyr) );
    rateVert_theta = repmat( rateTheta_pyr', length(rateTheta_pyr), 1 );
    lHor = repmat( L_Rat, 1, length(L_Rat) );
    lVert = repmat( L_Rat', length(L_Rat), 1 );
    rateHor_rip = repmat( rateRipple_pyr, 1, length(rateRipple_pyr) );
    rateVert_rip = repmat( rateRipple_pyr', length(rateRipple_pyr), 1 );
    % Label duplicate pairs (e.g., if neuron 5 and 10 are tagged,  
    % pairs (5,10) and (10,5) will have the same correlation)
    if length(pyr_opto_id)>1
        flag = true;
        toremove = nchoosek( pyr_opto_id, 2 );
        
        converge_tagtag = [ converge_tagtag ; convergence_mat( sub2ind( size(convergence_mat), toremove(:,1), toremove(:,2) ) ) ];
        
        ripcorr_tagtag = [ ripcorr_tagtag ; corrmat_ripple( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) ];
        thetacorr_tagtag = [ thetacorr_tagtag ; corrmat_theta( sub2ind( size(corrmat_theta), toremove(:,1), toremove(:,2) ) ) ];
        theta_tagtag_rates = [theta_tagtag_rates ; [rateHor_theta(sub2ind( size(rateHor_theta), toremove(:,1), toremove(:,2) ))...
                                              rateVert_theta(sub2ind( size(rateVert_theta), toremove(:,1), toremove(:,2) ))  ] ];
        ripple_tagtag_rates = [ripple_tagtag_rates ; [rateHor_rip(sub2ind( size(rateHor_rip), toremove(:,1), toremove(:,2) ))...
                                      rateVert_rip(sub2ind( size(rateVert_rip), toremove(:,1), toremove(:,2) ))  ] ];
        id_tag = [ id_tag ; repmat(tmp(1), size( toremove, 1), 1)];
        tagtag_dist = [ tagtag_dist ; dist_mat( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) ];
        tagtag_shankdist = [ tagtag_shankdist ; shankdist_mat( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) ];
        tagtag_lRat = [tagtag_lRat ; [lHor(sub2ind( size(lHor), toremove(:,1), toremove(:,2) ))...
                                              lVert(sub2ind( size(lVert), toremove(:,1), toremove(:,2) ))  ]];
        tagtag_wfoverlap = [tagtag_wfoverlap ; wf_overlap( sub2ind( size( wf_overlap ), toremove(:,1), toremove(:,2) ) ) ];
        animal_tag = [animal_tag ; repmat({mouse}, size( toremove, 1), 1)];
        id_session = [id_session ; tmp(1)];
        ripcorr_tagtag_session = [ripcorr_tagtag_session ; nanmedian(corrmat_ripple( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) )];
        propotion_tagged = [propotion_tagged ; length(pyr_opto_id) ./ length(pyr_id)];
        
        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
        convergence_mat( sub2ind( size(convergence_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        convergence_mat( sub2ind( size(convergence_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;
        
        
        rateHor_rip( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        rateHor_rip( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        rateVert_rip( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        rateVert_rip( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        
        rateHor_theta( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        rateHor_theta( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        rateVert_theta( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        rateVert_theta( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        
        lHor( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        lHor( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        lVert( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        lVert( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        
        wf_overlap( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        wf_overlap( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        
        dist_mat( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        dist_mat( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        shankdist_mat( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        shankdist_mat( sub2ind( size(corrmat_ripple), toremove(:,2), toremove(:,1) ) ) = randLabel;
        corrmat_ripple( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        corrmat_ripple( sub2ind( size(corrmat_ripple), toremove(:,2),  toremove(:,1) ) ) = randLabel;
        corrmat_theta( sub2ind( size(corrmat_ripple), toremove(:,1), toremove(:,2) ) ) = randLabel;
        corrmat_theta( sub2ind( size(corrmat_ripple), toremove(:,2),  toremove(:,1) ) ) = randLabel;
    end
    dist_mat(1:length(pyr_id)+1:end) = randLabel;
    shankdist_mat(1:length(pyr_id)+1:end) = randLabel;
    corrmat_ripple(1:length(pyr_id)+1:end) = randLabel;
    corrmat_theta(1:length(pyr_id)+1:end) = randLabel;
    wf_overlap(1:length(pyr_id)+1:end) = randLabel;
    rateHor_rip(1:length(pyr_id)+1:end) = randLabel;
    rateVert_rip(1:length(pyr_id)+1:end) = randLabel;
    rateHor_theta(1:length(pyr_id)+1:end) = randLabel;
    rateVert_theta(1:length(pyr_id)+1:end) = randLabel;
    lHor(1:length(pyr_id)+1:end) = randLabel;
    lVert(1:length(pyr_id)+1:end) = randLabel;
    convergence_mat(1:length(pyr_id)+1:end) = randLabel;
    
    convergence_vec = convergence_mat(pyr_opto_id,:);   convergence_vec = convergence_vec(:); 
    wf_overlap_vec = wf_overlap(pyr_opto_id,:);         wf_overlap_vec = wf_overlap_vec(:);
    rateHor_rip_vec = rateHor_rip(pyr_opto_id,:);       rateHor_rip_vec = rateHor_rip_vec(:);
    rateVert_rip_vec = rateVert_rip(pyr_opto_id,:);     rateVert_rip_vec = rateVert_rip_vec(:);
    rateHor_theta_vec = rateHor_theta(pyr_opto_id,:);   rateHor_theta_vec = rateHor_theta_vec(:);
    rateVert_theta_vec = rateVert_theta(pyr_opto_id,:); rateVert_theta_vec = rateVert_theta_vec(:);
    lHor_vec = lHor(pyr_opto_id,:);                     lHor_vec = lHor_vec(:);
    lVert_vec = lVert(pyr_opto_id,:);                   lVert_vec = lVert_vec(:);
    corrvec_ripple = corrmat_ripple(pyr_opto_id,:);     corrvec_ripple = corrvec_ripple(:);
    corrvec_theta = corrmat_theta(pyr_opto_id,:);       corrvec_theta = corrvec_theta(:); 
    distvec = dist_mat(pyr_opto_id,:);                  distvec = distvec(:);
    shankdistvec = shankdist_mat(pyr_opto_id,:);        shankdistvec = shankdistvec(:);
    
    % Create an index intp upper triangular matrix
    idx = triu( true( size(convergence_mat) ), 1);
    % Label all pairs involving tagged cells - both sbd and dbd
    convergence_mat(pyr_opto_id,:) = randLabel; convergence_mat(:,pyr_opto_id) = randLabel; 
    wf_overlap(pyr_opto_id,:) = randLabel;      wf_overlap(:,pyr_opto_id) = randLabel;
    rateHor_rip(pyr_opto_id,:) = randLabel;     rateHor_rip(:,pyr_opto_id) = randLabel;
    rateVert_rip(pyr_opto_id,:) = randLabel;    rateVert_rip(:,pyr_opto_id) = randLabel;
    rateHor_theta(pyr_opto_id,:) = randLabel;   rateHor_theta(:,pyr_opto_id) = randLabel;
    rateVert_theta(pyr_opto_id,:) = randLabel;  rateVert_theta(:,pyr_opto_id) = randLabel;
    lHor(pyr_opto_id,:) = randLabel;            lHor(:,pyr_opto_id) = randLabel;
    lVert(pyr_opto_id,:) = randLabel;           lVert(:,pyr_opto_id) = randLabel;
    corrmat_ripple(pyr_opto_id,:) = randLabel;  corrmat_ripple(:,pyr_opto_id) = randLabel;
    corrmat_theta(pyr_opto_id,:) = randLabel;   corrmat_theta(:,pyr_opto_id) = randLabel;
    dist_mat(pyr_opto_id,:) = randLabel;        dist_mat(:,pyr_opto_id) = randLabel;
    shankdist_mat(pyr_opto_id,:) = randLabel;   shankdist_mat(:,pyr_opto_id) = randLabel;
    % Extract the vectorized pairs, these have no labels
    convergence_vec1 = convergence_mat(idx);
    wf_overlap_vec1 = wf_overlap(idx);
    rateHor_rip_vec1 = rateHor_rip(idx);
    rateVert_rip_vec1 = rateVert_rip(idx);
    rateHor_theta_vec1 = rateHor_theta(idx);
    rateVert_theta_vec1 = rateVert_theta(idx);
    lHor_vec1 = lHor(idx);
    lVert_vec1 = lVert(idx);
    corrvec1_ripple = corrmat_ripple(idx);
    corrvec1_theta = corrmat_theta(idx);
    distvec1 = dist_mat(idx);
    shankdistvec1 = shankdist_mat(idx);
    
    % Remove the labeled pairs - dbd
    convergence_vec( convergence_vec == randLabel ) = [];
    wf_overlap_vec(wf_overlap_vec == randLabel) = [];
    corrvec_ripple( corrvec_ripple == randLabel ) = [];
    corrvec_theta( corrvec_theta == randLabel ) = [];
    distvec( distvec == randLabel ) = [];
    shankdistvec( shankdistvec == randLabel ) = [];
    rateHor_rip_vec(rateHor_rip_vec == randLabel) = [];
    rateVert_rip_vec(rateVert_rip_vec == randLabel) = [];
    rateHor_theta_vec(rateHor_theta_vec == randLabel) = [];
    rateVert_theta_vec(rateVert_theta_vec == randLabel) = [];
    lHor_vec(lHor_vec == randLabel) = [];
    lVert_vec(lVert_vec == randLabel) = [];
    % Remove the labeled pairs - untagged
    convergence_vec1( convergence_vec1 == randLabel ) = [];
    wf_overlap_vec1(wf_overlap_vec1 == randLabel) = [];
    corrvec1_ripple( corrvec1_ripple == randLabel ) = [];
    corrvec1_theta( corrvec1_theta == randLabel ) = [];
    distvec1( distvec1 == randLabel ) = [];
    shankdistvec1( shankdistvec1 == randLabel ) = [];
    rateHor_rip_vec1(rateHor_rip_vec1 == randLabel) = [];
    rateVert_rip_vec1(rateVert_rip_vec1 == randLabel) = [];
    rateHor_theta_vec1(rateHor_theta_vec1 == randLabel) = [];
    rateVert_theta_vec1(rateVert_theta_vec1 == randLabel) = [];
    lHor_vec1(lHor_vec1 == randLabel) = [];
    lVert_vec1(lVert_vec1 == randLabel) = [];
    
    % Store dbd as ntag
    converge_ntag = [ converge_ntag ; convergence_vec ];
    ntag_wfoverlap = [ntag_wfoverlap ; wf_overlap_vec];
    ripple_ntag_rates = [ripple_ntag_rates ; [rateHor_rip_vec rateVert_rip_vec]];
    theta_ntag_rates = [theta_ntag_rates ; [rateHor_theta_vec rateVert_theta_vec]];
    ntag_lRat = [ntag_lRat ; [lHor_vec lVert_vec]];
    ntag_dist = [ ntag_dist ; distvec ];
    ntag_shankdist = [ ntag_shankdist ; shankdistvec ];
    ripcorr_ntag = [ ripcorr_ntag ; corrvec_ripple ];
    thetacorr_ntag = [ thetacorr_ntag ; corrvec_theta ];
    id_ntag = [ id_ntag ; repmat(tmp(1), length(corrvec_ripple), 1)];
    animal_ntag = [animal_ntag ; repmat({mouse}, length(corrvec_ripple), 1)];
    
    % Store unlabelled as nolab
    converge_nolab = [ converge_nolab ; convergence_vec1 ];
    nolab_wfoverlap = [ nolab_wfoverlap ; wf_overlap_vec1 ];
    ripple_nolab_rates = [ ripple_nolab_rates ; [rateHor_rip_vec1 rateVert_rip_vec1] ];
    theta_nolab_rates = [ theta_nolab_rates ; [rateHor_theta_vec1 rateVert_theta_vec1] ];
    nolab_lRat = [ nolab_lRat ; [lHor_vec1 lVert_vec1] ];
    nolab_dist = [ nolab_dist ; distvec1 ];
    nolab_shankdist = [ nolab_shankdist ; shankdistvec1 ];
    ripcorr_nolab = [ ripcorr_nolab ; corrvec1_ripple ];
    thetacorr_nolab = [ thetacorr_nolab ; corrvec1_theta ];
    id_nolab = [ id_nolab ; repmat(tmp(1), length(corrvec1_ripple), 1) ];
    animal_nolab = [animal_nolab ; repmat({mouse}, length(corrvec1_ripple), 1) ];
    
    if flag
        ripcorr_nolab_session = [ripcorr_nolab_session ; nanmedian(corrvec1_ripple )];
    end

end

save('corrHome_v7.mat', 'theta_ntag_rates','ripple_ntag_rates', 'ntag_dist', 'ntag_shankdist','ripcorr_ntag', 'thetacorr_ntag', 'id_ntag', ...
     'ripcorr_tagtag', 'thetacorr_tagtag', 'theta_tagtag_rates', 'ripple_tagtag_rates', 'id_tag', 'tagtag_dist', 'tagtag_shankdist', ...
     'ntag_lRat', 'tagtag_lRat', 'animal_ntag', 'animal_tag', 'id_session', 'ripcorr_tagtag_session', 'propotion_tagged', 'ripcorr_ntag_session',...
     'ntag_wfoverlap', 'tagtag_wfoverlap', 'converge_nolab', 'nolab_wfoverlap', 'ripple_nolab_rates', 'theta_nolab_rates', 'nolab_lRat', ...
     'nolab_dist', 'nolab_shankdist', 'ripcorr_nolab', 'thetacorr_nolab', 'id_nolab', 'animal_nolab', 'ripcorr_nolab_session')

 

%%