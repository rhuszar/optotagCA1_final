


basepaths_famnov = [{'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e13\e13_26m1\e13_26m1_211019'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e13\e13_26m1\e13_26m1_211116'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e13\e13_26m1\e13_26m1_211119'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e13\e13_26m1\e13_26m1_211124'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211210'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211211'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211212'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e16\e16_3m2\e16_3m2_211213'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220117'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220118'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220119'} ; ...
                    {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\e15\e15_13f1\e15_13f1_220120'}];

                
%%
tagged = [];
pyr = [];
for kp = 1:length(basepaths_famnov)
    basepath = alterPath( basepaths_famnov{kp}, true);
    basename = bz_BasenameFromBasepath(basepath);
    load( fullfile(basepath, [ basename '.cell_metrics.cellinfo.mat']) )
    
    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id & pyr_id;
    tagged(kp) = sum(pyr_opto_id);
    pyr(kp) = sum(pyr_id);
end
%% Update the spatial information metrics for 2D exploration

for kp = 1:length(basepaths_famnov)
    basepath = alterPath( basepaths_famnov{kp}, true);
    basename = bz_BasenameFromBasepath(basepath);
    fprintf('%d/%d\n', kp, length(basepaths_famnov))
    ratePath = fullfile(basepath, [basename '.firingMapsAvg2D.cellinfo.mat']);
    if exist(ratePath)
        load(ratePath);
    else
        continue
    end
   
    occ = firingMaps.occupancy{1}(:)'; p = occ ./ sum(occ);
    Isec1 = [];
    Isec2 = [];
    Ispike = [];
    for jp = 1:length( firingMaps.UID )
        rmap = firingMaps.rateMaps{jp}; rmap = rmap(:)';
        Isec1(jp) = nanmean( rmap  .* log2( rmap  ./ nanmean( rmap ) ) );
        Ispike(jp) = Isec1(jp) ./ nanmean( rmap );
%         tmp = rmap  ./ nanmean( rmap ); tmp( tmp == 0 ) = .000000001;
%         Isec2(jp) = nansum( p' .* [ rmap  .* log2( tmp ) ]' );
    end
    firingMaps.Ispike = Ispike;
    firingMaps.Isec = Isec1;
    save(ratePath, 'firingMaps');
end

%% Correlations between spatial firing maps of optotagged neurons!
% theta cycle firing rate - 
% compute separately for novel and familiar load the waveforms

% same birthdate versus different birthdate
% Do we wish to load data from local hard drive? 
isLocal = true;

% Do we consider the stem of the maze, or not
consider_stem = false;
% Indices of what's after the maze stem
post_stem_inds_fam = 75:205;
post_stem_inds_nov = 51:164;

wfoverlap_tag = [];
pfcorr_tag_all = [];
theta_tag_rates_byMaze = cell(1,3);
theta_tag_rates_all = [];
pfcorr_tag_dist = [];
pfcorr_tag_shankdist = [];
lRat_tag = [];
id_tag = [];
animal_tag = [];
mazeid_pairs_tag = [];
basename_pairs_tag = [];

wfoverlap_ntag = [];
pfcorr_ntag_all = [];
theta_ntag_rates_byMaze = cell(1,3);
theta_ntag_rates_all = [];
pfcorr_ntag_dist = [];
pfcorr_ntag_shankdist = [];
lRat_ntag = [];
id_ntag = [];
animal_ntag = [];
mazeid_pairs_ntag = [];
basename_pairs_ntag = [];

% Theta cycle firing rates of individual neurons
rateTheta_all_tag = [];
rateTheta_byMaze_tag = [];
rateTheta_nov_tag = [];
rateTheta_fam_tag= [];
basename_cell_tag = [];
mazeid_cell_tag = [];
rateTheta_all_ntag = [];
rateTheta_byMaze_ntag  = [];
rateTheta_nov_ntag = [];
rateTheta_fam_ntag = [];
basename_cell_ntag = [];
mazeid_cell_ntag = [];
isec_byMaze_tag = [];
ispike_byMaze_tag = [];
isec_byMaze_ntag = [];
ispike_byMaze_ntag = [];

randLabel = 2000*rand;

%%
for bp = 1:length(basepaths_famnov)
   

    fprintf('%d/%d\n', bp, length(basepaths_famnov))
    basepath = alterPath( basepaths_famnov{bp}, isLocal );
    %basepath = basepaths_famnov{bp};
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    mouse = [ tmp{1} '_' tmp{2} ];
    clear v1 v2
    v1 = load(fullfile(basepath, [basename '.firingMapsAvg_multSess.cellinfo.mat']));
    v1_len = length( v1.firingMaps );
    try
        v2 = load(fullfile(basepath, [basename '.firingMapsAvg2D.cellinfo.mat']));
        v2_len = length( v2.firingMaps );
    catch
        v2_len = 0;
    end
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.thetaCycles.events.mat']))
    load(fullfile(basepath, [basename '.SleepState.states.mat']))
    load(fullfile(basepath, [basename '.Tracking.behavior.mat']))
    cluQual = load(fullfile(basepath, 'PYR_clusterQuality.mat'));
    load( fullfile(basepath, 'wf1.mat') )
    load(fullfile(basepath, 'chanMap.mat'))
    load(fullfile(basepath, 'maze_info.mat'))
    v.chanMap = [xcoords ; ycoords];
    v.kcoords = kcoords;
    nMazes = size( tracking.events.subSessions,1 );
    if nMazes < 3
        maze_info(3) = {'n/a'};
    end
    
    opto_id_log = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id_log = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id_log & pyr_id_log;
    pyr_id = find(pyr_id_log);
    
    % Get the L ratio
    [~, loc] = ismember(pyr_id, cluQual.pyr_UID);
    L_Rat = cluQual.L_Rat(loc); L_Rat = L_Rat(:);
    
    % Prepare theta intervals
    thetaCycle_mazeIndicator = arrayfun(@(x) InIntervals( Theta.cycles(:,1), tracking.events.subSessions(x,:) ) & ...
         InIntervals( Theta.cycles(:,2), tracking.events.subSessions(x,:) ),  1:nMazes, 'UniformOutput', false );
    thetaCycle_byMaze = cellfun(@(x) Theta.cycles(x,:), thetaCycle_mazeIndicator, 'UniformOutput', false);
    thetaCycle_fam = thetaCycle_byMaze{1};
    thetaCycle_nov = cell2mat( thetaCycle_byMaze(2:end)' );
    thetaCycle_all = cell2mat( thetaCycle_byMaze(1:end)' );
    
    theta_byMaze_dur = cellfun(@(x) sum( diff( x' ) ), thetaCycle_byMaze );
    theta_fam_dur = theta_byMaze_dur(1);
    theta_nov_dur = sum( theta_byMaze_dur(2:end) );
    theta_all_dur = sum( theta_byMaze_dur(1:end) );
    
    % There are indices into pyr_id - i.e., if 5th index in pyr_id is optotagged, you want to store it 
    pyr_opto_id = arrayfun(@(x) find(x==pyr_id),  find(pyr_opto_id) );
    pyr_nopto_id = setdiff(1:length(pyr_id), pyr_opto_id);
    
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

    % Collect ratemaps from linearizable mazes (tmaze & linear track)
    ratemaps_all = {};
    ispike_all = {};
    isec_all = {};
    for jp = 1:v1_len
        % Skip nonlinearizable session
        if isempty( v1.firingMaps(jp).rateMaps )
            continue
        end
        ispike_all{jp} = v1.firingMaps(jp).Ispike(:);
        isec_all{jp} = v1.firingMaps(jp).Isec(:);
        for ip = 1:length(pyr_id) 
            % Familiar linearizable track
            if jp == 1
                if consider_stem
                    ratemaps_all{jp}(:,ip) = cell2mat( v1.firingMaps(jp).rateMaps{pyr_id(ip)} );
                else
                    ratemaps_all{jp}(:,ip) = [v1.firingMaps(jp).rateMaps{pyr_id(ip)}{1}( post_stem_inds_fam ) ...
                                                        v1.firingMaps(jp).rateMaps{pyr_id(ip)}{2}( post_stem_inds_fam )];
                end
            end
            % Novel linearizable track - either tmaze, or linear track
            if jp>1 && strcmp(maze_info{jp}, 'tmaze_ipshita')
                if consider_stem
                    ratemaps_all{jp}(:,ip) = cell2mat( v1.firingMaps(jp).rateMaps{pyr_id(ip)} );
                else
                    ratemaps_all{jp}(:,ip) = [v1.firingMaps(jp).rateMaps{pyr_id(ip)}{1}( post_stem_inds_nov ) ...
                                                        v1.firingMaps(jp).rateMaps{pyr_id(ip)}{2}( post_stem_inds_nov )];
                end
            elseif jp>1
                ratemaps_all{jp}(:,ip) = cell2mat( v1.firingMaps(jp).rateMaps{pyr_id(ip)} );
            end
        end
    end
    
    index = setdiff( 1:nMazes, find( ~cellfun(@isempty, ratemaps_all ) ) );
    % Collect ratemaps from open fields, if we have any
    if v2_len ~= 0
        ispike_all{index} = v2.firingMaps.Ispike(:);
        isec_all{index} = v2.firingMaps.Isec(:);
        for ip = 1:length(pyr_id) 
            ratemaps_all{ index }(:,ip) = v2.firingMaps.rateMaps_nan{pyr_id(ip)}(:);
        end
    end

    thetaCycle_famnov_temp = nan(length(pyr_id), 2);
    thetaCycle_byMaze_temp = nan(length(pyr_id), 3);
    thetaCycle_all_temp = nan(length(pyr_id), 1);
    % Number of spikes per theta cycle
    for jp = 1:length(pyr_id)
        % Spike counts in theta cycles
        thetaSpk_fam = InIntervals( spikes.times{ pyr_id(jp) }, thetaCycle_fam );
        thetaSpk_nov = InIntervals( spikes.times{ pyr_id(jp) }, thetaCycle_nov );
        thetaSpk_all = InIntervals( spikes.times{ pyr_id(jp) }, thetaCycle_all );
        thetaSpk_byMaze = cellfun(@(x) InIntervals( spikes.times{ pyr_id(jp) }, x), thetaCycle_byMaze, 'UniformOutput', false);
        % Theta cycle rates
        thetaCycle_famnov_temp(jp,:) = [sum(thetaSpk_fam) ./ theta_fam_dur  sum(thetaSpk_nov) ./ theta_nov_dur ];
        thetaCycle_byMaze_temp(jp,1:nMazes) = cellfun(@sum, thetaSpk_byMaze ) ./ theta_byMaze_dur;
        thetaCycle_all_temp(jp) = sum(thetaSpk_all) ./ theta_all_dur;
    end

    basename_cell_tag = [basename_cell_tag ; repmat({basename}, length( pyr_opto_id ), 1)];
    basename_cell_ntag = [basename_cell_ntag ; repmat({basename}, length( pyr_nopto_id ), 1)];
    mazeid_cell_tag = [mazeid_cell_tag ; repmat( maze_info, length( pyr_opto_id ), 1)];
    mazeid_cell_ntag = [mazeid_cell_ntag ; repmat( maze_info, length( pyr_nopto_id ), 1)];
    rateTheta_nov_tag = [rateTheta_nov_tag ; thetaCycle_famnov_temp(pyr_opto_id,2)];
    rateTheta_nov_ntag = [rateTheta_nov_ntag ; thetaCycle_famnov_temp(pyr_nopto_id,2)];
    rateTheta_fam_tag= [rateTheta_fam_tag ; thetaCycle_famnov_temp(pyr_opto_id,1)];
    rateTheta_fam_ntag = [rateTheta_fam_ntag ; thetaCycle_famnov_temp(pyr_nopto_id,1)];
    rateTheta_byMaze_tag = [rateTheta_byMaze_tag ; thetaCycle_byMaze_temp(pyr_opto_id,:)];
    rateTheta_byMaze_ntag = [rateTheta_byMaze_ntag ; thetaCycle_byMaze_temp(pyr_nopto_id,:)];
    rateTheta_all_tag = [ rateTheta_all_tag ; thetaCycle_all_temp(pyr_opto_id)];
    rateTheta_all_ntag = [ rateTheta_all_ntag ; thetaCycle_all_temp(pyr_nopto_id)];
    
    
    % correlation + theta rate matrices
    corrmat_all = {}; rateHor_theta_byMaze = {}; rateVert_theta_byMaze = {};
    for ip = 1:nMazes
        corrmat_all{ip} = corr( ratemaps_all{ip}, 'rows', 'complete' );
        rateHor_theta_byMaze{ip} = repmat( thetaCycle_byMaze_temp(:,ip), 1, length( thetaCycle_byMaze_temp(:,ip) ) );
        rateVert_theta_byMaze{ip} = repmat( thetaCycle_byMaze_temp(:,ip)', length( thetaCycle_byMaze_temp(:,ip)), 1 );
    end
    % put nans in place of a matrices for the 3rd maze, in case it's missing
    if nMazes < 3
        corrmat_all{3} = nan( size(corrmat_all{ip}) );
        rateHor_theta_byMaze{3} = nan( size(corrmat_all{ip}) );
        rateVert_theta_byMaze{3} = nan( size(corrmat_all{ip}) );
        ispike_all{3} = nan(size(ispike_all{2}));
        isec_all{3} = nan(size(isec_all{2}));
    end
    
    ispike_all = cell2mat( ispike_all );
    isec_all = cell2mat( isec_all );
    
    isec_byMaze_tag = [isec_byMaze_tag ; isec_all(pyr_opto_id,:) ];
    ispike_byMaze_tag = [ispike_byMaze_tag ; ispike_all(pyr_opto_id,:) ];
    isec_byMaze_ntag = [isec_byMaze_ntag ; isec_all(pyr_nopto_id,:) ];
    ispike_byMaze_ntag = [ispike_byMaze_ntag ; ispike_all(pyr_nopto_id,:)  ];
    
    rateHor_theta_all = repmat( thetaCycle_all_temp, 1, length( thetaCycle_all_temp ) );
    rateVert_theta_all = repmat( thetaCycle_all_temp', length( thetaCycle_all_temp ), 1 );
    lHor = repmat( L_Rat, 1, length(L_Rat) );
    lVert = repmat( L_Rat', length(L_Rat), 1 );
    
    % Label duplicate pairs (e.g., if neuron 5 and 10 are tagged,  
    % pairs (5,10) and (10,5) will have the same correlation)
    if length(pyr_opto_id)>1
        toremove = nchoosek( pyr_opto_id, 2 );
        
        % Store ratemap correlations associated with each maze - last column might be nans
        pfcorr_tag_all = [ pfcorr_tag_all ; [ corrmat_all{1}( sub2ind( size(corrmat_all{ip}), toremove(:,1), toremove(:,2) ) ) ...
                                              corrmat_all{2}( sub2ind( size(corrmat_all{ip}), toremove(:,1), toremove(:,2) ) ) ...
                                        corrmat_all{3}( sub2ind( size(corrmat_all{ip}), toremove(:,1), toremove(:,2) ) )] ];
        
        % Store the pairwise distances - accounting for the number of mazes 
        pfcorr_tag_dist = [ pfcorr_tag_dist ; dist_mat( sub2ind( size(dist_mat), toremove(:,1), toremove(:,2) ) ) ];
        pfcorr_tag_shankdist = [ pfcorr_tag_shankdist ; shankdist_mat( sub2ind( size(shankdist_mat), toremove(:,1), toremove(:,2) ) ) ];
        wfoverlap_tag = [wfoverlap_tag ; wf_overlap( sub2ind( size( wf_overlap ), toremove(:,1), toremove(:,2) ) ) ];
        lRat_tag = [lRat_tag ; [lHor(sub2ind( size(lHor), toremove(:,1), toremove(:,2) ))...
                                      lVert(sub2ind( size(lVert), toremove(:,1), toremove(:,2) ))  ]];
        theta_tag_rates_all = [theta_tag_rates_all ; [rateHor_theta_all(sub2ind( size(rateHor_theta_all), toremove(:,1), toremove(:,2) ))...
                                      rateVert_theta_all(sub2ind( size(rateVert_theta_all), toremove(:,1), toremove(:,2) ))  ]];
        % Identifiers
        animal_tag = [animal_tag ; repmat({mouse}, size( toremove, 1), 1)];
        id_tag = [ id_tag ; repmat(tmp(1), size( toremove, 1), 1) ];
        mazeid_pairs_tag = [ mazeid_pairs_tag ; repmat(maze_info, size( toremove, 1), 1) ];
        basename_pairs_tag = [ basename_pairs_tag ; repmat( {basename}, size( toremove, 1), 1) ];
        
        % Store joint firing rates (saved as cell array) & label to-be-removed tagged pairs
        for ip = 1:3
            theta_tag_rates_byMaze{ip} = [theta_tag_rates_byMaze{ip} ; [rateHor_theta_byMaze{ip}(sub2ind( size(rateHor_theta_byMaze{ip}), toremove(:,1), toremove(:,2) ))...
                                              rateVert_theta_byMaze{ip}(sub2ind( size(rateVert_theta_byMaze{ip}), toremove(:,1), toremove(:,2) ))  ] ];
            
            corrmat_all{ip}( sub2ind( size(corrmat_all{ip}), toremove(:,1), toremove(:,2) ) ) = randLabel;
            corrmat_all{ip}( sub2ind( size(corrmat_all{ip}), toremove(:,2),  toremove(:,1) ) ) = randLabel;
            rateHor_theta_byMaze{ip}( sub2ind( size(corrmat_all{ip}), toremove(:,1), toremove(:,2) ) ) = randLabel;
            rateHor_theta_byMaze{ip}( sub2ind( size(corrmat_all{ip}), toremove(:,2), toremove(:,1) ) ) = randLabel;
            rateVert_theta_byMaze{ip}( sub2ind( size(corrmat_all{ip}), toremove(:,1), toremove(:,2) ) ) = randLabel;
            rateVert_theta_byMaze{ip}( sub2ind( size(corrmat_all{ip}), toremove(:,2), toremove(:,1) ) ) = randLabel;
        end
        
        dist_mat( sub2ind( size(dist_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        dist_mat( sub2ind( size(dist_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;
        shankdist_mat( sub2ind( size(shankdist_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        shankdist_mat( sub2ind( size(shankdist_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;
        wf_overlap( sub2ind( size(wf_overlap), toremove(:,1), toremove(:,2) ) ) = randLabel;
        wf_overlap( sub2ind( size(wf_overlap), toremove(:,2), toremove(:,1) ) ) = randLabel;
        lHor( sub2ind( size(lHor), toremove(:,1), toremove(:,2) ) ) = randLabel;
        lHor( sub2ind( size(lHor), toremove(:,2), toremove(:,1) ) ) = randLabel;
        lVert( sub2ind( size(lVert), toremove(:,1), toremove(:,2) ) ) = randLabel;
        lVert( sub2ind( size(lVert), toremove(:,2), toremove(:,1) ) ) = randLabel;
        rateHor_theta_all( sub2ind( size(rateHor_theta_all), toremove(:,1), toremove(:,2) ) ) = randLabel;
        rateHor_theta_all( sub2ind( size(rateHor_theta_all), toremove(:,2), toremove(:,1) ) ) = randLabel;
        rateVert_theta_all( sub2ind( size(rateVert_theta_all), toremove(:,1), toremove(:,2) ) ) = randLabel;
        rateVert_theta_all( sub2ind( size(rateVert_theta_all), toremove(:,2), toremove(:,1) ) ) = randLabel;
   
    end
    
    % Label the diagonal, pull out pairs w/ a tagged cell, exclude tag-tag pairs 
    dist_mat(1:length(pyr_id)+1:end) = randLabel;
    shankdist_mat(1:length(pyr_id)+1:end) = randLabel;
    wf_overlap(1:length(pyr_id)+1:end) = randLabel;
    lHor(1:length(pyr_id)+1:end) = randLabel;
    lVert(1:length(pyr_id)+1:end) = randLabel;
    rateHor_theta_all(1:length(pyr_id)+1:end) = randLabel;
    rateVert_theta_all(1:length(pyr_id)+1:end) = randLabel;
 
    distvec = dist_mat(pyr_opto_id,:);                        distvec = distvec(:);
    shankdistvec = shankdist_mat(pyr_opto_id,:);              shankdistvec = shankdistvec(:);
    wf_overlap_vec = wf_overlap(pyr_opto_id,:);               wf_overlap_vec = wf_overlap_vec(:);
    lHor = lHor(pyr_opto_id,:);                               lHor = lHor(:);
    lVert = lVert(pyr_opto_id,:);                             lVert = lVert(:);
    rateHor_theta_all = rateHor_theta_all(pyr_opto_id,:);     rateHor_theta_all = rateHor_theta_all(:);
    rateVert_theta_all = rateVert_theta_all(pyr_opto_id,:);   rateVert_theta_all = rateVert_theta_all(:);
    
    distvec( distvec == randLabel ) = [];
    shankdistvec( shankdistvec == randLabel ) = [];
    wf_overlap_vec(wf_overlap_vec == randLabel) = [];
    lHor(lHor == randLabel) = [];
    lVert(lVert == randLabel) = [];
    rateHor_theta_all(rateHor_theta_all == randLabel) = [];
    rateVert_theta_all(rateVert_theta_all == randLabel) = [];
    
    % Same as above, for correlations / joint rates on each maze
    corrvec = cell(1,3);
    rateHor_vec = cell(1,3);
    rateVert_vec = cell(1,3);
    for ip = 1:3
        % Label diagonal
        corrmat_all{ip}(1:length(pyr_id)+1:end) = randLabel;
        rateHor_theta_byMaze{ip}(1:length(pyr_id)+1:end) = randLabel;
        rateVert_theta_byMaze{ip}(1:length(pyr_id)+1:end) = randLabel;
        % Pull out pairs involving tagged cells
        corrvec{ip} = corrmat_all{ip}(pyr_opto_id,:);   corrvec{ip} = corrvec{ip}(:);
        rateHor_vec{ip} = rateHor_theta_byMaze{ip}(pyr_opto_id,:);   rateHor_vec{ip} = rateHor_vec{ip}(:);
        rateVert_vec{ip} = rateVert_theta_byMaze{ip}(pyr_opto_id,:);   rateVert_vec{ip} = rateVert_vec{ip}(:);
        % Remove tag-tag pairs
        corrvec{ip}( corrvec{ip} == randLabel ) = [];
        rateHor_vec{ip}(rateHor_vec{ip} == randLabel) = [];
        rateVert_vec{ip}(rateVert_vec{ip} == randLabel) = [];
    end
   
    pfcorr_ntag_all = [ pfcorr_ntag_all ; cell2mat( corrvec ) ];
    pfcorr_ntag_dist = [ pfcorr_ntag_dist ; distvec ];
    pfcorr_ntag_shankdist = [ pfcorr_ntag_shankdist ; shankdistvec ];
    wfoverlap_ntag = [wfoverlap_ntag ; wf_overlap_vec ];
    id_ntag = [ id_ntag ; repmat(tmp(1), length(distvec), 1) ];
    animal_ntag = [ animal_ntag ; repmat({mouse}, length(distvec), 1) ];
    mazeid_pairs_ntag = [ mazeid_pairs_ntag ; repmat(maze_info, length(distvec), 1) ];
    basename_pairs_ntag = [ basename_pairs_ntag ; repmat( {basename}, length(distvec), 1) ];
    lRat_ntag = [lRat_ntag ; [lHor lVert]];
    theta_ntag_rates_all = [theta_ntag_rates_all ; [rateHor_theta_all rateVert_theta_all]];

    for ip = 1:3
        theta_ntag_rates_byMaze{ip} = [theta_ntag_rates_byMaze{ip} ; [rateHor_vec{ip} rateVert_vec{ip}]];
    end
   
end

%%

save('corrBeh_nov_v1.mat', 'pfcorr_ntag_all', 'pfcorr_ntag_dist', 'pfcorr_ntag_shankdist', 'wfoverlap_ntag', 'id_ntag',...
    'animal_ntag', 'mazeid_pairs_ntag', 'basename_pairs_ntag', 'lRat_ntag', 'theta_ntag_rates', 'pfcorr_tag_all',...
    'pfcorr_tag_dist', 'pfcorr_tag_shankdist', 'wfoverlap_tag', 'lRat_tag', 'animal_tag', 'id_tag', 'mazeid_pairs_tag',...
    'basename_pairs_tag', 'theta_tag_rates', 'basename_cell_tag', 'basename_cell_ntag', 'mazeid_cell_tag', 'mazeid_cell_ntag', ...
    'rateTheta_nov_tag', 'rateTheta_nov_ntag', 'rateTheta_fam_tag', 'rateTheta_fam_ntag', 'rateTheta_byMaze_tag', 'rateTheta_byMaze_ntag', ...
    'theta_ntag_rates_all', 'theta_tag_rates_all', 'rateTheta_all_tag', 'rateTheta_all_ntag')

%% pfcorr by shank distance on a session by session basis

uniq_bp = unique( basename_pairs_tag );
pfcorr_by_shankdist_tag = [];
pfcorr_by_shankdist_ntag = [];
for kp = 1:length(uniq_bp)
    indic_tag = cellfun(@(x) strcmp(x, uniq_bp{kp}), basename_pairs_tag);
    indic_ntag = cellfun(@(x) strcmp(x, uniq_bp{kp}), basename_pairs_ntag);
    cc_tag = [ pfcorr_tag_all(indic_tag,2) ; pfcorr_tag_all(indic_tag,3) ];
    cc_ntag = [ pfcorr_ntag_all(indic_ntag,2) ; pfcorr_ntag_all(indic_ntag,3) ];
    shankd_tag = [ pfcorr_tag_shankdist(indic_tag) ; pfcorr_tag_shankdist(indic_tag) ];
    shankd_ntag = [ pfcorr_ntag_shankdist(indic_ntag) ; pfcorr_ntag_shankdist(indic_ntag) ];
    pfcorr_by_shankdist_tag = [pfcorr_by_shankdist_tag ; [ nanmean( cc_tag(shankd_tag==0) ) nanmean( cc_tag(shankd_tag==1) ) nanmean( cc_tag(shankd_tag==2) ) ] ];
    pfcorr_by_shankdist_ntag = [pfcorr_by_shankdist_ntag ; [ nanmean( cc_ntag(shankd_ntag==0) ) nanmean( cc_ntag(shankd_ntag==1) ) nanmean( cc_ntag(shankd_ntag==2) ) ] ];
end

%%
figure
set(gcf,'Position', [455   289   282   326])
plot( pfcorr_by_shankdist_tag', 'b' )
hold on
plot( nanmean( pfcorr_by_shankdist_tag ), 'k', 'LineWidth', 2 )
plot( pfcorr_by_shankdist_ntag', 'r' )
plot( nanmean( pfcorr_by_shankdist_ntag ), 'm', 'LineWidth', 2 )


%% Ratemap remap of SBD v. DBD between familiar and novel environments

minvals = 100;
step = 33;
vals = [ pfcorr_tag_all(:) ; pfcorr_ntag_all(:) ];
vals(isnan(vals)) = [];
edges = prctile( vals,0:step:100);
edges(1) = -inf; edges(end) = inf;
% edges = prctile([ pf_ntag_left ; pf_ntag_right ] ,0:step:100);
% Familiar
[cnts_tt_fam, ~, ib_tt_fam] = histcounts( pfcorr_tag_all(:,1), edges );
[cnts_nt_fam, ~, ib_nt_fam] = histcounts( pfcorr_ntag_all(:,1), edges );
% Novel - session 1
[cnts_tt_nov1, ~, ib_tt_nov1] = histcounts( pfcorr_tag_all(:,2), edges );
[cnts_nt_nov1, ~, ib_nt_nov1] = histcounts( pfcorr_ntag_all(:,2), edges );
% Novel - session 2
[cnts_tt_nov2, ~, ib_tt_nov2] = histcounts( pfcorr_tag_all(:,3), edges );
[cnts_nt_nov2, ~, ib_nt_nov2] = histcounts( pfcorr_ntag_all(:,3), edges );

good_bins = find( cnts_tt_fam > minvals );

tt_mu = []; tt_sem  = [];
nt_mu = []; nt_sem = [];
tt = [];
nt = [];
for kp = good_bins
   tt = [tt {[ pfcorr_tag_all( ib_tt_fam == kp,2 ) ; pfcorr_tag_all( ib_tt_fam == kp,3 )]}];
   nt = [nt {[ pfcorr_ntag_all( ib_nt_fam == kp, 2 ) ; pfcorr_ntag_all( ib_nt_fam == kp, 3 ) ]}]; 
   tt_mu = [ tt_mu ; nanmean( tt{kp} ) ]; tt_sem = [tt_sem ; nansem( tt{kp} )];
   nt_mu = [ nt_mu ; nanmean( nt{kp} )  ]; nt_sem = [nt_sem ; nansem( nt{kp} )];
end

xv = [1:length(good_bins)]';
ftt = [tt_mu+tt_sem ; flip( tt_mu-tt_sem )];
fnt = [nt_mu+nt_sem ; flip( nt_mu-nt_sem )];

figure
set(gcf,'Position',[440   377   306   420])
fill( [xv ; flip(xv)], ftt, 'b', 'facealpha', 0.1, 'linestyle', 'none' )
hold on
fill( [xv ; flip(xv)], fnt, 'r', 'facealpha', 0.1, 'linestyle', 'none' )
plot(tt_mu, '-b', 'linewidth', 1)
plot(nt_mu, '-r', 'linewidth', 1)
% errorbar(xv,tt_mu,tt_sem, '-s', 'MarkerSize',8,...
% 'MarkerEdgeColor','blue','MarkerFaceColor','blue','CapSize',18')
% hold on
% errorbar(xv,nt_mu,nt_sem, '-s', 'MarkerSize',8,...
% 'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',18')
% axis square
% xticks(1:length(good_bins))

xlabel('Ratemap corr. on Trial A (percentiles)', 'fontsize', 14)
ylabel('Ratemap corr. on Trial B', 'fontsize', 14)

%%

y = [ cell2mat( tt') ; cell2mat(nt') ];
gp1 = [repmat({'tag'}, size(cell2mat( tt'))) ; repmat({'ntag'}, size(cell2mat( nt'))) ];
gp2 = [ ones(size( tt{1} ) ) ; 2.*ones(size( tt{2} )) ; 3.*ones(size( tt{3} )) ; ...];
        ones(size( nt{1} ) ) ; 2.*ones(size( nt{2} )) ; 3.*ones(size( nt{3} ))  ];

[p,tbl,stats] = anovan(y,{gp1,gp2},'model','interaction',...
    'varnames',{'g1','g2'});
[c,~,~,gnames] = multcompare(stats,'Dimension',[1 2]);

[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))]


%% familiar_novel_ratemapCorr

ind1_t = ~isnan( pfcorr_tag_all(:,1) );      ind2_t = ~isnan( pfcorr_tag_all(:,2) );       ind3_t = ~isnan( pfcorr_tag_all(:,3) );
ind1_nt = ~isnan( pfcorr_ntag_all(:,1) );    ind2_nt = ~isnan( pfcorr_ntag_all(:,2) );     ind3_nt = ~isnan( pfcorr_ntag_all(:,3) );

tt = [ pfcorr_tag_all(ind1_t & ind2_t,1) pfcorr_tag_all(ind1_t & ind2_t,2)  ;...
       pfcorr_tag_all(ind1_t & ind3_t,1) pfcorr_tag_all(ind1_t & ind3_t,3) ];
nt = [ pfcorr_ntag_all(ind1_nt & ind2_nt,1) pfcorr_ntag_all(ind1_nt & ind2_nt,2)  ;...
       pfcorr_ntag_all(ind1_nt & ind3_nt,1) pfcorr_ntag_all(ind1_nt & ind3_nt,3) ];

   plot(nt(:,1), nt(:,2), '.r')
   
   hold on
plot(tt(:,1), tt(:,2), '.b')
lsline
xlabel('Spatial ratemap correlation in FAM')
ylabel('Spatial ratemap correlation in NOV')
legend('different birthdate', 'same birthdate','location', 'best')

%% 

ip1 = cellfun(@(x) strcmp(x, basename_pairs_tag(ind1_t & ind2_t)), unique( basename_pairs_tag ) , 'UniformOutput', false);
ip2 = cellfun(@(x) strcmp(x, basename_pairs_tag(ind1_t & ind3_t)), unique( basename_pairs_tag ) , 'UniformOutput', false);
v = [ ip1 ip2 ];
v_tt = arrayfun(@(x) cell2mat( v(x,:)' ), 1:length(unique( basename_pairs_tag )), 'UniformOutput', false);
ip1 = cellfun(@(x) strcmp(x, basename_pairs_ntag(ind1_nt & ind2_nt)), unique( basename_pairs_ntag ) , 'UniformOutput', false);
ip2 = cellfun(@(x) strcmp(x, basename_pairs_ntag(ind1_nt & ind3_nt)), unique( basename_pairs_ntag ) , 'UniformOutput', false);
v = [ ip1 ip2 ];
v_nt = arrayfun(@(x) cell2mat( v(x,:)' ), 1:length(unique( basename_pairs_ntag )), 'UniformOutput', false);

tt_mu = cell2mat( cellfun(@(x) [ nanmean(tt(x,1)) nanmean(tt(x,2)) ], v_tt, 'UniformOutput', false)' );
tt_sem = cell2mat( cellfun(@(x) [ nansem(tt(x,1)) nansem(tt(x,2)) ], v_tt, 'UniformOutput', false)' );
nt_mu = cell2mat( cellfun(@(x) [ nanmean(nt(x,1)) nanmean(nt(x,2)) ], v_nt, 'UniformOutput', false)' );
nt_sem = cell2mat( cellfun(@(x) [ nansem(nt(x,1)) nansem(nt(x,2)) ], v_nt, 'UniformOutput', false)' );

errorbar(tt_mu(:,1),tt_mu(:,2),tt_sem(:,2),tt_sem(:,2),tt_sem(:,1),tt_sem(:,1),'.b')
hold on
errorbar(nt_mu(:,1),nt_mu(:,2),nt_sem(:,2),nt_sem(:,2),nt_sem(:,1),nt_sem(:,1),'.r')
%%

vals = [ pfcorr_tag_all(:) ; pfcorr_ntag_all(:) ];
nanvals = isnan(vals);
pd = fitdist(vals(~nanvals),'kernel','Width',0.05);
x = min( vals(~nanvals) ):0.01:max( vals(~nanvals) );
y = pdf(pd,x);

figure
set(gcf,'Position', [1000         918         560         281])
histogram(vals(~nanvals),25, 'Normalization', 'pdf', 'DisplayStyle', 'stairs' );
hold on

plot(x,y  ,'Color','k','LineStyle','--')
xlim([min(x) max(x)])

ylim_vals = ylim;
xlim_vals = xlim;

boundaries = edges(2:end-1); % we use edges defined by the plot for visualization
b1 = rectangle('Position', [ xlim_vals(1), ylim_vals(1), boundaries(1) - xlim_vals(1), ylim_vals(2) ], 'FaceColor',[0.9290 0.6940 0.1250], 'linestyle', 'none');
b1.FaceColor(4) = 0.4;
b2 = rectangle('Position', [ boundaries(1), ylim_vals(1), boundaries(2) - boundaries(1), ylim_vals(2) ], 'FaceColor',[0.9290 0.6940 0.1250], 'linestyle', 'none');
b2.FaceColor(4) = 0.6;
b3 = rectangle('Position', [ boundaries(2), ylim_vals(1), xlim_vals(2) - boundaries(2), ylim_vals(2) ], 'FaceColor',[0.9290 0.6940 0.1250], 'linestyle', 'none');
b3.FaceColor(4) = 0.8;

% Now divide this distribution into percentiles
xlabel('Ratemap corr. in FAM', 'fontsize', 14)
ylabel('Probability density', 'fontsize', 14)

%% FAMILIAR VERSUS NOVEL - pool ratemap correlations across conditions

fam_corr = [ pfcorr_tag_all(:,1) ; pfcorr_ntag_all(:,1) ];
fam_corr(isnan(fam_corr)) = [];
nov_corr = [reshape( pfcorr_tag_all(:,2:3), prod( size( pfcorr_tag_all(:,2:3) ) ), 1 ) ; ...
               reshape( pfcorr_ntag_all(:,2:3), prod( size( pfcorr_ntag_all(:,2:3) ) ), 1 ) ];
nov_corr(isnan(nov_corr)) = [];         

figure
set(gcf,'Position', [521   874   409   398])
plot(sort(fam_corr), [ 1:length(fam_corr) ] ./length(fam_corr), 'b' )
hold on
plot(sort(nov_corr), [ 1:length(nov_corr) ] ./length(nov_corr), 'r' )
xlim([-0.5 0.9])
xlabel('Spatial ratemap correlation','fontsize', 16)
ylabel('Cumulative proportion','fontsize', 16)

figure
set(gcf,'Position', [440   377   396   420])
boxplot([fam_corr ; nov_corr], [ones(size(fam_corr)) ; 2.*ones(size(nov_corr))], 'whisker', inf, 'notch', 'on')
ylim([-0.2  0.35])
xticklabels({'FAM', 'NOV'})


%% FAMILIAR  ratemap correlation - all data
tt_fam = pfcorr_tag_all(:,1);
nt_fam = pfcorr_ntag_all(:,1);
figure
set(gcf,'Position', [521   874   409   398])
plot(sort(tt_fam), [ 1:length(tt_fam) ] ./length(tt_fam), 'b' )
hold on
plot(sort(nt_fam), [ 1:length(nt_fam) ] ./length(nt_fam), 'r' )
xlim([-0.4 0.6])
title('Familiar environment')
legend('same birthdate', 'different birthdate', 'location', 'best')

%% FAMILIAR ratemap correlation - individual animals

uniq_animals = unique(animal_tag);
figure
set(gcf,'Position', [779         580        1015         321])
for kp = 1:length( uniq_animals )
    
    subplot(1,length( uniq_animals ), kp)
    tt_fam = pfcorr_tag_all(cellfun(@(x) strcmp(x, uniq_animals{kp}), animal_tag),1);
    nt_fam = pfcorr_ntag_all(cellfun(@(x) strcmp(x, uniq_animals{kp}), animal_ntag),1);
    
    plot(sort(tt_fam), [ 1:length(tt_fam) ] ./length(tt_fam), 'b' )
    hold on
    plot(sort(nt_fam), [ 1:length(nt_fam) ] ./length(nt_fam), 'r' )
    xlim([-0.4 0.6])
    title(sprintf('%s, n=%d SBD', uniq_animals{kp}, length(tt_fam)),  'Interpreter', 'none')
    legend('same birthdate', 'different birthdate', 'location', 'best')
    
end

%% NOVEL ratemap correlation - all data

tt_nov = pfcorr_tag_all(:,2:3); tt_nov = tt_nov(:); tt_nov( isnan( tt_nov ) ) = [];
nt_nov = pfcorr_ntag_all(:, 2:3); nt_nov = nt_nov(:); nt_nov( isnan( nt_nov ) ) = [];
figure
set(gcf,'Position', [521   874   409   398])
plot(sort(tt_nov), [ 1:length(tt_nov) ] ./length(tt_nov), 'b' )
hold on
plot(sort(nt_nov), [ 1:length(nt_nov) ] ./length(nt_nov), 'r' )
xlim([-0.4 0.8])
%title('Novel environment')
xlabel('Spatial ratemap correlation','fontsize', 16)
ylabel('Cumulative proportion','fontsize', 16)

figure
boxplot([ nt_nov ; tt_nov ], [ ones(size( nt_nov )) ; 2.*ones(size(tt_nov)) ], 'notch', 'on', 'whisker', inf)
ylim([-0.2    0.4])

%% NOVEL ratemap correlation - maze type


maze_type = {'linear', 'tmaze', 'openField'};

% Maze type indicator
mazeid_pairs_tag_lin = mazeid_pairs_tag(:,2:3); mazeid_pairs_tag_lin = mazeid_pairs_tag_lin(:);     % tag
mazeid_pairs_ntag_lin = mazeid_pairs_ntag(:,2:3); mazeid_pairs_ntag_lin = mazeid_pairs_ntag_lin(:); % ntag
% Spatial ratemap correlation
pfcorr_tag_all_lin = pfcorr_tag_all(:,2:3); pfcorr_tag_all_lin = pfcorr_tag_all_lin(:);         % tag
pfcorr_ntag_all_lin = pfcorr_ntag_all(:,2:3); pfcorr_ntag_all_lin = pfcorr_ntag_all_lin(:);     % ntag

set(gcf,'Position', [779         580        1015         321])
ranksum_p = [];
for kp = 1:length(maze_type)
    tt_nov = pfcorr_tag_all_lin( cellfun(@(x) ~isempty( regexp(x, maze_type{kp}, 'once') ), mazeid_pairs_tag_lin) );
    nt_nov = pfcorr_ntag_all_lin( cellfun(@(x) ~isempty( regexp(x, maze_type{kp}, 'once') ), mazeid_pairs_ntag_lin) );
    
    subplot(1,length( uniq_animals ), kp)
    plot(sort(tt_nov), [ 1:length(tt_nov) ] ./length(tt_nov), 'b' )
    hold on
    plot(sort(nt_nov), [ 1:length(nt_nov) ] ./length(nt_nov), 'r' )
    xlim([-0.4 0.8])
    
    title(sprintf('%s, n=%d SBD', maze_type{kp}, length( tt_nov )))
    legend('same birthdate', 'different birthdate', 'location', 'best')
    ranksum_p(kp) = ranksum(tt_nov, nt_nov);
end



%% Both fam and nov - run after running both of the above

figure
set(gcf,'Position', [521   874   409   398])
plot(sort(tt_fam), [ 1:length(tt_fam) ] ./length(tt_fam), '--b', 'linewidth', 2 )
hold on
plot(sort(tt_nov), [ 1:length(tt_nov) ] ./length(tt_nov), 'b', 'linewidth', 2 )

plot(sort(nt_fam), [ 1:length(nt_fam) ] ./length(nt_fam), '--r', 'linewidth', 2 )
plot(sort(nt_nov), [ 1:length(nt_nov) ] ./length(nt_nov), 'r', 'linewidth', 2 )
xlim([-0.4 0.6])
title('Novel environment')
legend('SBD familiar', 'SBD novel','DBD familiar','DBD novel', 'location', 'best')


%%

tt_nov_shankdist = pfcorr_tag_shankdist(mazeid_tagtag == 1);
nt_nov_shankdist = pfcorr_ntag_shankdist(mazeid_pairs_ntag == 1);

tagtag_mu = [ nanmean( tt_nov( tt_nov_shankdist == 0 ) );
      nanmean( tt_nov( tt_nov_shankdist == 1  ) ) ;
      nanmean( tt_nov(tt_nov_shankdist == 2 ) ) ];
tagtag_sem = [ sem( tt_nov( tt_nov_shankdist == 0  )) ; ...
               sem( tt_nov( tt_nov_shankdist == 1  ) ) ;
               sem( tt_nov( tt_nov_shankdist == 2  ) ) ];

ntag_mu = [ nanmean( nt_nov( nt_nov_shankdist == 0 ) ) ;
            nanmean( nt_nov( nt_nov_shankdist == 1  ) );
             nanmean( nt_nov( nt_nov_shankdist == 2 ) ) ];
ntag_sem = [ sem( nt_nov( nt_nov_shankdist == 0) ) ; ...
             sem( nt_nov( nt_nov_shankdist == 1) ) ; 
             sem( nt_nov( nt_nov_shankdist == 2) ) ];

         
 errorbar(1:3,tagtag_mu,tagtag_sem, '-s', 'MarkerSize',10,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','CapSize',18')
hold on
errorbar(1:3,ntag_mu,ntag_sem,'-s', 'MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',18)
xticks([1 2 3])
xticklabels({'0', '1', '2'})

set(gcf,'Position',[440 378 404 420])
legend('same birthdate', 'different birthdate')
xlabel('Shank distance', 'fontsize', 20)
ylabel('Ripple cofiring (\rho)', 'fontsize', 20)




%% Firing rates - familiar versus novel environment

% figure
% set(gcf,'Position', [540   377   367   420])
% boxplot([rateTheta_fam_ntag ; rateTheta_fam_tag ; rateTheta_nov_ntag ; rateTheta_nov_tag], ...
%     [ones(size( rateTheta_fam_ntag)) ; 2.*ones(size(rateTheta_fam_tag)) ; 3.*ones(size(rateTheta_nov_ntag)) ; 4.*ones(size(rateTheta_nov_tag))], ...
%       'notch', 'on', 'whisker', inf)
% ylim([0.15 3.5])
% set(gca,'yscale', 'log')
% xticklabels({'FAM', 'FAM', 'NOV', 'NOV'})


a = [rateTheta_fam_ntag ; rateTheta_fam_tag ; rateTheta_nov_ntag ; rateTheta_nov_tag];
b = [ ones(size([rateTheta_fam_ntag ; rateTheta_fam_tag])) ; 2.*ones(size( [rateTheta_nov_ntag ; rateTheta_nov_tag] )) ]; 

figure
set(gcf,'Position', [540   377   367   420])
boxplot(a, b, 'notch', 'on', 'whisker', inf)
ylim([0.2 2])
set(gca,'yscale', 'log')
xticklabels({'FAM', 'NOV'})
ylabel('Firing rate (Hz)', 'fontsize', 16)
ranksum([rateTheta_fam_ntag ; rateTheta_fam_tag ], [rateTheta_nov_ntag ; rateTheta_nov_tag])


%% Ispike - familiar versus novel

% ispike_byMaze_tag = [ispike_byMaze_tag ; ispike_all(pyr_opto_id,:) ];
% ispike_byMaze_ntag = [ispike_byMaze_ntag ; ispike_all(pyr_nopto_id,:)  ];
    
isec_fam_ntag = isec_byMaze_ntag(:,1); %isec_fam_ntag( isnan(isec_fam_ntag) ) = [];
isec_fam_tag = isec_byMaze_tag(:,1);    %isec_fam_tag(isnan(isec_fam_tag)) = [];
isec_nov_ntag = isec_byMaze_ntag(:,2:3); isec_nov_ntag = isec_nov_ntag(:); %isec_nov_ntag(isnan(isec_nov_ntag)) = [];
isec_nov_tag = isec_byMaze_tag(:,2:3);   isec_nov_tag = isec_nov_tag(:);    %isec_nov_tag(isnan(isec_nov_tag)) = [];

% figure
% boxplot([isec_fam_ntag ; isec_fam_tag ; isec_nov_ntag ; isec_nov_tag], ...
%     [ones(size( isec_fam_ntag)) ; 2.*ones(size(isec_fam_tag)) ; 3.*ones(size(isec_nov_ntag)) ; 4.*ones(size(isec_nov_tag))], ...
%             'notch', 'on', 'whisker', inf)
% title('isec')

figure
set(gcf,'Position', [540   377   367   420])
a = [isec_fam_ntag ; isec_fam_tag ; isec_nov_ntag ; isec_nov_tag];
b = [ones(size( [ isec_fam_ntag ;  isec_fam_tag])) ; 2.*ones(size([isec_nov_ntag ; isec_nov_tag]))];
boxplot(a, b, 'notch', 'on', 'whisker', inf)
set(gca,'yscale', 'log')
ylim([0.08    3])
xticklabels({'FAM', 'NOV'})
ylabel('Bits per second', 'fontsize', 18)
ranksum([isec_fam_ntag ; isec_fam_tag ], [isec_nov_ntag ; isec_nov_tag])


%% 
ispike_fam_ntag = ispike_byMaze_ntag(:,1); %ispike_fam_ntag( isnan(ispike_fam_ntag) ) = [];
ispike_fam_tag = ispike_byMaze_tag(:,1);    %ispike_fam_tag(isnan(ispike_fam_tag)) = [];
ispike_nov_ntag = ispike_byMaze_ntag(:,2:3); ispike_nov_ntag = ispike_nov_ntag(:); %isec_nov_ntag(isnan(isec_nov_ntag)) = [];
ispike_nov_tag = ispike_byMaze_tag(:,2:3);   ispike_nov_tag = ispike_nov_tag(:);    %isec_nov_tag(isnan(isec_nov_tag)) = [];
% 
% figure
% boxplot([ispike_fam_ntag ; ispike_fam_tag ; ispike_nov_ntag ; ispike_nov_tag], ...
%     [ones(size( ispike_fam_ntag)) ; 2.*ones(size(ispike_fam_tag)) ; 3.*ones(size(ispike_nov_ntag)) ; 4.*ones(size(ispike_nov_tag))], ...
%             'notch', 'on', 'whisker', inf)
% title('ispike')

% figure
% a = [ispike_fam_ntag ; ispike_fam_tag ; ispike_nov_ntag ; ispike_nov_tag];
% b = [ones(size( [ ispike_fam_ntag ;  ispike_fam_tag])) ; 2.*ones(size([ispike_nov_ntag ; ispike_nov_tag]))];
% boxplot(a, b, 'notch', 'on', 'whisker', inf)
% set(gca,'yscale', 'log')
% %ylim([0.08    3])



%% Firing rates - familiar versus maze type (tmaze, linear, open field)

% Pull out firing rates associated with novel environments

% Tagged
rateTheta_byMaze_tag_lin = rateTheta_byMaze_tag(:,2:3); rateTheta_byMaze_tag_lin = rateTheta_byMaze_tag_lin(:);
mazeid_cell_tag_lin = mazeid_cell_tag(:,2:3); mazeid_cell_tag_lin = mazeid_cell_tag_lin(:);

rateTheta_nov_tag_linear = rateTheta_byMaze_tag_lin( cellfun(@(x) ~isempty( regexp(x, 'linear', 'once') ), mazeid_cell_tag_lin) );
rateTheta_nov_tag_tmaze = rateTheta_byMaze_tag_lin( cellfun(@(x) ~isempty( regexp(x, 'tmaze', 'once') ), mazeid_cell_tag_lin) );
rateTheta_nov_tag_openfield = rateTheta_byMaze_tag_lin( cellfun(@(x) ~isempty( regexp(x, 'open', 'once') ), mazeid_cell_tag_lin) );

% Untagged
rateTheta_byMaze_ntag_lin = rateTheta_byMaze_ntag(:,2:3); rateTheta_byMaze_ntag_lin = rateTheta_byMaze_ntag_lin(:);
mazeid_cell_ntag_lin = mazeid_cell_ntag(:,2:3); mazeid_cell_ntag_lin = mazeid_cell_ntag_lin(:);

rateTheta_nov_ntag_linear = rateTheta_byMaze_ntag_lin( cellfun(@(x) ~isempty( regexp(x, 'linear', 'once') ), mazeid_cell_ntag_lin) );
rateTheta_nov_ntag_tmaze = rateTheta_byMaze_ntag_lin( cellfun(@(x) ~isempty( regexp(x, 'tmaze', 'once') ), mazeid_cell_ntag_lin) );
rateTheta_nov_ntag_openfield = rateTheta_byMaze_ntag_lin( cellfun(@(x) ~isempty( regexp(x, 'open', 'once') ), mazeid_cell_ntag_lin) );

boxplot([rateTheta_fam_ntag ; rateTheta_fam_tag ; rateTheta_nov_ntag_linear ; rateTheta_nov_tag_linear ; ...
    rateTheta_nov_ntag_tmaze ; rateTheta_nov_tag_tmaze ; rateTheta_nov_ntag_openfield ; rateTheta_nov_tag_openfield],...
    [ones(size( rateTheta_fam_ntag)) ; 2.*ones(size(rateTheta_fam_tag)) ; 3.*ones(size(rateTheta_nov_ntag_linear)) ; 4.*ones(size(rateTheta_nov_tag_linear)) ; ...
    5.*ones(size( rateTheta_nov_ntag_tmaze)) ; 6.*ones(size(rateTheta_nov_tag_tmaze)) ; 7.*ones(size(rateTheta_nov_ntag_openfield)) ; 8.*ones(size(rateTheta_nov_tag_openfield))],...
      'notch', 'on', 'whisker', inf)
  
ylim([0.125 3.5])
set(gca,'yscale', 'log')

%% Correlation in firing rate across familiar / novel

loglog( rateTheta_fam_ntag, rateTheta_nov_ntag, '.r')
hold on
loglog( rateTheta_fam_tag, rateTheta_nov_tag, '.b', 'markersize', 10)

xlabel('familiar', 'fontsize', 20)
ylabel('novel', 'fontsize', 20)


lims = [ 0.001 50 ];
xlim(lims)
ylim(lims)
axis square

logx = log10([rateTheta_fam_ntag ; rateTheta_fam_tag]);
logy = log10([rateTheta_nov_ntag ; rateTheta_nov_tag]);
exclude = isnan(logx) | isnan(logy) | isinf(logx) | isinf(logy);
logx( exclude ) = [];
logy( exclude ) = [];

xvec = lims(1):0.025:lims(2);
const = polyfit(logx,logy, 1);
m = const(1);
k = const(2);
hold on
plot( xvec, xvec.^m.*10.^k, 'k' )
xlim(lims)
ylim(lims)

a = [ rateTheta_fam_ntag ; rateTheta_fam_tag ];
b = [ rateTheta_nov_ntag ; rateTheta_nov_tag ];
corr(a, b, 'type', 'Spearman')

%%

save('pfcorr_v6.mat', 'pfcorr_ntag_all', 'pfcorr_ntag_left', 'pfcorr_ntag_right', 'pfcorr_ntag_dist', 'pfcorr_ntag_shankdist', ...
             'theta_ntag_rates', 'thetacorr_ntag', 'id_ntag', 'pfcorr_tagtag_all', 'pfcorr_tagtag_left', 'pfcorr_tagtag_right', ...
          'pfcorr_tagtag_dist', 'pfcorr_tagtag_shankdist', 'theta_tagtag_rates', 'thetacorr_tagtag', 'id_tag', ...
       'ttcorr_tag', 'ttcorr_ntag', 'rateWake_tag', 'rateWake_ntag', 'rateTheta_tag', 'rateTheta_ntag',...
     'id_tag_cells', 'id_ntag_cells', 'ntag_lRat', 'tagtag_lRat', 'animal_ntag', 'animal_tag', 'ntag_wfoverlap', 'tagtag_wfoverlap')