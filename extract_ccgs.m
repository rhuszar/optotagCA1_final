
%{

Get the CCGs associated with tagged / untagged pyramidal neurons

%}


% same birthdate versus different birthdate

%load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
% load('C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac.mat')

% Do we wish to load data from local hard drive? 
isLocal = true;

binsize = 0.002;
duration = 0.6;

randLabel = 0;
while randLabel <= 1   
    randLabel = 2000*rand;
end

for bp = 1:length(basepaths_all)
    
    % Distance estimates
    ntag_dist = [];
    ntag_shankdist = [];
    % Convergence index
    ccg_ntag = [];

    % Distance estimates
    tagtag_dist = [];
    tagtag_shankdist = [];
    % Convergence index
    ccg_tagtag = [];
    
    fprintf('%d/%d\n', bp, length(basepaths_all))
    basepath = alterPath( basepaths_all{bp}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.pulses.events.mat']))
    %load(fullfile(basepath, [basename '.thetaCycles.events.mat']))
    %load(fullfile(basepath, [basename '.SleepState.states.mat']))
    
    T = max( cellfun(@max, spikes.times ) ) - sum( diff( pulses.intsPeriods' ) );
    
    load(fullfile(basepath, 'chanMap.mat'))
    v.chanMap = [xcoords ; ycoords];
    v.kcoords = kcoords;
    

    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    int_id = cellfun(@(x) ~isempty(regexp(x, 'Interneuron', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id & pyr_id;
    
    
    pyr_id = find(pyr_id);
    pyr_opto_id = arrayfun(@(x) find(x==pyr_id),  find(pyr_opto_id) );
    
%     if length(pyr_opto_id) >= 5
%         continue
%     end
    
    dist_mat = nan(length(pyr_id), length(pyr_id));
    shankdist_mat = nan(length(pyr_id), length(pyr_id));
    indicator_mat = zeros(length(pyr_id), length(pyr_id));
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
   
    
     % Label duplicate pairs (e.g., if neuron 5 and 10 are tagged,  
    % pairs (5,10) and (10,5) will have the same correlation)
    indicator_mat(pyr_opto_id,:) = 1;
    if length(pyr_opto_id)>1
        toremove = nchoosek( pyr_opto_id, 2 );
        
        for ii = 1:size(toremove,1)
            ip = pyr_id( toremove(ii,1) ); jp = pyr_id( toremove(ii,2) );
            n1 = spikes.times{ip}( ~InIntervals( spikes.times{ip}, pulses.intsPeriods ) );
            n2 = spikes.times{jp}( ~InIntervals( spikes.times{jp}, pulses.intsPeriods ) ); fr_n2 = length(n2) / T;
            cc = CCG([n1 ; n2], [ones(size(n1)) ; 2*ones(size(n2))], 'binsize', binsize, 'duration', duration );
            ccg_tagtag = [ ccg_tagtag ; ( ( cc(:,1,2)' / length(n1) ) / binsize ) ];
        end

        
        %id_tag = [ id_tag ; repmat(tmp(1), size( toremove, 1), 1)];
        tagtag_dist = [ tagtag_dist ; dist_mat( sub2ind( size(indicator_mat), toremove(:,1), toremove(:,2) ) ) ];
        tagtag_shankdist = [ tagtag_shankdist ; shankdist_mat( sub2ind( size(indicator_mat), toremove(:,1), toremove(:,2) ) ) ];
        
        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
        indicator_mat( sub2ind( size(indicator_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        indicator_mat( sub2ind( size(indicator_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;
        
        dist_mat( sub2ind( size(indicator_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        dist_mat( sub2ind( size(indicator_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;
        shankdist_mat( sub2ind( size(indicator_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        shankdist_mat( sub2ind( size(indicator_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;

    end
    indicator_mat(1:length(pyr_id)+1:end) = randLabel;
    dist_mat(1:length(pyr_id)+1:end) = randLabel;
    shankdist_mat(1:length(pyr_id)+1:end) = randLabel;
    
    distvec = dist_mat(pyr_opto_id,:);                  distvec = distvec(:);
    shankdistvec = shankdist_mat(pyr_opto_id,:);        shankdistvec = shankdistvec(:);
    
    distvec( distvec == randLabel ) = [];
    shankdistvec( shankdistvec == randLabel ) = [];
    
    ntag_dist = [ ntag_dist ; distvec ];
    ntag_shankdist = [ ntag_shankdist ; shankdistvec ];
    
    [ip,kp] = ind2sub(size(indicator_mat), find( indicator_mat == 1 ) );
    % Compute the CCGs between tagged and untagged cells
    for ii = 1:length(ip)
        n1_ind = pyr_id( ip(ii) ); n2_ind = pyr_id( kp(ii) );
        n1 = spikes.times{n1_ind}( ~InIntervals( spikes.times{n1_ind}, pulses.intsPeriods ) );
        n2 = spikes.times{n2_ind}( ~InIntervals( spikes.times{n2_ind}, pulses.intsPeriods ) ); fr_n2 = length(n2) / T;
        cc = CCG([n1 ; n2], [ones(size(n1)) ; 2*ones(size(n2))], 'binsize', binsize, 'duration', duration );
        ccg_ntag = [ ccg_ntag ; ( ( cc(:,1,2)' / length(n1) ) / binsize ) ];
    end

    save( fullfile(basepath, 'ccgs.mat' ), 'ccg_ntag', 'ntag_shankdist', 'ntag_dist',...
                                     'ccg_tagtag', 'tagtag_dist', 'tagtag_shankdist' )
    
end