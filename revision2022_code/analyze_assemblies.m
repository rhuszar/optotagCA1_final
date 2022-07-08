saveMat = false;
% cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV3_OUT')
cd("C:\Users\Roman\Documents\DATA\runAssembly_Vrev_OUT")

% Interested in everything
load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac_rev1.mat")
basepaths = basepaths_all;

pre_path = 'C:\Users\Roman\Documents\DATA\runAssembly_trackV3_OUT_rev';
prepost_path = 'C:\Users\Roman\Documents\DATA\runAssembly_Vrev_OUT';
fils_pre = dir(pre_path);
fils_pre = {fils_pre.name}; fils_pre(1:2) = [];
fils_prepost = dir(prepost_path);
fils_prepost = {fils_prepost.name}; fils_prepost(1:2) = [];

clear assembly_dynamics
assembly_dynamics(4) = struct('nHC',[], 'ripple_ar',[],'nPYR', [],...
                              'tagged',[], 'basename', []);
 
sbd_assembly_membership = [];
dbd_assembly_membership = [];
nsd = 2;    
% Each file is a single session 
% Single out assemblies with at least one high contributing opto tagged
% neuron
exclude_animal = {'e15_13f1', 'e16_3m2'};
for ii = 1:length( basepaths )
    
    basepath = alterPath( basepaths{ii}, true );
    basename = bz_BasenameFromBasepath( basepath );
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    disp(basename)
    % If we have familiar environment exploration, consider assemblies in PRE and POST
    if exist(fullfile(basepath, [basename '.behavior.mat']))
        ind = cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils_prepost);
        this_fil = fils_prepost{ ind };
        prefix = prepost_path;
    % If we have only homecage, consider 2h of data
    else
        ind = cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils_pre);
        this_fil = fils_pre{ ind };
        prefix = pre_path;
    end
    
    % What optotag session are we in?
    tmp = strsplit( this_fil, '_');
    id = tmp{1}; 
    
    % Load assembly output
    v = load( fullfile( prefix, this_fil ) );  
    if strcmp(id, 'e13')
        index = 1;
    elseif strcmp(id, 'e14')
        index = 2;
    elseif strcmp(id, 'e15')
        index = 3;
    elseif strcmp(id, 'e16')
        index = 4;
    end
    
    if ~isfield(v, 'assemblyTemplates')
        assemblyTemplates = v.assemblyTemplates_pre;
        assembly_act_spont = v.assembly_act_pre;
        ripple_int = v.ripple_int_pre;
    else
        assemblyTemplates = v.assemblyTemplates;
        assembly_act_spont = v.assembly_act_spont;
        ripple_int = v.ripple_int;
    end
    pyrInds_opto = arrayfun(@(x) find( x == v.pyrInds), v.pyrInds_opto);
    % Number of "high contributing" neurons (weights > 2SD)
    HC_indicator = zscore(assemblyTemplates) > repmat( nsd, size(assemblyTemplates));
    % Find assemblies with at least one tagged neuron as a high contributor
    taggedAssembly_id = logical( sum( HC_indicator( pyrInds_opto, : ),1 ))';
    
    % SBD assembly membership probability
    if length( pyrInds_opto ) > 1
        sbd_pairs = nchoosek(pyrInds_opto, 2);
        sbd_am = nan( size(sbd_pairs, 1), 1 );
        %size(sbd_pairs, 1)
        for pair_id = 1:size(sbd_pairs, 1)
            sbd_am(pair_id) = sum( HC_indicator(sbd_pairs( pair_id, 1), taggedAssembly_id) & HC_indicator(sbd_pairs( pair_id, 2), taggedAssembly_id) ) ./ sum(taggedAssembly_id);
        end
        sbd_assembly_membership = [sbd_assembly_membership ; sbd_am];
    end
    
    % DBD assembly membership probability
    nonresponsive = setdiff( 1:length(v.pyrInds), pyrInds_opto );
    [m,n] = ndgrid(pyrInds_opto, nonresponsive);
    dbd_pairs = [m(:), n(:)];
    dbd_am = nan( size(dbd_pairs, 1), 1 );
    for pair_id = 1:size(dbd_pairs, 1)
        dbd_am(pair_id) = sum( HC_indicator(dbd_pairs( pair_id, 1), taggedAssembly_id) & HC_indicator(dbd_pairs( pair_id, 2), taggedAssembly_id) ) ./ sum(taggedAssembly_id);
    end
    dbd_assembly_membership = [dbd_assembly_membership ; dbd_am];
    
    % For each assembly, pull out cross correlograms with tagged / untagged
    % nonassembly members
    
    % Activity rates in the different pre sleep states
    ripple_ar = cellfun( @(x) sum(InIntervals(x, ripple_int)), assembly_act_spont ) ./ sum( diff(ripple_int, [], 2) );
    
    assembly_dynamics(index).ripple_ar = [assembly_dynamics(index).ripple_ar ; ripple_ar];
    assembly_dynamics(index).nPYR = [ assembly_dynamics(index).nPYR ; size(assemblyTemplates,1) ];
    assembly_dynamics(index).tagged = logical( [ assembly_dynamics(index).tagged ; taggedAssembly_id ] );
    assembly_dynamics(index).basename = [ assembly_dynamics(index).basename ; repmat({basename}, size(assemblyTemplates,2), 1  ) ];
    
    

end

% save('assemblies_homecage.mat', 'assembly_dynamics')

%%
%{
this is not the right analysis

%}



%%

e13_tag = logical( assembly_dynamics(1).tagged );% & ~cellfun(@(x) ~isempty( regexp(x, 'e13_26m1', 'once') ), assembly_dynamics(1).basename );
e14_tag = logical( assembly_dynamics(2).tagged );
e15_tag = logical( assembly_dynamics(3).tagged ) & ~cellfun(@(x) ~isempty( regexp(x, 'e15_13f1', 'once') ), assembly_dynamics(3).basename ) ;
e16_tag = logical( assembly_dynamics(4).tagged ) & ~cellfun(@(x) ~isempty( regexp(x, 'e16_3m2', 'once') ), assembly_dynamics(4).basename );

mus = [ mean(assembly_dynamics(1).ripple_ar(e13_tag)) mean(assembly_dynamics(2).ripple_ar(e14_tag)) ...
        mean(assembly_dynamics(3).ripple_ar(e15_tag)) mean(assembly_dynamics(4).ripple_ar(e16_tag)) ];
% errhigh = [std(m_e13) std(m_e14) std(m_e15) std(m_e16)];
% errlow  = [std(m_e13) std(m_e14) std(m_e15) std(m_e16)];

errhigh = [std(assembly_dynamics(1).ripple_ar(e13_tag)) std(assembly_dynamics(2).ripple_ar(e14_tag)) ...
           std(assembly_dynamics(3).ripple_ar(e15_tag)) std(assembly_dynamics(4).ripple_ar(e16_tag))];
errlow  = [std(assembly_dynamics(1).ripple_ar(e13_tag)) std(assembly_dynamics(2).ripple_ar(e14_tag)) ...
           std(assembly_dynamics(3).ripple_ar(e15_tag)) std(assembly_dynamics(4).ripple_ar(e16_tag))];

%bar(13.5:1:16.5,diag(mus),'stacked'); 
figure
set(gcf,'Position',[440   377   369   420])
er = errorbar(13.5:1:16.5,mus,errlow);%,errhigh);
hold on
er.Color = [0 0 0]; 
%er.LineStyle = 'none';
ylabel('Ripple activity rate (Hz)', 'fontsize', 14)

hold on
plot( 13.5 + 0.1.*randn( size( assembly_dynamics(1).ripple_ar(e13_tag) ) ), assembly_dynamics(1).ripple_ar(e13_tag), '.r', 'markersize', 10 )
plot( 14.5 + 0.1.*randn( size( assembly_dynamics(2).ripple_ar(e14_tag) ) ), assembly_dynamics(2).ripple_ar(e14_tag), '.b', 'markersize', 10 )
plot( 15.5 + 0.1.*randn( size( assembly_dynamics(3).ripple_ar(e15_tag) ) ), assembly_dynamics(1).ripple_ar(e15_tag), '.m', 'markersize', 10 )
plot( 16.5 + 0.1.*randn( size( assembly_dynamics(4).ripple_ar(e16_tag ) ) ), assembly_dynamics(1).ripple_ar(e16_tag ), '.y', 'markersize', 10 )
xticks([13.5 14.5 15.5 16.5])
% birthdates = [ repmat(13,sum(e13_tag),1) ;
% repmat(14,sum(e14_tag),1) ;
% repmat(15,sum(e15_tag),1) ;
% repmat(16,sum(e16_tag),1)];
a = [ assembly_dynamics(1).ripple_ar(e13_tag) ; assembly_dynamics(2).ripple_ar(e14_tag) ; assembly_dynamics(3).ripple_ar(e15_tag) ; assembly_dynamics(4).ripple_ar(e16_tag ) ];
b = [ ones(sum( e13_tag ), 1 ) ; 2 .* ones(sum( e14_tag ), 1 ) ; 3 .* ones(sum( e15_tag ), 1 ) ; 4 .* ones(sum( e16_tag  ), 1 ) ];
stats = boot_anova1(a,b, 'classical', true);
stats

%% For a non-member neuron - both tagged and untagged - quantify contribution to assembly 
%  expression depending on whether that assembly has at least one tagged neuron or not

% tagged non-members
ccg_tagNmem_all = [];
tagNmem_assemblyIndicator_all = [];
tagNmem_fr_all = [];
tagNmem_weights_all = [];
% untagged non-members
ccg_ntagNmem_all = [];
ntagNmem_assemblyIndicator_all = [];
ntagNmem_fr_all = [];
ntagNmem_weights_all = [];
 
% ccg_tagNmem_nonrip_all = [];
% ccg_ntagNmem_nonrip_all = [];

saveMat = false;
% cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV3_OUT')
cd("C:\Users\Roman\Documents\DATA\runAssembly_Vrev_OUT")

% Interested in everything
load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac_rev1.mat")
basepaths = basepaths_all;

pre_path = 'C:\Users\Roman\Documents\DATA\runAssembly_trackV3_OUT_rev';
prepost_path = 'C:\Users\Roman\Documents\DATA\runAssembly_Vrev_OUT';
fils_pre = dir(pre_path);
fils_pre = {fils_pre.name}; fils_pre(1:2) = [];
fils_prepost = dir(prepost_path);
fils_prepost = {fils_prepost.name}; fils_prepost(1:2) = [];


nsd = 2;
% Each file is a single session 
% Single out assemblies with at least one high contributing opto tagged
% neuron
exclude_animal = {'e15_13f1', 'e16_3m2'};
for ii = 1:length( basepaths )
    
    
    basepath = alterPath( basepaths{ii}, true );
    basename = bz_BasenameFromBasepath( basepath );
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    disp(basename)
    if exist(fullfile(basepath, [basename '.behavior.mat']))
        ind = cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils_prepost);
        this_fil = fils_prepost{ ind };
        prefix = prepost_path;
    else
        ind = cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils_pre);
        this_fil = fils_pre{ ind };
        prefix = pre_path;
    end
    
    % What optotag session are we in?
    tmp = strsplit( this_fil, '_');
    id = tmp{1}; 
    
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    % Load assembly output
    v = load( fullfile( prefix, this_fil ) );  
    if strcmp(id, 'e13')
        index = 1;
    elseif strcmp(id, 'e14')
        index = 2;
    elseif strcmp(id, 'e15')
        index = 3;
    elseif strcmp(id, 'e16')
        index = 4;
    end
    
    if ~isfield(v, 'assemblyTemplates')
        load(fullfile(basepath, [basename '.pulses.events.mat']))
        assemblyTemplates = v.assemblyTemplates_pre;
        assembly_act_spont = v.assembly_act_pre;
        % Can do this analysis inside and outside ripples
        ripple_int = v.ripple_int_pre;
        nonpulse_period = SubtractIntervals( [ 0 v.pre_end ], pulses.intsPeriods );
        nonpulse_nonrip_period = SubtractIntervals(nonpulse_period, ripple_int);
    else
        assemblyTemplates = v.assemblyTemplates;
        assembly_act_spont = v.assembly_act_spont;
        % Can do this analysis inside and outside ripples
        ripple_int = v.ripple_int;
        nonpulse_period = v.nonpulse_period;
        nonpulse_nonrip_period = SubtractIntervals(nonpulse_period, ripple_int);
    end
    %pyrInds_opto = logical(sum( cell2mat( [ cellfun(@(x) x == v.pyrInds, num2cell(v.pyrInds_opto), 'UniformOutput', false) ]') ));
     pyrInds_opto = arrayfun(@(x) find( x == v.pyrInds ), v.pyrInds_opto);
     tmp = false( size( v.pyrInds ) ); tmp(pyrInds_opto) = true;
     pyrInds_opto = tmp;
    % Number of "high contributing" neurons (weights > 2SD)
    HC_indicator = zscore(assemblyTemplates) > repmat( nsd, size(assemblyTemplates));
    % Find assemblies with at least one tagged neuron as a high contributor
    taggedAssembly_id = logical( sum( HC_indicator( pyrInds_opto, : ),1 ))';
    
    ccg_ntagNmem = [];
    ccg_tagNmem = [];
    ntagNmem_assemblyIndicator = [];
    tagNmem_assemblyIndicator = [];
    ntagNmem_fr = [];
    tagNmem_fr = [];
    tagNmem_weights = [];
    ntagNmem_weights = [];
    
    % ccg_tagNmem_nonrip = [];
    % ccg_ntagNmem_nonrip = [];
    
    % For a given assembly - it's tagged or untagged ; a nonmember can be
    % tagged or untagged ; store all of this
    for assemblyN = 1:size(assemblyTemplates,2)
        % tagged nonmembers
        % this is utterly incorrect !!!!! we are pulling out members from
        % nonmembers - we want only nonm
        tagNmem_id = ~HC_indicator(:,assemblyN) & pyrInds_opto';
        tagNmem_assemblyIndicator = [tagNmem_assemblyIndicator ; repmat( taggedAssembly_id(assemblyN), sum(tagNmem_id), 1) ];
        tagNmem_weights = [tagNmem_weights ; assemblyTemplates(tagNmem_id, assemblyN)];
        for tagNmem = v.pyrInds( tagNmem_id )
            currSpk = spikes.times{tagNmem}( InIntervals(spikes.times{tagNmem}, nonpulse_period) );
            % currSpk_nonrip = spikes.times{tagNmem}( InIntervals(spikes.times{tagNmem}, nonpulse_nonrip_period) );
            % all
            ccg = CCG([ currSpk ; assembly_act_spont{assemblyN} ], [ ones(size( currSpk )) ; 2*ones(size( assembly_act_spont{assemblyN} )) ], 'binsize', .002, 'duration', .2 );
            ccg_tagNmem = [ ccg_tagNmem ; ccg(:,1,2)' ./ length(currSpk) ];
            tagNmem_fr = [tagNmem_fr ; length(currSpk) ./ sum( diff( nonpulse_period' ) )];
            % exclude ripples
            % ccg = CCG([ currSpk_nonrip ; assembly_act_spont{assemblyN} ], [ ones(size( currSpk_nonrip )) ; 2*ones(size( assembly_act_spont{assemblyN} )) ], 'binsize', .002, 'duration', .2 );
            % ccg_tagNmem_nonrip = [ ccg_tagNmem_nonrip ; ccg(:,1,2)' ./ length(currSpk_nonrip) ];
        end
        % untagged nonmembers
        ntagNmem_id = ~HC_indicator(:,assemblyN) & ~pyrInds_opto';
        ntagNmem_assemblyIndicator = [ntagNmem_assemblyIndicator ; repmat( taggedAssembly_id(assemblyN), sum(ntagNmem_id), 1) ];
        ntagNmem_weights = [ntagNmem_weights ; assemblyTemplates(ntagNmem_id, assemblyN)];
        for ntagNmem = v.pyrInds( ntagNmem_id )
            currSpk = spikes.times{ntagNmem}( InIntervals(spikes.times{ntagNmem}, nonpulse_period) );
            % currSpk_nonrip = spikes.times{ntagNmem}( InIntervals(spikes.times{ntagNmem}, nonpulse_nonrip_period) );
            % all 
            ccg = CCG([ currSpk ; assembly_act_spont{assemblyN} ], [ ones(size( currSpk )) ; 2*ones(size( assembly_act_spont{assemblyN} )) ], 'binsize', .002, 'duration', .2 );
            ccg_ntagNmem = [ ccg_ntagNmem ; ccg(:,1,2)' ./ length(currSpk) ];
            ntagNmem_fr = [ntagNmem_fr ; length(currSpk) ./ sum( diff( nonpulse_period' ) )];
            % exclude ripples
            % ccg = CCG([ currSpk_nonrip ; assembly_act_spont{assemblyN} ], [ ones(size( currSpk_nonrip )) ; 2*ones(size( assembly_act_spont{assemblyN} )) ], 'binsize', .002, 'duration', .2 );
            % ccg_ntagNmem_nonrip = [ ccg_ntagNmem_nonrip ; ccg(:,1,2)' ./ length(currSpk_nonrip) ];
        end
    end
    ccg_tagNmem_all = [ccg_tagNmem_all ; ccg_tagNmem];
    ccg_ntagNmem_all = [ccg_ntagNmem_all ; ccg_ntagNmem];
    ntagNmem_assemblyIndicator_all = [ntagNmem_assemblyIndicator_all ; ntagNmem_assemblyIndicator];
    tagNmem_assemblyIndicator_all = [tagNmem_assemblyIndicator_all ; tagNmem_assemblyIndicator];
    tagNmem_fr_all = [tagNmem_fr_all ; tagNmem_fr];
    ntagNmem_fr_all = [ntagNmem_fr_all ; ntagNmem_fr];
    ntagNmem_weights_all = [ntagNmem_weights_all ; ntagNmem_weights];
    tagNmem_weights_all = [tagNmem_weights_all ; tagNmem_weights];
    
%     ccg_tagNmem_nonrip_all = [ccg_tagNmem_nonrip ; ccg_tagNmem_nonrip];
%     ccg_ntagNmem_nonrip_all = [ccg_ntagNmem_nonrip ; ccg_ntagNmem_nonrip];


end

%%

plot( nanmean( ccg_tagNmem( logical(tagNmem_assemblyIndicator), : ) ) )
hold on
plot( nanmean( ccg_tagNmem( ~logical(tagNmem_assemblyIndicator), : ) ) )

plot( nanmean( ccg_ntagNmem( logical(ntagNmem_assemblyIndicator), : ) ) )
plot( nanmean( ccg_ntagNmem( ~logical(ntagNmem_assemblyIndicator), : ) ) )

%%


saveMat = true;
% Interested in everything
%load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac.mat")
basepaths = basepaths_all;
% basepaths = basepaths_beh;

% Interested in behavior
% load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_beh.mat')
% basepaths = basepaths_beh;

pre_path = 'C:\Users\Roman\Documents\DATA\runAssembly_trackV3_OUT_rev';
prepost_path = 'C:\Users\Roman\Documents\DATA\runAssembly_Vrev_OUT';
fils_pre = dir(pre_path);
fils_pre = {fils_pre.name}; fils_pre(1:2) = [];
fils_prepost = dir(prepost_path);
fils_prepost = {fils_prepost.name}; fils_prepost(1:2) = [];
                   
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
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    
    disp(basename)
    if exist(fullfile(basepath, [basename '.behavior.mat']))
        ind = cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils_prepost);
        this_fil = fils_prepost{ ind };
        prefix = prepost_path;
    else
        ind = cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils_pre);
        this_fil = fils_pre{ ind };
        prefix = pre_path;
    end
    
    % What optotag session are we in?
    tmp = strsplit( this_fil, '_');
    id = tmp{1}; 
 
    % Load assembly output
    v = load( fullfile( prefix, this_fil ) );
    if isfield(v, 'assemblyTemplates_pre')
        assemblyTemplates = v.assemblyTemplates_pre;
    else
        assemblyTemplates = v.assemblyTemplates;
    end
    
    pyrInds_opto = arrayfun(@(x) find( x == v.pyrInds) ,  v.pyrInds_opto);
    % Number of "high contributing" neurons (weights > 2SD)
    HC_indicator = zscore(assemblyTemplates) > repmat( nsd, size(assemblyTemplates));
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

save('assemblies_homecage.mat', 'tagtag_dist', 'tagtag_shankdist', 'ripcorr_tagtag', 'ripple_tagtag_rates', 'id_tag', 'animal_tag', ...
    'basenames_tag', 'ripple_ntag_rates', 'ntag_dist', 'ntag_shankdist', 'ripcorr_ntag', 'id_ntag', 'animal_ntag', 'basenames_ntag', '-append');

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

%%
a = [ripcorr_ntag(e13_ntag) ; ripcorr_ntag(e14_ntag) ; ripcorr_ntag(e15_ntag) ; ripcorr_ntag(e16_ntag)];
b = [ones(sum( e13_ntag), 1) ; 2.*ones(sum( e14_ntag), 1) ; 3.*ones(sum( e15_ntag), 1) ; 4.*ones(sum( e16_ntag), 1)];
stats = boot_anova1(a,b,'classical', true);


%% Activity k neurons - pool across animals
% Express each neuron's probability of expression as a zscore with respect
% to the control distribution
%cd('C:\Users\Roman\Documents\DATA\estimate_k_OUT')

nex = 20;

fils = dir;
fils = {fils.name}; fils(1:2) = [];

prob_zval_all = {};
k_all = {};
for kp = 12:length(fils)
    fils{kp}
    v = load(fils{kp});
    k_nonzero = find( v.k_tag > 0 );
    prob_zval = [];
    k = k_nonzero+1;
    for jp = 1:length(k_nonzero)
        % Control distribution - log transform
%         d = log10( v.k_rand(:,k_nonzero(jp)) ) ;
%         d( isinf(d) ) = -30;
        d = v.k_rand(:,k_nonzero(jp)) ;
        mu = mean( d ); sd = std( d );
         prob_zval(jp) = (  v.k_tag( k_nonzero(jp) ) - mu ) ./ sd;
 %       prob_zval(jp) = (  log10( v.k_tag( k_nonzero(jp) ) ) - mu ) ./ sd;
    end
    k( isnan( prob_zval ) ) = [];
    prob_zval( isnan( prob_zval ) ) = [];

    prob_zval_all{kp} = prob_zval;
    k_all{kp} = k;
    
    plot(2:length(v.k_tag)+1, v.k_tag, '-b' )
    hold on
    h = {};
    for jp = 1:nex
        h{jp} = plot(2:length(v.k_tag)+1, v.k_rand(randi(1000),:), '-k' );
        h{jp}.Color(4) = 0.2;
    end
    shg
    uiwait
    

end

prob_zval_mat = nan( length(prob_zval_all), max( cellfun(@length, prob_zval_all )) );
for kp = 1:length(prob_zval_all)
    prob_zval_mat(kp, 1:length( prob_zval_all{kp} ) ) = prob_zval_all{kp};
end

%%

se = nansem(prob_zval_mat);
mu = nanmean(prob_zval_mat);
nanind = isnan(mu) | isinf(mu) | isnan(se) | isinf(se);
mu(nanind) = [];
se(nanind) = [];

xvec = 2:length(mu)+1;
f = [ mu+se flip(mu-se) ];
xf = [ xvec flip( xvec ) ];

fill(xf, f, 'b', 'FaceAlpha',0.2, 'linestyle', 'none');
hold on
plot( xvec, mu, '-b', 'linewidth', 2 )
yline(0)

%%
for kp = 1:size(k_tag_all,2)
    figure
    a = k_tag_all_mat(kp,:);
    b = k_tag_rand_mat(kp,:);
    ind = ~isnan(a) & ~isnan(b);
    semilogy( a(ind) )
    hold on
    semilogy( b(ind) )
    uiwait
    clc
end

%%

ind = ~( sum( ~isnan( k_tag_all_mat ), 2) == 1 );
semilogy( nanmean( k_tag_all_mat(ind,:) ) )
hold on
semilogy( nanmean( k_tag_rand_mat(ind,:) ) )
%%

nan( length(k_tag_rand),  max( cellfun(@length, k_tag_all ) ) )

%%

cd('C:\Users\Roman\Documents\DATA\runAssembly_hardcoded_OUT')
fils = dir;
fils = {fils.name}; fils(1:2) = [];
ripz=[];
animal_id = {};
for kp = 1:length(fils)
    v = load( fils{kp} );
    tmp = strsplit( v.basename, '_' );
%     histogram( v.ripple_rate_null, 10, 'Normalization', 'probability','linestyle', 'none' )
%     title(sprintf('Session ID: %s', v.basename), 'Interpreter', 'none')
%     xlabel('Assembly expression rate (Hz) in SPW-Rs')
%     xline( mean( v.ripple_rate_tag ) )
    ripz(kp) = ( v.ripple_rate_tag - mean( v.ripple_rate_null ) ) ./ std( v.ripple_rate_null );
    animal_id(kp) = {[ tmp{1} ]};
%     uiwait
end

%%

histogram( ripz, 13, 'normalization', 'probability', 'linestyle', 'none' )
xlabel('Assembly expression rate (Z)')
xline(0)








