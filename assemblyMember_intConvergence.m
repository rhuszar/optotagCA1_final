
%
% High contributing - look at all pairs  

saveMat = true;
cd('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/DATA/runAssembly_trackV3_OUT')
load('/Users/romanhuszar/Google Drive (not syncing)/buzz_workspace/optotag/basepaths_mac.mat')
% cd("C:\Users\Roman\Google Drive\buzz_workspace\optotag\DATA\runAssembly_trackV3_OUT")
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

% Distance estimates
tagtag_dist = [];
tagtag_shankdist = [];
% Correlations
ripcorr_tagtag = [];
% Ripple rates
ripple_tagtag_rates = [];
id_tag = [];

members = [];
members_sbd = [];
convergence = [];
convergence_sbd = [];

randLabel = 2000*rand;

% Each file is a single session 
% Single out assemblies with at least one high contributing opto tagged
% neuron
for ii = 1:length( basepaths )
    
    basepath = alterPath( basepaths{ii}, true );
    basename = bz_BasenameFromBasepath( basepath );
    
    disp(basename)
    fils_id = find( cellfun(@(x)~isempty(regexp(x, basename, 'once')), fils) );
    if ~any( fils_id )
        continue
    end
    
    % Load assembly output
    v = load( fils{fils_id} );
    % What optotag session are we in?
    tmp = strsplit( fils{fils_id}, '_');    
    
    % Number of "high contributing" neurons (weights > 2SD)
    HC_indicator = zscore(v.assemblyTemplates_pre) > repmat( nsd, size(v.assemblyTemplates_pre));

    load(fullfile(basepath, [basename '.mono_res.cellinfo.mat']))

    load(fullfile(basepath, 'chanMap.mat'))
    v.chanMap = [xcoords ; ycoords];
    v.kcoords = kcoords;
    
    % For each pair of cells, calculate the number of assemblies they are a
    % member, as well as their convergence index
    convergence_mat = nan(length(v.pyrInds ),length(v.pyrInds ));
    assemblyMember_mat = nan(length(v.pyrInds ),length(v.pyrInds ));
    for ip = 1:length( v.pyrInds )
        for jp = 1:length( v.pyrInds )
            assemblyMember_mat(ip,jp) = sum( HC_indicator(ip,:) & HC_indicator(jp,:) );
            n1 = unique( mono_res.pyr2int( mono_res.pyr2int(:,1) == v.pyrInds(ip), 2) );
            n2 = unique( mono_res.pyr2int( mono_res.pyr2int(:,1) == v.pyrInds(jp), 2) );
            convergence_mat(ip,jp) = length( intersect(n1, n2) ) / length( unique( [n1 ; n2] ) );
        end
    end
    % upper triangle indicator
    indicator_triu = triu( true(size(assemblyMember_mat)), 1);
    % same birthdate indicator
    sbd = double( arrayfun(@(x) any( x == v.pyrInds_opto), v.pyrInds) );
    sbd_indicator = logical( sbd' * sbd );
    
    members = [ members ; assemblyMember_mat(indicator_triu & ~sbd_indicator) ];
    convergence = [ convergence ; convergence_mat(indicator_triu & ~sbd_indicator) ];
    members_sbd = [ members_sbd ; assemblyMember_mat(indicator_triu & sbd_indicator) ];
    convergence_sbd = [ convergence_sbd ; convergence_mat(indicator_triu & sbd_indicator) ];
    
    
end

%%

conv_mu = arrayfun(@(x) nanmean( convergence( x == members) ), 0:2);
conv_ci = arrayfun(@(x) 1.96*sem( convergence( x == members) ), 0:2);
conv_sbd_mu = arrayfun(@(x) nanmean( convergence_sbd( x == members_sbd) ), 0:1);
conv_sbd_ci = arrayfun(@(x) 1.96*sem( convergence_sbd( x == members_sbd) ), 0:1);

figure
set(gcf,'Position',[440   377   364   420])
f = [conv_mu+conv_ci flip(conv_mu-conv_ci)];
fill([0:2 flip(0:2)], f, [7 7 7]/8,'linestyle', 'none')
hold on; plot(0:2, conv_mu, '-k', 'linewidth', 1.5)

f = [conv_sbd_mu+conv_sbd_ci flip(conv_sbd_mu-conv_sbd_ci)];
fill([0:1 flip(0:1)], f, 'b','facealpha', 0.2,'linestyle', 'none')
hold on; plot(0:1, conv_sbd_mu, '-b','linewidth', 1.5)
xticks([0 1 2])
%legend('other', 'SBD')

xlabel('Shared assemblies per neuron pair', 'fontsize', 18)
ylabel('Interneuron convergence per neuron pair', 'fontsize', 18)

