



% same birthdate versus different birthdate
% This code should be part of the get_corrBeh_v2 pipeline > 

% load('/Users/romanhuszar/Google Drive/buzz_workspace/optotag/basepaths_mac.mat')
load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac_rev1.mat")
% load('C:\Users\Roman\Google Drive\buzz_workspace\optotag\basepaths_mac.mat')

% Do we wish to load data from local hard drive? 
isLocal = true;

% Distance estimates
ntag_dist = [];
ntag_shankdist = [];
% Convergence index
converge_ntag = [];
id_ntag = [];
basename_ntag = [];

% Distance estimates
tagtag_dist = [];
tagtag_shankdist = [];
% Convergence index
converge_tagtag = [];
id_tag = [];
basename_tag = [];

randLabel = 2000*rand;

for ip = 1:length(basepaths_all)
    
    fprintf('%d/%d\n', ip, length(basepaths_all))
    basepath = alterPath( basepaths_all{ip}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.mono_res.cellinfo.mat']))
    %load(fullfile(basepath, [basename '.pulses.events.mat']))
    %load(fullfile(basepath, [basename '.thetaCycles.events.mat']))
    %load(fullfile(basepath, [basename '.SleepState.states.mat']))
    
    load(fullfile(basepath, 'chanMap.mat'))
    v.chanMap = [xcoords ; ycoords];
    v.kcoords = kcoords;
    

    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    int_id = cellfun(@(x) ~isempty(regexp(x, 'Interneuron', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id & pyr_id;
    
    
    pyr_id = find(pyr_id);
    pyr_opto_id = arrayfun(@(x) find(x==pyr_id),  find(pyr_opto_id) );
    
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
    
     % Label duplicate pairs (e.g., if neuron 5 and 10 are tagged,  
    % pairs (5,10) and (10,5) will have the same correlation)
    if length(pyr_opto_id)>1
        toremove = nchoosek( pyr_opto_id, 2 );
        
        converge_tagtag = [ converge_tagtag ; convergence_mat( sub2ind( size(convergence_mat), toremove(:,1), toremove(:,2) ) ) ];
        id_tag = [ id_tag ; repmat(tmp(1), size( toremove, 1), 1)];
        basename_tag = [ basename_tag ; repmat({basename}, size( toremove, 1), 1) ];
        tagtag_dist = [ tagtag_dist ; dist_mat( sub2ind( size(convergence_mat), toremove(:,1), toremove(:,2) ) ) ];
        tagtag_shankdist = [ tagtag_shankdist ; shankdist_mat( sub2ind( size(convergence_mat), toremove(:,1), toremove(:,2) ) ) ];
        
        % Label all optotagged pairs - i.e., (3,6) as well as (6,3)
        convergence_mat( sub2ind( size(convergence_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        convergence_mat( sub2ind( size(convergence_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;
        
        dist_mat( sub2ind( size(convergence_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        dist_mat( sub2ind( size(convergence_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;
        shankdist_mat( sub2ind( size(convergence_mat), toremove(:,1), toremove(:,2) ) ) = randLabel;
        shankdist_mat( sub2ind( size(convergence_mat), toremove(:,2), toremove(:,1) ) ) = randLabel;

    end
    
    dist_mat(1:length(pyr_id)+1:end) = randLabel;
    shankdist_mat(1:length(pyr_id)+1:end) = randLabel;
    convergence_mat(1:length(pyr_id)+1:end) = randLabel;
    
    convergence_vec = convergence_mat(pyr_opto_id,:);   convergence_vec = convergence_vec(:); 
    distvec = dist_mat(pyr_opto_id,:);                  distvec = distvec(:);
    shankdistvec = shankdist_mat(pyr_opto_id,:);        shankdistvec = shankdistvec(:);
    
    convergence_vec( convergence_vec == randLabel ) = [];
    distvec( distvec == randLabel ) = [];
    shankdistvec( shankdistvec == randLabel ) = [];
    
    ntag_dist = [ ntag_dist ; distvec ];
    ntag_shankdist = [ ntag_shankdist ; shankdistvec ];
    converge_ntag = [ converge_ntag ; convergence_vec ];
    id_ntag = [ id_ntag ; repmat(tmp(1), length(convergence_vec), 1)];
    basename_ntag = [ basename_ntag ; repmat({basename},length(convergence_vec), 1) ];
    
end

%%


figure
set(gcf,'Position',[440 378 396 420])
plot(sort( converge_ntag), [ 1:length(converge_ntag) ] ./ length(converge_ntag), '-r');
hold on
plot(sort( converge_tagtag), [ 1:length(converge_tagtag) ] ./ length(converge_tagtag), '-b')
xlim( [-0.1 0.15] )
xlabel('Ripple cofiring (\rho)', 'fontsize', 20)
ylabel('Cumulative proportion', 'fontsize', 20)

figure
boxplot([converge_ntag ; converge_tagtag], [ones(size(converge_ntag)) ; 2*ones(size(converge_tagtag))],'notch','on', 'whisker',inf,...
        'labels',{'different birthdate','same birthdate'})
ylim([-0.05 0.5])
ylabel('Convergence onto I cells', 'fontsize', 20)

%%

% e13_ntag = cellfun(@(x) strcmp(x, 'e13'), id_ntag );  e14_ntag = cellfun(@(x) strcmp(x, 'e14'), id_ntag ); e15_ntag = cellfun(@(x) strcmp(x, 'e15'), id_ntag ); e16_ntag = cellfun(@(x) strcmp(x, 'e16'), id_ntag );
% e13_tag = cellfun(@(x) strcmp(x, 'e13'), id_tag );  e14_tag = cellfun(@(x) strcmp(x, 'e14'), id_tag ); e15_tag = cellfun(@(x) strcmp(x, 'e15'), id_tag ); e16_tag = cellfun(@(x) strcmp(x, 'e16'), id_tag );
ntag_median = [nanmedian(converge_ntag(e13_ntag)) ; nanmedian(converge_ntag(e14_ntag)) ; nanmedian(converge_ntag(e15_ntag)) ; nanmedian(converge_ntag(e16_ntag)) ];
figure
boxplot([converge_tagtag(e13_tag)-ntag_median(1) ; converge_tagtag(e14_tag)-ntag_median(2) ; converge_tagtag(e15_tag)-ntag_median(3) ; converge_tagtag(e16_tag)-ntag_median(4)], ...
         [ones(sum( e13_tag ),1) ; 2*ones(sum( e14_tag ),1) ; 3*ones(sum( e15_tag ),1) ; 4*ones(sum( e16_tag ),1)],'notch','on', 'whisker',inf,...
        'labels',{'E13.5','E14.5', 'E15.5', 'E16.5'})
    
ylim([-0.25 0.4]) 
ylabel('Convergence onto I cells (corrected)', 'fontsize', 20)

%%

tagtag_mu = [ nanmean( converge_tagtag( tagtag_shankdist == 0 ) );
      nanmean( converge_tagtag( tagtag_shankdist == 1  ) ) ;
      nanmean( converge_tagtag(tagtag_shankdist == 2 ) ) ];
tagtag_sem = [ sem( converge_tagtag( tagtag_shankdist == 0  )) ; ...
               sem( converge_tagtag( tagtag_shankdist == 1  ) ) ;
               sem( converge_tagtag( tagtag_shankdist == 2  ) ) ];

ntag_mu = [ nanmean( converge_ntag( ntag_shankdist == 0 ) ) ;
            nanmean( converge_ntag( ntag_shankdist == 1  ) );
             nanmean( converge_ntag( ntag_shankdist == 2 ) ) ];
ntag_sem = [ sem( converge_ntag( ntag_shankdist == 0) ) ; ...
             sem( converge_ntag( ntag_shankdist == 1) ) ; 
             sem( converge_ntag( ntag_shankdist == 2) ) ];

         
 errorbar(1:3,tagtag_mu,tagtag_sem, '-s', 'MarkerSize',10,...
    'MarkerEdgeColor','blue','MarkerFaceColor','blue','CapSize',18')
hold on
errorbar(1:3,ntag_mu,ntag_sem,'-s', 'MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red','CapSize',18)
xticks([1 2 3])
xticklabels({'0', '1', '2'})

set(gcf,'Position',[440 378 404 420])
legend('PRE pairs of the same birthdate', 'PRE pairs of different birthdate')
xlabel('Shank distance', 'fontsize', 20)
ylabel('Convergence onto I cells', 'fontsize', 20)

%% Extract the rate modulation in interneurons modultion
pulse_excl = 0.02;
isLocal = true;
for ip = 1:length(basepaths_all)
    
    fprintf('%d/%d\n', ip, length(basepaths_all))
    basepath = alterPath( basepaths_all{ip}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    % Check if we've extracted interneuron modulation already
    if isfield( cell_metrics.optoTag, 'int_ratemod' ) && ...
            isfield( cell_metrics.optoTag, 'int_ratemod_pval' )
        continue
    end
    
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    % load(fullfile(basepath, [basename '.mono_res.cellinfo.mat']))
    load(fullfile(basepath, [basename '.pulses.events.mat']))
    
    pulses_sh_ints = pulses.timestamps( diff( pulses.timestamps' ) <0.005,: );
    post_pulse_sh_ints = [ pulses_sh_ints(:,2) pulses_sh_ints(:,2) + pulse_excl ];
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

id_pretag = [];
spktrans_intdriven = [];
spktrans_splitISI_intdriven = [];
isLocal = true;

for ip = 1:length(basepaths_all)
    
    fprintf('%d/%d\n', ip, length(basepaths_all))
    basepath = alterPath( basepaths_all{ip}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.mono_res.cellinfo.mat']))
    load(fullfile(basepath, [basename '.pulses.events.mat']))
    
    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    int_id = cellfun(@(x) ~isempty(regexp(x, 'Interneuron', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id & pyr_id;
    int_opto_id = opto_id & int_id;
    
    driven_int_ids = find( int_id & cell_metrics.optoTag.int_ratemod > 0 & cell_metrics.optoTag.int_ratemod_pval < 0.05  );
    % driven_int_ids = find(int_opto_id);
    % Compute splitISI spike transmission of connections between optotagged
    % units and driven postsynaptic interneurons
    presyn_tagged = cellfun(@(x) x == mono_res.pyr2int(:,1), num2cell( find( pyr_opto_id) ), 'UniformOutput', false );
    % Loop over connections where presynaptic neuron is tagged
    for jj = find( sum( cell2mat( presyn_tagged ), 2) )'
        pyr2int = mono_res.pyr2int(jj,:);
        if ~any( pyr2int(2) == driven_int_ids )
            continue
        end
        % Presynaptic spikes outside pulses
        pre = spikes.times{pyr2int(1)}( ~InIntervals( spikes.times{pyr2int(1)}, pulses.intsPeriods ) );
        post = spikes.times{pyr2int(2)};
        ses = ShortTermCCG(pre, post,.0008,.2,'time', [ logspace(log10(5),log10(1000),15)/1000 inf]);
        spktrans_splitISI_intdriven = [ spktrans_splitISI_intdriven ; ses.trans ];
        
        % Also get the average spike transmission probability
%         [ccg,~] = CCG([pre ; post], [ones(size(pre)) ; 2*ones(size(post))], 'binsize', .0008, 'duration', 0.2);
        gammaSpk = find( diff( pre ) >= 0.0156 &  diff( pre ) <= 0.0484 )+1;
        [ccg,~] = CCG([pre(gammaSpk) ; post], [ones(length(gammaSpk),1) ; 2*ones(size(post))], 'binsize', .0008, 'duration', 0.2);
        
        [ccgSpkTrans,~,~,~] = GetTransProb(ccg(:,1,2),length(gammaSpk),.0008);
        spktrans_intdriven = [spktrans_intdriven ; ccgSpkTrans]; 
        id_pretag = [ id_pretag ; tmp(1) ];
    end
    
end

%%


e13_id = cellfun(@(x) strcmp(x, 'e13'), id_pretag); 
e14_id = cellfun(@(x) strcmp(x, 'e14'), id_pretag);
e15_id = cellfun(@(x) strcmp(x, 'e15'), id_pretag); 
e16_id = cellfun(@(x) strcmp(x, 'e16'), id_pretag);

isi = [ logspace(log10(5),log10(1000),15)/1000 ];


figure
hold on

spktrans_e13 = spktrans_splitISI_intdriven(e13_id,:);
mu_e13 = nanmean(spktrans_e13 );
se = nanstd(spktrans_e13) ./ sqrt( sum( ~isnan( spktrans_e13 ) ) );
fe13 = [mu_e13+se flip( mu_e13-se )];

spktrans_e14 = spktrans_splitISI_intdriven(e14_id,:);
mu_e14 = nanmean(spktrans_e14);
se = nanstd(spktrans_e14) ./ sqrt( sum( ~isnan( spktrans_e14 ) ) );
fe14 = [mu_e14+se flip( mu_e14-se )];

spktrans_e15 = spktrans_splitISI_intdriven(e15_id,:);
mu_e15 = nanmean(spktrans_e15);
se = nanstd(spktrans_e15) ./ sqrt( sum( ~isnan( spktrans_e15 ) ) );
fe15 = [mu_e15+se flip( mu_e15-se )];

spktrans_e16 =spktrans_splitISI_intdriven(e16_id,:);
mu_e16 = nanmean(spktrans_e16);
se = nanstd(spktrans_e16) ./ sqrt( sum( ~isnan( spktrans_e16 ) ) );
fe16 = [mu_e16+se flip( mu_e16-se )];



plot(isi, mu_e13, 'Color', brown, 'linewidth', 2)
plot(isi, mu_e14, 'r', 'linewidth', 2)
plot(isi, mu_e15, 'b', 'linewidth', 2)
plot(isi, mu_e16, 'k', 'linewidth', 2)
fill( [isi flip(isi)], fe13, brown, 'facealpha', 0.1, 'linestyle', 'none' )
fill( [isi flip(isi)], fe14, 'r', 'facealpha', 0.1, 'linestyle', 'none' )
fill( [isi flip(isi)], fe15, 'b', 'facealpha', 0.1, 'linestyle', 'none' )
fill( [isi flip(isi)], fe16, 'k', 'facealpha', 0.1, 'linestyle', 'none' )

legend('E13.5', 'E14.5', 'E15.5', 'E16.5')


set(gca,'xscale', 'log')
set(gca,'yscale', 'log')

xlabel('presynaptic ISI (s)', 'fontsize', 20)
ylabel('Spike transmission probability', 'fontsize', 20)

%%
e13_spktrans = spktrans_splitISI_intdriven(e13_id,5);
e14_spktrans = spktrans_splitISI_intdriven(e14_id,5);
e15_spktrans = spktrans_splitISI_intdriven(e15_id,5);
e16_spktrans = spktrans_splitISI_intdriven(e16_id,5);

figure
boxplot([e13_spktrans ; e14_spktrans ; e15_spktrans ; e16_spktrans], ...
         [ones(size(e13_spktrans)) ; 2*ones(size(e14_spktrans)) ; 3*ones(size(e15_spktrans)) ; 4*ones(size(e16_spktrans))],'notch','on', 'whisker',inf,...
        'labels',{'E13.5','E14.5', 'E15.5', 'E16.5'})
    
ylabel('Spike transmission probability', 'fontsize', 20)
xticks([1 2 3 4])
xticklabels([{'E13', 'E14', 'E15', 'E16'}])
set(gca,'yscale', 'log')
ylim([-.01 0.05]) 

%% Spike transmission, split by birthdate

e13_id = cellfun(@(x) strcmp(x, 'e13'), id_pretag); 
e14_id = cellfun(@(x) strcmp(x, 'e14'), id_pretag);
e15_id = cellfun(@(x) strcmp(x, 'e15'), id_pretag); 
e16_id = cellfun(@(x) strcmp(x, 'e16'), id_pretag);

tagtag_mu = [nanmean(spktrans_intdriven(e13_id)) ; nanmean(spktrans_intdriven(e14_id)) ; nanmean(spktrans_intdriven(e15_id)) ; nanmean(spktrans_intdriven(e16_id)) ; ];
tagtag_sem = [sem(spktrans_intdriven(e13_id)) ; sem(spktrans_intdriven(e14_id)) ; sem(spktrans_intdriven(e15_id)) ; sem(spktrans_intdriven(e16_id)); ];


errorbar(1:4,tagtag_mu,tagtag_sem, '-sk', 'MarkerSize',10,...
    'MarkerEdgeColor','k','MarkerFaceColor','k','CapSize',18')

set(gcf,'Position',[440 378 404 420])


ylabel('Spike transmission probability', 'fontsize', 20)
xlabel('Birthdate','fontsize', 20)
xticks([1 2 3 4])
xticklabels([{'E13', 'E14', 'E15', 'E16'}])


figure
boxplot([spktrans_intdriven(e13_id) ; spktrans_intdriven(e14_id) ; spktrans_intdriven(e15_id) ; spktrans_intdriven(e16_id)], ...
         [ones(sum( e13_id ),1) ; 2*ones(sum( e14_id ),1) ; 3*ones(sum( e15_id ),1) ; 4*ones(sum( e16_id ),1)],'notch','on', 'whisker',inf,...
        'labels',{'E13.5','E14.5', 'E15.5', 'E16.5'})
    
ylim([0 0.04]) 
ylabel('Spike transmission probability', 'fontsize', 20)
xticks([1 2 3 4])
xticklabels([{'E13', 'E14', 'E15', 'E16'}])

%% This is for visualization purposes
%  Pick out a single driven INT - all we care is that more birthdate
%  neurons project to it than non birthdated neurons

for ip = 57:-1:1%length(basepaths_all):-1:1
    
    fprintf('%d/%d\n', ip, length(basepaths_all))
    basepath = alterPath( basepaths_all{ip}, isLocal);
    basename = bz_BasenameFromBasepath(basepath);
    tmp = strsplit( basename, '_' );
    
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    
    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    int_id = cellfun(@(x) ~isempty(regexp(x, 'Interneuron', 'once')), cell_metrics.putativeCellType);
    pyr_opto_id = opto_id & pyr_id;
    int_opto_id = opto_id & int_id;
    
    if sum(pyr_opto_id) < 5
        continue
    end
    
    load(fullfile(basepath, [basename '.spikes.cellinfo.mat']))
    load(fullfile(basepath, [basename '.mono_res.cellinfo.mat']))
    load(fullfile(basepath, [basename '.pulses.events.mat']))

    driven_int_ids = find( int_id & cell_metrics.optoTag.int_ratemod > 0 & cell_metrics.optoTag.int_ratemod_pval < 0.05  );
    % driven_int_ids = find(int_opto_id);
    % Compute splitISI spike transmission of connections between optotagged
    % units and driven postsynaptic interneurons
    % presyn_tagged = cellfun(@(x) x == mono_res.pyr2int(:,1), num2cell( find( pyr_opto_id) ), 'UniformOutput', false );
    
    % Loop over connections where presynaptic neuron is tagged
    for jj = 1:length( driven_int_ids )
        
        inttag = intersect( mono_res.pyr2int(mono_res.pyr2int(:,2) == driven_int_ids(jj),1), find( pyr_opto_id ));
        intntag = intersect( mono_res.pyr2int(mono_res.pyr2int(:,2) == driven_int_ids(jj),1), find( ~pyr_opto_id ));
        
        if ( length(inttag) <= 1 && length(intntag) <= 1) || ( length(inttag) <= length(intntag))
            continue
        end
        
        %% store UIDs of the relevant units
        % pyr_tagged = [ 11    14    15];
        % pyr_untagged = 33
        % int_driven that's targeted by tagged and untagged = 31
        % int that is driven by untagged alone = [22 35]
        post = spikes.times{driven_int_ids(jj)};
        
        figure
        curr = 1;
        for kp = [5,7,8]
             subplot(1,3, curr)
             pre = spikes.times{inttag(kp)}( ~InIntervals( spikes.times{ inttag(kp) }, pulses.intsPeriods ) );
            % Also get the average spike transmission probability
            [ccg,~] = CCG([pre ; post], [ones(size(pre)) ; 2*ones(size(post))], 'binsize', .0008, 'duration', 0.2);
            [ccgSpkTrans,ccgProb,~,~] = GetTransProb(ccg(:,1,2),length(pre),.0008);
            plot(ccgProb)
            ylim([-0.001 0.012])
            curr = curr+1;
        end
        
        figure
        for kp = 1
             pre = spikes.times{intntag(kp)}( ~InIntervals( spikes.times{ intntag(kp) }, pulses.intsPeriods ) );
            % Also get the average spike transmission probability
            [ccg,~] = CCG([pre ; post], [ones(size(pre)) ; 2*ones(size(post))], 'binsize', .0008, 'duration', 0.2);
            [ccgSpkTrans,ccgProb,~,~] = GetTransProb(ccg(:,1,2),length(pre),.0008);
            plot(ccgProb)
            ylim([-0.001 0.012])
        end
        
        %%
        pyr_tagged = [ 11    14    15];
        % Which other interneurons does my untagged cell target;
        int_untagged = [ 22    34    35    51    64    67    68    69 ];
        
        % Which among these ints receives no connections from my tagged
        % PYRs
        mono_res.pyr2int(:,2)
        %%
        
        %% store UIDs of the relevant units
        % pyr_tagged = [ 11    14    15];
        % pyr_untagged = 33
        % int_driven that's targeted by tagged and untagged = 31
        % int that is driven by untagged alone = [22 35]
        
        UIDs = [11 14 15  33 31 35 ];
        spks = []; ids = [];
        for kp = 1:length(UIDs)
            spks = [spks ; spikes.times{UIDs(kp)}];
            ids = [ids ; kp.*ones(size(spikes.times{UIDs(kp)}))];
        end
        
        [ccgs, t] = CCG(spks, ids, 'binsize', 0.0008, 'duration', 0.4);
        
        curr = 1;
        for kp = 1:length(UIDs)
            for jp = kp:length(UIDs)
                subplot(6,6,curr)
                if kp ~= jp
                    [ccgSpkTrans,ccgProb,~,~] = GetTransProb(ccgs(:,kp,jp) ,length( spikes.times{UIDs(kp)} ),.0008 );
                    plot(t,ccgProb )
                    ylim([-0.001 0.012])
                else
                    plot(t,ccgs(:,kp,jp) )
                end
                curr = curr+1;
                xlim([-.1 0.1])
            end
            curr = curr+kp;
            disp('foo')
        end
        
        %%
        % Presynaptic spikes outside pulses
        disp('something')
    end
    
end

