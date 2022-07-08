

inds = toremove;
%inds = [[4 6 9] ; [2 20 26] ]';
im = []
for kp = 1:length(inds)
    im = [ im ; squeeze( cc_mat( inds(kp,1), inds(kp,2), :) )' ] ;
end



imagesc( reshape( cc_mat([ 5 6 7], [3  5 8],:), 9, 121) )

%%
% 53, 73, 87
%
Nlines = 30;

write_files = false;

outpath = '/Users/romanhuszar/Google Drive/buzz_workspace/optotag/analyze/figs/may2021/ccgs';


for kp = [53, 63, 64, 73, 87]
%for kp = 1:length(basepaths_all)
    
    basepath = alterPath(basepaths_all{kp}, true);
    basename = bz_BasenameFromBasepath(basepath);
    ccg_path = fullfile(basepath, 'ccgs.mat');
    if ~exist(ccg_path)
        continue
    end
    
    v = load(ccg_path);
    
    tvec=-.3:0.002:0.3;
    tails = [ v.ccg_ntag(:,1:75) v.ccg_ntag(:,end-75:end) ];
    ccf_ntag_z = ( v.ccg_ntag - nanmean(tails,2) ) ./ std(tails,[],2);
    to_flip = find( rand(size( ccf_ntag_z,1 ), 1) < 0.5 );
    ccf_ntag_z(to_flip,:) = flip( ccf_ntag_z(to_flip,:), 2);
    
    figure
    set(gcf,'Position', [440   378   419   420])
    hold on
    h = {};
    cc = 1;
    
    if size(v.ccg_tagtag,1) > Nlines
        inds_ntag = randi(size(v.ccg_ntag,1), 1, Nlines);
        inds_tag = randi(size(v.ccg_tagtag,1), 1, Nlines);
    else
        inds_ntag = randi(size(v.ccg_ntag,1), 1, size(v.ccg_tagtag,1));
        inds_tag = 1:size(v.ccg_tagtag,1);
    end
    
    for jp = inds_ntag
        h{cc} = plot(tvec,smooth( ccf_ntag_z(jp,:) ), 'r');
        h{cc}.Color(4) = 0.15;
        cc = cc+1;
    end
    line1 = plot( tvec,smooth( nanmean( ccf_ntag_z ) ), 'r', 'linewidth', 3 );
    
    
    tails = [ v.ccg_tagtag(:,1:75) v.ccg_tagtag(:,end-75:end) ];
    ccf_tagtag_z = ( v.ccg_tagtag - nanmean(tails,2) ) ./ std(tails,[],2);
    to_flip = find( rand(size( ccf_tagtag_z,1 ), 1) < 0.5 );
    ccf_tagtag_z(to_flip,:) = flip( ccf_tagtag_z(to_flip,:), 2);
    h1 = {};
    cc = 1;
    for jp = inds_tag
        h1{cc} = plot(tvec,smooth( ccf_tagtag_z(jp,:) ), 'b');
        h1{cc}.Color(4) = 0.15;
        cc = cc+1;
    end
    line2 = plot( tvec,smooth( nanmean( ccf_tagtag_z ) ), 'b','linewidth', 3 );
    line_h = [line1 line2];
    legend(line_h, 'different birthdate', 'same birthdate')
    xlim([-0.2 0.2])
    ylim([-3.5 18.5])
    xlabel('Time lag (s)', 'fontsize', 20)
    ylabel('Cross correlation (z-scored)', 'fontsize', 20)
    basename_dot = strrep( basename, '_', '.');
    title(sprintf('%s', basename_dot), 'fontsize', 20)
    %xline(0)
    % uiwait
    if write_files
        export_fig( fullfile(outpath, sprintf( '%s_ccgs', basename_dot) ), '-png')
        close all
    else
        uiwait
    end
end


%%

sig = []; n = [];
for kp = 1:length(basepaths_all)
    
    basepath = alterPath(basepaths_all{kp}, true);
    basename = bz_BasenameFromBasepath(basepath);
    ccg_path = fullfile(basepath, 'ccgs.mat');
    if ~exist(ccg_path)
        continue
    end
    
    v = load(ccg_path);
    
    tvec=-.3:0.002:0.3;
    tails = [ v.ccg_ntag(:,1:75) v.ccg_ntag(:,end-75:end) ];
    ccg_ntag_z = ( v.ccg_ntag - nanmean(tails,2) ) ./ std(tails,[],2);
    
    tails = [ v.ccg_tagtag(:,1:75) v.ccg_tagtag(:,end-75:end) ];
    ccg_tagtag_z = ( v.ccg_tagtag - nanmean(tails,2) ) ./ std(tails,[],2);
    
%     ntag = nanmean(v.ccg_ntag(:,tvec >= -.025 & tvec <= .025), 2);
%     tag = nanmean(v.ccg_tagtag(:,tvec >= -.025 & tvec <= .025), 2);
    ntag = nanmean(ccg_ntag_z(:,tvec >= -.05 & tvec <= .05), 2);
    tag = nanmean(ccg_tagtag_z(:,tvec >= -.05 & tvec <= .05), 2);
    %sig = [ sig ; ttest2(ntag, tag(randi(size(tag,1 ), 1,length(ntag))) ) ];
    ff = [];
    for jp = 1:1000
        ff(jp) = mean( ntag(randi(size(tag,1 ), 1,length(tag))) );
    end
    sig = [sig ; sum( mean(tag) > ff ) / 1000];
    
    n = [n ; size(ccg_tagtag_z,1)];
    

end


%% Here we take average z scored ccgs in all sessions 

Nlines = 30;

write_files = false;

outpath = '/Users/romanhuszar/Google Drive/buzz_workspace/optotag/analyze/figs/may2021/ccgs';
ntag = [];
tag = [];
tvec=-.3:0.002:0.3;
%for kp = [53, 63, 64, 73, 87]
for kp = 1:length(basepaths_all)
    
    basepath = alterPath(basepaths_all{kp}, true);
    basename = bz_BasenameFromBasepath(basepath);
    ccg_path = fullfile(basepath, 'ccgs.mat');
    if ~exist(ccg_path)
        continue
    end
    
    v = load(ccg_path);
    
    tails = [ v.ccg_ntag(:,1:75) v.ccg_ntag(:,end-75:end) ];
    ccf_ntag_z = ( v.ccg_ntag - nanmean(tails,2) ) ./ std(tails,[],2);
    to_flip = find( rand(size( ccf_ntag_z,1 ), 1) < 0.5 );
    ccf_ntag_z(to_flip,:) = flip( ccf_ntag_z(to_flip,:), 2);
    
    ntag = [ntag; nanmean( ccf_ntag_z ) ];
    
    tails = [ v.ccg_tagtag(:,1:75) v.ccg_tagtag(:,end-75:end) ];
    ccf_tagtag_z = ( v.ccg_tagtag - nanmean(tails,2) ) ./ std(tails,[],2);
    to_flip = find( rand(size( ccf_tagtag_z,1 ), 1) < 0.5 );
    ccf_tagtag_z(to_flip,:) = flip( ccf_tagtag_z(to_flip,:), 2);

    tag = [tag ; nanmean( ccf_tagtag_z )];
    

end


%% Look at the z-score CCG by animal

outpath = '/Users/romanhuszar/Google Drive/buzz_workspace/optotag/analyze/figs/may2021/ccgs';
ntag = []; animal_ntag = [];
tag = []; animal_tag = [];
rip_med = [];
tvec=-.3:0.002:0.3;
%for kp = [53, 63, 64, 73, 87]
exclude_animal = {'e16_3m2', 'e15_13f1'};
for kp = 1:length(basepaths_all)
    
    fprintf('%d/%d\n', kp, length(basepaths_all))
    basepath = alterPath(basepaths_all{kp}, true);
    basename = bz_BasenameFromBasepath(basepath);
    
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), exclude_animal) )
        continue
    end
    
    ccg_path = fullfile(basepath, 'ccgs.mat');
    load( fullfile(basepath, [basename '.ripples.events.mat']) )
    if ~exist(ccg_path)
        continue
    end
    
    rip_med = [ rip_med ; median( ripples.data.duration ) ];
    
    ff = strsplit( basename, '_' );
    an = [ ff{1} '_' ff{2} ];
    
    v = load(ccg_path);
    if isempty( v.ccg_tagtag )
        continue
    end
    
    tails = [ v.ccg_ntag(:,1:75) v.ccg_ntag(:,end-75:end) ];
    ccf_ntag_z = ( v.ccg_ntag - nanmean(tails,2) ) ./ std(tails,[],2);
    to_flip = find( rand(size( ccf_ntag_z,1 ), 1) < 0.5 );
    ccf_ntag_z(to_flip,:) = flip( ccf_ntag_z(to_flip,:), 2);
    
    ntag = [ntag; ccf_ntag_z ];
    animal_ntag = [animal_ntag ; repmat({an},size( ccf_ntag_z, 1),1)];
    
    tails = [ v.ccg_tagtag(:,1:75) v.ccg_tagtag(:,end-75:end) ];
    ccf_tagtag_z = ( v.ccg_tagtag - nanmean(tails,2) ) ./ std(tails,[],2);
    to_flip = find( rand(size( ccf_tagtag_z,1 ), 1) < 0.5 );
    ccf_tagtag_z(to_flip,:) = flip( ccf_tagtag_z(to_flip,:), 2);

    tag = [tag ; ccf_tagtag_z];
    animal_tag = [animal_tag ; repmat({an},size( ccf_tagtag_z, 1),1)];
    

end

%%
uniq_tag = unique( animal_tag );
animal_ccgResid = [];
for kp = 1:length( uniq_tag )
    t_ind = cellfun(@(x) ~isempty( regexp(x, uniq_tag{kp}, 'once') ), animal_tag);
    nt_ind =  cellfun(@(x) ~isempty( regexp(x, uniq_tag{kp}, 'once') ), animal_ntag);
    
    animal_ccgResid = [animal_ccgResid ; smooth( nanmean( tag(t_ind, : ) ) - nanmean( ntag(nt_ind, : ) ) )' ];

end
animal_ccgResid(2:3,:) = [];
%%

imagesc(tvec, 1:size(animal_ccgResid,1), animal_ccgResid)
%set(gca,'colorscale', 'log')
set(gca,'clim', [-0.5 2])
colormap(cool)
xline([-mean( rip_med )/2 mean( rip_med )/2])


%%
close all
n = 1;
t_ind = find( cellfun(@(x) ~isempty( regexp(x, uniq_tag{n}, 'once') ), animal_tag) );
nt_ind =  find( cellfun(@(x) ~isempty( regexp(x, uniq_tag{n}, 'once') ), animal_ntag) );

nlines = 10;

t_ind1 = t_ind( randi(length(t_ind), 1,Nlines) );
nt_ind1 = nt_ind( randi(length(nt_ind), 1, Nlines) );

%%
figure
hold on


h1 = {};
cc = 1;
for jp = nt_ind1'
    h1{cc} = plot(tvec,smooth( ntag(jp,:) ), 'r');
    h1{cc}.Color(4) = 0.15;
    cc = cc+1;
end
line1 = plot( tvec,smooth( nanmean( ntag(nt_ind,:) ) ), 'r','linewidth', 3 );

h2 = {};
cc = 1;
for jp = t_ind1'
    h2{cc} = plot(tvec,smooth( tag(jp,:) ), 'b');
    h2{cc}.Color(4) = 0.15;
    cc = cc+1;
end

line2 = plot( tvec,smooth( nanmean( tag(t_ind,:) ) ), 'b','linewidth', 3 );
line_h = [line1 line2];
legend(line_h, 'different birthdate', 'same birthdate')

xlim([-0.2 0.2])
ylim([-3 14])
xlabel('Time lag (s)', 'fontsize', 20)
ylabel('Cross correlation (z-scored)', 'fontsize', 20)
title(sprintf('%s', uniq_tag{n}), 'fontsize', 20)
%%

figure
plot( tvec,smooth( nanmean( tag(t_ind,:)) - nanmean( ntag(nt_ind,:) ) ), 'k','linewidth', 3 )
xlabel('Time lag (s)', 'fontsize', 20)
ylabel('Cross correlation difference', 'fontsize', 20)
title('Same birthdate - different birthdate','fontsize', 20)
xlim([-0.2 0.2])

