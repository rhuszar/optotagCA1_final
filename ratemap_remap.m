

%%
cc_all = [];

animal_exclude = [{'e16_3m2'} {'e15_13f1'}];
for kp = 1:length(basepaths_beh)
    fprintf('%d/%d\n',kp,length(basepaths_beh))
    basepath = alterPath( basepaths_beh{kp}, true );
    basename = bz_BasenameFromBasepath(basepath);
    
    if any( cellfun(@(x) ~isempty(regexp(basename, x, 'once')), animal_exclude) )
        continue
    end
    
    load(fullfile(basepath, [basename '.cell_metrics.cellinfo.mat']))
    try
        load(fullfile(basepath, [basename '.firingMapsAvg_v2.cellinfo.mat']))
    catch
        load(fullfile(basepath, [basename '.firingMapsAvg_multSess.cellinfo.mat']))
    end
    pyr_id_log = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    pyrRatemaps = firingMaps(1).rateMaps(pyr_id_log);
    
    % pyrRatemaps = firingMaps.rateMaps( setdiff( firingMaps.UID, firingMaps.pyr_ind ) );
    % Combine the ratemaps
    ratemaps_left = [];
    ratemaps_right = [];
    for kp = 1:length(pyrRatemaps)
        ratemaps_left = [ ratemaps_left ; pyrRatemaps{kp}{1}  ];
        ratemaps_right = [ ratemaps_right ; pyrRatemaps{kp}{2}  ];
    end
    cc = nan(1,205);
    for kp = 1:205
        t = corr( [ ratemaps_left(:,kp) ratemaps_right(:,kp)] );
        cc(kp) = t(1,2);
    end
    cc_all = [cc_all ; cc];
        
%     [coeff,score,latent] = pca(ratemaps');
%     plot3(score(1:205,1), score(1:205,2), score(1:205,3), 'r')
%     hold on
%     plot3(score(206:end,1), score(206:end,2), score(206:end,3), 'b')
%     plot3(score(74,1), score(74,2), score(74,3), 'xk', 'markersize', 15)
%     plot3(score(20,1), score(20,2), score(20,3), 'xm', 'markersize', 15)
%     grid on
%     uiwait
end

%%

subplot(3,1,1)
imagesc( cc_all )
colormap(cool)
xlabel('Linearized position (cm)', 'fontsize', 16)
ylabel('Session #', 'fontsize', 16)
axis xy
colorbar
title('Population vector decorrelation')

subplot(3,1,2)


xv = 1:205;
fnt = [nanmean( cc_all )+std( cc_all )  flip( nanmean( cc_all )-std( cc_all ) )];
fill( [xv  flip(xv)], fnt, 'k', 'facealpha', 0.1, 'linestyle', 'none' )
hold on
%plot(tt_mu, '-b', 'linewidth', 1)
plot( nanmean( cc_all ), 'k', 'linewidth', 2 )
xlim([1 205])
colorbar

xline([74 111 185 222])

subplot(3,1,3)
cc_stem = cc_all(:,1:74); cc_stem = cc_stem(:);
cc_turn1 = cc_all(:,75:111); cc_turn1 = cc_turn1(:);
cc_arm = cc_all(:,112:185); cc_arm = cc_arm(:);
cc_turn2 = cc_all(:,186:205); cc_turn2 = cc_turn2(:);
boxplot([cc_stem ; cc_turn1 ; cc_arm ; cc_turn2], [ones(size(cc_stem)) ; 2*ones(size(cc_turn1)) ; 3*ones(size(cc_arm)) ; 4*ones(size(cc_turn2))],...
    'notch','on','whisker',3,'labels',{'stem','turn1', 'arm', 'turn2'})


%%

