

cellpops = [{'basepaths_e13'}, {'basepaths_e14'}, {'basepaths_e15'}, {'basepaths_e16'} ];

    
%%
animal = 'e16_3m2';
indicator = find( cellfun(@(x) ~isempty(regexp(x, animal, 'once') ), basepaths_all) );

npyr = [];
nopto_pyr = [];
for kp = 1:length(indicator)
    fprintf('%d/%d\n', kp, length(indicator))
    basepath = alterPath( basepaths_all{indicator(kp)}, true );
    basename = bz_BasenameFromBasepath(basepath);
    disp(basename)
    
    load( fullfile( basepath, [basename '.cell_metrics.cellinfo.mat'] ) )
    
    pyr_id = cellfun(@(x) ~isempty(regexp(x, 'Pyramidal', 'once')), cell_metrics.putativeCellType);
    opto_id = ( cell_metrics.optoTag.p_salt < 0.001 ) & cell_metrics.optoTag.h_reliable;
    npyr = [npyr ; sum(pyr_id)];
    nopto_pyr = [nopto_pyr ; sum( pyr_id & opto_id )];
end
%%
fprintf( '%.3f out of %.3f\n',  sum(nopto_pyr), sum(npyr) )
mean( nopto_pyr ./ npyr )
