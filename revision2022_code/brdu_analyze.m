


name = 'e15_brdu_minus3_J_section1_z1.czi';

C1_roi_table = readtable([ 'C1-' name '.csv'] );
%C2_roi_table = readtable([ 'C2-' name '.csv'] );

ind_1 = cellfun(@(x) ~isempty(regexp(x, name,'once')), C1_roi_table.Label );
%ind_2 = cellfun(@(x) ~isempty(regexp(x, name,'once')), C2_roi_table.Label );

C1_hist_table = readtable([ 'hist_C1-' name '.csv'] );
% Retrieve the values
c1_vals_all = [];
for kp = 1:length( C1_hist_table.Values )
    c1_vals_all = [ c1_vals_all ; repmat( C1_hist_table.Values(kp), C1_hist_table.Counts(kp), 1) ];
end

%%

histogram(c1_vals_all, 256)
h = xline(C1_roi_table.Mean(ind_1));
for jp = [2 4 5 9 12 13 15]
    h{jp}.Color = 'r';
end

%%

paths_to_czi = [ {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e15_brdu_plus24_A1'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e15_brdu_plus24_B1'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e15_brdu_plus24_C1'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e15_brdu_minus6_N'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e15_brdu_minus6_D'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e15_brdu_minus3_J'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e15_brdu_minus3_I'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e15_brdu_minus3_H'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e14_brdu_zero_G'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e14_brdu_zero_E'} ; ...
                {'Z:\Buzsakilabspace\LabShare\RomanHuszar\DATA\other\brduIUE\confocal\e14_brdu_plus3_L'} ];
all_fils = dir('C:\Users\Roman\Documents\DATA\imageJ_test');
all_fil_names = {all_fils.name}; all_fil_names(1:2) = [];

czi_to_process = [];
for kp = 1:length(paths_to_czi)
   fils = dir( paths_to_czi{kp} ); fils = {fils.name};
   fils(1:2) = [];
   
   for jp = 1:length(fils)
       if strcmp( fils{jp}, 'Thumbs.db'); continue; end
       if ~any( cellfun(@(x) ~isempty(regexp(x, fils{jp}, 'once')), all_fil_names ) )
           czi_to_process = [czi_to_process ; fils(jp)];
       end
       
   end
end


%% Get all the unique files

fils = dir;
fils(1:2) = [];
fils = {fils.name};

roi_fils = find( cellfun(@(x) ~isempty(regexp(x, 'zip','once')), fils) );
z_stacks = {};
for kp = 1:length(roi_fils)
    tmp = strsplit( fils{ roi_fils(kp) }, '.' );
    z_stacks(kp) = {[ tmp{1} '.czi' ]};
end

tdTom_rank = [];
brdu_indicator = [];
brdu_roi_mean = [];
hist_mean = [];
hist_std = [];
timepoint = [];
identifier = [];

% loop through all the stacks and do some basic analysis
for jp = 1:length(z_stacks)
    fprintf('%d/%d\n', jp, length(z_stacks))
    z_stack_name =  z_stacks{jp};
    tmp = strsplit( z_stack_name, '_' );
    C1_roi_table = readtable([ 'C1-' z_stack_name '.csv'] );
    C2_roi_table = readtable([ 'C2-' z_stack_name '.csv'] );
    ind_1 = cellfun(@(x) ~isempty(regexp(x, z_stack_name,'once')), C1_roi_table.Label );
    ind_2 = cellfun(@(x) ~isempty(regexp(x, z_stack_name,'once')), C2_roi_table.Label );
    C1_hist_table = readtable([ 'hist_C1-' z_stack_name '.csv'] );
    
    c1_vals_all = [];
    for kp = 1:length( C1_hist_table.Values )
        c1_vals_all = [ c1_vals_all ; repmat( C1_hist_table.Values(kp), C1_hist_table.Counts(kp), 1) ];
    end
    
    brdu_means = C1_roi_table.Mean(ind_1);
    tdTom_means = C2_roi_table.Mean(ind_2);
    [ tdTom_means_sorted, ip ] = sort(tdTom_means);
    
    rank_perc = [ 1:length( tdTom_means_sorted) ] ./ length(tdTom_means_sorted);
    tdTom_rank = [tdTom_rank ; rank_perc']; 
    brdu_indicator = [ brdu_indicator ; brdu_means(ip) > ( mean(c1_vals_all) + 2*std(c1_vals_all) ) ];
    brdu_roi_mean = [brdu_roi_mean ; brdu_means(ip)];
    hist_mean = [hist_mean ; repmat(mean(c1_vals_all), length( brdu_means ), 1)];
    hist_std = [hist_std ; repmat( std(c1_vals_all), length( brdu_means ), 1)];
    timepoint = [timepoint ; repmat(tmp(3), length( brdu_means ), 1)];
    identifier = [identifier ; repmat({[ tmp{3} '_' tmp{4} ]}, length( brdu_means ), 1)];
end

%%

tp_uniq = [  {'minus6', ''} ; ...
            {'minus3', ''} ; ...
            {'zero','plus3' } ; ...
            {'plus24', ''}];
        
se_doubleLab = [];
mu_doubleLab = [];
mu_tdTomRank = [];
se_tdTomRank = [];
Ncells = [];
Nanimals = [];
% Cycle through timepoints, calculate a boostrapped standard error of the
% mean
for jp = 1:size(tp_uniq,1)
    ind = cellfun(@(x) strcmp(x, tp_uniq{jp,1}) || strcmp(x, tp_uniq{jp,2}) , timepoint);
    bootstat = bootstrp(100,@mean,brdu_indicator(ind));
    se_doubleLab = [se_doubleLab ; std(bootstat)];
    mu_doubleLab = [mu_doubleLab ; mean(brdu_indicator(ind))];
    mu_tdTomRank = [mu_tdTomRank ; median( tdTom_rank( logical( brdu_indicator ) & ind ) )];
    se_tdTomRank = [se_tdTomRank ; sem( tdTom_rank( logical( brdu_indicator ) & ind ) )];
    Ncells = [Ncells ; sum(ind)];
    Nanimals = [Nanimals ; length( unique( identifier(ind) ) )]; 
end


x = 1:length(tp_uniq);

figure
set(gcf,'Position', [836   567   448   537])
bar(x,mu_doubleLab)                
hold on
er = errorbar(x,mu_doubleLab,se_doubleLab,se_doubleLab);  
xticklabels({'-6', '-3', '0h', '+24'})
xlabel('Time from electroporation (h)', 'fontsize', 14)
ylabel('Fraction BrdU+ among electroporated neurons', 'fontsize', 14)

% figure
% set(gcf,'Position',[1466         549         448         537])
% bar(x,mu_tdTomRank)                
% hold on
% er = errorbar(x,mu_tdTomRank,se_tdTomRank,se_tdTomRank);
% xticklabels({'-6', '-3', '0h', '+24'})
% xlabel('Time from electroporation (h)', 'fontsize', 14)
% ylabel('tdTomato fluorescence (normalized)', 'fontsize', 14)