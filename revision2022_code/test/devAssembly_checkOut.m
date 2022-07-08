

outpath = 'C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\final_code\revision2022_code\test';
path = 'C:\Users\Roman\Documents\DATA\devAssembly_simVrev_null';
cd(path)
load("C:\Users\Roman\Google Drive\buzz_workspace\optotag\analyze\final_code\revision2022_code\test\fast_inds_to_process.mat")

fils = dir;
fils = {fils.name};
fils(1:2) = []; 

res = [];
for kp = 1:length(fils)
    curr_fil = fils{kp};
    path2ind = fullfile(path, curr_fil);
    tmp = strsplit(curr_fil, '_');
    index = str2num( tmp{2} );
    nfils = dir(path2ind);
    nfils = {nfils.name}; nfils(1:2) = [];
    res = [res ; index length( nfils )];
end

inds_to_process = setdiff( 1:420, res(res(:,2) == 100, 1 )  );
save( fullfile(outpath, 'fast_inds_to_process.mat'), 'inds_to_process' ) 
%%
fils = dir;
fils = {fils.name};
fils(1:2) = [];
err = [];
for kp = 1:length(fils)
    curr_fil = fils{kp};
    v = load(curr_fil);
    err = [ err ; -sum( v.lp1 ) ];
end
histogram(err)

%%

fast_inds_to_process = find( A>=9 );

