% Extract most possible pathogenic sites(bits) using chi-square test.
% Number of extracted sites is influenced by 'threshold',
% i.e. 'p' in chi-square test.

%%
% Configurations.
% 1: use 0~2 mode. data from genotype;
% 2: use 3 bit mode. data from genotype_3x.
% 0: Pass chi2 test.
% !!! 3 bit mode gets better result than 0~2 mode!!!
encoding_mode = 2;
% save 'possible_pathogenic_idx' to .mat file.
save_to_mat = 1;
global threshold
threshold = 0.001;   % probability threshold in chi2 test.

%%
% Do Chi-square test.
% crosstab(x,y): cross-tabulation of vectors.
% x: phenotype; y: each column of genotype(_3x).
x = phenotype;
possible_pathogenic_idx = [];   %　可能的致病位点index
possibility = [];

%%
if encoding_mode == 1
%------------------------------------------------------------
tic
for i = 1 : num_sites
    [~,~,p] = crosstab(x,genotype(:,i));
    % return params 'table' and 'chi2' are useless so replaced by '~'.
    if p < threshold
        possible_pathogenic_idx = [possible_pathogenic_idx i];
        possibility = [possibility p];
        % display new index immediately so we know how far we've gone...
        fprintf('possible_pathogenic_idx = %d.\n',i);
    end
end
toc
% Takes around 210s, gets:
%       73 sites when p < 0.01;
%        9 sites when p < 0.001.
%------------------------------------------------------------
end

%%
if encoding_mode == 2
%------------------------------------------------------------
tic
for i = 1 : num_sites*3
    [~,~,p] = crosstab(x,genotype_3x(:,i));
    % return params 'table' and 'chi2' are useless so replaced by '~'.
    if p < threshold
        possible_pathogenic_idx = [possible_pathogenic_idx; i];
        possibility = [possibility; p];
        % display new index immediately so we know how far we've gone...
        fprintf('possible_pathogenic_idx = %d.\n',i);
    end
end

% sort
[sorted_psb,sorted_psb_idx_temp] = sort(possibility,'ascend');
% sorted_psb_idx_temp is the index in 'possibility',
% we need to calculate the index in 'genotype_3x'.
temp = possible_pathogenic_idx(sorted_psb_idx_temp);
sorted_psb_idx = temp;
toc
% Takes around 600s. Smallest p is 1e-7. 
% There are 276 bits satisfying p < 0.01;
%            24 bits satisfying p < 0.001;
%             5 bits satisfying p < 0.0001.
%------------------------------------------------------------
end
clear encoding_mode temp sorted_psb_idx_temp

%%
if save_to_mat
    switch threshold
        case 0.01
            str = '_0_01';     % extract 276 sites.
        case 0.001
            str = '_0_001';    % extract  73 sites.
        case 0.0001
            str = '_0_0001';   % extract   5 sites.
    end
    mat_name_str = ['p2_chi2_pathogenic_idx_3x_p' str '.mat'];
    save(mat_name_str,'possible_pathogenic_idx');
    save(mat_name_str,'sorted_psb','sorted_psb_idx','-append');
end
clear table chi2 p save_to_mat x y mat_name_str

%%









