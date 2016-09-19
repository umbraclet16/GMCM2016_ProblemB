
% Synthesize possible pathogenic bits obtained from both chi-square test
% and infinite norm sorting, and give the final results.

%%
num_sites = 9455;
% Only evaluate top 10 bits in each result.
max = 10;

%%
for i = 1 : 4
    % Load files.
    switch i
        case 1
            filename = 'p2_chi2_pathogenic_idx_3x_p_0_01.mat';
        case 2
            filename = 'p2_chi2_pathogenic_idx_3x_p_0_001.mat';
        case 3
            filename = 'p2_chi2_pathogenic_idx_3x_p_0_0001.mat';
        case 4
            % We only need top 10, so using top 50 file is enough.
            % Other files(100/200/300) give same results.
            filename = 'p2_inf_norm_pathogenic_idx_3x_50.mat';
    end
    
    load(filename);
    fprintf('Loaded file: %s.\n',filename)
    
    % Deal with variable name inconsistency from 2 methods...
    % When using chi-square test, the variable we need is called
    % 'sorted_psb_idx';
    % When using infinite norm sorting, the variable we need is called
    % 'sorted_dif_bit_idx'.
    if i > 3
        sorted_psb_idx = sorted_dif_bit_idx;
    end
    
    %%
    n = max;
    if length(sorted_psb_idx) < n
        % Results less than 10.
        n = length(sorted_psb_idx);
    end
    
    % Bits to sites.
    result = floor((sorted_psb_idx(1:n) - ones(n,1)) / 3) + 1;
    
    
    fprintf('Results from %s:\n',filename)
    for i = 1 : n
        fprintf('%d\t',result(i));
    end
    fprintf('\n');
    

    clear sorted_psb_idx
end

%%
% Transform site numbers to site names.
load('data.mat','site_name');

% Results from p2_chi2_pathogenic_idx_3x_p_0_01.mat:
% 2938	3398	4526	3398	9196	5944	8380	2938	3341	4526	
%Results from p2_inf_norm_pathogenic_idx_3x_50.mat:
%2938	4526	3398	3398	9196	8380	2938	4526	250	5588	

% OVERALL RESULT:
result = [2938; 3398; 4526; 9196; 8380];

disp('Most likely pathogenic sites are:')
for i = 1 : 5
    result_site_name = site_name{result(i)};
    fprintf('%s\n',result_site_name)
end


%%











