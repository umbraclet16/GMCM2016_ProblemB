
% Use Fisher's combination test to find genes most relevant to the disease.

%%
% Configurations.
% If calculate P? Pass to save time if we can load P from file.
calc_P = 0;
% Method: 1: Fisher's combination test; 2: Most-significant SNP.
method = 2;
% Both methods come from <全基因组关联分析对银屑病遗传学研究的启示_张学军>.

%%
% Load data: base_combination; genotype; genotype_3x; site_name.
load('data.mat');

num_samples = 1000;
num_sites = 9445;
num_genes = 300;
phenotype = [zeros(500,1);ones(500,1)];
x = [phenotype(1:400);phenotype(501:900)];
P = zeros(num_sites,1);

%%
% Calculate p for each site.
if calc_P
%------------------------------------------------------------
tic
for i = 1 : num_sites
    [~,~,p] = crosstab(x,[genotype(1:400,i);genotype(501:900,i)]);
    % return params 'table' and 'chi2' are useless so replaced by '~'.

    P(i) = p;
end
toc
% Takes around 200s.

save('p3_method1.mat','P');
%------------------------------------------------------------
else
    load('p3_method1.mat','P')
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method 1: Fisher's combination test.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if method == 1
%------------------------------------------------------------
% Calculate statistic(统计量)
% Stat_X = -2ln(product_{i=1}^{n}(p_i)) = sigma_{i=1}^{n}[-2ln(p_i)].
% It obeys the chi square distribution with degree of freedom 2n.
Stat_X = zeros(num_genes,1);
for i = 1 : num_genes
    for j = 1 : num_sites_in_gene(i)
    Stat_X(i) = Stat_X(i) - 2 * log(P(index_gene(i) + j - 1));
    end
end

deg_of_freedom = 2 * num_sites_in_gene; % degree of freedom = 2 * n.

save('p3_method1.mat','Stat_X','deg_of_freedom','-append')

%%
% Calculate P_gene.
P_gene = zeros(num_genes,1);
for i = 1 : num_genes
    P_gene(i) = 1 - chi2cdf(Stat_X(i),deg_of_freedom(i));
end

save('p3_method1.mat','P_gene','-append')

%%
% Find pathogenic genes(satisfying p < 0.05).
pathogenic_genes = find(P_gene < 0.05);
num_pathogenic_genes = length(pathogenic_genes);
% Sort all genes by P in ascending order.
[sorted_P_gene,sorted_P_gene_idx] = sort(P_gene,'ascend');
% Sort pathogenic genes by P in ascending order.
sorted_pathogenic_genes = sorted_P_gene_idx(1:length(pathogenic_genes));
sorted_P_pathogenic_genes = sorted_P_gene(1:length(pathogenic_genes));

save('p3_method1.mat','sorted_pathogenic_genes','sorted_P_pathogenic_genes','-append')

%%
% Display result.
disp('By means of Fisher''s combination test,')
fprintf('found %d genes possibly relevant to the diseases,\n',num_pathogenic_genes)
disp('satisfying p < 0.05. They are:')
for i = 1 : length(sorted_pathogenic_genes)
    fprintf('gene NO: %d, p = %g.\n',sorted_pathogenic_genes(i),sorted_P_pathogenic_genes(i))
end
fprintf('\n')
%------------------------------------------------------------
end         % method 1: Fisher's combination test.

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% method 2: Most-significant SNP.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Take smallest p of sites in each gene as the p of that gene.
if method == 2
%------------------------------------------------------------
P_gene = zeros(num_genes,1);
for i = 1 : num_genes
    P_gene(i) = P(index_gene(i));
    for j = 2 : num_sites_in_gene(i)
        if P_gene(i) > P(index_gene(i) + j - 1)
            P_gene(i) = P(index_gene(i) + j - 1);
        end
    end
end

save('p3_method2.mat','P_gene')

%%
% Find pathogenic genes(top 10 smallest p).
num_pathogenic_genes = 10;
% Sort all genes by P in ascending order.
[sorted_P_gene,sorted_P_gene_idx] = sort(P_gene,'ascend');
% Sort pathogenic genes by P in ascending order.
sorted_pathogenic_genes = sorted_P_gene_idx(1:num_pathogenic_genes);
sorted_P_pathogenic_genes = sorted_P_gene(1:num_pathogenic_genes);

save('p3_method2.mat','sorted_pathogenic_genes','sorted_P_pathogenic_genes','-append')

%%
% Display result.
disp('By means of Most-significant SNP,')
fprintf('Top %d genes possibly relevant to the diseases are:\n',num_pathogenic_genes)
for i = 1 : length(sorted_pathogenic_genes)
    fprintf('gene NO: %d, p = %g.\n',sorted_pathogenic_genes(i),sorted_P_pathogenic_genes(i))
end
fprintf('\n')

%------------------------------------------------------------
end     % method 2: Most-significant SNP.

%%


