
%

%%
training_samples = genotype_3x;
% training_samples(701:1000,:) = [];
% training_samples(201:500,:) = [];
training_samples(901:1000,:) = [];
training_samples(401:500,:) = [];
X  = training_samples;
Y = phenotype;
% Y(701:1000,:) = [];
% Y(201:500,:) = [];
Y(901:1000,:) = [];
Y(401:500,:) = [];

correlation = corr(X,Y);
[sorted_correlation,sorted_correlation_idx] = sort(abs(correlation),'descend');

possible_pathogenic_idx = sorted_correlation_idx(1:200);
possible_pathogenic_idx = sort(possible_pathogenic_idx,'ascend');

