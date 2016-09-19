%

%%
healthy = genotype_3x(1:500,possible_pathogenic_idx);
ill = genotype_3x(501:1000,possible_pathogenic_idx);
% ill = genotype_3x(501:1000,:);

X = ill;

threshold = 0.05;
similar = [];
for i = 5 : 500
    [~,~,p] = crosstab(X(4,:),X(i,:));
    % return params 'table' and 'chi2' are useless so replaced by '~'.
    if p > threshold
        similar = [similar i];
        % display new index immediately so we know how far we've gone...
        fprintf('similar_idx = %d.\n',i);
    end
    fprintf('.')
end










%%
if 0
% cidx = kmeans(X, 5, 'distance', 'sqeuclid');
cidx = kmeans(X, 20, 'distance', 'correlation');

% silhouette(X, cidx, 'sqeuclid');
silhouette(X, cidx, 'correlation');
% S = silhouette(X, CLUST);

end

%%
if 0
X2 = zscore(X);
Y2 = pdist(X2);

Z2 = linkage(Y2);

C2 = cophenet(Z2,Y2);

T = cluster(Z2,6);
H = dendrogram(Z2);

end

%%
if 0
D = pdist(X);
Z = squareform(D);

end
%%