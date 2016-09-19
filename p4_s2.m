
%
iterate = 0;
min_coef = 5;

%%
% Use all samples. Keep no samples for test.
x = genotype_3x;
y = multi_phenos; % [1000*10]
% Extract top 150 sites for each character(性状).
num_sites_each_character = 150;

% Use to save 150*10=1500 possible pathogenic sites(bits).
% There could be repetitions! Do we need to deal with it? Guess so...
possible_pathogenic_idx = zeros(num_sites_each_character,num_characters);

for i = 1 : num_characters
    
correlation = corr(x,y(:,i));
[sorted_correlation,sorted_correlation_idx] = sort(abs(correlation),'descend');

possible_pathogenic_idx(:,i) = sorted_correlation_idx(1:num_sites_each_character);
% possible_pathogenic_idx(:,i) = sort(possible_pathogenic_idx(:,i),'ascend');

% sorted_psb_idx = floor((sorted_correlation_idx - 1) / 3) + 1;


end

%%
% Remove repeated bits in possible_pathogenic_idx.
possible_pathogenic_idx = ...
    reshape(possible_pathogenic_idx,numel(possible_pathogenic_idx),1);
possible_pathogenic_idx = sort(possible_pathogenic_idx,'ascend');
possible_pathogenic_idx = fun_delete_duplicate(possible_pathogenic_idx);
% Now there are 869 elements in possible_pathogenic_idx.

%%
Y = multi_phenos_sum.^2; % (Y(1)+...+Y(10))^2.

X = [];
for i = 1 : length(possible_pathogenic_idx)
    X = [X, genotype_3x(:,possible_pathogenic_idx(i))];
end

%%
% B = regress(Y,X) returns the vector B of regression coefficients
% in the linear model Y = X * B.
% NOTICE: must add constant term to the left of X!
X = [ones(num_samples,1),X];    % add constant term to X.
alpha = 0.05; % def: 0.05
[B,BINT,R,RINT,STATS] = regress(Y,X,alpha);

Yhat = X * B;
disp('Regress fit finished.')
fprintf('Originally X has %d columns.\n',size(X,2))

% 残差与残差区间的杠杆图
rcoplot(R,RINT)
str_title = ['site number:' num2str(size(X,2)) ...
    ' alpha:' num2str(alpha)];
title(str_title)

%%
if iterate
%-------------------------------------------
last_B_dim = 0;

fprintf('Start iteration to reduce dimension, min_coef = %g:\n',min_coef)
while 1
    for j = length(B) : -1 : 2 % Ignore column 1(constant term).
        if abs(B(j)) < min_coef
            % Remove terms with coef less than min_coef.
            X(:,j) = [];
        end
    end
    fprintf('Now X has %d columns.\n',size(X,2)) % DEBUG
    [B,BINT,R,RINT,STATS] = regress(Y,X,alpha);
    Yhat = X * B;
    
    % break condition: the dimension of stops decreasing.
    if last_B_dim == length(B)
        break;
    end
    last_B_dim = length(B);
end
% 残差与残差区间的杠杆图
figure(2)
rcoplot(R,RINT)
str_title = ['After iteration. site number:' num2str(size(X,2)) ...
    ' alpha:' num2str(alpha)];
title(str_title)
disp('Iteration finished.')
clear last_B_dim
%-------------------------------------------
end

%%

%TODO: give final results.


