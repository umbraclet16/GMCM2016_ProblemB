
%
iterate = 0;
min_coef = 5;
num_psb_sites = 286; % p = 0.001: 28;  p = 0.01: 286; p = 0.05: 1481.

%%
% Use only samples with all 10 characters(10 '1' in multi_phenos,num=294),
% and samples with no characters(10 '0' in multi_phenos,num=300).
% Then the dimension of y is reduced from [1000*10] to [(294+300)*1].
% Model in P2 can be used.
y10 = find(multi_phenos_sum == 10); % [294*1]
y0 = find(multi_phenos_sum == 0); % [300*1]
y = [ones(length(y10),1);zeros(length(y0),1)];
x = zeros(length(y10)+length(y0),num_sites*3);
for i = 1 : length(y10)
    x(i,:) = genotype_3x(y10(i),:);
end
for i = 1 : length(y0)
   x(length(y10)+i,:) =  genotype_3x(y0(i),:);
end

[correlation,pval] = corr(x,y);

[sorted_psb,sorted_psb_idx] = sort(pval,'ascend');

possible_pathogenic_idx = sorted_psb_idx(1:num_psb_sites);
possible_pathogenic_idx = sort(possible_pathogenic_idx,'ascend');
sorted_psb_idx = floor((possible_pathogenic_idx - 1) / 3) + 1;
% Remove duplicate sites.
sorted_psb_idx = fun_delete_duplicate(sorted_psb_idx);

save('p4_result','sorted_psb_idx')

%%%%%%%%%%%%%%%%%%%%%%THE END%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
%%
Y = multi_phenos_sum;
Y(find(multi_phenos_sum ~= 10 && multi_phenos_sum ~= 0)) = [];

X = [];
for i = 1 : length(possible_pathogenic_idx)
    X = [X, genotype_3x(:,possible_pathogenic_idx(i))];
end
X(find(multi_phenos_sum ~= 10 && multi_phenos_sum ~= 0,:)) = [];

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

end

