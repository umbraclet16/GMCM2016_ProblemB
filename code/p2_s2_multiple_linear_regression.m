
%%
% Configurations.
% Choose encoding mode, 1 for 3 bits; 0 for 0~2.
use_genotype_3x = 1;
% Choose multiple linear regression method.
% 1: regress; 2: stepwisefit; 3: robustfit(logistic).
% TODO: REMOVE 2!!!
reg_method = 1;
% Method used to extract possible pathogenic sites(bits).
% 1: chi-square test; 2: infinite norm.
p2_extract_method = 2;
% Save result to .mat file?
p2_save_result = 0;
% Iterate regression to reduce dimension of final result?
iterate = 1;
% Terms with oefficients smaller than 'min_coef' will be removed
% in the iteration to reduce dimension of final result.
% try: 0.05; 0.1; 0.15; 0.2.
% (0.2 is already too much. Dimension is reduced to 4,
%  elements in Yhat(Y_test) are almost same.)
min_coef = 0.05;
%------------------------------------------------------------
% Threshold(p) in chi-square test.
% Takes value in [0.01;0.001;0.0001]
threshold = 0.01;
%------------------------------------------------------------
% Amount of bits extracted using infinite norm.
% Takes value in [1000;300;200;100;50]
amount_of_bits_extracted = 300;
%------------------------------------------------------------
%%
% Load possible_pathogenic_idx from .mat file.
if p2_extract_method == 1   % chi-square test.
    str1 = 'chi2';
    switch threshold
        case 0.01
            str2 = 'p_0_01';   % extract 276 sites.
        case 0.001
            str2 = 'p_0_001';   % extract 24 sites.
        case 0.0001
            str2 = 'p_0_0001';   % extract 5 sites.
        case 0.05
            str2 = 'p_0_05';
    end
%------------------------------------------------------------
elseif p2_extract_method == 2   % infinite norm.
    str1 = 'inf_norm';
    str2 = num2str(amount_of_bits_extracted);
else
    disp('p2_extract_method is invalid!')
end

mat_name_str = ['p2_' str1 '_pathogenic_idx_3x_' str2 '.mat'];

% ALWAYS load data from file!!!
% if ~exist('possible_pathogenic_idx','var')
    load(mat_name_str)
    fprintf('Loading from %s...\n',mat_name_str)
% end

clear p2_extract_method mat_name_str

%%
switch reg_method
    case 1
        reg = 1; stepwise_reg = 0; logistic = 0;
    case 2
        reg = 0; stepwise_reg = 1; logistic = 0;
    case 3
        reg = 0; stepwise_reg = 0; logistic = 1;
end

%%
% Multiple linear regression analysis on the 
% possible pathogenic sites(bits) obtained above.
%------------------------------------------------------------
% Solve for the vector B of regression coefficients
% in the linear model Y = X * B.
% Y: phenotype;
% X: most distinct n columns in genotype_3x corresponding to each of the 
%       possible pathogenic sites obtained above.
Y = phenotype;
% Remove 200 testing samples. Cannot change order of next 2 lines
% because it causes 'Matrix index is out of range for deletion.'
Y(901:1000,:) = [];
Y(401:500,:) = [];
X = [];

if use_genotype_3x
%------------------------------------------------------------
for i = 1 : length(possible_pathogenic_idx)
    X = [X, genotype_3x(:,possible_pathogenic_idx(i))];
end
%------------------------------------------------------------
else % use genotype(0~2)
%------------------------------------------------------------
for i = 1 : length(possible_pathogenic_idx)
    X = [X, genotype(:,possible_pathogenic_idx(i))];
end
%------------------------------------------------------------
end
% Remove 200 testing samples.
X(901:1000,:) = [];
X(401:500,:) = [];

%%
%=============================================
% Prepare test samples 'X_test' for calculating model accuracy
% later in 'p2_s3_calc_correct_pct.m'.
% This part is moved from 'p2_s3_calc_correct_pct.m'.
% Very UGLY but it works...
% TODO: figure out how to deal with this problem in p2_s3.
X_test = [genotype_3x(401:500,possible_pathogenic_idx); ...
          genotype_3x(901:1000,possible_pathogenic_idx)];

% Add constant term to X_test if using regress() or robustfit()!!!
if reg || logistic
    X_test = [ones(200,1),X_test];
end
%=============================================

%%
if reg
%------------------------------------------------------------
% B = regress(Y,X) returns the vector B of regression coefficients
% in the linear model Y = X * B.
% NOTICE: must add constant term to the left of X!
X = [ones(num_samples - 200,1),X];    % add constant term to X.
alpha = 0.05; % def: 0.05
[B,BINT,R,RINT,STATS] = regress(Y,X,alpha);

% STATS

% STATS contains the R-square statistic, the F statistic and p value
% for the full model, and an estimate of the error variance.
if STATS(1) >= 0.95     % R-square > 0.95
    fprintf('Multiple linear fitting succeeded!\n');
end

Yhat = X * B;
disp('Regress fit finished.')
fprintf('Originally X has %d columns.\n',size(X,2))

% 残差与残差区间的杠杆图
rcoplot(R,RINT)
str_title = ['site number:' num2str(size(X,2)) ...
    ' alpha:' num2str(alpha)];
title(str_title)
%------------------------------------------------------------
if iterate
%-------------------------------------------
% Already added constant term above!
% X = [ones(num_samples - 200,1) X];

% When threshold = 0.01:
% (筛除系数小于0.05的项，可降维至111.)
% Iterate to remove terms with coefficients < 0.05, B dimension -> 111;
% Iterate to remove terms with coefficients <  0.1, B dimension ->  31;
% Iterate to remove terms with coefficients < 0.15, B dimension ->   8;
% Iterate to remove terms with coefficients <  0.2, B dimension ->   4.
last_B_dim = 0;
% Used to record pathogenic sites(bits) indices during the iteration,
% due to that the dimension of X is dynamically decreased.
% -1 is for occupying the correspondent col of constant term in X.
psb_pathogenic_idx_iterated = [-1;possible_pathogenic_idx];

fprintf('Start iteration to reduce dimension, min_coef = %g:\n',min_coef)
while 1
    for j = length(B) : -1 : 2 % Ignore column 1(constant term).
        if abs(B(j)) < min_coef
            % Remove terms with coef less than min_coef.
            X(:,j) = [];
            % Remove corresponding columns in X_test.
            X_test(:,j) = [];
            % Remove corresponding index in record array.
            psb_pathogenic_idx_iterated(j) = [];
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
disp('Iteration finished.')

% Bits to sites. Also remove fake idx for constant term.
num_result = length(psb_pathogenic_idx_iterated);
result = floor((psb_pathogenic_idx_iterated(2:num_result) - ones(num_result-1,1)) / 3) + 1;
% Remove duplicates.
result = fun_delete_duplicate(result);
num_result = length(result);

fprintf('After the iteration, we get %d most likely pathogenic sites:\n',num_result)
for i = 1 : num_result
    fprintf('NO.%d: %d, coef = %g\n',i,result(i),B(i+1))
end
fprintf('Coefficient of constant term is %g.\n',B(1))

% 残差与残差区间的杠杆图
figure(2)
rcoplot(R,RINT)
str_title = ['After iteration. site number:' num2str(size(X,2)) ...
    ' alpha:' num2str(alpha)];
title(str_title)
clear last_B_dim
%-------------------------------------------
end
%------------------------------------------------------------
end
%%
% stepwise regression
if stepwise_reg
%------------------------------------------------------------
PENTER = 0.05;
PREMOVE = 0.10;
[B,SE,PVAL,INMODEL,STATS,NEXTSTEP,HISTORY] = stepwisefit(X,Y,'penter',PENTER,'premove',PREMOVE);

Yhat = X * B;
disp('Stepwise fit finished.')
clear SE PVAL INMODEL NEXTSTEP HISTORY % Don't need them now.
%------------------------------------------------------------
end

%%
if logistic
%------------------------------------------------------------
[B,STATS] = robustfit(X,Y,'logistic');
Yhat = [ones(num_samples - 200,1) X] * B;
disp('Logistic fit finished.')
fprintf('Originally X has %d columns.\n',size(X,2))

if iterate
%-------------------------------------------
% Add constant term to the left of X.
X = [ones(num_samples - 200,1),X];
% When threshold = 0.01:
% (筛除系数小于0.05的项，可降维至109.)
% Iterate to remove terms with coefficients < 0.05, B dimension -> 109;
% Iterate to remove terms with coefficients <  0.1, B dimension ->  31;
% Iterate to remove terms with coefficients < 0.15, B dimension ->  13;
% Iterate to remove terms with coefficients <  0.2, B dimension ->   4.
last_B_dim = 0;
% Used to record pathogenic sites(bits) indices during the iteration,
% due to that the dimension of X is dynamically decreased.
% -1 is for occupying the correspondent col of constant term in X.
psb_pathogenic_idx_iterated = [-1;possible_pathogenic_idx];

fprintf('Start iteration to reduce dimension, min_coef = %g:\n',min_coef)
while 1
    for j = length(B) : -1 : 2 % Ignore column 1(constant term).
        if abs(B(j)) < min_coef
            % Remove terms with coef less than min_coef.
            X(:,j) = [];
            % Remove corresponding columns in X_test.
            X_test(:,j) = [];
            % Remove corresponding index in record array.
            psb_pathogenic_idx_iterated(j) = [];
        end
    end
    fprintf('Now X has %d columns.\n',size(X,2)) % DEBUG
%--------------------
%     [B,~,~,~,~] = regress(Y,X);
    X(:,1) = [];
    [B,STATS] = robustfit(X,Y,'logistic');
    X = [ones(num_samples - 200,1),X];
%--------------------
    Yhat = X * B;
    
    % break condition: the dimension of stops decreasing.
    if last_B_dim == length(B)
        break;
    end
    last_B_dim = length(B);
end
disp('Iteration finished.')

% Bits to sites. Also remove fake idx for constant term.
num_result = length(psb_pathogenic_idx_iterated);
result = floor((psb_pathogenic_idx_iterated(2:num_result) - ones(num_result-1,1)) / 3) + 1;
% Remove duplicates.
result = fun_delete_duplicate(result);
num_result = length(result);

fprintf('After the iteration, we get %d most likely pathogenic sites:\n',num_result)
for i = 1 : num_result
    fprintf('NO.%d: %d, coef = %g\n',i,result(i),B(i+1))
end
fprintf('Coefficient of constant term is %g.\n',B(1))

clear last_B_dim
%-------------------------------------------
end
%------------------------------------------------------------
end

%%
% Calculate accuracy of discrimination using 200 test samples.
p2_s3_calc_correct_pct;

%%
% Save to .mat file.
if p2_save_result
%------------------------------------------------------------
% str1:extract_method; str2: threshold
switch reg_method
    case 1
        str_method = '_regress_';
    case 2
        str_method  = '_stepwise_';
    case 3
        str_method = '_logistic_';
end

if iterate
    str_iterate = '_iterate';
else
    str_iterate = [];
end

filename = ['p2_result_mat/p2_result_' str1 str_method str2 str_iterate '.mat'];
save(filename,'threshold','B','healthy_test','healthy_test_correct_pct',...
    'ill_test','ill_test_correct_pct','healthy_training','healthy_training_correct_pct',...
    'ill_training','ill_training_correct_pct','STATS')

fprintf('Save to file %s.\n',filename)
disp('========================================')
%------------------------------------------------------------
end
clear p2_save_result

%%
clear str1 str2 filename str_method reg_method iterate


