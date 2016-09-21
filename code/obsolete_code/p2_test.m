method = 2;

%%
% manual chi2 test.
if method == 1
%------------------------------------------------------------
chi2 = [];
chi2_index = [];
for i = 1 : num_sites*3
temp = fun_calc_chi2(genotype_3x(:,i));
chi2 = [chi2 temp];
chi2_index = [chi2_index i];
end

% [chi2_max,max_idx] = max(chi2);
[chi2_max2min,I] = sort(chi2,'descend');
%------------------------------------------------------------
end

%%
% infinite norm.
if method == 2
%------------------------------------------------------------
sum_1 = zeros(num_sites*3,1);   % healthy.
sum_2 = zeros(num_sites*3,1);   % ill.
% sign_of_sum = zeros(num_sites*3,1);  % 前５００中０多则取０；１多则取１

sum_1(:) = sum(genotype_3x(1:500,:));
sum_2(:) = sum(genotype_3x(501:1000,:));
dif = abs(sum_1 - sum_2);
% sign_of_sum(sum1 - sum2 > 0) = 1;

[sort_dif,dif_idx] = sort(dif,'descend');

inf_norm = zeros(num_sites,1);
inf_norm_idx_in_site = zeros(num_sites,1);
for i = 1 : num_sites
    [inf_norm(i),inf_norm_idx_in_site(i)] = max(dif(3*(i-1)+1:3*i));
%     inf_norm_idx(i) = inf_norm_idx_in_site(i) + 3 * (i-1);
end

[sort_inf_norm,inf_norm_idx] = sort(inf_norm,'descend');

%------------------------------------------------------------
end

%%
% interval constraint.
if method == 3
%------------------------------------------------------------
possible = zeros(9445,1);
cnt0 = 0; cnt1 = 0; cnt2 = 0;

for j = 1 : 9445
    for i = 1 : 1000
        switch genotype(i,j)
            case 0
                cnt0 = cnt0 + 1;
            case 1
                cnt1 = cnt1 + 1;
            case 2
                cnt2 = cnt2 + 1;
        end
    end
    cnt = max(cnt0,cnt1);
    cnt = max(cnt,cnt2);
    cnt0 = 0; cnt1 = 0; cnt2 = 0;
    if cnt > 450 && cnt < 550       % still over 5000 possible...
        possible(j) = 1;
    end
end

num_possible = 0;
for i = 1 : 9445
    if possible(i) == 1
        num_possible = num_possible + 1;
    end
end
%------------------------------------------------------------
end

%%


