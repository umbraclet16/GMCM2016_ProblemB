
% Extract most possible pathogenic sites(bits) using
% infinite norm to extract possible pathogenic sites(bits).

%%
% Configuration.
% Amount of bits to be extracted.
amount_of_bits_extracted = 300;
% Save to .mat file?
save_to_mat = 1;

%%
%------------------------------------------------------------
% Sort all bits by the difference of sum of healthy and sum of ill
% samples in descending order.
sum_1 = zeros(num_sites*3,1);   % healthy.
sum_2 = zeros(num_sites*3,1);   % ill.

sum_1(:) = sum(genotype_3x(1:400,:));
sum_2(:) = sum(genotype_3x(501:900,:));
dif = abs(sum_1 - sum_2);

[sorted_dif_bit,sorted_dif_bit_idx] = sort(dif,'descend');

%%
% Bits to Sites.
inf_norm_sites = zeros(num_sites,1);
inf_norm_idx_in_site = zeros(num_sites,1);
for i = 1 : num_sites
    [inf_norm_sites(i),inf_norm_idx_in_site(i)] = max(dif(3*(i-1)+1:3*i));
end

[sorted_inf_norm_sites,sorted_inf_norm_site_idx] = sort(inf_norm_sites,'descend');

%%
% I use sorted bits rather than sites, because it's easier to handle.
% But it's actually NOT infinite norm...
possible_pathogenic_idx = sorted_dif_bit_idx(1 : amount_of_bits_extracted);

%%
if save_to_mat
%------------------------------------------------------------
str = num2str(amount_of_bits_extracted);
    
filename = ['p2_inf_norm_pathogenic_idx_3x_' str];
save(filename,'possible_pathogenic_idx', ...
    'sorted_dif_bit','sorted_dif_bit_idx', ...
    'sorted_inf_norm_sites','sorted_inf_norm_site_idx')
fprintf('Saved to file: %s.\n',filename)
%------------------------------------------------------------
end

%%
clear str filename save_to_mat amount_of_bits_extracted

