% Convert encoding mode from 'AGCT' to numerical code(0~2).

%%
% get base combination(碱基组合) of each column.
% AT = 1; AC = 2; AG = 3; TC = 4; TG = 5; CG = 6.
% Takes less than 1s.
tic
base_combination = zeros(n,1);
for i = 1 : n
    base_combination(i) = fun_get_base_combination(genotype_cell(:,i));
end
toc

save('base_combination.mat','base_combination');

%%
% transfer to numerical coding and save to file.
% Takes around 6s.
tic
% genotype[1000*9445]
genotype = fun_genotype_encoding(genotype_cell,base_combination);
toc

%%








