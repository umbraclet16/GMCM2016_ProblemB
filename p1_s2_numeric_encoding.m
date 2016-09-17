% Convert encoding mode from 'AGCT' to numerical code(0~2).

%%
% Configurations.
% Save encoding results(base_combination,genotype,genotype_3x) 
% to .mat files?
save_data_to_mat = 1;
% Clear genotype_cell to save memory?
save_memory = 0;

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

if save_data_to_mat
    save('data.mat','base_combination','-append');
    % '-append' is neccessary, or former data in .mat file will be lost!!!
end
    
%%
% transfer to numerical coding and save to file.
% Takes around 7s.
tic
% genotype[1000*9445]
genotype = fun_genotype_encoding(genotype_cell,base_combination);
if save_data_to_mat
    save('data.mat','genotype','-append');
end
toc

%%
% Remove genotype_cell to save memory.
%------------------------------------------------------------
if save_memory
    clear genotype_cell
end
clear save_memory
%------------------------------------------------------------

%%
% Transform encoding from 1 digit[0/1/2] to 3 bit[0/1].
% Takes less than 1s.
encoding_3x = 1;
%------------------------------------------------------------
if encoding_3x
    tic
    num_samples = 1000;
    n_3x = n * 3;
    genotype_3x = zeros(num_samples,n_3x);
    for i = 1 : num_samples
        for j = 1 : n
            if genotype(i,j) == 0
                genotype_3x(i,3*(j-1)+1) = 1;
            else if genotype(i,j) == 1
                    genotype_3x(i,3*(j-1)+2) = 1;
                else if genotype(i,j) == 2
                        genotype_3x(i,3*j) = 1;
                    else        % should never happen.
                        fprintf('Error: genotype(%d,%d) is beyond 0 ~ 2.\n',i,j)
                        break;
                    end
                end
            end
        end
    end
    
    if save_data_to_mat
        save('data.mat','genotype_3x','-append');
    end
    toc
end
%------------------------------------------------------------
clear encoding_3x

%%
clear i j




