function genotype_num = fun_genotype_encoding( cell, base_combination )
% Convert genotype coding mode from cell consisting of 'ATCG'
% to numerical code.
% Input:
%       cell            : genotype_cell[1000*9445];
%       base_combination: the base combination array[9455*1];
%                         (AT = 1; AC = 2; AG = 3; TC = 4; TG = 5; CG = 6.)
% Output:
%       the numerical genotype array[1000*9445].
%       
% TODO: temporarily element value is 0~2. We'll look for better way.

m = 1000;
n = 9445;
genotype_num = zeros(m,n);

for j = 1 : n
    switch base_combination(j)
%%
        case 1              % AT
            for i = 1 : m
                if (cell{i,j}(1) == 'A' && cell{i,j}(2) == 'A')
                    genotype_num(i,j) = 0;
                else if (cell{i,j}(1) == 'T' && cell{i,j}(2) == 'T')
                        genotype_num(i,j) = 2;
                    else
                        genotype_num(i,j) = 1;
                    end
                end
            end
            
%%
        case 2              % AC
            for i = 1 : m
                if (cell{i,j}(1) == 'A' && cell{i,j}(2) == 'A')
                    genotype_num(i,j) = 0;
                else if (cell{i,j}(1) == 'C' && cell{i,j}(2) == 'C')
                        genotype_num(i,j) = 2;
                    else
                        genotype_num(i,j) = 1;
                    end
                end
            end
            
%%
        case 3              % AG
            for i = 1 : m
                if (cell{i,j}(1) == 'A' && cell{i,j}(2) == 'A')
                    genotype_num(i,j) = 0;
                else if (cell{i,j}(1) == 'G' && cell{i,j}(2) == 'G')
                        genotype_num(i,j) = 2;
                    else
                        genotype_num(i,j) = 1;
                    end
                end
            end
            
%%
        case 4              % TC
            for i = 1 : m
                if (cell{i,j}(1) == 'T' && cell{i,j}(2) == 'T')
                    genotype_num(i,j) = 0;
                else if (cell{i,j}(1) == 'C' && cell{i,j}(2) == 'C')
                        genotype_num(i,j) = 2;
                    else
                        genotype_num(i,j) = 1;
                    end
                end
            end
            
%%
        case 5              % TG
            for i = 1 : m
                if (cell{i,j}(1) == 'T' && cell{i,j}(2) == 'T')
                    genotype_num(i,j) = 0;
                else if (cell{i,j}(1) == 'G' && cell{i,j}(2) == 'G')
                        genotype_num(i,j) = 2;
                    else
                        genotype_num(i,j) = 1;
                    end
                end
            end
            
%%
        case 6              % CG
            for i = 1 : m
                if (cell{i,j}(1) == 'C' && cell{i,j}(2) == 'C')
                    genotype_num(i,j) = 0;
                else if (cell{i,j}(1) == 'G' && cell{i,j}(2) == 'G')
                        genotype_num(i,j) = 2;
                    else
                        genotype_num(i,j) = 1;
                    end
                end
            end
%%    
    end
end












