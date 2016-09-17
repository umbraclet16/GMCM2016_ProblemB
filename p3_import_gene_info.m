
%%
% This is convenient, but filenames are sorted like:
% gene_1.dat, gene_10.dat, gene_100.dat, gene_101.dat, gene_102.dat...
read_in_a_batch = 0;    % DO NOT USE
%------------------------------------------------------------
if read_in_a_batch
    filename = dir('gene_info/*.dat');
    filename(1).name;
end
%------------------------------------------------------------
clear read_in_a_batch

%%
% Read data from 300 files in an inconvinient way.
% Consumed time can be ignored.
num_genes = 300;    % 300 genes
line_max = 60;          % Each .dat file has no more than 60 lines.
% There are totally 9445 lines in all 300 files.

% Preallocate space for gene_info[60*300]
% which consists of site names in each gene.(每个基因包含的位点名称)
% Each COLUMN corresponds to one gene.
% TODO: this wastes spaces. Should I use 300 individual cells???
gene_info = cell(line_max,num_genes);

for j = 1 : num_genes
    filename = ['gene_info/gene_' num2str(j) '.dat'];
    temp = textread(filename,'%s');
    
    n_site = size(temp,1);  % number of sites in the current gene.
    for i = 1 : n_site
        gene_info{i,j} = temp{i};
    end
end

%%
% check if gene_info involves all 9445 sites.
check_gene_info = 0;        % Already CHECKED!
%------------------------------------------------------------
if check_gene_info
    cnt = 0;
    for i = 1 : line_max
        for j = 1 : num_genes
            if ~isempty(gene_info{i,j})
                cnt = cnt + 1;
            end
        end
    end
end
%------------------------------------------------------------
clear check_gene_info

%%
% check if the order of sites in gene_info accords with
% the order in site_name(i.e. genotype).
check_order = 0;        % Already CHECKED!
%------------------------------------------------------------
if check_order
    orders_accord_cnt = 0;
    index_genotype = 0;
    for j = 1 : num_genes
        for i = 1 : line_max
            if ~isempty(gene_info{i,j})
                index_genotype = index_genotype + 1;
                % if string in gene_info and site_name are same, cnt++
                if strcmp(gene_info{i,j},site_name{index_genotype})
%                     fprintf('index = %d.\n',index_genotype);
                    orders_accord_cnt = orders_accord_cnt + 1;
                end
            end
        end
    end
    
% Or take gene_info as a 1-D cell:

%     for i = 1 : line_max * num_shape
%         if ~isempty(gene_info{i})
%             index_genotype = index_genotype + 1;
%             
%             if strcmp(gene_info{i},site_name{index_genotype})
% %                 fprintf('index = %d.\n',index_genotype);
%                 orders_accord_cnt = orders_accord_cnt + 1;
%             end
%         end
%     end
    
    if orders_accord_cnt == n
        disp('The order of sites in gene_info accords with genotype.')
    end

end
%------------------------------------------------------------
clear check_order

%%
clear i j temp










