function combination = fun_get_base_combination( cell_col )
%%
% which 2 bases are in the current column?
% All 6 possibilities:
% AT = 1; AC = 2; AG = 3; TC = 4; TG = 5; CG = 6.
debug = 0;
%%
ATCG = 0;
i = 1;
while bitget(ATCG,1)+bitget(ATCG,2)+bitget(ATCG,3)+bitget(ATCG,4) < 2
    
    if debug
    fprintf('size of cell_col(%d) is %d.\n',i,size(cell_col,1));
    end
    
    if cell_col{i}(1) == 'A' || cell_col{i}(2) == 'A'
        ATCG = bitor(ATCG,8);         % set bit 4.
%         bitset(ATCG,4);   % This does not change ATCG! Only return value.
        
    end

    if cell_col{i}(1) == 'T' || cell_col{i}(2) == 'T'
        ATCG = bitor(ATCG,4);         % set bit 3.
%         bitset(ATCG,3);
    end

    if cell_col{i}(1) == 'C' || cell_col{i}(2) == 'C'
        ATCG = bitor(ATCG,2);         % set bit 3.
%         bitset(ATCG,2);
    end
    
    if cell_col{i}(1) == 'G' || cell_col{i}(2) == 'G'
        ATCG = bitor(ATCG,1);         % set bit 3.
        bitset(ATCG,1);
    end
    
    i = i + 1;
end

switch ATCG
    case 12        % 8 + 4 = A + T
        combination = 1;
    case 10        % 8 + 2 = A + C
        combination = 2;
    case 9         % 8 + 1 = A + G
        combination = 3;
    case 6         % 4 + 2 = T + C
        combination = 4;
    case 5         % 4 + 1 = T + G
        combination = 5;
    case 3         % 2 + 1 = C + G
        combination = 6;
    otherwise
        combination = 0;
end

end

