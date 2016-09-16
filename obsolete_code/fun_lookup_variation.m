function variation = fun_lookup_variation( cell_col )
%%
% loopup variation(变异), i.e. is there a third base in a column?
% return 1 if variation happens, otherwise 0.
%%
variation = 0;
ATCG = 0;
for i = 1 : 1000
    if cell_col{i}(1) == 'A' || cell_col{i}(2) == 'A'
        bitset(ATCG,4);         % set bit 4.
%         ATCG = bitor(ATCG,8);
    end

    if cell_col{i}(1) == 'T' || cell_col{i}(2) == 'T'
        bitset(ATCG,3);         % set bit 3.
    end

    if cell_col{i}(1) == 'C' || cell_col{i}(2) == 'C'
        bitset(ATCG,2);         % set bit 2.
    end
    
    if cell_col{i}(1) == 'G' || cell_col{i}(2) == 'G'
        bitset(ATCG,1);         % set bit 1.
    end
    
    i = i + 1;
end

if bitget(ATCG,1)+bitget(ATCG,2)+bitget(ATCG,3)+bitget(ATCG,4) > 2
    variation = 1;
end

end

