
% loopup variation(变异), i.e. is there a third base in a column?

%%
for i = 1 : n
    variation_happened = fun_lookup_variation(genotype(:,i));
    if variation_happened == 1
        disp('variation happened in col i.');
        break;
    end
end
