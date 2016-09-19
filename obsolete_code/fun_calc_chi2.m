function chi2 = fun_calc_chi2( genotype_3x_col )

%%
a = 0; b = 0; c = 0; d = 0;
n = 1000;
for i = 1 : 500 % healthy
    if genotype_3x_col(i) == 1
        a = a + 1;
    end
end
c = 500 - a;

for i = 501 : 1000 % ill
    if genotype_3x_col(i) == 1
        b = b + 1;
    end
end
d = 500 - b;

% fprintf('a=%d,b=%d,c=%d,d=%d.\n',a,b,c,d)

chi2 = n*(a*d-b*c)*(a*d-b*c)/((a+b)*(c+d)*(a+c)*(b+d));

end

