% Read data from file 'genotype.dat',
% and save to cell 'genotype_cell'(基因型数据) and 'site_name'(点位名称).

%%
clc
clear
close all

%%
tic
% read genotype(基因型) from file genotype.dat[1001*9445]. Takes around 4.5s.
m = 1001;
n = 9445;
data_raw = textread('genotype.dat','%s',m * n);    % col vector.
% reshape. Takes another 3s.
genotype_cell = reshape(data_raw,n,m);
genotype_cell = genotype_cell';
clear data_raw;
toc

% mat file will be very large, causing MATLAB to crash!!!...
% save to mat file to save time.
% save('genotype.mat', 'genotype');

%%
% separate site(位点) name and allelic gene(等位基因) info.
site_name = cell(1,n);
for i = 1 : n
    site_name{i} = genotype_cell{1,i};
end

genotype_cell(1,:) = []; % remove first line.

%%
clear i












