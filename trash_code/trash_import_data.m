
%%
clc
clear
close all

%%
% read file
[fid, msg] = fopen('genotype.dat'); % 1001 * 9445
if fid == -1
    disp(msg);  % display error info
end
fclose(fid);

% index = 1;
% A = fread(fid,1001*9445);     % ASCII!
A= textread('genotype.dat','%s',1001*9445);
