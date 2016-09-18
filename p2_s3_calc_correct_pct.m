% Calculate accuracy of genetic disease discrimination using 
% 200 test samples(100 healthy and 100 ill).

%%
% Calculate model accuracy using training samples:
healthy_training = size(find(Yhat(1:400) < 0.5),1);
healthy_training_correct_pct = healthy_training / 400;

% NOTICE: test samples are removed, so indices of ill samples
% are 401:800!!!
ill_training = size(find(Yhat(401:800) > 0.5),1);
ill_training_correct_pct = ill_training / 400;

fprintf('healthy_training = %d, ill_training = %d.\n',healthy_training,ill_training)
fprintf('healthy_training_correct_pct = %g, ill_training_correct_pct = %g.\n',...
    healthy_training_correct_pct,ill_training_correct_pct)

%%
% Calculate model accuracy using test samples:
%=============================================
% 'X_test' has to be updated together with 'X' in the iterate procedure
% in p2_s2_multiple_linear_regression.m, because deleting cols causes
% the indices change constantly which I don't know how to handle.
% Therefore I move this part of work to p2_s2, so jump it here.
if 0
%-------------------------------------------
X_test = [genotype_3x(401:500,possible_pathogenic_idx); ...
          genotype_3x(901:1000,possible_pathogenic_idx)];

% Add constant term to X_test if using regress() or robustfit()!!!
if reg || logistic
    X_test = [ones(200,1),X_test];
end

% Delete corresponding columns that were deleted in the 
% training samples(X) in the iterate procedure.
for i = length(record_deleted_col) : -1 : 1 % Cannot omit '-1' in descending order.
        X_test(:,record_deleted_col(i)) = [];
end
%-------------------------------------------
end
%=============================================
Y_test = X_test * B;

healthy_test = size(find(Y_test(1:100) < 0.5),1);
healthy_test_correct_pct = healthy_test / 100;

ill_test = size(find(Y_test(101:200) > 0.5),1);
ill_test_correct_pct = ill_test / 100;

fprintf('healthy_test = %d, ill_test = %d.\n',healthy_test,ill_test)
fprintf('healthy_test_correct_pct = %g, ill_test_correct_pct = %g.\n',...
    healthy_training_correct_pct,ill_training_correct_pct)

%%