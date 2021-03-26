function RMSE = calc_RMS_resid(error)
%%% Calculate the RMS error for a set of data
%%% Assumes true error = 0.
n = length(error(:,1));
RMSE = zeros(n,1);
for k = 1:n
    RMSE_comp = sqrt(nanmean(error(k,:)).^2);
    RMSE(k) = RMSE_comp;
end
end