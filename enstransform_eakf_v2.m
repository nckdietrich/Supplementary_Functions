function post_ensemble = enstransform_eakf_v2(prior_ensemble, prior_obs_ensemble, y, R)
    % Lewis Redner: 6/4/20

    % PURPOSE: perform an ensemble adjustment kalman filter transform of
    % incoming tiegcm ensemble data
    
    % INPUTS:
        % prior_ensemble - mass density data in an nxm matrix, where n is
        % the densities from a 3d grid arranged into a single column vector
        % and m is the number of ensembles used
        
        % prior_obs_ensemble - observed electron density in same config as
        % mass density
        
        % y - observations of electron density by cosmic1, 2x no obs
        
        % R - error covariance value of each observation, fed in as a
        % vector

    %% Find key dimensions
    % find number of states in prior ensembles
    n = size(prior_ensemble,1);
    % find number of ensembles in prior ensembles
    m = size(prior_ensemble,2);
    % find number of entries in each observation
    yi = size(y,1);
    % find number of observations total
    no_obs = size(y,2);
    % initialise array for post ensemble
    post_ensemble = zeros(n,m);
    
    %% Set up Data Matrices and R
    % find the mean of the prior ensemble along each row, i.e. x bar f
    prior_ens_mean = mean(prior_ensemble,2);
    % calculate the state data matrix
    Dx = 1/sqrt(m-1) * (prior_ensemble-prior_ens_mean);
    % calculate the mean of the prior observations along each row
    prior_obs_mean = mean(prior_obs_ensemble,2);
    % calculate the observation data matrix
    Dy = 1/sqrt(m-1)*(prior_obs_ensemble-prior_obs_mean);
    
    % R is fine as is because we just want diagonal elements. No point
    % setting up a matrix just to pick off the diagonals we created it with
    
    %% Performing the transform
    % NOTE: yf is prior obs and y is actual obs. Double check with Nick
    
    % initialise the prior_ensemble which is changed in each iteration
    for i = 1:1
        % extract the observation data matrix values across all ensembles for the ith
        % observation
        Dyi = Dy(i,:);
        % extract the R value corresponding to the ith observation
        Ri = R(i);
        % update the mean of prior_obs_ensemble with the y observation
        post_obs_mean_i = prior_obs_mean(i) + (Dyi*Dyi')/(Dyi*Dyi'+Ri)*(y(1,i)-prior_obs_mean(i));
        % update the prior ensemble in observation space
        post_obs_i = sqrt(Ri/(Ri+Dyi*Dyi'))*(prior_obs_ensemble(i,:) - prior_obs_mean(i))+post_obs_mean_i;
        % update the prior ensemble member with the post observations in
        % state space
        post_ensemble = prior_ensemble + (Dx*Dyi')/(Dyi*Dyi')*(post_obs_i - prior_obs_ensemble(i,:));
        % update the incoming prior ensemble with the post ensemble for the
        % next iteration
        prior_ensemble = post_ensemble;
        % update the state data matrix Dx
        prior_ens_mean = mean(prior_ensemble,2);
        Dx = 1/sqrt(m-1) * (prior_ensemble-prior_ens_mean);
    end
end
