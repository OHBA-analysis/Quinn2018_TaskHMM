function [cope,varcope,tstats] = run_glm( data, design_matrix, contrasts )
%%function [tstat, thresh] = run_glm( data, contrasts )
%
%

% Preallocate arrays
cope = zeros(size(contrasts,1),size(data,2));
varcope = zeros(size(contrasts,1),size(data,2));
tstats = zeros(size(contrasts,1),size(data,2));

% reject bad_trials containing nans
good_inds = ~isnan(sum(data,2));

% Precompute design matrix invs
pinvx = pinv(design_matrix(good_inds,:));
pinvxtx=pinv(design_matrix(good_inds,:)'*design_matrix(good_inds,:));

% Main loop
for ii = 1:size(data,2)
    
    % Compute per-time point GLM across the good trials
    [cope(:,ii),varcope(:,ii)] = glm_fast_for_meg(squeeze(data(good_inds,ii)),...
                                                  design_matrix(good_inds,:),...
                                                  pinvxtx,pinvx,contrasts');
end


% Create gaussian variance smoothing function
gaussdist =  @(x,mu,sd) 1/(2*pi*sd)*exp(-(x-mu).^2/(2*sd^2));
variance_smoothing = gaussdist((1:size(data,2))',size(data,2)/2,25);
% Normalise so that area under Gaussian is 1
variance_smoothing = variance_smoothing ./ sum(variance_smoothing);

for ii = 1:size(cope,1)
    tstats(ii,:) = cope(ii,:) ./ fftconv(varcope(ii,:),variance_smoothing);
end