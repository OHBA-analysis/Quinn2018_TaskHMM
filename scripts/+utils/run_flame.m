function [cope,varcope,tstats] = run_flame( data, vars, design_matrix, contrasts )
%%function [tstat, thresh] = run_glm( data, contrasts )
%
%

% Preallocate arrays
cope = zeros(size(contrasts,1),size(data,2));
varcope = zeros(size(contrasts,1),size(data,2));
tstats = zeros(size(contrasts,1),size(data,2));

% Main loop
for ii = 1:size(data,2)
                                           
    % fit GLM with varcopes as random effects variables in error term
    [b, beta, covgam]=flame1(data(:,ii),design_matrix,vars(:,ii));
    cope(:,ii) = contrasts*b;
    varcope(:,ii) = sqrt(contrasts*covgam*contrasts');
         
end


% Create gaussian variance smoothing function
gaussdist =  @(x,mu,sd) 1/(2*pi*sd)*exp(-(x-mu).^2/(2*sd^2));
variance_smoothing = gaussdist((1:size(data,2))',size(data,2)/2,25);
% Normalise so that area under Gaussian is 1
variance_smoothing = variance_smoothing ./ sum(variance_smoothing);

for ii = 1:size(cope,1)
    tstats(ii,:) = cope(ii,:) ./ fftconv(varcope(ii,:),variance_smoothing);
end