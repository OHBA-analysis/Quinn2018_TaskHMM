function [nnmf_res,ss] = run_nnmf( S, niterations, summary_plot )
% function [nnmf_res,ss] = wh.run_nnmf( S, niters, summary_plot )
%
%

if nargin < 3 || isempty( summary_plot )
    summary_plot = false;
end

if nargin < 2 || isempty( niterations )
    niterations = 10;
end

if summary_plot == true && niterations > 10
    warning('Summary plot is likely to be crowded with more than 10 iterations!');
end

S.do_plots = 0;

% Preallocate for SumSquare of residuls
ncomps = S.maxP;
nsamples = size( S.psds,3 );
ss = zeros( niterations, ncomps);

% Specify fit function, a unimodal gaussian
gauss_func = @(x,f) f.a1.*exp(-((x-f.b1)/f.c1).^2);

% Default fit options
options = fitoptions('gauss1');

% constrain lower and upper bounds
options.Lower = [0,1,0];
options.Upper = [Inf,nsamples,nsamples];

% Main loop
winning_value = Inf;
if summary_plot == true
    specs = zeros( ncomps, nsamples, niterations);
end

for ii = 1:niterations

    next_nnmf = teh_spectral_nnmf( S );

    for jj = 1:ncomps
        f = fit( linspace(1,nsamples,nsamples)',next_nnmf.nnmf_coh_specs(jj,:)', 'gauss1',options);
        resid = next_nnmf.nnmf_coh_specs(jj,:) - gauss_func(1:nsamples,f);
        ss(ii,jj) = sum( resid.^2 );
    end

    if sum(ss(ii,:)) < winning_value
        nnmf_res = next_nnmf;
        winning_value = sum(ss(ii,:));
    end

    if summary_plot == true
        specs(:,:,ii) = next_nnmf.nnmf_coh_specs;
    end

end


if summary_plot
    nrows = ceil( niterations/5 );
    winning_ind = find(winning_value == sum(ss,2));

    figure('Position',[100 100 1536 768])
    for ii = 1:niterations
        subplot( nrows,5, ii);
        plot( specs(:,:,ii)','linewidth',2);grid on;
        title_text = [ num2str(ii) '- SS: ' num2str(sum(ss(ii,:)))];
        if ii == winning_ind
            title_text = [title_text ' - WINNER'];
        end
        title(title_text,'FontSize',14);
    end

    figure
    x_vect = 1:niterations;
    h = bar(x_vect, sum(ss,2) );
    grid on;hold on
    bar(x_vect(winning_ind),sum(ss(winning_ind,:)),'r')
    xlabel('Iteration')
    ylabel('Residual Sum Squares')
    set(gca,'FontSize',14);

end

