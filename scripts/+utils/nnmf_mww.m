function [wbest,hbest,normbest] = nnmf(a,k,varargin)
%NNMF Non-negative matrix factorization.
%   [W,H] = NNMF(A,K) factors the N-by-M matrix A into non-negative factors
%   W (N-by-K) and H (K-by-M).  The result is not an exact factorization,
%   but W*H is a lower-rank approximation to the original matrix A.  The W
%   and H matrices are chosen to minimize the objective function that is
%   defined as the root mean squared residual between A and the
%   approximation W*H.  This is equivalent to
%
%          D = norm(A-W*H,'fro')/sqrt(N*M)
%
%   The factorization uses an iterative method starting with random initial
%   values for W and H.  Because the objective function often has local
%   minima, repeated factorizations may yield different W and H values.
%   Sometimes the algorithm converges to solutions of lower rank than K,
%   and this is often an indication that the result is not optimal.
%
%   [W,H,D] = NNMF(...) also returns D, the root mean square residual.
%
%   [W,H,D] = NNMF(A,K,'PARAM1',val1,'PARAM2',val2,...) specifies one or
%   more of the following parameter name/value pairs:
% 
%      Parameter    Value 
%      'algorithm'  Either 'als' (default) to use an alternating least
%                   squares algorithm, or 'mult' to use a multiplicative
%                   update algorithm.
%      'w0'         An N-by-K matrix to be used as the initial value for W.
%      'h0'         A K-by-M matrix to be used as the initial value for H.
%      'replicates' The number of times to repeat the factorization, using
%                   new random starting values for W and H, except at the
%                   first replication if w0 and h0 are given (default 1).
%                   This tends to be most beneficial with the 'mult'
%                   algorithm. 
%        'options'  An options structure as created by the STATSET
%                   function.  NNMF uses the following fields:
%
%            'Display'   Level of display output.  Choices are 'off'
%                        (the default), 'final', and 'iter'.
%            'MaxIter'   Maximum number of steps allowed. The default
%                        is 100.  Unlike in optimization settings,
%                        reaching MaxIter is regarded as convergence.
%            'TolFun'    Positive number giving the termination tolerance
%                        for the criterion.  The default is 1e-4.
%            'TolX'      Positive number giving the convergence threshold
%                        for relative change in the elements of W and H.
%                        The default is 1e-4.
%            'UseParallel'
%            'UseSubstreams'
%            'Streams'   These fields specify whether to perform multiple
%                        replicates in parallel, and how to use random 
%                        numbers when generating the starting points for
%                        the replicates. For information on these fields 
%                        see PARALLELSTATS.
%                        NOTE: If 'UseParallel' is TRUE and 'UseSubstreams' is FALSE,
%                        then the length of 'Streams' must equal the number of workers 
%                        used by NNMF.  If a parallel pool is already open, this 
%                        will be the size of the parallel pool.  If a parallel pool 
%                        is not already open, then MATLAB may try to open a pool for 
%                        you (depending on your installation and preferences).
%                        To ensure more predictable results, it is best to use 
%                        the PARPOOL command and explicitly create a parallel pool 
%                        prior to invoking NNMF with 'UseParallel' set to TRUE.
%
%    Examples:
%        % Non-negative rank-2 approximation of the Fisher iris measurements
%        load fisheriris
%        [w,h] = nnmf(meas,2);
%        gscatter(w(:,1),w(:,2),species);
%        hold on; biplot(max(w(:))*h','VarLabels',{'sl' 'sw' 'pl' 'pw'},'positive',true); hold off;
%        axis([0 12 0 12]);
%
%        % Try a few iterations at several replicates using the
%        % multiplicative algorithm, then continue with more iterations
%        % from the best of these results using alternating least squares
%        x = rand(100,20)*rand(20,50);
%        opt = statset('maxiter',5,'display','final');
%        [w,h] = nnmf(x,5,'rep',10,'opt',opt,'alg','mult');
%        opt = statset('maxiter',1000,'display','final');
%        [w,h] = nnmf(x,5,'w0',w,'h0',h,'opt',opt,'alg','als');
%
%    See also BIPLOT, PCA, STATSET, PARALLELSTATS.

%   Copyright 2007-2013 The MathWorks, Inc.

% Reference:
%   M.W. Berry et al. (2007), "Algorithms and Applications for Approximate
%     Nonnegative Matrix Factorization," Computational Statistics and Data
%     Analysis, vol. 52, no. 1, pp. 155-173.

% The factorization is not uniquely defined.  This function normalizes H so
% that its rows have unit length. It orders the columns of W and rows of H
% so that the columns of W have decreasing length.

% Check required arguments
narginchk(2,Inf);
[n,m] = size(a);
if ~isscalar(k) || ~isnumeric(k) || k<1 || k>min(m,n) || k~=round(k)
    error(message('stats:nnmf:BadK'));
end

% Process optional arguments
pnames = {'algorithm' 'w0' 'h0' 'replicates' 'options'};
dflts  = {'als'       []   []   1            []       };
[alg,w0,h0,tries,options] = ...
        internal.stats.parseArgs(pnames,dflts,varargin{:});
    
% Check optional arguments
alg = internal.stats.getParamVal(alg,{'mult' 'als'},'ALGORITHM');
ismult = strncmp('mult',alg,numel(alg));
checkmatrices(a,w0,h0,k);
if ~isscalar(tries) || ~isnumeric(tries) || tries<1 || tries~=round(tries)
    error(message('stats:nnmf:BadReplicates'));
end

defaultopt = statset('nnmf');
tolx = statget(options,'TolX',defaultopt,'fast');
tolfun = statget(options,'TolFun',defaultopt,'fast');
maxiter = statget(options,'MaxIter',defaultopt,'fast');
dispopt = statget(options,'Display',defaultopt,'fast');

[~,dispnum] = internal.stats.getParamVal(dispopt, {'off','notify','final','iter'},'Display');
dispnum = dispnum - 1;

[useParallel, RNGscheme, poolsz] = ...
    internal.stats.parallel.processParallelAndStreamOptions(options,true);

usePool = useParallel && poolsz>0;

% Special case, if K is full rank we know the answer
if isempty(w0) && isempty(h0)
    if k==m
        w0 = a;
        h0 = eye(k);
    elseif k==n
        w0 = eye(k);
        h0 = a;
    end
end


% Define the function that will perform one iteration of the
% loop inside smartFor
loopbody = @loopBody;

% Suppress undesired warnings.
if usePool
    % On workers and client
    pctRunOnAll internal.stats.parallel.muteParallelStore('rankDeficientMatrix', ...
        warning('off','MATLAB:rankDeficientMatrix') );
else
    % On client
    ws = warning('off','MATLAB:rankDeficientMatrix');
end

% Prepare for in-progress 
if dispnum > 1 % 'iter' or 'final'
    if usePool
       % If we are running on a parallel pool, each worker will generate
       % a separate periodic report.  Before starting the loop, we
       % seed the parallel pool so that each worker will have an
       % identifying label (eg, index) for its report.
       internal.stats.parallel.distributeToPool( ...
           'workerID', num2cell(1:poolsz) );

       % Periodic reports behave differently in parallel than they do
       % in serial computation (which is the baseline).
       % We advise the user of the difference.
       warning(message('stats:nnmf:displayParallel2'));
       
       % Leave formatted by \t UI strings untranslated. 8/17/2011
       fprintf('    worker\t      rep\t   iteration\t     rms resid\t     |delta x|\n' );
    else
       if useParallel
          warning(message('stats:nnmf:displayParallel'));
       end
       fprintf('    rep\t   iteration\t   rms resid\t  |delta x|\n');
    end
end

try
    whbest = internal.stats.parallel.smartForReduce(...
        tries, loopbody, useParallel, RNGscheme, 'argmin');
catch ME
    % Revert warning setting for rankDeficientMatrix to value prior to nnmf.
    if usePool
        % On workers and on client
        pctRunOnAll warning(internal.stats.parallel.statParallelStore('rankDeficientMatrix').state,'MATLAB:rankDeficientMatrix');
    else
        % On client
        warning(ws);
    end
    rethrow(ME);
end

normbest = whbest{1};
wbest = whbest{3};
hbest = whbest{4};
% whbest{2} contains the iteration chosen for the best factorization,
% but it has no meaning except as a "reproducible" tie-breaker, and
% is not supplied as a return value.
    
if dispnum > 1   % 'final' or 'iter'
    fprintf('%s\n',getString(message('stats:nnmf:FinalRMSResidual',sprintf('%g',normbest))));
end

% Revert warning setting for rankDeficientMatrix to value prior to nnmf.
if usePool
    % On workers and on client
    pctRunOnAll warning(internal.stats.parallel.statParallelStore('rankDeficientMatrix').state,'MATLAB:rankDeficientMatrix');
else
    % On client
    warning(ws);
end

if normbest==Inf
    error(message('stats:nnmf:NoSolution'))
end

% MWW Put the outputs in a standard form - first normalize h
% note w is freq, h is space
if 1
    % original code - puts amplitude in w
    hlen = sqrt(sum(hbest.^2,2));
else
    % puts amplitude in h
    hlen = 1./sqrt(sum(wbest.^2,1))';
end

if any(hlen==0)
    warning(message('stats:nnmf:LowRank', k - sum( hlen==0 ), k));
    hlen(hlen==0) = 1;
end

wbest = bsxfun(@times,wbest,hlen');
hbest = bsxfun(@times,hbest,1./hlen);

% Then order by w
[~,idx] = sort(sum(wbest.^2,1),'descend');
wbest = wbest(:,idx);
hbest = hbest(idx,:);

% ---- Nested functions ----

    function cellout = loopBody(iter,S)
        if isempty(S)
            S = RandStream.getGlobalStream;
        end

        % whtry is a "temporary variable" and hence needs to be
        % reinitialized at start of each loop.  
        whtry = cell(4,1); % whtry{1} = norm of error
                           % whtry{3} = w
                           % whtry{4} = h

        % Get random starting values if required
        if( ~isempty(w0) && iter ==1 ) 
            whtry{3} = w0;
        else
            whtry{3} = rand(S,n,k); 
        end
        if( ~isempty(h0) && iter ==1 )
            whtry{4} = h0;
        else
            whtry{4} = rand(S,k,m); 
        end

        % Perform a factorization
        [whtry{3},whtry{4},whtry{1}] = ...
            nnmf1(a,whtry{3},whtry{4},ismult,maxiter,tolfun,tolx,...
                                  dispnum,iter,usePool);
        whtry{2} = iter;

        cellout = whtry;
end

end  % of nnmf

% -------------------
function [w,h,dnorm] = nnmf1(a,w0,h0,ismult,maxiter,tolfun,tolx,...
                                   dispnum,repnum,usePool)
% Single non-negative matrix factorization
nm = numel(a);
sqrteps = sqrt(eps);


% Display progress.  For parallel computing, the replicate number will be
% displayed under the worker performing the replicate.
if dispnum>1 % 'final' or 'iter' 
    if usePool 
        labindx = internal.stats.parallel.workerGetValue('workerID');
        dispfmt = '%8d\t%8d\t%8d\t%14g\t%14g\n';
    else
        dispfmt = '%7d\t%8d\t%12g\t%12g\n';
    end   
end    
    
for j=1:maxiter
    % note w is freq, h is space

    if ismult
        % Multiplicative update formula
        numer = w0'*a;
        h = max(0,h0 .* (numer ./ ((w0'*w0)*h0 + eps(numer))));
        numer = a*h';
        w = max(0,w0 .* (numer ./ (w0*(h*h') + eps(numer))));
    else
        % Alternating least squares
        h = max(0, w0\a);
        if maxiter>1
            w = max(0, a/h);
        else
            w=w0;
        end
    end
    
    % Get norm of difference and max change in factors
    d = a - w*h;
    dnorm = sqrt(sum(sum(d.^2))/nm);
    dw = max(max(abs(w-w0) / (sqrteps+max(max(abs(w0))))));
    dh = max(max(abs(h-h0) / (sqrteps+max(max(abs(h0))))));
    delta = max(dw,dh);
    
    % Check for convergence
    if j>1
        if delta <= tolx
            break;
        elseif dnorm0-dnorm <= tolfun*max(1,dnorm0)
            break;
        elseif j==maxiter
            break
        end
    end

    if dispnum>2 % 'iter'
       if usePool 
           fprintf(dispfmt,labindx,repnum,j,dnorm,delta);
       else
           fprintf(dispfmt,repnum,j,dnorm,delta);
       end
    end

    % Remember previous iteration results
    dnorm0 = dnorm;
    w0 = w;
    h0 = h;
end

if dispnum>1   % 'final' or 'iter'
    if usePool 
       fprintf(dispfmt,labindx,repnum,j,dnorm,delta);
   else
       fprintf(dispfmt,repnum,j,dnorm,delta);
   end
end

end

% ---------------------------
function checkmatrices(a,w,h,k)
% check for non-negative matrices of the proper size

if ~ismatrix(a) || ~isnumeric(a) || ~isreal(a) || any(any(~isfinite(a)))
    error(message('stats:nnmf:BadA'))
end
[n,m] = size(a);
if ~isempty(w)
    if ~ismatrix(w) || ~isnumeric(w)|| ~isreal(w) || any(any(w<0)) || any(any(~isfinite(w)))
        error(message('stats:nnmf:BadWNegativeValues'))
    elseif ~isequal(size(w),[n k])
        error(message('stats:nnmf:BadWSizeIsWrong', sprintf( '%d', n ), sprintf( '%d', k )));
    end
end
if ~isempty(h)
    if ~ismatrix(h) || ~isnumeric(h)|| ~isreal(h) || any(any(h<0)) || any(any(~isfinite(h)))
        error(message('stats:nnmf:BadHNegativeValues'))
    elseif ~isequal(size(h),[k m])
        error(message('stats:nnmf:BadHSizeIsWrong', sprintf( '%d', k ), sprintf( '%d', m )));
    end
end
end % checkmatrices

