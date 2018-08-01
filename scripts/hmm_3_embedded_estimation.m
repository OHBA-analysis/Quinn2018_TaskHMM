%% Overview
%
% Here we load in our source parcellated MEG data from the preprocessing stage,
% perform some normalisation and compute an Time-Delay-Embedded HMM.
%

config = utils.get_studydetails;

%%

% Find preprocessed data files
datapath = fullfile( config.analysisdir,'spm_sss_processed' );
s = study(datapath,'Giles');

% Load in epoch info and initialise parcellation object
load( fullfile( datapath,'epochinfo.mat') );
p = parcellation('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm');

% We'll need to collect the data, T and epoch information
data = [];            % HMM ready dataset
T = [];               % Length of continuous good segments
R = [];               % Indices of single run within data
B = cell(s.n,1);      % Indices of bad samples per session
trl = cell(s.n,1);    % epoch info per segment
runlen = zeros(s.n,1);          % Length of run per good segment

% Preallocate array to store ERF
nsamples = diff(epochinfo{1}.trl(1,1:2)) + 1;
ntrials = 148;
erf = nan(nsamples,p.n_parcels,ntrials,s.n);

% main normalisation loop
for ind = 1:s.n

    fprintf('Processing %s\n',s.fnames{ind});


    %-------------------------------
    % continuous file
    D = s.read(ind);
    D_orig = D.montage('switch',0);

    runlen(ind) = size(D,2);

    %-------------------------------
	% get data and orthogonalise
	dat = D(:,:,1);

	dat = ROInets.remove_source_leakage(dat, 'symmetric');

	%-------------------------------
    % Get badsamples
    runlen(ind) = size(dat,2);
    bs = ~good_samples( D );

    % find single good samples - bug when we have consecutive bad segments
    xx = find(diff(diff(bs)) == 2)+1;
    if ~isempty(xx)
        bs(xx) = 1;
    end

    % store bad samples
    B{ind}=find(bs);

    % indices of good samples
    good_inds=setdiff(1:runlen(ind),B{ind});

    % remove bad samples,
    % replace with >> a = zeros(44,size(D,2))*nan;a(:,inds) = dat;

    dat = dat(:,good_inds);

    if any(bs)

        t_good = ~bs;
        db = find(diff([0; t_good(:); 0]));
        onset = db(1:2:end);
        offset = db(2:2:end);
        t = offset-onset;

        % sanity check
        if size(dat,2) ~= sum(t)
            disp('Mismatch between Data and T!!');
        end
    else
        t = size(dat,2);
    end

    %--------------------------------
    % Store info

    offset = sum(T);

    R = cat(1,R,[offset+1 offset+size(dat,2)]);

    T = cat(1,T,t);

    data = cat(2,data,dat);

    %--------------------------------
    % Check evoked result

    % get trial info
    trl{ind} = epochinfo{ind}.trl;

    % replace bad samples
    subj_data = nan(p.n_parcels,size(D,2));
    subj_data(:,good_inds) = data(:,R(ind,1):R(ind,2));

    for kk = 1:size(trl{ind},1)
        start = trl{ind}(kk,1);
        stop = trl{ind}(kk,2);

        if stop > size(subj_data,2) || start > size(subj_data,2)
            disp('skipping');
            continue
        else
           erf(:,:,kk,ind) = subj_data(:,start:stop)';
        end
    end
end

%%

% Define HMM folder and save HMM-ready data
hmm_folder = fullfile( config.analysisdir, 'embedded_hmm' );
if ~exist( hmm_folder )
    mkdir( hmm_folder );
end
outfile = fullfile( hmm_folder, 'embedded_hmm_data' );

save( outfile, 'data', 'R', 'T', 'B', 'runlen', '-v7.3' );

%% Check ERF
%
% Here we plot the ERF created above as a final check that the input data
% to the HMM is properly aligned with respect to our triggeres

time_vect = linspace(-1,2,751) - .032;
figure;
plot(time_vect, nanmean(nanmean(erf,4),3))
grid on;
xlabel('Time (secs)')
ylabel('Amplitude Envelope')


%%
run_flip = true;
if run_flip
    % need to compute erp above
    x = squeeze(nanmean(erf,3));
    x = reshape(permute(x,[2 1 3]),39,[]);
    T2 = repmat(751,114,1);

    options_sf = struct();
    options_sf.maxlag = 5;
    options_sf.noruns = 20;
    options_sf.nbatch = 3;
    options_sf.verbose = 1;
    flips = findflip(x',T2',options_sf);

    T3 = [R(1,2) sum(diff(R(:,2)))]; % Single scan sessions
    data = flipdata(data',T3',flips);

    % Define HMM folder and save HMM-ready data
    outfile = fullfile(hmm_folder, 'embedded_hmm_data_flipped.mat' );

    save( outfile, 'data', 'R', 'T', 'B', 'runlen', '-v7.3' );
end


%% HMM inference
%
% Here we infer the HMM itself, a detailed description of the HMM-MAR toolbox
% can be found on https://github.com/OHBA-analysis/HMM-MAR/wiki
%

load(outfile);

% Prepare options structure
options = struct();
options.verbose = 1;

% These options specify the data and preprocessing that hmmmar might perform. Further options are discussed here
options.onpower = 0;
options.standardise = 0;
options.Fs = 250;

% Here we specify the HMM parameters
options.K = 6;  	         % The number of states to infer
options.order = 0; 	         % The lag used, this is only relevant when using MAR observations
options.zeromean = 1; 	     % We do not want to model the mean, so zeromean is set on
options.covtype = 'full';    % We want to model the full covariance matrix
options.embeddedlags = -7:7; % 15 lags are used from -7 to 7
options.pca = 39*4;          % The PCA dimensionality reduction is 4 times the number of ROIs

% These options specify parameters relevant for the Stochastic inference. They
% may be omitted to run a standard inference, but this will greatly increase
% the memory and CPU demands during processing. A detailed description of the
% Stochastic options and their usage can be found here:
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#stochastic
options.BIGNinitbatch = 15;
options.BIGNbatch = 15;
options.BIGtol = 1e-7;
options.BIGcyc = 500;
options.BIGundertol_tostop = 5;
options.BIGdelay = 5;
options.BIGforgetrate = 0.7;
options.BIGbase_weights = 0.9;

% The following loop performs the main HMM inference. We start by
% estimating a 6 state HMM as used in the manuscript.
states_to_infer = [6];

% Optionally, we can explore a wider range of values for K by looping through
% several values. This can be done by uncommenting the line below.
% Warning: This is likely to be extremely time-consuming to infer!

%states_to_infer = 2:2:12; % uncomment this line to explore different numbers of states

% The HMM inference is repeated a number of times and the results based on
% the iteration with the lowest free energy. Note that this can be
% extremely time-consuming for large datasets. For a quick exploration of
% results, nrepeats can be set to a smaller value or even 1. The full inference
% is run over 10 repeats.
nrepeats = 1;

for kk = states_to_infer
    best_freeenergy = nan;
    options.K = kk;

    for irep = 1:nrepeats
        % Run the HMM, note we only store a subset of the outputs
        % more details can be found here: https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#estimation
        [hmm_iter, Gamma_iter, ~, vpath_iter, ~, ~, ~, ~, fehist] = hmmmar (data,T',options);

        if isnan(best_freeenergy) || fehist(end) < best_freeenergy
            hmm = hmm_iter;
            Gamma = Gamma_iter;
            vpath = vpath_iter;
        end
    end
    % Save the HMM outputs
    hmm_outfile = fullfile( config.analysisdir, 'embedded_hmm', sprintf('embedded_HMM_K%d',options.K));
    save( hmm_outfile ,'hmm','Gamma','vpath','T')
end


%% Load 6 state HMM
% The state-wise spectra will be computed for 6 states

nstates = 6;
hmm_outfile = fullfile( config.analysisdir, 'embedded_hmm', sprintf('embedded_HMM_K%d',nstates) );
load( hmm_outfile ,'hmm','Gamma','vpath','T')

%% Statewise-Spectra
% Next we estimate the state-wise multitaper for each subject and each parcel.

% account for delay embedding in state gammas
pad_options = struct;
pad_options.embeddedlags = -7:7;
Gamma = padGamma(Gamma, T, pad_options);

if size(Gamma,1) ~= size(data,1)
	warning('The size of data and Gamma do not match');
end

% These options specify the how the spectra will be computed. Full details can
% be found here:
% https://github.com/OHBA-analysis/HMM-MAR/wiki/User-Guide#spectra

spec_options = struct();
spec_options.fpass = [1 40];
spec_options.p = 0; % no confidence intervals
spec_options.to_do = [1 0]; % no pdc
spec_options.win = 256;
spec_options.embeddedlags = -7:7;
spec_options.Fs = D.Fsample;

psd = zeros(114,nstates,39,39,39);
coh = zeros(114,nstates,39,39,39);
for ind = 1:114
    disp(ind);
    subj_data = data(R(ind,1):R(ind,2),:);

    fit = hmmspectramt(subj_data,R(ind,2)-R(ind,1),Gamma(R(ind,1):R(ind,2),:),spec_options);
    for jj = 1:6
        psd(ind,jj,:,:,:) = fit.state(jj).psd;
        coh(ind,jj,:,:,:) = fit.state(jj).coh;
    end
    clear fit subj_data
end

% Save the MT outputs
mt_outfile = fullfile( config.analysisdir, 'embedded_hmm', sprintf('embedded_HMM_K%d_spectra',options.K));
save( mt_outfile ,'psd','coh')

