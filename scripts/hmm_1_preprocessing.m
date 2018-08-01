%% Overview
%
% Our preprocessing will take place within a specified analysis directory and
% will process data using three sub-directories.
% raw_data - containing the fif files (after Maxfilter SSS has been applied)
% spm_sss - contains a the raw data imported into an SPM object with minimal processing (artefact channel labelling and coregistration)
% spm_sss_processed - contains the spm objects with the preprocessing pipeline applied.
%
% We begin with SSS fif files from downloaded from
% https://openfmri.org/dataset/ds000117/ these should be present in within the
% raw_data subfolder in the main analysis directory. The folder hierarchy from
% the openfmri.org download should be preserved.
%
% raw_data
% |--> sub001
%   |--> run_01_sss.fif
%   |--> run_02_sss.fif
%   |--> run_03_sss.fif
%   |--> run_04_sss.fif
%   |--> run_05_sss.fif
%   |--> run_06_sss.fif
%   |--> highres001.nii
% |--> sub002
%   |--> run_01_sss.fif
%   |--> run_02_sss.fif
%   |--> run_03_sss.fif
%   |--> run_04_sss.fif
%   |--> run_05_sss.fif
%   |--> run_06_sss.fif
%   |--> highres001.nii
% |
% |--> sub019
%   |--> run_01_sss.fif
%   |--> run_02_sss.fif
%   |--> run_03_sss.fif
%   |--> run_04_sss.fif
%   |--> run_05_sss.fif
%   |--> run_06_sss.fif
%   |--> highres001.nii
%
% These files will be converted into SPM12 format into the spm_sss directory
% before computing the coregistration. All SPM files are copied into the top
% level of spm_sss and named to reflect the relevant subject and session
% numbers.
%
% spm_sss
% |--> spm_sub01_run01.mat
% |--> spm_sub01_run01.dat
% |--> spm_sub01_run02.mat
% |--> spm_sub01_run02.dat
% |
% |--> spm_sub19_run06.mat
% |--> spm_sub19_run06.dat
%
% Finally the spm_sss_processed folder contains the files during and after MEEG
% data preprocessing. These files are copied from spm_sss with the same naming
% conventions. These files will be the inputs to the HMM data preperation
% scripts.
%
% The script will check for pre-existing files before running any
% preprocessing to avoid repeating stages, only missing data files will be
% processed.

%check raw data is in place
missing_runs = utils.check_raw_data;

%% IMPORT RAW FIFS WITH SSS
% Here, we read in the raw fif object and conver the data into SPM12 format and
% save a copy in our spm_sss directory.

indir = fullfile( config.datadir, 'raw_data' );
outdir = fullfile( config.analysisdir, 'spm_sss' );

if ~exist( outdir )
    mkdir(outdir);
end


subjects = 1:19;
sessions = 1:6;

% loop through sessions and subjects
for j = 1:length(subjects)
    for k = 1:length(sessions)
        if exist(getfullpath(fullfile(outdir,sprintf('spm_sub%02d_run_%02d.mat',subjects(j),sessions(k)))),'file')
            % We don't need to copy existing files
            fprintf('Already imported spm_sub%02d_run_%02d.mat\n',subjects(j),sessions(k));
        else
            % Run the conversion to SPM12 format
            fif_in = fullfile(indir,sprintf('sub%03d/run_%02d_sss.fif',subjects(j),sessions(k)));
            D = osl_import(fif_in,'outfile',getfullpath(fullfile(outdir,sprintf('spm_sub%02d_run_%02d',subjects(j),sessions(k)))));
        end
    end
end

%% Perform coregistration
% Two processing steps are carried out in spm_sss. Firstly, we relabel the
% relevant artefact (ECG/EOG) channels for later use and secondly, we run the
% coregistration. Neither of these steps interact with the MEEG data themselves
% so we keep them as a separate stage.
%

 % If true, |use_existing| will prevent rerunning the coregstration on any files which already contain a coregistration
use_existing = true;
subjects = 1:19;
sessions = 1:6;

% Main loop through subjects and sessions
for j = 1:length(subjects)
    for k = 1:length(sessions)

        % Load data in from spm_sss, note that any changes in this loop are saved into the same file.
        D = spm_eeg_load(getfullpath(fullfile(outdir,sprintf('spm_sub%02d_run_%02d',subjects(j),sessions(k)))));

        % Next we re-label the artefact channels
        D = D.chantype(find(strcmp(D.chanlabels,'EEG062')),'EOG');
        D = D.chanlabels(find(strcmp(D.chanlabels,'EEG062')),'VEOG');

        D = D.chantype(find(strcmp(D.chanlabels,'EEG061')),'EOG');
        D = D.chanlabels(find(strcmp(D.chanlabels,'EEG061')),'HEOG');

        D = D.chantype(find(strcmp(D.chanlabels,'EEG063')),'ECG');
        D = D.chanlabels(find(strcmp(D.chanlabels,'EEG063')),'ECG');

        D.save();

        % Skip coreg if needed
        if use_existing && isfield(D,'inv')
            fprintf('Already coregistered %s\n',D.fname);
            continue
        end

        % Coregistration is carried out using a call to osl_headmodel. This
        % function takes a struct detailing how the coregistration should be
        % carried out. Importantly, we specify a D object from spm_sss and a
        % strutrual MRI scan from raw_data.
        % Note that useheadshape is set to false, typically we would run this
        % stage with useheadshape set to true. As the MRI scans included in the
        % download have been defaced to ensure participant anonymity the
        % headshape points on the face and nose can cause rotational errors in
        % the coreg. To avoid this we do not include the headshape points and
        % rely only on the fiducials for the coregistration.
        coreg_settings = struct;
        coreg_settings.D = D.fullfile;
        coreg_settings.mri = sprintf('%sraw_data/sub%03d/highres001.nii',config.datadir,subjects(j))
        coreg_settings.useheadshape = false;
        coreg_settings.forward_meg = 'Single Shell';
        coreg_settings.use_rhino = true;
        coreg_settings.fid.label.nasion='Nasion';
        coreg_settings.fid.label.lpa='LPA';
        coreg_settings.fid.label.rpa='RPA';
        D = osl_headmodel(coreg_settings);

        % Next we generate and save a summary image so we can check each
        % coregistration has been carried out sucessfully.
        h = report.coreg(D);
        report.save_figs(h,outdir,D.fname);
        close(h);

    end
end

%% MEEG Data Preprocessing
%
% All the processing including the electrophysiological data is carried out
% within this loop. Here we apply a range of data preprocessing steps to remove
% artefacts from our data, project the data into source-space and apply a
% parcellation.
%
% Downsampling - Downsample from 1000Hz to 250Hz
% Filtering - We apply a single passband filter to isolate data between 1 and 45Hz
% Bad Segment Detection - An automatic algorithm is used to identify noisy data segments which are removed from subsequent analysis
% Independant Components Analysis - ICA is used to identify artefactual componets by correlation with the EOG and ECG, these are removed from subsequent analysis
% Sensor Normalisation - The Magnetometers and Gradiometers within each dataset are normalised to make their variances comparable.
% Beamforming - An LCMV Beamformer is used to project the sensor data into an 8mm grid in source space
% Parcellation - The source space data is reduced into a network of parcels defined by a NIFTI file
%
% The analysis is carried out on files copied across from spm_sss.  Where
% possible, the processsing is carried out on the files in spm_sss_processed
% 'in place', meaning that they do not generate a new SPM object. The
% downsampling and filtering overwrite the old data file, but the ICA,
% beamforming and parcellation are applied by adding online montages to the SPM
% object.

% Intialise parcellation and study objects
p = parcellation( 'fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm' );
s = study( fullfile(config.analysisdir,'spm_sss'),0);

% Set output and working directories
outdir = fullfile( config.analysisdir,'spm_sss_processed' );
wd = getfullpath('tempdir');
mkdir(wd)
if ~exist( outdir )
    mkdir(outdir);
end

% Preallocate array for ICA sessions that might need manual checking
sessions_to_check = [];

% Main loop through files within the study object
for j = 1:s.n

    % We don't need to repeat files which are already preprocessed
    if exist(fullfile(outdir,s.fnames{j}),'file')
        fprintf('Done: %s\n',fullfile(outdir,s.fnames{j}))
        continue
    else
        fprintf('Todo: %s\n',fullfile(outdir,s.fnames{j}))
    end

    % Load in MEEG data as an SPM object
    D = s.read(j);

    % Downsample and copy, or just copy
    if D.fsample > 250
        D = spm_eeg_downsample(struct('D',D,'fsample_new',250,'prefix',[wd '/'])); % Note - downsampling cannot be done in-place using prefix='', it just fails
    else
        D = D.copy(getfullpath(fullfile(pwd,wd,D.fname))); % Copy into working directory
    end

    % Apply a 1-45Hz passband filter
    D = osl_filter(D,[1 45],'prefix','');

    % Apply automatric bad segment detection
    D = osl_detect_artefacts(D,'badchannels',false);

    % Run ICA artefact detection. This will automatically reject components
    % which have correlations larger than .5 with either of the artefact
    % channels.
    D = osl_africa(D,'used_maxfilter',true);

    % Though the automatic correlations generally work well, we should be
    % careful to check for unsusual datasets and possibly manually correct the
    % automatic assessment. This is particularly important for relatively noisy
    % data or when analysing a new dataset for the first time.
    %
    % Here we will manally inspect the ICA component rejections for any dataset meeting the following criteria
    % # More than 4 ICs rejected
    % # Zero ICs rejected
    % # No component rejected due to correlation with EOG
    % # No component rejected due to correlation with ECG
    artefact_chan_corr_thresh = .5;
    if isempty(D.ica.bad_components) || length(D.ica.bad_components) > 4
        disp('%s components rejected, recommend checking session', length(D.ica.bad_components));
        sessions_to_check = cat(1,sessions_to_check,j);
    elseif max(D.ica.metrics.corr_chan_367_EOG.value) < artefact_chan_corr_thresh || ...
            max(D.ica.metrics.corr_chan_368_EOG.value) < artefact_chan_corr_thresh
        disp('no candidate components for either EOG, recommend checking session');
        sessions_to_check = cat(1,sessions_to_check,j);
    elseif max(D.ica.metrics.corr_chan_369_ECG.value) < artefact_chan_corr_thresh
              disp('no candidate components for  ECG, recommend checking session');
        sessions_to_check = cat(1,sessions_to_check,j);
    end

    % This is where manual Africa would go
    % D = D.montage('remove',1:D.montage('getnumber'));
    % D = osl_africa(D,'do_ident','manual');

    % Normalise sensor types
    S = [];
    S.D = D;
    S.modalities = {'MEGMAG','MEGPLANAR'};
    S.do_plots = 0;
    S.samples2use = good_samples(D,D.indchantype(S.modalities,'GOOD'));
    S.trials = 1;
    S.pca_dim = 99;
    S.force_pca_dim = 0;
    S.normalise_method = 'min_eig';
    D = normalise_sensor_data( S );

    % Run LCMV Beamformer
    D = osl_inverse_model(D,p.template_coordinates,'pca_order',50);

    % Do parcellation
    D = ROInets.get_node_tcs(D,p.parcelflag,'spatialBasis','Giles');

    % Save out
    D = D.montage('switch',0);
    D.copy(fullfile(outdir,s.fnames{j}));

end


%% Epoch information
%
% Finally, we extract the epoch information from the continuous SPM12 files. We
% do not apply the epoching here as the HMM will be run on the continuous data,
% without knowledge of any task structure. The epoch definitions here will
% instead be applied to the HMM state time-courses.

% Preallocate results array
epochinfo = cell(114,1);

% Load in preprocessed data
s = study( fullfile(config.analysisdir,'spm_sss_processed'),0);

for j = 1:s.n

    % Load in MEEG data as an SPM object
    D = s.read(j);

    S2 = struct;
    S2.D = D;

    % We define a wide window to include the participant responses
    pretrig = -1000; % epoch start in ms
    posttrig = 2000; % epoch end in ms
    S2.timewin = [pretrig posttrig];

    % define the trials we want from the event information
    S2.trialdef(1).conditionlabel = 'Famous_first';
    S2.trialdef(1).eventtype = 'STI101_up';
    S2.trialdef(1).eventvalue = [5];
    S2.trialdef(2).conditionlabel = 'Famous_imm';
    S2.trialdef(2).eventtype = 'STI101_up';
    S2.trialdef(2).eventvalue = [6];
    S2.trialdef(3).conditionlabel = 'Famous_last';
    S2.trialdef(3).eventtype = 'STI101_up';
    S2.trialdef(3).eventvalue = [7];

    S2.trialdef(4).conditionlabel = 'Unfamiliar_first';
    S2.trialdef(4).eventtype = 'STI101_up';
    S2.trialdef(4).eventvalue = [13];
    S2.trialdef(5).conditionlabel = 'Unfamiliar_imm';
    S2.trialdef(5).eventtype = 'STI101_up';
    S2.trialdef(5).eventvalue = [14];
    S2.trialdef(6).conditionlabel = 'Unfamiliar_last';
    S2.trialdef(6).eventtype = 'STI101_up';
    S2.trialdef(6).eventvalue = [15];

    S2.trialdef(7).conditionlabel = 'Scrambled_first';
    S2.trialdef(7).eventtype = 'STI101_up';
    S2.trialdef(7).eventvalue = [17];
    S2.trialdef(8).conditionlabel = 'Scrambled_imm';
    S2.trialdef(8).eventtype = 'STI101_up';
    S2.trialdef(8).eventvalue = [18];
    S2.trialdef(9).conditionlabel = 'Scrambled_last';
    S2.trialdef(9).eventtype = 'STI101_up';
    S2.trialdef(9).eventvalue = [19];

    S2.reviewtrials = 0;
    S2.save = 0;
    S2.epochinfo.padding = 0;
    S2.event = D.events;
    S2.fsample = D.fsample;
    S2.timeonset = D.timeonset;

    [epochinfo{j}.trl, epochinfo{j}.conditionlabels] = spm_eeg_definetrial(S2);

end

% Save epoch information for use after HMM inference
save( fullfile(config.analysisdir,'spm_sss_processed','epochinfo.mat'), 'epochinfo', '-v7.3');

%% Run a final sanity check to make sure everything is in place
%

missing_runs = utils.check_preproc_data;
