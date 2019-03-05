function out  = check_preproc_data()
%%function [missing_files,missing_proc_files] = check_preproc_data()
%
% Helper function to check that the preprocessing has completed successfully

    % Get study information
    config = utils.get_studydetails();

    % Preallocate bad file arrays
    out = [];
    out.missing_files = {};
    out.missing_proc_files = {};
    out.missing_epochinfo = {};

    % Main loop
    for subj = 1:19
        fprintf('Checking subject: %d\n',subj);
        for run = 1:6

            % Generate run name
            runname = sprintf( 'spm_sub%02d_run_%02d.mat', subj,run );

            %% Check spm_sss directory
            rawpath = fullfile( config.analysisdir, 'spm_sss',runname );
            if ~exist( rawpath )
                out.missing_files{end+1} = rawpath;
                continue
            end
            D = spm_eeg_load(rawpath);

            % Check coregistration has run
            if ~isfield( D, 'inv' )
                fprintf( 'Coreg out.missing for file: %s\n',runpath);
            end

            %% Check spm_sss_preprocessed directory
            preprocpath = fullfile( config.analysisdir, 'spm_sss_processed',runname );
            if ~exist( preprocpath )
                out.missing_files{end+1} = preprocpath;
                continue
            end
            D = spm_eeg_load(preprocpath);

            % Check that we have 5 montages
            if D.montage('getnumber') ~= 5
                out.missing_proc_files{end+1} = preprocpat;
            end

            % Check that the fifth montage has 39 channels
            D = D.montage('switch',5);
            if D.nchannels ~= 39
                out.missing_proc_files{end+1} = preprocpat;
            end

        end
    end

    fprintf('\n\n');

    %% Check epoch information
    load( fullfile( config.analysisdir,'spm_sss_processed','epochinfo.mat') );
    missing_epochs = cellfun( @isempty, epochinfo );
    if sum( missing_epochs ) == 0
        disp('All epochinfo present in epochinfo file');
    elseif sum( missing_epochs ) > 0
        msg = sprintf( '%d subjects epochinfo out.missing', length(find(missing_epochs)) );
        out.missing_epochs = find(missing_epochs);
    end

    % Inform user utils.ther files are out.missing.
    if isempty( out.missing_files)
        disp('All coregistered data found in spm_sss dir');
    else
        msg = sprintf( '%d raw data files out.missing', length(out.missing_files) );
        warning(msg);
	    for ii = 1:length(out.missing_files)
	        disp(out.missing_files{ii});
	    end
    end

    % Inform user utils.ther files are out.missing.
    if isempty( out.missing_proc_files)
        disp('All parcel-network data found in spm_sss_processed dir');
    else
        msg = sprintf( '%d raw data files out.missing', length(out.missing_proc_files) );
        warning(msg);
	    for ii = 1:length(out.missing_proc_files)
	        disp(out.missing_proc_files{ii});
	    end
    end
