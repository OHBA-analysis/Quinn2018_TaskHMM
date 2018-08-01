function missing_files = check_raw_data()

    config = wh.get_studydetails();

    missing_files = {};
    for subj = 1:19

        subjname = sprintf( 'sub0%02d', subj );
        subjpath = fullfile(config.datadir,'raw_data',subjname);

        % Check fif locations
        for run = 1:6
            runpath = fullfile( subjpath, sprintf( 'run_0%d_sss.fif', run ) );
            if ~exist( runpath )
                missing_files{end+1} = runpath;
            end
        end

        % Check structural location
        runpath = fullfile( subjpath, 'highres001.nii' );
        if ~exist( runpath )
            missing_files{end+1} =  runpath;
        end

    end

    % Inform user whether files are missing.
    if isempty( missing_files)
        disp('All raw data found in studydir');
    else
        msg = sprintf( '%d raw data files missing', length(missing_files) );
        warning(msg);
	for ii = 1:length(missing_files)
	    disp(missing_files{ii});
	end
    end
