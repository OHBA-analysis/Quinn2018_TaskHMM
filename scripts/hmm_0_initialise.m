% specify the location of this download in the line below:
% eg download_path = '/home/disk3/ajquinn/HMM_frontiers_download';
download_path = ''; % <- change this line

if isempty(download_path)
    error('Please specify the location of the download path inside hmm_0_initialise.m');
end

if exist( download_path, 'dir' ) == 0
    error('Download path in hmm_0_initialise.m not found: %s', download_path);
end

% Add study utilities
addpath( fullfile( download_path , 'scripts' ) );

% Get study paths
config = utils.get_studydetails;

% Initialise OSL
addpath( fullfile( config.scriptdir,'toolboxes','osl','osl-core') );
osl_startup

% Add netlab and fmt
addpath( fullfile(osldir,'ohba-external','netlab3.3','netlab') );
addpath( fullfile(osldir,'ohba-external','fmt') );

% Add HMM-MAR to path
addpath(genpath( fullfile( config.scriptdir,'toolboxes','HMM-MAR-master') ));

% Add distibutionplot to path
addpath(genpath( fullfile( config.scriptdir,'toolboxes','distributionPlot') ));
