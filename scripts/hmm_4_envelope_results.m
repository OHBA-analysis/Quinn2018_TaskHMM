%% Overview
%
% Here we summarise the results of the HMM inference and compute the
% task-evoked GLM statistics

config = utils.get_studydetails;

% Define colours to use in state plots
set1_cols = utils.set1_cols;

% Define sample rate
sample_rate = 250;

%% Load in results from envelope data

% Meta data
method = 'envelope';
K = 6;

% Find HMM directory
base = fullfile( config.analysisdir, 'envelope_hmm');
mkdir( fullfile( base, 'figures' ) );

% Load in HMM results
load( fullfile(base, sprintf('envelope_HMM_K%s.mat',num2str(K))) );

% Load in run indices
load( fullfile(base, 'envelope_hmm_data.mat'), 'R','B','runlen' );

% Load in epoch info
load( fullfile( config.analysisdir, 'spm_sss_processed','epochinfo.mat' ) );

% Create basepath for saving results
savebase = fullfile( config.analysisdir, 'envelope_hmm','figures','envelope_HMM_K6');

%% Temporal statistics
%
% Here we compute the global temporal statistics from the Amplitude Envelope
% HMM. These are computed per subject and visualised as violin plots

scan_T = [R(1,2) diff(R(:,2))']; % Indexing individual scan sessions
subj_T = sum(reshape(scan_T,6,[])); % Indexing individal subjects

% Compute temporal stats

% Fractional Occupancy is the proportion of time spent in each state
FO = getFractionalOccupancy( Gamma, subj_T, 2);
% Interval Time is the time between subsequent visits to a state
IT = getStateIntervalTimes( Gamma, subj_T, []);
ITmerged = cellfun(@mean,IT);clear IT
% Life Times (or Dwell Times) is the duration of visits to a state
LT = getStateLifeTimes( Gamma, subj_T, []);
LTmerged = cellfun(@mean,LT); clear LT

% Make summary figures
fontsize = 18;

figure;subplot(111);
distributionPlot(FO,'showMM',2,'color',{set1_cols{1:size(FO,2)}});
set(gca,'YLim',[0 1],'FontSize',fontsize)
title('Fractional Occupancy');xlabel('State');ylabel('Proportion');grid on;
print([savebase '_temporalstats_FO'],'-depsc')

figure;subplot(111);
distributionPlot(LTmerged ./ sample_rate * 1000,'showMM',2,'color',{set1_cols{1:size(FO,2)}})
title('Life Times');xlabel('State');ylabel('Time (ms)');grid on;
set(gca,'YLim',[0 300],'FontSize',fontsize,'FontSize',fontsize)
print([savebase '_temporalstats_LT'],'-depsc')

figure;subplot(111);
distributionPlot(ITmerged ./ sample_rate,'showMM',2,'color',{set1_cols{1:size(FO,2)}})
title('Interval Times');xlabel('State');ylabel('Time (secs)');grid on
set(gca,'YLim',[0 3],'FontSize',fontsize)
print([savebase '_temporalstats_IT'],'-depsc')

%% Extract Event Related Field and Event Related Gamma

% Create time vector and adjust for projector lag
time_vect = linspace(-1,2,751);
time_vect_adj = time_vect - .032; % account for projector lag

% Apply epoching to posterior probabilities, merge subjects and baseline
% correct
erg = utils.load_epoch_results( Gamma', epochinfo, runlen, B, R);
erg = utils.merge_sessions( erg, repelem(1:19,6) );
erg_bl = utils.baseline_correct( erg, 225:250 );

%% Compute two-level GLM
%
% This cell computes the two level GLM from the baseline corrected task-evoked state probabilities.

% extract conditions for each subject and session.
conds={'Famous','Unfamiliar','Scrambled'};
cond_inds = cell(19,3);
for ii = 1:19
    for jj = 1:6
         session = ((ii-1)*6)+jj;
         condlabels = epochinfo{session}.conditionlabels;
         for kk = 1:3
             cond_inds{ii,kk} = cat(2,cond_inds{ii,kk}, find(~cellfun(@isempty, strfind(condlabels,conds{kk}))) + ((jj-1)*150));
         end
    end
end

% Build first level design matrix for each subject
first_level_design_matrix = zeros(size(erg_bl,3),4,19);
first_level_design_matrix(:,1,:) = 1; % Mean term
for ii = 1:19
    first_level_design_matrix(cond_inds{ii,1},2,ii) = 1; % Famous Faces
    first_level_design_matrix(cond_inds{ii,2},3,ii) = 1; % Unfamiliar Faces
    first_level_design_matrix(cond_inds{ii,3},4,ii) = 1; % Scrambled Faces

    % Remove the mean from non-constant regressors as we have an explicit
    % constant regressor and contrast for the mean
    for jj = 2:4
        first_level_design_matrix(:,jj,ii) = demean(first_level_design_matrix(:,jj,ii));
    end
end

% Define first level contrasts, these are constant across subjects
first_level_contrasts = zeros(3,4);
first_level_contrasts(1,:) = [1 0 0 0]; % Grand Mean
first_level_contrasts(2,:) = [0 1 1 -2]; % Faces>Non-Faces
first_level_contrasts(3,:) = [0 1 -1 0]; % Famous>Unfamiliar

% Group level design matrix and contrasts - just the mean of the first levels.
group_level_design_matrix = ones(19,1);
group_level_contrasts = 1;

% Estimate GLM
[copes,thresh_glm] = utils.run_group_glm( erg_bl(:,250:650,:,:),...
                            first_level_design_matrix,first_level_contrasts,...
                            group_level_design_matrix,group_level_contrasts,...
                            1000,0,0);

%% Plot GLM results

% Grand mean figure
figure;subplot(111);hold on;grid on
t = time_vect_adj(250:650);

for ii = 1:6
    plot(t,squeeze(copes(ii,1,:))','Color',set1_cols{ii},'linewidth',2)

    if sum(copes(ii,1,:)>thresh_glm(1)) > 0
        sig_inds = find(squeeze(copes(ii,1,:))>thresh_glm(1));
        t2 = ones(size(t))*nan;t2(sig_inds) = t(sig_inds);
        x2 = ones(size(t))*nan;x2(sig_inds) = -.08;
        plot(t2,x2,'Color',set1_cols{ii},'linewidth',5);
    end

end

plot([0 0],ylim,'k--')
annotation('textbox',...
        [.16 .89 .1 .03],...
        'String','Stimulus Onset',...
        'FontSize',12,...
        'FontWeight','bold',...
        'EdgeColor',[1 1 1],...
        'LineWidth',3,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
plot([.932 .932],ylim,'k--')
annotation('textbox',...
        [.61 .89 .2 .03],...
        'String','Average Reaction Time',...
        'FontSize',12,...
        'FontWeight','bold',...
        'EdgeColor',[1 1 1],...
        'LineWidth',3,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
xlim([time_vect_adj(250) time_vect_adj(650)])
set(gca,'XTick',0:.2:1.5,'FontSize',18)
ylim([-.2 .2])
xlabel('Time (seconds)')
title('Grand Average')
ylabel('Task Evoked Occupancy')
print([savebase '_FOglm_mean'],'-depsc')

% Faces>Scrambled Faces figure
figure;
subplot(111);hold on;grid on
t = time_vect_adj(250:650);

for ii = 1:6
    plot(t,squeeze(copes(ii,2,:))','Color',set1_cols{ii},'linewidth',2)

    if sum(copes(ii,2,:)>thresh_glm(2)) > 0
        sig_inds = find(squeeze(copes(ii,2,:))>thresh_glm(2));
        t2 = ones(size(t))*nan;t2(sig_inds) = t(sig_inds);
        x2 = ones(size(t))*nan;x2(sig_inds) = -.08;
        plot(t2,x2,'Color',set1_cols{ii},'linewidth',5);
    end

end
ylim([-.1 .125])

plot([0 0],ylim,'k--')
annotation('textbox',...
        [.16 .89 .1 .03],...
        'String','Stimulus Onset',...
        'FontSize',12,...
        'FontWeight','bold',...
        'EdgeColor',[1 1 1],...
        'LineWidth',3,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
plot([.932 .932],ylim,'k--')
annotation('textbox',...
        [.61 .89 .2 .03],...
        'String','Average Reaction Time',...
        'FontSize',12,...
        'FontWeight','bold',...
        'EdgeColor',[1 1 1],...
        'LineWidth',3,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
xlim([time_vect_adj(250) time_vect_adj(650)])
set(gca,'XTick',0:.2:1.5,'FontSize',18)
xlabel('Time (seconds)')
ylabel('Task Evoked Occupancy')
title('Faces > Scrambled Faces');
print([savebase '_FOglm_face'],'-depsc')

% Famous>Unfamiliar Figure
figure;
subplot(111);hold on;grid on
t = time_vect_adj(250:650);

for ii = 1:6
    plot(t,squeeze(copes(ii,3,:))','Color',set1_cols{ii},'linewidth',2)

    if sum(copes(ii,3,:)>thresh_glm(3)) > 0
        sig_inds = find(squeeze(copes(ii,3,:))>thresh_glm(3));
        t2 = ones(size(t))*nan;t2(sig_inds) = t(sig_inds);
        x2 = ones(size(t))*nan;x2(sig_inds) = -.08;
        plot(t2,x2,'Color',set1_cols{ii},'linewidth',5);
    end

end
ylim([-.1 .125])

plot([0 0],ylim,'k--')
annotation('textbox',...
        [.16 .89 .1 .03],...
        'String','Stimulus Onset',...
        'FontSize',12,...
        'FontWeight','bold',...
        'EdgeColor',[1 1 1],...
        'LineWidth',3,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
plot([.932 .932],ylim,'k--')
annotation('textbox',...
        [.61 .89 .2 .03],...
        'String','Average Reaction Time',...
        'FontSize',12,...
        'FontWeight','bold',...
        'EdgeColor',[1 1 1],...
        'LineWidth',3,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
xlim([time_vect_adj(250) time_vect_adj(650)])
set(gca,'XTick',0:.2:1.5,'FontSize',18)
xlabel('Time (seconds)')
ylabel('Task Evoked Occupancy')
title('Famous > Unfamiliar')
print([savebase '_FOglm_famous'],'-depsc')

%% Mean Activation Maps
%
% For the envelope HMM, we take the summary of each state directly from the
% observtion models stored in hmm.state, this provides the information for both
% the mean and functional connectivity results.

% load parcellation
parc = parcellation('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm');

% Activation maps are normalised within each state to allow for simple visualisation of the states topology
net_mean = zeros(39,size(Gamma,2));
for k = 1:size(Gamma,2)

    net_mean(:,k) = zscore( diag(hmm.state(k).Omega.Gam_rate) ./ hmm.state(k).Omega.Gam_shape );

end

% visualise state in OSLEYES
parc.osleyes(net_mean);

% Optionally save a nifti of the results, these are used to generate the
% figures in the paper via HCP Workbench
parc.savenii( net_mean, [savebase '_meanactivations']);


%% Node Weight Maps

% As with the mean activation maps, the node weights are normalised to aid visualisation
net_nw = zeros(39,size(Gamma,2));
thresh_mean = zeros(size(Gamma,2),1);
for k = 1:size(Gamma,2)

    G = hmm.state(k).Omega.Gam_rate ./ hmm.state(k).Omega.Gam_shape;
    G = G - diag(diag(G));
    nw = sum(G,1)' + sum(G,2);
    net_nw(:,k) = zscore(nw);

end

% visualise state in OSLEYES
parc.osleyes(net_nw);

% Optionally save a nifti of the results
parc.savenii( net_nw, [savebase '_nodeweights']);
