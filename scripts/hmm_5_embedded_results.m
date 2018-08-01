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
method = 'embedded';
K = 6;

% Find HMM directory
base = fullfile( config.analysisdir, sprintf('%s_hmm',method));
mkdir( fullfile( base, 'figures' ) );

% Load in HMM results
load( fullfile(base, sprintf('%s_HMM_K%s.mat',method,num2str(K))) );

% Load in run indices
load( fullfile(base, sprintf('%s_hmm_data_flipped.mat',method)), 'R','B','runlen' );

% Load in epoch info
load( fullfile( config.analysisdir, 'spm_sss_processed','epochinfo.mat' ) );

% Create basepath for saving results
savebase = fullfile( config.analysisdir, sprintf('%s_hmm',method),'figures',sprintf('%s_HMM_K6',method));
if ~exist( savebase )
    mkdir(savebase);
end

% account for delay embedding in state gammas
pad_options.embeddedlags = -7:7;
Gamma = padGamma(Gamma, T, pad_options);

%% Temporal statistics
%
% Here we compute the global temporal statistics from the Time-Delay-Embedded
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


% Plot temporal stats
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

% Build first level design matrix
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

% Define first level contrasts
first_level_contrasts = zeros(3,4);
first_level_contrasts(1,:) = [1 0 0 0]; % Grand Mean
first_level_contrasts(2,:) = [0 1 1 -2]; % Faces>Non-Faces
first_level_contrasts(3,:) = [0 1 -1 0]; % Famous>Unfamiliar

% Group level design matrix and contrasts - just the mean of the first
% levels
group_level_design_matrix = ones(19,1);
group_level_contrasts = 1;

% Estimate GLM
[copes,thresh_glm] = utils.run_group_glm( erg_bl(:,250:650,:,:),...
                            first_level_design_matrix,first_level_contrasts,...
                            group_level_design_matrix,group_level_contrasts,...
                            1000);

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

%% State descriptions
%
% The state descriptions for the embedded HMM are taken from the state-wise
% multitaper estimation (in contrast to the envelope, where we made direct use of
% the observation model).
%
% To aid visualisation, we compute two Non-Negative Matrix Factorisations so we
% can avoid setting arbitrary frequency bands. This does not change the results
% of the HMM, but simply provides a data-driven method for splitting the
% continuous spectrum into a set of peaks.

parc = parcellation('fmri_d100_parcellation_with_PCC_tighterMay15_v2_8mm');

%% Broadband power plots
mt_outfile = fullfile( config.analysisdir, 'embedded_hmm', sprintf('embedded_HMM_K%d_spectra',K));
load( mt_outfile )

net_mean = zeros(39,size(psd,2));
for kk = 1:size(psd,2)
	tmp = squeeze( mean(mean(abs(psd(:,kk,1:29,:,:)),3),1) );
    net_mean(:,kk) = zscore(diag(tmp));
end

% visualise state in OSLEYES
parc.osleyes(net_mean);

% Optionally save a nifti of the results, these are used to generate the
% figures in the paper via HCP Workbench
parc.savenii( net_mean, [savebase '_meanactivations_broadband']);

%% Broadband glass brain networks

% Compute GMM Threshold
H = [];
for kk = 1:6
    G = squeeze( mean(mean(abs(psd(:,kk,1:29,:,:)),3),1) );
    G = G-diag(diag(G));
    H = [H, reshape(G(G~=0),1,[])];
end

S2 = struct;
S2.do_fischer_xform = false;
S2.do_plots = true;
S2.data = H;
[ graphgmm_res ] = teh_graph_gmm_fit( S2 );

% Plot glass brain networks
for kk = 1:6

    G = squeeze( mean(mean(abs(psd(:,kk,1:29,:,:)),3),1) );
    G(G<graphgmm_res.orig_th) = 0;
    if sum(sum(G))>0

        h = parc.plot_network( G );
        h.Parent.Parent.Position(3) = 420;
        h.Parent.FontSize = 18;
        view([0 90]);
        zoom(1);

    end
end

%% Spectral Mode NNMF
%

% Compute the NNMF
mt_outfile = fullfile( config.analysisdir, 'embedded_hmm', 'embedded_HMM_K6_spectra');
load( mt_outfile, 'psd' );

S = [];
S.psds = psd(:,:,1:29,:,:);
S.maxP=4;
S.maxPcoh=4;
S.do_plots = true;
nnmf_res = nnmf_res = wh.run_nnmf( S, 20 );

nnmf_outfile = fullfile( config.analysisdir, 'embedded_hmm', 'embedded_HMM_K6_nnmf');
save(nnmf_outfile,'nnmf_res')

% Visualise the mode shapes
for ii = 1:3
    figure('Position',[100 100*ii 256 256])
    h = area(nnmf_res.nnmf_coh_specs(ii,:));
    h.FaceAlpha = .5;
    h.FaceColor = [.5 .5 .5];
    grid on;axis('tight');
    set(gca,'YTickLabel',[],'FontSize',14);
end


%% Spectral Mode Power Plots
nnmf_outfile = fullfile( config.analysisdir, 'embedded_hmm', 'embedded_HMM_K6_nnmf');
load( nnmf_outfile );

net_mean = zeros(39,4,size(Gamma,2));
thresh_mean = zeros(size(Gamma,2),1);
for k = 1:size(Gamma,2)
    net_mean(:,:,k) = squeeze(nnmf_res.nnmf_psd_maps(k,:,:))';
end

for k = 1:size(Gamma,2)
    dat = net_mean(:,1:3,k);
    parc.osleyes( dat )
end


%% Spectral Mode Network Plots

for kk = 1:6

    G = squeeze(sum(nnmf_res.nnmf_coh_maps(kk,3,:,:),2));

    S2 = struct;
    S2.do_fischer_xform = false;
    S2.do_plots = false;
    S2.data = nnmf_res.nnmf_coh_maps(kk,1:3,:,:);
    S2.data = S2.data(S2.data~=0);
    [ graphgmm_res ] = teh_graph_gmm_fit( S2 );

    for jj = 1:3
        G = squeeze(nnmf_res.nnmf_coh_maps(kk,jj,:,:));
        G(G<graphgmm_res.orig_th) = 0;

        if sum(G(:)) > 0
            h = parc.plot_network(G);
            h.Parent.Parent.Position(3) = 420;
            h.Parent.FontSize = 18;
            set(gca,'CLim',[graphgmm_res.orig_th max(S2.data)]);
            view([0 90]);
            zoom(1);
        else
            disp('No connections survived thresholding');
        end

    end
end


%% HMM regularised TF plots
%
% Finally we compute the HMM regularised TF plots utils.ch show the induced power
% responses as modelled by the HMM for a given parcel.

% nodes of interest
node = 9; % this can be changed to select another node from the parcellation

tf_hmm = zeros(K,39,751,19);

% Generate HMM Regularised TF Response
for ii = 1:19
    for jj = 1:6
        tf_hmm(jj,:,:,ii) = (abs(squeeze(psd(ii,jj,:,node,node)))*squeeze(nanmean(erg_bl(jj,:,:,ii),3))) .* FO(ii,jj);
    end
end

% Make a plot summarising the TF response, task-evoked gammas and state-wise spectra
freq_vect = 1:39;
figure('Position',[100 100 768 768])
% Spectra
ax(1) = axes('Position',[.1 .3 .2 .6]);hold on
for ii = 1:6
    plot(ax(1),abs(squeeze(nanmean(psd(:,ii,:,node,node),1))),freq_vect,...
                'linewidth',2,'Color',set1_cols{ii});
end
set(ax(1),'Xdir','reverse','XTickLabel',[]);grid on
ylabel('Frequency (Hz)')

% ERG
ax(2) = axes('Position',[.3 .1 .6 .2]);hold on
dat = nanmean(nanmean(erg_bl,4),3);% - repmat(erg_baselines,1,751);
for ii = 1:6
    plot(ax(2), time_vect_adj,dat(ii,:),'linewidth',2,'Color',set1_cols{ii})
end
axis('tight');grid on
hold on;plot([0 0],ylim,'k--','linewidth',1.5);
plot([.932 .932],ylim,'k--','linewidth',1.5);
xlim([-.1 2])
xlabel('Time (secs)')
ylabel('Relative Occupancy')

% Contour plot
ax(3) = axes('Position',[.3 .3 .6 .6]);
dat = squeeze(sum(mean(tf_hmm(:,:,:,:),4),1));
bc = squeeze(mean(dat(:,225:250),2));
dat = dat-repmat(bc,1,751);
contourf(ax(3), time_vect_adj,freq_vect,dat,24,'linestyle','none')
grid on
c = colorbar;c.Position(1) = .91;
set(ax(3),'XTickLabel',[],'YTickLabel',[])
hold on;plot([0 0],ylim,'k--','linewidth',1.5);
plot([.932 .932],ylim,'k--','linewidth',1.5);
xlim([-.1 2])
utils.set_redblue_colourmap(ax(3),min(min(dat)), max(max(dat)))


