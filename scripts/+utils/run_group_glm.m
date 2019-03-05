function [group_copes, thresh] = run_group_glm( data, ...
                            first_level_design_matrix,first_level_contrasts,...
                            group_level_design_matrix,group_level_contrasts,...
                            nperms, use_tstat,use_flame)
%%function [group_copes, thresh] = run_group_glm( data, ...
%                            first_level_design_matrix,first_level_contrasts,...
%                            group_level_design_matrix,group_level_contrasts,...
%                            nperms)
%
% data - [states,time,trials,subjects]
% first_level_design_matrix - [trials, subject predictors,subjects]
% first_level_contrasts - [subject predictors, subject contrasts]
% group_level_design_matrix - [subjects, group predictors]
% group_level_contrasts - [group predictors, group contrasts]
%

if nargin < 8 || isempty(use_flame)
    use_flame=1;
end

if nargin < 7 || isempty(use_tstat)
    use_tstat=0;
end

if nargin < 6 || isempty(nperms)
    nperms = 0;
end

% Meta info
[nstates,nsamples,ntrials,nsubjs] = size(data);
ncontrasts = size(first_level_contrasts,1);

% Preallocate data arrays
copes = zeros(nstates,ncontrasts,nsamples,nsubjs);
varcopes = zeros(nstates,ncontrasts,nsamples,nsubjs);
tstats = zeros(nstates,ncontrasts,nsamples,nsubjs);

group_copes = zeros(nstates,ncontrasts,nsamples);
group_tstats = zeros(nstates,ncontrasts,nsamples);

% Compute first and group level GLM per state.
fprintf('\nComputing GLM\n');
for ii = 1:nstates
    for jj = 1:nsubjs
        [copes(ii,:,:,jj),varcopes(ii,:,:,jj),tstats(ii,:,:,jj)] = utils.run_glm( squeeze(data(ii,:,:,jj))',...
                                                first_level_design_matrix(:,:,jj),first_level_contrasts);
    end

    if use_flame==1
        % Group level estimate
        for jj = 1:ncontrasts
            [group_copes(ii,jj,:),~,group_tstats(ii,jj,:)] = utils.run_flame( squeeze(copes(ii,jj,:,:))',...
                                                                       squeeze(varcopes(ii,jj,:,:))',...
                                                group_level_design_matrix,group_level_contrasts);
        end
    else
        for jj = 1:ncontrasts
            [group_copes(ii,jj,:),~,group_tstats(ii,jj,:)] = utils.run_glm( squeeze(copes(ii,jj,:,:))',...
                                                group_level_design_matrix,group_level_contrasts);
        end
    end
end

% Compute group level stats with sign-flipping permutations
if nperms > 0
    msg_base = 'Computing permutation %d of %d';
    msg = sprintf( msg_base, 1, nperms);
    fprintf(msg);

    % Preallocate null distributions and add observed vales
    null_dist = zeros(nperms,size(first_level_contrasts,1));
    null_dist(1,:) = squeeze(max(max(group_copes,[],3),[],1));

    % Run permutations
    for ii = 2:nperms
        fprintf(repmat('\b',1,length(msg)))
        msg = sprintf( msg_base, ii, nperms);
        fprintf(msg);

        perm_copes = zeros(nstates,ncontrasts,size(copes,3));
        perms = sign(randn( size(copes,4),1 ));
        for jj = 1:nstates% per state

            if use_flame==1
                % Group level estimate
                for kk = 1:ncontrasts
                    [perm_copes(jj,kk,:),~,~] = utils.run_flame( squeeze(copes(jj,kk,:,:))',...
                        squeeze(varcopes(jj,kk,:,:))',...
                        group_level_design_matrix.*perms,group_level_contrasts);
                end
            else
                for kk = 1:ncontrasts
                    [perm_copes(jj,kk,:),~,~] = utils.run_glm( squeeze(copes(jj,kk,:,:))',...
                        group_level_design_matrix.*perms,group_level_contrasts);
                end
            end

        end
        null_dist(ii,:) = squeeze(max(max(perm_copes,[],3),[],1));
    end

    fprintf('\n')
    % Estimate threshold
    thresh = prctile(null_dist,95,1);
else
    thresh = [];
end
