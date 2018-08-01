function [ trial_data ] = load_epoch_results( data, epochinfo, runlen, B, R)

% meta information
nchannels = size(data,1); % this could also be states
nsamples = epochinfo{1}.trl(1,2) - epochinfo{1}.trl(1,1) + 1;

% return array
trial_data = zeros(nchannels,nsamples,148,114)*nan;

% main loop
for ii = 1:19
    for jj = 1:6
        ind = ((ii-1)*6)+jj;

        % get bad samples
        good_inds=setdiff(1:runlen(ind),B{ind});

        % extract subject data, accounting for bad samples
        subj_data = zeros( nchannels,runlen(ind) )*nan;
        subj_data(:,good_inds) = data(:,R(ind,1):R(ind,2));

        % main epoch loop
        for kk = 1:size(epochinfo{ind}.trl,1)
            % get trial start and stop samples
            start = epochinfo{ind}.trl(kk,1);
            stop = epochinfo{ind}.trl(kk,2);

            if stop > size(subj_data,2) || start > size(subj_data,2)
                continue % not all sessions have 148 trials
            else
               trial_data(:,:,kk,ind) = subj_data(:,start:stop);
            end
        end
    end
end
