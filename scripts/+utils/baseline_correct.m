function [data] = baseline_correct( data, baseline_inds )
%%function [data] = baseline_correct( data, baseline_inds )
%
% baseline corrects the second dimenion

for ii = 1:size(data,1)
    for jj = 1:size(data,4)
        bl = squeeze(nanmean(nanmean(data(ii,baseline_inds,:,jj),2),3));
        data(ii,:,:,jj) = data(ii,:,:,jj) - bl;
    end
end