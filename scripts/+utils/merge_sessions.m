function data_merged = merge_sessions( data, session_subj_mapping )
%function data = merge_sessions( data, session_subj_mapping )
%
% last dimension of data should be sessions
% session_subj_mapping is repelem(1:19,6);

nsubjs = length(unique( session_subj_mapping ));
nsess = sum( session_subj_mapping==1 );
ntrials = size(data,3)*nsess;

data_merged = zeros( size(data,1), size(data,2), ntrials, nsubjs );

for ii = 1:nsubjs

    data_merged(:,:,:,ii) = reshape( data(:,:,:,session_subj_mapping==ii), size(data,1),size(data,2),ntrials );

end



