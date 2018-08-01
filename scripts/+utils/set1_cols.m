function set1 = set1_cols

% Colourscheme based on the ColorBrewer Set 1

% http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=8
set1 = { [228,26,28], [55,126,184], [77,175,74], [152,78,163],...
              [255,127,0], [247,129,191], [255,255,51],[166,86,40] };
          
% normalise and convert back to cell array
set1 = cellfun(@(x) x./255, set1,'UniformOutput', false );