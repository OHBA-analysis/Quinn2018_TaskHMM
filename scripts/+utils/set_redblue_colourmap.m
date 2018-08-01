function cm = set_redblue_colourmap( ax, clmin, clmax, N )
%%function cm = set_redblue_colourmap( ax, clmin, clmax, N )
%
%

if nargin < 4 || isempty(N)
    N = 64;
end

% scale to sensible order of magnitude
min_scale = floor(log10(abs(clmin)));
clmin = clmin*10.^-min_scale;

max_scale = floor(log10(abs(clmax)));
clmax = clmax*10.^-max_scale;

% Round up to nearest 1
clmin = -round(abs(clmin),1);
clmax = round(abs(clmax),1);

% make largest order of magnitude the ones unit
if abs(max_scale) < abs(min_scale)
    clmin = clmin ./ 10*(max_scale-min_scale);
elseif abs(max_scale) > abs(min_scale)
    clmax = clmax ./ 10*(min_scale-max_scale);
end

% get positive values
max_ratio = abs(clmax)/ ( abs(clmin) + abs(clmax) );
max_lims = round(N*max_ratio);
c1 = cat(1,linspace(.5,1,max_lims),linspace(0,1,max_lims),linspace(0,1,max_lims));

% get negative values
min_ratio = abs(clmin)/( abs(clmax)+abs(clmin) );
min_lims = round(N*min_ratio);
c2 = cat(1,linspace(1,0,min_lims),linspace(1,0,min_lims),linspace(1,.5,min_lims));

% set map
cm = flipud(cat(2,c1,c2)');
colormap(ax,cm);
if abs(max_scale) < abs(min_scale)
    caxis(ax,[clmin.*10^(1+min_scale) clmax.*10^max_scale]);
elseif abs(max_scale) > abs(min_scale)
    caxis(ax,[clmin.*10^min_scale clmax.*10^(1+max_scale)]);
else
    caxis(ax,[clmin.*10^min_scale clmax.*10^(max_scale)]);
end