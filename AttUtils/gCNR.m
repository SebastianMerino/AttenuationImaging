function out = gCNR(data1,data2, varargin)
%GCNR Summary of this function goes here
%   Detailed explanation goes here
    if nargin == 2
        [~,edg] = histcounts([data1;data2],32);
    else
        edg = varargin{1};
    end
    f = histcounts(data1,edg, 'Normalization','probability');
    g = histcounts(data2,edg, 'Normalization','probability');
    out = 1 - sum(min(f, g));
end

