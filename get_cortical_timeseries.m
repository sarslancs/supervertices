function [ timeseries ] = get_cortical_timeseries( dtseries, hem )
%GET_CORTICAL_TIMESERIES Get timeseries for the given hemisphere or
%cerebellum
%   [ TIMESERIES ] = GET_CORTICAL_TIMESERIES( DTSERIES, HEM ) returns a 
%   n x d timeseries matrix, in which n is the number of vertices and d is
%   the length of timeseries. HEM is a char that specifies the area of the 
%   brain that has been requested ('L' = Left, 'R' = Right, 'C' = 
%   Cerebellum). DTSERIES Should be a dense cifti matrix in the size of 
%   91282 x d.

if strcmp(hem, 'L')
    timeseries = dtseries(1:29696,:);
elseif strcmp(hem, 'R')
    timeseries = dtseries(29697:59412,:);
else
    timeseries = dtseries(59413:end,:);
end



