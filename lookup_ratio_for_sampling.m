function [ keepRatio ] = lookup_ratio_for_sampling( n, hem )
%LOOKUP_RATIO_FOR_SAMPLING Keepratio for sampling
%   Returns KEEPRATIO, the decimation rate for the mesh subsampling
%   process. It is a number less than 1, indicates the percentage of the 
%   remaining vertices after the subsampling. The values below have been 
%   precisely tuned for the HCP surfaces, for the given number of vertices, 
%   N and the hemisphere, HEM. For other Ns, please work out the exact 
%   number with the following formula: 
%
%   keepratio = N / 29696   if HEM is 'L' 
%   keepratio = N / 29716   if HEM is 'R'
%       
%                  Left        Right
%   n =  500      .01683      .01668
%   n = 1000      .03353      .03366
%   n = 1500      .05065      .05060
%   n = 2000      .06725      .06718


if strcmp(hem, 'L')
    if n == 500
        keepRatio = .01683;
    elseif n == 1000
        keepRatio = .03353;
    elseif n == 1500
        keepRatio = .05065;
    elseif n == 2000;
        keepRatio = .06725;
    else
        keepRatio = n / 29696;
    end
elseif strcmp(hem, 'R')
    if n == 500
        keepRatio = .01668;
    elseif n == 1000
        keepRatio = .03366;
    elseif n == 1500
        keepRatio = .05060;
    elseif n == 2000;
        keepRatio = .06718;
    else
        keepRatio = n / 29696;
    end
else
    error('Unrecognized hemisphere...')
end


