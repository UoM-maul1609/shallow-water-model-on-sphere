function [lon_plot,data_plot,lon_native,data_native] = cyclic_lon(lon,data)
%CYCLIC_LON Prepare periodic longitude data for plotting.
%   Handles both the current model output (N distinct samples on [0,2*pi))
%   and legacy output containing both 0 and 2*pi. DATA must have longitude
%   along its second dimension.

lon = lon(:);
if size(data,2) ~= numel(lon)
    error('cyclic_lon:dimensionMismatch', ...
        'DATA must have longitude along dimension 2 (%d columns expected, got %d).', ...
        numel(lon), size(data,2));
end

isdup = false;
if numel(lon) > 1
    dlon = median(abs(diff(lon)));
    tol = max(1.e-10, dlon.*1.e-6);
    isdup = abs(abs(lon(end)-lon(1)) - 2*pi) <= tol;
end

if isdup
    lon_native = lon(1:end-1);
    data_native = data(:,1:end-1);
else
    lon_native = lon;
    data_native = data;
end

if numel(lon_native) > 1
    direction = sign(median(diff(lon_native)));
    if direction == 0
        direction = 1;
    end
else
    direction = 1;
end

lon_plot = [lon_native; lon_native(1) + direction.*2*pi];
data_plot = [data_native data_native(:,1)];
end
