function sp = music_1d(R, n, design, wavelength, grid_size, varargin)
%MUSIC_1D 1D MUSIC.
%Syntax:
%   sp = MUSIC_1D(R, n, design, wavelength, grid_size, ...);
%   sp = MUSIC_1D(R, n, f_steering, [], grid_size, ...);
%Inputs:
%   R - Sample covariance matrix.
%   n - Number of sources.
%   design - Array design. Can also be a function handle that generates
%            a steering matrix. This function must take two arguments,
%            wavelength and the doa vector.
%   wavelength - Wavelength.
%   grid_size - Number of grid points used.
%   ... - Options:
%           'Unit' - Can be 'radian', 'degree', or 'sin'. Default value is
%                   'radian'.
%           'RefineEstimates' - If set to true, will refine the estimated
%                               direction of arrivals around the grid.
%Output:
%   sp - Spectrum structure with the following fields:
%           x - An 1 x grid_size vector.
%           y - An 1 x grid_size vector. Calling `plot(x, y)` will plot the
%               spectrum.
%           x_est - An 1 x n vector storing the estimated DOAs. May not
%                   fall on the grid if 'RefineEstimates' is set to true.
%           x_unit - The same as the unit specified by 'Unit'.
%           resolved - True if the number of peaks in the spectrum is
%                      greater or equal to the number of sources.
%           discrete - Constant value false.
%Reference:
%   [1] R. Schmidt, "Multiple emitter location and signal parameter
%       estimation," IEEE Transactions on Antennas and Propagation,
%       vol. 34, no. 3, pp. 276-280, Mar. 1986.
unit = 'radian';
refine_estimates = false;
for ii = 1:2:nargin-5
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'unit'
            unit = lower(option_value);
        case 'refineestimates'
            refine_estimates = true;
        otherwise
            error('Unknown option ''%s''.', option_name);
    end
end
m = size(R, 1);
if n >= m
    error('Too many sources.');
end
% discretize and create the corresponding steering matrix
[doa_grid_rad, doa_grid_display, ~] = default_doa_grid(grid_size, unit, 1);
% find noise subspace
[U, D] = eig(0.5*(R + R'), 'vector');
% possible asymmetry due to floating point error
if ~isreal(D)
    eig_values = abs(D);
    [~, I] = sort(eig_values);
    Un = U(:, I(1:end-n));
else
    Un = U(:, 1:end-n);
end
% compute spectrum
sp_intl = 1./compute_inv_spectrum(Un, design, wavelength, doa_grid_rad);
[x_est, x_est_idx, resolved] = find_doa_from_spectrum_1d(doa_grid_display, sp_intl, n);
% refine
if resolved && refine_estimates
    switch unit
        case 'radian'
            f_obj = @(x) compute_inv_spectrum(Un, design, wavelength, x);
        case 'degree'
            f_obj = @(x) compute_inv_spectrum(Un, design, wavelength, deg2rad(x));
        case 'sin'
            f_obj = @(x) compute_inv_spectrum(Un, design, wavelength, asin(x));
        otherwise
            error('Invalid unit ''%s''.', unit);
    end
    x_est = refine_grid_estimates(f_obj, doa_grid_display, x_est_idx);
end
% return
sp = struct();
sp.x = doa_grid_display;
sp.x_est = x_est;
sp.x_unit = unit;
sp.y = sp_intl;
sp.resolved = resolved;
sp.discrete = false;
end

function v = compute_inv_spectrum(Un, design, wavelength, theta)
if ishandle(design)
    A = design(wavelength, theta);
else
    A = steering_matrix(design, wavelength, theta);
end
v = Un' * A;
v = real(sum(conj(v) .* v, 1));
end

