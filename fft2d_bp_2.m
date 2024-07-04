
function [data_f, data_fft2, filter] = fft2d_bp_2(data_0, cutoff_x_min, cutoff_x_max, ...
    cutoff_y_min, cutoff_y_max, w_min, w_max)
    % include both distance based methods and gradual change at the boundary in
    % this method
    
    % data_0: original image
    % cutoff_x_min: min x-cutoff, # of frequency points
    % cutoff_x_min: max x-cutoff, # of frequency points
    % cutoff_y_min: min y-cutoff, # of frequency points
    % cutoff_y_min: max y-cutoff, # of frequency points
    % w_min: factor of min gradual change, range: [0, 1), outside this
    % range will be ignored by the function
    % w_max: factor of min gradual change, range: (1, +infty], outside
    % ignored

    % output:
    % data_f: data after filtering. Should be a real number.
    % data_fft2: data after 2D fft and shift the zero frequency to center
    % filter: frequency domain filter

    % check the dimension before performing the filtering, zero components
    % at the center
    data_fft2 = fftshift(fft2(data_0));

    % Dimensions of the data
    [rows, cols] = size(data_0);
    
    % Define the filter
    one_arr = ones([rows, cols]);
    filter = zeros([rows, cols]);
    mask_grad_min = zeros([rows, cols]);
    mask_grad_max = zeros([rows, cols]);

    % Distance matrix
    % Even and odd cases are different. 
    % When having odd number of points, define the zero frequency as the 
    % center;
    % When even number of points, define the first negative frequency 
    % as artificial frequency, zero frequency is at num/2.
    if mod(cols, 2) == 0
        x_arr = linspace(-cols/2, cols/2 - 1, cols);
    else
        x_arr = linspace(-(cols-1)/2, (cols-1)/2, cols);
    end

    if mod(rows, 2) == 0
        y_arr = linspace(-rows/2, rows/2 - 1, rows);
    else
        y_arr = linspace(-(rows-1)/2, (rows-1)/2, rows);
    end

    
    [xv, yv] = meshgrid(x_arr, y_arr);

    % consider the case for low pass
    if cutoff_x_min == 0 || cutoff_y_min == 0
        elips_2 = (xv/cutoff_x_max).^2 + (yv/cutoff_y_max).^2;
        elips_4 = (xv/cutoff_x_max/w_max).^2 + (yv/cutoff_y_max/w_max).^2;
        filter(elips_2 <= one_arr) = 1;
        if w_max > 1
            mask_grad_max(elips_4 <= one_arr & elips_2 > one_arr) = 1;
            gradual_region_max = (sqrt(elips_2) - w_max)/(1 - w_max);
            filter = filter + mask_grad_max .* gradual_region_max;
        end

    % high pass
    elseif cutoff_x_max == 0 || cutoff_y_max == 0
        elips_1 = (xv/cutoff_x_min).^2 + (yv/cutoff_y_min).^2;
        elips_3 = (xv/cutoff_x_min/w_min).^2 + (yv/cutoff_y_min/w_min).^2;
        filter(elips_1 >= one_arr) = 1;
        if w_min < 1 && w_min >=0
            mask_grad_min(elips_3 >= one_arr & elips_1 < one_arr) = 1;
            gradual_region_min = (sqrt(elips_1) - w_min)/(1 - w_min);
            filter = filter + mask_grad_min .* gradual_region_min;
        end

    % band pass
    else
        % The elipses of cutoff
        elips_1 = (xv/cutoff_x_min).^2 + (yv/cutoff_y_min).^2;
        elips_2 = (xv/cutoff_x_max).^2 + (yv/cutoff_y_max).^2;
    
        % add non-zero components
        filter(elips_1 >= one_arr & elips_2 <= one_arr) = 1;
        
        % The elipses of gradual change
        if w_min < 1 && w_min >=0
            elips_3 = (xv/cutoff_x_min/w_min).^2 + (yv/cutoff_y_min/w_min).^2;
            mask_grad_min(elips_1 < one_arr & elips_3 >= one_arr) = 1;
            % Compute the expression with a linear gradual change
            gradual_region_min = (sqrt(elips_1) - w_min)/(1 - w_min);
            filter = filter + mask_grad_min .* gradual_region_min;
        end

        if w_max > 1
            elips_4 = (xv/cutoff_x_max/w_max).^2 + (yv/cutoff_y_max/w_max).^2;
            mask_grad_max(elips_4 <= one_arr & elips_2 > one_arr) = 1;
            gradual_region_max = (sqrt(elips_2) - w_max)/(1 - w_max);
            filter = filter + mask_grad_max .* gradual_region_max;
        end
    
    end
    
    data_fft2_f = data_fft2 .* filter;

    data_f=(ifft2(ifftshift(data_fft2_f)));

end