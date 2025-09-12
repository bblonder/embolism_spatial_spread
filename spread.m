function [result_focal_normalized, result_random_normalized] = spread(sample_name, filename_stack_1, filename_stack_2, filename_mask, filename_veins, damage_r, damage_c, damage_radius, focal_distance_max, n_samples_random)
    % damage_* are in original coordinates on vein image
    
    mkdir('results');

    focal_distance_step = 10; % not really analyzed
    
    fprintf('load in stacks\n');
    stack_1 = tiffreadVolume(filename_stack_1);

    if (~strcmp(filename_stack_2,''))
        stack_2 = tiffreadVolume(filename_stack_2);
    end
    
    fprintf('mask the stacks\n');
    mask = logical(imread(filename_mask));
    
    if (~strcmp(filename_stack_2,''))
        stack = logical(cat(3, stack_1, stack_2));
        clear stack_1;
        clear stack_2;
    else
        stack = stack_1;
        clear stack_1;
    end
    
    mask_scaled = imresize(mask, size(stack,[1 2]));
    
    mask_array = repmat(mask_scaled, [1 1 size(stack, 3)]);
    
    stack_masked = stack & ~mask_array;
    
    fprintf('get all embolisms\n');
    stack_max = sum(stack_masked, 3);
    
    fprintf('resize veins\n');
    veins = imread(filename_veins);
    veins_scaled = imresize(veins, size(stack,[1 2]));
    
    % estimate scale factor
    scale_factor_veins = mean(size(veins) ./ size(stack,[1 2]));
    %scale_factor_mask = mean(size(mask) ./ size(stack,[1 2]));
    
    fprintf('show damage and veins\n');
    f0 = figure;
    img = veins_scaled;
    if (~isnan(damage_radius))
        img = insertShape(img,'filledcircle',[damage_c/scale_factor_veins damage_r/scale_factor_veins, damage_radius/scale_factor_veins],ShapeColor="blue");
    end
    imshow(img);
    saveas(f0, sprintf('results/%s_damage',sample_name),'png');
    close(f0);
    
    % save memory
    clear stack;
    clear mask_array;
    
    fprintf('pick random coordinates\n');
    % pick random coordinates within the buffer area (so any focal regions do
    % not hit image boundary)
    mask_border = true(size(mask_scaled));
    mask_border(1:focal_distance_max,:) = 0;
    mask_border(:,1:focal_distance_max) = 0;
    mask_border = mask_border .* flip(mask_border,1) .* flip(mask_border,2);
    
    [idx_r,idx_c] = find(mask_scaled==0 & mask_border==1);
    ordering = randperm(length(idx_r),n_samples_random);
    idx_r_ss = idx_r(ordering);
    idx_c_ss = idx_c(ordering);
    
    fprintf('show damage and random locations\n');
    f1 = figure;
    img = ind2rgb(8*uint8(stack_max) + 2*uint8(mask_scaled) + 4*uint8(mask_border),hsv(64));
    img = insertShape(img,'filled-rectangle',[idx_c_ss - focal_distance_max, idx_r_ss - focal_distance_max, repmat(2*focal_distance_max, [n_samples_random 1]) repmat(2*focal_distance_max, [n_samples_random 1])],ShapeColor="black",Opacity=0.2);
    if (~isnan(damage_radius))
        img = insertShape(img,'filledcircle',[damage_c/scale_factor_veins damage_r/scale_factor_veins, damage_radius/scale_factor_veins],ShapeColor="blue",Opacity=0.6);
        img = insertShape(img,'filled-rectangle',[damage_c/scale_factor_veins - focal_distance_max, damage_r/scale_factor_veins - focal_distance_max, 2*focal_distance_max, 2*focal_distance_max],ShapeColor="magenta",Opacity=0.6);
    end
    
    imshow(img);
    saveas(f1, sprintf('results/%s_image',sample_name),'png');
    close(f1);
    
    fprintf('extract random damage location data\n');

    result_random = NaN([size(stack_masked,3) n_samples_random]);
    for i=1:n_samples_random
        focal_r = idx_r_ss(i);
        focal_c = idx_c_ss(i);
    
        fprintf('%d/%d x=%d y=%d\n',i, n_samples_random, focal_r, focal_c);
    
        stack_subset = stack_masked((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max), :);
        values_focal_time_series = sum(stack_subset, 1:2);
    
        result_random(:,i) = values_focal_time_series;
    end
    % normalize by the area being considered
    result_random_normalized = result_random ./ (2*focal_distance_max)^2;
    
    fprintf('extract focal damage location data\n');
    focal_distance_series = focal_distance_step:focal_distance_step:focal_distance_max;
    result_focal = NaN([size(stack_masked,3) length(focal_distance_series)]);
    if (~isnan(damage_radius))
        % try at different distance scales
        for i=1:length(focal_distance_series)
            focal_r = damage_r/scale_factor_veins;
            focal_c = damage_c/scale_factor_veins;
        
            fprintf('%d/%d r=%d c=%d distance=%d\n',i, length(focal_distance_series), focal_r, focal_c, focal_distance_series(i));
        
            stack_subset = stack_masked((focal_r-focal_distance_series(i)):(focal_r+focal_distance_series(i)), (focal_c-focal_distance_series(i)):(focal_c+focal_distance_series(i)), :);
            values_focal_time_series = sum(stack_subset, 1:2);
        
            result_focal(:,i) = values_focal_time_series;
        end
    end
    result_focal_normalized = result_focal ./ (2*repmat(focal_distance_series, [size(stack_masked, 3) 1])).^2;
    
    fprintf('show observed data normalized\n');
    f2 = figure;
    plot(cumsum(result_random_normalized),'-k');
    hold on;
    % plot mean of observed data
    plot(mean(cumsum(result_random_normalized),2),'-r','LineWidth',3);
    % plot focal damage data
    plot(cumsum(result_focal_normalized),'-b','LineWidth',3)
    saveas(f2, sprintf('results/%s_time_series',sample_name),'png');
    close(f2);

    fprintf('write out files\n');
    writematrix(result_focal_normalized,sprintf('results/%s_focal_normalized.csv',sample_name));
    writematrix(result_random_normalized,sprintf('results/%s_random_normalized.csv',sample_name));
end