function [result_focal_normalized, result_random_big_normalized, result_random_small_normalized] = spread(sample_name, filename_stack_1, filename_stack_2, filename_stack_3, filename_mask, filename_veins, damage_r, damage_c, damage_radius, focal_distance_max, n_samples_random)
    % damage_* are in original coordinates on vein image

    SIZE_BIG = 15;
    SIZE_SMALL = 5;
    
    mkdir('results');
    mkdir('temp');

    focal_distance_step = 10; % not really analyzed
    
    fprintf('load in stack 1\n');
    fn_1 = unzip(filename_stack_1,'temp');
    stack_1 = tiffreadVolume(fn_1{1});
    delete(fn_1{1});

    if (~strcmp(filename_stack_2,''))
        fprintf('load in stack 2\n');
        fn_2 = unzip(filename_stack_2,'temp');
        stack_2 = tiffreadVolume(fn_2{1});
        delete(fn_2{1});
    end

    if (~strcmp(filename_stack_3,''))
        fprintf('load in stack 3\n');
        fn_3 = unzip(filename_stack_3,'temp');
        stack_3 = tiffreadVolume(fn_3{1});
        delete(fn_3{1});
    end

    % remove temp files
    try
        rmdir('temp');
    end
    
    fprintf('mask the stacks\n');
    mask = logical(imread(filename_mask));
    
    fprintf('concatenate the stacks\n');
    stack = stack_1;
    if (~strcmp(filename_stack_2,''))
        stack = logical(cat(3, stack, stack_2));
        clear stack_2;

        if (~strcmp(filename_stack_3,''))
            stack = logical(cat(3, stack, stack_3));
            clear stack_3;
        end
    end
    clear stack_1;
    
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

    % make mask
    fprintf('making border mask\n');
    mask_border = true(size(mask_scaled));
    mask_border(1:focal_distance_max,:) = 0;
    mask_border(:,1:focal_distance_max) = 0;
    mask_border = mask_border .* flip(mask_border,1) .* flip(mask_border,2);

    % identify big veins
    fprintf('finding big and small veins\n');
    veins_scaled(mask_scaled~=0) = 0;
    veins_scaled_binary = imbinarize(veins_scaled, graythresh(veins_scaled));
    veins_dist = bwdist(veins_scaled_binary==0);

    veins_dist_transformed = veins_dist;
    veins_dist_transformed(veins_scaled_binary==0) = 0;

    veins_big = (veins_dist_transformed > SIZE_BIG);
    veins_small = (veins_dist_transformed < SIZE_SMALL & veins_dist_transformed > 2);
    
    fprintf('pick random coordinates\n');
    % pick random coordinates within the buffer area at a certain size class 
    % (so any focal regions do not hit image boundary)
    
    [idx_r_big,idx_c_big] = find(mask_scaled==0 & mask_border==1 & veins_big==1);
    ordering_big = randperm(length(idx_r_big),n_samples_random);
    idx_r_big_ss = idx_r_big(ordering_big);
    idx_c_big_ss = idx_c_big(ordering_big);

    [idx_r_small,idx_c_small] = find(mask_scaled==0 & mask_border==1 & veins_small==1);
    ordering_small = randperm(length(idx_r_small),n_samples_random);
    idx_r_small_ss = idx_r_small(ordering_small);
    idx_c_small_ss = idx_c_small(ordering_small);
    
    fprintf('show damage and random locations\n');
    f1 = figure;
    img = ind2rgb(8*uint8(stack_max) + 2*uint8(mask_scaled) + 4*uint8(mask_border),gray(32));
    img = insertShape(img,'filled-rectangle',[idx_c_big_ss - focal_distance_max, idx_r_big_ss - focal_distance_max, repmat(2*focal_distance_max, [n_samples_random 1]) repmat(2*focal_distance_max, [n_samples_random 1])],ShapeColor="blue",Opacity=0.2);
    img = insertShape(img,'filled-rectangle',[idx_c_small_ss - focal_distance_max, idx_r_small_ss - focal_distance_max, repmat(2*focal_distance_max, [n_samples_random 1]) repmat(2*focal_distance_max, [n_samples_random 1])],ShapeColor="cyan",Opacity=0.2);
    if (~isnan(damage_radius))
        img = insertShape(img,'filledcircle',[damage_c/scale_factor_veins damage_r/scale_factor_veins, damage_radius/scale_factor_veins],ShapeColor="black",Opacity=0.6);
        img = insertShape(img,'filled-rectangle',[damage_c/scale_factor_veins - focal_distance_max, damage_r/scale_factor_veins - focal_distance_max, 2*focal_distance_max, 2*focal_distance_max],ShapeColor="red",Opacity=0.6);
    end
    
    imshow(img);
    saveas(f1, sprintf('results/%s_image',sample_name),'png');
    close(f1);
    





    fprintf('extract random damage location data\n');

    result_random_big = NaN([size(stack_masked,3) n_samples_random]);
    for i=1:n_samples_random
        focal_r = idx_r_big_ss(i);
        focal_c = idx_c_big_ss(i);
    
        fprintf('%d/%d x=%d y=%d\n',i, n_samples_random, focal_r, focal_c);
    
        stack_subset = stack_masked((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max), :);
        values_focal_time_series = sum(stack_subset, 1:2);
    
        result_random_big(:,i) = values_focal_time_series;
    end
    % normalize by the area being considered
    result_random_big_normalized = result_random_big ./ (2*focal_distance_max)^2;

    result_random_small = NaN([size(stack_masked,3) n_samples_random]);
    for i=1:n_samples_random
        focal_r = idx_r_small_ss(i);
        focal_c = idx_c_small_ss(i);
    
        fprintf('%d/%d x=%d y=%d\n',i, n_samples_random, focal_r, focal_c);
    
        stack_subset = stack_masked((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max), :);
        values_focal_time_series = sum(stack_subset, 1:2);
    
        result_random_small(:,i) = values_focal_time_series;
    end
    % normalize by the area being considered
    result_random_small_normalized = result_random_small ./ (2*focal_distance_max)^2;
    
    fprintf('extract focal damage location data\n');
    focal_distance_series = focal_distance_step:focal_distance_step:focal_distance_max;
    result_focal = NaN([size(stack_masked,3) length(focal_distance_series)]);
    if (~isnan(damage_radius))
        % try at different distance scales
        for i=1:length(focal_distance_series)
            focal_r = floor(damage_r/scale_factor_veins);
            focal_c = floor(damage_c/scale_factor_veins);
        
            fprintf('%d/%d r=%d c=%d distance=%d\n',i, length(focal_distance_series), focal_r, focal_c, focal_distance_series(i));
        
            stack_subset = stack_masked((focal_r-focal_distance_series(i)):(focal_r+focal_distance_series(i)), (focal_c-focal_distance_series(i)):(focal_c+focal_distance_series(i)), :);
            values_focal_time_series = sum(stack_subset, 1:2);
        
            result_focal(:,i) = values_focal_time_series;
        end
    end
    result_focal_normalized = result_focal ./ (2*repmat(focal_distance_series, [size(stack_masked, 3) 1])).^2;
    
    fprintf('show observed data normalized\n');
    f2 = figure;
    plot(cumsum(result_random_big_normalized),'-b');
    hold on;
    plot(cumsum(result_random_small_normalized),'-c');
    % plot mean of observed data
    plot(mean(cumsum(result_random_big_normalized),2),'-b','LineWidth',3);
    plot(mean(cumsum(result_random_small_normalized),2),'-c','LineWidth',3);
    % plot focal damage data
    plot(cumsum(result_focal_normalized),'-r','LineWidth',3);
    saveas(f2, sprintf('results/%s_time_series',sample_name),'png');
    close(f2);

    fprintf('write out files\n');
    writematrix(result_focal_normalized,sprintf('results/%s_focal_normalized.csv',sample_name));
    writematrix(result_random_big_normalized,sprintf('results/%s_random_big_normalized.csv',sample_name));
    writematrix(result_random_small_normalized,sprintf('results/%s_random_small_normalized.csv',sample_name));
end