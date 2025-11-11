function [] = spread(sample_name, filename_stack_1, filename_stack_2, filename_stack_3, filename_mask, filename_veins, damage_r, damage_c, damage_radius, focal_distance_max, n_samples_random)
    % damage_* are in original coordinates on vein image

    SIZE_BIG = 12;
    SIZE_SMALL = 4;
    
    mkdir('results');
    mkdir('temp');

    % load in the mask
    mask = logical(imread(filename_mask));

    % make a mask image with the damage cutout and then inverted
    fprintf('make damage mask\n');
    if (~isnan(damage_c))
        mask_with_damage = insertShape(uint8(mask), 'filledcircle',[damage_c damage_r, damage_radius],ShapeColor='white');
        mask_with_damage = ~logical(mask_with_damage(:,:,1)>0);
    else
        mask_with_damage = ~mask;
    end
    f000=figure;
    imshow(mask_with_damage);
    saveas(f000, sprintf('results/%s_%d_mask_with_damage',sample_name, focal_distance_max),'png');
    close(f000);
    
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
    catch
        fprintf('*** temp dir could not be removed\n')
    end
    


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
    
    % get the mask and use it
    fprintf('mask the stacks\n');
    mask_scaled = imresize(mask, size(stack,[1 2]));
    
    mask_array = repmat(mask_scaled, [1 1 size(stack, 3)]);
    
    stack_masked = stack & ~mask_array;

    % also scale the damage mask
    mask_with_damage_scaled = imresize(mask_with_damage, size(stack,[1 2]));
    
    % get the embolism time series
    fprintf('get time series of all embolisms across leaf\n');
    ts = squeeze(sum(stack_masked, [1 2]));
    
    f00 = figure;
    plot(ts); hold on; plot(cumsum(ts));
    saveas(f00, sprintf('results/%s_%d_ts',sample_name, focal_distance_max),'png');
    close(f00);

    fprintf('summing all embolisms\n');
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
        img = insertShape(img,'filledcircle',[damage_c/scale_factor_veins damage_r/scale_factor_veins, damage_radius/scale_factor_veins],ShapeColor="red");
    end
    imshow(img);
    saveas(f0, sprintf('results/%s_%d_damage',sample_name, focal_distance_max),'png');
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


    % get total lamina pixels and total vein pixels
    px_count_mask = sum(mask_scaled(:) == 0);
    px_count_veins = sum(veins_scaled_binary(:) ~= 0);
    px_count_total = prod(size(mask_scaled));
    writematrix([px_count_mask, px_count_veins, px_count_total], sprintf('results/%s_%d_px_count.csv',sample_name, focal_distance_max));


    
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
    saveas(f1, sprintf('results/%s_%d_image',sample_name, focal_distance_max),'png');
    close(f1);
    





    fprintf('extract random damage location data\n');
    result_random_big = NaN([size(stack_masked,3) n_samples_random]);
    result_random_big_mask = NaN(1, n_samples_random);
    for i=1:n_samples_random
        focal_r = idx_r_big_ss(i);
        focal_c = idx_c_big_ss(i);
    
        stack_subset = stack_masked((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max), :);
        values_focal_time_series = sum(stack_subset, 1:2);
    
        result_random_big(:,i) = values_focal_time_series;

        values_mask_with_damage_big = mask_with_damage_scaled((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max));

        result_random_big_mask(:,i) = sum(values_mask_with_damage_big(:));

        fprintf('big %d/%d x=%d y=%d %d\n', i, n_samples_random, focal_r, focal_c, result_random_big_mask(:,i));
    end

    result_random_small = NaN([size(stack_masked,3) n_samples_random]);
    result_random_small_mask = NaN(1, n_samples_random);
    for i=1:n_samples_random
        focal_r = idx_r_small_ss(i);
        focal_c = idx_c_small_ss(i);
    
        stack_subset = stack_masked((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max), :);
        values_focal_time_series = sum(stack_subset, 1:2);
    
        result_random_small(:,i) = values_focal_time_series;

        values_mask_with_damage_small = mask_with_damage_scaled((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max));

        result_random_small_mask(:,i) = sum(values_mask_with_damage_small(:));

        fprintf('small %d/%d x=%d y=%d %d\n',i, n_samples_random, focal_r, focal_c, result_random_small_mask(:,i));
    end
    
    fprintf('extract focal damage location data\n');
    result_focal = NaN([size(stack_masked,3) 1]);
    result_focal_mask = NaN;

    if (~isnan(damage_radius))
        focal_r = floor(damage_r/scale_factor_veins);
        focal_c = floor(damage_c/scale_factor_veins);

        stack_subset = stack_masked((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max), :);
        values_focal_time_series = sum(stack_subset, 1:2);
    
        result_focal = squeeze(values_focal_time_series);

        values_mask_with_damage_focal = mask_with_damage_scaled((focal_r-focal_distance_max):(focal_r+focal_distance_max), (focal_c-focal_distance_max):(focal_c+focal_distance_max));
        result_focal_mask = sum(values_mask_with_damage_focal(:));

        fprintf('focal r=%d c=%d %d\n', focal_r, focal_c, result_focal_mask);
    end
    
    fprintf('show observed data normalized\n');
    f2 = figure;
    plot(cumsum(result_random_big),'-b');
    hold on;
    plot(cumsum(result_random_small),'-c');
    % plot mean of observed data
    plot(mean(cumsum(result_random_big),2),'-b','LineWidth',3);
    plot(mean(cumsum(result_random_small),2),'-c','LineWidth',3);
    % plot focal damage data
    
    plot(cumsum(result_focal,1),'-r','LineWidth',3);
    saveas(f2, sprintf('results/%s_%d_time_series',sample_name, focal_distance_max),'png');
    close(f2);

    fprintf('write out files\n');
    writematrix(ts,sprintf('results/%s_%d_ts.csv',sample_name, focal_distance_max));
    writematrix([result_focal_mask; result_focal],sprintf('results/%s_%d_focal.csv',sample_name, focal_distance_max));
    size(result_random_big)
    size(result_random_big_mask)
    writematrix([squeeze(result_random_big_mask); result_random_big ],sprintf('results/%s_%d_random_big.csv',sample_name, focal_distance_max));
    writematrix([squeeze(result_random_small_mask); result_random_small ],sprintf('results/%s_%d_random_small.csv',sample_name, focal_distance_max));
end