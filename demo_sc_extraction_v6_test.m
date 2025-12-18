%% Demo: SAR 散射中心提取测试 (ASC-NOMP V6)
% 测试改进后的算法：
% 1. 通过 L 判断散射中心类型（L=0 局域式，L>0 分布式）
% 2. 允许多个局域式散射中心
% 3. 分布式 gamma=0，局域式 gamma>0
% 4. 改进的参数估计顺序

function demo_sc_extraction_v6_test()
    fprintf('=== SAR Scattering Center Extraction Test (V6) ===\n\n');
    
    % 生成测试图像
    [Img, ground_truth] = generate_test_image();
    
    % 配置参数
    user_cfg = struct();
    user_cfg.verbose = true;
    user_cfg.max_atoms = 50;
    user_cfg.max_localized = Inf;  % 允许多个局域式散射中心
    user_cfg.stop_explained_ratio = 0.85;
    
    % 运行提取算法
    fprintf('\n--- Running ASC-NOMP extraction ---\n');
    [theta_list, residual_norms, info] = main_ascm_nomp_v6(Img, user_cfg);
    
    % 显示结果
    fprintf('\n--- Extraction Results ---\n');
    fprintf('Total atoms extracted: %d\n', numel(theta_list));
    fprintf('Energy explained: %.2f%%\n', 100*info.explained_ratio);
    
    % 统计分析
    print_statistics(theta_list);
    
    % 验收测试
    run_validation_tests(theta_list);
    
    % 可视化
    visualize_results(Img, theta_list, info);
    
    fprintf('\n=== Test Complete ===\n');
end

%% 生成测试图像
function [Img, ground_truth] = generate_test_image()
    % 创建一个 64x64 的测试图像，包含分布式和局域式散射中心
    H = 64; W = 64;
    Img = zeros(H, W);
    
    % Ground truth
    ground_truth = struct('localized', [], 'distributed', []);
    
    % 添加分布式散射中心（线状目标）
    % 斜线 1
    for i = 1:20
        r = 20 + i;
        c = 15 + i;
        if r <= H && c <= W
            Img(r, c) = 0.8 * exp(-i*0.05);
        end
    end
    ground_truth.distributed(1).type = 'line';
    ground_truth.distributed(1).center = [30, 25];
    
    % 斜线 2
    for i = 1:15
        r = 40 - i;
        c = 40 + i;
        if r > 0 && r <= H && c <= W
            Img(r, c) = 0.7 * exp(-i*0.05);
        end
    end
    ground_truth.distributed(2).type = 'line';
    ground_truth.distributed(2).center = [33, 47];
    
    % 添加局域式散射中心（点状目标）
    % 点 1
    r1 = 25; c1 = 50;
    for dr = -3:3
        for dc = -3:3
            r = r1 + dr; c = c1 + dc;
            if r > 0 && r <= H && c > 0 && c <= W
                dist = sqrt(dr^2 + dc^2);
                Img(r, c) = Img(r, c) + 0.9 * exp(-dist * 0.5);
            end
        end
    end
    ground_truth.localized(1).center = [r1, c1];
    
    % 点 2
    r2 = 45; c2 = 20;
    for dr = -2:2
        for dc = -2:2
            r = r2 + dr; c = c2 + dc;
            if r > 0 && r <= H && c > 0 && c <= W
                dist = sqrt(dr^2 + dc^2);
                Img(r, c) = Img(r, c) + 0.85 * exp(-dist * 0.6);
            end
        end
    end
    ground_truth.localized(2).center = [r2, c2];
    
    % 添加噪声
    Img = Img + 0.05 * randn(H, W);
    Img(Img < 0) = 0;
    
    fprintf('Generated test image: %dx%d\n', H, W);
    fprintf('  - %d distributed scattering centers (lines)\n', numel(ground_truth.distributed));
    fprintf('  - %d localized scattering centers (points)\n', numel(ground_truth.localized));
end

%% 统计分析
function print_statistics(theta_list)
    if isempty(theta_list)
        fprintf('\nNo atoms extracted.\n');
        return;
    end
    
    fprintf('\n--- Parameter Statistics ---\n');
    
    % 分类统计
    is_localized = [theta_list.L] == 0;
    is_distributed = [theta_list.L] > 0;
    
    n_loc = sum(is_localized);
    n_dist = sum(is_distributed);
    
    fprintf('\nScattering center types:\n');
    fprintf('  Localized (L=0):     %d atoms\n', n_loc);
    fprintf('  Distributed (L>0):   %d atoms\n', n_dist);
    
    % 局域式散射中心参数
    if n_loc > 0
        loc_atoms = theta_list(is_localized);
        gammas = [loc_atoms.gamma];
        alphas = [loc_atoms.alpha];
        
        fprintf('\nLocalized scattering centers:\n');
        fprintf('  Gamma range: [%.3f, %.3f]\n', min(gammas), max(gammas));
        fprintf('  Gamma mean:  %.3f\n', mean(gammas));
        fprintf('  Alpha range: [%.3f, %.3f]\n', min(alphas), max(alphas));
        fprintf('  Alpha mean:  %.3f\n', mean(alphas));
    end
    
    % 分布式散射中心参数
    if n_dist > 0
        dist_atoms = theta_list(is_distributed);
        Ls = [dist_atoms.L];
        gammas = [dist_atoms.gamma];
        alphas = [dist_atoms.alpha];
        
        fprintf('\nDistributed scattering centers:\n');
        fprintf('  L range:     [%.3f, %.3f]\n', min(Ls), max(Ls));
        fprintf('  L mean:      %.3f\n', mean(Ls));
        fprintf('  Gamma range: [%.3f, %.3f]\n', min(gammas), max(gammas));
        fprintf('  Gamma mean:  %.3f (should be 0)\n', mean(gammas));
        fprintf('  Alpha range: [%.3f, %.3f]\n', min(alphas), max(alphas));
        fprintf('  Alpha mean:  %.3f\n', mean(alphas));
    end
end

%% 验收测试
function run_validation_tests(theta_list)
    fprintf('\n--- Validation Tests ---\n');
    
    if isempty(theta_list)
        fprintf('FAIL: No atoms extracted\n');
        return;
    end
    
    n_pass = 0;
    n_total = 0;
    
    % Test 1: 分布式原子 L > 0
    n_total = n_total + 1;
    dist_atoms = theta_list([theta_list.L] > 0);
    if ~isempty(dist_atoms)
        Ls = [dist_atoms.L];
        if all(Ls > 0)
            fprintf('PASS: Test 1 - Distributed atoms have L > 0\n');
            n_pass = n_pass + 1;
        else
            fprintf('FAIL: Test 1 - Some distributed atoms have L <= 0\n');
        end
    else
        fprintf('SKIP: Test 1 - No distributed atoms found\n');
    end
    
    % Test 2: 分布式原子 gamma = 0
    n_total = n_total + 1;
    if ~isempty(dist_atoms)
        gammas = [dist_atoms.gamma];
        if all(abs(gammas) < 1e-6)
            fprintf('PASS: Test 2 - Distributed atoms have gamma = 0\n');
            n_pass = n_pass + 1;
        else
            fprintf('FAIL: Test 2 - Some distributed atoms have gamma != 0 (max=%.3f)\n', max(abs(gammas)));
        end
    else
        fprintf('SKIP: Test 2 - No distributed atoms found\n');
    end
    
    % Test 3: 局域式原子 gamma > 0
    n_total = n_total + 1;
    loc_atoms = theta_list([theta_list.L] == 0);
    if ~isempty(loc_atoms)
        gammas = [loc_atoms.gamma];
        if all(gammas > 0)
            fprintf('PASS: Test 3 - Localized atoms have gamma > 0\n');
            n_pass = n_pass + 1;
        else
            fprintf('FAIL: Test 3 - Some localized atoms have gamma <= 0\n');
        end
    else
        fprintf('SKIP: Test 3 - No localized atoms found\n');
    end
    
    % Test 4: Alpha 值有效
    n_total = n_total + 1;
    alphas = [theta_list.alpha];
    if all(alphas > 0) && all(isfinite(alphas))
        fprintf('PASS: Test 4 - All alpha values are valid (> 0)\n');
        n_pass = n_pass + 1;
    else
        fprintf('FAIL: Test 4 - Some alpha values are invalid\n');
    end
    
    % Test 5: 成功提取原子
    n_total = n_total + 1;
    if numel(theta_list) > 0
        fprintf('PASS: Test 5 - Successfully extracted %d atoms\n', numel(theta_list));
        n_pass = n_pass + 1;
    else
        fprintf('FAIL: Test 5 - No atoms extracted\n');
    end
    
    fprintf('\nTest Summary: %d/%d tests passed (%.1f%%)\n', n_pass, n_total, 100*n_pass/n_total);
end

%% 可视化结果
function visualize_results(Img, theta_list, info)
    fprintf('\n--- Generating Visualization ---\n');
    
    figure('Position', [100, 100, 1200, 400]);
    
    % 子图 1: 原始图像
    subplot(1, 3, 1);
    imagesc(Img);
    colormap('gray');
    colorbar;
    title('Original Image');
    axis image;
    
    % 子图 2: 提取的散射中心（按类型标记）
    subplot(1, 3, 2);
    imagesc(Img);
    colormap('gray');
    hold on;
    
    if ~isempty(theta_list)
        % 分离局域式和分布式
        is_localized = [theta_list.L] == 0;
        is_distributed = [theta_list.L] > 0;
        
        % 绘制局域式（红色圆圈）
        if any(is_localized)
            loc_atoms = theta_list(is_localized);
            rows = [loc_atoms.row];
            cols = [loc_atoms.col];
            plot(cols, rows, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        end
        
        % 绘制分布式（蓝色方块）
        if any(is_distributed)
            dist_atoms = theta_list(is_distributed);
            rows = [dist_atoms.row];
            cols = [dist_atoms.col];
            plot(cols, rows, 'bs', 'MarkerSize', 8, 'LineWidth', 2);
            
            % 绘制方向
            for i = 1:numel(dist_atoms)
                phi = dist_atoms(i).phi;
                L = dist_atoms(i).L;
                len = min(L, 10);
                dx = len * cos(phi);
                dy = -len * sin(phi);
                quiver(cols(i), rows(i), dx, dy, 0, 'b', 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
            end
        end
    end
    
    hold off;
    title(sprintf('Extracted Centers (%d total)', numel(theta_list)));
    legend({'Localized (L=0)', 'Distributed (L>0)'}, 'Location', 'best');
    axis image;
    
    % 子图 3: ROI 和统计
    subplot(1, 3, 3);
    imagesc(info.roi_mask);
    colormap(gca, 'parula');
    title(sprintf('ROI (Explained: %.1f%%)', 100*info.explained_ratio));
    axis image;
    
    fprintf('Visualization complete.\n');
end
