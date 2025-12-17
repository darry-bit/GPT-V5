%% demo_sc_extraction_v6_test.m
% 测试 ASC 提取的 L、alpha、gamma 参数估计功能
%
% 此脚本演示：
% 1. 创建合成散射体图像（包含不同长度的散射体）
% 2. 使用增强的 main_ascm_nomp_v6 提取 ASC 参数
% 3. 显示和验证提取的 L、alpha、gamma 参数
% 4. 可视化结果

clear; close all; clc;

fprintf('=== ASC Parameter Estimation Test ===\n\n');

%% 1. 创建测试图像
% 由于没有实际 MSTAR 数据，我们创建一个合成图像进行测试
fprintf('1. Creating synthetic test image...\n');

im_size = 128;
Img = zeros(im_size, im_size);

% 添加背景噪声
Img = Img + 0.05 * randn(im_size, im_size);

% 创建几个不同长度的散射体
% 散射体 1: 中心位置，长度约 15 像素，角度 30 度
center1 = [64, 64];
angle1 = 30 * pi / 180;
length1 = 15;
for i = -length1/2:length1/2
    x = round(center1(1) + i * cos(angle1));
    y = round(center1(2) + i * sin(angle1));
    if x > 0 && x <= im_size && y > 0 && y <= im_size
        Img(y, x) = Img(y, x) + 1.0 * exp(-0.1 * abs(i));
    end
end

% 散射体 2: 偏移位置，长度约 8 像素，角度 -45 度
center2 = [80, 50];
angle2 = -45 * pi / 180;
length2 = 8;
for i = -length2/2:length2/2
    x = round(center2(1) + i * cos(angle2));
    y = round(center2(2) + i * sin(angle2));
    if x > 0 && x <= im_size && y > 0 && y <= im_size
        Img(y, x) = Img(y, x) + 0.8 * exp(-0.15 * abs(i));
    end
end

% 散射体 3: 短散射体，长度约 3 像素
center3 = [40, 70];
Img(center3(2)-1:center3(2)+1, center3(1)) = 0.6;

% 平滑图像使其更真实
Img = imgaussfilt(Img, 1.0);
Img = max(0, Img);

fprintf('   Image size: %dx%d\n', size(Img,1), size(Img,2));
fprintf('   Number of synthetic scatterers: 3\n');
fprintf('   Expected lengths: ~15, ~8, ~3 pixels\n\n');

%% 2. 配置参数
fprintf('2. Configuring parameters...\n');

user_cfg = struct();
user_cfg.verbose = true;
user_cfg.max_atoms = 30;
user_cfg.aggregate_atoms = true;
user_cfg.refine_parameters = true;
user_cfg.use_gamma_decay = true;

fprintf('   aggregate_atoms: %d\n', user_cfg.aggregate_atoms);
fprintf('   refine_parameters: %d\n', user_cfg.refine_parameters);
fprintf('   use_gamma_decay: %d\n\n', user_cfg.use_gamma_decay);

%% 3. 执行 ASC 提取
fprintf('3. Running ASC extraction...\n');
tic;
[theta_list, residual_norms, info] = main_ascm_nomp_v6(Img, user_cfg);
elapsed = toc;
fprintf('   Extraction completed in %.3f seconds\n\n', elapsed);

%% 4. 显示结果
fprintf('4. Extraction Results:\n');
fprintf('   ========================================\n');
fprintf('   Total atoms extracted: %d\n', numel(theta_list));
fprintf('   Explained energy: %.2f%%\n', info.explained_ratio * 100);

if isfield(info, 'groups')
    fprintf('   Aggregated groups: %d\n\n', numel(info.groups));
else
    fprintf('   No aggregated groups\n\n');
end

% 显示每个原子的参数
fprintf('   Individual Atoms:\n');
fprintf('   --------------------------------------------------\n');
fprintf('   # |  Type  |  Pos(r,c)  |   A   |  L  | alpha | gamma |  phi(deg)\n');
fprintf('   --------------------------------------------------\n');

for i = 1:numel(theta_list)
    t = theta_list(i);
    fprintf('   %2d| %-7s| (%3d,%3d) | %5.2f | %4.1f| %5.2f | %5.3f | %6.1f\n', ...
        i, t.class(1:min(7,end)), t.row, t.col, t.A, t.L, t.alpha, t.gamma, t.phi*180/pi);
end
fprintf('   --------------------------------------------------\n\n');

%% 5. 参数统计分析
fprintf('5. Parameter Statistics:\n');
fprintf('   ========================================\n');

% 分离 localized 和 distributed 原子
localized_idx = find(strcmp({theta_list.class}, 'localized'));
distributed_idx = find(strcmp({theta_list.class}, 'distributed'));

fprintf('   Localized atoms: %d\n', numel(localized_idx));
if ~isempty(distributed_idx)
    L_vals = [theta_list(distributed_idx).L];
    alpha_vals = [theta_list(distributed_idx).alpha];
    gamma_vals = [theta_list(distributed_idx).gamma];
    
    fprintf('   Distributed atoms: %d\n', numel(distributed_idx));
    fprintf('   L statistics:\n');
    fprintf('     - Range: [%.2f, %.2f]\n', min(L_vals), max(L_vals));
    fprintf('     - Mean: %.2f\n', mean(L_vals));
    fprintf('     - Unique values: %d\n', numel(unique(L_vals)));
    
    fprintf('   Alpha statistics:\n');
    fprintf('     - Range: [%.3f, %.3f]\n', min(alpha_vals), max(alpha_vals));
    fprintf('     - Mean: %.3f\n', mean(alpha_vals));
    fprintf('     - All different from -1.0: %s\n', ...
        mat2str(all(alpha_vals ~= -1.0)));
    
    fprintf('   Gamma statistics:\n');
    fprintf('     - Range: [%.4f, %.4f]\n', min(gamma_vals), max(gamma_vals));
    fprintf('     - Mean: %.4f\n', mean(gamma_vals));
    fprintf('     - All positive: %s\n', mat2str(all(gamma_vals > 0)));
else
    fprintf('   No distributed atoms found\n');
end
fprintf('\n');

%% 6. 验收标准检查
fprintf('6. Acceptance Criteria Check:\n');
fprintf('   ========================================\n');

pass_count = 0;
total_tests = 4;

% 测试 1: L 不再固定为 1
if ~isempty(distributed_idx)
    L_vals = [theta_list(distributed_idx).L];
    test1_pass = any(L_vals ~= 1);
    fprintf('   [%s] Test 1: L values not fixed at 1\n', ...
        iif(test1_pass, 'PASS', 'FAIL'));
    if test1_pass, pass_count = pass_count + 1; end
else
    fprintf('   [SKIP] Test 1: No distributed atoms\n');
    total_tests = total_tests - 1;
end

% 测试 2: alpha 是有效数值（不是 -1.0）
if ~isempty(distributed_idx)
    alpha_vals = [theta_list(distributed_idx).alpha];
    test2_pass = all(alpha_vals ~= -1.0) && all(isfinite(alpha_vals));
    fprintf('   [%s] Test 2: Alpha values are valid (not -1.0)\n', ...
        iif(test2_pass, 'PASS', 'FAIL'));
    if test2_pass, pass_count = pass_count + 1; end
else
    fprintf('   [SKIP] Test 2: No distributed atoms\n');
    total_tests = total_tests - 1;
end

% 测试 3: gamma 是正数
if ~isempty(distributed_idx)
    gamma_vals = [theta_list(distributed_idx).gamma];
    test3_pass = all(gamma_vals > 0);
    fprintf('   [%s] Test 3: Gamma values are positive\n', ...
        iif(test3_pass, 'PASS', 'FAIL'));
    if test3_pass, pass_count = pass_count + 1; end
else
    fprintf('   [SKIP] Test 3: No distributed atoms\n');
    total_tests = total_tests - 1;
end

% 测试 4: 提取到了原子
test4_pass = numel(theta_list) > 0;
fprintf('   [%s] Test 4: Atoms extracted successfully\n', ...
    iif(test4_pass, 'PASS', 'FAIL'));
if test4_pass, pass_count = pass_count + 1; end

fprintf('\n   Overall: %d/%d tests passed\n', pass_count, total_tests);
fprintf('   ========================================\n\n');

%% 7. 可视化
fprintf('7. Generating visualization...\n');

figure('Name', 'ASC Extraction Results', 'Position', [100, 100, 1200, 400]);

% 子图 1: 原始图像
subplot(1,3,1);
imagesc(Img); axis image; colormap(gca, 'gray'); colorbar;
title('Input Image');
xlabel('Column'); ylabel('Row');

% 子图 2: 提取的散射中心位置
subplot(1,3,2);
imagesc(Img); axis image; colormap(gca, 'gray'); colorbar;
hold on;
for i = 1:numel(theta_list)
    t = theta_list(i);
    if strcmp(t.class, 'localized')
        plot(t.col, t.row, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
    else
        % 绘制线段表示方向和长度
        dx = t.L/2 * cos(t.phi);
        dy = -t.L/2 * sin(t.phi);
        plot([t.col-dx, t.col+dx], [t.row-dy, t.row+dy], ...
            'g-', 'LineWidth', 2);
        plot(t.col, t.row, 'go', 'MarkerSize', 6, 'MarkerFaceColor', 'g');
    end
end
hold off;
title(sprintf('Extracted Scatterers (N=%d)', numel(theta_list)));
xlabel('Column'); ylabel('Row');
legend('Localized', 'Distributed');

% 子图 3: 参数分布
subplot(1,3,3);
if ~isempty(distributed_idx)
    L_vals = [theta_list(distributed_idx).L];
    alpha_vals = [theta_list(distributed_idx).alpha];
    gamma_vals = [theta_list(distributed_idx).gamma];
    
    yyaxis left;
    plot(distributed_idx, L_vals, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
    ylabel('L (Length)');
    ylim([0, max(L_vals)*1.2]);
    
    yyaxis right;
    plot(distributed_idx, alpha_vals, 'rs-', 'LineWidth', 2, 'MarkerSize', 8);
    hold on;
    plot(distributed_idx, gamma_vals*10, 'g^-', 'LineWidth', 2, 'MarkerSize', 8);
    hold off;
    ylabel('Alpha / Gamma×10');
    
    xlabel('Atom Index');
    title('Parameter Distribution');
    legend('L', 'Alpha', 'Gamma×10', 'Location', 'best');
    grid on;
else
    text(0.5, 0.5, 'No distributed atoms', ...
        'HorizontalAlignment', 'center', 'FontSize', 14);
    axis off;
end

fprintf('   Visualization complete.\n\n');

fprintf('=== Test Complete ===\n');

%% Helper function
function out = iif(condition, true_val, false_val)
    if condition
        out = true_val;
    else
        out = false_val;
    end
end
