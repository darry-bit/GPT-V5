 function [theta_list, residual_norms, info] = main_ascm_nomp_v6(Img, user_cfg)
      if nargin < 2 || isempty(user_cfg), user_cfg = struct(); end
      cfg = v6_default_cfg();
      fn  = fieldnames(user_cfg);
      for i = 1:numel(fn), cfg.(fn{i}) = user_cfg.(fn{i}); end

      % 预处理
      I = abs(double(Img));
      if ~all(isfinite(I(:))), error('Input Img contains NaN/Inf.'); end
      if cfg.do_log, I = log1p(I); end
      I = I - min(I(:)); if max(I(:))>0, I = I./max(I(:)); end
      if cfg.do_smooth
          % 创建高斯滤波器（不依赖工具箱）
          kernel_size = ceil(3*cfg.sigma_smooth)*2 + 1;
          half = floor(kernel_size/2);
          [x, y] = meshgrid(-half:half, -half:half);
          kernel = exp(-(x.^2 + y.^2)/(2*cfg.sigma_smooth^2));
          kernel = kernel / sum(kernel(:));
          I = conv2(I, kernel, 'same');
      end
      cfg = v6_validate_cfg(cfg, size(I));

      roi_mask = v6_build_roi(I, cfg);
      dict = v6_build_dict(cfg);
      [theta_list, residual_norms] = v6_run_nomp(I, roi_mask, dict, cfg);

      % 原子聚合和参数估计
      do_aggregate = isfield(cfg, 'aggregate_atoms') && cfg.aggregate_atoms;
      if do_aggregate && ~isempty(theta_list)
          groups = v6_aggregate_atoms(theta_list, cfg);
          if cfg.verbose && ~isempty(groups)
              fprintf('[AGGREGATE] Found %d groups:\n', numel(groups));
              for g = 1:numel(groups)
                  fprintf('  Group %d: atoms [%s]\n', g, num2str(groups{g}));
              end
          end
          do_refine = isfield(cfg, 'refine_parameters') && cfg.refine_parameters;
          if do_refine
              theta_list = v6_refine_L_alpha_gamma(I, roi_mask, theta_list, groups, dict, cfg);
          end
      else
          groups = {};
      end

      % 解释能量
      R0 = I .* roi_mask; E0 = sum(R0(:).^2);
      if isempty(theta_list)
          E_res = E0;
      else
          R = R0; [H,W] = size(I);
          for k = 1:numel(theta_list)
              atom = dict.atom(:,:,theta_list(k).dict_idx);
              patch_full = zeros(H,W);
              patch_full = v6_paste_patch(patch_full, atom, theta_list(k).row, theta_list(k).col);
              R = R - theta_list(k).A * patch_full;
          end
          E_res = sum(R(:).^2);
      end
      info = struct();
      info.explained_ratio = max(0, 1 - E_res / max(E0, eps));
      info.cfg             = cfg;
      info.num_atoms       = numel(theta_list);
      info.roi_mask        = roi_mask;
      info.dict            = dict;
      info.groups          = groups;
      if cfg.verbose
          fprintf('--- ASC summary ---\n');
          fprintf('  atoms       = %d\n', info.num_atoms);
          fprintf('  explained   = %.2f %%\n', 100*info.explained_ratio);
      end
  end

  %% 默认参数：CornerDiff + 最小线元
  function cfg = v6_default_cfg()
      cfg = struct();
      cfg.do_log       = false;
      cfg.do_smooth    = true;
      cfg.sigma_smooth = 0.8;

      cfg.auto_roi          = true;
      cfg.roi_thresh_scale  = 1.0;
      cfg.roi_dilate        = 2;
      cfg.roi_margin        = 4;

      cfg.patch_size    = 17;           % 奇数
      cfg.L_list        = 1;            % 最小线元
      cfg.phi_list_deg  = -80:5:80;     % 粗角度网格
      cfg.sigma_loc     = 1.0;
      cfg.sigma_short_thick = 1.0;
      cfg.sigma_short_thin  = 0.6;      % 未用但保留

      % 阈值与上限（分开控制）
      cfg.min_ncc_localized   = 0.50;   % CornerDiff 更严格
      cfg.min_ncc_distributed = 0.20;   % 线元放宽
      cfg.min_gain_ratio      = 0.002;  % 单线元增益下限
      cfg.stop_explained_ratio= 0.80;
      cfg.max_atoms           = 24;
      cfg.max_localized       = 1;      % CornerDiff 最多 1 个
      cfg.max_iters           = 128;

      cfg.do_phi_refine        = true;
      cfg.phi_refine_range_deg = 15;
      cfg.phi_refine_step_deg  = 1;

      % Gamma 衰减参数
      cfg.use_gamma_decay      = true;  % 是否使用 gamma 衰减建模
      cfg.gamma_default        = 0.1;   % 默认 gamma 值
      cfg.gamma_range          = [0.01, 1.0]; % gamma 搜索范围
      
      % 原子聚合参数
      cfg.aggregate_atoms      = true;  % 是否进行原子聚合
      cfg.aggregate_dist_thresh = 3.0;  % 聚合距离阈值（像素）
      cfg.aggregate_angle_thresh = 15;  % 聚合角度阈值（度）
      
      % 参数精化
      cfg.refine_parameters    = true;  % 是否进行 L, alpha, gamma 精化

      cfg.verbose = true;
      cfg.dx=1; cfg.dy=1; cfg.x0=0; cfg.y0=0;
  end

  function cfg = v6_validate_cfg(cfg, im_size)
      if mod(cfg.patch_size,2)==0 || cfg.patch_size<=0, error('patch_size must be positive odd.'); end
      if any(cfg.L_list<=0), error('L_list must be positive.'); end
      if max(cfg.L_list) >= cfg.patch_size, error('max(L_list) must be < patch_size.'); end
      if isempty(cfg.phi_list_deg) || any(~isfinite(cfg.phi_list_deg)), error('phi_list_deg must be finite.'); end
      cfg.max_atoms = max(1, floor(cfg.max_atoms));
      cfg.max_iters = max(1, floor(cfg.max_iters));
      cfg.max_localized = max(0, floor(cfg.max_localized));
      cfg.min_gain_ratio = max(cfg.min_gain_ratio,0);
      cfg.roi_dilate = max(0, round(cfg.roi_dilate));
      cfg.roi_margin = max(0, round(cfg.roi_margin));
      if numel(im_size)~=2 || any(im_size<cfg.patch_size), error('Image too small for patch_size.'); end
  end

  function roi_mask = v6_build_roi(I, cfg)
      [H,W] = size(I);
      if ~cfg.auto_roi, roi_mask = true(H,W); return; end
      if max(I(:)) <= min(I(:)), roi_mask = true(H,W); return; end
      T  = v6_graythresh(I) * cfg.roi_thresh_scale;
      roi = I > T;
      roi = v6_imclose(roi, 2);
      roi = v6_imfill(roi);
      if ~any(roi(:)), roi_mask = true(H,W); return; end
      if cfg.roi_dilate>0, roi = v6_imdilate(roi, cfg.roi_dilate); end
      CC = v6_bwconncomp(roi);
      if CC.NumObjects > 1
          lens = cellfun(@numel, CC.PixelIdxList); [~, idxMax] = max(lens);
          roi = false(H,W); roi(CC.PixelIdxList{idxMax}) = true;
      end
      [r,c] = find(roi);
      r1 = max(1, min(r)-cfg.roi_margin); r2 = min(H, max(r)+cfg.roi_margin);
      c1 = max(1, min(c)-cfg.roi_margin); c2 = min(W, max(c)+cfg.roi_margin);
      roi_mask = false(H,W); roi_mask(r1:r2, c1:c2) = true;
  end

  function dict = v6_build_dict(cfg)
      p  = cfg.patch_size; h  = (p-1)/2; [xx,yy] = meshgrid(-h:h, -h:h);
      atoms = []; geom  = {}; cls = {}; Ls=[]; phis=[]; idx=0;

      % CornerDiff (L=0)
      sigma = cfg.sigma_loc;
      g = exp(-(xx.^2 + yy.^2)/(2*sigma^2));
      g = g / max(norm(g(:)), eps);
      idx=idx+1; atoms(:,:,idx)=g; geom{idx}='CornerDiff'; cls{idx}='localized'; Ls(idx)=0; phis(idx)=0;

      % 最小线元（ED_thick，L=1）
      for ang = cfg.phi_list_deg(:).'
          phi = deg2rad(ang);
          patch = v6_build_ed_patch(cfg, 1, phi, 'ED_thick', xx, yy);
          idx=idx+1; atoms(:,:,idx)=patch; geom{idx}='ED_thick'; cls{idx}='distributed'; Ls(idx)=1; phis(idx)=phi;
      end

      dict = struct('atom',atoms,'geom',{geom},'class',{cls},'L',Ls,'phi',phis,'num',size(atoms,3));
      if cfg.verbose
          fprintf('[DICT-line] patch_size=%d, entries=%d (%d localized + %d distributed)\n', ...
              cfg.patch_size, dict.num, 1, dict.num-1);
      end
  end

  function patch = v6_build_ed_patch(cfg, L, phi, geom_type, xx, yy, gamma)
      if nargin < 5
          p = cfg.patch_size; h = (p-1)/2; [xx,yy] = meshgrid(-h:h, -h:h);
      end
      if nargin < 7
          gamma = cfg.gamma_default; % 使用默认 gamma 值
      end
      L_eff = max(L, 1); % 保底 1 像素代表 0.1
      switch geom_type
          case 'ED_thick'
              sigma_long  = max(L_eff/3, 0.8);
              sigma_short = cfg.sigma_short_thick;
          case 'ED_thin'
              sigma_long  = max(L_eff/2, 0.8);
              sigma_short = cfg.sigma_short_thin;
          otherwise
              error('Unknown geom_type: %s', geom_type);
      end
      x_rot =  xx*cos(phi) + yy*sin(phi);
      y_rot = -xx*sin(phi) + yy*cos(phi);
      patch = exp(-(x_rot.^2/(2*sigma_long^2) + y_rot.^2/(2*sigma_short^2)));
      
      % 应用 gamma 衰减（沿散射体长度方向）
      if cfg.use_gamma_decay && gamma > 0
          % gamma 衰减模型: exp(-gamma * |x_rot|)
          decay = exp(-gamma * abs(x_rot));
          patch = patch .* decay;
      end
      
      patch(abs(x_rot) > (L_eff/2 + 1)) = 0;
      patch = patch / max(norm(patch(:)), eps);
  end

  function [theta_list, residual_norms] = v6_run_nomp(I, roi_mask, dict, cfg)
      [H,W] = size(I);
      R  = I .* roi_mask; E0 = sum(R(:).^2);
      residual_norms = [];
      theta_list = struct('row',{},'col',{},'A',{},'geom',{},'class',{},'L',{},'alpha',{},'gamma',{},'phi',{},'dict_idx',{},'x',{},'y',{});
      if E0 <= 0, return; end
      n_localized = 0;

      for it = 1:cfg.max_iters
          best = struct('ncc', -Inf);
          for j = 1:dict.num
              class_type = dict.class{j};
              geom_type  = dict.geom{j};
              patch = dict.atom(:,:,j);

              % 局域原子数量限制
              if strcmp(class_type,'localized') && n_localized >= cfg.max_localized
                  continue;
              end

              C = conv2(R, rot90(patch,2), 'same');
              C(~roi_mask) = -Inf;
              [val, idx_lin] = max(C(:));
              if ~isfinite(val), continue; end
              [r, c] = ind2sub([H,W], idx_lin);
              [R_patch, mask_patch] = v6_get_patch(R, roi_mask, r, c, cfg.patch_size);
              P_patch = patch .* mask_patch;
              num = sum(R_patch(:) .* P_patch(:));
              den = sqrt(sum(R_patch(:).^2) * sum(P_patch(:).^2));
              if den <= eps, continue; end
              ncc = num / den;

              % 分开阈值
              min_ncc_req = cfg.min_ncc_distributed;
              if strcmp(class_type,'localized'), min_ncc_req = cfg.min_ncc_localized; end
              if ncc < min_ncc_req, continue; end

              if ncc > best.ncc
                  best.ncc=ncc; best.val=num; best.r=r; best.c=c; best.idx=j;
              end
          end
          if ~isfield(best,'idx')
              if cfg.verbose, fprintf('[NOMP] stop: no candidate above NCC thresholds\n'); end
              break;
          end

          geom_type  = dict.geom{best.idx};
          class_type = dict.class{best.idx};
          L0         = dict.L(best.idx);
          phi0       = dict.phi(best.idx);

          if cfg.do_phi_refine && strcmp(class_type,'distributed')
              [phi_ref, patch_ref] = v6_refine_phi_joint(R, roi_mask, best.r, best.c, geom_type, L0, phi0, cfg);
          else
              phi_ref  = phi0; patch_ref = dict.atom(:,:,best.idx);
          end

          % 局部 LS 幅度估计
          [R_patch_full, mask_patch_full] = v6_get_patch(R, roi_mask, best.r, best.c, cfg.patch_size);
          P = patch_ref .* mask_patch_full;
          den_loc = sum(P(:).^2);
          if den_loc <= eps
              if cfg.verbose, fprintf('[NOMP] reject: zero-energy atom at (%d,%d)\n', best.r, best.c); end
              break;
          end
          num_loc = sum(R_patch_full(:) .* P(:));
          A = num_loc / den_loc;

          atom_full = zeros(H,W);
          atom_full = v6_paste_patch(atom_full, patch_ref, best.r, best.c);

          E_before = sum(R(:).^2);
          R_new    = R - A * atom_full;
          E_after  = sum(R_new(:).^2);
          gain     = max(0, (E_before - E_after) / max(E0, eps));
          if gain < cfg.min_gain_ratio || ~isfinite(gain)
              if cfg.verbose
                  fprintf('[NOMP] reject: gain=%.3f%% < %.3f%% (atom=%s)\n', gain*100, cfg.min_gain_ratio*100,geom_type);
              end
              break;
          end

          R = R_new;
          residual_norms(end+1) = sqrt(E_after); %#ok<AGROW>
          [x_m, y_m] = v6_rowcol_to_xy(best.r, best.c, H, W, cfg);
          t = struct('row',best.r,'col',best.c,'A',A,'geom',geom_type,'class',class_type,...
                     'L',L0,'alpha',-1.0,'gamma',cfg.gamma_default,'phi',phi_ref,'dict_idx',best.idx,'x',x_m,'y',y_m);
          theta_list(end+1) = t; %#ok<AGROW>
          if strcmp(class_type,'localized'), n_localized = n_localized + 1; end

          explained = max(0, 1 - E_after / max(E0, eps));
          if cfg.verbose
              fprintf('[NOMP] #%d accept: %s NCC=%.3f A=%.3f gain=%.2f%% @(%.0f,%.0f) phi=%.1f\n', ...
                  numel(theta_list), geom_type, best.ncc, A, gain*100, t.row, t.col, t.phi*180/pi);
          end
          if numel(theta_list) >= cfg.max_atoms || explained >= cfg.stop_explained_ratio
              if cfg.verbose, fprintf('[NOMP] stop: atoms=%d explained=%.2f%%\n', numel(theta_list), explained*100); end
              break;
          end
      end
  end

  function [phi_best, patch_best] = v6_refine_phi_joint(R, roi_mask, r, c, geom_type, L, phi_init, cfg)
      range_deg = cfg.phi_refine_range_deg; step_deg  = cfg.phi_refine_step_deg;
      phi0_deg = phi_init*180/pi;
      phi_cands = deg2rad(phi0_deg - range_deg : step_deg : phi0_deg + range_deg);
      [R_patch, mask_patch] = v6_get_patch(R, roi_mask, r, c, cfg.patch_size);
      best_cost = Inf; best_idx=1; costs = inf(size(phi_cands));
      for ii = 1:numel(phi_cands)
          phi = phi_cands(ii);
          patch0 = v6_build_ed_patch(cfg, L, phi, geom_type);
          P = patch0 .* mask_patch;
          den = sum(P(:).^2); if den <= eps, continue; end
          num = sum(R_patch(:).*P(:));
          A_phi = num / den;
          Res  = R_patch - A_phi*P;
          cost = sum(Res(:).^2);
          costs(ii)=cost;
          if cost < best_cost, best_cost=cost; best_idx=ii; end
      end
      phi_best = phi_cands(best_idx);
      if numel(phi_cands)>=3 && best_idx>1 && best_idx<numel(phi_cands)
          c1=costs(best_idx-1); c2=costs(best_idx); c3=costs(best_idx+1);
          denom = (c1 - 2*c2 + c3);
          if abs(denom)>eps
              delta = 0.25*(c1 - c3)/denom; delta = max(-0.75, min(0.75, delta));
              phi_best = phi_cands(best_idx) + delta*deg2rad(step_deg);
          end
      end
      patch_best = v6_build_ed_patch(cfg, L, phi_best, geom_type);
      if ~isfinite(best_cost)
          phi_best = phi_init;
          patch_best = v6_build_ed_patch(cfg, L, phi_init, geom_type);
      end
  end

  function [patch, mask_patch] = v6_get_patch(Im, roi_mask, rC, cC, patch_size)
      [H,W] = size(Im); h = (patch_size-1)/2;
      r1 = max(1, rC-h); r2 = min(H, rC+h);
      c1 = max(1, cC-h); c2 = min(W, cC+h);
      patch = zeros(patch_size, patch_size);
      mask_patch = zeros(patch_size, patch_size);
      rr = r1:r2; cc = c1:c2; pr = rr - rC + h + 1; pc = cc - cC + h + 1;
      patch(pr,pc)      = Im(rr,cc);
      mask_patch(pr,pc) = roi_mask(rr,cc);
  end

  function Im = v6_paste_patch(Im, patch, rC, cC)
      [H,W] = size(Im); [ph,pw] = size(patch); h = (ph-1)/2;
      r1 = max(1, rC-h); r2 = min(H, rC+h); c1 = max(1, cC-h); c2 = min(W, cC+h);
      rr = r1:r2; cc = c1:c2; pr = rr - rC + h + 1; pc = cc - cC + h + 1;
      Im(rr,cc) = Im(rr,cc) + patch(pr,pc);
  end

  function [x, y] = v6_rowcol_to_xy(r, c, H, W, cfg)
      if ~isfield(cfg,'dx'), cfg.dx=1; end
      if ~isfield(cfg,'dy'), cfg.dy=1; end
      if ~isfield(cfg,'x0'), cfg.x0=0; end
      if ~isfield(cfg,'y0'), cfg.y0=0; end
      c0 = (W + 1) / 2; r0 = (H + 1) / 2;
      x = (c - c0) * cfg.dx + cfg.x0;
      y = (r0 - r) * cfg.dy + cfg.y0;
  end

  %% 原子聚合函数
  function groups = v6_aggregate_atoms(theta_list, cfg)
      % 将属于同一散射体的原子聚合成组
      % 输入：theta_list - 原子列表
      %       cfg - 配置参数
      % 输出：groups - 聚合组，每组包含原子索引列表
      
      if isempty(theta_list)
          groups = {};
          return;
      end
      
      n = numel(theta_list);
      groups = {};
      assigned = false(1, n);
      
      % 只聚合 distributed 类型的原子
      for i = 1:n
          if assigned(i) || strcmp(theta_list(i).class, 'localized')
              continue;
          end
          
          % 创建新组
          group = i;
          assigned(i) = true;
          
          % 查找相邻且方向相近的原子
          for j = i+1:n
              if assigned(j) || strcmp(theta_list(j).class, 'localized')
                  continue;
              end
              
              % 检查是否应该聚合
              should_aggregate = false;
              for k = group
                  % 计算距离
                  dr = theta_list(j).row - theta_list(k).row;
                  dc = theta_list(j).col - theta_list(k).col;
                  dist = sqrt(dr^2 + dc^2);
                  
                  % 计算角度差
                  angle_diff = abs(theta_list(j).phi - theta_list(k).phi) * 180 / pi;
                  if angle_diff > 180
                      angle_diff = 360 - angle_diff;
                  end
                  
                  % 如果距离和角度都满足阈值，则聚合
                  if dist <= cfg.aggregate_dist_thresh && angle_diff <= cfg.aggregate_angle_thresh
                      should_aggregate = true;
                      break;
                  end
              end
              
              if should_aggregate
                  group(end+1) = j; %#ok<AGROW>
                  assigned(j) = true;
              end
          end
          
          % 只有多个原子的组才保存（单原子不需要聚合）
          if numel(group) > 1
              groups{end+1} = group; %#ok<AGROW>
          end
      end
  end

  %% 估计散射系数 alpha
  function alpha = v6_estimate_alpha(I, roi_mask, atom_indices, theta_list, dict, cfg)
      % 基于散射模型和幅度信息估计散射系数 alpha
      % 输入：I - 输入图像
      %       roi_mask - ROI 掩码
      %       atom_indices - 原子索引列表
      %       theta_list - 原子参数列表
      %       dict - 字典
      %       cfg - 配置
      % 输出：alpha - 估计的散射系数
      
      if isempty(atom_indices)
          alpha = -1.0;
          return;
      end
      
      % 计算组的平均幅度和 NCC
      total_A = 0;
      total_ncc = 0;
      for idx = atom_indices
          total_A = total_A + abs(theta_list(idx).A);
          
          % 重新计算 NCC
          r = theta_list(idx).row;
          c = theta_list(idx).col;
          patch = dict.atom(:,:,theta_list(idx).dict_idx);
          [R_patch, mask_patch] = v6_get_patch(I, roi_mask, r, c, cfg.patch_size);
          P_patch = patch .* mask_patch;
          num = sum(R_patch(:) .* P_patch(:));
          den = sqrt(sum(R_patch(:).^2) * sum(P_patch(:).^2));
          if den > eps
              total_ncc = total_ncc + (num / den);
          end
      end
      
      avg_A = total_A / numel(atom_indices);
      avg_ncc = total_ncc / numel(atom_indices);
      
      % 基于幅度和 NCC 估计 alpha
      % alpha 与幅度和匹配质量相关
      alpha = avg_A * avg_ncc;
      
      % 确保 alpha 在合理范围内
      alpha = max(0.01, min(10.0, alpha));
  end

  %% 估计衰减因子 gamma
  function gamma = v6_estimate_gamma(I, roi_mask, atom_indices, theta_list, dict, cfg)
      % 通过幅度沿散射体方向的衰减特性估计 gamma
      % 输入：I - 输入图像
      %       roi_mask - ROI 掩码
      %       atom_indices - 原子索引列表
      %       theta_list - 原子参数列表
      %       dict - 字典
      %       cfg - 配置
      % 输出：gamma - 估计的衰减因子
      
      if isempty(atom_indices) || numel(atom_indices) < 2
          gamma = cfg.gamma_default;
          return;
      end
      
      % 提取原子位置和幅度
      positions = zeros(numel(atom_indices), 2);
      amplitudes = zeros(numel(atom_indices), 1);
      for i = 1:numel(atom_indices)
          idx = atom_indices(i);
          positions(i,:) = [theta_list(idx).row, theta_list(idx).col];
          amplitudes(i) = abs(theta_list(idx).A);
      end
      
      % 计算平均角度
      avg_phi = 0;
      for idx = atom_indices
          avg_phi = avg_phi + theta_list(idx).phi;
      end
      avg_phi = avg_phi / numel(atom_indices);
      
      % 沿主方向对原子排序
      proj = positions(:,1) * cos(avg_phi) + positions(:,2) * sin(avg_phi);
      [~, sort_idx] = sort(proj);
      
      % 计算相邻原子的幅度衰减
      decay_rates = [];
      for i = 1:numel(sort_idx)-1
          idx1 = sort_idx(i);
          idx2 = sort_idx(i+1);
          A1 = amplitudes(idx1);
          A2 = amplitudes(idx2);
          dist = norm(positions(idx2,:) - positions(idx1,:));
          
          if A1 > eps && A2 > eps && dist > eps
              % gamma 估计: A2/A1 = exp(-gamma*dist)
              decay_rate = -log(A2/A1 + eps) / (dist + eps);
              if decay_rate > 0 && decay_rate < 2.0 % 合理范围
                  decay_rates(end+1) = decay_rate; %#ok<AGROW>
              end
          end
      end
      
      if ~isempty(decay_rates)
          gamma = median(decay_rates);
          gamma = max(cfg.gamma_range(1), min(cfg.gamma_range(2), gamma));
      else
          gamma = cfg.gamma_default;
      end
  end

  %% 参数联合精化
  function theta_list = v6_refine_L_alpha_gamma(I, roi_mask, theta_list, groups, dict, cfg)
      % 对聚合后的原子组进行 (L, alpha, gamma) 联合优化
      % 输入：I - 输入图像
      %       roi_mask - ROI 掩码
      %       theta_list - 原子参数列表
      %       groups - 聚合组
      %       dict - 字典
      %       cfg - 配置
      % 输出：theta_list - 精化后的原子参数列表
      
      if isempty(groups)
          return;
      end
      
      for g = 1:numel(groups)
          atom_indices = groups{g};
          if isempty(atom_indices)
              continue;
          end
          
          % 估计等效长度 L
          % 基于原子数量和空间分布
          if numel(atom_indices) > 1
              positions = zeros(numel(atom_indices), 2);
              for i = 1:numel(atom_indices)
                  idx = atom_indices(i);
                  positions(i,:) = [theta_list(idx).row, theta_list(idx).col];
              end
              
              % 计算主方向上的跨度
              avg_phi = mean([theta_list(atom_indices).phi]);
              proj = positions(:,1) * cos(avg_phi) + positions(:,2) * sin(avg_phi);
              L_eff = max(proj) - min(proj) + 1; % +1 因为包含端点
              L_eff = max(1, L_eff); % 至少为 1
          else
              L_eff = 1;
          end
          
          % 估计 alpha
          alpha = v6_estimate_alpha(I, roi_mask, atom_indices, theta_list, dict, cfg);
          
          % 估计 gamma
          gamma = v6_estimate_gamma(I, roi_mask, atom_indices, theta_list, dict, cfg);
          
          % 更新组内所有原子的参数
          for idx = atom_indices
              theta_list(idx).L = L_eff;
              theta_list(idx).alpha = alpha;
              theta_list(idx).gamma = gamma;
          end
      end
  end

  %% 简化的 Otsu 阈值算法（不依赖 image 包）
  function thresh = v6_graythresh(I)
      % 归一化到 [0, 1]
      I = double(I);
      I = (I - min(I(:))) / (max(I(:)) - min(I(:)) + eps);
      
      % 计算直方图
      nbins = 256;
      counts = histc(I(:), linspace(0, 1, nbins));
      p = counts / sum(counts);
      
      % Otsu 方法寻找最佳阈值
      omega = cumsum(p);
      mu = cumsum(p .* (1:nbins)');
      mu_t = mu(end);
      
      sigma_b_squared = (mu_t * omega - mu).^2 ./ (omega .* (1 - omega) + eps);
      [~, idx] = max(sigma_b_squared);
      thresh = (idx - 1) / (nbins - 1);
  end

  %% 简化的形态学膨胀（不依赖 image 包）
  function BW = v6_imdilate(BW, radius)
      if radius <= 0, return; end
      % 创建圆形结构元素
      [x, y] = meshgrid(-radius:radius, -radius:radius);
      se = (x.^2 + y.^2) <= radius^2;
      % 膨胀操作
      BW = double(BW);
      BW = conv2(BW, double(se), 'same') > 0;
  end

  %% 简化的形态学闭运算（不依赖 image 包）
  function BW = v6_imclose(BW, radius)
      % 先膨胀后腐蚀
      BW = v6_imdilate(BW, radius);
      BW = v6_imerode(BW, radius);
  end

  %% 简化的形态学腐蚀（不依赖 image 包）
  function BW = v6_imerode(BW, radius)
      if radius <= 0, return; end
      % 创建圆形结构元素
      [x, y] = meshgrid(-radius:radius, -radius:radius);
      se = (x.^2 + y.^2) <= radius^2;
      % 腐蚀操作
      BW = double(BW);
      kernel_sum = sum(se(:));
      result = conv2(BW, double(se), 'same');
      BW = result >= kernel_sum;
  end

  %% 简化的孔洞填充（不依赖 image 包）
  function BW = v6_imfill(BW)
      % 使用泛洪填充算法
      [H, W] = size(BW);
      filled = BW;
      
      % 从边界开始标记外部区域
      border = false(H, W);
      border(1,:) = ~BW(1,:);
      border(end,:) = ~BW(end,:);
      border(:,1) = border(:,1) | ~BW(:,1);
      border(:,W) = border(:,W) | ~BW(:,W);
      
      % 简单的4连通泛洪填充
      old_sum = 0;
      new_sum = sum(border(:));
      iter = 0;
      max_iter = H * W;
      
      while new_sum ~= old_sum && iter < max_iter
          old_sum = new_sum;
          % 膨胀标记但限制在非对象区域
          border_dilated = v6_imdilate(border, 1);
          border = border_dilated & ~BW;
          new_sum = sum(border(:));
          iter = iter + 1;
      end
      
      % 填充：所有未标记的非对象像素都是孔洞
      filled = BW | ~border;
  end

  %% 简化的连通分量标记（不依赖 image 包）
  function CC = v6_bwconncomp(BW)
      % 4连通分量标记
      [H, W] = size(BW);
      labeled = zeros(H, W);
      label = 0;
      
      for i = 1:H
          for j = 1:W
              if BW(i,j) && labeled(i,j) == 0
                  % 新的连通分量
                  label = label + 1;
                  % 使用栈进行深度优先搜索
                  stack = [i, j];
                  pixels = [];
                  
                  while ~isempty(stack)
                      r = stack(1,1);
                      c = stack(1,2);
                      stack(1,:) = [];
                      
                      if r < 1 || r > H || c < 1 || c > W
                          continue;
                      end
                      if ~BW(r,c) || labeled(r,c) > 0
                          continue;
                      end
                      
                      labeled(r,c) = label;
                      pixels = [pixels; sub2ind([H,W], r, c)]; %#ok<AGROW>
                      
                      % 添加4连通邻居
                      stack = [stack; r-1, c; r+1, c; r, c-1; r, c+1]; %#ok<AGROW>
                  end
                  
                  CC.PixelIdxList{label} = pixels;
              end
          end
      end
      
      if label == 0
          CC.NumObjects = 0;
          CC.PixelIdxList = {};
      else
          CC.NumObjects = label;
      end
  end