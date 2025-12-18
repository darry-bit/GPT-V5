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
      if cfg.do_smooth, I = imgaussfilt(I, cfg.sigma_smooth); end
      cfg = v6_validate_cfg(cfg, size(I));

      roi_mask = v6_build_roi(I, cfg);
      dict = v6_build_dict(cfg);
      [theta_list, residual_norms] = v6_run_nomp(I, roi_mask, dict, cfg);
      
      % 聚合分布式原子并细化参数
      if ~isempty(theta_list)
          theta_list = v6_aggregate_atoms(theta_list, cfg);
          theta_list = v6_refine_L_alpha_gamma(theta_list, I, roi_mask, cfg);
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
      cfg.max_localized       = Inf;    % 允许多个局域式散射中心
      cfg.max_iters           = 128;
      
      % Gamma 估计参数
      cfg.gamma_localized_range = [0.1, 2.0];  % 局域式 gamma 范围
      cfg.gamma_distributed     = 0;            % 分布式 gamma 固定为 0

      cfg.do_phi_refine        = true;
      cfg.phi_refine_range_deg = 15;
      cfg.phi_refine_step_deg  = 1;
      
      % 聚合参数
      cfg.agg_dist_thresh  = 6;   % 放宽距离阈值
      cfg.agg_angle_thresh = 30;  % 放宽角度阈值（度）

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
      T  = graythresh(I) * cfg.roi_thresh_scale;
      roi = I > T;
      roi = imclose(roi, strel('disk',2));
      roi = imfill(roi, 'holes');
      if ~any(roi(:)), roi_mask = true(H,W); return; end
      if cfg.roi_dilate>0, roi = imdilate(roi, strel('disk', cfg.roi_dilate)); end
      CC = bwconncomp(roi);
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

  function patch = v6_build_ed_patch(cfg, L, phi, geom_type, xx, yy)
      if nargin < 5
          p = cfg.patch_size; h = (p-1)/2; [xx,yy] = meshgrid(-h:h, -h:h);
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
          
          % 初始 gamma：distributed 用 0，localized 用默认值（后续细化）
          gamma_init = 0.0;
          if strcmp(class_type,'localized')
              gamma_init = 0.5;
          end
          
          t = struct('row',best.r,'col',best.c,'A',A,'geom',geom_type,'class',class_type,...
                     'L',L0,'alpha',-1.0,'gamma',gamma_init,'phi',phi_ref,'dict_idx',best.idx,'x',x_m,'y',y_m);
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

  %% 聚合分布式原子（Union-Find 算法）
  function theta_list = v6_aggregate_atoms(theta_list, cfg)
      N = numel(theta_list);
      if N == 0, return; end
      
      % 分离局域式和分布式原子
      is_distributed = false(1, N);
      for i = 1:N
          % 通过 L 判断类型：L=0 为局域式，L>0 为分布式
          if theta_list(i).L > 0
              is_distributed(i) = true;
          end
      end
      
      dist_idx = find(is_distributed);
      if isempty(dist_idx)
          return; % 没有分布式原子，无需聚合
      end
      
      % Union-Find 初始化
      parent = 1:numel(dist_idx);
      
      function root = find_root(x)
          if parent(x) ~= x
              parent(x) = find_root(parent(x));
          end
          root = parent(x);
      end
      
      function union(x, y)
          rx = find_root(x);
          ry = find_root(y);
          if rx ~= ry
              parent(ry) = rx;
          end
      end
      
      % 聚合条件：距离和角度
      for i = 1:numel(dist_idx)
          for j = i+1:numel(dist_idx)
              idx_i = dist_idx(i);
              idx_j = dist_idx(j);
              
              % 计算距离
              dx = theta_list(idx_i).x - theta_list(idx_j).x;
              dy = theta_list(idx_i).y - theta_list(idx_j).y;
              dist = sqrt(dx^2 + dy^2);
              
              if dist > cfg.agg_dist_thresh
                  continue;
              end
              
              % 计算角度差
              phi_i = theta_list(idx_i).phi;
              phi_j = theta_list(idx_j).phi;
              angle_diff = abs(mod(phi_i - phi_j + pi, 2*pi) - pi) * 180 / pi;
              
              if angle_diff > cfg.agg_angle_thresh
                  continue;
              end
              
              % 检查连接方向与散射体方向的一致性
              conn_angle = atan2(dy, dx);
              phi_avg = atan2(sin(phi_i) + sin(phi_j), cos(phi_i) + cos(phi_j));
              conn_diff = abs(mod(conn_angle - phi_avg + pi, 2*pi) - pi) * 180 / pi;
              
              % 连接方向应与散射体方向接近（允许较大偏差）
              if conn_diff < 60 || conn_diff > 120
                  union(i, j);
              end
          end
      end
      
      % 聚合成组
      groups = containers.Map('KeyType', 'double', 'ValueType', 'any');
      for i = 1:numel(dist_idx)
          root = find_root(i);
          if ~groups.isKey(root)
              groups(root) = [];
          end
          groups(root) = [groups(root), dist_idx(i)];
      end
      
      % 更新 theta_list
      group_keys = cell2mat(groups.keys);
      for k = 1:numel(group_keys)
          g_idx = groups(group_keys(k));
          if numel(g_idx) > 1
              % 计算聚合后的 L（所有线元长度之和）
              L_total = 0;
              for i = 1:numel(g_idx)
                  L_total = L_total + theta_list(g_idx(i)).L;
              end
              % 更新组内所有原子的 L
              for i = 1:numel(g_idx)
                  theta_list(g_idx(i)).L = L_total;
              end
          end
      end
      
      if cfg.verbose
          fprintf('[AGG] aggregated %d distributed atoms into %d groups\n', ...
              numel(dist_idx), numel(group_keys));
      end
  end

  %% 细化 L、alpha、gamma 参数
  function theta_list = v6_refine_L_alpha_gamma(theta_list, I, roi_mask, cfg)
      N = numel(theta_list);
      if N == 0, return; end
      
      for i = 1:N
          % 通过 L 判断类型
          if theta_list(i).L > 0
              % 分布式散射中心：gamma = 0
              theta_list(i).gamma = cfg.gamma_distributed;
              % L 已在聚合中确定
          else
              % 局域式散射中心：L = 0，估计 gamma > 0
              theta_list(i).gamma = v6_fit_localized_gamma(theta_list(i), I, roi_mask, cfg);
          end
          
          % 估计 alpha（散射强度）
          theta_list(i).alpha = theta_list(i).A; % 简化：使用幅度作为 alpha
      end
      
      if cfg.verbose
          n_loc = sum([theta_list.L] == 0);
          n_dist = sum([theta_list.L] > 0);
          fprintf('[REFINE] refined %d localized (L=0, gamma>0) and %d distributed (L>0, gamma=0)\n', ...
              n_loc, n_dist);
      end
  end

  %% 拟合局域式散射中心的 gamma
  function gamma = v6_fit_localized_gamma(theta, I, roi_mask, cfg)
      % 对局域式散射中心拟合径向衰减模型 I(r) = A * exp(-gamma * r)
      r = theta.row; c = theta.col;
      [H, W] = size(I);
      
      % 提取局部区域
      half_win = min(10, floor(cfg.patch_size / 2));
      r1 = max(1, r - half_win); r2 = min(H, r + half_win);
      c1 = max(1, c - half_win); c2 = min(W, c + half_win);
      
      I_patch = I(r1:r2, c1:c2);
      mask_patch = roi_mask(r1:r2, c1:c2);
      
      % 计算径向距离
      [rows, cols] = size(I_patch);
      [xx, yy] = meshgrid(1:cols, 1:rows);
      cx = c - c1 + 1; cy = r - r1 + 1;
      r_vals = sqrt((xx - cx).^2 + (yy - cy).^2);
      
      % 提取有效像素
      valid = mask_patch & (I_patch > 0) & (r_vals > 0.5);
      if sum(valid(:)) < 5
          gamma = 0.5; % 默认值
          return;
      end
      
      r_data = r_vals(valid);
      I_data = I_patch(valid);
      
      % 加权最小二乘拟合：log(I) = log(A) - gamma * r
      % 权重：距离中心越近权重越大
      weights = exp(-r_data / 3);
      
      % 构建线性系统：[1, r] * [log(A); -gamma] = log(I)
      X = [ones(size(r_data)), r_data];
      y = log(I_data + eps);
      W = diag(weights);
      
      % 加权最小二乘
      try
          beta = (X' * W * X) \ (X' * W * y);
          gamma_est = -beta(2);
      catch
          gamma_est = 0.5;
      end
      
      % 限制在合理范围内
      gamma = max(cfg.gamma_localized_range(1), min(cfg.gamma_localized_range(2), gamma_est));
  end