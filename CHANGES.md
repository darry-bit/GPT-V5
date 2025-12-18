# Changes to ASC-NOMP V6

## Summary

This update improves the scattering center extraction algorithm by implementing proper type detection based on physical parameters and improving gamma estimation.

## Key Changes

### 1. Scattering Center Type Detection (Based on L Parameter)

**Previous:** Type determined by dictionary's 'class' field ('localized' or 'distributed')

**New:** Type determined by L parameter value:
- **L = 0**: Localized scattering center (point target)
- **L > 0**: Distributed scattering center (extended target)

### 2. Multiple Localized Centers Support

**Previous:** `cfg.max_localized = 1` limited extraction to single localized center

**New:** `cfg.max_localized = Inf` allows multiple localized scattering centers

### 3. Gamma Estimation by Type

**Previous:** No distinction between scattering center types for gamma estimation

**New:** Type-specific gamma estimation:
- **Distributed (L > 0)**: `gamma = 0` (no attenuation, uniform energy distribution)
- **Localized (L = 0)**: `gamma > 0` (estimated via radial decay fitting)

### 4. Improved Parameter Estimation Order

**Correct order implemented:**
1. OMP matching determines position and initial L
2. Aggregation determines final L value (for distributed centers)
3. Type determination based on L
4. Gamma estimation based on type
5. Alpha (scattering intensity) estimation

## New Functions

### `v6_aggregate_atoms(theta_list, cfg)`
- Implements Union-Find algorithm for chain-based aggregation
- Groups nearby distributed atoms with similar orientations
- Considers connection direction consistency with scattering direction
- Relaxed thresholds: `dist_thresh=6`, `angle_thresh=30°`

### `v6_refine_L_alpha_gamma(theta_list, I, roi_mask, cfg)`
- Determines scattering center type based on L value
- Sets gamma appropriately for each type
- Estimates alpha (scattering intensity)

### `v6_fit_localized_gamma(theta, I, roi_mask, cfg)`
- Fits radial decay model: `I(r) = A * exp(-gamma * r)`
- Uses weighted least squares fitting
- Constrains gamma to range `[0.1, 2.0]`

## Configuration Changes

### New Parameters

```matlab
cfg.max_localized = Inf;  % Allow multiple localized centers
cfg.gamma_localized_range = [0.1, 2.0];  % Valid range for localized gamma
cfg.gamma_distributed = 0;  % Fixed gamma for distributed centers
cfg.agg_dist_thresh = 6;  % Aggregation distance threshold
cfg.agg_angle_thresh = 30;  % Aggregation angle threshold (degrees)
```

## Physical Model Reference

| Scattering Type | L | gamma | Physical Meaning |
|----------------|---|-------|------------------|
| Distributed | L > 0 | γ = 0 | Energy uniformly distributed along length, no attenuation |
| Localized | L = 0 | γ > 0 | Energy concentrated at point, radial attenuation from center |

## Testing

Run `demo_sc_extraction_v6_test.m` to test the implementation:

```matlab
demo_sc_extraction_v6_test()
```

### Test Coverage

1. Distributed atoms have L > 0
2. Distributed atoms have gamma = 0
3. Localized atoms have gamma > 0
4. All alpha values are valid (> 0)
5. Successful atom extraction

### Visualization

The demo script creates visualizations showing:
- Original image
- Extracted scattering centers (color-coded by type)
  - Red circles: Localized (L=0)
  - Blue squares: Distributed (L>0)
- Direction indicators for distributed centers
- ROI and energy statistics

## Compatibility

All changes are backward compatible. Existing code will continue to work with default parameters.
