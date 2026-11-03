# Photometry pipeline：后续 notebook tutorial 备忘录

更新日期：2026-11-03。

本备忘录记录当前已经实际跑通的 SARA 光学 aperture photometry 流程，
供后续编写正式 notebook tutorial 使用。从已完成 plate solving、FWHM
测量的 Master Light 开始，结束于 science targets 的 calibrated magnitudes。
本次分析已经运行到最终定标结果，并由用户确认通过。

详细 API 行为见 [Photometry.md](Photometry.md)；source/result 数据模型见
[Sources.md](Sources.md)。本文件按照实际运行顺序组织，每段代码可作为一个 cell。

## 1. 前置对象与 imports

本节接在之前的处理流程后面，不重新执行 reduction、registration、plate solving
或 FWHM 测量。前面的 notebook 应明确创建下列对象，避免依赖 kernel 残留变量：

| 对象 | 用途 |
| --- | --- |
| `analysis_dir: Path` | 当前晚上的 analysis 根目录 |
| `analysis_images` | 已有的 `ImageCollectionGroup`，本例使用第一个 child collection |
| `DLS_field: ObservationField` | 标准场，包含 field name 和 Telescope |
| `sources: Sources` | 已定义的 science targets，各自具有坐标和 observation field |

实际分析 notebook 为 `../20240215_has_perfect_flats/Photometry_Pipeline.ipynb`
（相对于仓库根目录），属于分析数据目录，不随本仓库提交。

Master Light 应具有 WCS、`FIELD`、`FILTER` 和以像素为单位的 `FWHM`。
测光还读取 `BUNIT`（ADU）、`GAIN`、`RDNOISE`、`NCOMBINE`。
这些 metadata 来自之前的处理步骤，不在这里重复手工填写。

```python
from pathlib import Path

from photozpy.photometry import PhotometryRegions
from photozpy.photometry.photometry import Photometry
from photozpy.standard_stars import StandardStarCatalog
```

## 2. 从标准场选择一张 catalog query 图像

`ObservationField` 不保存 FITS path。用 field name 和一个 reference filter
筛选已有 collection；这张图像仅用于 catalog footprint query 和初始位置 preview。
后面的 regions 仍按所有可用 filter 分别生成。

```python
# Select the observation field containing the standard stars.
standard_star_field = DLS_field

# Choose a reference filter from the field's Telescope configuration.
standard_reference_filter = standard_star_field.telescope.filters()[-1]

# Get this field/filter's Master Light from the existing analysis collection.
standard_star_path = Path(
    analysis_images[0].files_filtered(
        imtype="Master Light",
        field=standard_star_field.field_name,
        filter=standard_reference_filter,
        include_path=True,
    )[0]
)
```

如需指定 reference filter，把对应变量改为
`standard_star_field.telescope.resolve_filter_alias("SDSS g")` 即可。
本例假定所选 field/filter 对应一张 Master Light，不需要手工查文件名。

## 3. 读取 standard catalog，创建标准星集合

```python
# Locate the compiled standard-star catalog used for this analysis.
catalog_path = Path(
    "/Users/yong/Desktop/Research/Photozpy/"
    "Compiled_Standard_Stars/Smith_Southern_Sky.h5"
)

# Initialize the catalog reader and inspect its metadata.
standard_star = StandardStarCatalog(catalog_path=catalog_path)
standard_star.print_metadata()

# Query catalog stars in the image and attach their reference magnitudes.
standard_sources = standard_star.find_standard_sources(
    fits_path=standard_star_path,
    observation_field=standard_star_field,
    magnitude_system="AB",
    edge_margin_pix=100.0,
    mag_limit=23.0,
)

# Inspect all matched catalog stars and measurements in one DataFrame.
standard_star.matched_stars
```

正式 tutorial 需要说明：

- 这里关联 catalog 坐标和 WCS footprint，还不是图像中的自动 source detection。
- 每颗星创建一个 `SourceInfo`，catalog 的 `id_name` 是稳定 source ID。
- 各星的 catalog magnitudes/errors 写入 `Sources` 的 `reference` 槽；
  此时 `instrumental` 和 `calibrated` 仍为 `None`。
- filter/magnitude/error 列通过 field 的 Telescope aliases 自动关联，
  正常使用不需要逐个填写 mapping dictionary。
- `mag_limit=23.0` 要求所有自动匹配的 catalog bands 满足限制。
  SARA 的 aliases 可匹配 u/g/r/i/z，因此 catalog 有 u 时也参与筛选，
  即使本晚没有 u 图像。本次 DLS 样例返回 6 颗候选标准星。
- `magnitude_system="AB"` 来自本例 catalog 的 metadata，不执行 AB/Vega
  或 passband conversion。
- `matched_stars` 是经过所有筛选后的表格副本。修改它不会修改已创建 source
  的 reference magnitudes。

## 4. 用 interactive DS9 预览 catalog 位置

```python
# Preview matched stars at their initial catalog coordinates in DS9.
standard_viewer = standard_star.preview_standard_stars(
    fits_path=standard_star_path,
    sources=standard_sources,
)
```

蓝色固定 cross 表示初始 catalog 位置。这里仅 preview；可编辑的 aperture、
annulus 以及接受/删除决定在下一步处理。正常 pipeline 不用旧的 image plot。

## 5. Centroid 并生成标准场所有 filter 的 regions

```python
# Generate image-specific regions for all standards and available filters.
standard_regions = PhotometryRegions.generate_regions(
    image_collection=analysis_images,
    sources=standard_sources,
    root=analysis_dir,
    aperture_scaling=3.0,
    inner_annulus_scaling=5.0,
    outer_annulus_scaling=8.0,
    box_size=19,
)

# Review and save the standard-star regions in interactive DS9.
standard_regions.review()
```

每张图像使用自己的 header FWHM 和 centroid。默认 aperture radius 为
`3 × FWHM`，annulus 内外半径为 `5 × FWHM`、`8 × FWHM`；三个 scaling
factors 可修改。`box_size` 是 centroiding 的局部像素窗口宽度，不是 aperture 半径。

旧 centroiding 主线保留：背景估计、两次 COM，然后 quadratic centroid；
失败时尝试 COM，再回退初始 WCS 坐标并记录原因。原始 catalog/source sky
coordinates 不会被手工 region 调整覆盖。

## 6. 批量处理全部 science targets

```python
# Generate per-source, per-filter regions for the science-source collection.
science_regions = PhotometryRegions.generate_regions(
    image_collection=analysis_images,
    sources=sources,
    root=analysis_dir,
    aperture_scaling=3.0,
    inner_annulus_scaling=5.0,
    outer_annulus_scaling=8.0,
    box_size=25,
)

# Review every science target and filter in one queue.
science_regions.review()
```

不需要逐颗调用生成函数。同一 target 的 filters 排在一起；切换同一 source/group
的 filter 时继承 zoom 和 image-coordinate view center，换 target/field 时
恢复 zoom-to-fit。

### DS9 review 操作约定

- 绿色 circle 是 aperture，红色 annulus 是 background，蓝色固定 cross 是初始
  坐标。仅 aperture 显示 source name，避免双重标签。
- 手工移动或修改后点击 **Save & Next**，保存并覆写对应文件；最后一张同样
  要点击此按钮保存。
- aperture 与 annulus 可以不同心，以避开附近污染源。偏心或重叠会显示红色
  warning，但不强行移回原位，也不阻止保存。
- 不采用的标准星同时删除 aperture 和 annulus；只对当前图像生效，不删除
  其他 filter 的记录或 catalog reference。
- 保留稳定 ID/role tags。只删一个 shape、破坏配对或非法半径会阻止保存。
- **Previous** 读取已保存状态；**Quit** 不保存当前尚未确认的修改。
- 靠近边缘、stack 覆盖不足或污染严重的标准星，在这里决定是否保留。
  几何 coverage diagnostics 不代表每个像素实际贡献的曝光张数。

### 文件组织和坐标记录

```text
analysis/
    photometry/regions/
        DLS_field_standards/
            SDSS_g'_aperture.reg
            SDSS_g'_background.reg
            SDSS_g'_regions.json
            SDSS_r'_aperture.reg
            SDSS_r'_background.reg
            SDSS_r'_regions.json
            ...
        4FGL_J0746.6-4754/
            fgl_field_SDSS_g'_aperture.reg
            fgl_field_SDSS_g'_background.reg
            fgl_field_SDSS_g'_regions.json
            ...
```

标准星每个 field/filter 共用一组文件；science targets 每个 source/field/filter
一组。不同 filter 放在同一 source/group 文件夹下，不再各建一层目录。
注释记录 field、source、filter 和对应 image。图像/region 关联使用相对于 analysis
根目录的 paths；移动整套 analysis 时保留目录结构及内部记录文件。
JSON 由程序管理，正常调用无需手工定位。

`SourceInfo` 保存身份和 sky coordinate，不保存 region paths。最终像素位置和
文件关联保存在 `PhotometryRegions` records、测光结果表中；各阶段 magnitudes
保存在 `Sources` 的结果槽中。

## 7. 统一做 aperture photometry

```python
# Use the existing image selection for this night's measurements.
photometry = Photometry(image_collection=analysis_images)

# Measure standards and science targets together with their reviewed regions.
photometry_table = photometry.run_photometry(
    sources=(standard_sources, sources),
    regions=(standard_regions, science_regions),
    bkg_clip_sigma=3.0,
    src_detection_sigma=3.0,
    verbose=False,
)

# Inspect counts, instrumental magnitudes, significance, and measurement status.
photometry_table
```

这个入口统一处理两类 source，也支持它们位于同一 science frame。
测光使用最终保存的 aperture/annulus，包括手工中心偏移和半径，
不会重新 centroid 或再次乘 FWHM scaling。

这里只更新 `instrumental` 槽，保留 reference/calibrated 槽。如果修改 regions
后重新测光，应重新执行后面的零点计算和定标；旧 calibrated 结果不会自动刷新。

### 科学定义与保留的误差方法

1. 背景沿用旧的 sigma-clipped annulus **mean**，不自动改成 median。
2. source aperture 做 exact sum，再扣除平均背景乘 aperture area。
3. `MASK` 和非有限像素无效；`CMASK` 只记录 cosmic-ray 检测，修补后的数值
   仍参加测光。
4. source aperture 有无效/越界像素时报告 `invalid_aperture`，不补估缺失 flux；
   annulus 无效像素排除，沿用旧的有效面积定义。
5. 误差调用原来的 `estimate_source_region_error`，手动使用 gain、read noise、
   aperture/annulus 面积及 `NCOMBINE`；不用自动生成的 CCDData uncertainty。
6. 当前使用 average-stack **ADU counts**，不是 count rate，不在这里除以曝光
   时间。本轮保留原噪声模型，不新增配准 covariance 或 calibration-frame
   uncertainty 传播项。

对 detection：

$$
m_{\mathrm{inst}}=-2.5\log_{10}(s),\qquad
\sigma_m=\frac{2.5}{\ln 10}\frac{\sigma_s}{s},\qquad
\mathrm{significance}=\frac{s}{\sigma_s}.
$$

significance 小于阈值或 net counts 非正时，沿用旧 background-only error
计算 limiting magnitude，默认阈值为 3。测光 QTable 的 magnitude error 保留
NaN，写入 `FilterMagnitude` 时用现有 `-99` error flag 标识 upper limit；
这个 flag 不是物理误差。背景像素 sigma clipping 与后面的零点 z-score 筛选
属于不同步骤，正式 notebook 不要混为一谈。

## 8. 在 standard sources 内计算各 filter 的 zero point

```python
# Calculate zero points using the original single-pass z-score clipping recipe.
zero_point_table = standard_sources.calculate_zero_points(sigma=3.0)

# Inspect the per-filter mean, scatter, and retained/clipped sample counts.
zero_point_table
```

```python
# Inspect every standard/filter offset and its inclusion or exclusion reason.
standard_zero_point_table = standard_sources.zero_point_details

# Display the complete standard-star zero-point overview.
standard_zero_point_table
```

每个 filter 独立使用有效、非 upper-limit、significance ≥ 3 的标准星：

$$
ZP_i=m_{\mathrm{ref},i}-m_{\mathrm{inst},i},\qquad
ZP_{\mathrm{filter}}=\operatorname{mean}(ZP_{i,\mathrm{retained}}).
$$

保留旧的一次 `scipy.stats.zscore` 筛选和不加权平均，**不是迭代 sigma clipping**。
样本少时 3σ 可能保留视觉上偏离的点。单颗有效星或所有零点相同的情况会保留
样本，避免 undefined z-score 全部排除。

summary 的 `scatter` 是保留样本的 population standard deviation，不是
zero-point uncertainty。无有效标准星的 filter 为 NaN。`n_clipped` 只计算
z-score 剔除，不包含 significance 等条件的排除。

## 9. 查看 zero point statistics

```python
# Plot all finite offsets and distinguish retained and rejected measurements.
# Red lines are average +/- 0.1 mag visual guides, not uncertainty bounds.
zero_point_figure = standard_sources.plot_zero_point_statistics()
```

每个有有效零点的 filter 一个 panel：蓝色圆点为 retained，橙色叉号为 rejected，
绿色点为平均零点。标题明确显示数量；无法计算零点的记录另列 unavailable。
所有有效 rejected 点也纳入坐标范围，包括因 significance 或 clipping 未参与
平均的点。默认只显示图，不自动写文件。

图用于检查标准星定标的一致性。scatter 大不能单独确定是天气造成的，
还需结合污染、catalog 和测光情况判断。

## 10. 统一定标 science targets，显示紧凑摘要

```python
# Apply the inspected standard-star zero points to all science targets.
# Preserve instrumental errors, significance, and upper-limit flags.
calibrated_table = sources.calibrate_targets(
    standard_sources=standard_sources,
)

# Display one filter/magnitude/error/significance summary per source.
sources.show_magnitudes(kind="calibrated")
```

最终使用：

$$
m_{\mathrm{cal}}=m_{\mathrm{inst}}+ZP_{\mathrm{filter}},\qquad
\sigma_{m,\mathrm{cal}}=\sigma_{m,\mathrm{inst}}.
$$

本轮有意**不合并 catalog error 或 zero point uncertainty**，这两个误差项
留作增强功能。输出系统沿用 reference 的 AB/Vega 标签，不隐式转换。

定标复用已经检查过的零点，更新 science sources 的 `calibrated` 槽；
重复调用从 instrumental 开始，不会重复加零点。标准星结果变化后要明确重新
计算零点，程序不会静默继续用旧值。如果未算过零点，`calibrate_targets` 会
先按默认参数计算；正式 notebook 仍推荐保留前面的计算与检查 cell。

`show_magnitudes` 按 source 显示名字、field、system 和四列紧凑结果，复用
`FilterMagnitude.to_table()`。缺失值为 `--`，upper-limit error 为 `UL (-99)`。
这个方法只展示，不修改或重新计算。完整 DataFrame 保留在 `calibrated_table`，
供详细检查和导出。

```python
# Inspect instrumental measurements with the same layout if needed.
sources.show_magnitudes(kind="instrumental")

# Inspect catalog reference magnitudes for the standards if needed.
standard_sources.show_magnitudes(kind="reference")
```

## 11. Kernel 重启后读取已有 regions

主流程目标是一次完成整晚分析，不要求配置复杂的断点恢复。
如果 kernel 重启，先重新定义 fields、sources、collection，并通过 catalog
query 恢复 `standard_sources` 及 reference magnitudes。然后用下面两段**替代**
`generate_regions`，不要为了恢复对象覆盖手工调整。

```python
# Load saved standard regions and review decisions without centroiding again.
standard_regions = PhotometryRegions.load_regions(
    image_collection=analysis_images,
    sources=standard_sources,
    root=analysis_dir,
)

# Load saved science regions, retaining all manually edited geometry.
science_regions = PhotometryRegions.load_regions(
    image_collection=analysis_images,
    sources=sources,
    root=analysis_dir,
)
```

随后重跑 photometry、zero points 和 calibration。此入口只恢复 region 关联，
不序列化整个 `Sources`/`ObservationField`，也不恢复 kernel 的数值结果。
缺文件会明确报错，不自动生成替代文件，不需要手工查 JSON。
对象还在时，可直接再次调用 `standard_regions.review()` 或 `science_regions.review()`。

## 12. 正式 notebook 编写与验收记录

- 主流程拆分成短 cell，每个小步骤保留英文 comment。
- 用 `matched_stars` 总览 candidates；interactive DS9 检查图像和 regions；
  statistics plot 检查零点；`show_magnitudes` 展示最终结果。
- 当前实际样例：6 颗标准星 × 4 filters，加 2 个 science targets × 4 filters，
  共 32 项测光。没有 u-band 测量时，该结果槽仍为 NaN。
- 用户已确认 region 手工调整保存、review 的视图继承正常，并取得最终
  calibrated photometry。
- 本轮自动测试共 500 项，覆盖 catalog/result 关联、region/DS9 往返、原误差
  公式、unified photometry、zero points、upper limits、定标和展示。

后续可以独立讨论 catalog/zero-point uncertainty、aperture correction、自动
aperture 优化、迭代 clipping 或加权平均，以及 science FOV 内自动搜索标准星
的 notebook 入口。当前统一测光/定标已允许两类 source 共处一个 field；
这些增强项不作为本次流程的完成条件。
