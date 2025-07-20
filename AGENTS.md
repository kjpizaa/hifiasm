### hifiasm有参考基因组组装 - AI开发助手快速指南

> **受众:** GitHub Copilot、ChatGPT、Claude、Cursor AI等AI编程助手  
> **人类开发者:** 详细设计文档请参考 `docs/ref_guided_assembly_design.md`

-----
## 8 · Quick-API cheat-sheet (Keep these stable!)

| 分类 | 核心函数 / 结构 | 所在文件 (commit ec9a8b) | 你改动时要注意 |
|------|----------------|-------------------------|----------------|
| **参考弧注入** | `load_paths_and_inject(asg_t *g, const char *tsv, int64_t span_bp, int check_hifi)` | `src/ref_arc.c` | 解析 `paths.tsv` → 对称地加 `ARC_FLAG_REF` 弧；**勿在此调用 `asg_cleanup()`** |
| | `asg_path_within(const asg_t *g, uint32_t v0, uint32_t v1, int64_t max_bp, int max_step)` | `src/sg_utils.c` | 小 BFS 判 HiFi≤6 kb 是否已连通；忽略 `ARC_FLAG_REF` |
| | `ARC_FLAG_REF` (`1u<<30`)、`ARC_REF_W` (默认 400) | `Overlaps.h` | 高位标志唯一；权重越大参考弧越少被剪 |
| **Overlap → String-graph** | `ma_hit2arc()` | `Overlaps.h` §366-442 / `Overlaps.cpp` 896-972 | 把 `ma_hit_t` 转换成对称弧；保持 rev/方向一致性 |
| | `normalize_ma_hit_t_single_side_advance[_mult]` | `Overlaps.cpp` 1139-1586 | 确保 A↔B 都有 overlap；别误删 `rev` 位 |
| | `ma_hit_contained_advance` | `Overlaps.cpp` 1782-1862 | contained read 处理；参考弧对 contained 逻辑无影响 |
| | `ma_hit_flt` | `Overlaps.cpp` 1865-1922 | 改阈值时同步更新 `--ref-path` 教程中的示例 |
| **基础 graph 清洗** | `asg_arc_del_multi` | `Overlaps.cpp` 980-1005 | 保留最长弧；**切勿删除 `ARC_FLAG_REF`** |
| | `asg_arc_del_asymm` | `Overlaps.cpp` 1007-1035 | 保证弧对称；同上 |
| | `asg_cleanup()` | 多处调用 | 注入参考弧后仅需调用一次 |
| **高级 graph 算法** | `asg_bub_pop1_primary_trio`, `asg_pop_bubble_primary_trio` | `Overlaps.cpp` 7963-8146 / 11054-11200 | 泡 >50 kb 时**不要**自动替换为参考弧路径 |
| | `check_tip`, `asg_arc_del_simple_circle_untig` | `Overlaps.cpp` 3215-3290 / 3163-3214 | 在删除 tip 时：`if(e->ul & ARC_FLAG_REF) continue;` |
| | `resolve_tangles` | `Overlaps.cpp` 11546-11725 | 参考弧参与评分，但不能压倒 HiFi 证据 |
| **Unitig 生成 / 遍历** | `ma_ug_gen()` | `Overlaps.cpp` 36310-36479 | 压缩 1-进-1-出链；`ARC_FLAG_REF` 只影响弧权，不改链判定 |
| | `ma_utg_t` / `ma_ug_t` | `Overlaps.h` 214, 266-271 | `a[j] = (vertex<<32 | edgeLen)` 格式不能破 |
| **UL 复用管线** | `hifi_unitigs_map_to_reference()` | `inter.cpp` (新增) | 把参考染色体当 “虚拟 UL”；产生 `uc_block_t` 带 `BLOCK_REF` |
| | `ul_resolve()` / `ul_refine_alignment()` | `inter.cpp` 881-960 等 | 签名别改；只读 `BLOCK_REF` |
| | `extend_coordinates()` | `inter.cpp` 912-960 | 需支持 800 Mb 染色体坐标，避免 `uint32_t` 溢出 |
| **Overlap 检测内部结构** | `ma_hit_t` / `ma_hit_t_alloc` | `Overlaps.h` 40-45 / 116-132 | `el==1` 精确；`del` 删边标记；参考弧注入不触碰这些 |
| | `overlap_region`、`window_list` | `Hash_Table.h` 39-106 / 160-189 | 调 k-mer 长度时同步修改测试脚本 |
| **性能 / 调试** | `radix_sort_arch64` | `Overlaps.cpp` 24-45 | 新字段别破 64-bit 键格式 `(weight<<32|id)` |
| | `fprintf(stderr,"[REF-ARC] ...")` 宏 | `src/ref_arc.c` | 标准 debug tag，方便 grep |



## 10 · “一层参考 vs. 多层真 UL” - 权重与清洗隔离

| 组件 | 对 **真 UL / HiFi** 的行为 | 对 **参考虚拟 UL** 的特殊处理 | 所在文件 / 常量 |
|------|---------------------------|--------------------------------|-----------------|
| **Overlap→Arc 权重** | `weight = ol` (overlap 长) × coverage | 固定 `ARC_REF_W ≈ 400`<br>不乘深度 | `Overlaps.h` |
| **Unitig-occ 统计** | `ug_occ_w()` 按弧权 * 深度叠加 | 直接写 `UL_REF_WEIGHT = 900`<br>(只一条) | `gfa_ut.h:ug_occ_w` |
| **深度投票解歧 (`ul_resolve`)** | 多条 UL 比票数 | 参考 block 仅做端点合法性校验<br>**不**参加投票 | `inter.h:ul_resolve` |
| **Bubble / Tip 剪枝** | 可删除弱弧 | `if(e->ul & ARC_FLAG_REF) continue;`<br>参考弧豁免 | 多处 *pop/clip* 例外分支 |
| **BFS 封顶距离** | `max_bp = 6 kb` | 参考弧用 `REF_MAX_SPAN (50 Mb)` | `sg_utils.c:asg_path_within` |
| **坐标类型** | 大多 `uint32_t` | 参考 chr 800 Mb → **全部 `uint64_t`** | `inter.cpp:extend_coordinates` |

### 调参须知  
* **调弱参考** → 降 `ARC_REF_W` 或 `UL_REF_WEIGHT`  
* **调强参考** → 升上述常量，或放宽 `REF_MAX_SPAN`  
* 任何时候 **保持弧对称**：加 `v→w` 弧必须同步 `w^1→v^1`.

当AI需要基于这个设计实现代码时，应该：

1. **优先理解现有UL流程**：重点关注uc_block_t数据结构和UL_INF全局管理
2. **严格遵循标记系统**：正确使用BLOCK_REF标记
3. **复用现有函数**：尽量调用现有比对和处理函数
4. **注重错误处理**：每个步骤都要有失败回退机制
5. **保持接口兼容**：不修改任何现有函数签名
6. **🔧 注意类型转换**：正确处理`ref_ul_idx_t`与`ul_idx_t`的转换
7. **🔧 避免内存泄漏**：参考基因组索引使用引用模式，避免浅拷贝导致double free

**特别注意编译修复**：
- 使用`ref_ul_idx_t`而非`ul_idx_t`
- 正确的HPC开关：`!(asm_opt.flag & HA_F_NO_HPC)`
- 每条染色体作为独立unitig处理
- 内存优化选项：支持可选的merged_seq释放

这样可以确保实现既功能完整又稳定可靠，同时避免编译错误和运行时问题。# AGENTS.md - HiFiasm参考基因组增强实现指南

## 1. 项目概述

本项目是hifiasm的增强版本，在Ultra-Long (UL) 组装流程中集成参考基因组支持。核心理念是将参考基因组转换为**"虚拟UL读段"**，通过现有UL处理基础设施实现参考指导的组装优化。

### 核心设计原则 🎯
- **98%代码复用**：最大化利用现有UL处理流程
- **最小侵入性**：仅在关键位置添加功能，不修改现有函数
- **向后兼容**：不传入`--ref-fasta`时行为完全不变
- **优雅降级**：参考基因组处理失败时自动回退到标准hifiasm

---

## 2. 关键文件和修改点

### 粘贴内容分析：文件位置 `inter.cpp`

从粘贴的代码可以看出，这是在`ul_realignment`函数中添加的参考基因组集成功能：

```cpp
ma_ug_t *ul_realignment(const ug_opt_t *uopt, asg_t *sg, uint32_t double_check_cache, const char *bin_file)
{
    // ... 原有初始化代码 ...
    ma_ug_t *ug = gen_polished_ug(uopt, sg);
    
#ifdef ENABLE_REF_GENOME_V4
    // 🔧 关键集成点：在UL管线内直接集成参考基因组
    extern ref_genome_t *global_ref_genome;
    extern hifiasm_opt_t asm_opt;

    if (global_ref_genome && asm_opt.ref_fasta && asm_opt.ref_fasta[0] && ug && ug->u.n > 0) {
        fprintf(stderr, "[M::%s] Integrating reference genome into UL pipeline\n", __func__);

        int result = integrate_reference_blocks_to_existing_ul_pipeline(
            ug,
            (const ul_idx_t*)global_ref_genome->ul_index,
            &asm_opt
        );

        if (result == 0) {
            fprintf(stderr, "[M::%s] Reference integration successful\n", __func__);
        } else {
            fprintf(stderr, "[WARNING] Reference integration failed, continuing with standard UL\n");
        }
    }
#endif
    
    // ... 原有UL处理流程继续 ...
}
```

---

## 3. 核心技术架构

### 3.1 数据流转换链路

```
参考基因组FASTA → ref_genome_t → ul_idx_t → unitig比对 → overlap_region → uc_block_t → UL_INF全局结构
```

### 3.2 关键数据结构映射

| 原始数据 | 转换目标 | 作用 |
|---------|---------|------|
| `染色体序列` | `All_reads.seq[i]` | 虚拟UL读段 |
| `参考基因组索引` | `ul_idx_t` | 支持现有比对流程 |
| `unitig-参考比对结果` | `uc_block_t` | 统一数据格式 |
| `uc_block_t.el高位` | `BLOCK_REF标记` | 区分参考基因组blocks |

### 3.3 集成时机和位置

**最佳集成点**：`ul_realignment`函数内，`gen_polished_ug`调用之后
- **原因**：此时unitig图已构建完成，可以进行参考比对
- **优势**：不影响原有HiFi数据处理流程
- **安全性**：条件编译保护，失败时优雅降级

---

## 4. 关键函数解析

### 4.1 核心新增函数（来自粘贴内容）

#### `ensure_unitig_seq(ma_ug_t* ug, uint32_t uid)`
- **作用**：确保unitig序列可用，处理序列缺失情况
- **返回**：unitig序列指针或NULL

#### `overlap_to_uc_block_ref_mode()`
- **作用**：将overlap结果转换为uc_block_t格式
- **特点**：设置`BLOCK_REF`标记，区分参考基因组blocks
- **质量控制**：过滤短比对（<500bp）和高错误率比对（>15%）

#### `unitigs_map_to_reference_batch()`
- **作用**：批量比对所有unitigs到参考基因组
- **核心逻辑**：
  ```cpp
  // 对每个unitig
  for (uint32_t uid = 0; uid < unitigs->u.n; uid++) {
      // 调用现有比对函数
      ha_get_ul_candidates_interface();
      // 转换结果为uc_block_t
      overlap_to_uc_block_ref_mode();
  }
  ```

#### `integrate_reference_blocks_to_existing_ul_pipeline()`
- **作用**：主集成函数，将参考blocks注入UL处理流程
- **关键步骤**：
  1. 生成参考基因组blocks
  2. 扩展`UL_INF`全局结构容量
  3. 添加blocks到对应unitig的block列表
  4. 调用现有UL处理函数完成集成

### 4.2 现有函数复用

从项目知识可以看出，以下现有函数被完全复用：
- `filter_ul_ug(ug)` - UL unitig过滤
- `gen_ul_vec_rid_t(&UL_INF, NULL, ug)` - 生成UL向量索引
- `update_ug_arch_ul_mul(ug)` - 更新unitig架构

---

## 5. 数据结构标记系统

### 5.1 Block标记机制
```cpp
#define BLOCK_REF             (1u<<15)
#define BLOCK_SET_REF(block)   ((block)->el |= BLOCK_REF)
#define BLOCK_IS_REF(block)    ((block)->el & BLOCK_REF)
#define BLOCK_CLEAR_REF(block) ((block)->el &= ~BLOCK_REF)
```

### 5.2 标记的作用
- **区分数据源**：参考基因组blocks vs 真实UL数据
- **处理差异化**：参考blocks在某些清洗步骤中豁免删除
- **调试支持**：便于追踪参考基因组数据流

---

## 6. 集成流程详解

### 6.1 集成触发条件
```cpp
if (global_ref_genome && asm_opt.ref_fasta && asm_opt.ref_fasta[0] && ug && ug->u.n > 0)
```
- 全局参考基因组已初始化
- CLI指定了参考基因组文件
- unitig图构建成功

### 6.2 集成执行流程

1. **准备阶段**：检查参数有效性
2. **比对阶段**：`unitigs_map_to_reference_batch()` - 批量比对
3. **转换阶段**：overlap_region → uc_block_t转换
4. **注入阶段**：扩展UL_INF结构，添加参考blocks
5. **处理阶段**：调用现有UL函数完成处理

### 6.3 错误处理策略
- **比对失败**：跳过该unitig，继续处理其他
- **内存不足**：释放已分配资源，返回错误
- **集成失败**：输出警告，继续标准UL流程

---

## 7. 内存管理和资源清理

### 7.1 内存分配模式
```cpp
// 两阶段分配：临时存储 + 最终合并
uc_block_t **all_results = (uc_block_t**)calloc(unitigs->u.n, sizeof(uc_block_t*));
uint64_t *all_counts = (uint64_t*)calloc(unitigs->u.n, sizeof(uint64_t));

// 最终合并为单一数组
uc_block_t *final_blocks = (uc_block_t*)malloc(total_blocks * sizeof(uc_block_t));
```

### 7.2 清理责任链
1. **overlap临时结构**：每次比对后立即清理
2. **临时block数组**：合并后释放
3. **最终block数组**：集成完成后释放
4. **全局UL_INF**：由现有流程管理

---

## 8. 与现有UL流程的协同

### 8.1 完美匹配的设计
- **数据格式兼容**：uc_block_t是UL流程的标准格式
- **标记透明传递**：BLOCK_REF标记在整个流程中保持
- **函数接口不变**：现有UL函数无需修改

### 8.2 处理流程融合
```cpp
// 参考blocks注入后，标准UL流程正常运行：
filter_ul_ug(ug);                    // 过滤（参考blocks豁免）
gen_ul_vec_rid_t(&UL_INF, NULL, ug); // 生成索引
update_ug_arch_ul_mul(ug);          // 更新架构
```

---

## 9. 调试和验证要点

### 9.1 关键日志信息
```cpp
"[M::%s] Processing %zu unitigs → reference mapping"
"[M::%s] Generated %lu reference blocks"  
"[M::%s] Added %lu reference blocks to UL_INF"
"[M::%s] Reference-guided UL pipeline completed"
```

### 9.2 验证检查点
- 参考blocks数量是否合理
- BLOCK_REF标记是否正确设置
- UL_INF结构是否正确扩展
- 现有UL函数是否正常执行

---

## 10. 修改原则和约束

### 10.1 代码修改约束
- **禁止修改现有函数签名**：保持向后兼容
- **最小化全局变量**：仅使用已有的global_ref_genome
- **条件编译保护**：所有新代码用`#ifdef ENABLE_REF_GENOME_V4`包围
- **错误优雅处理**：失败时必须能回退到标准流程


## 12. 重要编译和类型修复

### 12.1 关键类型冲突解决方案

项目知识显示存在类型冲突问题，必须修复：

```c
// ❌ 错误：直接使用ul_idx_t会与现有类型冲突
// ul_idx_t *ref_index;

// ✅ 正确：使用专门的参考基因组索引类型
typedef struct {
    void *flt_tab;       // ha_ft_t* filter table  
    ha_pt_t *pt_idx;     // ha_pt_t* position table
} ref_ul_idx_t;
```

### 12.2 现有函数接口适配

在调用现有hifiasm函数时需要类型转换：

```c
// 在unitigs_map_to_reference_batch中：
int unitigs_map_to_reference_batch(ma_ug_t *unitigs,
                                   const ul_idx_t *ref_index,  // 接受ul_idx_t*
                                   // ...其他参数);

// 调用时的正确转换：
result = unitigs_map_to_reference_batch(
    unitigs,
    (const ul_idx_t*)ref->ul_index,  // ref_ul_idx_t* → ul_idx_t*
    &blocks, &count, opt
);
```

### 12.3 内存管理注意事项

```c
// ✅ 正确的UL索引构建（基于项目知识）
ref->ul_index = (ref_ul_idx_t*)malloc(sizeof(ref_ul_idx_t));
ref->ul_index->flt_tab = ha_ft_ul_gen(&asm_opt, &ref_unitigs, 19, 10, 5);
ref->ul_index->pt_idx = ha_pt_ul_gen(&asm_opt, ref->ul_index->flt_tab, &ref_unitigs, 19, 10, 5);

// ✅ 正确的清理顺序
void ref_genome_destroy(ref_genome_t *ref) {
    if (ref->ul_index) {
        if (ref->ul_index->pt_idx) ha_pt_destroy(ref->ul_index->pt_idx);
        if (ref->ul_index->flt_tab) ha_ft_destroy(ref->ul_index->flt_tab);
        free(ref->ul_index);
    }
    // ... 其他清理
}
```

当AI需要基于这个设计实现代码时，应该：

1. **优先理解现有UL流程**：重点关注uc_block_t数据结构和UL_INF全局管理
2. **严格遵循标记系统**：正确使用BLOCK_REF标记
3. **复用现有函数**：尽量调用现有比对和处理函数
4. **注重错误处理**：每个步骤都要有失败回退机制
5. **保持接口兼容**：不修改任何现有函数签名


