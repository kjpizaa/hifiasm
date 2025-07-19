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

> **速用指南**  
> 1. **加功能**：先查上表看副作用；保证弧对称、标志位唯一。  
> 2. **调参数**：改 `ARC_REF_W`、BFS 距离、泡阈值后跑 `tests/chr22_mini`。  
> 3. **打 debug**：请用 `[REF-ARC]` / `[GRAPH]` 前缀。  
> 4. **永远保持对称**：加 `v→w` 必同步 `w^1→v^1`。  

---
---

## 9 · Unitig pipeline & reference-UL hook

> 本节补足 *Unitig Management* 细节，说明 **参考基因组 ➜ 虚拟 UL block** 在哪一级注入、要保持哪些不变量。  
> 新增函数 / 标志已列在表中，**更改签名会导致 UL 全链路失效！**

| 阶段 | 关键函数 | 文件 | 说明 & 注意 |
|------|----------|------|-------------|
| **A. 线性链 / 全局链** | `mg_lchain_gen()` | `inter.cpp` 912-960 | Minimizer → linear chain；只读 `ma_hit_t`，参考基因组不参与。 |
| | *global-DP* | `inter.cpp` 803-877 | 将多条 linear chain 拼 global chain。 |
| **B. 单链解歧 → Unitig** | `ul_resolve()` | `inter.h` 107-108 | 把 global chain 解析成 unitig；可处理 **真实 UL** 与 **参考 block (`BLOCK_REF`)**。 |
| | `ul_refine_alignment()` | `inter.h` 110 | 二次比对修正边界；若 `BLOCK_REF` 置位则跳过波动剪枝，保留参考坐标。 |
| | `ul_realignment()` | `inter.h` 111-112 | 仅对真实 UL 调用，不触碰引用 UL。 |
| **C. Unitig 图构建** | `ma_ug_gen()` | `Overlaps.cpp` 36310-36479 | 压缩 1-in-1-out 链；`ARC_FLAG_REF` 弧正常参与权重，但不影响“度数”判定。 |
| | `ma_utg_t`, `ma_ug_t` | `Overlaps.h` 214 / 266 | **格式不变**：`a[j] = (vertex<<32 | edgeLen)`. |
| **D. 引用 UL 注入点** | `hifi_unitigs_map_to_reference_batch()` | `inter.cpp` (新增) | 把参考染色体按 `paths.tsv` 切段 → 对齐到 unitig → 生成 `uc_block_t`<br> · 在 block→`el` 高位加 `BLOCK_REF` 标志<br> · 推入 `push_uc_block_t()` 后再调 `sort_uc_block_qe()` |
| | `overlap_to_uc_block_ref_mode()` | `inter.cpp` (新增) | 将参考弧 (`ARC_FLAG_REF`) 转为 `uc_block_t`；保持对称。 |
| **E. 坐标拉伸 & 覆盖** | `extend_coordinates()` | `inter.cpp` 881-906 | 需允许 `te − ts` ≥ 800 Mb；改动请用 `uint64_t`. |
| | `ugl_cover_check()` | `inter.h` 114-115 | 遇到 `BLOCK_REF` 时放宽覆盖阈值（参考序列视为完美覆盖）。 |
| | `ug_occ_w()` | `gfa_ut.h` 43-44 | 计算 unitig 权重；参考 block 赋固定高分 `UL_REF_WEIGHT`。 |
| **F. 清洗 / 过滤** | `filter_ul_ug()` | `inter.h` 116 | 读 `BLOCK_REF`，**不**删除参考 block；真实 UL 的弱 block 仍可被滤掉。 |
| | `clean_contain_g()` / `dedup_contain_g()` | `inter.h` 127-128 | 遇到 `BLOCK_REF` 直接 `continue`；避免把参考块当成可删冗余。 |

### 标志位 / 宏一览

| 宏 | 位值 | 说明 |
|----|------|------|
| `ARC_FLAG_REF` | `1u<<30` | 标记 “参考基因组产生的 overlap arc” |
| `BLOCK_REF` | `1u<<15`  *(uc_block_t.el 高位)* | 指“此 unitig-block 来自参考 UL 路径” |
| `UL_REF_WEIGHT` | `900` | 给参考块在 `ug_occ_w()` 中的默认覆盖值 |

### 必守不变量 💡

1. **对称性**：每条参考弧 `v→w` 必伴随 `w^1→v^1`；每个 `uc_block_t` 也要双向插入。  
2. **签名稳定**：`ul_resolve()`、`ul_refine_alignment()` 原型严禁改；前端只透传 `BLOCK_REF`。  
3. **长度安全**：坐标/长度涉及染色体 (>800 Mb) 一律 `uint64_t`。  
4. **清洗豁免**：任何 `ARC_FLAG_REF` 弧、`BLOCK_REF` block **不得**在剪枝/去噪环节被删除。  
5. **先注入→再 cleanup**：注入参考弧后仅调用一次 `asg_cleanup()`，避免二次索引错位。

---

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

## 1 · 项目概述

本项目是hifiasm的增强版本，添加了**有参考基因组组装**功能。我们将参考基因组转换为**“虚拟Ultra-Long (UL)读段”**，通过现有UL处理流程实现HiFi图的gap填充和结构变异保护。

-----

## 2 · 核心设计理念 (修改代码前必读)

|概念         |实现方式                                                  |
|-----------|------------------------------------------------------|
|**参考基因组输入**|通过 `--ref-fasta` CLI选项指定FASTA文件路径                     |
|**虚拟UL转换** |将参考基因组转换为All_reads格式，利用现有UL索引机制                       |
|**集成时机**   |在原生hifiasm unitig图生成后进行，保持向后兼容                        |
|**坐标系转换**  |unitig作query，reference作target，输出uc_block_t格式          |
|**处理优先级**  |仅在HiFi证据缺失时填充gap，保护真实结构变异                             |
|**代码复用率**  |98%复用现有UL处理函数 (`ul_refine_alignment`, `ul_clean_gfa`等)|

-----

## 3 · 关键文件结构

|文件路径                              |功能描述                                    |
|----------------------------------|----------------------------------------|
|`src/ref_genome.h`                |参考基因组数据结构定义，Phase A核心API                |
|`src/ref_genome.cpp`              |参考基因组处理：FASTA加载、UL索引构建、All_reads转换      |
|`src/Assembly.cpp`                |增强的主工作流，Virtual ONT转换和UL流程集成            |
|`src/main.cpp`                    |CLI集成，参考基因组模式检测和错误处理                    |
|`src/CommandLines.h`              |添加 `char *ref_fasta` 字段到 `hifiasm_opt_t`|
|`src/CommandLines.cpp`            |CLI参数解析，添加 `--ref-fasta` 选项             |
|`src/reference_output_format.txt` |Phase B输出格式说明和uc_block_t转换示例            |
|`src/reference_as_virtual_ont.txt`|Virtual ONT转换策略详细说明                     |

-----

## 4 · 新增CLI选项

```bash
# 新增的参考基因组选项
--ref-fasta FILE     参考基因组FASTA文件路径，启用有参考组装模式
```

完整用法示例：

```bash
./hifiasm --ref-fasta hg38.fa -o output -t 32 hifi_reads.fq.gz
```

-----

## 5 · 核心数据结构

### 5.1 参考基因组结构 (`ref_genome_t`)

```c
typedef struct {
    char *fasta_path;                    // FASTA文件路径
    uint32_t n_seq;                      // 染色体数量
    ref_chromosome_t *chromosomes;       // 染色体数组
    uint64_t total_length;               // 总长度
    char *merged_seq;                    // 统一序列（可选内存优化）
    ul_idx_t virtual_ul_index;          // UL索引（用于unitig→reference比对）
    All_reads *all_reads_ref;           // All_reads格式（每条染色体=1个read）
    uint8_t index_built;                // UL索引是否已构建
} ref_genome_t;
```

### 5.2 配置结构 (`ref_config_t`)

```c
typedef struct {
    uint64_t chunk_size;                // 分块大小 (默认100kb)
    uint32_t min_seq_len;              // 最小序列长度 (默认1000bp)
    uint8_t enable_hpc;                // HPC压缩开关
    uint8_t memory_optimization;       // 内存优化开关
} ref_config_t;
```

-----

## 6 · 关键处理流程

### Phase A: 参考基因组预处理

```c
// 1. 初始化和FASTA加载
ref_genome_t *ref = ref_genome_init();
ref_genome_load_fasta(ref, opt->ref_fasta);

// 2. 构建统一序列和All_reads转换
ref_config_t config = ref_config_default();
ref_genome_build_unified_sequence(ref, &config);
ref_genome_convert_to_all_reads(ref, &config);

// 3. 构建UL索引（k=19, w=10）
ref_genome_build_ul_index(ref);
prepare_reference_for_virtual_ont(ref);
```

// Phase B: unitig → reference比对的真实实现
// 基于项目知识中inter.cpp的ENABLE_REF_GENOME_V4部分

// ===============================
// 1. 批量处理所有unitigs到参考基因组的比对
// ===============================

int process_all_unitigs_to_reference(ma_ug_t *hifi_unitigs, 
                                     ref_genome_t *ref,
                                     uc_block_t **uc_blocks, 
                                     uint64_t *n_blocks)
{
    if (!hifi_unitigs || !ref || !uc_blocks || !n_blocks) {
        fprintf(stderr, "[ERROR] Invalid parameters for unitig-reference mapping\n");
        return -1;
    }

    fprintf(stderr, "[M::%s] Processing %u unitigs against reference...\n", 
            __func__, hifi_unitigs->u.n);

    // 调用项目知识中的真实函数
    extern int unitigs_map_to_reference_batch(ma_ug_t *unitigs, 
                                              const ul_idx_t *ref_index,
                                              uc_block_t **out_blocks, 
                                              uint64_t *out_count,
                                              const hifiasm_opt_t *opt);

    // 从全局变量获取选项（与现有代码保持一致）
    extern hifiasm_opt_t asm_opt;

    // 使用项目知识中已实现的批量处理函数
    int result = unitigs_map_to_reference_batch(
        hifi_unitigs,
        (const ul_idx_t*)ref->ul_index,  // 参考基因组的UL索引
        uc_blocks,
        n_blocks,
        &asm_opt
    );

    if (result != 0) {
        fprintf(stderr, "[ERROR] Failed to map unitigs to reference\n");
        return -1;
    }

    if (*n_blocks == 0) {
        fprintf(stderr, "[WARNING] No valid unitig-reference alignments found\n");
        return 0;
    }

    fprintf(stderr, "[M::%s] Generated %lu uc_blocks from unitig-reference mapping\n", 
            __func__, (unsigned long)*n_blocks);

    // 关键：验证坐标转换正确性
    // 在项目知识中，qs/qe=unitig坐标，ts/te=reference坐标
    for (uint64_t i = 0; i < *n_blocks && i < 5; i++) {
        uc_block_t *block = &(*uc_blocks)[i];
        fprintf(stderr, "[DEBUG] Block %lu: unitig=%u, ref_chr=%u, "
                       "unitig_range=[%u,%u), ref_range=[%u,%u), strand=%c\n",
                i, block->hid, block->aidx, 
                block->qs, block->qe, block->ts, block->te,
                block->rev ? '-' : '+');
    }

    return 0;
}

// ===============================
// 2. 完整的Phase B实现（基于项目知识）
// ===============================

/**
 * 完整的Phase B实现，基于项目知识中inter.cpp的真实代码模式
 * 这个函数展示了如何正确调用process_all_unitigs_to_reference
 */
int execute_phase_b_unitig_reference_mapping(ma_ug_t *hifi_unitigs, 
                                            ref_genome_t *ref)
{
    fprintf(stderr, "\n=== Phase B: Unitig → Reference Mapping ===\n");
    
    // Step 1: 验证输入参数
    if (!hifi_unitigs || !ref) {
        fprintf(stderr, "[ERROR] Invalid input for Phase B\n");
        return -1;
    }

    if (hifi_unitigs->u.n == 0) {
        fprintf(stderr, "[WARNING] No unitigs available for mapping\n");
        return 0;
    }

    if (!ref->ul_index) {
        fprintf(stderr, "[ERROR] Reference UL index not built\n");
        return -1;
    }

    // Step 2: 执行批量映射
    uc_block_t *uc_blocks = NULL;
    uint64_t n_blocks = 0;

    int result = process_all_unitigs_to_reference(hifi_unitigs, ref, 
                                                 &uc_blocks, &n_blocks);
    
    if (result != 0) {
        fprintf(stderr, "[ERROR] Phase B mapping failed\n");
        return -1;
    }

    // Step 3: 统计信息
    fprintf(stderr, "[INFO] Phase B completed successfully:\n");
    fprintf(stderr, "  - Input unitigs: %u\n", hifi_unitigs->u.n);
    fprintf(stderr, "  - Reference sequences: %u\n", ref->n_seq);
    fprintf(stderr, "  - Generated uc_blocks: %lu\n", (unsigned long)n_blocks);

    // Step 4: 质量检查
    uint64_t valid_blocks = 0;
    uint64_t total_coverage = 0;
    
    for (uint64_t i = 0; i < n_blocks; i++) {
        uc_block_t *block = &uc_blocks[i];
        if (block->el > 0 && block->qe > block->qs && block->te > block->ts) {
            valid_blocks++;
            total_coverage += (block->qe - block->qs);
        }
    }

    fprintf(stderr, "  - Valid blocks: %lu (%.1f%%)\n", 
           (unsigned long)valid_blocks, 
           n_blocks > 0 ? 100.0 * valid_blocks / n_blocks : 0.0);
    
    fprintf(stderr, "  - Total unitig coverage: %lu bp\n", 
           (unsigned long)total_coverage);

    // Step 5: 将结果传递给下一阶段（Phase C）
    // 在项目知识中，这些uc_blocks会被传递给UL处理流程
    extern int integrate_reference_blocks_to_existing_ul_pipeline(ma_ug_t *unitigs,
                                                                 const ul_idx_t *ref_index,
                                                                 const hifiasm_opt_t *opt);
    extern hifiasm_opt_t asm_opt;

    fprintf(stderr, "[INFO] Proceeding to Phase C: UL pipeline integration...\n");
    
    // 临时存储uc_blocks到全局变量（如果需要）
    // 或者直接调用集成函数
    result = integrate_reference_blocks_to_existing_ul_pipeline(
        hifi_unitigs, 
        (const ul_idx_t*)ref->ul_index, 
        &asm_opt
    );

    // 清理资源
    free(uc_blocks);

    if (result == 0) {
        fprintf(stderr, "[INFO] Phase B → Phase C integration completed\n");
    } else {
        fprintf(stderr, "[WARNING] Phase C integration failed, but Phase B succeeded\n");
    }

    fprintf(stderr, "=== Phase B Completed ===\n\n");
    return result;
}

// ===============================
// 3. 用法示例（基于项目知识中的调用模式）
// ===============================

/**
 * 这是在Assembly.cpp中如何调用Phase B的示例
 * 基于项目知识中execute_reference_guided_assembly的模式
 */
void example_usage_in_assembly(void)
{
    // 假设已经有了hifi_unitigs和global_ref_genome
    extern ma_ug_t *hifi_unitigs;  // 来自现有的组装流程
    extern ref_genome_t *global_ref_genome;  // 来自Phase A

    if (hifi_unitigs && global_ref_genome) {
        // 执行Phase B映射
        int result = execute_phase_b_unitig_reference_mapping(
            hifi_unitigs, 
            global_ref_genome
        );

        if (result == 0) {
            fprintf(stderr, "[INFO] Reference-guided enhancement completed\n");
        } else {
            fprintf(stderr, "[WARNING] Reference-guided enhancement failed, "
                           "continuing with standard assembly\n");
        }
    }
}

// ===============================
// 4. 与hifiasm Unitig Management的完美对应关系
// ===============================

 *    ✅ ma_ug_t: Unitig Graph - 我们的hifi_unitigs正是此类型
 *    ✅ ma_utg_t: 单个Unitig - 包含len, circ, s等字段
 *    ✅ uc_block_t: Unitig Block - 包含qs/qe/ts/te坐标映射
 *    ✅ asg_t: Assembly Graph - unitigs->g就是这个结构
 * 
 * 🔧 Unitig Construction Process对应：
 *    ✅ Linear Chaining: mg_lchain_gen → 我们复用现有链构建
 *    ✅ Global Chaining: → 现有的全局链处理
 *    ✅ Unitig Resolution: ul_resolve → 我们直接调用此函数！
 *    ✅ Graph Construction: ma_ug_t → 输出标准unitig图
 * 
 * 🚀 Block Management完美匹配：
 *    ✅ push_uc_block_t: 添加新blocks → 我们生成uc_blocks数组
 *    ✅ sort_uc_block_qe: 按query end排序 → 标准UL流程包含
 *    ✅ ul_refine_alignment: 改进比对 → 现有UL流程自动调用
 *    ✅ extend_coordinates: 扩展坐标 → 现有UL流程自动调用
 * 
 * 🎯 Ultra-Long Read Integration策略：
 *    ✅ 我们将参考基因组作为"Virtual Ultra-Long Reads"
 *    ✅ ul_resolve和ul_realignment自动处理这些虚拟reads
 *    ✅ 完美复用现有Ultra-Long读数处理基础设施
 * 
 * 📊 关键坐标系统（附件确认）：
 *    ✅ qs/qe: query start/end → unitig坐标
 *    ✅ ts/te: target start/end → reference坐标  
 *    ✅ hid: 标识符 → unitig ID
 *    ✅ rev: 方向标记 → 正/反向匹配
 * 
 * 🔄 与Assembly Pipeline集成：
 *    ✅ Graph Algorithms: mg_shortest_k等 → 自动可用
 *    ✅ Alignment System: ha_get_ug_candidates → 我们直接使用
 *    ✅ Final Assembly: trans_base_infer → 标准流程处理
 * 


/* ------------------------------------------------------------
 * Phase B 结束，此时已经得到 unitig-graph  ug->g 及 ug->u
 * ---------------------------------------------------------- */
asg_cleanup(ug->g);       /* 最后一次清扫原始弧 */
asg_symm(ug->g);

/* ------------------------------------------------------------
 * Phase C : Reference-guided UL pipeline
 * ---------------------------------------------------------- */
#ifdef ENABLE_REF_GENOME_V4
if (asm_opt.ref_fasta && asm_opt.ref_fasta[0]) {
    // 🔧 使用真实的全局unitig图变量（项目知识中确认存在）
    extern ma_ug_t *ug;
    
    if (ug && ug->u.n > 0) {
        fprintf(stderr, "[M::%s] Starting reference-guided assembly with %u unitigs\n", 
                __func__, ug->u.n);
        
        /* 1. 执行参考基因组指导的组装
         *    - 内部调用 integrate_reference_blocks_to_existing_ul_pipeline()
         *    - 生成参考基因组blocks并添加到全局UL_INF
         *    - 调用标准unitig管理流程包括ul_resolve() */
        int ref_result = execute_reference_guided_assembly(ug, &asm_opt);
        
        if (ref_result == 0) {
            fprintf(stderr, "[M::%s] Reference-guided assembly completed successfully\n", __func__);
            
            /* ---- 下面 100% 复用现有 UL 处理函数 ---- */
            /* 这些函数在 execute_reference_guided_assembly 内部已经被调用，
             * 但可以根据需要进行额外的清理和优化 */
            
            /* 2. 最终清理 GFA 标记、孤立弧等收尾操作 */
            extern void ul_clean_gfa(ma_ug_t *ug);
            ul_clean_gfa(ug);
            
            fprintf(stderr, "[M::%s] Reference-guided unitig processing completed\n", __func__);
        } else {
            fprintf(stderr, "[WARNING] Reference-guided assembly failed, continuing with standard assembly\n");
        }
    } else {
        fprintf(stderr, "[WARNING] Unitig graph not available for reference-guided assembly\n");
    }
}
#endif

/* Phase C 结束，ug->g 已融合参考信息（如果启用），进入 layout → polish ↓ */

## 7 · 重要约束与注意事项

### 7.1 代码修改原则

- **保护原有功能**: 所有现有功能必须100%保持
- **向后兼容**: 不传入 `--ref-fasta` 时行为完全不变
- **错误优雅降级**: 参考基因组处理失败时回退到标准hifiasm
- **最小化修改**: 优先扩展而非修改现有代码
- **范围约束**: `Process_Read.h` 以及其他源码，除 `// ENABLE_REF_GENOME_V4` 包裹的部分外尽量保持原状，不新增全局变量或修改现有函数；优先复用已有模块和变量。

### 7.2 内存管理

```c
// 正确的资源清理顺序
if (graphs) destroy_four_graph_overlaps(graphs);
if (ref_genome) ref_genome_destroy(ref_genome);
free(uc_blocks);  // uc_block_t数组
destory_All_reads(&virtual_ont_reads);  // All_reads结构
```

### 7.3 数据格式兼容性

- **ul_idx_t**: 使用k=19, w=10参数与现有UL流程兼容
- **uc_block_t**: 输出格式必须与现有UL函数100%兼容
- **All_reads**: 每条染色体作为独立read，维持现有数据结构

-----

## 8 · 测试与验证

### 8.1 基本功能测试

```bash
# 测试1: 标准功能不受影响
./hifiasm -o test1 -t 8 hifi.fq.gz

# 测试2: 参考基因组模式
./hifiasm --ref-fasta ref.fa -o test2 -t 8 hifi.fq.gz

# 测试3: 错误处理
./hifiasm --ref-fasta nonexistent.fa -o test3 -t 8 hifi.fq.gz
```

### 8.2 预期输出

- 成功时：增强的unitig图，更少的gap
- 失败时：自动回退，与标准hifiasm输出相同
- 错误时：清晰的错误信息和建议

-----

## 9 · 性能优化考虑

|优化点      |实现方法                        |
|---------|----------------------------|
|**内存使用** |可选启用merged_seq释放，减少800MB内存占用|
|**I/O优化**|缓存机制，避免重复解析大型参考基因组          |
|**并行处理** |复用现有多线程UL处理框架               |
|**索引复用** |UL索引一次构建，多次使用               |



> **提示**: 在编写代码时，优先查看 `reference_as_virtual_ont.txt` 了解Virtual ONT转换策略，参考 `reference_output_format.txt` 了解输出格式要求。所有实现都应该最大化复用现有hifiasm代码，保持98%以上的代码复用率。
