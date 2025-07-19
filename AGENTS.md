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

> 以上行号基于 `ec9a8b` commit；若 upstream 更新请用 `grep -n` 校准。  
> 将此表格追加后：`git add AGENTS.md && git commit -m "docs: expand cheat-sheet with overlap & graph APIs"`。
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

### Phase B: unitig → reference比对

```c
// 使用现有UL比对函数，注意坐标转换
uc_block_t *uc_blocks;
uint64_t n_blocks;
process_all_unitigs_to_reference(hifi_unitigs, ref, &uc_blocks, &n_blocks);

// 关键：qs/qe=unitig坐标，ts/te=reference坐标
```

### Phase C: UL流程集成

```c
// 100%复用现有UL处理函数
ul_refine_alignment(&uopt, unitigs->g);
extend_coordinates(unitigs, uc_blocks, n_blocks);
update_ug_arch_ul_mul(unitigs);
filter_ul_ug(unitigs);
ul_clean_gfa(unitigs);
```

主流程会在生成unitig图完成后调用 `execute_reference_guided_assembly(ug, &asm_opt)`，
将unitig与参考基因组对齐，并把生成的 `uc_block_t` 结果纳入上述UL处理管线。

-----

## 7 · 重要约束与注意事项

### 7.1 代码修改原则

- **保护原有功能**: 所有现有功能必须100%保持
- **向后兼容**: 不传入 `--ref-fasta` 时行为完全不变
- **错误优雅降级**: 参考基因组处理失败时回退到标准hifiasm
- **最小化修改**: 优先扩展而非修改现有代码

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

-----

## 10 · 常见问题排查

### 10.1 编译错误

- 确保所有新增的.h文件都包含了必要的include
- 检查CommandLines.h中是否正确添加了ref_fasta字段
- 验证函数声明与实现的一致性

### 10.2 运行时错误

- 参考基因组文件权限和存在性检查
- 内存分配失败的优雅处理
- UL索引构建失败的回退机制

### 10.3 输出异常

- 验证uc_block_t格式的正确性
- 检查坐标转换 (query↔target) 是否正确
- 确认UL流程集成点的数据一致性

-----

## 11 · 开发状态与里程碑

- [x] **Phase A**: 参考基因组预处理 (ref_genome.cpp)
- [x] **Phase B**: Virtual ONT转换策略 (Assembly.cpp)
- [x] **Phase C**: CLI集成与错误处理 (main.cpp, CommandLines.cpp)
- [x] **Phase D**: 测试框架与文档完善
- [ ] **Phase E**: 性能优化与边界案例处理

-----

> **提示**: 在编写代码时，优先查看 `reference_as_virtual_ont.txt` 了解Virtual ONT转换策略，参考 `reference_output_format.txt` 了解输出格式要求。所有实现都应该最大化复用现有hifiasm代码，保持98%以上的代码复用率。
