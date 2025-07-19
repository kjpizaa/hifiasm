### hifiasm有参考基因组组装 - AI开发助手快速指南

> **受众:** GitHub Copilot、ChatGPT、Claude、Cursor AI等AI编程助手  
> **人类开发者:** 详细设计文档请参考 `docs/ref_guided_assembly_design.md`

-----

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
