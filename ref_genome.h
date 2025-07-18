#ifndef __REF_GENOME_H__
#define __REF_GENOME_H__

#include <stdint.h>
#include <stdio.h>
#include "Process_Read.h"

#ifdef __cplusplus
extern "C" {
#endif

// 前向声明，避免循环依赖
struct ha_pt_s;

/**
 * @brief 参考基因组染色体结构 (完整版)
 */
typedef struct {
    char *name;              // 染色体名称
    char *seq;               // 染色体序列（压缩后可能释放以节省内存）
    uint64_t length;         // 染色体长度
    uint64_t offset;         // 在unified sequence中的偏移量
} ref_chromosome_t;

/**
 * @brief UL索引结构（修正版 - 重命名避免冲突）
 * 
 * ⚠️ 修复说明：
 * - 原项目中Overlaps.h已经定义了ul_idx_t结构
 * - 为了避免"conflicting declaration"错误，这里重命名为ref_ul_idx_t
 * - 保持功能完全一致，只是名称不同
 */
typedef struct {
    void *flt_tab;          // filter table指针
    struct ha_pt_s *pt_idx; // position table指针
} ref_ul_idx_t;

/**
 * @brief 虚拟ONT数据结构 (完整版，保持项目兼容性)
 */
typedef struct {
    int ref_as_ont_read;          // 是否将参考作为ONT read
    uint64_t ont_read_length;     // ONT read长度
    uint32_t chunk_size;          // chunk大小（保持兼容性）
    uint64_t n_chunks;            // chunk数量（保持兼容性）
    All_reads *virtual_reads;     // 指向虚拟reads数据（引用）
    ref_ul_idx_t *virtual_ul_index;   // 虚拟UL索引（引用，避免double free）
} virtual_ont_data_t;

/**
 * @brief 参考基因组配置结构 (完整版)
 */
typedef struct {
    uint32_t min_chromosome_length;     // 最小染色体长度
    int use_homopolymer_compression;    // 是否使用HPC（从全局asm_opt读取）
    uint32_t chunk_size;                // chunk大小（保持兼容性）
    uint32_t overlap_threshold;         // overlap阈值
    int enable_caching;                 // 是否启用缓存
} ref_config_t;

/**
 * @brief 参考基因组主结构 (完整版，基于项目知识但修正了问题)
 * 
 * 设计思路：
 * - 保留项目知识中的所有原有功能和兼容性
 * - 修正了内存优化、HPC支持、错误检查等问题
 * - 每条染色体可作为独立read（修正版）
 * - 支持B阶段unitig→reference比对
 */
typedef struct {
    char *fasta_path;                   // FASTA文件路径
    ref_chromosome_t *chromosomes;      // 染色体数组
    uint32_t n_seq;                     // 序列数量
    uint64_t total_length;              // 总长度
    char *merged_seq;                   // 合并序列（保持兼容性，可选内存优化）
    uint64_t n_bases;                   // 碱基总数
    uint64_t n_chunks;                  // chunk数量（保持兼容性）
    int index_built;                    // 索引是否已构建
    ref_ul_idx_t *ul_index;             // UL索引（使用重命名后的类型）
    virtual_ont_data_t *virtual_ont_data; // 虚拟ONT数据
    All_reads *all_reads_ref;           // All_reads格式：每条染色体作为独立read
} ref_genome_t;

// ===== Phase A: 参考基因组索引化处理核心API (完整版) =====

/**
 * @brief 初始化参考基因组结构
 * @return 成功返回初始化的结构体指针，失败返回NULL
 */
ref_genome_t* ref_genome_init(void);

/**
 * @brief 销毁参考基因组结构（修正版：避免double free）
 * @param ref 参考基因组结构指针
 */
void ref_genome_destroy(ref_genome_t *ref);

/**
 * @brief 加载参考基因组FASTA文件（修正版：一次遍历+优化）
 * 修正：一次遍历+kvec动态扩展，快速大写转换（seq[j] &= 0xDF）
 * @param ref 参考基因组结构指针
 * @param fasta_path FASTA文件路径
 * @return 成功返回0，失败返回非0
 */
int ref_genome_load_fasta(ref_genome_t *ref, const char *fasta_path);

/**
 * @brief 构建统一序列（修正版：内存优化但保留原功能）
 * 保留原功能但可选内存优化：避免800MB*2内存占用
 * @param ref 参考基因组结构指针
 * @param config 配置参数
 * @return 成功返回0，失败返回非0
 */
int ref_genome_build_unified_sequence(ref_genome_t *ref, const ref_config_t *config);

/**
 * @brief 转换为All_reads格式（修正版：每条染色体作为独立read）
 * 修正：7条染色体→7个reads，添加错误检查，可选内存优化
 * @param ref 参考基因组结构指针
 * @param config 配置参数
 * @return 成功返回0，失败返回非0
 */
int ref_genome_convert_to_all_reads(ref_genome_t *ref, const ref_config_t *config);

/**
 * @brief 构建UL索引（修正版：支持多个参考序列，正确的HPC支持）
 * 修正：每条染色体对应unitig，正确HPC开关，为B阶段unitig→reference做准备
 * @param ref 参考基因组结构指针
 * @return 成功返回0，失败返回非0
 */
int ref_genome_build_ul_index(ref_genome_t *ref);

/**
 * @brief 准备虚拟ONT数据（修正版：避免浅拷贝但保持兼容性）
 * 修正：避免double free，引用模式，保持All_reads兼容性
 * @param ref 参考基因组结构指针
 * @return 成功返回0，失败返回非0
 */
int prepare_reference_for_virtual_ont(ref_genome_t *ref);

// ===== Phase B: 缓存和持久化功能 (完整版) =====

/**
 * @brief 保存缓存文件（项目知识原有功能 - 完整版）
 * @param ref 参考基因组结构指针
 * @param cache_prefix 缓存文件前缀
 * @return 成功返回0，失败返回非0
 */
int ref_genome_save_cache(const ref_genome_t *ref, const char *cache_prefix);

/**
 * @brief 加载缓存文件（项目知识原有功能 - 完整版）
 * @param ref 参考基因组结构指针
 * @param cache_prefix 缓存文件前缀
 * @return 成功返回0，失败返回非0
 */
int ref_genome_load_cache(ref_genome_t *ref, const char *cache_prefix);

/**
 * @brief 检查缓存有效性（项目知识原有功能 - 增强版）
 * 增强：添加文件修改时间检查
 * @param cache_prefix 缓存文件前缀
 * @param fasta_path FASTA文件路径
 * @return 缓存有效返回1，否则返回0
 */
int ref_genome_cache_is_valid(const char *cache_prefix, const char *fasta_path);

// ===== Phase C: 增强功能API (完整保留项目知识中的功能) =====

/**
 * @brief 加载缓存并构建UL索引（项目知识原有功能）
 * @param ref 参考基因组结构指针
 * @param cache_prefix 缓存文件前缀
 * @return 成功返回0，失败返回非0
 */
int ref_genome_load_cache_with_ul_index(ref_genome_t *ref, const char *cache_prefix);

/**
 * @brief 创建参考映射缓存（项目知识原有功能）
 * @param ref 参考基因组结构指针
 * @param cache_prefix 缓存文件前缀
 * @return 成功返回0，失败返回非0
 */
int create_reference_mapping_cache(ref_genome_t *ref, const char *cache_prefix);

/**
 * @brief 验证虚拟ONT质量（项目知识原有功能）
 * @param vont 虚拟ONT数据结构
 * @return 成功返回0，失败返回非0
 */
int validate_virtual_ont_quality(const virtual_ont_data_t *vont);

/**
 * @brief 构建虚拟ONT UL索引（项目知识原有功能）
 * @param vont 虚拟ONT数据结构
 * @return 成功返回0，失败返回非0
 */
int build_virtual_ont_ul_index(virtual_ont_data_t *vont);

// ===== Phase D: 配置和工具函数API (完整版) =====

/**
 * @brief 获取默认配置（修正版：基于全局HPC设置）
 * @return 默认配置结构体
 */
ref_config_t ref_config_default(void);

/**
 * @brief 验证参考基因组完整性（项目知识原有功能）
 * @param ref 参考基因组结构指针
 * @return 成功返回0，失败返回非0
 */
int validate_reference_genome(const ref_genome_t *ref);

/**
 * @brief 打印参考基因组统计信息（项目知识原有功能 - 完整增强版）
 * 增强：包含虚拟ONT详细信息，内存使用统计
 * @param ref 参考基因组结构指针
 */
void ref_genome_print_stats(const ref_genome_t *ref);

/**
 * @brief 获取染色体偏移量（项目知识原有功能）
 * @param ref 参考基因组结构指针
 * @param chr_name 染色体名称
 * @param offset 输出偏移量
 * @param length 输出长度
 * @return 成功返回0，失败返回非0
 */
int get_chromosome_offset(const ref_genome_t *ref, const char *chr_name, 
                         uint64_t *offset, uint64_t *length);

/**
 * @brief 清理临时数据（项目知识原有功能）
 * @param ref 参考基因组结构指针
 * @param keep_sequences 是否保留序列数据
 */
void ref_genome_cleanup_temp_data(ref_genome_t *ref, int keep_sequences);

// ===== Phase E: 高级功能API (完整保留项目知识中的功能) =====

/**
 * @brief 优化参考基因组内存使用（项目知识原有功能）
 * @param ref 参考基因组结构指针
 * @param optimization_level 优化级别 (0=无优化, 1=基础优化, 2=激进优化)
 * @return 成功返回0，失败返回非0
 */
int ref_genome_optimize_memory(ref_genome_t *ref, int optimization_level);

/**
 * @brief 重建染色体序列（从合并序列恢复）（项目知识原有功能）
 * @param ref 参考基因组结构指针
 * @return 成功返回0，失败返回非0
 */
int ref_genome_restore_chromosome_sequences(ref_genome_t *ref);

/**
 * @brief 创建参考基因组摘要（项目知识原有功能）
 * @param ref 参考基因组结构指针
 * @param output_file 输出文件路径
 * @return 成功返回0，失败返回非0
 */
int ref_genome_create_summary(const ref_genome_t *ref, const char *output_file);

// ===== 常量定义 (完整版) =====

#define REF_DEFAULT_CHUNK_SIZE      50000     // 保持兼容性
#define REF_DEFAULT_MIN_CHR_LEN     1000
#define REF_DEFAULT_OVERLAP_THR     500
#define REF_DEFAULT_K_SIZE          19        // UL索引k-mer大小
#define REF_DEFAULT_W_SIZE          10        // UL索引窗口大小  
#define REF_DEFAULT_CUTOFF          5         // UL索引cutoff

#define REF_CACHE_VERSION           1
#define REF_MAX_CHROMOSOMES         1000
#define REF_MAX_NAME_LEN            256

// 内存优化级别常量
#define REF_MEMORY_OPT_NONE         0
#define REF_MEMORY_OPT_BASIC        1
#define REF_MEMORY_OPT_AGGRESSIVE   2

// ===== 错误码定义 =====

#define REF_SUCCESS                 0
#define REF_ERROR_NULL_POINTER      -1
#define REF_ERROR_FILE_NOT_FOUND    -2
#define REF_ERROR_MEMORY_ALLOC      -3
#define REF_ERROR_INVALID_FORMAT    -4
#define REF_ERROR_INDEX_BUILD       -5
#define REF_ERROR_CACHE_CORRUPT     -6
#define REF_ERROR_INVALID_CONFIG    -7
#define REF_ERROR_OPTIMIZATION      -8

#ifdef __cplusplus
}
#endif

#endif /* __REF_GENOME_H__ */

/*
 * 📋 完整版总结：
 * 
 * ✅ 完全保留了项目知识中的所有原有功能：
 * 
 * 🔧 Phase A - 核心API（修正版）：
 * - ref_genome_init/destroy：完整内存管理 + double free修复
 * - ref_genome_load_fasta：一次遍历优化 + 快速大写转换
 * - ref_genome_build_unified_sequence：保留原功能 + 内存优化选项
 * - ref_genome_convert_to_all_reads：每条染色体独立read + 错误检查
 * - ref_genome_build_ul_index：支持多序列 + 正确HPC开关
 * - prepare_reference_for_virtual_ont：避免浅拷贝 + 保留兼容性
 * 
 * ✅ Phase B - 缓存管理（完整版）：
 * - ref_genome_save/load_cache：完整二进制缓存系统
 * - ref_genome_cache_is_valid：缓存有效性验证（增强：文件时间检查）
 * 
 * ✅ Phase C - 增强功能（完整保留）：
 * - ref_genome_load_cache_with_ul_index：带索引的缓存加载
 * - create_reference_mapping_cache：映射缓存创建
 * - validate_virtual_ont_quality：虚拟ONT质量验证
 * - build_virtual_ont_ul_index：虚拟ONT索引构建
 * 
 * ✅ Phase D - 配置和工具（完整版+增强）：
 * - ref_config_default：默认配置生成（修正：基于全局HPC）
 * - validate_reference_genome：完整性验证
 * - ref_genome_print_stats：详细统计（增强：虚拟ONT + 内存信息）
 * - get_chromosome_offset：染色体偏移查询
 * - ref_genome_cleanup_temp_data：临时数据清理
 * 
 * ✅ Phase E - 高级功能（完整保留）：
 * - ref_genome_optimize_memory：多级别内存优化
 * - ref_genome_restore_chromosome_sequences：序列恢复
 * - ref_genome_create_summary：摘要生成
 * 
 * 🔧 关键修复：
 * - ✅ ul_idx_t → ref_ul_idx_t（避免类型冲突）
 * - ✅ 所有相关结构体类型更新一致
 * - ✅ 保留所有原有API接口
 * - ✅ 增强错误码和常量定义
 * - ✅ 完整的功能分组和文档
 * 
 * 💡 兼容性保证：
 * - 保持与现有代码100%接口兼容
 * - 保留所有chunk、merged_seq等兼容性概念
 * - 只修复编译错误，不改变功能行为
 * - 向后兼容所有现有调用代码
 */