# AGENTS.md - Hifiasm RSP集成项目指南

## 项目概述：RSP集成到Hifiasm基因组组装器

本项目旨在将参考引导序列放置（Reference-guided Sequence Placement, RSP）功能集成到Hifiasm基因组组装器中。RSP系统通过利用参考基因组比对信息来指导unitig构建过程，实现reads的正确组装和unitig的智能扩展合并。

## 核心文件和目录结构

### 主要实现文件
- `paste.txt` - 包含原始`ma_ug_gen_primary`函数实现（已存在）
- `rsp_integration_solution.txt` - 完整的RSP集成规范和实现指南（已存在）
- `Assembly.cpp` - 需修改的主要组装逻辑文件
- `CommandLines.h` - 需添加RSP选项的头文件
- `CommandLines.cpp` - 需添加RSP参数处理的实现文件

### RSP集成组件结构
```
src/
├── rsp_integration.cpp        # RSP核心实现（新增）
├── rsp_integration.h          # RSP头文件定义（新增）
├── Assembly.cpp               # 修改ma_ug_gen_primary函数
├── CommandLines.{h,cpp}       # 添加RSP命令行选项
└── Overlaps.cpp              # 复用现有的overlap处理函数
```

## 与现有模块的兼容性检查

### 重用现有数据结构
```cpp
// 1. 重用现有的ma_ug_t和ma_utg_t结构（来自项目知识）
// 位置：Assembly.h中定义，无需新增
typedef struct {
    ma_utg_v u;           // unitig向量
    asg_t *g;            // unitig图
} ma_ug_t;

// 2. 重用现有的asg_t结构（来自项目知识）
// 位置：已在Assembly.h中定义，用于overlap图

// 3. 重用现有的hifiasm_opt_t结构（来自CommandLines.h）
// 需要扩展添加RSP字段
```

### 重用现有函数
```cpp
// 从项目知识中确认的现有函数，直接复用：
asg_arc_n(g, v)              // 获取顶点v的弧数量
asg_arc_a(g, v)              // 获取顶点v的弧数组
asg_arc_len(arc)             // 获取弧的长度
kv_pushp(ma_utg_t, ug->u, &p) // 向unitig向量添加元素
kdq_init(uint64_t)           // 初始化队列
kdq_push/kdq_unshift        // 队列操作
asg_cleanup(g)               // 图清理
```

## RSP集成架构

### 核心概念
RSP集成通过以下步骤实现智能unitig构建和扩展：

1. **预处理阶段**：从参考比对信息构建RSP段
2. **预标记阶段**：标记RSP段内部reads，防止错误起点
3. **预构建阶段**：创建RSP指导的unitigs
4. **扩展检测阶段**：在标准算法中检测并扩展RSP unitigs
5. **合并优化阶段**：实现跨区域的unitig合并

### 新增数据结构（必须新增）

#### RSP段结构定义
```cpp
// 位置：新增文件 rsp_integration.h
// 作用：存储参考序列上连续比对的reads段信息
typedef struct {
    uint32_t *read_ids;      // 输入：段内reads ID数组
    uint8_t *directions;     // 输入：reads方向数组(0=正向,1=反向)
    uint32_t n_reads;        // 输入：reads数量
    uint32_t chr_id;         // 输入：染色体ID
} rsp_segment_t;

// 作用：管理所有RSP段的集合
typedef struct {
    rsp_segment_t *segments; // 输入：RSP段数组
    uint32_t n_segments;     // 输入：段数量
} rsp_segments_t;

// 作用：跟踪RSP unitig的端点信息，支持扩展
typedef struct {
    uint32_t segment_id;     // 输入：所属RSP段ID
    uint8_t is_tail;         // 输入：是否为尾端点(0=头,1=尾)
} rsp_endpoint_info_t;
```

### 命令行参数扩展

#### 在hifiasm_opt_t中添加字段
```cpp
// 位置：CommandLines.h中的hifiasm_opt_t结构体
// 作用：添加RSP功能的配置参数
typedef struct {
    // ... 现有字段保持不变 ...
    
    // 新增RSP相关字段
    char *ref_alignments;     // 输入：参考比对文件路径(.paf格式)
    int32_t rsp_min_mapq;     // 输入：最小映射质量阈值(默认:10)
    int32_t rsp_max_gap;      // 输入：最大允许间隙大小(默认:50000)
    uint8_t enable_rsp;       // 输入：是否启用RSP功能(0/1)
} hifiasm_opt_t;
```

## 核心函数实现

### 1. RSP段构建函数
```cpp
// 位置：rsp_integration.cpp
// 作用：从参考比对信息解析并构建RSP段
// 输入：segments(输出容器), g(overlap图)
// 输出：成功构建的段数量
int rsp_build_segments(rsp_segments_t *segments, asg_t *g);
```

### 2. RSP预标记函数
```cpp
// 位置：rsp_integration.cpp  
// 作用：标记RSP段内部reads，防止标准算法错误处理
// 输入：segments(RSP段), mark(标记数组)
// 输出：void，通过mark数组标记内部reads
void rsp_pre_mark_internal_reads(rsp_segments_t *segments, int32_t *mark);
```

### 3. RSP预构建函数
```cpp
// 位置：rsp_integration.cpp
// 作用：根据RSP段预先构建unitigs
// 输入：segments(RSP段), ug(unitig图), g(overlap图)  
// 输出：void，向ug中添加RSP unitigs
void rsp_pre_build_unitigs(rsp_segments_t *segments, ma_ug_t *ug, asg_t *g);
```

### 4. RSP端点扩展函数
```cpp
// 位置：rsp_integration.cpp
// 作用：检测并双向扩展RSP unitigs到非参考区域
// 输入：ug(unitig图), g(overlap图), read_id(端点read)
// 输出：void，就地扩展RSP unitigs
// 重要：支持头端点向上游扩展和尾端点向下游扩展
void extend_rsp_unitigs_if_possible(ma_ug_t *ug, asg_t *g, uint32_t read_id) {
    rsp_endpoint_info_t *info = get_rsp_endpoint_info(read_id);
    if (!info) return;
    
    uint32_t unitig_id = info->segment_id;
    
    if (info->is_tail == 0) { // 头端点，检查向上游扩展
        uint32_t vertex = read_id << 1;  // 正向顶点
        if (arc_cnt(g, vertex^1) == 1) {  // 入度=1，满足扩展条件
            extend_rsp_unitig_backward(ug, unitig_id, g, read_id);
        }
    } else { // 尾端点，检查向下游扩展
        uint32_t vertex = read_id << 1;  // 正向顶点  
        if (arc_cnt(g, vertex) == 1) {   // 出度=1，满足扩展条件
            extend_rsp_unitig_forward(ug, unitig_id, g, read_id);
        }
    }
}
```

### 5. RSP端点检测函数
```cpp
// 位置：rsp_integration.cpp
// 作用：判断给定read是否为RSP unitig的端点
// 输入：read_id(read标识符)
// 输出：true表示是端点，false表示不是
bool is_rsp_endpoint(uint32_t read_id);
```

## 修改现有函数

### ma_ug_gen_primary函数修改
```cpp
// 位置：Assembly.cpp中的ma_ug_gen_primary函数
// 修改内容：在原函数基础上集成RSP预处理和扩展逻辑

ma_ug_t *ma_ug_gen_primary(asg_t *g, uint8_t flag)
{
    // ========== 现有代码保持不变 ==========
    asg_cleanup(g);
    int32_t *mark;
    uint32_t i, v, n_vtx = g->n_seq * 2;
    kdq_t(uint64_t) *q;
    ma_ug_t *ug;

    ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
    ug->g = asg_init();
    mark = (int32_t*)calloc(n_vtx, 4);

    // ========== 新增：RSP预处理步骤 ==========
    rsp_segments_t rsp_segments;
    memset(&rsp_segments, 0, sizeof(rsp_segments));
    
    // 检查是否启用RSP功能
    if (asm_opt.ref_alignments && asm_opt.ref_alignments[0]) {
        fprintf(stderr, "[RSP] Building RSP segments from %s...\n", asm_opt.ref_alignments);
        
        // 构建RSP段
        if (rsp_build_segments(&rsp_segments, g) > 0) {
            // 预标记RSP段内部reads
            rsp_pre_mark_internal_reads(&rsp_segments, mark);
            
            // 预构建RSP段为unitigs
            rsp_pre_build_unitigs(&rsp_segments, ug, g);
            
            // 初始化端点追踪
            init_rsp_endpoint_tracking(&rsp_segments);
            
            fprintf(stderr, "[RSP] Pre-built %u RSP unitigs from %u segments\n", 
                    (uint32_t)ug->u.n, rsp_segments.n_segments);
        }
    }

    // ========== 修改：标准算法集成RSP扩展检测 ==========
    q = kdq_init(uint64_t);
    for (v = 0; v < n_vtx; ++v) {
        uint32_t w, x, l, start, end, len;
        ma_utg_t *p;
        
        // 原有过滤条件保持不变
        if (g->seq[v>>1].del || mark[v]) continue;
        if (arc_cnt(g, v) == 0 && arc_cnt(g, (v^1)) != 0) continue;
        if(flag == PRIMARY_LABLE && g->seq[v>>1].c == ALTER_LABLE) continue;
        if(flag == ALTER_LABLE && g->seq[v>>1].c != ALTER_LABLE) continue;
        
        // 新增：检查是否为RSP端点，如果是则扩展
        uint32_t read_id = v >> 1;
        if (is_rsp_endpoint(read_id)) {
            extend_rsp_unitigs_if_possible(ug, g, read_id);
            mark[v] = 1;  // 标记为已处理
            continue;
        }

        // ========== 原有标准算法逻辑保持不变 ==========
        mark[v] = 1;
        q->count = 0, start = v, end = v^1, len = 0;
        
        // forward延伸逻辑不变
        w = v;
        while (1) {
            if (arc_cnt(g, w) != 1) break;
            x = arc_first(g, w).v;
            if (arc_cnt(g, x^1) != 1) break;
            mark[x] = mark[w^1] = 1;
            l = asg_arc_len(arc_first(g, w));
            kdq_push(uint64_t, q, (uint64_t)w<<32 | l);
            end = x^1, len += l;
            w = x;
            if (x == v) break;
        }
        
        // 线性/环形unitig处理逻辑不变
        if (start != (end^1) || kdq_size(q) == 0) {
            l = g->seq[end>>1].len;
            kdq_push(uint64_t, q, (uint64_t)(end^1)<<32 | l);
            len += l;
        } else {
            start = end = UINT32_MAX;
            goto add_unitig;
        }
        
        // backward延伸逻辑不变
        x = v;
        while (1) {
            if (arc_cnt(g, x^1) != 1) break;
            w = arc_first(g, x^1).v ^ 1;
            if (arc_cnt(g, w) != 1) break;
            mark[x] = mark[w^1] = 1;
            l = asg_arc_len(arc_first(g, w));
            kdq_unshift(uint64_t, q, (uint64_t)w<<32 | l);
            start = w, len += l;
            x = w;
        }

add_unitig:
        // unitig创建逻辑不变
        if (start != UINT32_MAX) mark[start] = mark[end] = 1;
        kv_pushp(ma_utg_t, ug->u, &p);
        p->s = 0, p->start = start, p->end = end, p->len = len, p->n = kdq_size(q), p->circ = (start == UINT32_MAX);
        p->m = p->n;
        kv_roundup32(p->m);
        p->a = (uint64_t*)malloc(8 * p->m);
        for (i = 0; i < kdq_size(q); ++i)
            p->a[i] = kdq_at(q, i);
    }
    kdq_destroy(uint64_t, q);

    // ========== 现有unitig图构建逻辑保持不变 ==========
    for (v = 0; v < n_vtx; ++v) mark[v] = -1;
    for (i = 0; i < ug->u.n; ++i) {
        if (ug->u.a[i].circ) continue;
        mark[ug->u.a[i].start] = i<<1 | 0;
        mark[ug->u.a[i].end] = i<<1 | 1;
    }
    
    for (i = 0; i < g->n_arc; ++i) {
        asg_arc_t *p = &g->arc[i];
        if (p->del) continue;
        if (mark[p->ul>>32^1] >= 0 && mark[p->v] >= 0) {
            asg_arc_t *q;
            uint32_t u = mark[p->ul>>32^1]^1;
            int l = ug->u.a[u>>1].len - p->ol;
            if (l < 0) l = 1;
            q = asg_arc_pushp(ug->g);
            q->ol = p->ol, q->del = 0;
            q->ul = (uint64_t)u<<32 | l;
            q->v = mark[p->v]; q->ou = 0;
            q->el = p->el;
        }
    }
    
    for (i = 0; i < ug->u.n; ++i)
        asg_seq_set(ug->g, i, ug->u.a[i].len, 0);
    asg_cleanup(ug->g);
    
    // ========== 新增：RSP清理 ==========
    if (rsp_segments.n_segments > 0) {
        rsp_cleanup_segments(&rsp_segments);
        cleanup_rsp_endpoint_tracking();
    }
    
    free(mark);
    return ug;
}
```

## PAF信息集成 - 复用现有hifiasm模块

### 重用现有的比对数据结构和函数

根据项目知识分析，hifiasm已有完整的比对处理系统，我们可以直接复用：

#### 1. 重用ma_hit_t结构存储PAF信息

**现有的ma_hit_t结构（来自Overlaps.h）**：
```cpp
// 位置：Overlaps.h - 已存在，无需新增
typedef struct {
    uint64_t qns;              // query name + start position
    uint32_t qe, tn, ts, te;   // query end, target name, target start/end
    uint32_t cc:30, ml:1, rev:1; // coverage, match length, reverse
    uint32_t bl:31, del:1;     // block length, delete flag
    uint8_t el;                // exact length
    uint8_t no_l_indel;        // no large indel
} ma_hit_t;

// 比对分配结构 - 已存在
typedef struct {
    ma_hit_t* buffer;          // 比对记录缓冲区
    uint32_t size;             // 缓冲区大小
    uint32_t length;           // 当前记录数量
    uint8_t is_fully_corrected;
    uint8_t is_abnormal;
} ma_hit_t_alloc;
```

#### 2. 重用现有的比对处理函数

**比对分配和处理函数（已存在）**：
```cpp
// 位置：Assembly.cpp - 已存在，直接复用
void add_ma_hit_t_alloc(ma_hit_t_alloc* x, ma_hit_t* element);  // 添加比对记录
void ma_hit_sort_tn(ma_hit_t *a, long long n);                 // 按target排序
void ma_hit_sort_qns(ma_hit_t *a, long long n);                // 按query排序

// 比对访问宏 - 已存在，直接复用
#define Get_qn(hit) ((uint32_t)((hit).qns>>32))                // 获取query ID
#define Get_qs(hit) ((uint32_t)(hit).qns)                      // 获取query start
#define Get_qe(hit) ((hit).qe)                                 // 获取query end
#define Get_tn(hit) ((hit).tn)                                 // 获取target ID
#define Get_ts(hit) ((hit).ts)                                 // 获取target start
#define Get_te(hit) ((hit).te)                                 // 获取target end
```

### RSP段构建的具体实现（复用现有系统）

```cpp
// 位置：rsp_integration.cpp
// 作用：将PAF信息转换为ma_hit_t格式并构建RSP段
// 重用：ma_hit_t_alloc, add_ma_hit_t_alloc, ma_hit_sort_tn等现有函数
int rsp_build_segments(rsp_segments_t *segments, asg_t *g) {
    if (!asm_opt.ref_alignments || !asm_opt.ref_alignments[0]) {
        segments->segments = NULL;
        segments->n_segments = 0;
        return 0;
    }
    
    // 1. 重用现有的ma_hit_t_alloc系统加载PAF
    ma_hit_t_alloc rsp_alignments;
    init_ma_hit_t_alloc(&rsp_alignments);  // 重用现有函数
    
    if (rsp_load_paf_to_ma_hit_t(asm_opt.ref_alignments, &rsp_alignments) <= 0) {
        fprintf(stderr, "[RSP] Failed to load PAF file: %s\n", asm_opt.ref_alignments);
        clear_ma_hit_t_alloc(&rsp_alignments);  // 重用现有函数
        return 0;
    }
    
    // 2. 重用现有的排序和过滤函数
    ma_hit_sort_tn(rsp_alignments.buffer, rsp_alignments.length);  // 按参考序列排序
    rsp_filter_alignments_by_mapq(&rsp_alignments, asm_opt.rsp_min_mapq);
    
    // 3. 基于已排序的ma_hit_t数据构建RSP段
    rsp_identify_segments_from_ma_hit_t(&rsp_alignments, segments, asm_opt.rsp_max_gap);
    
    // 4. 验证段内reads在overlap图中的连通性（重用asg_arc_n等现有函数）
    rsp_validate_segments_connectivity(segments, g);
    
    clear_ma_hit_t_alloc(&rsp_alignments);  // 重用现有函数
    return segments->n_segments;
}
```

#### 3. PAF到ma_hit_t的转换函数

```cpp
// 位置：rsp_integration.cpp
// 作用：将PAF格式转换为hifiasm内部的ma_hit_t格式
// 输入：paf_file(PAF文件路径), alns(输出的ma_hit_t_alloc结构)
// 输出：成功加载的比对记录数量
int rsp_load_paf_to_ma_hit_t(const char *paf_file, ma_hit_t_alloc *alns) {
    FILE *fp = fopen(paf_file, "r");
    if (!fp) {
        fprintf(stderr, "[RSP] Cannot open PAF file: %s\n", paf_file);
        return -1;
    }
    
    char line[4096];
    int count = 0;
    
    while (fgets(line, sizeof(line), fp)) {
        ma_hit_t hit;
        if (rsp_parse_paf_line_to_ma_hit_t(line, &hit) == 0) {
            add_ma_hit_t_alloc(alns, &hit);  // 重用现有函数
            count++;
        }
    }
    
    fclose(fp);
    fprintf(stderr, "[RSP] Loaded %d alignments from %s\n", count, paf_file);
    return count;
}

// PAF行解析到ma_hit_t格式
int rsp_parse_paf_line_to_ma_hit_t(char *line, ma_hit_t *hit) {
    char read_name[256], ref_name[256], strand[8];
    uint32_t read_len, read_start, read_end, ref_len, ref_start, ref_end;
    uint32_t matches, align_len, mapq = 0;
    
    // 解析PAF格式：qname qlen qstart qend strand tname tlen tstart tend matches alen mapq
    int n = sscanf(line, "%s %u %u %u %s %s %u %u %u %u %u %u",
                   read_name, &read_len, &read_start, &read_end, strand,
                   ref_name, &ref_len, &ref_start, &ref_end, &matches, &align_len, &mapq);
    
    if (n < 11) return -1;  // PAF格式错误
    
    // 过滤低质量比对
    if (mapq < asm_opt.rsp_min_mapq) return -1;
    
    // 转换为ma_hit_t格式
    uint32_t read_id = rsp_get_read_id_by_name(read_name);  // 获取内部read ID
    if (read_id == (uint32_t)-1) return -1;  // read不存在
    
    // 构造ma_hit_t结构（复用现有格式）
    hit->qns = ((uint64_t)read_id << 32) | read_start;
    hit->qe = read_end;
    hit->tn = rsp_get_ref_id_by_name(ref_name);  // 参考序列ID
    hit->ts = ref_start;
    hit->te = ref_end;
    hit->rev = (strand[0] == '-') ? 1 : 0;
    hit->ml = 1;  // 标记为参考比对
    hit->el = 1;  // 精确比对
    hit->del = 0; // 未删除
    hit->cc = mapq; // 使用cc字段存储映射质量
    
    return 0;
}
```

#### 4. 基于ma_hit_t构建RSP段

```cpp
// 位置：rsp_integration.cpp  
// 作用：从已排序的ma_hit_t数据识别连续的RSP段
// 输入：alns(已排序的比对数据), segments(输出段), max_gap(最大间隙)
void rsp_identify_segments_from_ma_hit_t(ma_hit_t_alloc *alns, rsp_segments_t *segments, int64_t max_gap) {
    if (alns->length == 0) return;
    
    // 按参考序列和位置分组
    uint32_t current_ref = Get_tn(alns->buffer[0]);
    uint32_t seg_start = 0;
    kvec_t(uint32_t) seg_reads;
    kv_init(seg_reads);
    
    for (uint32_t i = 0; i < alns->length; i++) {
        uint32_t ref_id = Get_tn(alns->buffer[i]);
        uint32_t read_id = Get_qn(alns->buffer[i]);
        
        // 检查是否需要开始新段
        if (ref_id != current_ref || 
            (i > 0 && Get_ts(alns->buffer[i]) - Get_te(alns->buffer[i-1]) > max_gap)) {
            
            // 完成当前段
            if (seg_reads.n >= 2) {  // 至少2个reads才构成段
                rsp_create_segment_from_reads(&seg_reads, alns, segments, current_ref);
            }
            
            // 开始新段
            seg_reads.n = 0;
            current_ref = ref_id;
        }
        
        kv_push(uint32_t, seg_reads, i);  // 添加read索引
    }
    
    // 处理最后一个段
    if (seg_reads.n >= 2) {
        rsp_create_segment_from_reads(&seg_reads, alns, segments, current_ref);
    }
    
    kv_destroy(seg_reads);
}
```

### 利用现有的Read管理系统

```cpp
// 位置：rsp_integration.cpp
// 作用：将PAF中的read名称映射到hifiasm内部的read ID
// 重用：现有的R_INF (All_reads)系统和Get_NAME宏
uint32_t rsp_get_read_id_by_name(const char *read_name) {
    // 重用现有的All_reads (R_INF)系统进行read查找
    for (uint32_t i = 0; i < R_INF.total_reads; i++) {
        if (strcmp(Get_NAME(R_INF, i), read_name) == 0) {  // 重用现有宏
            return i;
        }
    }
    return (uint32_t)-1;
}
```

### 过滤函数重用现有逻辑

```cpp
// 位置：rsp_integration.cpp
// 作用：基于映射质量过滤比对记录  
// 重用：现有的ma_hit_t结构和del字段标记删除
void rsp_filter_alignments_by_mapq(ma_hit_t_alloc *alns, int32_t min_mapq) {
    uint32_t keep_count = 0;
    
    for (uint32_t i = 0; i < alns->length; i++) {
        // 检查映射质量（存储在cc字段中）
        if (alns->buffer[i].cc >= min_mapq) {
            if (keep_count != i) {
                alns->buffer[keep_count] = alns->buffer[i];
            }
            keep_count++;
        }
    }
    
    alns->length = keep_count;
    fprintf(stderr, "[RSP] Filtered to %u high-quality alignments (MAPQ >= %d)\n", 
            keep_count, min_mapq);
}
```

## 关键优势

### 1. 最大化代码重用
- **ma_hit_t结构**：直接使用现有的比对记录格式
- **内存管理**：重用`ma_hit_t_alloc`的分配和释放函数
- **排序过滤**：重用现有的`ma_hit_sort_tn`、`ma_hit_sort_qns`等函数
- **Read管理**：重用`R_INF (All_reads)`系统和相关宏

### 2. 保持一致性
- **数据格式**：PAF信息转换为hifiasm内部标准格式
- **内存模式**：使用相同的内存分配和释放模式
- **错误处理**：遵循现有的错误处理约定

### 3. 性能优化
- **零拷贝**：直接在ma_hit_t格式上操作，避免格式转换开销
- **批量处理**：利用现有的批量排序和过滤函数
- **索引重用**：利用现有的read名称到ID的映射机制

### 命令行参数集成

#### 修改CommandLines.h和CommandLines.cpp

**在hifiasm_opt_t中添加**：
```cpp
// 位置：CommandLines.h
typedef struct {
    // ... 现有字段 ...
    
    // RSP相关字段
    char *ref_alignments;     // PAF文件路径
    int32_t rsp_min_mapq;     // 最小映射质量
    int32_t rsp_max_gap;      // 最大间隙大小
    uint8_t enable_rsp;       // 启用RSP功能
} hifiasm_opt_t;
```

**在CommandLines.cpp中添加**：
```cpp
// 位置：CommandLines.cpp的CommandLine_process函数
static const ketopt_t long_options[] = {
    // ... 现有选项 ...
    { "ref-alignments", ko_required_argument, 400 },
    { "rsp-min-mapq", ko_required_argument, 401 },
    { "rsp-max-gap", ko_required_argument, 402 },
    { 0, 0, 0 }
};

// 在参数处理循环中添加
case 400: // --ref-alignments
    asm_opt->ref_alignments = opt.arg;
    asm_opt->enable_rsp = 1;
    break;
case 401: // --rsp-min-mapq  
    asm_opt->rsp_min_mapq = atoi(opt.arg);
    break;
case 402: // --rsp-max-gap
    asm_opt->rsp_max_gap = atoi(opt.arg);
    break;
```

**在init_opt函数中初始化**：
```cpp
// 位置：CommandLines.cpp的init_opt函数
void init_opt(hifiasm_opt_t* asm_opt) {
    // ... 现有初始化 ...
    
    // RSP参数初始化
    asm_opt->ref_alignments = NULL;
    asm_opt->rsp_min_mapq = 10;       // 默认最小映射质量
    asm_opt->rsp_max_gap = 50000;     // 默认最大间隙50kb
    asm_opt->enable_rsp = 0;          // 默认关闭RSP
}
```

### 使用示例和工作流

#### 完整的RSP工作流

```bash
# 1. 准备参考基因组和HiFi reads
minimap2 -ax map-hifi -t 32 reference.fa hifi_reads.fq.gz > alignments.paf

# 2. 运行带RSP的hifiasm
./hifiasm -o output --ref-alignments alignments.paf --rsp-min-mapq 20 --rsp-max-gap 100000 -t 32 hifi_reads.fq.gz

# 3. 预期日志输出
[RSP] Loading PAF alignments from alignments.paf...
[RSP] Loaded 45000 alignments, filtered to 42000 (MAPQ >= 20)
[RSP] Built 120 RSP segments covering 38000 reads
[RSP] Pre-marked 35000 internal reads from 120 segments
[RSP] Pre-built 120 RSP unitigs from 120 segments
[RSP] Extended 45 RSP unitigs, merged 23 segments
[RSP] Final result: 89 unitigs (67 RSP-derived, 22 standard)
```

## 预期结果验证

### 目标输出验证
根据项目需求，实现后应产生以下unitigs：

```
最终unitigs（ug->u.n = 2）：

扩展RSP Unitig 0: Read0 → Read1 → Read2 → Read3 → Read6 → Read7 → Read8 → Read4 → Read5
- 起始：RSP段 Read0→Read1→Read2→Read3（Chr1前段）
- 扩展：Read3端点检测到连接Read6，向前扩展
- 继续：Read6→Read7→Read8（无参考区域）  
- 合并：Read8→Read4→Read5（Chr1后段，覆盖原RSP Unitig 1）
- 特点：跨越两个Chr1参考区域和中间无参考区域

分支 Unitig 1: Read9 → Read10
- 类型：标准算法处理的独立分支
- 原因：Read0和Read3不满足单入单出条件，分支无法连接主链
- 特点：无参考信息的独立路径
```

### PAF文件示例

**测试用PAF文件格式**：
```
read_000    12000   0     11800   +   chr1    248956422   1000    12800   11500   11800   60
read_001    12200   0     12000   +   chr1    248956422   5000    17000   11800   12000   60  
read_002    11800   0     11600   +   chr1    248956422   10000   21600   11400   11600   60
read_003    12500   0     12300   +   chr1    248956422   16000   28300   12100   12300   60
read_004    13000   0     12800   -   chr1    248956422   100000  112800  12600   12800   60
read_005    12800   0     12600   +   chr1    248956422   129000  141600  12400   12600   60
```

这个PAF文件将指导RSP系统：
- Read0-3形成Chr1前段的RSP段
- Read4-5形成Chr1后段的RSP段（Read4反向比对）
- Read6-10无参考信息，由标准算法处理

### 关键实现点
1. **RSP段预构建**：Read0-3和Read4-5分别构建为初始RSP unitigs
2. **内部reads标记**：Read1和Read2被标记，防止错误起点
3. **端点扩展检测（双向）**：
   - **Read0（头端点）**：可以向上游扩展，如果存在Read10→Read0连接且满足单入单出条件
   - **Read3（尾端点）**：可以向下游扩展，如果存在Read3→Read6连接且满足单入单出条件
   - **扩展条件**：关键是满足单入单出条件，而不是端点类型
4. **跨区域合并**：通过Read6-8连接两个Chr1区域
5. **分支独立处理**：Read9-10因连接条件不满足而独立成unitig

## 向后兼容性保证

### 不提供PAF文件时的正常运行

**关键设计原则**：RSP功能完全可选，不影响原有hifiasm功能。

```cpp
// 位置：Assembly.cpp中的ma_ug_gen_primary函数
ma_ug_t *ma_ug_gen_primary(asg_t *g, uint8_t flag)
{
    // ========== 现有代码保持不变 ==========
    asg_cleanup(g);
    int32_t *mark;
    uint32_t i, v, n_vtx = g->n_seq * 2;
    kdq_t(uint64_t) *q;
    ma_ug_t *ug;

    ug = (ma_ug_t*)calloc(1, sizeof(ma_ug_t));
    ug->g = asg_init();
    mark = (int32_t*)calloc(n_vtx, 4);

    // ========== RSP预处理（仅在提供PAF时执行）==========
    rsp_segments_t rsp_segments;
    memset(&rsp_segments, 0, sizeof(rsp_segments));
    
    // 关键检查：只有在提供参考比对文件时才启用RSP
    if (asm_opt.ref_alignments && asm_opt.ref_alignments[0]) {
        fprintf(stderr, "[RSP] Building RSP segments from %s...\n", asm_opt.ref_alignments);
        
        if (rsp_build_segments(&rsp_segments, g) > 0) {
            rsp_pre_mark_internal_reads(&rsp_segments, mark);
            rsp_pre_build_unitigs(&rsp_segments, ug, g);
            init_rsp_endpoint_tracking(&rsp_segments);
            fprintf(stderr, "[RSP] Pre-built %u RSP unitigs\n", (uint32_t)ug->u.n);
        }
    }
    // 如果没有提供PAF文件，rsp_segments.n_segments = 0，不会执行任何RSP逻辑

    // ========== 标准算法（与原版本完全一致）==========
    q = kdq_init(uint64_t);
    for (v = 0; v < n_vtx; ++v) {
        uint32_t w, x, l, start, end, len;
        ma_utg_t *p;
        
        // 原有过滤条件保持不变
        if (g->seq[v>>1].del || mark[v]) continue;
        if (arc_cnt(g, v) == 0 && arc_cnt(g, (v^1)) != 0) continue;
        if(flag == PRIMARY_LABLE && g->seq[v>>1].c == ALTER_LABLE) continue;
        if(flag == ALTER_LABLE && g->seq[v>>1].c != ALTER_LABLE) continue;
        
        // RSP端点检查（仅在有RSP段时执行）
        uint32_t read_id = v >> 1;
        if (rsp_segments.n_segments > 0 && is_rsp_endpoint(read_id)) {
            extend_rsp_unitigs_if_possible(ug, g, read_id);
            mark[v] = 1;
            continue;
        }

        // ========== 原有标准算法逻辑完全不变 ==========
        mark[v] = 1;
        q->count = 0, start = v, end = v^1, len = 0;
        
        // forward延伸逻辑...（保持原样）
        // backward延伸逻辑...（保持原样）
        // unitig创建逻辑...（保持原样）
    }
    
    // ========== 原有unitig图构建逻辑保持不变 ==========
    // ...所有原有代码保持不变...
    
    // ========== RSP清理（仅在有RSP段时执行）==========
    if (rsp_segments.n_segments > 0) {
        rsp_cleanup_segments(&rsp_segments);
        cleanup_rsp_endpoint_tracking();
    }
    
    free(mark);
    return ug;
}
```

### 兼容性验证测试

```bash
# 1. 原有功能测试（无RSP）
./hifiasm -o test_original -t 4 reads.fq.gz
# 预期：与原版本行为完全一致，无RSP相关日志

# 2. RSP功能测试（有PAF文件）
./hifiasm -o test_rsp --ref-alignments alignments.paf -t 4 reads.fq.gz
# 预期：额外的RSP日志，但核心算法保持一致

# 3. 结果对比验证
# 在相同reads输入下：
# - 无PAF：产生标准unitig图
# - 有PAF：产生RSP优化的unitig图，但分支结构保持逻辑一致性
```

### 关键安全机制

1. **条件编译**：RSP代码仅在检测到PAF文件时执行
2. **内存安全**：RSP相关内存仅在需要时分配，总是正确释放
3. **状态隔离**：RSP状态不影响原有算法的任何全局变量
4. **错误恢复**：PAF文件解析失败时，回退到标准算法

## 项目交付清单

### 新增文件
1. **rsp_integration.h** - RSP数据结构和函数声明
2. **rsp_integration.cpp** - RSP核心功能实现

### 修改文件  
1. **CommandLines.h** - 添加RSP配置字段到hifiasm_opt_t
2. **CommandLines.cpp** - 添加RSP命令行参数处理
3. **Assembly.cpp** - 修改ma_ug_gen_primary函数集成RSP

### 测试文件
1. **test.paf** - 测试用参考比对文件
2. **test_reads.fa** - 测试用reads文件
3. **validate_rsp.sh** - RSP功能验证脚本

## 重要实现注意事项

1. **内存管理**：所有RSP相关内存必须在函数结束前正确释放
2. **线程安全**：RSP功能必须与现有多线程框架兼容
3. **错误处理**：参考比对文件格式错误时的容错处理
4. **性能优化**：RSP预处理不应显著增加总体组装时间
5. **向后兼容**：不启用RSP时，行为与原版本完全一致

通过以上设计，RSP集成将实现智能的reads组装和unitig扩展，既保证参考区域的正确组装，又允许跨区域的自然连接，最终产生更长、更准确的unitigs。
