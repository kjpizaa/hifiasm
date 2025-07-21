#ifndef __INTER__
#define __INTER__
#include "Overlaps.h"
#include "Process_Read.h"
#include "hic.h"

#define G_CHAIN_BW 16//128
#define FLANK_M (0x7fffU)
#define P_CHAIN_COV 0.985
#define P_FRAGEMENT_CHAIN_COV 0.20
#define P_FRAGEMENT_PRIMARY_CHAIN_COV 0.70
#define P_FRAGEMENT_PRIMARY_SECOND_COV 0.25
#define P_CHAIN_SCORE 0.6
#define G_CHAIN_GAP 0.1
#define UG_SKIP 5
#define RG_SKIP 25
#define UG_SKIP_GRAPH_N 72
#define UG_SKIP_N 100
#define UG_ITER_N 5000
#define UG_DIS_N 50000
// #define UG_TRANS_W 2
#define UG_TRANS_W 2
// #define UG_TRANS_ERR_W 512
#define UG_TRANS_ERR_W 64
#define G_CHAIN_TRANS_RATE 0.25
#define G_CHAIN_TRANS_WEIGHT -1
#define G_CHAIN_INDEL 128
#define W_CHN_PEN_GAP 0.1
#define N_GCHAIN_RATE 0.04
#define PRIMARY_UL_CHAIN_MIN 75000

typedef struct {
	int w, k, bw, max_gap, is_HPC, hap_n, occ_weight, max_gap_pre, max_gc_seq_ext, seed;
	int max_lc_skip, max_lc_iter, min_lc_cnt, min_lc_score, max_gc_skip, ref_bonus;
	int min_gc_cnt, min_gc_score, sub_diff, best_n;
    float chn_pen_gap, mask_level, pri_ratio;
	///base-alignment
	double bw_thres, diff_ec_ul, diff_ec_ul_low, diff_ec_ul_hpc; int max_n_chain, ec_ul_round;
} mg_idxopt_t;

struct mg_tbuf_s {
	void *km;
	int frag_gap;
};
typedef struct mg_tbuf_s mg_tbuf_t;


mg_tbuf_t *mg_tbuf_init(void);

void mg_tbuf_destroy(mg_tbuf_t *b);

void *mg_tbuf_get_km(mg_tbuf_t *b);

typedef struct {
	FILE *fp;
	ul_vec_t u;
	uint64_t flag;
} ucr_file_t;

typedef struct {
	int32_t off, cnt;
	uint32_t v;
	int32_t score;
} mg_llchain_t;

typedef struct {
	int32_t id, parent;
	int32_t off, cnt;
	int32_t n_anchor, score;
	int32_t qs, qe;
	int32_t plen, ps, pe;
	int32_t blen, mlen;
	float div;
	uint32_t hash;
	int32_t subsc, n_sub;
	uint32_t mapq:8, flt:1, dummy:23;
} mg_gchain_t;

typedef struct {
	size_t n,m;
    uint64_t *a, tl;
	kvec_t(char) cc;
} mg_dbn_t;

typedef struct {
	int32_t cnt;
	uint32_t v;
	int32_t score;
	uint32_t qs, qe, ts, te;
} mg_lres_t;

typedef struct {
    int32_t n_gc, n_lc;
    mg_gchain_t *gc;///g_chain; idx in l_chains
    mg_lres_t *lc;///l_chain
    uint64_t qid, qlen;
} mg_gres_t;

typedef struct {
	size_t n,m;
    mg_gres_t *a;
	uint64_t total_base;
    uint64_t total_pair;
} mg_gres_a;

void push_uc_block_t(const ug_opt_t *uopt, kv_ul_ov_t *z, char **seq, uint64_t *len, uint64_t b_id);
void ul_resolve(ma_ug_t *ug, const asg_t *rg, const ug_opt_t *uopt, int hap_n);
void ul_load(const ug_opt_t *uopt);
uint64_t* get_hifi2ul_list(all_ul_t *x, uint64_t hid, uint64_t* a_n);
uint64_t ul_refine_alignment(const ug_opt_t *uopt, asg_t *sg);
ma_ug_t *ul_realignment(const ug_opt_t *uopt, asg_t *sg, uint32_t double_check_cache, const char *bin_file);
int32_t write_all_ul_t(all_ul_t *x, char* file_name, ma_ug_t *ug);
int32_t load_all_ul_t(all_ul_t *x, char* file_name, All_reads *hR, ma_ug_t *ug);
uint32_t ugl_cover_check(uint64_t is, uint64_t ie, ma_utg_t *u);
void filter_ul_ug(ma_ug_t *ug);
void gen_ul_vec_rid_t(all_ul_t *x, All_reads *rdb, ma_ug_t *ug);
void update_ug_arch_ul_mul(ma_ug_t *ug);
void print_ul_alignment(ma_ug_t *ug, all_ul_t *aln, uint32_t id, const char* cmd);
void clear_all_ul_t(all_ul_t *x);
void trans_base_infer(ma_ug_t *ug, asg_t *sg, ug_opt_t *uopt, kv_u_trans_t *res, bubble_type *bub);
hpc_re_t *gen_hpc_re_t(ma_ug_t *ug);
idx_emask_t* graph_ovlp_binning(ma_ug_t *ug, asg_t *sg, const ug_opt_t *uopt);
uint32_t gen_src_shared_interval_simple(uint32_t src, ma_ug_t *ug, uint64_t *flt, uint64_t flt_n, kv_ul_ov_t *res);
uint64_t check_ul_ov_t_consist(ul_ov_t *x, ul_ov_t *y, int64_t ql, int64_t tl, double diff);
uint32_t infer_se(uint32_t qs, uint32_t qe, uint32_t ts, uint32_t te, uint32_t rev, 
uint32_t rqs, uint32_t rqe, uint32_t *rts, uint32_t *rte);
uint32_t clean_contain_g(const ug_opt_t *uopt, asg_t *sg, uint32_t push_trans);
void dedup_contain_g(const ug_opt_t *uopt, asg_t *sg);
void trans_base_mmhap_infer(ma_ug_t *ug, asg_t *sg, ug_opt_t *uopt, kv_u_trans_t *res);
scaf_res_t *gen_contig_path(const ug_opt_t *uopt, asg_t *sg, ma_ug_t *ctg, ma_ug_t *ref);
void gen_contig_trans(const ug_opt_t *uopt, asg_t *sg, ma_ug_t *qry, scaf_res_t *qry_sc, ma_ug_t *ref, scaf_res_t *ref_sc, ma_ug_t *gfa, kv_u_trans_t *ta, uint32_t qoff, uint32_t toff, bubble_type *bu, kv_u_trans_t *res);
void gen_contig_self(const ug_opt_t *uopt, asg_t *sg, ma_ug_t *db, scaf_res_t *db_sc, ma_ug_t *gfa, kv_u_trans_t *ta, uint64_t soff, bubble_type *bu, kv_u_trans_t *res, uint32_t is_exact);
void order_contig_trans(kv_u_trans_t *in);
void sort_uc_block_qe(uc_block_t* a, uint64_t a_n);
// inter.h 文件末尾添加的完整代码块

#ifdef ENABLE_REF_GENOME_V4
#include "ref_genome.h"
#include "Hash_Table.h"
// 参考基因组Block标记宏 (uc_block_t.el 高位)
#define BLOCK_REF             (1u<<15)
#define BLOCK_SET_REF(block)   ((block)->el |= BLOCK_REF)
#define BLOCK_IS_REF(block)    ((block)->el & BLOCK_REF)
#define BLOCK_CLEAR_REF(block) ((block)->el &= ~BLOCK_REF)

// 🔧 修正：移除全局变量声明（因为已经在Assembly.cpp中有了）
// extern ma_ug_t *ug; /* global unitig graph for reference-guided pipeline */

// 🔧 修正：函数声明需要正确的条件编译包围
const char* ensure_unitig_seq(ma_ug_t* ug, uint32_t uid);

int overlap_to_uc_block_ref_mode(overlap_region_alloc *overlap_list,
                                uint32_t query_unitig_id,
                                uc_block_t **out_blocks,
                                uint64_t *out_count);

int unitigs_map_to_reference_batch(ma_ug_t *unitigs,
                                  const ul_idx_t *ref_index,
                                  uc_block_t **out_blocks,
                                  uint64_t *out_count,
                                  const hifiasm_opt_t *opt);

int integrate_reference_blocks_to_existing_ul_pipeline(ma_ug_t *unitigs,
                                                      const ul_idx_t *ref_index,
                                                      const hifiasm_opt_t *opt);

#endif // ENABLE_REF_GENOME_V4
#ifdef ENABLE_REF_GENOME_V4

// Ensure unitig sequence is allocated. When missing, allocate an N-filled
// string of the desired length. This is a lightweight fallback used when
// the sequence is not generated during earlier steps.
const char* ensure_unitig_seq(ma_ug_t* ug, uint32_t uid)
{
    if (!ug || uid >= ug->u.n) return NULL;
    ma_utg_t *u = &ug->u.a[uid];
    if (u->s) return u->s;
    u->s = (char*)calloc(u->len + 1, 1);
    if (!u->s) return NULL;
    memset(u->s, 'N', u->len);
    u->s[u->len] = '\0';
    return u->s;
}

// Convert overlap regions to uc_block_t list. Only long and high-quality
// overlaps are retained. BLOCK_REF is set on all generated blocks.
int overlap_to_uc_block_ref_mode(overlap_region_alloc *overlap_list,
                                uint32_t query_unitig_id,
                                uc_block_t **out_blocks,
                                uint64_t *out_count)
{
    if (!overlap_list || !out_blocks || !out_count) return -1;
    *out_blocks = NULL; *out_count = 0;
    if (overlap_list->length == 0) return 0;

    uc_block_t *res = (uc_block_t*)calloc(overlap_list->length, sizeof(uc_block_t));
    if (!res) return -1;

    uint64_t cnt = 0;
    for (uint64_t i = 0; i < overlap_list->length; i++) {
        overlap_region *ov = &overlap_list->list[i];
        if (ov->overlapLen < 500) continue;
        double er = 1.0 - (double)ov->align_length / (double)(ov->overlapLen?ov->overlapLen:1);
        if (er > 0.15) continue;

        uc_block_t *b = &res[cnt++];
        b->hid = ov->y_id;      // reference id
        b->qs = ov->x_pos_s; b->qe = ov->x_pos_e;
        b->ts = ov->y_pos_s; b->te = ov->y_pos_e;
        b->rev = ov->y_pos_strand;
        b->pchain = 0; b->base = 0; b->pidx = b->pdis = b->aidx = (uint32_t)-1;
        b->el = 0;             // normal block
        BLOCK_SET_REF(b);      // mark as reference-derived
    }

    if (cnt == 0) { free(res); return 0; }
    *out_blocks = res; *out_count = cnt;
    return 0;
}

// Map all unitigs to the reference index. Results of each unitig are
// converted to uc_block_t via overlap_to_uc_block_ref_mode and concatenated
// into a single array returned in out_blocks/out_count.
int unitigs_map_to_reference_batch(ma_ug_t *unitigs,
                                  const ul_idx_t *ref_index,
                                  uc_block_t **out_blocks,
                                  uint64_t *out_count,
                                  const hifiasm_opt_t *opt)
{
    if (!unitigs || !ref_index || !out_blocks || !out_count) return -1;
    *out_blocks = NULL; *out_count = 0;

    uc_block_t **all_results = (uc_block_t**)calloc(unitigs->u.n, sizeof(uc_block_t*));
    uint64_t *all_counts = (uint64_t*)calloc(unitigs->u.n, sizeof(uint64_t));
    if (!all_results || !all_counts) { free(all_results); free(all_counts); return -1; }

    uint64_t total = 0;
    for (uint32_t i = 0; i < unitigs->u.n; i++) {
        const char *seq = ensure_unitig_seq(unitigs, i);
        if (!seq) continue;
        uint64_t len = unitigs->u.a[i].len;

        overlap_region_alloc ov; init_overlap_region_alloc(&ov);
        overlap_region_alloc ov_hp; init_overlap_region_alloc(&ov_hp);
        Candidates_list cl; init_Candidates_list(&cl);
        ha_abufl_t *ab = ha_abufl_init();

        ha_get_ul_candidates_interface(ab, i, (char*)seq, len,
                opt->ul_mz_win ? opt->ul_mz_win : opt->mz_win,
                opt->ul_mer_length ? opt->ul_mer_length : opt->k_mer_length,
                ref_index, &ov, &ov_hp, &cl, 0, opt->max_n_chain, 1,
                NULL, NULL, NULL, NULL, 0, NULL);

        ha_abufl_destroy(ab);

        overlap_to_uc_block_ref_mode(&ov, i, &all_results[i], &all_counts[i]);
        total += all_counts[i];

        destory_overlap_region_alloc(&ov); destory_overlap_region_alloc(&ov_hp);
        destory_Candidates_list(&cl);
    }

    if (total == 0) { free(all_results); free(all_counts); return 0; }
    uc_block_t *final_blocks = (uc_block_t*)malloc(total * sizeof(uc_block_t));
    uint64_t off = 0;
    for (uint32_t i = 0; i < unitigs->u.n; i++) {
        if (all_counts[i] == 0) continue;
        memcpy(final_blocks + off, all_results[i], all_counts[i] * sizeof(uc_block_t));
        off += all_counts[i];
        free(all_results[i]);
    }
    free(all_results); free(all_counts);

    *out_blocks = final_blocks; *out_count = total;
    return 0;
}

// Integrate reference blocks into UL_INF and trigger existing UL update
// routines. Reference blocks are appended to the corresponding unitig
// entry and marked with BLOCK_REF.
int integrate_reference_blocks_to_existing_ul_pipeline(ma_ug_t *unitigs,
                                                      const ul_idx_t *ref_index,
                                                      const hifiasm_opt_t *opt)
{
    if (!unitigs || !ref_index) return -1;

    uc_block_t *blocks = NULL; uint64_t n_block = 0;
    if (unitigs_map_to_reference_batch(unitigs, ref_index, &blocks, &n_block, opt) != 0)
        return -1;

    fprintf(stderr, "[M::%s] Generated %lu reference blocks\n", __func__, (unsigned long)n_block);

    if (n_block == 0) { free(blocks); return 0; }

    if (unitigs->u.n > UL_INF.n) {
        kv_resize(ul_vec_t, UL_INF, unitigs->u.n);
        while (UL_INF.n < unitigs->u.n) {
            memset(&UL_INF.a[UL_INF.n], 0, sizeof(ul_vec_t));
            UL_INF.n++;
        }
    }

    for (uint64_t i = 0; i < n_block; i++) {
        uc_block_t *b = &blocks[i];
        uint32_t uid = b->hid; // here hid stores query unitig id
        if (uid >= UL_INF.n) continue;
        ul_vec_t *p = &UL_INF.a[uid];
        kv_push(uc_block_t, p->bb, *b);
        BLOCK_SET_REF(&p->bb.a[p->bb.n-1]);
    }

    free(blocks);

    filter_ul_ug(unitigs);
    gen_ul_vec_rid_t(&UL_INF, NULL, unitigs);
    update_ug_arch_ul_mul(unitigs);

    fprintf(stderr, "[M::%s] Added %lu reference blocks to UL_INF\n", __func__, (unsigned long)n_block);
    return 0;
}

#endif // ENABLE_REF_GENOME_V4
#endif // __INTER__ (这个必须是文件的最后一行)
