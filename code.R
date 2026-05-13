######################################################################################
options("repos"= c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")
#install.packages("epitools")
rm(list = ls())
library(epitools)
library(ComplexHeatmap)
library(circlize)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
setwd("/hpcdisk1/tianchx_group/panyt/61.LG_PanCancer_BM/02.reAnalysis/12.Cluster_combined")
# 这里读取例文作者提供的注释数据
# 数值为特定组织中特定细胞类型的实际细胞数目
scRNA_harmony <- readRDS("/hpcdisk1/tianchx_group/panyt/61.LG_PanCancer_BM/02.reAnalysis/19.Sub_UMAP_combined/Total_used.rds")
table(scRNA_harmony@meta.data$CellMajor)
table(scRNA_harmony@meta.data$CellMin1)

Plasma<-readRDS("/hpcdisk1/tianchx_group/panyt/61.LG_PanCancer_BM/02.reAnalysis/15.Plasma_sub/Plasma_Sub_test_PC_50.rds")
#Plasma<-subset(Plasma,CellMin %in% setdiff(levels(factor(Plasma@meta.data$CellMin)),"unknown"))
table(Plasma@meta.data$CellMin)
for(i in levels(factor(Plasma@meta.data$CellMin))){
	scRNA_harmony@meta.data[rownames(subset(Plasma@meta.data,CellMin %in% c(i))),"CellMin1"]=i
}
table(scRNA_harmony@meta.data$CellMin1)
scRNA_harmony<-subset(scRNA_harmony,CellMin1 %in% setdiff(levels(factor(scRNA_harmony@meta.data$CellMin1)),"unknown"))

scRNA_harmony@meta.data[1:4,]
table(scRNA_harmony@meta.data$Tissue)
length(levels(factor(scRNA_harmony@meta.data$Sample)))

print(scRNA_harmony)
print(names(scRNA_harmony@meta.data))

# ╔══════════════════════════════════════════════════════════════════╗
# ║   BoM-Atlas  |  Web Data Extraction Script  v2.0               ║
# ║   前提: scRNA_harmony 已在当前 R 环境中加载                      ║
# ║   建议内存: ≥ 32 GB                                              ║
# ╚══════════════════════════════════════════════════════════════════╝
#
#  输出文件结构:
#  web_data/
#  ├─ umap_data.tsv                ← 全量细胞坐标 (整数编码, 节省60%空间)
#  ├─ umap_labels.json             ← 类别标签字典
#  ├─ gene_list.csv                ← 过滤后所有可查询基因
#  ├─ cell_hierarchy.csv           ← 大类→小类 层级
#  ├─ major_name_map.csv           ← 大类名→安全文件名 映射
#  ├─ dp_major_avg.csv             ← Dotplot平均表达 (基因×大类)
#  ├─ dp_major_pct.csv             ← Dotplot表达比例 (基因×大类)
#  ├─ dp_minor_{name}_avg.csv      ← 各大类内 小类Dotplot (按大类分文件)
#  ├─ dp_minor_{name}_pct.csv      ← 各大类内 小类Dotplot
#  ├─ comp_cancer_major/minor.csv  ← 按癌种组成统计
#  └─ comp_tissue_major/minor.csv  ← 按组织组成统计

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(jsonlite)
})

# ════════════════════════════════════════════════════════════════
#  ⚙️  用户配置（只需修改这里）
# ════════════════════════════════════════════════════════════════
obj        <- scRNA_harmony   # 直接使用已加载的 Seurat 对象
OUTPUT_DIR <- "/hpcdisk1/tianchx_group/panyt/61.LG_PanCancer_BM/02.reAnalysis/67.Web_Prepare/web_data"    # 输出目录
COL_MAJOR  <- "CellMajor"    # 细胞大类列名
COL_MINOR  <- "CellMin1"     # 细胞小类列名
COL_TISSUE <- "Tissue"       # 组织来源列名 (PT/BoM/hBM)
COL_CANCER <- "Cancer"       # 癌种列名
COL_SAMPLE <- "Sample"       # 样本列名

# ════════════════════════════════════════════════════════════════
#  初始化 & 验证
# ════════════════════════════════════════════════════════════════
banner <- function(msg) cat(sprintf("\n%s\n  %s\n%s\n",
  paste(rep("═", 52), collapse=""), msg, paste(rep("═", 52), collapse="")))
step   <- function(n, msg) cat(sprintf("\n── STEP %d/5 │ %s\n", n, msg))

banner("BoM-Atlas Web Data Extraction v2.0")
cat(sprintf("  Start: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

dir.create(OUTPUT_DIR, showWarnings=FALSE, recursive=TRUE)
t_global <- proc.time()

# 检查必要列
req_cols <- c(COL_MAJOR, COL_MINOR, COL_TISSUE, COL_CANCER, COL_SAMPLE)
miss     <- setdiff(req_cols, names(obj@meta.data))
if (length(miss) > 0) stop("\n❌ Missing columns: ", paste(miss, collapse=", "))

cat(sprintf("✓ Object: %s cells × %s genes\n",
            formatC(ncol(obj), big.mark=","), formatC(nrow(obj), big.mark=",")))

# 提取并清理元数据
meta <- data.frame(
  cell_id   = colnames(obj),
  CellMajor = as.character(obj@meta.data[[COL_MAJOR]]),
  CellMin1  = as.character(obj@meta.data[[COL_MINOR]]),
  Tissue    = as.character(obj@meta.data[[COL_TISSUE]]),
  Cancer    = as.character(obj@meta.data[[COL_CANCER]]),
  stringsAsFactors = FALSE, row.names = NULL
)
meta[is.na(meta)] <- "Unknown"    # 用 "Unknown" 填充 NA

# 打印概况
cat(sprintf("  CellMajor (%d): %s\n",
            n_distinct(meta$CellMajor),
            paste(sort(unique(meta$CellMajor)), collapse=" | ")))
cat(sprintf("  CellMin1  (%d types)\n", n_distinct(meta$CellMin1)))
cat(sprintf("  Tissue    : %s\n",  paste(sort(unique(meta$Tissue)),  collapse=" | ")))
cat(sprintf("  Cancer    : %d types\n", n_distinct(meta$Cancer)))


# ════════════════════════════════════════════════════════════════
#  STEP 1 │ UMAP 全量细胞（~80万）
#
#  设计说明：将字符类别编码为 0-based 整数，文件体积缩小约 60%。
#  浏览器通过 umap_labels.json 还原标签。
# ════════════════════════════════════════════════════════════════
step(1, sprintf("UMAP Export — ALL %s cells", formatC(ncol(obj), big.mark=",")))
t0 <- proc.time()

# 自动检测 UMAP reduction
umap_key <- grep("umap", names(obj@reductions), ignore.case=TRUE, value=TRUE)[1]
if (is.na(umap_key)) {
  stop("❌ No UMAP found. Available reductions: ",
       paste(names(obj@reductions), collapse=", "))
}
cat(sprintf("  Using reduction: [[ %s ]]\n", umap_key))

coords <- obj@reductions[[umap_key]]@cell.embeddings[, 1:2]

# 构建整数编码映射
major_lv  <- sort(unique(meta$CellMajor))
minor_lv  <- sort(unique(meta$CellMin1))
tissue_lv <- sort(unique(meta$Tissue))
cancer_lv <- sort(unique(meta$Cancer))

labels_dict <- list(
  major  = major_lv,
  minor  = minor_lv,
  tissue = tissue_lv,
  cancer = cancer_lv
)

umap_out <- data.frame(
  x  = round(coords[, 1], 2),                        # 保留2位小数
  y  = round(coords[, 2], 2),
  ma = match(meta$CellMajor, major_lv)  - 1L,        # 0-based → 节省字节
  mi = match(meta$CellMin1,  minor_lv)  - 1L,
  ti = match(meta$Tissue,    tissue_lv) - 1L,
  ca = match(meta$Cancer,    cancer_lv) - 1L
)

# TSV 格式（无引号）：800K行 × 6列 ≈ 20-25 MB
write.table(umap_out,
            file.path(OUTPUT_DIR, "umap_data.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)

write_json(labels_dict, file.path(OUTPUT_DIR, "umap_labels.json"), auto_unbox=TRUE)

sz_umap <- file.size(file.path(OUTPUT_DIR, "umap_data.tsv")) / 1024^2
cat(sprintf("✓  umap_data.tsv    : %.1f MB\n", sz_umap))
cat(sprintf("✓  umap_labels.json : %d major | %d minor | %d tissue | %d cancer\n",
            length(major_lv), length(minor_lv), length(tissue_lv), length(cancer_lv)))
cat(sprintf("   Elapsed: %.0f s\n", (proc.time()-t0)[3]))


# ════════════════════════════════════════════════════════════════
#  STEP 2 │ 基因过滤
#
#  去除: 核糖体蛋白 / 线粒体 / LINC / 含点号 / LOC / ENSG ID
#  保留: 所有其余基因（通常 15,000–20,000 个）
# ════════════════════════════════════════════════════════════════
step(2, "Gene Filtering")
t0 <- proc.time()

gn    <- rownames(obj)
valid <- rep(TRUE, length(gn))

valid[grepl("^RP[SL]\\d+", gn)]                   <- FALSE   # 核糖体蛋白
valid[grepl("^MT-",         gn, ignore.case=TRUE)] <- FALSE   # 线粒体基因
valid[grepl("^LINC",        gn, ignore.case=TRUE)] <- FALSE   # 长链非编码RNA
valid[grepl("\\.",          gn)]                   <- FALSE   # 含点号（无效符号）
valid[grepl("^LOC",         gn, ignore.case=TRUE)] <- FALSE   # 未注释基因座
valid[grepl("^ENSG",        gn, ignore.case=TRUE)] <- FALSE   # Ensembl ID

fg <- gn[valid]   # filtered genes

# 打印过滤明细
cat(sprintf("  Total genes   : %6d\n", length(gn)))
cat(sprintf("  ─ RPS/RPL     : %6d\n", sum(grepl("^RP[SL]\\d+", gn))))
cat(sprintf("  ─ MT-         : %6d\n", sum(grepl("^MT-",  gn, ignore.case=TRUE))))
cat(sprintf("  ─ LINC        : %6d\n", sum(grepl("^LINC", gn, ignore.case=TRUE))))
cat(sprintf("  ─ dot (.)     : %6d\n", sum(grepl("\\.", gn))))
cat(sprintf("  ─ LOC         : %6d\n", sum(grepl("^LOC", gn, ignore.case=TRUE))))
cat(sprintf("  ─ ENSG        : %6d\n", sum(grepl("^ENSG", gn, ignore.case=TRUE))))
cat(sprintf("  ═ Kept        : %6d genes ✓\n", length(fg)))

write.csv(data.frame(gene = fg),
          file.path(OUTPUT_DIR, "gene_list.csv"), row.names=FALSE)
cat(sprintf("✓  gene_list.csv (%d genes)\n", length(fg)))
cat(sprintf("   Elapsed: %.0f s\n", (proc.time()-t0)[3]))


# ════════════════════════════════════════════════════════════════
#  STEP 3 │ Dotplot 统计量
#
#  为每个基因 × 每个细胞类型计算：
#    avg_expr = 平均归一化表达（全细胞，含0）
#    pct_expr = 表达该基因的细胞百分比（非零比例）
#
#  输出格式（宽矩阵 CSV）：
#    行 = 基因,  列 = 细胞类型
#    前端加载后建立 gene→行 索引，支持任意基因即时查询
#
#  ⚡ 速度关键：rowNNZ() 直接读取稀疏矩阵内部结构，
#              避免创建中间逻辑矩阵，速度快 3-5×
# ════════════════════════════════════════════════════════════════
step(3, "Dotplot Statistics (avg_expr + pct_expr)")
cat(sprintf("  %s genes × %s cells\n",
            formatC(length(fg), big.mark=","),
            formatC(ncol(obj), big.mark=",")))
cat("  ⏳ Estimated 5–20 min depending on RAM & CPU...\n\n")
t0 <- proc.time()

# ── 加载稀疏表达矩阵 ─────────────────────────────────────────────
cat("  Loading RNA data slot...\n")
expr_mat <- tryCatch(
  GetAssayData(obj, assay="RNA", layer="data"),  # Seurat v5
  error = function(e)
    GetAssayData(obj, assay="RNA", slot="data")  # Seurat v4
)
expr_mat <- expr_mat[fg, ]     # 只保留过滤后的基因
gc()

nnz       <- nnzero(expr_mat)
density   <- nnz / (nrow(expr_mat) * ncol(expr_mat)) * 100
cat(sprintf("  Sparse matrix: %s genes × %s cells\n",
            formatC(nrow(expr_mat), big.mark=","),
            formatC(ncol(expr_mat), big.mark=",")))
cat(sprintf("  Density: %.1f%% | Non-zeros: %s\n\n",
            density, formatC(nnz, format="d", big.mark=",")))

# ── 工具函数 ─────────────────────────────────────────────────────

# 最快的非零行计数：直接读 CsparseMatrix 的 @i slot（行索引数组）
row_nnz <- function(mat) {
  if (inherits(mat, "CsparseMatrix")) {
    tabulate(mat@i + 1L, nbins = nrow(mat))   # 直接统计，无需创建逻辑矩阵
  } else {
    as.integer(Matrix::rowSums(mat != 0))
  }
}

# 计算一组 group 的 avg + pct 统计
dot_stats <- function(mat, groups, step_label = "") {
  glevels  <- sort(unique(as.character(groups)))
  ng       <- length(glevels)
  avg_m    <- matrix(0.0, nrow(mat), ng, dimnames = list(rownames(mat), glevels))
  pct_m    <- matrix(0.0, nrow(mat), ng, dimnames = list(rownames(mat), glevels))

  t_in <- proc.time()
  for (k in seq_along(glevels)) {
    grp    <- glevels[k]
    cidx   <- which(groups == grp)
    nc     <- length(cidx)
    sub    <- mat[, cidx, drop = FALSE]

    avg_m[, k] <- Matrix::rowMeans(sub)           # 含0的平均值
    pct_m[, k] <- row_nnz(sub) / nc * 100         # % 非零细胞

    ela    <- (proc.time() - t_in)[3]
    eta    <- if (k > 1) sprintf("  ETA ≈ %.0f min", ela / k * (ng - k) / 60) else ""
    cat(sprintf("    [%2d/%2d] %-28s  n = %s%s\n",
                k, ng, grp, formatC(nc, big.mark=","), eta))

    rm(sub)
    if (k %% 4 == 0) gc()
  }

  list(avg = round(avg_m, 4), pct = round(pct_m, 2))
}

# 写出宽矩阵 CSV（行=基因，列=细胞类型）
write_dp <- function(st, prefix, genes) {
  a_df <- cbind(data.frame(gene = genes, stringsAsFactors = FALSE),
                as.data.frame(st$avg))
  p_df <- cbind(data.frame(gene = genes, stringsAsFactors = FALSE),
                as.data.frame(st$pct))
  f_avg <- file.path(OUTPUT_DIR, paste0(prefix, "_avg.csv"))
  f_pct <- file.path(OUTPUT_DIR, paste0(prefix, "_pct.csv"))
  write.csv(a_df, f_avg, row.names = FALSE)
  write.csv(p_df, f_pct, row.names = FALSE)
  cat(sprintf("  ✓ %s_avg.csv  (%.1f MB)\n",
              prefix, file.size(f_avg) / 1024^2))
  cat(sprintf("  ✓ %s_pct.csv  (%.1f MB)\n",
              prefix, file.size(f_pct) / 1024^2))
}

# ── 3a: 按 CellMajor ─────────────────────────────────────────────
cat("[3a] CellMajor level dotplot stats\n")
st_major <- dot_stats(expr_mat, meta$CellMajor, "CellMajor")
write_dp(st_major, "dp_major", fg)
rm(st_major); gc()


# ── 3b: 按 CellMin1（每个大类单独一对文件，前端按需加载）──────────
cat("\n[3b] CellMin1 level dotplot stats (one file pair per major type)\n")

major_types  <- sort(unique(meta$CellMajor))

# 生成唯一的安全文件名（只保留字母数字，用 make.unique 保证不重复）
safe_fn    <- function(x) gsub("[^A-Za-z0-9]", "_", trimws(x))
safe_names <- make.unique(safe_fn(major_types), sep = "_v")

major_map_df <- data.frame(
  CellMajor = major_types,
  safe_name = safe_names,
  stringsAsFactors = FALSE
)
write.csv(major_map_df, file.path(OUTPUT_DIR, "major_name_map.csv"), row.names = FALSE)

for (j in seq_along(major_types)) {
  maj   <- major_types[j]
  sname <- safe_names[j]
  idx   <- which(meta$CellMajor == maj)
  n_sub <- n_distinct(meta$CellMin1[idx])

  cat(sprintf("\n  [%d/%d] %s  (%s cells, %d subtypes)\n",
              j, length(major_types), maj,
              formatC(length(idx), big.mark=","), n_sub))

  sub_mat <- expr_mat[, idx, drop = FALSE]
  sub_grp <- meta$CellMin1[idx]

  st <- dot_stats(sub_mat, sub_grp, sname)
  write_dp(st, paste0("dp_minor_", sname), fg)

  rm(sub_mat, sub_grp, st, idx)
  gc()
}

cat(sprintf("\n✓ Step 3 total: %.1f min\n", (proc.time() - t0)[3] / 60))


# ════════════════════════════════════════════════════════════════
#  STEP 4 │ 组成分析（癌种 & 组织）
# ════════════════════════════════════════════════════════════════
step(4, "Composition Stats")
t0 <- proc.time()

comp_list <- list(
  cancer_major = meta %>%
    count(Cancer, CellMajor, name = "n") %>%
    group_by(Cancer) %>%
    mutate(pct = round(n / sum(n) * 100, 3)) %>% ungroup(),

  cancer_minor = meta %>%
    count(Cancer, CellMajor, CellMin1, name = "n") %>%
    group_by(Cancer, CellMajor) %>%
    mutate(pct = round(n / sum(n) * 100, 3)) %>% ungroup(),

  tissue_major = meta %>%
    count(Tissue, CellMajor, name = "n") %>%
    group_by(Tissue) %>%
    mutate(pct = round(n / sum(n) * 100, 3)) %>% ungroup(),

  tissue_minor = meta %>%
    count(Tissue, CellMajor, CellMin1, name = "n") %>%
    group_by(Tissue, CellMajor) %>%
    mutate(pct = round(n / sum(n) * 100, 3)) %>% ungroup()
)

for (nm in names(comp_list)) {
  fpath <- file.path(OUTPUT_DIR, paste0("comp_", nm, ".csv"))
  write.csv(comp_list[[nm]], fpath, row.names = FALSE)
  cat(sprintf("  ✓ comp_%s.csv  (%d rows)\n", nm, nrow(comp_list[[nm]])))
}
cat(sprintf("   Elapsed: %.0f s\n", (proc.time() - t0)[3]))


# ════════════════════════════════════════════════════════════════
#  STEP 5 │ 细胞层级关系
# ════════════════════════════════════════════════════════════════
step(5, "Cell Type Hierarchy")

hierarchy <- meta %>%
  distinct(CellMajor, CellMin1) %>%
  arrange(CellMajor, CellMin1)

write.csv(hierarchy, file.path(OUTPUT_DIR, "cell_hierarchy.csv"), row.names = FALSE)
cat(sprintf("  ✓ cell_hierarchy.csv  (%d major → %d minor)\n",
            n_distinct(hierarchy$CellMajor), n_distinct(hierarchy$CellMin1)))


# ════════════════════════════════════════════════════════════════
#  汇总报告
# ════════════════════════════════════════════════════════════════
banner("🎉  Extraction Complete!")
total_min <- (proc.time() - t_global)[3] / 60
cat(sprintf("\n  Total runtime : %.1f min\n", total_min))
cat(sprintf("  Finish time   : %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# 文件清单
all_files <- list.files(OUTPUT_DIR, full.names = TRUE, recursive = FALSE)
cat(sprintf("  %-48s  %8s\n", "File", "Size(MB)"))
cat(sprintf("  %s\n", paste(rep("─", 58), collapse="")))

tot_mb <- 0
for (f in sort(all_files)) {
  mb     <- file.size(f) / 1024^2
  tot_mb <- tot_mb + mb
  cat(sprintf("  %-48s  %7.1f\n", basename(f), mb))
}
cat(sprintf("  %s\n", paste(rep("─", 58), collapse="")))
cat(sprintf("  %-48s  %7.1f\n", "TOTAL", tot_mb))

cat("\n💡 Next steps:\n")
cat("   1. 将 index.html 与 web_data/ 放在同一目录\n")
cat("   2. 本地预览: python -m http.server 8000\n")
cat("      浏览器打开: http://localhost:8000\n")
cat("   3. 上线: 将二者上传至 GitHub → GitHub Pages\n\n")

#####################################################################
python -m http.server 8000
#成功的样子应该是：
#Serving HTTP on :: port 8000 (http://[::]:8000/) ...

#用网页打开下面这个
http://localhost:8000/