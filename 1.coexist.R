library(tidyverse)
library(readxl)

dat1=read_excel("./DNR-clinical.xlsx")
dat2=read_excel("./DR-clinical.xlsx")

dat1$group='0'
dat2$group='1'

dat=rbind(dat1,dat2)
fill_na_with_mode <- function(df) {
  # 定义一个辅助函数来计算众数
  get_mode <- function(x) {
    uniq_x <- unique(x)
    uniq_x[which.max(tabulate(match(x, uniq_x)))]
  }
  
  # 对以 "rs" 开头的列进行处理
  df <- df %>%
    mutate(across(starts_with("rs"), ~ ifelse(is.na(.), get_mode(na.omit(.)), .)))
  
  return(df)
}
dat=fill_na_with_mode(dat)
##----gene alele
str(dat)
df_expanded <- dat%>%
# 对每个以 rs 开头的列应用 separate_rows
mutate(across(starts_with("rs"), ~ str_split(., ""))) %>%
  unnest(cols = starts_with("rs"))
# 查看处理后的数据
table(df_expanded$rs7754561)

library(dplyr)
library(tidyr)
# 假设你的数据框名为 df_expanded
# 第一部分：处理 rs 列
for (rs_col in names(df_expanded)[grepl("^rs", names(df_expanded))]) {
  # 计算频次
  allele_freqs <- df_expanded %>%
    count(.data[[rs_col]]) %>%
    arrange(desc(n))
  
  # 找到高频和低频等位基因
  if (nrow(allele_freqs) > 1) {
    high_freq <- allele_freqs[[1, 1]]
    low_freq <- allele_freqs[[nrow(allele_freqs), 1]]
    
    # 更新等位基因
    df_expanded[[rs_col]] <- ifelse(df_expanded[[rs_col]] == high_freq, 0,
                                    ifelse(df_expanded[[rs_col]] == low_freq, 1, NA))
  }
}


# 第二部分：处理其他列
df_expanded <- df_expanded %>%
  mutate(across(9:12, ~ case_when(
    . == "A" ~ 0,
    . == "B" ~ 1
  )))

df_expanded$s=ifelse(df_expanded$s=='M'|df_expanded$s=='男',1,0)
# 查看处理后的数据
print(df_expanded)
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)

# 假设你的数据框名为 df_expanded

# 获取所有 rs 列名
rs_columns <- names(df_expanded)[grepl("^rs", names(df_expanded))]

# 初始化一个空列表来存储结果
results_list <- list()

# 函数：进行卡方检验并保存结果
perform_chi_square_test <- function(columns) {
  # 创建一个新列，表示这些 rs 列是否全部为 1
  df_expanded <- df_expanded %>%
    mutate(all_one = rowSums(select(., all_of(columns)) == 1) == length(columns))
  
  # 创建列联表
  contingency_table <- table(df_expanded$all_one, df_expanded$group)
  
  # 进行卡方检验
  chi_square_result <- chisq.test(contingency_table)
  
  # 提取统计量和 p 值
  result <- data.frame(
    combination = paste(columns, collapse = "_"),
    statistic = chi_square_result$statistic,
    p_value = chi_square_result$p.value,
    stringsAsFactors = FALSE
  )
  
  # 将列联表转换为数据框并添加到结果中
  contingency_df <- as.data.frame(contingency_table)
  colnames(contingency_df) <- c("all_one", "group", "count")
  
  # 合并结果
  final_result <- result %>%
    left_join(contingency_df, by = character())  # 不用按任何列连接，直接合并
  
  return(final_result)
}

# 对所有可能的组合进行卡方检验
for (n in 2:length(rs_columns)) {
  combn(rs_columns, n, FUN = function(cols) {
    results_list[[length(results_list) + 1]] <<- perform_chi_square_test(cols)
  }, simplify = FALSE)
}

# 将所有结果合并成一个数据框
results_df <- bind_rows(results_list)

# 查看结果
print(results_df)

# 保存结果到 CSV 文件
write.csv(results_df, "chi_square_results_summary_with_contingency.csv", row.names = FALSE)

###----互作


