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
str(dat)
#install.packages('SNPassoc')
library(SNPassoc)
#https://cran.r-project.org/web/packages/SNPassoc/vignettes/SNPassoc.html
help(package='SNPassoc')
colnames(dat)=gsub(x = colnames(dat),pattern ='.x',replacement = '' )
dat <- dat %>%
  mutate(across(9:12, ~ case_when(
    . == "A" ~ 0,
    . == "B" ~ 1
  )))

dat$s=ifelse(dat$s=='M'|dat$s=='男',1,0)
colnames(dat)[7]='Gender'
dat$Age=ifelse(dat$Age>median(dat$Age),1,0)
str(dat)
#------
idx <- grep("^rs", colnames(dat))
asthma.s <- setupSNP(data=dat, colSNPs=idx, sep="")
head(asthma.s$rs7754561)
plotMissing(asthma.s, print.labels.SNPs = FALSE)
##---
genotype_Result1=as.data.frame(summary(asthma.s, print=FALSE))

hwe <- tableHWE(asthma.s)
head(hwe)
association(group ~ rs3759890, data = asthma.s)
table(asthma.s$rs7754561)
str(asthma.s)
genotype_Result2=WGassociation(group, data=asthma.s)
write_csv(genotype_Result1,file = "./genotype_Result1.csv")
write_csv(genotype_Result2,file = "./genotype_Result2.csv")
association(group ~ dominant(rs10061133)*factor(Glucose), 
            data=asthma.s)
###---

##----gene alele
str(dat)
df_expanded <- dat%>%
  # 对每个以 rs 开头的列应用 separate_rows
  mutate(across(starts_with("rs"), ~ str_split(., ""))) %>%
  unnest(cols = starts_with("rs"))
# 查看处理后的数据
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
df <- df_expanded  # 替换为您的数据框名称
# 获取所有以 rs 开头的列名
rs_columns <- names(df)[grepl("^rs", names(df))]
# 获取环境变量列名
env_columns <- c("Gender", "Age", "Duration", "Glucose", "HbA1c", "Treatment")
# 初始化一个空列表来存储结果
results_list <- list()
# 函数：进行卡方检验并保存结果
perform_chi_square_test <- function(rs_col, env_col) {
  # 创建列联表
  contingency_table <- table(df[[rs_col]], df[[env_col]])
  
  # 进行卡方检验
  chi_square_result <- chisq.test(contingency_table)
  
  # 提取统计量和 p 值
  result <- data.frame(
    rs_column = rs_col,
    env_variable = env_col,
    statistic = chi_square_result$statistic,
    p_value = chi_square_result$p.value,
    stringsAsFactors = FALSE
  )
  
  # 将列联表转换为数据框并添加到结果中
  contingency_df <- as.data.frame(contingency_table)
  colnames(contingency_df) <- c("rs_value", "env_value", "count")
  
  # 合并结果
  final_result <- result %>%
    left_join(contingency_df, by = character())  # 不用按任何列连接，直接合并
  
  return(final_result)
}

# 对每一个 rs 列与每一个环境变量进行卡方检验
for (rs_col in rs_columns) {
  for (env_col in env_columns) {
    results_list[[length(results_list) + 1]] <- perform_chi_square_test(rs_col, env_col)
  }
}

# 将所有结果合并成一个数据框
results_df <- bind_rows(results_list)

# 将特定列转换为因子，设置 levels 为 0 和 1
results_df <- results_df %>%
  mutate(across(c(rs_value, env_value), ~ factor(., levels = c(0, 1))))

# 查看结果
print(results_df)

# 保存结果到 CSV 文件
write.csv(results_df, "chi_square_results_rs_env.csv", row.names = FALSE)

###----互作-yu group

df <- df_expanded  # 替换为您的数据框名称

# 获取所有以 rs 开头的列名
rs_columns <- names(df)[grepl("^rs", names(df))]

# 获取环境变量列名
env_columns <- c("Gender", "Age", "Duration", "Glucose", "HbA1c", "Treatment")
# 初始化一个空列表来存储结果
results_list <- list()
# 函数：进行卡方检验并保存结果
perform_chi_square_test <- function(rs_col, env_col) {
  # 创建新的列，表示共存情况
  df_expanded <- df %>%
    mutate(coexist = rowSums(cbind(df[[rs_col]], df[[env_col]])))
           # 标记共存为 2 的情况
           df_expanded <- df_expanded %>%
             mutate(coexist = ifelse(coexist == 2, 1, 0))
           # 创建列联表与 group
           contingency_table <- table(df_expanded$coexist, df_expanded$group)
           
           # 进行卡方检验
           chi_square_result <- chisq.test(contingency_table)
           
           # 提取统计量和 p 值
           result <- data.frame(
             rs_column = rs_col,
             env_variable = env_col,
             statistic = chi_square_result$statistic,
             p_value = chi_square_result$p.value,
             stringsAsFactors = FALSE)
           
           # 将列联表转换为数据框并添加到结果中
           contingency_df <- as.data.frame(contingency_table)
           colnames(contingency_df) <- c("coexist", "group", "count")
           
           # 合并结果
           final_result <- result %>%
             left_join(contingency_df, by = character())  # 不用按任何列连接，直接合并
           return(final_result)
}

# 对每一个 rs 列与每一个环境变量进行共存分析与卡方检验
for (rs_col in rs_columns) {
  for (env_col in env_columns) {
    results_list[[length(results_list) + 1]] <- perform_chi_square_test(rs_col, env_col)
  }
}

# 将所有结果合并成一个数据框
results_df <- bind_rows(results_list)

# 将特定列转换为因子，设置 levels 为 0 和 1
results_df <- results_df %>%
  mutate(across(c(coexist, group), ~ factor(., levels = c(0, 1))))

# 查看结果
print(results_df)

# 保存结果到 CSV 文件
write.csv(results_df, "chi_square_results_coexist_group.csv", row.names = FALSE)

