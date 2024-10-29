---
title: "Reviewing the past, one can understand the new！"
categories: 
 - Code Optimization
tags: 
 - R
---

# Introduction

Here are the general requirements for writing code:

1. **Code Readability**:
   - Use meaningful variable and function names so that others can understand what the code does.
   - Utilize comments appropriately to explain complex algorithms or logic.
   - Follow style guidelines such as indentation and spacing to keep the code neat and consistent.

2. **Modular Design**:
   - Break code down into small, reusable functions or modules, with each module responsible for a specific task.
   - Avoid duplicate code by using functions or classes to reuse logic.

3. **Error Handling**:
   - Implement error handling mechanisms to ensure that the code can gracefully handle exceptions rather than crashing.
   - Use exception handling to catch and manage potential errors.

4. **Performance Optimization**:
   - For applications with high performance requirements, pay attention to the time complexity and space complexity of algorithms.
   - Collect and analyze performance data, and optimize as necessary.

5. **Testing**:
   - Write unit tests to ensure the correctness of each functional module.
   - Use testing frameworks to automate the testing process and ensure that functionality is not broken after code changes.

6. **Documentation**:
   - Write project documentation to explain how to use, set up, and understand the functionalities of the code.
   - Provide API documentation to help developers understand how to utilize your codebase.

7. **Version Control**:
   - Use version control tools (such as Git) to manage the history of code changes.
   - Make clear commits with meaningful messages.

8. **Adhering to Best Practices**:
   - Always keep an eye on the best practices for the language and framework you are using, updating code to conform to the latest standards.
   - Engage in code reviews to learn from colleagues and improve code quality.

## Example One: Venn Diagram

> Detailed content [Click][click-veen]

- modified code

{% highlight r %}
#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Time      : 20241029
# @File      : veen.r

library(readr)
library(tidyverse)
library(ggVennDiagram)
library(R6)

# 创建一个R6类
VennDiagramAnalysis <- R6Class("VennDiagramAnalysis",
  public = list(
    orthogroups_file = NULL,
    unassigned_genes_file = NULL,
    output_prefix = "venn_diagram",  # 输出前缀
    plot_width = 20,                  # 图片宽度
    plot_height = 18,                 # 图片高度
    Orthogroup_all = NULL,
    orth_gen = NULL,
    result = NULL,

    initialize = function(orthogroups_file, unassigned_genes_file, output_prefix = "venn_diagram", plot_width = 20, plot_height = 18) {
      self$orthogroups_file <- orthogroups_file
      self$unassigned_genes_file <- unassigned_genes_file
      self$output_prefix <- output_prefix
      self$plot_width <- plot_width
      self$plot_height <- plot_height
      self$load_data()
      self$process_data()
    },

    load_data = function() {
      Orthogroups <- read_tsv(self$orthogroups_file)
      Orthogroups_UnassignedGenes <- read_tsv(self$unassigned_genes_file)
      self$Orthogroup_all <- bind_rows(Orthogroups, Orthogroups_UnassignedGenes) %>%
        select(-"Hartmannula_sinica_185.fna_remove")
    },

    first_obse = function(data, n = 1, sep = ",") {
      list <- as.list(data)
      for (v in 2:length(colnames(data))) {
        for (o in 1:length(rownames(data))) {
          list[[v]][o] <- list[[v]][[o]] %>% strsplit(sep)
          list[[v]][o] <- list[[v]][[o]][n]
        }
      }
      for (v in 2:length(colnames(data))) {
        list[[v]] <- as.character(list[[v]])
      }
      as.data.frame(list)
    },

    process_data = function() {
      Orthogroup_all <- self$first_obse(self$Orthogroup_all)
      col_names <- colnames(Orthogroup_all)[-1]

      unique_symbols <- unique(gsub("[A-Za-z0-9]", "", col_names))
      split_content_list <- lapply(col_names, function(name) {
        strsplit(name, paste0("[", paste(unique_symbols, collapse = ""), "]"))[[1]]
      })

      col_renames <- sapply(split_content_list, function(split_content) {
        paste(head(split_content, 2), collapse = "_")
      })

      self$orth_gen <- lapply(seq_along(col_renames), function(i) {
        Orthogroup_all %>%
          filter(!!rlang::sym(col_names[i]) != "NA") %>%
          select(Orthogroup, !!rlang::sym(col_names[i])) %>%
          rename(!!col_renames[i] := "Orthogroup")
      })

      result <- Orthogroup_all %>% select(-Orthogroup)

      for (i in 1:length(col_names)) {
        result <- result %>%
          left_join(self$orth_gen[[i]], by = setNames(col_names[i], colnames(self$orth_gen[[i]])[2]))
      }

      self$result <- result %>% select((length(col_names) + 1):length(result))
    },

    plot_venn_diagram = function() {
      p1 <- ggVennDiagram(self$result, label = "count", edge_size = 1, label_alpha = 0) +
        scale_color_brewer(palette = "Set1") +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(palette = "Blues", direction = 1)

      p2 <- ggVennDiagram(self$result, label = "count", edge_size = 1, label_alpha = 0, 
                          set_color = c("blue", "red", "green", "purple", "orange", "yellow", "cyan")) +
        scale_color_brewer(palette = "Set1") +
        scale_x_continuous(expand = expansion(mult = .2)) +
        scale_fill_distiller(palette = "Blues", direction = 1)

      # 使用输出前缀和尺寸进行保存
      ggsave(filename = paste0(self$output_prefix, "_p1.pdf"), plot = p1, width = self$plot_width, height = self$plot_height, units = "cm")
      ggsave(filename = paste0(self$output_prefix, "_p2.pdf"), plot = p2, width = self$plot_width, height = self$plot_height, units = "cm")
    }
  )
)
{% endhighlight %}

**Save the above code to a file, e.g. veen.r**

### Now, let's test the code:

{% highlight r %}
source(veen.r)

venn_analysis <- VennDiagramAnalysis$new(
  orthogroups_file = "./OrthoFinder/OrthoFinder/Results_May31/Orthogroups/Orthogroups.tsv",
  unassigned_genes_file = "./OrthoFinder/OrthoFinder/Results_May31/Orthogroups/Orthogroups_UnassignedGenes.tsv",
  output_prefix = "my_venn_diagram",  # 自定义输出前缀
  plot_width = 25,                     # 自定义宽度
  plot_height = 20                      # 自定义高度
)

venn_analysis$plot_venn_diagram()
{% endhighlight %}

**The output will be two pdf files: my_venn_diagram_p1.pdf and my_venn_diagram_p2.pdf.**

- The first change is to remove the redundant code.
- The second change is to use R6 class.
- The third change is that it is simpler to use.

> Email me with more questions!
> 584338215@qq.com

[click-veen]: https://mengqy2022.github.io/comparative%20genomics/Veen-plot/

<script src="https://giscus.app/client.js"
        data-repo="mengqy2022/mengqy2022.github.io"
        data-repo-id="R_kgDONFQ-nw"
        data-category="Announcements"
        data-category-id="DIC_kwDONFQ-n84CjtiY"
        data-mapping="pathname"
        data-strict="0"
        data-reactions-enabled="1"
        data-emit-metadata="0"
        data-input-position="bottom"
        data-theme="dark_high_contrast"
        data-lang="zh-CN"
        crossorigin="anonymous"
        async>
</script>