---
title: "ggtree: KEGG Pathway Heatmap！"
categories: 
 - Comparative genomics
tags: 
 - KEGG
 - R
---

# Introduction

<div style="text-align: center; margin-bottom: 20px;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-28-tree-heatmap-1.png"/>
</div>

KEGG pathway and module correspondence information. [Click][click-1]

## Data preparation

- tools_RStudio_result_all
<div style="text-align: center; margin-bottom: 20px;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-28-barplot-2.png"/>
</div>

- Maps_Moules_2023_11_16
<div style="text-align: center; margin-bottom: 20px;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-28-tree-heatmap-2.png"/>
</div>

- tools_RStudio_映射模块上的K号
<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-28-tree-heatmap-3.png"/>
</div>

## R script for tree + heatmap

{% highlight r %}
#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Time      : 20231116

library(readr)
library(tidyfst)

tools_RStudio_result_all <- read_delim("tools_RStudio_result_all.txt", 
                                       delim = "\t", escape_double = FALSE, 
                                       trim_ws = TRUE)

Maps_Moules_2023_11_16 <- read_delim("F:/database/KEGG/KEGG_General_Documentation_module/Maps_Moules_2023_11_16.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     col_names = FALSE, trim_ws = TRUE)

tools_RStudio_映射模块上的K号 <- read_delim("tools_RStudio_映射模块上的K号.txt", 
                                     delim = "\t", escape_double = FALSE, 
                                     col_names = FALSE, trim_ws = TRUE)


map01230 <- Maps_Moules_2023_11_16 %>% filter_dt(X2 == "map01230")
map01230 <- map01230 %>% rename_dt(m_id = X1)
map01230_1 <- map01230 %>% inner_join_dt(tools_RStudio_result_all,by = "m_id")
write.table(map01230_1,file = "map01230.txt",sep = "\t",quote = F,row.names = F)
{% endhighlight %}

- map01230.txt
<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-28-tree-heatmap-4.png"/>
</div>

{% highlight r %}
#!/usr/bin/env Rscript
# -*- coding: utf-8 -*-
# @Author    : mengqingyao
# @Time      : 20231118

library(ggtree)
library(ggtreeExtra)
library(tidyfst)

map01230 <- read_delim("map01230.txt", delim = "\t", 
                       escape_double = FALSE, trim_ws = TRUE)

map01230 <- map01230 %>% select(-X2)

head_data <- map01230 %>% longer_dt(name = "species",value = "complet")

a <- map01230 %>% column_to_rownames(var="m_name")

tree <- map01230 %>% column_to_rownames(var="m_name") %>%
  #  该函数使用指定的距离度量计算数据矩阵行间的距离，计算并返回距离矩阵。
  dist() %>% 
  #  Tree Estimation Based on an Improved Version of the NJ Algorithm
  ape::bionj()

ggtree(tree,branch.length = "none", layout = "circular",
       linetype = 1,size = 0.5, ladderize = T)+
  layout_fan(angle =180)+ # 调整开口
  geom_tiplab(offset=11,show.legend=FALSE,size=2,
              color = "black",starstroke = 0)+
  geom_fruit(data=head_data,
                       geom=geom_tile,
                       mapping=aes(y=m_name,x=species,fill=complet),
                       pwidth=0.8,
                       offset=0.04,
                       color = "white",
                       axis.params=list(axis="x",text.angle=-90,text.size=1,hjust=0))+
  scale_fill_gradientn(colours = RColorBrewer::brewer.pal(100,"Blues")) +
  theme(plot.margin=margin(8,-6,0,0,"cm"),# 上下,左右,远近
        legend.position = c(1.3,0.8),
        legend.background = element_rect(fill = NA), 
        )
{% endhighlight %}

**[Rsudio][rstudio-doc] is recommended for running this script.**

<div class="notice">
  <h4>I wrote this script a long time ago, and there is a lot of redundant code, but this is the original processing logic.<br><br>It's better to understand.</h4>
</div>

## Result

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-28-tree-heatmap-5.png"/>
</div>

## Beautification

1. I would suggest saving the image as a pdf.
2. Evolutionary landscaping using Adobe Illustrator.
3. Save images as needed.

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-28-tree-heatmap-6.png"/>
</div>

## Quote 

> - [RStudio][rstudio-doc]
> - [tidyfst][tidyfst-doc]
> - [ggtree][ggtree-doc]
> - [ggtreeExtra][ggtreeExtra-doc]

> Email me with more questions!
> 584338215@qq.com

[click-1]: https://rest.kegg.jp/link/module%20/pathway
[rstudio-doc]: https://posit.co/
[tidyfst-doc]: https://hope-data-science.github.io/tidyfst/
[ggtree-doc]: https://guangchuangyu.github.io/software/ggtree/
[ggtreeExtra-doc]: https://bioconductor.org/packages/devel/bioc/vignettes/ggtreeExtra/inst/doc/ggtreeExtra.html

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