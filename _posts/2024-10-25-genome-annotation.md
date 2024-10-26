---
title: "Bacterial Genome Annotation!"
categories: 
 - Genomics
tags: 
 - Bacteria
 - bash
 - python
---

# Introduction

Now, we will genome annotate the bacterial genome, obtain protein coding sequences.

These data can then be used for functional gene annotation, phylogenetic analysis, etc.

## Data preparation

### Obtain the Genome Sequence

- First, we can download the bacterial genome sequence from NCBI or other sources.
- Secondly, through high-throughput sequencing data (HTS) technology, we can obtain the sequence of the bacterial genome.

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-25-genome-annotation-1.png"/>
</div>

### Preparing genome name files

{% highlight bash %}
ls Genome/ | cut -d '.' -f 1 > head.txt
{% endhighlight %}

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-25-genome-annotation-2.png"/>
</div>

## Genome Annotation workflow

{% highlight bash %}
#!/bin/bash/
# **************************************************
# 共生菌基因预测以及去除假基因
# Date   : 2023-09-25
# Author : 孟庆瑶
# **************************************************

#  发生错误停止脚本
set -e

#  获取脚本名称并输出
name=$(basename $0) 
echo 
echo "   脚本名称: $name"

#  设置参数以及处理
OUTPREFIX="Bacterial_results"

while getopts hi:g:o: opt
do
    case "$opt" in 
    h)
        echo
        echo "脚本说明: 基因预测、假基因预测和获得去除假基因的CDs序列，(*￣︶￣)，东西不多，一字一字看！"
        echo
        echo "使用说明: bash $name -i genome_file_name -g genome_fasta -o out_prefix"
        echo "------------------------------------------------------------------------------------"
        echo -e "\t-i: 输入文本文件、包含基因文件名称、不带后缀、单列文件，一列为一个基因组；"
        echo -e "\t-g: 基因文件夹名称，文件夹中可以包含一个或者多个基因组；"
        echo -e "\t-o: 文件前缀, 默认:去除假基因序列；"
        echo "`date '+Date: %D %T'`"
        echo -e "\n             >>>> mqy <<<< \n "
        echo "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo "     !!!!   保证所有文件都在一个工作目录下   !!!!"
        echo "     !!!!          -g:为文件夹名称           !!!!"
        echo "     !!!!    -i内容和-g中的文件名称要相同    !!!!"
        echo "     !!!!-g文件夹种存放的物种核苷酸基因组文件!!!!"
        echo "     !!!!，文件名称与-i文件中相同。          !!!!"
        echo "     !!!! 输出结果中包含:                    !!!!"
        echo "     !!!!               基因预测文件;        !!!!"
        echo "     !!!!               假基因预测文件;      !!!!"
        echo "     !!!!               去除假基因序列文件。 !!!!"
        echo "     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
        echo
        exit 0
        ;;
    i) 
        echo "输入-i文件为: $OPTARG"
        FILENAME=$OPTARG
        ;;
    g) 
        echo "输入-g文件为: $OPTARG"
        FOLDERNAME=$OPTARG
        ;;
    o) 
        echo "设置输出文件名前缀为: $OPTARG"
        OUTPREFIX=$OPTARG
        ;;
    *) 
        echo "未知参数: $opt"
        ;;
    esac
done

#####################-z 判断字符串是否为0

if [ -z $FILENAME ]
then
    echo -e "\n请输入-h，查看帮助文档！"
    exit 1
fi

############设置单个物种文件夹
for i in `cat ./$FILENAME` ;do mkdir ./$FOLDERNAME/${i}  &&  mv ./$FOLDERNAME/${i}.* ./$FOLDERNAME/${i} ;done

############Prokka
source /home/mengqingyao/miniconda3/bin/activate prokka

if [ -d "prokka" ]; then
    rm -rf prokka
fi

mkdir prokka && cd prokka

for i in `cat ../$FILENAME`; do prokka \
 ../$FOLDERNAME/${i}/${i}.* --cpus 0 --force --outdir ${i}_prokka --prefix ${i} --addmrna --compliant 200; done 

############pseudofinder

conda activate pseudofinder

if [ -d "../pseudogene" ]; then
    rm -rf ../pseudogene
fi

mkdir ../pseudogene && cd ../pseudogene

for i in `cat ../$FILENAME` ; do /home/mengqingyao/biosoftware/pseudofinder-master/pseudofinder.py \
annotate --genome ../prokka/${i}_prokka/${i}.gb* --outprefix ${i} \
-di --database /data/mengqy/database/diamond_nr/diamond_makedb.dmnd --threads 16 -e 1e-15 ;done

#############获得去除假基因的序

if [ -d "../remove_pseudogene" ]; then
    rm -rf ../remove_pseudogene
fi

mkdir ../remove_pseudogene && cd ../remove_pseudogene

for i in `cat ../$FILENAME` ; do cat ../prokka/${i}_prokka/${i}.ffn | grep ">" | cut -c 2- | cut -d " " -f 1 > ${i}_cds.ids && \
 awk '{if(match($0,"old_locus_tag=")) {print substr($0,RSTART+RLENGTH) }}' ../pseudogene/${i}_pseudos.gff \
 | sed 's/,/\t/g' | awk '{print $1}' | sed '/^$/d' >> ${i}_pseudofene.ids && comm -23 <(sort ${i}_cds.ids) <(sort ${i}_pseudofene.ids) > ${i}_remove.ids && \
 python3 /home/mengqingyao/Script_mqy/get_ids_fasta.py ../prokka/${i}_prokka/${i}.faa ${i}_remove.ids ${i}_remove_ids.faa && \
 python3 /home/mengqingyao/Script_mqy/get_ids_fasta.py ../prokka/${i}_prokka/${i}.ffn ${i}_remove.ids ${i}_remove_ids.ffn;done

if [ -d "../$OUTPREFIX" ]; then
    rm -rf ../$OUTPREFIX
fi

mkdir ../$OUTPREFIX 

for i in `cat ../$FILENAME` ; do mv ${i}_remove_ids.faa ../$OUTPREFIX && mv ${i}_remove_ids.ffn ../$OUTPREFIX; done

#  运行结束
echo "   运行结束！"
echo "   感谢使用本脚本！"
echo "   版本信息: v1.0"
echo "   日期: `date '+%Y-%m-%d %H:%M:%S'`"
echo "   作者: 孟庆瑶"
echo "   如有疑问请联系: <15877464851@163.com>"
{% endhighlight %}

> **If you want to use this script, you need to have the following software installed:**
> - [Prokka][prokka-doc]
> - [Pseudofinder][pseudofinder-doc]

<div class="notice">
  <h4>Since the software is installed in a different way and in a different location, we need to make changes.</h4>
</div>

**A python script in the process. get_ids_fasta.py**

{% highlight python %}
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys

script_name = sys.argv[0]

def usage():
    print('Usage: python {} [fasta_file] [namelist_file] [outfile_name]'.format(script_name))


def main():
    outf = open(sys.argv[3], 'w')
    dict = {}
    with open(sys.argv[1], 'r') as fastaf:
        for line in fastaf:
            if line.startswith('>'):
                name = line
                dict[name] = ''
            else:
                dict[name] += line.replace('\n', '')  # 读取整个fasta文件构成字典

    with open(sys.argv[2], 'r') as listf:
        for row in listf:
            row = row.strip()
            for key in dict.keys():  # 选取包含所需名称的基因名和序列
                if row in key:
                    outf.write(key)
                    outf.write(dict[key] + '\n')
    outf.close()


try:
    main()
except IndexError:
    usage()
{% endhighlight %}

**The script is used to extract the sequences of the specified gene names from the fasta file.**

## Running the workflow

{% highlight bash %}
nohup bash /home/mengqy/ws/No_pseudogene_faa.sh -i head.txt -g Genome/ & 
{% endhighlight %}

**Wait for the completion of the job.**

## Result

<div style="text-align: center; margin-bottom: 20px;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-25-genome-annotation-3.png"/>
</div>

We won't explain the results of the different software in detail here, but if you're interested, you can visit [Prokka][prokka-doc] and [Pseudofinder][pseudofinder-doc] to find out.

## Quote

> - [Prokka][prokka-doc]
> - [Pseudofinder][pseudofinder-doc]

> Email me with more questions!
> 584338215@qq.com

[prokka-doc]: https://github.com/tseemann/prokka
[pseudofinder-doc]: https://help.rc.ufl.edu/doc/Pseudofinder
