---
title: "Eggnog-mapper Gene Function Annotation！"
categories: 
 - Genomics
tags: 
 - KEGG
 - COG
 - GO
---

# Introduction

Gene function annotation, in simple terms, is the comparison of protein sequences extracted from the genome based on existing databases to obtain the corresponding information.

- [eggNOG][eggNOG-doc]
1. The whole genome protein sequences of 5090 organisms (477 eukaryotes, 4445 representative bacteria and 168 archaea) and 2502 viruses were collected in the eggNOG public database version V5.0. 
2. These species were grouped into 379 categories (taxonomic levels), with the number of each category indicated by the NCBI classification number. 
3. Homologous gene classification of the whole genome protein sequences of species in each category yielded a total of 4.4M homologous gene classes. eggNOG database is an extension of NCBI's COG database, which collects a more comprehensive set of species and a larger amount of protein sequence data.
4. The eggNOG database provides a comprehensive set of gene function annotations, including the functional categories of each gene (COG), the gene ontology (GO) terms, and the KEGG pathway information.

Previously we obtained protein coding sequences through the bacterial genome [Click][ga-doc], and now we want to construct a single-copy genome phylogenetic tree to analyze the evolutionary relationships among bacteria.

## Go to the next step

### Obtain the protein coding sequences

- First, we can download the bacterial protein coding sequences from any public databases.
- Secondly, genome annotation was performed to obtain protein coding sequences.

##  eggNOG-mapper Gene function annotation

### Online Annotation

- [eggNOG][eggNOG-doc]

<div style="text-align: center; margin-bottom: 20px">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-26-eggnog-mapper-1.png"/>
</div>

Check if eggnog-mapper is installed.
1. Select sequence type
2. Upload the protein coding sequence file
3. Enter your e-mail address
4. `Search filters` and `Annotation options` can be adjusted as needed.
5. Click `Submit` to start the annotation process.

Wait for the completion of the job.

Download the result file for your email.

### Local Annotation

1. Need download the [eggNOG-mapper][eggNOG-mapper-doc] software.
2. Download the [eggNOG database][eggNOG-database].
3. Run the eggNOG-mapper software.

- Check if eggnog-mapper is installed.

{% highlight bash %}
emapper.py -v
{% endhighlight %}

<div style="text-align: center">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-26-eggnog-mapper-2.png"/>
</div>

- Download the [eggNOG database][eggNOG-database].

<div style="text-align: center; margin-bottom: 20px">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-26-eggnog-mapper-3.png"/>
</div>

You need set the database path in the `emapper.py` configuration file.

- The software requires [diamond][diamond-doc] support.

{% highlight bash %}
diamond help
{% endhighlight %}

<div style="text-align: center; margin-bottom: 20px">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-26-eggnog-mapper-4.png"/>
</div>

Diamond database

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-26-eggnog-mapper-5.png"/>
</div>

- Run the eggNOG-mapper software.

{% highlight bash %}
nohup python3 emapper.py -m diamond -i input.fa -o output.txt --tax_scope all --cpu 0 &
{% endhighlight %}

## results

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/2024-10-26-eggnog-mapper-6.png"/>
</div>

## Quote 

> - [eggNOG][eggNOG-doc]
> - [diamond][diamond-doc]

> Email me with more questions!
> 584338215@qq.com

[eggNOG-doc]: http://eggnog-mapper.embl.de/
[ga-doc]: https://mengqy2022.github.io/genomics/genome-annotation/
[eggNOG-mapper-doc]: https://github.com/eggnogdb/eggnog-mapper
[eggNOG-database]: http://eggnog5.embl.de/download/emapperdb-5.0.2/
[diamond-doc]: https://github.com/bbuchfink/diamond