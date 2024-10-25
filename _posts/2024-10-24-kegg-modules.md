---
title: "KEGG Metabolic Modules！"
categories: 
 - Gene annotation
tags: 
 - KEGG
---

# Introduction

I'll use one of my previous studies as an example:

Now in the same species of bacteria with both symbiotic and free-living modes of life, we wanted to analyze them histologically to see if there were any differences in their metabolism.

## Just do it.

We need to use [KEGG][KEGG-docs] and [KEGG API][KAAS-docs] to collect some of the data we need.

### [KEGG]

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-1.png"/>
</div>

Collect all data from the site, save it in a file and modify it.

`"Sublime Text" is recommended`

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-2.png"/>
</div>

Next, we need to convert it to a dataframe format. Use this script:

{% highlight python %}
#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @Author    : smallmeng
# @Email     : 15877464851@163.com
# @Time      : 2023/4/20 8:57

import sys

class KeggModuleProcessor:
    def __init__(self, input_file, output_file):
        self.input_file = input_file
        self.output_file = output_file
        self.B_class = ""
        self.C_class = ""

    def descriptive(self):
        print('\n描述: 将KEGG网页的模块分类自制文件输入，变成矩阵形式。\n')

    def usage(self):
        print('Usage: python3 Kegg_module_name_class_one.py [input_file] [outfile]')

    def process_file(self):
        with open(self.input_file, 'rt') as kaas, open(self.output_file, 'wt') as moules_anno:
            moules_anno.write('B_class\tC_class\tm_id\tm_name\n')

            for line in kaas:
                line = line.strip()
                if line[0] == 'A' and len(line) > 1:
                    self.B_class = line[4:len(line)]
                elif line[0] == 'B' and len(line) > 1:
                    self.C_class = line[5:len(line)]
                elif line[0] == 'C' and len(line) > 1:
                    m_id = line[6:12]
                    m_name = line[13:len(line)]
                    moules_anno.write(f'{self.B_class}\t{self.C_class}\t{m_id}\t{m_name}\n')

def main():
    try:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
        processor = KeggModuleProcessor(input_file, output_file)
        processor.process_file()
    except IndexError:
        processor = KeggModuleProcessor("", "")
        processor.descriptive()
        processor.usage()

if __name__ == '__main__':
    main()
{% endhighlight %}


{% highlight bash %}
python Kegg_module_name_class_one.py Kegg_module_name_class.txt Kegg_module_name_class_one.txt
{% endhighlight %}

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-3.png"/>
</div>

### KEGG API

The information corresponding to the KO number and the module can be obtained with [KEGG API][API-docs] and saved to a file.

**Entrances** [Click here][KO-moudles]

Very simple processing operation.

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-4.png"/>
</div>

### Download Genome Annotation File

The previous work was as follows:

> [KEGG ANN][KEGG-ann]

Wait for results and then download the file.

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-5.png"/>
</div>

- Save the first file.

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-6.png"/>
</div>

- Second document.

<div style="text-align: center; margin-bottom: 20px;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-7.png" title="one"/>
</div>

<div style="text-align: center; margin-bottom: 20px;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-8.png" title="two"/>
</div>

<div style="text-align: center;">
  <img src="https://mengqy2022.github.io/assets/images/20241025-9.png"  title="three"/>
</div>

<div class="notice">
  <h4>The required documents are ready.</h4>
</div>

> Email me with more questions!
> 584338215@qq.com

[KEGG-docs]: https://www.kegg.jp/brite/ko00002a
[API-docs]: https://www.kegg.jp/kegg/rest/
[KO-moudles]: https://rest.kegg.jp/link/module/ko
[KEGG-ann]: https://mengqy2022.github.io/gene%20annotation/kegg-ann/