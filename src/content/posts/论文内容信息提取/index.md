---
title: 论文内容提取
published: 2024-09-18
description: "使用一个python代码对pdf论文进行粗略的不同信息提取提取"
image: "./121631681_p0.jpg"
tags: [python,论文]
category: "论文"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 准备好的文章 `pdf`

## 信息提取

abstracts = []
introductions = []
methods = []
others = []

```python 
import os
import re
from PyPDF2 import PdfReader

def extract_text_by_section(pdf_text, start_pattern, end_pattern):
    start = re.search(start_pattern, pdf_text, re.IGNORECASE)
    end = re.search(end_pattern, pdf_text, re.IGNORECASE)
    if start and end:
        return pdf_text[start.start():end.start()].strip()
    elif start:
        return pdf_text[start.start():].strip()
    return ""

def process_pdf_files(folder_path, output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)  # 创建输出目录

    abstracts = []
    introductions = []
    methods = []
    others = []
    
    for filename in os.listdir(folder_path):
        if filename.endswith('.pdf'):
            file_path = os.path.join(folder_path, filename)
            reader = PdfReader(file_path)
            full_text = "\n".join([page.extract_text() for page in reader.pages if page.extract_text()])

            title = full_text.split('\n')[0]
            abstract = extract_text_by_section(full_text, "abstract", "introduction")
            introduction = extract_text_by_section(full_text, "introduction", "method")
            method = extract_text_by_section(full_text, "method", "results|discussion|conclusion")

            # Calculating other parts excluding references
            references_start = re.search("references|bibliography", full_text, re.IGNORECASE)
            if references_start:
                rest_text = full_text[:references_start.start()].strip()
            else:
                rest_text = full_text
                
            other_parts = rest_text.replace(abstract, '').replace(introduction, '').replace(method, '')
            
            # Append with titles
            abstracts.append(f"Title: {title}\n{abstract}")
            introductions.append(f"Title: {title}\n{introduction}")
            methods.append(f"Title: {title}\n{method}")
            others.append(f"Title: {title}\n{other_parts}")

    # Save results to txt files in the specified output folder
    with open(os.path.join(output_folder, "abstracts.txt"), "w") as f:
        f.write("\n\n".join(abstracts))
    with open(os.path.join(output_folder, "introductions.txt"), "w") as f:
        f.write("\n\n".join(introductions))
    with open(os.path.join(output_folder, "methods.txt"), "w") as f:
        f.write("\n\n".join(methods))
    with open(os.path.join(output_folder, "others.txt"), "w") as f:
        f.write("\n\n".join(others))

# Specify the path to the folder
folder_path = '/home/tim/hmt/文章/del_pdf'
output_folder = '/home/tim/hmt/文章/del_output'
process_pdf_files(folder_path, output_folder)

```