---
title: Selenium获取cyclicpepedia信息
published: 2025-01-13
description: "如何使用无头动态网页Selenium+ChromeDriver，通过输入smiles，爬取获从clicpepeycdia获得的氨基酸的序列信息"
image: "./124366833_p0.png"
tags: [Linux,氨基酸,web，python]
category: "指南"
draft: false
lang: ""
---

## 前言

:::tip
准备工作 chrome网页端 `ChromeDriver `，动态网页 `Selenium`\

:::

## 如何安装ChromeDriver

chrome 下载 ChromeDriver 

打开 `https://googlechromelabs.github.io/chrome-for-testing/#stable`


移动其到chrome安装位置 C:\Program Files (x86)\Google\Chrome\Application

设置:此电脑→右击属性→高级系统设置→环境变量→用户变量→Path→编辑→新建

path C:\Program Files (x86)\Google\Chrome\Application\

## 如何安装Selenium

```bash
pip install selenium
```

## 爬取代码
```python
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options

def get_sequence_from_smiles(smiles_value):
    """
    根据输入的 SMILES 值，获取对应的 Amino Acid Sequence（无头模式）。
    
    Args:
        smiles_value (str): 输入的 SMILES 字符串。
    
    Returns:
        str: 提取的 Amino Acid Sequence 值。如果出错，返回 None。
    """
    # 配置无头模式的 Chrome 浏览器选项
    chrome_options = Options()
    chrome_options.add_argument("--headless")  # 无头模式
    chrome_options.add_argument("--disable-gpu")  # 禁用 GPU（可选）
    chrome_options.add_argument("--no-sandbox")  # 防止沙盒问题（推荐服务器上使用）
    chrome_options.add_argument("--disable-dev-shm-usage")  # 防止资源限制错误（推荐服务器上使用）

    # 指定 chromedriver.exe 的路径
    service = Service("C:/Program Files/Google/Chrome/Application/chromedriver.exe")  # 替换为你的实际路径

    # 创建无头模式的 WebDriver 实例
    driver = webdriver.Chrome(service=service, options=chrome_options)

    try:
        # 打开目标网页
        url = "https://www.biosino.org/iMAC/cyclicpepedia/stru2seq"
        driver.get(url)

        # 等待 textarea 元素加载（最长等待 10 秒）
        textarea = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.XPATH, '/html/body/div[1]/div[3]/div/div/div/div/div/div/div/div[1]/div[2]/form/div[1]/div/textarea'))
        )

        # 输入 SMILES 数据
        textarea.clear()  # 清空输入框（如果需要）
        textarea.send_keys(smiles_value)  # 输入数据

        # 定位并点击提交按钮
        submit_button = WebDriverWait(driver, 10).until(
            EC.element_to_be_clickable((By.XPATH, '/html/body/div[1]/div[3]/div/div/div/div/div/div/div/div[1]/div[2]/form/div[2]/div[2]/button'))
        )
        submit_button.click()  # 点击按钮

        # 等待结果生成
        # 定位父级 <p> 标签
        result_element = WebDriverWait(driver, 10).until(
            EC.presence_of_element_located((By.XPATH, '/html/body/div[1]/div[3]/div/div/div/div/div/div/div/div[1]/div[3]/p[6]'))
        )

        # 提取目标值
        result_text = result_element.text

        # 去掉 <b> 标签的内容，只保留实际的值
        cleaned_result = result_text.replace("Amino acid sequence :", "").strip()

        # 返回结果
        return cleaned_result

    except Exception as e:
        print(f"出现错误: {e}")
        return None

    finally:
        # 关闭 WebDriver
        driver.quit()
```

## 封装使用
```python

import pandas as pd
import time
from concurrent.futures import ThreadPoolExecutor

# 定义批量处理函数
def process_smiles_with_threads(input_file, output_file, max_threads=8):
    """
    使用多线程处理 SMILES 列，并获取对应的 Sequence。
    
    Args:
        input_file (str): 输入的 CSV 文件路径。
        output_file (str): 输出的 CSV 文件路径。
        max_threads (int): 最大线程数。
    """
    # 读取 CSV 文件
    df = pd.read_csv(input_file)

    # 确保存在 'smiles' 列
    if 'smiles' not in df.columns:
        raise ValueError("输入文件中必须包含 'smiles' 列！")

    # 如果不存在 'sequence' 列，则添加该列
    if 'sequence' not in df.columns:
        df['sequence'] = None

    # 过滤未处理的行
    rows_to_process = df[df['sequence'].isna()]

    if rows_to_process.empty:
        print("所有行均已处理，无需重复操作。")
        return

    # 定义处理单行的函数
    def process_row(index, smiles_value):
        try:
            print(f"正在处理第 {index + 1} 行的 SMILES: {smiles_value}")
            sequence = get_sequence_from_smiles(smiles_value)  # 调用你的函数
            print(f"第 {index + 1} 行处理完成，结果为: {sequence}")
            return index, sequence
        except Exception as e:
            print(f"第 {index + 1} 行处理出错: {e}")
            return index, None

    # 使用线程池并行处理
    with ThreadPoolExecutor(max_threads) as executor:
        futures = {
            executor.submit(process_row, index, row['smiles']): index
            for index, row in rows_to_process.iterrows()
        }

        # 获取每个线程的执行结果
        for future in futures:
            index, sequence = future.result()
            df.at[index, 'sequence'] = sequence  # 将结果写入 DataFrame

            # 实时保存中间结果到文件
            df.to_csv(output_file, index=False)

    print(f"处理完成，结果已保存到 {output_file}")


if __name__ == "__main__":
    # 输入和输出文件路径
    input_csv = "a.csv"  # 输入文件路径
    output_csv = "output.csv"  # 输出文件路径

    # 调用多线程批量处理函数
    process_smiles_with_threads(input_csv, output_csv, max_threads=4)
```

