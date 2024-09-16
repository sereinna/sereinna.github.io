---
title: 如何筛选一个gostar数据库并从中提取蛋白-配体作用信息
published: 2024-09-17
description: "Linux系统，Del数据，docking"
image: "./7.jpg"
tags: [Linux,del,docking]
category: "代码"
draft: false
lang: ""
---

## 筛选代码
```python 
import pymysql
import pandas as pd


class MysqlSave:
    def __init__(self):
        self.content = pymysql.Connect(
            host='127.0.0.1',  # mysql的主机ip
            port=3306,  # 端口
            user='root',  # 用户名
            passwd='serein',  # 数据库密码
            charset='utf8',  # 使用字符集
        )
        self.cursor = self.content.cursor()

    def search_and_save(self,sql,csv_file):
        """
        导出为csv的函数
        :param sql: 要执行的mysql指令
        :param csv_file: 导出的csv文件名
        :return: 
        """
        # 执行sql语句
        self.cursor.execute(sql)
        # 拿到表头
        des = self.cursor.description
        title = [each[0] for each in des]

        # 拿到数据库查询的内容
        result_list = []
        for each in  self.cursor.fetchall():
            result_list.append(list(each))

        # 保存成dataframe
        df_dealed = pd.DataFrame(result_list, columns=title)
        # 保存成csv 这个编码是为了防止中文没法保存，index=None的意思是没有行号
        df_dealed.to_csv(csv_file, index = None, encoding ='utf_8_sig')
        
        
if __name__ == '__main__':
    mysql = MysqlSave()
    sql = '''
        SELECT DISTINCT
            act_id,
            std_activity_type,
            activity_uom,
            activity_value,
            assay_type,
            sub_smiles,
            mol_weight,
            all_activity_gostar.pdb_id,
            uniprot_id
        FROM
            gostar.all_activity_gostar,
            gostar.structure_details,
            gostar.drug_targets,
            gostar.standard_name_master,
            gostar.target_protein_master
        WHERE
            gostar.all_activity_gostar.gvk_id = gostar.structure_details.gvk_id
            AND gostar.all_activity_gostar.stdname_id = gostar.standard_name_master.stdname_id
            AND gostar.all_activity_gostar.target_id = gostar.target_protein_master.target_id
            AND activity_type = 'IC50'
            AND activity_type = 'Ki' or 'Kd'
            AND activity_prefix = "="
            AND assay_type = "B"
            AND activity_uom = "nM"
            AND sub_smiles IS NOT NULL
            AND uniprot_id IS NOT NULL
            AND all_activity_gostar.pdb_id IS NOT NULL
            AND mol_weight < 800
            #AND activity_value <= 100

        #LIMIT 50;
    '''
    
    #file_path = r'/mnt/c/Users/86173/rmsd/mysql3.csv'
    mysql.search_and_save(sql, 'mysql_kikd_100.csv')

    ```