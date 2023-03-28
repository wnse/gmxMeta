
# gmxMeta说明文档


## 命令行
```
docker run --rm -v {数据目录}:{数据目录} -v{报告结果目录}:{报告结果目录} -v /home/:/data_for_tmp -v /home/hdd/gmxMeta/:/data/ gmxmeta python /gmxMeta/gmxMeta.py -t {CPU数量} --tmp /data_for_tmp/tmp -db /data/db/ -DB /data/metaDB/ -i1 {R1端fastq} -i2 {R2端fastq} -o {样本报告结果目录}  -n {样本ID} -N {会员ID} 
```

举例：shell_cmd.xlsx

## 目录结构

数据分析及报告结果产生过程中，将使用多个目录，如下表所示。

|目录|作用|内容|涉及参数|
|-|-|-|-|
|db|所有样本结果|`old`、`ID.csv`、`sample_dir`|-db|
|metaDB|参考数据库|`CARD`、`chocophlan`、`metaphlan`、`uniref`、`VFDB`|-DB|
|tmp|分析过程文件暂存位置||--tmp|
|out|报告结果文件|`1.info.xlsx`、`4.alph.png`|-o|

## 正式使用前：
1. 更新db目录中的old目录(以防丢失样本结果)
2. 删除db目录中的`ID.csv`

