
# gmxMeta说明文档


## 命令行
```
docker run --rm -v {工作目录}:/data/ gmxmeta python /gmxMeta/gmxMeta.py -i1 {R1端fastq} -i2 {R2端fastq} -n {样本ID} -N {会员ID} -db {样本最终结果目录} -o {样本报告结果目录} --tmp {文件暂存目录} -DB {参考数据库目录} -t {CPU数量}
```

举例：

```
docker run --rm -v /home/hdd/yangk_test/:/data/ gmxmeta python /gmxMeta/gmxMeta.py -i1 /data/test_data/M202302160067_test_1.fastq.gz -i2 /data/test_data/M202302160067_test_2.fastq.gz -n M202302160067_test -N nlwushubin -db /data/db/ -o /data/test_out/M202302160067_test -DB /data/metaDB/ --tmp /data/tmp --debug -t 16 1>/home/hdd/yangk_test/M202302160067_test.log 2>/home/hdd/yangk_test/M202302160067_test.log
```

## 目录结构

数据分析及报告结果产生过程中，将使用多个目录，如下表所示。

|目录|作用|内容|涉及参数|
|-|-|-|-|
|db|所有样本结果|`old`、`ID.csv`、`sample_dir`|-db|
|metaDB|参考数据库|`CARD`、`chocophlan`、`metaphlan`、`uniref`、`VFDB`|-DB|
|tmp|分析过程文件暂存位置||--tmp|
|out|报告结果文件|`1.info.xlsx`、`4.alph.png`|-o|




