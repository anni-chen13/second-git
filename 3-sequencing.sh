####################################################################
#                基因学苑VIP课程（第4季）:测序原理                   #
#                            王通 2024/11                          #
####################################################################

#####################################################
#                      下载数据                 #
#####################################################
#=========================
#      aws下载数据       #
#=========================
#1 查看虚拟环境，*号表示当前环境
conda env list

#创建虚拟环境
mamba create -n aws
#激活虚拟环境
mamba activate aws
#安装软件
mamba install -y awscli
mamba install -y sra-tools
#退出虚拟环境
mamba deactivate

# https://www.ncbi.nlm.nih.gov/bioproject/PRJNA491657
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR14143424/SRR14143424 ./
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR14092160/SRR14092160
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR14092310/SRR14092310

#https://www.ncbi.nlm.nih.gov/bioproject/PRJNA381365
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486971/SRR12486971
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486972/SRR12486972
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486974/SRR12486974
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486978/SRR12486978
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486979/SRR12486979
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486981/SRR12486981
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486983/SRR12486983
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486988/SRR12486988
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486989/SRR12486989
aws s3 cp --no-sign-request s3://sra-pub-run-odp/sra/SRR12486990/SRR12486990

# 下载果蝇测序数据
aws s3 --no-sign-request sync s3://ont-open-data/contrib/melanogaster_bkim_2023.01/flowcells/D.melanogaster.R1041.400bps/D_melanogaster_1/20221217_1251_MN20261_FAV70669_117da01a/ D.melanogaster.R1041.400bps/ 

#=========================
#     使用循环下载      #
#=========================
#修改扩展名
ls -1 | while read i;do mv $i ${i}.sra;done;

#数据转换
ls -1 | while read i;do fasterq-dump ${i};done;

#压缩
pigz -p 24 *.fastq



#####################################################
#          三、illumina数据质控及过滤               #
#####################################################

#fastqc质控
mkdir illumina_qc
fastqc -f fastq -o illumina_qc/ illumina_1.fastq.gz illumina_2.fastq.gz

#fastp 数据过滤
fastp -i illumina_1.fastq.gz  -I illumina_2.fastq.gz -o clean.1.fq.gz -O clean.2.fq.gz -z 4 -q 20 -u 40 -n 10 -f 15 -t 15 -F 15 -T 15 -h fastp.html

#过滤完质控
mkdir illumina_clean
fastqc -f fastq -o illumina_clean/ clean.1.fq.gz  clean.2.fq.gz


#####################################################
#            五、pacbio数据质控及过滤               #
#####################################################
#fastqc质控
mkdir pacbio_qc/
fastqc -f fastq -o pacbio_qc/ Dros_ULI.CL_14Aug2020.trimmed.unique.good.fastq.gz

#过滤数据
filtlong --min_length 300 --min_mean_q 90 Dros_ULI.CL_14Aug2020.trimmed.unique.good.fastq |  gzip >pacbio.filtlong.fq.gz

#质控完过滤
mkdir pacbio_clean
fastqc -f fastq -o pacbio_clean/ pacbio.filtlong.fq.gz

#将bam转换为fastq
samtools fastq m84006_221229_002525_s1.hifi_reads.bam > m84006_221229_002525_s1.hifi_reads.fastq

#压缩文件
pigz -p 32 m84006_221229_002525_s1.hifi_reads.fastq

fastqc -f fastq -o revio/ -t 24 m84006_221229_002525_s1.hifi_reads.fastq.gz

fastqc -f fastq -o onso -t 24 HG002_Onso_R1.fastq.gz  HG002_Onso_R2.fastq.gz


#####################################################
#                  八、纳米孔碱基识别           #
#####################################################
#链接数据到当前目录下
ln -s /ifs1/Vip4Data/3.sequencing/data/fast5_files/ ./
ll
#1 Basecalling
guppy_basecaller -i fast5_files/ -s fastq_files --config dna_r9.4.1_450bps_fast.cfg -r 
#合并全部fastq并压缩
cat fastq_files/*.fastq | gzip >lambda.fastq.gz

# dorado碱基识别
dorado basecaller fast pod5/dna_r10.4.1_e8.2_400bps_4khz/ --emit-fastq --output-dir output

#basecalling同时拆分barcode
guppy_basecaller -i /ifs1/Vip4Data/3.sequencing/data/16s_fast5/ -s fastq --config \
    dna_r9.4.1_450bps_hac.cfg -r  -x cuda:all --trim_barcodes --barcode_kits SQK-RBK004 \
    --min_score 10

#basecalling与拆分barcode分布完成
#首先进行basecalling
guppy_basecaller -i /ifs1/Vip4Data/3.sequencing/data/16s_fast5/ -s fastq --config dna_r9.4.1_450bps_hac.cfg -r  -x cuda:all
#拆分barcode
guppy_barcoder -i  fastq -s barcode --trim_barcodes --barcode_kits SQK-RBK004 --min_score 10

#############################
#          数据质控与过滤    #
#############################
#激活nanoplot环境
conda activate nanoplot

#NanoPlot质控
NanoPlot --fastq nanopore.fastq.gz -o nanoplot
#过滤数据
filtlong --min_length 2000 --min_mean_q 90 nanopore.fastq.gz |  gzip >nanopore.filtlong.fq.gz
#过滤完质控
NanoPlot --fastq nanopore.filtlong.fq.gz -o clean
conda deactivate




