import os
import gzip
import sys
import shutil
import subprocess
import argparse
import logging

from MetaRes2Excel import check_samples
from MetaRep2Excel import Rep2Excel

def runMeta(**args):
    def check_file(d, k, t):
        if k in d:
            if t == 'value':
                return 0
            if t == 'file':
                if os.path.isfile(d[k]):
                    return 0
                else:
                    return f"{d[k]} not exists"
            if t == 'dir':
                if os.path.isdir(d[k]):
                    return 0
                else:
                    try:
                        os.makedirs(d[k])
                        return 0
                    except Exception as e:
                        return f'mkdir {d[k]}: {e}'
        else:
            return f'{k} not defined'
        
    check_args = {'infastq1': 'file', 'infastq2': 'file', 
                  'sample_name': 'value', 'member_name': 'value',
                  'db_dir': 'dir', 'meta_db': 'dir', 
                  'out_dir': 'dir', 'shell': 'value', 'tmp_dir': 'dir', 'threads':'value'}
    for k,t in check_args.items():
        check_res = check_file(args, k, t)
        if check_res:
            logging.error(check_res)
            return None
    logfile_err = os.path.join(args['tmp_dir'],'log.err')
    logfile_out = os.path.join(args['tmp_dir'],'log.out')
    shell_path = os.path.realpath(os.path.split(__name__)[0])
    cmd =  ["sh", os.path.join(shell_path, args['shell']), 
            "-i1", args['infastq1'], 
            "-i2", args['infastq2'],
            "-o", args['tmp_dir'], 
            "-t", str(args['threads']),
            "-n", os.path.join(args['meta_db'], 'chocophlan'), 
            "-p", os.path.join(args['meta_db'], 'uniref'),
            # "-u", os.path.join(args['meta_db'], 'utility_mapping'),
            "-m", os.path.join(args['meta_db'], 'metaphlan'),
            "-r", os.path.join(args['meta_db'], 'CARD', 'card.json'),
            "-vf", os.path.join(args['meta_db'], 'VFDB', 'VFDB_setB_pro.dmnd'),
            ]

    logging.info(f"meta analysis:{cmd}")
    logging.info(f"{' '.join(cmd)}")

    try:
        with open(logfile_err, 'w') as ERR, open(logfile_out,'w') as OUT:
            p = subprocess.Popen(cmd, stdout=OUT, stderr=ERR)
            logging.info(f"meta analysis pid:{p.pid}")
            p.wait()
            return p.returncode
    except Exception as e:
        logging.error(f'run meta error {e}')
        return None

def gz(sourcefile, targetfile):
    '''压缩文件'''
    gz_filename = targetfile #压缩后文件名
    f_ungz = open(sourcefile,'rb') 
    f_gz = gzip.open(gz_filename,'wb') 
    # f_gz = gzip.GzipFile(gz_filename,'wb') 使用GzipFile类创建压缩文件对象
    f_gz.writelines(f_ungz)
    f_ungz.close()
    f_gz.close()


def mvFiles(source_dir, target, sample_id, member_id):
    target_dir = os.path.join(target, sample_id[:5], sample_id[:8], sample_id)
    if not os.path.isdir(target_dir):
        os.makedirs(target_dir)
    
    logfile = 'log.err'
    fastp_json = 'fastp.json'
    fastp_html = 'fastp.html'
    rgi_res_file = 'card.txt'
    vf_res_file = 'vfdb.diamond.txt'
    raw2orf_file = 'raw2orf.bam.idxstats'
    metaphlan_res_file = 'final_metaphlan_bugs_list.tsv'
    pathabundance_res_file = 'final_pathabundance_unstratified.tsv'
    final_contig_file = 'final.contigs.fa'

    file_list = [logfile, fastp_json, fastp_html, rgi_res_file, vf_res_file, raw2orf_file]
    for f in file_list:
        try:
            shutil.copy(os.path.join(source_dir, f), os.path.join(target_dir, f))
        except Exception as e:
            logging.error(f'copy error {f}:{e}')
    try:
        shutil.copy(os.path.join(source_dir, 'humann', pathabundance_res_file), os.path.join(target_dir, pathabundance_res_file))
        shutil.copy(os.path.join(source_dir, 'humann', 'final_humann_temp', metaphlan_res_file), os.path.join(target_dir, metaphlan_res_file))
    except Exception as e:
        logging.error(f'copy error {metaphlan_res_file}:{e}')
    try:
        gz(os.path.join(source_dir, 'megahit_out', final_contig_file), os.path.join(target_dir, final_contig_file+'.gz'))
    except Exception as e:
        logging.error(f'gz error {final_contig_file}:{e}')
    logging.info(f"mv files done")
    try:
        os.system(f"echo {member_id},{sample_id}  >>{os.path.join(target, 'ID.csv')}")
    except Exception as e:
        logging.error(f'echo ID error {member_id} {sample_id}:{e}')

def runGmx(final_out_excel, config_file, member_name, sample_name, db_dir, old_db_dir=None):
    if not old_db_dir:
        old_db_dir = os.path.join(db_dir, 'old')
        if not os.path.isdir(old_db_dir):
            logging.error(f"{old_db_dir} not exists")
            old_db_dir = None
    sample_res_dir = os.path.join(db_dir, sample_name[:5], sample_name[:8], sample_name)
    rep_dict_list = check_samples(config_file, db_dir, member_name, outputdir=sample_res_dir, old_db_dir=old_db_dir)
    Rep2Excel(final_out_excel, rep_dict_list)
    logging.info(f"excel result done")


if __name__ == "__main__":

    shell_path = os.path.realpath(os.path.split(__file__)[0])
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parse.add_argument('-i1', '--infastq1', help=f'input fastq R1', required=True)
    parse.add_argument('-i2', '--infastq2', help=f'input fastq R2', required=True)
    parse.add_argument('-n', '--sample_name', help=f'sample name (样本ID)', required=True)
    parse.add_argument('-N', '--member_name', help=f'member name (会员ID)', required=True)
    parse.add_argument('-db', '--db_dir', help=f'sample results dir', default='/sample_db/')
    parse.add_argument('-o', '--out_dir', help=f'output dir, default to db_dir', default=None)
    parse.add_argument('-DB', '--meta_db', help=f'metageome analysis database dir, including humann, metaphlan, CARD, VFDB etc', default='/meta_db/')
    parse.add_argument('--shell', help=f'meta shell name in {shell_path}', default=os.path.join(shell_path, 'run_humann.sh'))
    parse.add_argument('--tmp', help=f'tmp dir for analysis', default='/data/tmp')
    parse.add_argument('--config', help=f'function config', default=f"{os.path.join(shell_path, 'function_config.xlsx')}")
    parse.add_argument('--debug', help=f'delete tmp result', action='store_true')
    parse.add_argument('-t', '--threads', help=f'threads for analysis', default=10)

    args = parse.parse_args()  
    args.tmp_dir = os.path.join(args.tmp, args.sample_name)

    config_excel = args.config
    if not os.path.isfile(config_excel):
        sys.exit(f"{config_excel} not exists")

    if not os.path.join(args.tmp):
        logging.info(f'{args.tmp} not exstis, mkdir...')
        os.mkdirs(args.tmp)
    if os.path.isdir(args.tmp_dir):
        logging.info(f'{args.tmp_dir} exstis, remove...')
        shutil.rmtree(args.tmp_dir)

    if not args.out_dir:
        args.out_dir = os.path.join(args.db_dir, args.sample_name)

    logging.basicConfig(level=logging.INFO, encoding='utf-8',
                        filename=args.out_dir + '.log',
                        format='%(asctime)s %(levelname)s %(message)s',
                        datefmt="%Y-%m-%d %H:%m:%S")

    try:
        res = runMeta(**args.__dict__)
        mvFiles(args.tmp_dir, args.db_dir, args.sample_name, args.member_name)
        final_out_excel = os.path.join(args.out_dir, '1.info.xlsx')
        runGmx(final_out_excel, config_excel, args.member_name, args.sample_name, args.db_dir, old_db_dir=None)
    except Exception as e:
        logging.error(f"LAST:{e}")
    if not args.debug:
        try:
            shutil.rmtree(args.tmp_dir)
            logging.info(f"remove tmp {args.tmp_dir} {e}")
        except Exception as e:
            logging.error(f"remove tmp {args.tmp_dir} {e}")
