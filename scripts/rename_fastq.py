import os
import glob

dir_name = 'fastq_files'
for filename in glob.glob(f'{dir_name}/*'):
    sra_code = filename.split('/')[1].split('_')[0]
    new_name = f'{dir_name}/{sra_code}.fastq.gz'
    os.rename(filename, new_name)
    print(sra_code)
