import pandas as pd

# load a list of house-keeping genes
hk_df = pd.read_csv('references/mm10.HouseKeepingGene.bed', sep='\t', header=None)
hk_df.columns = ['chr', 'start', 'end', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'ItemRgb', 'blockCount', 'blockSizes', 'blockStarts']
hk_df = hk_df.set_index('name')
refseq_id = list(hk_df.index.unique())

# load gtf file to get chr positions
gtf  = pd.read_csv("references/genome.gtf", sep='\t', header=None, skiprows=5, usecols=[0, 3, 4, 8], dtype='str')
gtf.columns = ['chr_ind', 'chrstart', 'chrend', 'Gene stable ID']
gtf['Gene stable ID'] = gtf['Gene stable ID'].map(lambda x: x.split(';')[0].replace('gene_id "', '').replace('"', '') if 'transcript_id' not in x else None)
gtf.dropna(inplace=True)

# load biomart annotation to get refseq ids
biomart_df = pd.read_csv("references/biomart_export.txt", sep='\t')
biomart_hk = biomart_df[biomart_df["RefSeq mRNA ID"].isin(refseq_id)]
biomart_hk = biomart_hk.merge(gtf, on='Gene stable ID').set_index('RefSeq mRNA ID')

# create a bed file
merge_col = ['Gene stable ID', 'chrstart', 'chrend', 'chr_ind']
merge_df = hk_df.merge(biomart_hk[merge_col], left_index=True, right_index=True).reset_index()
bed_file = merge_df[['chr_ind', 'chrstart', 'chrend', 'Gene stable ID', 'score', 'strand', 'thickStart', 'thickEnd', 'ItemRgb', 'blockCount', 'blockSizes', 'blockStarts']]
bed_file.to_csv('references/hk_genes_ensembl.bed', sep='\t', header=None, index=None)
