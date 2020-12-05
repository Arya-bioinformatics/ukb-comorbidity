../magma/magma_v1.07b/magma --annotate window=2,0.5 --snp-loc ../genome/varaint_info.txt --gene-loc ../magma/NCBI37.3/NCBI37.3.gene.loc --out ../genome/magma_file/snp

../magma/magma_v1.07b/magma --annotate window=2,0.5 --snp-loc ../genome/varaint_info_rmhla.txt --gene-loc ../magma/NCBI37.3/NCBI37.3.gene.loc --out ../genome/magma_file/snp_rmhla


# batch
'../magma/magma_v1.07b/magma ' \
               '--bfile ../magma/g1000_eur/g1000_eur ' \
               'synonyms=../magma/g1000_eur/g1000_eur.synonyms ' \
               '--gene-annot ../genome/magma_annotate/snp.genes.annot ' \
               '--pval ../genome/magma_file/' + each + \
               ' N=452264 duplicate=first --out ../genome/magma_annotate/' + key
