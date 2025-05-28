
TRAIT="spleen_original_gldm_SmallDependenceHighGrayLevelEmphasis"
python prioritize_nearest_genes.py -trait ${TRAIT}
python gene_highest_pops_by_locus.py -trait ${TRAIT}
python merge_prioritized_genes.py -trait ${TRAIT}