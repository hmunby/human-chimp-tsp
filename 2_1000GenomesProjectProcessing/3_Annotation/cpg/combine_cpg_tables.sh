out_file="cpg_annot_table.txt"

for i in {1..22}
do
    if [ $i -eq 1 ]
    then
        cp cpg_tables_per_chrom/cpg.annot_table.chr${i}.txt $out_file
    else
        cat cpg_tables_per_chrom/cpg.annot_table.chr${i}.txt >> $out_file
    fi
done

# bgzip the complete annotation table
bgzip $out_file

tabix -s1 -b2 -e2 ${out_file}.gz



