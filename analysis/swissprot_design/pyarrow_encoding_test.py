import pyarrow.parquet as pq, pyarrow as pa, os, sys
t = pq.read_table('v5_min_zstd.parquet')
ints = ['pident_x10','qcov_x10','scov_x10','aln_length','q_start','q_end','s_start','s_end']
enc = {c:'DELTA_BINARY_PACKED' for c in ints}
enc['rk']='DELTA_BINARY_PACKED'
enc['bitscore']='BYTE_STREAM_SPLIT'; enc['neg_log10_evalue']='BYTE_STREAM_SPLIT'
enc['protein_id']='DELTA_BYTE_ARRAY'
enc['sp_id']='DELTA_BINARY_PACKED'
variants = {
 'pa_zstd3_dict_default.parquet': dict(compression='zstd', compression_level=3, row_group_size=1_000_000),
 'pa_zstd3_delta_bss.parquet': dict(compression='zstd', compression_level=3, row_group_size=1_000_000, use_dictionary=False, column_encoding=enc),
 'pa_zstd9_delta_bss.parquet': dict(compression='zstd', compression_level=9, row_group_size=1_000_000, use_dictionary=False, column_encoding=enc),
}
for name, kw in variants.items():
    pq.write_table(t, name, **kw)
    md = pq.read_metadata(name)
    sizes = {}
    for rg in range(md.num_row_groups):
        r = md.row_group(rg)
        for i in range(r.num_columns):
            c = r.column(i); sizes[c.path_in_schema] = sizes.get(c.path_in_schema,0)+c.total_compressed_size
    print(name, os.path.getsize(name)/1e6, 'MB')
    for k,v in sorted(sizes.items(), key=lambda x:-x[1]): print('   ', k, round(v/1e6,2))
