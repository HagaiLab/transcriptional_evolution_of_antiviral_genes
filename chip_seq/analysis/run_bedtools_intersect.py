import os

species=['human','mouse']
promoters=['prom1000100','prom2000500']
histones=['H3K27ac','H3K4me3']
treat=['PIC4','UNST']

chip_seq_dir=r'/tzachi_storage/lilachs/my_project_from_22_11/chip_seq'

cmd_load='module load bedtools/bedtools2.25'
os.system(cmd_load)

cmd_cd='cd /tzachi_storage/lilachs/my_project_from_22_11/chip_seq/bedtools_intersect/'
os.system(cmd_cd)

for s in species:
    if s=='human':
        s_short='HS'
    elif s=='mouse':
        s_short='MM'
    for p in promoters:
        aft=int(p[-3:])
        bef=int(p[-7:-3])
        p_file_name=f'bef{bef}_aft{aft}'
        for h in histones:
            for t in treat:
                cmd00=f'bedtools intersect -a {chip_seq_dir}/create_promoters_files/{s}/BED_promoters_{s}_{p_file_name}_no_headers.bed -b {chip_seq_dir}/narrowPeak_files/merged_{s}/{s_short}_{h}_{t}_merged.bed > {s}/{p}/{h}/{t}peaks_intersect_00.txt'
                cmd03=f'bedtools intersect -a {chip_seq_dir}/create_promoters_files/{s}/BED_promoters_{s}_{p_file_name}_no_headers.bed -b {chip_seq_dir}/narrowPeak_files/merged_{s}/{s_short}_{h}_{t}_merged.bed > {s}/{p}/{h}/{t}peaks_intersect_03.txt'
                cmd05=f'bedtools intersect -a {chip_seq_dir}/create_promoters_files/{s}/BED_promoters_{s}_{p_file_name}_no_headers.bed -b {chip_seq_dir}/narrowPeak_files/merged_{s}/{s_short}_{h}_{t}_merged.bed > {s}/{p}/{h}/{t}peaks_intersect_05.txt'

                os.system(cmd00)
                os.system(cmd03)
                os.system(cmd05)




