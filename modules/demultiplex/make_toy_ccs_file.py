# %%

import pysam

def make_toy_file(input_filename, output_filename, toy_row_size):
    bam_filename = '../ULSC17_Gloria_Iso_031121.hifi_reads.bam'
    # touse = [line.strip() for line in open('zmws_to_whitelist_codethon_toy.txt')]
    reader = pysam.AlignmentFile(input_filename, 'rb', check_sq=False)
    f = pysam.AlignmentFile(output_filename, 'wb', header=reader.header)
    i = 0
    for r in reader:
        i += 1
        if i > toy_row_size:
            break
        f.write(r)
    f.close()

for i in range(1,8):
    print(i)
    make_toy_file(f'/Volumes/sheynkman/projects/bone_proteogenomics/pacbio_ccs_analysis/ccs_bams/osteo_cell{i}.ccs.bam', f'toy_files/cell{i}_toy.ccs.bam', 100)


# %%
