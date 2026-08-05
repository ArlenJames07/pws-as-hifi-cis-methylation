import os, subprocess

folder_bam_files="/mnt/diskrare/arlenb/08/hiphase_phased_bam/t2t"
reference_file="/home/rare/arlen/reference/chm13v22.fasta"
folder_snvs="/mnt/diskrare/arlenb/08/structural_variants/hiphase_phased_variants/t2t"
exclude_regions="/home/rare/arlen/reference/cnv.exclude.bed"

def hifi_cnv(folder_bam_files,reference_file, folder_snvs, exclude_regions):
    out=f"/home/rare/arlen/outputs/Structural_variants/hifi_cnv"
    os.makedirs(out,exist_ok=True)
    for x in os.listdir(folder_bam_files):
        if x.endswith(".bam"):
            file=os.path.join(folder_bam_files,x)
            basename=x.replace(".bam","")
            cmd=f"hificnv --bam {file} --ref {reference_file} --maf {folder_snvs}/{basename}.small.vcf.gz --exclude {exclude_regions} --threads 32 --output-prefix {out}/{basename}"
            subprocess.run(cmd, check=True, shell=True)
            #print(file)

hifi_cnv(folder_bam_files,reference_file, folder_snvs, exclude_regions)            


