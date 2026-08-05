import os, subprocess

output="/home/rare/arlen/outputs"
small_variant_folder="/home/rare/arlen/outputs/small_variants/SNVs_filtered"
bam_folder="/mnt/diskrare/arlenb/08/aligned_reads/t2t"
SV_folder="/home/rare/arlen/outputs/Structural_variants/08/structural_variants/vcfs_filtered/t2t"
mnt="/mnt/diskrare/arlenb/08/hiphase_results"
os.makedirs(mnt, exist_ok=True)
reference="/home/rare/arlen/reference/chm13v22.fasta"

def hiphase_phasing (output, small_variant_folder, bam_folder, SV_folder, mnt, reference):
    out_bam=f"{mnt}/bamfiles"
    os.makedirs(out_bam, exist_ok=True)
    out_vcf=f"{mnt}/variants"
    os.makedirs(out_vcf, exist_ok=True)
    for x in os.listdir(SV_folder):
        if x.endswith("A.vcf.gz"):
            SV_file=os.path.join(SV_folder,x)
            basename=os.path.basename(SV_file).replace(".vcf.gz", "")
            small_variant_file=os.path.join(small_variant_folder, x)    
            bamfile=f"{os.path.join(bam_folder, basename)}.bam"
            #print(SV_file,small_variant_file,bamfile)
            program="/home/rare/programs/hiphase-v1.4.5-x86_64-unknown-linux-gnu/hiphase"
            command=[f"{program}",
                     "--threads", "32",
                     "--reference", f"{reference}",
                     "--bam", f"{bamfile}",
                     "--output-bam", f"{out_bam}/{basename}.bam",
                     "--vcf", f"{small_variant_file}",
                     "--output-vcf", f"{out_vcf}/{basename}.small.vcf.gz",
                     "--vcf", f"{SV_file}",
                     "--output-vcf", f"{out_vcf}/{basename}.SV.vcf.gz",
                     "--stats-file", f"{out_vcf}/{basename}.stats.csv",
                     "--blocks-file", f"{out_vcf}/{basename}.blocks.tsv",
                     "--summary-file", f"{out_vcf}/{basename}.summary.tsv"]
            subprocess.run(command, check=True)
            
hiphase_phasing (output, small_variant_folder, bam_folder, SV_folder, mnt, reference)
