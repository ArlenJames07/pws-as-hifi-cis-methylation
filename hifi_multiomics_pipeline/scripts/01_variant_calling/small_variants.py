import os, subprocess

def small_variants(output, reference, genome, force=False):
    aligned_reads = os.path.abspath(f"{output}/08/aligned_reads/{genome}")
    out = os.path.abspath(f"{output}/08/small_variants/{genome}")
    os.makedirs(out, exist_ok=True)  # Ensure output directory exists

    reference_name = os.path.basename(reference)
    dir_reference = os.path.abspath(os.path.dirname(reference))
    num_cores = os.cpu_count()

    bam_files = [f for f in os.listdir(aligned_reads) if f.endswith(".bam")]

    for file in bam_files:
        basename = os.path.splitext(file)[0]
        vcf_path = os.path.join(out, f"{basename}.vcf.gz")

        if force or not os.path.exists(vcf_path):
            print(f"Calling variants for {file}...")

            BIN_VERSION = "1.8.0"
            command = [
                "docker", "run", "--rm",
                "-v", f"{aligned_reads}:/input",
                "-v", f"{dir_reference}:/reference",
                "-v", f"{out}:/output",
                f"google/deepvariant:{BIN_VERSION}",
                "/opt/deepvariant/bin/run_deepvariant",
                f"--model_type=PACBIO",
                f"--ref=/reference/{reference_name}",
                f"--reads=/input/{file}",
                f"--output_vcf=/output/{basename}.vcf.gz",
                f"--output_gvcf=/output/{basename}.g.vcf.gz",
                f"--num_shards={num_cores}"
            ]
            subprocess.run(command, check=True)
        else:
            print(f"Skipping {file}, VCF already exists.")

#small_variants(output, reference_t2t, genome_t2t, force=False)
#small_variants(output, reference_hg38, genome_hg38, force=False)

output="/home/rare/arlen/outputs"
mnt="/mnt/diskrare/arlenb"
t2t="t2t"
def filter_SNVs(mnt, output, genome):
    folder_SNVs=f"{mnt}/08/small_variants/deepvariant_fixed/{genome}"
    out=f"{output}/small_variants/SNVs_filtered"
    os.makedirs(out, exist_ok=True)
    for x in os.listdir(folder_SNVs):
        if x.endswith(".vcf.gz"):
            file=os.path.join(folder_SNVs, x)
            basename=os.path.basename(file)
            command=f"sudo bcftools view -f 'PASS' {file}  -o {out}/{basename}"
            subprocess.run(command, check=True, shell=True)
            

#filter_SNVs(mnt, output, t2t)
