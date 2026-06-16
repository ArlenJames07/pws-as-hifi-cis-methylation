import os, subprocess


folder_genomes="/home/rare/data/genomes/bamfiles"
reference_t2t="/home/rare/arlen/reference/chm13v22.fasta"
mnt_diskrare="/mnt/diskrare/arlenb"
smrtlink="/home/rare/programs/smrtlink/smrtlink/smrtcmds/bin"
genome_t2t="t2t"
genome_hg38="hg38"
output="/home/rare/arlen/outputs"
tandem_repeats="/home/rare/arlen/reference/human_T2T_CHM13v2.trf.bed"

def aligment(folder_genomes, reference, mnt_diskrare, smrtlink, genome):
    out=f"{mnt_diskrare}/08/aligned_reads/{genome}"
    os.makedirs(out, exist_ok=True)
    for x in os.listdir(folder_genomes):
        if x.endswith(".bam"):
            file = os.path.join(folder_genomes, x)
            basename = os.path.basename(file).replace(".bam", "")
            cmd=f"{smrtlink}/pbmm2 align {reference} {file} --sort --preset HIFI --log-level INFO -j 32 {out}/{basename}.bam"
            subprocess.run(cmd, check=True, shell=True)
#aligment(folder_genomes, reference_t2t, mnt_diskrare, smrtlink, genome_t2t)


def discover_signatures(output, mnt_diskrare, smrtlink, genome, tandem_repeats):
    aligned_reads = f"{mnt_diskrare}/08/aligned_reads/{genome}"
    out = f"{output}/Structural_variants/signatures/{genome}"

    os.makedirs(out, exist_ok=True)

    for x in os.listdir(aligned_reads):
        if x.endswith(".bam"):
            file = os.path.join(aligned_reads, x)
            basename = os.path.basename(file).replace(".bam", "")
            cmd=f"smrtlink/psbv discover {file} {out}/{basename}.svsig.gz --tandem-repeats {tandem_repeats}"

            command = (
                f"{smrtlink}/pbsv discover {file} {out_file} "
                f"--tandem-repeats /home/rare/arlen/reference/human_T2T_CHM13v2.trf.bed"
            )
            print(f"Running: {command}")
            subprocess.run(command, shell=True, check=True)
            
#discover_signatures(output, mnt_diskrare, smrtlink, genome_t2t, tandem_repeats)


def structural_variants_call(output, reference, smrtlink, genome):
    signatures = f"{output}/Structural_variants/08/structural_variants/signatures/{genome}"
    out = f"{output}/Structural_variants/pbsv_discover/{genome}"
    os.makedirs(out, exist_ok=True)
    lista=[]
    for x in os.listdir(signatures):
        if x.endswith("bam.svsig.gz"):
            file = os.path.join(signatures, x)
            basename = os.path.basename(file).replace("bam.svsig.gz", "")
            lista.append(file)
    files=" ".join(lista)
    cmd=f"{smrtlink}/pbsv call  {reference_t2t} {files}  -j 32 {out}/SVs_rare.vcf"
    subprocess.run(cmd, check=True, shell=True)
            
structural_variants_call(output, reference_t2t, smrtlink, genome_t2t)


def filter_SVs(output, genome):
    SVs_folder = f"{output}/08/structural_variants/vcfs/{genome}"
    out = f"{output}/08/structural_variants/vcfs_filtered/{genome}"

    os.makedirs(out, exist_ok=True)

    for x in os.listdir(SVs_folder):
        if x.endswith("vcf"):
            file = os.path.join(SVs_folder, x)
            basename = os.path.basename(file).replace(".vcf", "")
            filtered_file = f"{out}/{basename}.vcf"

            if os.path.exists(filtered_file):
                print(f"Skipping {basename}, filtered VCF already exists.")
                continue

            program = "../programs/svpack/svpack"
            command = f"{program} filter --pass-only --min-svlen 50 {file} > {filtered_file}"
            print(f"Running: {command}")
            subprocess.run(command, shell=True, check=True)

#filter_SVs(output,genome_t2t)