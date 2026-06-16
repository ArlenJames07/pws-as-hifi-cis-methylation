import os, subprocess

def pbcpgtools():
    ruta_bam_alineados="/mnt/diskrare/arlenb/08/hiphase_results/bamfiles"
    programa="/home/rare/programs/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin/aligned_bam_to_cpg_scores"
    modelo="/home/rare/programs/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/models/pileup_calling_model.v1.tflite"
    ruta_salida="/home/rare/arlen/outputs/methylation/genomes_2"
    os.makedirs(ruta_salida, exist_ok=True)
    for x in os.listdir(ruta_bam_alineados):
        if x.endswith("08_1_A01_bc2043_001P.bam"):
            file=os.path.join(ruta_bam_alineados,x)
            basename=os.path.basename(file).replace(".bam","")
            cmd=f"{programa} --bam {file} --output-prefix {ruta_salida}/{basename} --model {modelo} --threads 32"
            subprocess.run(cmd, check=True, shell=True)
            
pbcpgtools()            
