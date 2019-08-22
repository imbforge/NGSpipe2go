FastqScreen_vars=[
    outdir : QC + "/fastqscreen",
    threads: Integer.toString(ESSENTIAL_THREADS),
    //the fastqscreen_conf defines your references, with these we will create a fastqscreen conf script and then run the fastqscreen
    //this could be e.g.
    conf   : "Human::BowtieIndex/genome,rRNA::mouse_all_rRNA.fasta",
    //fastqscreen additional param e.g. subset or bowtie /bowtie 2 parameters
    extra  : "--nohits --subset 100000"
]
