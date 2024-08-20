import subprocess
import os
import docker


def run_docker(path_to_fastq,samples,path_to_shuffled,runDate,snp_output_mount,path_to_referance):
#https://forums.docker.com/t/how-is-var-run-docker-sock-created-on-mac-if-docker-daemon-runs-in-hyperkit-vm-on-linux/29090/2
    client = docker.from_env()

    #get command
    paired_end_reads = join_paired_end_reads(samples) 

    #start up continer to mix reads
    client.containers.run("staphb/lyveset:latest",volumes=[path_to_shuffled+":/data/CRAB",path_to_fastq+"/fastp:/data/FASTQ"],command="/bin/bash -c '"+paired_end_reads+"'")
    print("-"*50+"\n")
    print("Shuffled Reads Docker has exited")
    print("-"*50+"\n")
    #get snp  string
    #snp_command_string = run_SNPCreation(samples,runDate)
    snp_command_string="set_manage.pl --create /data/Output/"+runDate+" && set_manage.pl /data/Output/"+runDate+" --change-reference /data/referance/GCF_008632635.1_ASM863263v1_referance_genome.fna && cd /data/Output/"+runDate+"/reads && ln -sv /data/CRAB/*.fastq.gz .&& launch_set.pl --numcpus 20 /data/Output/"+runDate

    #startup SNP container to run analysis
    #client.containers.run("staphb/lyveset:latest",volumes=[path_to_shuffled+":/data/CRAB",path_to_referance+":/data/referance",snp_output_mount+":/data/Output"],command="/bin/bash -c '"+snp_command_string+"'")


def join_paired_end_reads(samples):
    #path to fastq should put us in the fastp folder created with each sample 
    #path_to_allreads/fastp/sample
    docker_string=""

    for sample in samples:
        docker_string+="cd /data/FASTQ/"+sample+" &&  shuffleSplitReads.pl --numcpus 15 "+sample+"_R* -o /data/CRAB/ && "
    
        #docker_string+="shuffleSplitReads.pl --numcpus 15 "+sample+"_R* -o /data/CRAB/ && "
    
    #print("docker run -v "+path_to_shuffled+":/data/CRAB -v "+path_to_fastq+"/fastp:/data/FASTQ staphb/lyveset:1.1.4f /bin/bash -c '"+docker_string+"'")
    #subprocess.run("docker run -v "+path_to_shuffled+":/data/CRAB -v "+path_to_fastq+"/fastp::/data/FASTQ staphb/lyveset:1.1.4f /bin/bash -c '"+docker_string+"'",shell =True)
    return docker_string[:-3]

#def run_SNPCreation(path_to_shuffled_reads,samples,output_mount,run_date,path_to_referance):

def run_SNPCreation(samples,run_date):
    #set_manage.pl --create /data/Output/031124 && set_manage.pl /data/Output/031124 
    #--change-reference /data/referance/GCF_008632635.1_ASM863263v1_referance_genome.fna 

    #&& cd /data/Output/031124/reads && ln -sv /data/CRAB/*.fastq.gz . 
    #&& launch_set.pl --numcpus 20 /data/Output/031124'
    docker_string="set_manage.pl --create /data/Output/"+run_date+" && set_manage.pl /data/Output/"+run_date+" --change-reference /data/referance/GCF_008632635.1_ASM863263v1_referance_genome.fna && "

    for sample in samples:
        docker_string+="set_manage.pl /data/Output/"+run_date+"  --add-reads /data/CRAB/"+sample+".fastq.gz && "

    #subprocess.run("docker run -v "+path_to_shuffled_reads+":/data/CRAB -v "+path_to_referance+":/data/referance -v "+output_mount+":/data/Output staphb/lyveset:1.1.4f /bin/bash -c '"+docker_string+"launch_set.pl /data/Output"+run_date+"'",shell=True)
    #print("docker run -v "+path_to_shuffled_reads+":/data/CRAB -v "+path_to_referance+":/data/referance -v "+output_mount+":/data/Output staphb/lyveset:1.1.4f /bin/bash -c '"+docker_string+"launch_set.pl --numcpus 8 /data/Output"+run_date+"'")
#launch_set.pl --numcpus 20 /data/Output/ -ref /data/referance/GCF_008632635.1_ASM863263v1_referance_genome.fna 
    return docker_string+"launch_set.pl --numcpus 20 /data/Output/"+run_date


if __name__ == "__main__":
    run_docker("something",['2278019', '2278016', '2281037', '2281793'],"/Users/adrian/Desktop/CRAB_TESTING/Assembled","062422","/Users/adrian/Desktop/CRAB_TESTING","/Users/adrian/Desktop/CRAB_TESTING/temp")


