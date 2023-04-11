#CRAB Pipeline
from WF_0_Assembler.WF_0_Assembler_runner import run_assembly
from WF_1_Annotate.WF_1_Annotate import run_annotate
from WF_2_FindAMR.WF_2_FindAMRs import find_AMR_genes
from WF_3_DB.WF_3_DB_push import run_DB_push
from WF_3_5_SNP_Phylo.WF_3_5_SNP_Phylo import run_WF_3_5
from WF_4_CreateReport.WF_4_helper import run_create_PDF
from WF_5_extra.WF_5_upload_ncbi import run_WF_5
import os
import sys
import json
import reader

#Main body of the scrip which will run the rest


class CRAB_pipeline_worker():

    def __init__(self, cache_path) :
        self.cache_path = cache_path
        #self.path_to_reads = path_to_reads
        #self.sample_sheet_p = sample_sheet_p
        #self.run_data = sample_sheet_p.split("/")[-1][:-4]

        demo_cahce= reader.read_json(cache_path+"/data/pipeline_variables.json")

        for item in [*demo_cahce] :
            setattr(self,item, demo_cahce[item])

    def run_pipeline(self,path_to_reads,run_date):
        
        #WF_0
        #Fastq pre proccessing, runs SPADES assembler, RETURNS list of HSN
        sample_HSN , Assembly_stats = run_assembly(self.cache_path,path_to_reads,self.assembly_output,self.busco_output,run_date)        

        print("Assembly Done")

        #WF_1
        #runs Prokka
        #runs MLST typing, RETURNS MLST TYPE in DICT {"HSH":[species,type, something, ...]} 
                                                    #{'2296669_manualy': ['2296669_manualy', 'abaumannii_2', '2']}
        self.assembly_output+="/"+run_date  
        self.prokka_output+="/"+run_date      
        mlst = run_annotate(self.assembly_output,self.prokka_output,sample_HSN)
        print("Annotation Done")
        #print(mlst)

        #WF_2
        #Runs Abricate, converts the output to something to be pushed to DB
        self.abricate_output+="/"+run_date 
        found_genes = find_AMR_genes(sample_HSN,self.assembly_output,self.abricate_output)
        #print(found_genes)
        print("found AMR genes")
   
                #found_genes DICT {HSN:[GENE,%COV,%IDENT,DB_Used,Accession_Seq,Gene_Product,Resistance]}

        #WF_3 DB push
        #demographical push
        #gene and anti-micorable data
        #MLST typing 
        run_DB_push(self.cache_path,sample_HSN,mlst,found_genes,Assembly_stats,run_date,self.csv_path) 
        print("Push data to DB")
        
        #3.5 workflow to pull contigs into assembled genome
        #then do snp stuff 
        #and phylogenetic things
        run_WF_3_5(path_to_reads,sample_HSN, self.path_to_shuffled_reads,run_date,self.path_to_referance_genome, self.path_to_snp_output )
        print("Sequences Aligned")
        
        #WF_4 report generation
        #Phylogentics Tree of all samples on run
        #SNP heat map of all samples
        #bring together all informationself.path_to_pdf_output
        run_create_PDF(sample_HSN,run_date, self.path_to_pdf_output ,self.cache_path,found_genes, mlst,self.path_to_snp_output)
        print("Report Generated!")

        
        run_WF_5(sample_HSN,run_date,self.assembly_output,self.path_to_pdf_output,self.StrainDB)



if __name__ == "__main__":
    
    dir_path = "/".join(os.path.dirname(os.path.realpath(__file__)).split("/")[:-1]) #path minus scripts 
            #TO DO
        #Add assembler stats into WF_
        #also create clean up function to delete all random grabo data
    print(sys.argv)
    input_path = sys.argv[1]
    rundate = sys.argv[2]
    print(input_path)
    print("-----------")
    print(rundate)
    print("-----------")
 
    CRAB_p = CRAB_pipeline_worker(dir_path)
       
    CRAB_p.run_pipeline(input_path,rundate)

   # CRAB_pipeline("/Users/adrian/Desktop/CRAB_DATA/062422","/Users/adrian/Desktop/CRAB_DATA/062422_samplesheet.csv")