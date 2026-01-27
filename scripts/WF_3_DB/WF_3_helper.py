import sys
sys.path.insert(0,'/epi/home/ashita.jawali@kdhe.state.ks.us/Documents/GitHub/CRAB_Analysis/scripts')
from ms_sql_handler import ms_sql_handler
import pandas as pd
import cx_Oracle as co
import reader 
import datetime
from other import add_cols
import json



class demographics_import():

    def __init__(self,cache_path) : #0
        #here need to import json file
        #and used that to store
       
        #demo_cahce= reader.read_json(cache_path+"/data/demographics.json")
        demo_cahce= json.load(open(cache_path+"/data/demographics.json"))
        for item in [*demo_cahce] :
            setattr(self,item, demo_cahce[item])
        

        #and import metric data needed
#        self.df_hsn = pd.read_json(cache_path+"/data/sample_metrics.json")
        #df2['year']=df2['year'].astype(int)
#        for df_column in self.df_hsn.columns:
#            self.df_hsn[df_column]=self.df_hsn[df_column].astype(int)

    
    def get_lims_demographics(self,hsn,date,csv_path): #1
        self.wgs_run_date = date[:2]+"/"+date[2:4]+"/20"+date[4:]
        unfound_hsn=[]
        conn = co.connect(self.lims_connection)      
        query="select * from crecrpa5demo where HSN in (" + ",".join(hsn) + ")"
        self.lims_df = pd.read_sql(query, conn)
        excel_df =pd.read_excel(csv_path + "/HAI_Metadata.xlsx", sheet_name='Sheet1', converters={'HSN':int})
        found_hsn = [str(i) for i in self.lims_df['HSN'].values.tolist() ]
        self.no_lims_hsn = pd.DataFrame(columns=excel_df.columns.values.tolist())
        for h in hsn: #if no demographical information add a blank line
            if h not in found_hsn :
                print(str(h)+" not found in lims df")
                try:
                    r= excel_df.query("HSN == " + str(h))
                    self.no_lims_hsn = pd.concat([r,self.no_lims_hsn])
                except Exception as e:
                    print("This error was found: \t"+e)
                    print("-"*10)
                    print(str(h)+" not found in CSV file")
                    #hsn.remove(h)
                    unfound_hsn.append(h) 
        conn.close()
        print("Not found in LIMS or CSV")
        print(unfound_hsn)
        return hsn

    def format_lims_df(self): #2
        # manipulate sql database to format accepted by the master EXCEL worksheet
        #self.log.write_log("format_lims_DF","Manipulating demographics to database format")

        self.lims_df = self.lims_df.rename(columns = self.demo_names)
        #print(self.lims_df.head().to_string())
        def safe_to_datetime(date_str):
            try:
                return pd.to_datetime(date_str)
            except ValueError:
                return pd.NaT

        self.no_lims_hsn = self.no_lims_hsn.loc[:, ~self.no_lims_hsn.columns.str.contains('^Unnamed')]
        self.no_lims_hsn['name']=self.no_lims_hsn['First Name'] + " " + self.no_lims_hsn['Last Name']
        self.no_lims_hsn['HSN']=self.no_lims_hsn['HSN'].astype(int)
        self.no_lims_hsn = self.no_lims_hsn.rename(columns={"Rec'd":"date_recd","WGS Rec'd": "pcr_run_date","HSN": "hsn","Collected": "doc","DOB": "dob","Sex": "sex","State": "state","Source Site": "source","CO":"county"})
        self.no_lims_hsn['pcr_run_date']= pd.to_datetime(self.no_lims_hsn['pcr_run_date'])
        self.no_lims_hsn['date_recd']= pd.to_datetime(self.no_lims_hsn['date_recd'])
        self.no_lims_hsn['dob']= self.no_lims_hsn['dob'].apply(safe_to_datetime)
        self.no_lims_hsn['doc']= pd.to_datetime(self.no_lims_hsn['doc'])
        self.no_lims_hsn.drop(columns=['Extracted','Sequenced','Last Name','First Name','Age','Source Type','Country','Comment','WGS serotype','coverage (calculated from workbook)','#total reads','Clusters passing filter',"HAI WGS ID"], axis=1, inplace=True)

        # if not self.no_lims_hsn.empty:
        self.lims_df = pd.concat([self.lims_df, self.no_lims_hsn])
        #self.log.write_log("format_lims_DF","Done!")


    def create_mlst_df(self,mlst_dict):
    #mlst  DICT {HSN:[HSN,species,overallType]}
        #{'2296669_manualy': ['2296669_manualy', 'abaumannii_2', '2']}

        self.mlst_df = pd.DataFrame.from_dict(mlst_dict,orient='index',columns=['hsn','species','mlst'])
        self.mlst_df['hsn']=self.mlst_df['hsn'].astype(int)
        # print(len(self.mlst_df))
        # print(self.mlst_df.head())

    def create_metrics_df(self,assembly_m,pipeline):
        #res[sample]=sample_buso_res["results"]["one_line_summary"]
        #assemblys  DICT {HSN:"C:98.4%[S:98.4%,D:0.0%],F:1.6%,M:0.0%,n:124"}
        # meaning of output "Complete": 98.4, "Single copy": 98.4,"Multi copy": 0.0,"Fragmented": 1.6, "Missing": 0.0, "n_markers": 124,
        for h in [*assembly_m]:
            if pipeline:
                assembly_m[h] =  assembly_m[h]
            else:
                assembly_m[h] =  assembly_m[h][2:6]
            
        self.metrics_df = pd.DataFrame.from_dict(assembly_m,orient='index',columns=['busco_out'])
        self.metrics_df['hsn']=self.metrics_df.index.astype(int)
        #print(len(self.metrics_df))
        #print(self.metrics_df.head())
    
    def create_genes_df(self,found_genes_dict):
        #found_genes DICT {HSN:[[GENE,%COV,%IDENT,DB_Used,Accession_Seq,Gene_Product,Resistance],[GENE2....]]}

        self.genes_df = pd.DataFrame.from_dict(found_genes_dict, orient='index',columns=['hsn','gene','coverage','identity','db_used','accession_seq','gene_product','resistance'])

        self.genes_df['hsn']=self.genes_df['hsn'].astype(int)
        #print(len(self.genes_df))
        #print(self.genes_df.head())

        

    def merge_dfs(self): #3
        #will use this to merge all the different DFs 
        #Demographics
        #Run Stats
        #MLST typing

        #self.log.write_log("merge_dfs","Merging dataframes")
        self.lims_df['hsn']=self.lims_df['hsn'].astype(int)
        
        self.df = pd.merge(self.lims_df, self.mlst_df, how="inner", on="hsn")

        self.df = pd.merge(self.df, self.metrics_df, how="inner", on="hsn")


       #self.log.write_log("merge_dfs","Done")
    
    def format_dfs(self): #4
        self.df = self.df.rename(columns = self.demo_names)
        #self.log.write_log("format_dfs","Starting")
        # format columns, insert necessary values
        #self.log.write_log("format_dfs","Adding/Formatting/Sorting columns")
        #print(self.df.columns.to_list())
        self.df = add_cols(obj=self, \
            df=self.df, \
            col_lst=self.add_col_lst, \
            col_func_map=self.col_func_map)
        # sort/remove columns to match list
        self.df = self.df[self.sample_data_col_order]

    
    def create_ncbi_csv(self,csv_out_path,rundate,demo_list):
        
        csv_f = open(csv_out_path+"/"+rundate+"_ncbi_pathogen.csv","w+")
        #header
        csv_f.write("sample_name,bioproject_accession,organism,isolate,collected_by,collection_date,geo_loc_name,host,host_disease,isolation_source,lat_lon,MLST\n")
        #YEARDG-0001
        for samp in demo_list:
            #loop through demolist and find the corret piece of info needed
            csv_f.write(samp[24]+',PRJNA288601,Acinetobacter baumannii,'+'ISOLATE'+',NA,'+samp[14]+',USA,Homo sapiens,NA,'+samp[17]+',Missing,MLST'+samp[18]+"_PubMLST\n")                                     
        
        csv_f.close()

    def assign_HAI_ID(self,formatted_df,year,max_id,path_to_excel): #could alos be used to add in run metrics
        
        run_metrics =pd.read_excel(path_to_excel+"/HAI_Metadata.xlsx",sheet_name='Sheet1',converters={'HSN':int})

        for i,row in formatted_df.iterrows():
            max_id+=1
            hsn = str(formatted_df['hsn'].iloc[i])
            res = run_metrics.query("HSN == "+hsn)
            
            formatted_df.at[i,'HAI_WGS_ID'] = year+"DG-"+str(max_id).rjust(5,'0') 

            if not(res.empty):              
                
                formatted_df.at[i,'PhiX174_Recovery'] = res["PhiX recovery"].iloc[0]
                formatted_df.at[i,'Q30'] = res["Q30%"].iloc[0]
                formatted_df.at[i,'Cluster_Density'] = res["Cluster density"].iloc[0]
                formatted_df.at[i,'pass cluster'] = res["Clusters passing filter"].iloc[0]
                formatted_df.at[i,'Total_Num_Reads'] = res["#total reads"].iloc[0]
            else:
                print("HSN "+hsn+" Not found on excel sheet")
            #calculated coverage need to make a field for this
            #seq coverage? (cause the other is assembly coverage)
            #formatted_df.at[i,'coveraaaaaa'] = res["coverage (calculated from workbook)"].iloc[0]

        return formatted_df

    def database_push(self, excel_path,runD): #5
        #self.log.write_log("database_push","Starting")
        self.setup_db()
   
        try:
            #change query to where HAI ID is like current year
            HAI_ID = self.db_handler.sub_read(query="SELECT MAX(HAI_WGS_ID) AS MAX_ID FROM dbo.Results where wgs_run_date >= cast('"+str(self.wgs_run_date[-4:])+"' as datetime)")
            HAI_ID_MAX=int(HAI_ID.to_string()[-5:])
            #2023DG-00080
            #2023DG00080
        except:
            print("No HAI ID")
            HAI_ID_MAX = 0


        self.df = self.assign_HAI_ID(self.df,self.wgs_run_date[-4:],HAI_ID_MAX,excel_path)
        #print(self.df.head())
        df_demo_lst = self.df.values.astype(str).tolist()
        # print(df_demo_lst[:5])  # Check the first few rows of data being pushed to SQL
        i=0
        while i < len(df_demo_lst):
        
            t= df_demo_lst[i][-1].split("%")[0][2:]
            
            df_demo_lst[i]= df_demo_lst[i][:-1]+[t]

            i+=1
    
        #df_table_col_query = "(" + ", ".join(self.df.columns.astype(str).tolist()) + ")"
        
        self.write_query_tbl1 = (" ").join(self.write_query_tbl1)
        # print(self.write_query_tbl1)
        self.create_ncbi_csv(excel_path,runD,df_demo_lst)

        try:
            self.db_handler.lst_ptr_push(df_lst=df_demo_lst, query=self.write_query_tbl1)
        except Exception as e:
            print(f"Error pushing data to SQL: {e}")

        #self.log.write_log("database_push","Done!`")

    def setup_db(self):
        self.db_handler = ms_sql_handler(self)
        self.db_handler.establish_db()
    
    #will need a function for resistance table push
    def database_push_genes(self):
        self.setup_db()
        df_demo_lst = self.genes_df.values.astype(str).tolist()
        self.write_query_genes = (" ").join(self.write_query_genes)
   
        self.db_handler.lst_ptr_push(df_lst=df_demo_lst, query=self.write_query_genes)
    


if __name__ == "__main__":
    
    import_demo = demographics_import("/epi/home/ashita.jawali@kdhe.state.ks.us/Documents/GitHub/CRAB_Analysis")
    sample_hsn = import_demo.get_lims_demographics(['2434975','2445821','2468507','2488768','2492075','2506355','2510743','2527973'],"111323","/epi/home/ashita.jawali@kdhe.state.ks.us/WGS_Drive/CRAB_WGS_Sequencing")
    
        
