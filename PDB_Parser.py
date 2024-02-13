from imports import *


class PDB_Parser():
    def __init__(self, ATOM_DATA_FILE_PATH = f'./AtomData_{int(time.time())}.csv'):
        # Sequence clusters at <identity> % sequence identity clustering
        self.CLUSTER_FILE_URL = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt"
        self.PDB_URL_TEMPLATE = Template("https://files.rcsb.org/download/${pdb_id}.pdb")
        self.FASTA_URL_TEMPLATE = Template("https://www.rcsb.org/fasta/entry/${pdb_id}")
        self.structQueue = []
        self.structQueueSIZE = 125
        self.STRUCT_FORMAT = re.compile('^.{4}_[0-9]$')
        self.CHAIN_FORMAT = re.compile('\|Chains? (.*?)\|')
        self.AUTH_CHAIN_FORMAT = re.compile('\[auth (.+)\]')
        self.BACKBONE_ATOMS = (' N  ',' CA ',' C  ')
        self.atomData = []
        self.ATOM_DATA_FILE_PATH = ATOM_DATA_FILE_PATH
        
    
    @log_method
    def getStruct(self):
        """
        loop over the largest 100 clusters (the first 100 lines), 
        and extract one random structure for each cluster.
        """
        r = requests.get(self.CLUSTER_FILE_URL, stream=True)
        self.structQueue.clear()
        
        for line in r.iter_lines():
            line = line.decode("utf-8")
            if len(self.structQueue) == self.structQueueSIZE:
                break
            if line:
                # randomly sample a structure 
                structs = line.split()             
                struct_id = sampleWithConstraints(structs, self.STRUCT_FORMAT)
                
                if not struct_id:
                    continue
                    
                self.structQueue.append(struct_id) 
              
            
    def getFASTA(self, struct_id):
        """
        A subroutine that fetches the FASTA File, and return the chain_id
            Note: struct_id (instances in structQueue) = {pdb_id}_{pe_id}
        """
        fasta_url = self.FASTA_URL_TEMPLATE.substitute(pdb_id=struct_id[:4])
        r = requests.get(fasta_url, stream=True)
        
        for line in r.iter_lines():
            line = line.decode("utf-8")
            # match the line with the struct_id
            if line[0] != ">":
                continue
                
            if line[1:7] == struct_id:
                chain_ids = self.CHAIN_FORMAT.findall(line)[0].split(", ")
                
                # Select ONE of the chain_id(s)
                while len(chain_ids) != 0:
                    chain_id = chain_ids.pop(-1)
                    auth_chain_id = self.AUTH_CHAIN_FORMAT.search(chain_id)
                    if auth_chain_id:
                        return auth_chain_id.group(1)
                    else:
                        return chain_id
                
        return False
      
        
    def getPDB(self, struct_id, chain_id):
        """
        Get the PDB File, and parse it to a clean format for further usage.
         Return a list of dicts containing atomData (segment of the whole)
            Note: 
                - struct_id (instances in structQueue) = {pdb_id}_{pe_id}
                - chain_id is from FASTA (use method getFASTA())
                - for the line(str) slicing indices, please referr to the PDB file format
        """
        pdb_url = self.PDB_URL_TEMPLATE.substitute(pdb_id=struct_id[:4])
        r = requests.get(pdb_url, stream=True)
        segAtomData = []
        
        for line in r.iter_lines():
            line = line.decode("utf-8")
            if line.startswith('ATOM'):
                if line[21] == chain_id and line[12:16] in self.BACKBONE_ATOMS:
                    segAtomData.append({
                            'atom_name' : line[12:16],
                            'residue_name' : line[17:20],
                            'x' : line[30:38],
                            'y' : line[38:46],
                            'z' : line[46:54]
                        })
                    
            elif len(segAtomData) != 0 and line[21] != chain_id:
                return segAtomData
        
        return False
   
    @log_method
    def getAtomData(self, save_to_csv=False):
        """
        Get the FULL Atom Data for further analysis, go through the subroutine:
            - getFASTA
            - getPDB
         and return a dataframe with columns: atom_name, residue_name, x, y, z.
        """
        cnt = 0
        for struct_id in self.structQueue:
            chain_id = self.getFASTA(struct_id)
            seg_atomData = self.getPDB(struct_id, chain_id)
            
            if seg_atomData:
                # skip the cif format ones
                # TODO: To improve, modify to handle the cif ones
                self.atomData += seg_atomData
                cnt += 1
                
        atom_data = pd.DataFrame.from_dict(self.atomData)
        
        if save_to_csv:
            atom_data.to_csv(self.ATOM_DATA_FILE_PATH)
            print(f'The Atom Data is saved at {self.ATOM_DATA_FILE_PATH}.')
            
        print(f'{cnt} Structures Included.')
        
        return atom_data
        
    
    def main(self, save_to_csv=False):
        self.getStruct()
        
        return self.getAtomData(save_to_csv=save_to_csv)
        
