from imports import *

class Ramachandran_Analysis():
    
    def __init__(self, AtomData=None):
        """
        AtomData is a pandas dataframe having the format (columns): atom_name, residue_name, x, y, z
        """
        self.AtomData = AtomData
        self.AngleData = []
        self.ANGLE_CSV_PATH = f'./data/Angles_{int(time.time())}.csv'
    
    @log_method
    def getAngles(self, save_to_csv=False): # TODO: split all chains and calculate
        """
        Use a sliding window to calculate the (phi, psi) angles, 
        and return (update) the self.AngleData (list of dict).
        
        - Phi Angle = arccos(n(C_i-1, N_i, CA_i)*n(N_i, CA_i, C_i))
        - Psi Angle = arccos(n(N_i, CA_i, C_i)*n(CA_i, C_i, N_i+1))
        
        Note that the sign of the phi angle is the sign of the inner product of 
            n(C_i-1, N_i, CA_i) and (C_i-CA_i), sign of psi in a similar fashion.
        """
        if self.AtomData is None:
            print("Atom Data is not imported! Please new a instance and initialize with one!")
            return
        
        # Find indices of 'CA' (discard the first and last CA since they have unexist angles)
        ca_indices = self.AtomData[self.AtomData['atom_name'] == ' CA '].index[1:-1]        
        self.AngleData.clear()
        
        for ca_idx in ca_indices:
            
            # get the plane norms and some vectors
            planeNorm_C_N_CA = normVecCross(
                                    self.AtomData.iloc[ca_idx-2, 2:5],    # C_i-1
                                    self.AtomData.iloc[ca_idx-1, 2:5],    # N_i
                                    self.AtomData.iloc[ca_idx, 2:5]       # CA_i
                                )
            
            planeNorm_N_CA_C = normVecCross(
                                    self.AtomData.iloc[ca_idx-1, 2:5],    # N_i
                                    self.AtomData.iloc[ca_idx, 2:5],      # CA_i
                                    self.AtomData.iloc[ca_idx+1, 2:5]     # C_i
                                )
            
            planeNorm_CA_C_N = normVecCross(
                                    self.AtomData.iloc[ca_idx, 2:5],      # CA_i
                                    self.AtomData.iloc[ca_idx+1, 2:5],    # C_i
                                    self.AtomData.iloc[ca_idx+2, 2:5]     # N_i+1
                                )
            
            vec_CA_C         =  (self.AtomData.iloc[ca_idx+1, 2:5].astype(float)       # C_i
                                     - self.AtomData.iloc[ca_idx, 2:5].astype(float)   # CA_i
                                ).to_numpy(dtype=float) 
            
            vec_C_N          = (self.AtomData.iloc[ca_idx+2, 2:5].astype(float)        # N_i+1
                                    - self.AtomData.iloc[ca_idx+1, 2:5].astype(float)  # C_i
                                ).to_numpy(dtype=float)
            
            # get the sign and cos value of phi angle
            phi_sign = np.sign(np.dot(planeNorm_C_N_CA, vec_CA_C))
            phi_cos = np.dot(planeNorm_C_N_CA, planeNorm_N_CA_C)
            
            # get the sign and cos value of psi angle
            psi_sign = np.sign(np.dot(planeNorm_N_CA_C, vec_C_N))
            psi_cos = np.dot(planeNorm_N_CA_C, planeNorm_CA_C_N)
            
            # get the signed (phi, psi) angles
            phi = np.degrees(phi_sign * np.arccos(phi_cos))
            psi = np.degrees(psi_sign * np.arccos(psi_cos))
            
            self.AngleData.append({
                            'residue_name' : self.AtomData.loc[ca_idx, 'residue_name'],
                            'phi'          : phi,
                            'psi'          : psi
                        })
        
        
        self.angle_df = pd.DataFrame.from_dict(self.AngleData)
        
        if save_to_csv:
            self.angle_df.to_csv(self.ANGLE_CSV_PATH)
            print(f'Angle Data is saved at: {self.ANGLE_CSV_PATH}.')
            
        return
    
    def subplot(self, angles_df, df_idx):
        # Create scatter plot
        sns.set(style="whitegrid")
        plt.figure(figsize=(8, 6))
        sns.scatterplot(x=angles_df['phi'],
                        y=angles_df['psi'],
                        alpha=0.6)

        # Set labels and title
        title = ['All Residues but Glycines and Prolines', 'Glycines', 'Prolines']
        
        plt.xlabel('Phi (degrees)')
        plt.ylabel('Psi (degrees)')
        plt.title(f'Ramachandran Plot (Scatter) - {title[df_idx]}')
        
        # Set x and y axis limits
        plt.xlim(-180, 180)
        plt.ylim(-180, 180)
        
        # Set x and y axis ticks
        plt.xticks(np.arange(-180, 181, 45))
        plt.yticks(np.arange(-180, 181, 45))
        
        # Draw x and y axes as dotted lines
        plt.axhline(0, color='black', linestyle='--', lw=2)  # Horizontal line for x axis
        plt.axvline(0, color='black', linestyle='--', lw=2)  # Vertical line for y axis

        # Show plot
        plt.show()
        
        
    def plot(self, angle_filepath=None):
        """
        Generate Ramachandran plots for:
         (a) all residues but glycines and prolines
         (b) all glycines
         (c) all prolines
         
        Note that the Angle DataFrame should have columns: residue_name, phi, psi
        """
        
        if angle_filepath is not None:
            angle_data = pd.read_csv(angle_filepath, index_col=0)
        else:
            try:
                angle_data = self.angle_df
                # will raise error if not executed self.getAngles()
            except:
                print('There is no Angle Data to be plot. \
                        It seems you have not executed self.getAngles()!')
            
        angle_rest = angle_data.loc[~angle_data['residue_name'].isin(('GLY', 'PRO'))]
        angle_gly  = angle_data.loc[angle_data['residue_name'] == 'GLY']
        angle_pro  = angle_data.loc[angle_data['residue_name'] == 'PRO']
        
        for df_idx, angles_df in enumerate([angle_rest, angle_gly, angle_pro], start=0):
            self.subplot(angles_df, df_idx)
            
    def main(self):
        self.getAngles(save_to_csv=True)
        self.plot()
        
