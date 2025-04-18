import pickle
import h5py
import numpy as np
import subprocess
import os
import glob
import bilby



class ligo_bilby_pe_fetch:
    def __init__(self, file_name,
                 waveform_name=None,
                 override_previous_run=False, 
                 user='albert.einstein',
                 psd_files_dict=None,
                 shift_trigger_time=None,
                 run_diagnoostics_samples=10):
        """
        Initialize the ReproducePE object.

        Parameters
        ----------
        file_name : str
            The name (or path) to the file that contains the data (e.g., a pickle file).
        waveform_name : str
            The name of the waveform to retrieve from the file.
        """
        self.file_name = file_name
        base = os.path.basename(file_name)   # returns 'file_name.ext'
        name, ext = os.path.splitext(base) 
        self.label=name
        self.ini_file_name = file_name+'_config.ini'
        self.data =None        
        self.waveform_name = waveform_name
        self.run_diagnoostics_samples=run_diagnoostics_samples
        self.likelihood_priors=None
        self.pre_wavefprm_text=''
        self.override_previous_run=override_previous_run
        self.psd_files_dict=psd_files_dict
        self.shift_trigger_time=shift_trigger_time
        self.user=user
        self.diag_res=None
        if self.shift_trigger_time is not None:
            print('skipping diagnostics a shifted trigger detected')
            self.run_diagnoostics_samples=0
        
        self._process_PE_file()
        
        self._run_bilby_pipe()
        self._run_data_generation_step()
        
        self._load_data() # load data into self.data

        if self.run_diagnoostics_samples>0:
            self._run_diagnostics()

    def get_diagnostics_result(self):
        return self.diag_res

    def _run_diagnostics(self):
        import pandas as pd
        import random
        import matplotlib.pyplot as plt
        import tqdm

        test_likelihood,likelihood_priors = self.get_likelihood(calibration_marginalization_in=False,
        calibration_number_of_response_curves=10)
        
        # get posterior samples 
        print('Running diagnostics',flush=True)
        f = h5py.File(self.file_name, 'r')
        
        PE_posterior =  f[f"{self.pre_wavefprm_text}{self.waveform_name}"]['posterior_samples']
        random_indices = np.random.choice(len(PE_posterior), size=int(self.run_diagnoostics_samples*2), replace=False)
        random_indices.sort()
        df = pd.DataFrame.from_records(PE_posterior[random_indices],columns=PE_posterior.dtype.names)
        df = df.reset_index(drop=True)
  
        print(f"retriving posterior samples with these columns: {df.columns}")
        dict_samples = [{key: sample[key] for key in df}
                        for _, sample in df.iterrows()]
    
        sample_nums = random.sample(range(0, len(dict_samples)), self.run_diagnoostics_samples)
        ll_original=[]
        ll_test=[]
        for k in tqdm.tqdm(sample_nums,desc='running diagnostics'):
            p=dict_samples[k]
            ll_original.append(p['log_likelihood'])
        
            if test_likelihood.distance_marginalization:            
                p["luminosity_distance"] = likelihood_priors["luminosity_distance"]
            
            if test_likelihood.phase_marginalization:
                p['phase']=likelihood_priors['phase']
        
            if test_likelihood.time_marginalization:
                p["geocent_time"] = likelihood_priors["geocent_time"]
                    
            test_likelihood.parameters.update(p)
            
            
            llr = test_likelihood.log_likelihood_ratio()
            ll_test.append(llr)

        ll_test = np.array(ll_test)
        ll_original = np.array(ll_original)
        plt.figure()
        plt.plot(ll_original,ll_test,'o')
        plt.plot(ll_original,ll_original,'-r')
        plt.ylabel('reconstructed log  likelihood')
        plt.xlabel('PE original likelihood')
        plt.grid(True)
        # save the output
        output_folder='PE_diagnostics'
        os.makedirs(output_folder, exist_ok=True)
        plt.savefig(f"{output_folder}/ReproducePE_{self.label}_diagnostics.png")    
        text='Passed'
        if np.mean(np.abs(ll_test-ll_original))> 0.1:
            text='Failed'
        print(f"diagnostics {text}: mean dlogL = {np.mean(np.abs(ll_test-ll_original))}, the max dlogL = {round(np.max(np.abs(ll_test-ll_original)),2)}.  created ReproducePE_diagnostics.png")
        if text=='Passed':
            self.diag_res= True 
        else:
            self.diag_res= False        



    def _process_PE_file(self):
        f = h5py.File(self.file_name, 'r')
        content= list(f) 
        print(f"This file contains the following: {content}")
        
        if self.waveform_name is None:
            tmp_name = content[0]
            if 'C01:' in tmp_name:
                self.waveform_name=tmp_name[4:]
            else:
                self.waveform_name=tmp_name    
            print(f"No waveform provided, defaulting to {self.waveform_name}")

        if f"C01:{self.waveform_name}" in content:
            self.pre_wavefprm_text='C01:'

        
        
        data_struct = f[f"{self.pre_wavefprm_text}{self.waveform_name}"]['config_file']['config']
             
        with open(self.ini_file_name, "w") as f:
            for key, value in data_struct.items():
                if key=='webdir' or key=='prior-file':
                    continue
                if key=='label':
                    f.write(f"label= {self.label}\n") 
                    continue

                if key=='outdir':
                    f.write(f"outdir = {self.label}_outdir\n")       
                    continue

                # Check if the value is an HDF5 dataset by testing for a 'shape' attribute.
                if hasattr(value, "shape"):
                    # Extract the stored value
                    val_extracted = value[()]  
                    # For single-element datasets, we might want to collapse the array:
                    if isinstance(val_extracted, np.ndarray) and val_extracted.size == 1:
                        val_extracted = val_extracted.item()
                    # If it's a byte string, decode it.
                    if isinstance(val_extracted, bytes):
                        val_extracted = val_extracted.decode('utf-8')
                    value_to_write = val_extracted
                else:
                    value_to_write = value
                
                if key=='trigger-time' and self.shift_trigger_time is not None:
                    value_to_write=+self.shift_trigger_time

                if key=='psd-dict': 
                    if self.psd_files_dict is not None:    
                        value_to_write=self.psd_files_dict
                    else:
                        value_to_write=self._create_local_psd_files()    
                    f.write("psd-dict = {")
                    for key, path in value_to_write.items():
                        f.write(f"{key}: {path},")
                    f.write("}\n")
                    continue

                
                f.write(f"{key} = {value_to_write}\n")
            #make sure teh accounting function is there 
            if 'accounting' not in data_struct.keys():
                 f.write(f"accounting = ligo.prod.o4.cbc.pe.bilby\n")
            
            if 'accounting_user' not in data_struct.keys():
                f.write(f"accounting_user = {self.user}\n")
                      
        print(f"creating {self.ini_file_name}")

    def _create_local_psd_files(self):
        f = h5py.File(self.file_name, 'r')
        psd_data = f[f"{self.pre_wavefprm_text}{self.waveform_name}"]['psds']

        ret_dict={}
    
        for ifo_name in list(psd_data):
            frq = psd_data[ifo_name][:,0]
            vals = psd_data[ifo_name][:,1]
            ret_dict[ifo_name] = np.vstack(list([frq,vals])).T

        output_folder = "psd_files"
        os.makedirs(output_folder, exist_ok=True)
         
        file_path_dict = {}

        # Loop through and save into files
        for key, array in ret_dict.items():
            filename = f"psd_{self.label}_{key}.dat"
            filepath = os.path.join(output_folder, filename)
            
            np.savetxt(filepath, array)  # Save as plain text
            file_path_dict[key] = filepath
        return file_path_dict


    def _run_bilby_pipe(self):

        if os.path.exists(f"{self.label}_outdir/submit/bash_{self.label}.sh") and self.override_previous_run==False:
            print('found a copy of bilby pipe run already, skipping ')
            return

        print('runnign ini bilby pipe')
        result = subprocess.run(f"bilby_pipe {self.ini_file_name}", shell=True,capture_output=True, text=True)
        if result.returncode != 0:
            print("bilby_pipe failed!",flush=True)
            print("--- STDOUT ---")
            print(result.stdout)
            print("--- STDERR ---")
            print(result.stderr)
        else:
            print("bilby_pipe ran successfully.",flush=True)
            print(result.stdout)
        
    def _run_data_generation_step(self):
        pickle_files = glob.glob(f"{self.label}_outdir/data/*.pickle")        
        if pickle_files and self.override_previous_run==False:
            print('found a copy of bilby pipe data generation step already, skipping ')
            return


        print('runnign data generation bilby pipe')
        result = subprocess.run(f"bash {self.label}_outdir/submit/bash_{self.label}.sh generation",shell=True, capture_output=True, text=True)
        if result.returncode != 0:
            print("bilby_pipe failed!",flush=True)
            print("--- STDOUT ---")
            print(result.stdout)
            print("--- STDERR ---")
            print(result.stderr)
        else:
            print("bilby_pipe ran successfully.",flush=True)
            print(result.stdout)

    
    def _load_data(self):

        pattern = f"{self.label}_outdir/data/{self.label}*.pickle"
        matching_files = glob.glob(pattern)
        if matching_files:
            # Return the first matching file (or modify as needed to handle multiple files)
            self.data_bilby_pipe_file_name =  matching_files[0]
        else:
            print('couldnt locate the data file ') 
            return   

       
        try:
            with open(self.data_bilby_pipe_file_name, 'rb') as f:
                self.data = pickle.load(f)
                
        except Exception as e:
            print(f"Error loading file '{self.data_bilby_pipe_file_name}': {e}")
            

    def get_interformeters(self):
       
        if self.data is not None:
            ifos = self.data.interferometers
            return ifos    
        else:
            print('get_interformeters::couldnt find data')

    def get_priors(self):


        f = h5py.File(self.file_name, 'r')
        
        priors_data= f[f"{self.pre_wavefprm_text}{self.waveform_name}"]['priors']['analytic']
        p_dict={}   
        for p in list(priors_data):    
            #if p== 'chirp_mass' or p== 'mass_ratio' or p=='zenith' or p=='azimuth': # skip a probelm 
            #    continue   

            if p== 'chirp_mass':
                cl = priors_data['chirp_mass'][0].decode("utf-8").split("(")
                if cl[0] == "UniformInComponentsChirpMass":
                    complete_cl = "bilby.gw.prior.UniformInComponentsChirpMass("
                    cl[0] = complete_cl
                    string = ''.join(cl)
                    p_dict[p]=string
                    continue
                    
            if p== 'mass_ratio': 
                cl = priors_data['mass_ratio'][0].decode("utf-8").split("(")
                if cl[0] == "UniformInComponentsMassRatio":
                    complete_cl = "bilby.gw.prior.UniformInComponentsMassRatio("
                    cl[0] = complete_cl
                    string = ''.join(cl)
                    p_dict[p]=string
                    continue

            p_dict[p]=priors_data[p][0].decode("utf-8")
        return bilby.core.prior.dict.PriorDict(p_dict) 
        
        
    

    def get_likelihood(self,calibration_marginalization_in=None,calibration_number_of_response_curves=None):
        from bilby.gw.likelihood.base import GravitationalWaveTransient

        if self.data is None:  
            print('get_likelihood::couldnt find data') 
            return None
        meta_data = self.data.meta_data['command_line_args'] 
        
        
        
        if meta_data['likelihood_type']  != 'GravitationalWaveTransient':
            print('get_likelihood:: likelihood-type is not GravitationalWaveTransient, return None') 
            return None 

        if any("recalib" in key for key in self.get_priors().keys()) and meta_data['calibration_marginalization']==False:
            # There is at least one key that contains 'recalib'
           print(f"get_likelihood:: WARNING note the calibration marginlization is set to {meta_data['calibration_marginalization']} but cliabration priors were found , condier setting it to True  ")
        
        calibration_marginalization = meta_data['calibration_marginalization']
        if calibration_marginalization_in is not None:
            calibration_marginalization = calibration_marginalization_in
            print(f"get_likelihood:: overrideing calibration_marginalization found ini settings {meta_data['calibration_marginalization']} and setting it to {calibration_marginalization}")
        

        # calibration settings in case it is acive 
        number_of_response_curves=meta_data['number_of_response_curves']
        if calibration_number_of_response_curves is not None:
            number_of_response_curves=calibration_number_of_response_curves 
           
        # the priors get changes inside of the likelihood, 
        # we want to keep a copy that will get overriden everytime we create the likelihood    
        self.likelihood_priors = self.get_priors()   
        
                  
        likelihood = GravitationalWaveTransient(
            
            interferometers=self.get_interformeters(),
            waveform_generator=self.get_waveform(),
            priors=self.likelihood_priors,

            reference_frame=meta_data['reference_frame'],
            time_reference=meta_data['time_reference'],     
            
            distance_marginalization=meta_data['distance_marginalization'],
            distance_marginalization_lookup_table=meta_data['distance_marginalization_lookup_table'],

            phase_marginalization=meta_data['phase_marginalization'],
            
            time_marginalization=meta_data['time_marginalization'],
            jitter_time=meta_data['jitter_time'],
            
            calibration_marginalization=calibration_marginalization,
            calibration_lookup_table = meta_data['calibration_lookup_table'],
            number_of_response_curves=number_of_response_curves,
            
            )

       
        return likelihood,self.likelihood_priors
        

    def get_waveform(self):
        import ast
        import re

        if self.data is not None:
            metadata = self.data.meta_data['command_line_args'] 
            #metadata = f[f"C01:{self.waveform_name}"]['meta_data']['meta_data']
        else:
            print('get_waveform::couldnt find data')
            return None    
        
        raw_string = metadata['minimum_frequency']
        # Add quotes around keys
        fixed_string = re.sub(r'(\w+)\s*:', r'"\1":', raw_string)
        data_dict = ast.literal_eval(fixed_string)
        minimum_frequency = float(min(data_dict.values()))

        raw_string = metadata['maximum_frequency']
        # Add quotes around keys
        fixed_string = re.sub(r'(\w+)\s*:', r'"\1":', raw_string)
        data_dict = ast.literal_eval(fixed_string)
        maximum_frequency = float(max(data_dict.values()))

        print(f"get_waveform:: metadata[minimum_frequency] = {metadata['minimum_frequency']}, ended up with {minimum_frequency}")
        print(f"get_waveform:: metadata[maximum_frequency] = {metadata['maximum_frequency']}, ended up with {maximum_frequency}")

        
        waveform_arguments = {
            'reference_frequency': float(metadata['reference_frequency']),
            'waveform_approximant': metadata['waveform_approximant'],
            'minimum_frequency': minimum_frequency,
            'maximum_frequency': maximum_frequency,
            'catch_waveform_errors': metadata['catch_waveform_errors'],
            'pn_spin_order': int(metadata['pn_spin_order']),
            'pn_tidal_order': int(metadata['pn_tidal_order']),
            'pn_phase_order': int(metadata['pn_phase_order']),
            'pn_amplitude_order': int(metadata['pn_amplitude_order']),            
        }
        print(type(metadata['waveform_arguments_dict']))
      
        if 'waveform_arguments_dict' in metadata.keys():
            print('get_waveform::found waveform arguments')
            if metadata['waveform_arguments_dict'] not in [None, 'None', '', {}]: 
                print('get_waveform::updating waveform arguments')
                waveform_arguments.update(ast.literal_eval(metadata['waveform_arguments_dict']))

            
        if 'lal_binary_black_hole' != metadata['frequency_domain_source_model']:
            print(f"get_waveform:: frequency_domain_source_model = {metadata['frequency_domain_source_model']}, \
                  currently supporting only lal_binary_black_hole ")
            return None
        
    
        print('get_waveform:: WARNING (not suitable for BNS/NSBH) using bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters')

        WaveformGenerator_class = bilby.gw.WaveformGenerator      
        if metadata['waveform_generator']=='bilby.gw.waveform_generator.LALCBCWaveformGenerator':
            WaveformGenerator_class =bilby.gw.waveform_generator.LALCBCWaveformGenerator
            print('workign with LALCBCWaveformGenerator')

        waveform_generator = WaveformGenerator_class(
            duration=metadata['duration'],
            frequency_domain_source_model=bilby.gw.source.lal_binary_black_hole,
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_black_hole_parameters,
            waveform_arguments=waveform_arguments,      
            sampling_frequency=metadata['sampling_frequency'],        
        )    

        return waveform_generator



    def __repr__(self):
        return f"ReproducePE(file_name='{self.file_name}', waveform_name='{self.waveform_name}')"


    



