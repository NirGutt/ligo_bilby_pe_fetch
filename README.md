# ligo_bilby_pe_fetch

A lightweight Python module to fetch and process Bilby parameter estimation (PE) objects, which might be needed for LIGO/Virgo data analysis workflows.

---

## Manual

###  `ligo_bilby_pe_fetch` Inputs:


ligo_bilby_pe_fetch(
    file_name,
    waveform_name=None,
    override_previous_run=False,
    user='albert.einstein',
    psd_files_dict=None,
    shift_trigger_time=None,
    run_diagnostics_samples=10
)

file_name               (str)   – Path to the Bilby PE result file (e.g., result.h5 or directory). Required.

waveform_name           (str)   – Optional waveform name of not given will take the first one in the file.

override_previous_run   (bool)  – If True, reruns setup even if cached/intermediate files exist. Default: False.

user                    (str)   – Username for internal logging or tagging. 'albert.einstein'.

psd_files_dict          (dict)  – Optional dict like {'H1': 'path/to/psd.dat'}. If None, PSDs are taken from the PE result file.

shift_trigger_time      (float) – If provided, shifts the trigger time. Might be useful.. Default: None.

run_diagnostics_samples (int)   – Number of posterior samples used for diagnostics. Default: 10.you just want to make sure you get x=y

##  Available Functions

Once we have the class,
RPE = ligo_bilby_pe_fetch(
    file_name,
    waveform_name=None,
    override_previous_run=False,
    user='albert.einstein',
    psd_files_dict=None,
    shift_trigger_time=None,
    run_diagnostics_samples=10
)

One can use it to call the following functions:

- `RPE.get_interferometers()`  
  → Returns the interferometer data used in the analysis.

- `RPE.get_waveform()`  
  → Retrieves the waveform generator used for the PE.

- `RPE.get_likelihood()`  
  → Returns the likelihood function for the analysis.

- `RPE.get_diagnostics_result()`  
  → Returns the diagnostics result for the PE (i.e. comparison between the fetched likelihood values and the ones stored in the file).
