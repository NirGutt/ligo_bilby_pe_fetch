# ligo_bilby_pe_fetch

A lightweight Python module to fetch and process Bilby parameter estimation (PE) results, particularly useful for LIGO/Virgo data analysis workflows.

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
##  Available Functions

- `RPE.get_interferometers()`  
  → Returns the interferometer data used in the analysis.

- `RPE.get_waveform()`  
  → Retrieves the waveform generator used for the PE.

- `RPE.get_likelihood()`  
  → Returns the likelihood function for the analysis.

- `RPE.get_diagnostics_result()`  
  → Returns diagnostics for the run (e.g. comparison between prior and posterior).
