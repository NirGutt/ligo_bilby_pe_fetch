import ligo_bilby_pe_fetch 


posterior_file='IGWN-GWTC2p1-v2-GW190929_012149_PEDataRelease_mixed_cosmo.h5'
RPE = ligo_bilby_pe_fetch.ligo_bilby_pe_fetch(posterior_file
                                              ,'IMRPhenomXPHM'
                                              ,run_diagnoostics_samples=100)

# to show what this can get back, # check the code on how to construct teh likelihood in case of marginlizations 
RPE.get_interformeters()
RPE.get_waveform()
RPE.get_likelihood()

print(f"Diagnositics results {RPE.get_diagnostics_result()}, True means good. The mean abs logL difference < 0.1")
