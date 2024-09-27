
papermill notebooks/AnoPrimer-long.ipynb qPCR_run.ipynb -k AnoPrimer -f tests/cDNA_Params_fun.json &&
papermill notebooks/AnoPrimer-long.ipynb qPCR2_run.ipynb -k AnoPrimer -f tests/cDNA_Params_2_fun.json &&
papermill notebooks/AnoPrimer-long.ipynb gDNA_run.ipynb -k AnoPrimer -f tests/gDNA_probe_Params_fun.json &&
papermill notebooks/AnoPrimer-long.ipynb qPCR_run.ipynb -k AnoPrimer -f tests/cDNA_Params.json &&
papermill notebooks/AnoPrimer-long.ipynb qPCR2_run.ipynb -k AnoPrimer -f tests/cDNA_Params_2.json &&
papermill notebooks/AnoPrimer-long.ipynb gDNA_run.ipynb -k AnoPrimer -f tests/gDNA_probe_Params.json &&
papermill notebooks/AnoPrimer-long.ipynb probe_run.ipynb -k AnoPrimer -f tests/probe_Params.json &&
papermill notebooks/AnoPrimer-short.ipynb short_run.ipynb -k AnoPrimer
