
papermill notebooks/AnoPrimer-long.ipynb tests/qPCR_run.ipynb -k AnoPrimer -f tests/cDNA_Params_fun.json &&
papermill notebooks/AnoPrimer-long.ipynb tests/qPCR2_run.ipynb -k AnoPrimer -f tests/cDNA_Params_2_fun.json &&
papermill notebooks/AnoPrimer-long.ipynb tests/gDNA_run.ipynb -k AnoPrimer -f tests/gDNA_probe_Params_fun.json &&
papermill notebooks/AnoPrimer-long.ipynb tests/qPCR_run.ipynb -k AnoPrimer -f tests/cDNA_Params.json &&
papermill notebooks/AnoPrimer-long.ipynb tests/qPCR2_run.ipynb -k AnoPrimer -f tests/cDNA_Params_2.json &&
papermill notebooks/AnoPrimer-long.ipynb tests/gDNA_run.ipynb -k AnoPrimer -f tests/gDNA_probe_Params.json &&
papermill notebooks/AnoPrimer-long.ipynb tests/probe_run.ipynb -k AnoPrimer -f tests/probe_Params.json &&
papermill notebooks/AnoPrimer-short.ipynb tests/short_run.ipynb -k AnoPrimer
