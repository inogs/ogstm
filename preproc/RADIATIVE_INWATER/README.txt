elena terzic @ 04.10.2019
comment to read_nc.py

before running the code save a version of bit.sea/instruments/matchup_manager.py

with a changed line 215
instead of ref_varname = self.reference_var(p, model_varname)
put ref_varname = 'IRR_380'
link this script to that version of bit.sea through PYTHONPATH=...
