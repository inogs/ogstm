#! /bin/bash
useintel
rm *.nc *.tmp *.png
./data_parser.sh TEST02
./crea_file.sh
pyload 
python plot_hovmoeller_full.py
