#!/usr/bin/env python3
import netCDF4 as nc
import numpy as np
from pathlib import Path

# Cartella di input/output
input_file  = "/leonardo_work/OGS_test2528_0/plazzari/OGSTM-FABM/OGSTM-BFM-qDeg/MED_BIODIV_T0/wrkdir/MODEL/BC/ATL_yyyy0630-00:00:00.nc"
output_file = "/leonardo_work/OGS_test2528_0/plazzari/OGSTM-FABM/OGSTM-BFM-qDeg/MED_BIODIV_T0/wrkdir/MODEL/BC/ATL98_yyyy0630-00:00:00.nc"

# Dizionario di mapping:
# chiave   = nome file originale
# valore   = lista dei nuovi file da creare
mapping = {
    "N1_p": ["N1_p"],
    "N3_n": ["N3_n"],
    "N4_n": ["N4_n"],
    "N5_s": ["N5_s"],
    "N6_r": ["N6_r"],
    "O2_o": ["O2_o"],
    "O3_c": ["O3_c"],
    "O3h_h": ["O3h_h"],
    "O4_n": ["O4_n"],
    "O5_c": ["O5_c"],
    "R1_c": ["R1_c"],
    "R1_n": ["R1_n"],
    "R1_p": ["R1_p"],
    "R1_s": ["R1_s"],
    "R2_c": ["R2_c"],
    "R3_c": ["R3_c"],
    "R6_c": ["R6_c"],
    "R6_n": ["R6_n"],
    "R6_p": ["R6_p"],
    "R6_s": ["R6_s"],
    "R8_c": ["R8_c"],
    "R8_n": ["R8_n"],
    "R8_p": ["R8_p"],
    "R8_s": ["R8_s"],
    "X1_c": ["X1_c"],
    "X2_c": ["X2_c"],
    "X3_c": ["X3_c"],
    "B1_c": ["B1_m5_c"],
    "B1_n": ["B1_m5_n"],
    "B1_p": ["B1_m5_p"],
    "P1_c": ["P1_11_c", "P1_12_c", "P1_13_c", "P1_14_c", "P1_15_c", "P1_16_c", "P1_17_c", "P1_18_c", "P1_19_c", "P1_20_c"],
    "P1_n": ["P1_11_n", "P1_12_n", "P1_13_n", "P1_14_n", "P1_15_n", "P1_16_n", "P1_17_n", "P1_18_n", "P1_19_n", "P1_20_n"],
    "P1_p": ["P1_11_p", "P1_12_p", "P1_13_p", "P1_14_p", "P1_15_p", "P1_16_p", "P1_17_p", "P1_18_p", "P1_19_p", "P1_20_p"],
    "P1_Chl": ["P1_11_Chl", "P1_12_Chl", "P1_13_Chl", "P1_14_Chl", "P1_15_Chl", "P1_16_Chl", "P1_17_Chl", "P1_18_Chl", "P1_19_Chl", "P1_20_Chl"],
    "P1_s": ["P1_11_s", "P1_12_s", "P1_13_s", "P1_14_s", "P1_15_s", "P1_16_s", "P1_17_s", "P1_18_s", "P1_19_s", "P1_20_s"],
    "P4_c": ["P4_11_c", "P4_12_c", "P4_13_c", "P4_14_c", "P4_15_c", "P4_16_c", "P4_17_c", "P4_18_c", "P4_19_c", "P4_20_c"],
    "P4_n": ["P4_11_n", "P4_12_n", "P4_13_n", "P4_14_n", "P4_15_n", "P4_16_n", "P4_17_n", "P4_18_n", "P4_19_n", "P4_20_n"],
    "P4_p": ["P4_11_p", "P4_12_p", "P4_13_p", "P4_14_p", "P4_15_p", "P4_16_p", "P4_17_p", "P4_18_p", "P4_19_p", "P4_20_p"],
    "P4_Chl": ["P4_11_Chl", "P4_12_Chl", "P4_13_Chl", "P4_14_Chl", "P4_15_Chl", "P4_16_Chl", "P4_17_Chl", "P4_18_Chl", "P4_19_Chl", "P4_20_Chl"],
    "P2_c": ["P2_4_c", "P2_5_c", "P2_6_c", "P2_7_c", "P2_8_c", "P2_9_c", "P2_10_c", "P2_11_c", "P2_12_c", "P2_13_c", "P5_6_c", "P5_7_c", "P5_8_c", "P5_9_c", "P5_10_c", "P5_11_c", "P5_12_c", "P5_13_c", "P5_14_c", "P5_15_c", "P7_2_c", "P7_3_c", "P7_4_c", "P7_5_c", "P7_6_c", "P7_7_c", "P8_2_c", "P8_3_c", "P8_4_c", "P8_5_c", "P8_6_c", "P8_7_c"],
    "P2_n": ["P2_4_n", "P2_5_n", "P2_6_n", "P2_7_n", "P2_8_n", "P2_9_n", "P2_10_n", "P2_11_n", "P2_12_n", "P2_13_n", "P5_6_n", "P5_7_n", "P5_8_n", "P5_9_n", "P5_10_n", "P5_11_n", "P5_12_n", "P5_13_n", "P5_14_n", "P5_15_n", "P7_2_n", "P7_3_n", "P7_4_n", "P7_5_n", "P7_6_n", "P7_7_n", "P8_2_n", "P8_3_n", "P8_4_n", "P8_5_n", "P8_6_n", "P8_7_n"],
    "P2_p": ["P2_4_p", "P2_5_p", "P2_6_p", "P2_7_p", "P2_8_p", "P2_9_p", "P2_10_p", "P2_11_p", "P2_12_p", "P2_13_p", "P5_6_p", "P5_7_p", "P5_8_p", "P5_9_p", "P5_10_p", "P5_11_p", "P5_12_p", "P5_13_p", "P5_14_p", "P5_15_p", "P7_2_p", "P7_3_p", "P7_4_p", "P7_5_p", "P7_6_p", "P7_7_p", "P8_2_p", "P8_3_p", "P8_4_p", "P8_5_p", "P8_6_p", "P8_7_p"],
    "P2_Chl": ["P2_4_Chl", "P2_5_Chl", "P2_6_Chl", "P2_7_Chl", "P2_8_Chl", "P2_9_Chl", "P2_10_Chl", "P2_11_Chl", "P2_12_Chl", "P2_13_Chl", "P5_6_Chl", "P5_7_Chl", "P5_8_Chl", "P5_9_Chl", "P5_10_Chl", "P5_11_Chl", "P5_12_Chl", "P5_13_Chl", "P5_14_Chl", "P5_15_Chl", "P7_2_Chl", "P7_3_Chl", "P7_4_Chl", "P7_5_Chl", "P7_6_Chl", "P7_7_Chl", "P8_2_Chl", "P8_3_Chl", "P8_4_Chl", "P8_5_Chl", "P8_6_Chl", "P8_7_Chl"],
    "P3_c": ["P3_0_c", "P3_1_c", "P3_2_c", "P9_m2_c", "P9_m1_c", "P9_0_c", "P6_m4_c", "P6_m3_c"],
    "P3_n": ["P3_0_n", "P3_1_n", "P3_2_n", "P9_m2_n", "P9_m1_n", "P9_0_n", "P6_m4_n", "P6_m3_n"],
    "P3_p": ["P3_0_p", "P3_1_p", "P3_2_p", "P9_m2_p", "P9_m1_p", "P9_0_p", "P6_m4_p", "P6_m3_p"],
    "P3_Chl": ["P3_0_Chl", "P3_1_Chl", "P3_2_Chl", "P9_m2_Chl", "P9_m1_Chl", "P9_0_Chl", "P6_m4_Chl", "P6_m3_Chl"],
    "Z6_c": ["Z6_5_c", "Z6_6_c", "Z6_7_c", "Z6_8_c", "Z6_9_c", "Z6_10_c", "Z6_11_c"],
    "Z6_n": ["Z6_5_n", "Z6_6_n", "Z6_7_n", "Z6_8_n", "Z6_9_n", "Z6_10_n", "Z6_11_n"],
    "Z6_p": ["Z6_5_p", "Z6_6_p", "Z6_7_p", "Z6_8_p", "Z6_9_p", "Z6_10_p", "Z6_11_p"],
    "Z5_c": ["Z5_9_c", "Z5_10_c", "Z5_11_c", "Z5_12_c", "Z5_13_c", "Z5_14_c", "Z5_15_c", "Z5_16_c", "Z5_17_c", "Z5_18_c"],
    "Z5_n": ["Z5_9_n", "Z5_10_n", "Z5_11_n", "Z5_12_n", "Z5_13_n", "Z5_14_n", "Z5_15_n", "Z5_16_n", "Z5_17_n", "Z5_18_n"],
    "Z5_p": ["Z5_9_p", "Z5_10_p", "Z5_11_p", "Z5_12_p", "Z5_13_p", "Z5_14_p", "Z5_15_p", "Z5_16_p", "Z5_17_p", "Z5_18_p"],
    "Z4_c": ["Z4_15_c", "Z4_16_c", "Z4_17_c", "Z4_18_c", "Z4_19_c", "Z4_20_c", "Z4_21_c", "Z4_22_c", "Z4_23_c", "Z4_24_c"],
    "Z4_n": ["Z4_15_n", "Z4_16_n", "Z4_17_n", "Z4_18_n", "Z4_19_n", "Z4_20_n", "Z4_21_n", "Z4_22_n", "Z4_23_n", "Z4_24_n"],
    "Z4_p": ["Z4_15_p", "Z4_16_p", "Z4_17_p", "Z4_18_p", "Z4_19_p", "Z4_20_p", "Z4_21_p", "Z4_22_p", "Z4_23_p", "Z4_24_p"],
    "Z3_c": ["Z3_22_c", "Z3_23_c", "Z3_24_c", "Z3_25_c", "Z3_26_c", "Z3_27_c", "Z3_28_c", "Z3_29_c", "Z3_30_c", "Z3_31_c"],
    "Z3_n": ["Z3_22_n", "Z3_23_n", "Z3_24_n", "Z3_25_n", "Z3_26_n", "Z3_27_n", "Z3_28_n", "Z3_29_n", "Z3_30_n", "Z3_31_n"],
    "Z3_p": ["Z3_22_p", "Z3_23_p", "Z3_24_p", "Z3_25_p", "Z3_26_p", "Z3_27_p", "Z3_28_p", "Z3_29_p", "Z3_30_p", "Z3_31_p"]
}

def main():
    with nc.Dataset(input_file, 'r') as src:
        with nc.Dataset(output_file, 'w', format=src.file_format) as dst:
            dst.setncatts(src.__dict__)
            for name, dimension in src.dimensions.items():
                dst.createDimension(name, (len(dimension) if not dimension.isunlimited() else None))

            for name, variable in src.variables.items():
                if name in mapping:
                    target_names = mapping[name]
                    n_val = len(target_names)
                    
                    print(f"Dividendo {name}")
                    
                    for t_name in target_names:

                        new_var = dst.createVariable(t_name, variable.datatype, variable.dimensions)
                        new_var.setncatts(variable.__dict__)
                        
                        # dividi i nomi per il numero di valori
                        new_var[:] = src[name][:] / n_val
                
                else:
                    if name not in dst.variables:
                        new_var = dst.createVariable(name, variable.datatype, variable.dimensions)
                        new_var.setncatts(variable.__dict__)
                        new_var[:] = src[name][:]

    print(f"Creato {output_file}")

if __name__ == "__main__":
    main()
    
