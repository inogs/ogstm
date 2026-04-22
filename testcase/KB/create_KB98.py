#!/usr/bin/env python3
import os
from pathlib import Path

# Cartella di input/output
input_dir = Path("INIT_NWM_KB_FABM")
output_dir = Path("INIT_NWM_KB_FABM98")
output_dir.mkdir(exist_ok=True)

# Dizionario di mapping:
# chiave   = nome file originale
# valore   = lista dei nuovi file da creare
mapping = {
    "INIT.depth": ["INIT.depth"],
    "INIT.N1_p": ["INIT.N1_p"],
    "INIT.N3_n": ["INIT.N3_n"],
    "INIT.N4_n": ["INIT.N4_n"],
    "INIT.N5_s": ["INIT.N5_s"],
    "INIT.N6_r": ["INIT.N6_r"],
    "INIT.O2_o": ["INIT.O2_o"],
    "INIT.O3_c": ["INIT.O3_c"],
    "INIT.O3h_h": ["INIT.O3h_h"],
    "INIT.O4_n": ["INIT.O4_n"],
    "INIT.O5_c": ["INIT.O5_c"],
    "INIT.R1_c": ["INIT.R1_c"],
    "INIT.R1_n": ["INIT.R1_n"],
    "INIT.R1_p": ["INIT.R1_p"],
    "INIT.R1_s": ["INIT.R1_s"],
    "INIT.R2_c": ["INIT.R2_c"],
    "INIT.R3_c": ["INIT.R3_c"],
    "INIT.R6_c": ["INIT.R6_c"],
    "INIT.R6_n": ["INIT.R6_n"],
    "INIT.R6_p": ["INIT.R6_p"],
    "INIT.R6_s": ["INIT.R6_s"],
    "INIT.R8_c": ["INIT.R8_c"],
    "INIT.R8_n": ["INIT.R8_n"],
    "INIT.R8_p": ["INIT.R8_p"],
    "INIT.R8_s": ["INIT.R8_s"],
    "INIT.X1_c": ["INIT.X1_c"],
    "INIT.X2_c": ["INIT.X2_c"],
    "INIT.X3_c": ["INIT.X3_c"],
    "INIT.B1_c":   ["INIT.B1_m5_c"],
    "INIT.B1_n":   ["INIT.B1_m5_n"],
    "INIT.B1_p":   ["INIT.B1_m5_p"],
    "INIT.P1_c":   ["INIT.P1_11_c", "INIT.P1_12_c", "INIT.P1_13_c", "INIT.P1_14_c", "INIT.P1_15_c", "INIT.P1_16_c", "INIT.P1_17_c", "INIT.P1_18_c", "INIT.P1_19_c", "INIT.P1_20_c"],
    "INIT.P1_n":   ["INIT.P1_11_n", "INIT.P1_12_n", "INIT.P1_13_n", "INIT.P1_14_n", "INIT.P1_15_n", "INIT.P1_16_n", "INIT.P1_17_n", "INIT.P1_18_n", "INIT.P1_19_n", "INIT.P1_20_n"],
    "INIT.P1_p":   ["INIT.P1_11_p", "INIT.P1_12_p", "INIT.P1_13_p", "INIT.P1_14_p", "INIT.P1_15_p", "INIT.P1_16_p", "INIT.P1_17_p", "INIT.P1_18_p", "INIT.P1_19_p", "INIT.P1_20_p"],
    "INIT.P1_Chl": ["INIT.P1_11_Chl", "INIT.P1_12_Chl", "INIT.P1_13_Chl", "INIT.P1_14_Chl", "INIT.P1_15_Chl", "INIT.P1_16_Chl", "INIT.P1_17_Chl", "INIT.P1_18_Chl", "INIT.P1_19_Chl", "INIT.P1_20_Chl"],
    "INIT.P1_s":   ["INIT.P1_11_s", "INIT.P1_12_s", "INIT.P1_13_s", "INIT.P1_14_s", "INIT.P1_15_s", "INIT.P1_16_s", "INIT.P1_17_s", "INIT.P1_18_s", "INIT.P1_19_s", "INIT.P1_20_s"],
    "INIT.P4_c":   ["INIT.P4_11_c", "INIT.P4_12_c", "INIT.P4_13_c", "INIT.P4_14_c", "INIT.P4_15_c", "INIT.P4_16_c", "INIT.P4_17_c", "INIT.P4_18_c", "INIT.P4_19_c", "INIT.P4_20_c"],
    "INIT.P4_n":   ["INIT.P4_11_n", "INIT.P4_12_n", "INIT.P4_13_n", "INIT.P4_14_n", "INIT.P4_15_n", "INIT.P4_16_n", "INIT.P4_17_n", "INIT.P4_18_n", "INIT.P4_19_n", "INIT.P4_20_n"],
    "INIT.P4_p":   ["INIT.P4_11_p", "INIT.P4_12_p", "INIT.P4_13_p", "INIT.P4_14_p", "INIT.P4_15_p", "INIT.P4_16_p", "INIT.P4_17_p", "INIT.P4_18_p", "INIT.P4_19_p", "INIT.P4_20_p"],
    "INIT.P4_Chl": ["INIT.P4_11_Chl", "INIT.P4_12_Chl", "INIT.P4_13_Chl", "INIT.P4_14_Chl", "INIT.P4_15_Chl", "INIT.P4_16_Chl", "INIT.P4_17_Chl", "INIT.P4_18_Chl", "INIT.P4_19_Chl", "INIT.P4_20_Chl"],
    "INIT.P2_c":   ["INIT.P2_4_c", "INIT.P2_5_c", "INIT.P2_6_c", "INIT.P2_7_c", "INIT.P2_8_c", "INIT.P2_9_c", "INIT.P2_10_c", "INIT.P2_11_c", "INIT.P2_12_c", "INIT.P2_13_c","INIT.P5_6_c", "INIT.P5_7_c", "INIT.P5_8_c", "INIT.P5_9_c", "INIT.P5_10_c", "INIT.P5_11_c", "INIT.P5_12_c", "INIT.P5_13_c", "INIT.P5_14_c", "INIT.P5_15_c","INIT.P7_2_c", "INIT.P7_3_c", "INIT.P7_4_c", "INIT.P7_5_c", "INIT.P7_6_c", "INIT.P7_7_c","INIT.P8_2_c", "INIT.P8_3_c", "INIT.P8_4_c", "INIT.P8_5_c", "INIT.P8_6_c", "INIT.P8_7_c"],
    "INIT.P2_n":   ["INIT.P2_4_n", "INIT.P2_5_n", "INIT.P2_6_n", "INIT.P2_7_n", "INIT.P2_8_n", "INIT.P2_9_n", "INIT.P2_10_n", "INIT.P2_11_n", "INIT.P2_12_n", "INIT.P2_13_n","INIT.P5_6_n", "INIT.P5_7_n", "INIT.P5_8_n", "INIT.P5_9_n", "INIT.P5_10_n", "INIT.P5_11_n", "INIT.P5_12_n", "INIT.P5_13_n", "INIT.P5_14_n", "INIT.P5_15_n","INIT.P7_2_n", "INIT.P7_3_n", "INIT.P7_4_n", "INIT.P7_5_n", "INIT.P7_6_n", "INIT.P7_7_n","INIT.P8_2_n", "INIT.P8_3_n", "INIT.P8_4_n", "INIT.P8_5_n", "INIT.P8_6_n", "INIT.P8_7_n"],
    "INIT.P2_p":   ["INIT.P2_4_p", "INIT.P2_5_p", "INIT.P2_6_p", "INIT.P2_7_p", "INIT.P2_8_p", "INIT.P2_9_p", "INIT.P2_10_p", "INIT.P2_11_p", "INIT.P2_12_p", "INIT.P2_13_p","INIT.P5_6_p", "INIT.P5_7_p", "INIT.P5_8_p", "INIT.P5_9_p", "INIT.P5_10_p", "INIT.P5_11_p", "INIT.P5_12_p", "INIT.P5_13_p", "INIT.P5_14_p", "INIT.P5_15_p","INIT.P7_2_p", "INIT.P7_3_p", "INIT.P7_4_p", "INIT.P7_5_p", "INIT.P7_6_p", "INIT.P7_7_p","INIT.P8_2_p", "INIT.P8_3_p", "INIT.P8_4_p", "INIT.P8_5_p", "INIT.P8_6_p", "INIT.P8_7_p"],
    "INIT.P2_Chl": ["INIT.P2_4_Chl", "INIT.P2_5_Chl", "INIT.P2_6_Chl", "INIT.P2_7_Chl", "INIT.P2_8_Chl", "INIT.P2_9_Chl", "INIT.P2_10_Chl", "INIT.P2_11_Chl", "INIT.P2_12_Chl", "INIT.P2_13_Chl","INIT.P5_6_Chl", "INIT.P5_7_Chl", "INIT.P5_8_Chl", "INIT.P5_9_Chl", "INIT.P5_10_Chl", "INIT.P5_11_Chl", "INIT.P5_12_Chl", "INIT.P5_13_Chl", "INIT.P5_14_Chl", "INIT.P5_15_Chl","INIT.P7_2_Chl", "INIT.P7_3_Chl", "INIT.P7_4_Chl", "INIT.P7_5_Chl", "INIT.P7_6_Chl", "INIT.P7_7_Chl","INIT.P8_2_Chl", "INIT.P8_3_Chl", "INIT.P8_4_Chl", "INIT.P8_5_Chl", "INIT.P8_6_Chl", "INIT.P8_7_Chl"],
    "INIT.P3_c":   ["INIT.P3_0_c", "INIT.P3_1_c", "INIT.P3_2_c","INIT.P9_m2_c", "INIT.P9_m1_c", "INIT.P9_0_c","INIT.P6_m4_c", "INIT.P6_m3_c"],
    "INIT.P3_n":   ["INIT.P3_0_n", "INIT.P3_1_n", "INIT.P3_2_n","INIT.P9_m2_n", "INIT.P9_m1_n", "INIT.P9_0_n","INIT.P6_m4_n", "INIT.P6_m3_n"],
    "INIT.P3_p":   ["INIT.P3_0_p", "INIT.P3_1_p", "INIT.P3_2_p","INIT.P9_m2_p", "INIT.P9_m1_p", "INIT.P9_0_p","INIT.P6_m4_p", "INIT.P6_m3_p"],
    "INIT.P3_Chl": ["INIT.P3_0_Chl", "INIT.P3_1_Chl", "INIT.P3_2_Chl","INIT.P9_m2_Chl", "INIT.P9_m1_Chl", "INIT.P9_0_Chl","INIT.P6_m4_Chl", "INIT.P6_m3_Chl"],
    "INIT.Z6_c":   ["INIT.Z6_5_c", "INIT.Z6_6_c", "INIT.Z6_7_c", "INIT.Z6_8_c", "INIT.Z6_9_c", "INIT.Z6_10_c", "INIT.Z6_11_c"],
    "INIT.Z6_n":   ["INIT.Z6_5_n", "INIT.Z6_6_n", "INIT.Z6_7_n", "INIT.Z6_8_n", "INIT.Z6_9_n", "INIT.Z6_10_n", "INIT.Z6_11_n"],
    "INIT.Z6_p":   ["INIT.Z6_5_p", "INIT.Z6_6_p", "INIT.Z6_7_p", "INIT.Z6_8_p", "INIT.Z6_9_p", "INIT.Z6_10_p", "INIT.Z6_11_p"],
    "INIT.Z5_c":   ["INIT.Z5_9_c", "INIT.Z5_10_c", "INIT.Z5_11_c", "INIT.Z5_12_c", "INIT.Z5_13_c", "INIT.Z5_14_c", "INIT.Z5_15_c", "INIT.Z5_16_c", "INIT.Z5_17_c", "INIT.Z5_18_c"],
    "INIT.Z5_n":   ["INIT.Z5_9_n", "INIT.Z5_10_n", "INIT.Z5_11_n", "INIT.Z5_12_n", "INIT.Z5_13_n", "INIT.Z5_14_n", "INIT.Z5_15_n", "INIT.Z5_16_n", "INIT.Z5_17_n", "INIT.Z5_18_n"],
    "INIT.Z5_p":   ["INIT.Z5_9_p", "INIT.Z5_10_p", "INIT.Z5_11_p", "INIT.Z5_12_p", "INIT.Z5_13_p", "INIT.Z5_14_p", "INIT.Z5_15_p", "INIT.Z5_16_p", "INIT.Z5_17_p", "INIT.Z5_18_p"],
    "INIT.Z4_c":   ["INIT.Z4_15_c", "INIT.Z4_16_c", "INIT.Z4_17_c", "INIT.Z4_18_c", "INIT.Z4_19_c", "INIT.Z4_20_c", "INIT.Z4_21_c", "INIT.Z4_22_c", "INIT.Z4_23_c", "INIT.Z4_24_c"],
    "INIT.Z4_n":   ["INIT.Z4_15_n", "INIT.Z4_16_n", "INIT.Z4_17_n", "INIT.Z4_18_n", "INIT.Z4_19_n", "INIT.Z4_20_n", "INIT.Z4_21_n", "INIT.Z4_22_n", "INIT.Z4_23_n", "INIT.Z4_24_n"],
    "INIT.Z4_p":   ["INIT.Z4_15_p", "INIT.Z4_16_p", "INIT.Z4_17_p", "INIT.Z4_18_p", "INIT.Z4_19_p", "INIT.Z4_20_p", "INIT.Z4_21_p", "INIT.Z4_22_p", "INIT.Z4_23_p", "INIT.Z4_24_p"],
    "INIT.Z3_c":   ["INIT.Z3_22_c", "INIT.Z3_23_c", "INIT.Z3_24_c", "INIT.Z3_25_c", "INIT.Z3_26_c", "INIT.Z3_27_c", "INIT.Z3_28_c", "INIT.Z3_29_c", "INIT.Z3_30_c", "INIT.Z3_31_c"],
    "INIT.Z3_n":   ["INIT.Z3_22_n", "INIT.Z3_23_n", "INIT.Z3_24_n", "INIT.Z3_25_n", "INIT.Z3_26_n", "INIT.Z3_27_n", "INIT.Z3_28_n", "INIT.Z3_29_n", "INIT.Z3_30_n", "INIT.Z3_31_n"],
    "INIT.Z3_p":   ["INIT.Z3_22_p", "INIT.Z3_23_p", "INIT.Z3_24_p", "INIT.Z3_25_p", "INIT.Z3_26_p", "INIT.Z3_27_p", "INIT.Z3_28_p", "INIT.Z3_29_p", "INIT.Z3_30_p", "INIT.Z3_31_p"]
}

def read_values(filepath):
    values = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            values.append(float(line))
    return values

def write_values(filepath, values):
    with open(filepath, "w") as f:
        for v in values:
            f.write(f"{v:.20f}\n")

def main():
    for infile in input_dir.iterdir():
        if not infile.is_file():
            continue

        old_name = infile.name

        # considera solo i file presenti nel dizionario
        if old_name not in mapping:
            print(f"Salto {old_name}: non presente nel dizionario")
            continue

        new_names = mapping[old_name]
        nout = len(new_names)

        if nout == 0:
            print(f"Salto {old_name}: lista di output vuota")
            continue

        try:
            values = read_values(infile)
        except Exception as e:
            print(f"Errore leggendo {old_name}: {e}")
            continue

        # dividi i valori per il numero di file di uscita
        new_values = [v / nout for v in values]

        for new_name in new_names:
            outfile = output_dir / new_name
            try:
                write_values(outfile, new_values)
                print(f"Creato {outfile}")
            except Exception as e:
                print(f"Errore scrivendo {outfile}: {e}")

if __name__ == "__main__":
    main()
