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
    "INIT.B1_c": ["INIT.B1_m5"],
    "INIT.P1_c": ["INIT.P1_11_c", "INIT.P1_12_c", "INIT.P1_13_c", "INIT.P1_14_c", "INIT.P1_15_c", "INIT.P1_16_c", "INIT.P1_17_c", "INIT.P1_18_c", "INIT.P1_19_c", "INIT.P1_20_c"],
    "INIT.P1_n": ["INIT.P1_11_n", "INIT.P1_12_n", "INIT.P1_13_n", "INIT.P1_14_n", "INIT.P1_15_n", "INIT.P1_16_n", "INIT.P1_17_n", "INIT.P1_18_n", "INIT.P1_19_n", "INIT.P1_20_n"],
    "INIT.P1_p": ["INIT.P1_11_p", "INIT.P1_12_p", "INIT.P1_13_p", "INIT.P1_14_p", "INIT.P1_15_p", "INIT.P1_16_p", "INIT.P1_17_p", "INIT.P1_18_p", "INIT.P1_19_p", "INIT.P1_20_p"],
    "INIT.P1_Chl":  ["INIT.P1_11_Chl", "INIT.P1_12_Chl", "INIT.P1_13_Chl", "INIT.P1_14_Chl", "INIT.P1_15_Chl", "INIT.P1_16_Chl", "INIT.P1_17_Chl", "INIT.P1_18_Chl", "INIT.P1_19_Chl", "INIT.P1_20_Chl"],
    "INIT.P1_s": ["INIT.P1_11_s", "INIT.P1_12_s", "INIT.P1_13_s", "INIT.P1_14_s", "INIT.P1_15_s", "INIT.P1_16_s", "INIT.P1_17_s", "INIT.P1_18_s", "INIT.P1_19_s", "INIT.P1_20_s"], 
    "INIT.P4_c": ["INIT.P4_11_c", "INIT.P4_12_c", "INIT.P4_13_c", "INIT.P4_14_c", "INIT.P4_15_c", "INIT.P4_16_c", "INIT.P4_17_c", "INIT.P4_18_c", "INIT.P4_19_c", "INIT.P4_20_c"],
    "INIT.P4_n": ["INIT.P4_11_n", "INIT.P4_12_n", "INIT.P4_13_n", "INIT.P4_14_n", "INIT.P4_15_n", "INIT.P4_16_n", "INIT.P4_17_n", "INIT.P4_18_n", "INIT.P4_19_n", "INIT.P4_20_n"],
    "INIT.P4_p": ["INIT.P4_11_p", "INIT.P4_12_p", "INIT.P4_13_p", "INIT.P4_14_p", "INIT.P4_15_p", "INIT.P4_16_p", "INIT.P4_17_p", "INIT.P4_18_p", "INIT.P4_19_p", "INIT.P4_20_p"],
    "INIT.P4_Chl":  ["INIT.P4_11_Chl", "INIT.P4_12_Chl", "INIT.P4_13_Chl", "INIT.P4_14_Chl", "INIT.P4_15_Chl", "INIT.P4_16_Chl", "INIT.P4_17_Chl", "INIT.P4_18_Chl", "INIT.P4_19_Chl", "INIT.P4_20_Chl"],
    "INIT.P2_c": ["INIT.P2_4_c", "INIT.P2_5_c", "INIT.P2_6_c", "INIT.P2_7_c", "INIT.P2_8_c", "INIT.P2_9_c", "INIT.P2_10_c", "INIT.P2_11_c", "INIT.P2_12_c", "INIT.P2_13_c"],
    "INIT.P2_n": ["INIT.P2_4_n", "INIT.P2_5_n", "INIT.P2_6_n", "INIT.P2_7_n", "INIT.P2_8_n", "INIT.P2_9_n", "INIT.P2_10_n", "INIT.P2_11_n", "INIT.P2_12_n", "INIT.P2_13_n"],
    "INIT.P2_p": ["INIT.P2_4_p", "INIT.P2_5_p", "INIT.P2_6_p", "INIT.P2_7_p", "INIT.P2_8_p", "INIT.P2_9_p", "INIT.P2_10_p", "INIT.P2_11_p", "INIT.P2_12_p", "INIT.P2_13_p"],
    "INIT.P2_Chl":  ["INIT.P2_4_Chl", "INIT.P2_5_Chl", "INIT.P2_6_Chl", "INIT.P2_7_Chl", "INIT.P2_8_Chl", "INIT.P2_9_Chl", "INIT.P2_10_Chl", "INIT.P2_11_Chl", "INIT.P2_12_Chl", "INIT.P2_13_Chl"],   
    "INIT.P5_c": ["INIT.P5_6_c", "INIT.P5_7_c", "INIT.P5_8_c", "INIT.P5_9_c", "INIT.P5_10_c", "INIT.P5_11_c", "INIT.P5_12_c", "INIT.P5_13_c", "INIT.P5_14_c", "INIT.P5_15_c"],
    "INIT.P5_n": ["INIT.P5_6_n", "INIT.P5_7_n", "INIT.P5_8_n", "INIT.P5_9_n", "INIT.P5_10_n", "INIT.P5_11_n", "INIT.P5_12_n", "INIT.P5_13_n", "INIT.P5_14_n", "INIT.P5_15_n"],
    "INIT.P5_p": ["INIT.P5_6_p", "INIT.P5_7_p", "INIT.P5_8_p", "INIT.P5_9_p", "INIT.P5_10_p", "INIT.P5_11_p", "INIT.P5_12_p", "INIT.P5_13_p", "INIT.P5_14_p", "INIT.P5_15_p"],
    "INIT.P5_Chl":  ["INIT.P5_6_Chl", "INIT.P5_7_Chl", "INIT.P5_8_Chl", "INIT.P5_9_Chl", "INIT.P5_10_Chl", "INIT.P5_11_Chl", "INIT.P5_12_Chl", "INIT.P5_13_Chl", "INIT.P5_14_Chl", "INIT.P5_15_Chl"],
    "INIT.P7_c": ["INIT.P7_2_c", "INIT.P7_3_c", "INIT.P7_4_c", "INIT.P7_5_c", "INIT.P7_6_c", "INIT.P7_7_c"],
    "INIT.P7_n": ["INIT.P7_2_n", "INIT.P7_3_n", "INIT.P7_4_n", "INIT.P7_5_n", "INIT.P7_6_n", "INIT.P7_7_n"],
    "INIT.P7_p": ["INIT.P7_2_p", "INIT.P7_3_p", "INIT.P7_4_p", "INIT.P7_5_p", "INIT.P7_6_p", "INIT.P7_7_p"],
    "INIT.P7_Chl":  ["INIT.P7_2_Chl", "INIT.P7_3_Chl", "INIT.P7_4_Chl", "INIT.P7_5_Chl", "INIT.P7_6_Chl", "INIT.P7_7_Chl"],
    "INIT.P8_c": ["INIT.P8_2_c", "INIT.P8_3_c", "INIT.P8_4_c", "INIT.P8_5_c", "INIT.P8_6_c", "INIT.P8_7_c"],
    "INIT.P8_n": ["INIT.P8_2_n", "INIT.P8_3_n", "INIT.P8_4_n", "INIT.P8_5_n", "INIT.P8_6_n", "INIT.P8_7_n"],
    "INIT.P8_p": ["INIT.P8_2_p", "INIT.P8_3_p", "INIT.P8_4_p", "INIT.P8_5_p", "INIT.P8_6_p", "INIT.P8_7_p"],
    "INIT.P8_Chl":  ["INIT.P8_2_Chl", "INIT.P8_3_Chl", "INIT.P8_4_Chl", "INIT.P8_5_Chl", "INIT.P8_6_Chl", "INIT.P8_7_Chl"],
   # aggiungi qui gli altri mapping
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