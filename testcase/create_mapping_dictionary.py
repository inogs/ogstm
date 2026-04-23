import re
import os
from collections import Counter

# Choose the file and path
yaml_file_path = "fabm.yaml"

############# FUNCTION ###################
def generate_mapping(file_path):
    if not os.path.exists(file_path):
      print(f"Error: File '{file_path}' not found.")
      return

    # 1. Extract all the environment variables and size classes (eg. N5, B1_m5, P2_4, Z3, etc.)
    pattern = re.compile(r' (O3h|[NORXBPZ][0-9](?:_[A-Za-z0-9]{1,2})?):')      # [WARNING#1] Use "O3h|" to catch the specific exception, duplicates are filtered below
    
    # Initialize the list and the set
    items = ["depth"]                                                          # [WARNING#2] Hardcoded "depth" since it is not defined in the .yaml files
    seen = set(["depth"])
    
    with open(file_path, 'r') as f:
      for line in f:
        match = pattern.search(line)
        if match:
          val = match.group(1)
          if val not in seen:     # Check for duplicate entries
            items.append(val)     # Append the new entry to the list
            seen.add(val)         # Add the new entry to the set

    if len(items) == 1: # Only "depth" is in there
      print("No matching identifiers found in the YAML file.")
      return

    summary = Counter(item[0] for item in items)
    
    # Print the summary
    print("--- Summary ---")
    for key in sorted(summary.keys()):
      print(f"{key}: {summary[key]}")
    print(f"TOTAL: {len(items)}")
    print("=============================")
    
    # Print all the founded types  
    print("--- Extracted Types ---")
    for item in items:
      print(f"- {item}")
    print("=============================")


    # 2. Group the environment variables and size classes by their main definition (e.g., N, B1, P4, Z3)
    groups = {}
    for item in items:
      family = item.split('_')[0]
      if family not in groups:
        groups[family] = []
      groups[family].append(item)
        
    print("--- Sorted Groups ---")
    for group in groups.items():
      print(f"- {group}")
    print("=============================")


    # 3. Build the mapping dictionary
    mapping = {}

    # Suffix map for environment variables
    env_suffixes = {
      'N1': ['p'], 'N3': ['n'], 'N4': ['n'], 'N5': ['s'], 'N6': ['r'],
      'O2': ['o'], 'O3': ['c'], 'O3h': ['h'], 'O4': ['n'], 'O5': ['c'],
      'R1': ['c', 'n', 'p', 's'],
      'R2': ['c'],
      'R3': ['c'],
      'R6': ['c', 'n', 'p', 's'],
      'R8': ['c', 'n', 'p', 's'],
      'X1': ['c'], 'X2': ['c'], 'X3': ['c']
    }
        
    # Define the phytoplankton exceptions (P3-P6-P9, P2-P5-P7-P8) for the dictionary's framework
    merger_map = {
      'P6': 'P3',
      'P9': 'P3',
      'P5': 'P2',
      'P7': 'P2',
      'P8': 'P2'
    }
    
    # Loop on all the state variables
    mapping["INIT.depth"] = ["INIT.depth"]
    for family, members in groups.items():
      # Environment state variables
      if family in env_suffixes:
        for suffix in env_suffixes[family]:
          key = f"INIT.{family}_{suffix}"
          if key not in mapping: mapping[key] = []
          mapping[key].extend([f"INIT.{m}_{suffix}" for m in members])
      else:
      
        # Organism state variables
        target_family = merger_map.get(family, family)
        if family.startswith('B'):
          suffixes = ['c', 'n', 'p']
        elif family.startswith('P'):
          suffixes = ['c', 'n', 'p', 'Chl']
          if family == 'P1' or target_family == 'P1':  # Apply the 's' (si) for P1 (diatoms)
            suffixes.append('s')
        elif family.startswith('Z'):
          suffixes = ['c', 'n', 'p']
        else:
          continue 
        
        for suffix in suffixes:
          key = f"INIT.{target_family}_{suffix}"
          if key not in mapping:
            mapping[key] = []
          
          formatted_members = [f"INIT.{m}_{suffix}" for m in members]
          mapping[key].extend(formatted_members)


    # 4. Print the mapping dictionary
    print("--- Mapping Dictionary ---")
    print("mapping = {")
    keys = list(mapping.keys())
    if keys:
      max_key_len = max(len(key) for key in keys) + 1
      for i, key in enumerate(keys):
        key_print = f'"{key}":'
        padding = " " * (max_key_len - len(key_print) + 3)
        value_print = ", ".join(f'"{v}"' for v in mapping[key])
        comma = "," if i < len(keys) - 1 else ""
        print(f'    {key_print}{padding}[{value_print}]{comma}')
    
    print("}")


############# RUN #############
generate_mapping(yaml_file_path)