import sys

def raise_bgc_error(received_type):
    
    valid_types = "[DEFAULT[Default,default], BFM[bfm], FABM-BFM[fabm-bfm]]"
    print(f"[ERROR] BGC_TYPE='{received_type}' is wrong/undefined. Expected one of: {valid_types}")
    sys.exit()