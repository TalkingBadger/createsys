
def get_pdbfilename_HCNH2O(seed: int):
    seed = int(seed)
    if seed not in range(20):
        raise ValueError("HCN/H2O system is only available for seeds between 0 and 19")
    return f"/home/abeckmann/ma/potqc/exps/40-hcn_h2o_configs/configs/hcn_h2o_{seed:02}.pdb"

def get_pdbfilename_fixed_26HCN10H2O(seed: int, qm_only=False):
    seed = int(seed)
    if seed not in range(19):
        raise ValueError("fixed ratio HCN/H2O system is only available for seeds between 0 and 18")
    return f"/home/abeckmann/ma/potqc/exps/47-HCN-H2O-fixed-ratio-configs/{'qm' if qm_only else ''}configs/{seed:02}.pdb"

def get_qmatoms_HCNH2O(seed: int):
    seed = int(seed)
    if seed not in range(20):
        raise ValueError("HCN/H2O system is only available for seeds between 0 and 19")
    return f"/home/abeckmann/ma/potqc/exps/40-hcn_h2o_configs/configs/qm_atoms_{seed:02}"

def get_pdbfilename_nucleobase(nucleobase: str, seed: int):
    seed = int(seed)
    if seed not in range(15):
        raise ValueError("Nucleobase systems are only available for seed between 0 and 14.")
    letter = nucleobase[0].upper()
    if letter not in "ACGT":
        raise ValueError(f"Nucleobase {nucleobase} is not one of (A)denine, (C)ytosin, (T)hymine, nor (G)uanine")
    id = f"{nucleobase[0].upper()}{seed:02}"
    return f"/home/abeckmann/ma/potqc/exps/43-nucleosides-seedconfs/{id}/F_{id}.pdb"


if __name__ == "__main__":
    print(get_pdbfilename_HCNH2O(19))
    print(get_qmatoms_HCNH2O(4))
    print(get_pdbfilename_nucleobase("thymine", 10))
    print(get_pdbfilename_fixed_26HCN10H2O(10))

