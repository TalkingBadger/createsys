
def get_pdbfilename_HCNH2O(seed: int):
    seed = int(seed)
    if seed not in range(20):
        raise ValueError("HCN/H2O system is only available for seeds between 0 and 19")
    return f"/home/abeckmann/ma/potqc/exps/40-hcn_h2o_configs/configs/hcn_h2o_{seed:02}.pdb"

def get_qmatoms_HCNH2O(seed: int):
    seed = int(seed)
    if seed not in range(20):
        raise ValueError("HCN/H2O system is only available for seeds between 0 and 19")
    return f"/home/abeckmann/ma/potqc/exps/40-hcn_h2o_configs/configs/qm_atoms_{seed:02}"

if __name__ == "__main__":
    print(get_pdbfilename_HCNH2O(19))
    print(get_qmatoms_HCNH2O(4))