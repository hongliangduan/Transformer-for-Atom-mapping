from rdkit import Chem

def mol_with_atom_index( mol ):
    max_map=max(a.GetAtomMapNum() for a in mol.GetAtoms())
    for a in mol.GetAtoms():
        if not a.HasProp('molAtomMapNumber'):
            a.SetProp( 'molAtomMapNumber',str(max_map+1))
            max_map += 1
    return mol
#
with open('file.txt','r') as f1,open('all-file.txt','w') as f2:
    t1 = f1.readlines()
    for x in range(len(t1)):
        source=t1[x].replace('\n','').split('>>')[0]
        target=t1[x].replace('\n','').split('>>')[1]
        # mol = Chem.MolFromSmiles('CC(=O)[c:5]1[cH:6][cH:7][c:8]2[c:9]([cH:10]1)[O:11][CH2:12][CH2:13][O:14][CH2:15][CH2:16][O:17][CH2:18][CH2:19][O:20][CH2:21][CH2:22][O:23]2.O([C:2]([CH3:1])=[O:3])[OH:4]>>[CH3:1][C:2](=[O:3])[O:4][c:5]1[cH:6][cH:7][c:8]2[c:9]([cH:10]1)[O:11][CH2:12][CH2:13][O:14][CH2:15][CH2:16][O:17][CH2:18][CH2:19][O:20][CH2:21][CH2:22][O:23]2')
        mol = Chem.MolFromSmiles(source)
        mol_M=mol_with_atom_index(mol)
        smi=Chem.MolToSmiles(mol_M)
        rxn=smi+'>>'+target
        f2.write(rxn+'\n')
    # print(smi)
f1.close()
f2.close()
