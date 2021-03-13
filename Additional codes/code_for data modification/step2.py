from rdkit import Chem
import random
with open('all-file.txt','r')as f1,open('file.source','w')as f2,open('file.target','w')as f3,open('random_reaction-train.txt','w')as f4:
    l=f1.readlines()
    for i in l:
        smi1=i.replace('\n', '').split('>>')[0]
        smi2=i.replace('\n', '').split('>>')[1]
        mol1=Chem.MolFromSmiles(smi1)
        mol2=Chem.MolFromSmiles(smi2)
        map1 = list((a.GetAtomMapNum() for a in mol1.GetAtoms() if a.HasProp('molAtomMapNumber')))
        map2 = list((a.GetAtomMapNum() for a in mol2.GetAtoms() if a.HasProp('molAtomMapNumber')))
        r = []
        rr=[]
        smap = [i for i in map1 if i in map2]
        dmap=[i for i in map1 if i not in map2]+[i for i in map2 if i not in map1]
        d = smap.copy()
        dd=dmap.copy()
        for i in smap:
            n = random.choice(d)
            smi1 = smi1.replace(':' + str(i) + ']', ':' + str(n+1000) + ']')
            smi2 = smi2.replace(':' + str(i) + ']', ':' + str(n+1000) + ']')
            r.append(n + 1000)
            d.remove(n)
        for i in dmap:
            n = random.choice(dd)
            smi1 = smi1.replace(':' + str(i) + ']', ':' + str(n+1000) + ']')
            smi2 = smi2.replace(':' + str(i) + ']', ':' + str(n+1000) + ']')
            rr.append(n + 1000)
            dd.remove(n)
        for a in r:
             smi1 = smi1.replace(':' + str(a) + ']', ':' + str(a - 1000) + ']')
             smi2 = smi2.replace(':' + str(a) + ']', ':' + str(a - 1000) + ']')
        for a in rr:
             smi1 = smi1.replace(':' + str(a) + ']', ':' + str(a - 1000) + ']')
             smi2 = smi2.replace(':' + str(a) + ']', ':' + str(a - 1000) + ']')
        rxn = smi1 + '>>' + smi2
        f2.write(smi1+'\n')
        f3.write(smi2 + '\n')
        f4.write(rxn + '\n')

