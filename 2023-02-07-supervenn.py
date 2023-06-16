from supervenn import supervenn
from tqdm import tqdm
import requests
import matplotlib.pyplot as plt
import time

verified = {}
def verify_drug(name):
    if name not in verified:
        res = requests.get("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/%s/property/title/JSON"%name)
        time.sleep(0.1)
        if not res.ok:
            verified[name] = name
        else:
            verified[name] = res.json()["PropertyTable"]["Properties"][0].get("Title", name)
    return verified[name]

drug_sets = {}
with open("data/resources.dmt") as o:
    for line in o:
        term, _, *drugs = line.strip().split("\t")
        drug_sets[term] = set(drugs)

mapped_drug_sets = {}
for k,v in drug_sets.items():
    mapped_drug_sets[k] = set()
    for i in tqdm(v):
         mapped_drug_sets[k].add(verify_drug(i))
    print(k, len(v), len(mapped_drug_sets[k]))

sets = []
labels = []
for k, v in mapped_drug_sets.items():
    if len(v) > 1 and not '20' in k:
        sets.append(v)
        labels.append(k)
plt.figure(figsize=(16, 8))
supervenn(sets, labels, sets_ordering= 'minimize gaps', widths_minmax_ratio=0.1 )
plt.savefig("supervenn.png", dpi=500)
plt.savefig("supervenn.svg", dpi=500)
