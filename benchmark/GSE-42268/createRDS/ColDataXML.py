from xml.dom.minidom import parse, parseString
import os

os.chdir("/Users/jackpeplinski/CycleMix/benchmark")
# ignore all phase when it's in a matrix
dom = parse("../benchmarkData/GSE42268/GSE42268_family.xml")
samples = dom.getElementsByTagName("Sample")
accs = []
orgs = []
cell_types = []
phases = []
for sample in samples:
    acc = sample.getElementsByTagName("Accession").item(0).firstChild.nodeValue
    org = sample.getElementsByTagName("Organism").item(0).firstChild.nodeValue
    characteristics = sample.getElementsByTagName("Characteristics")
    for characteristic in characteristics:
        # print(characteristic.getAttribute("tag") == "cell cycle phase")
        if characteristic.getAttribute("tag") == "cell type":
            cell_type = characteristic.firstChild.nodeValue.replace("\n", "").rstrip()
        if characteristic.getAttribute("tag") == "cell cycle phase":
            phase = characteristic.firstChild.nodeValue.replace("\n", "").rstrip()
    cell_types.append(cell_type)
    accs.append(acc)
    orgs.append(org)
    phases.append(phase.strip())

# Remove "All phase" from the list
indices = [i for i, x in enumerate(phases) if x == "All phase"]
for idx, index in enumerate(indices):
    accs.pop(index - idx)
    cell_types.pop(index - idx)
    phases.pop(index - idx)

print(accs)
# print("####")
# print(orgs)
# print("###")
# print(cell_types)
# print("###")
print(phases)
