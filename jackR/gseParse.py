from xml.dom.minidom import parse, parseString

# ignore all phase when it's in a matrix
dom = parse("../jackData/GSE42268/GSE42268_family.xml")
samples = dom.getElementsByTagName("Sample")
accs = []
orgs = []
phases = []
for sample in samples:
    acc = sample.getElementsByTagName("Accession").item(0).firstChild.nodeValue
    org = sample.getElementsByTagName("Organism").item(0).firstChild.nodeValue
    characteristics = sample.getElementsByTagName("Characteristics")
    for characteristic in characteristics:
        # print(characteristic.getAttribute("tag") == "cell cycle phase")
        if characteristic.getAttribute("tag") == "cell cycle phase":
            phase = characteristic.firstChild.nodeValue.replace("\n","").rstrip()
    accs.append(acc)
    orgs.append(org)
    phases.append(phase)

print(accs)
print("####")
print(orgs)
print("###")
print(phases)

