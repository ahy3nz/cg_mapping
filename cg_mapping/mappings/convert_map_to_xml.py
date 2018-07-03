from lxml import etree
import pdb

inputFile = "C16FFA.map"
outputFile = "C16FFA.xml"
root = etree.Element('Mapping')
mapping = etree.SubElement(root,"Molecule", name=outputFile[:-4])
beadtypes = etree.SubElement(mapping, 'Beads')
bonds = etree.SubElement(mapping, 'Bonds')
for line in open(inputFile, 'r').readlines():
    if 'bond' not in line:
        [index, beadtype, indices] = line.split(':')
        new_element = etree.SubElement(beadtypes, "Bead", index=index.strip(),
                                       beadtype=beadtype.strip(), map=indices.strip())
    else: 
        [first, second] = line.split(':')[1].split()
        new_element = etree.SubElement(bonds, "Bond", bead1=first.strip(), 
                                       bead2=second.strip())

bonds.set('n_bonds', str(len(bonds.getchildren())))
beadtypes.set('n_beads', str(len(beadtypes.getchildren())))
tree =  etree.ElementTree(root)
tree.write(outputFile, pretty_print=True)
