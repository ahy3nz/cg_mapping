from lxml import etree


outputFile = "chol.xml"
root = etree.Element('Mapping')
mapping = etree.SubElement(root,"Molecule", name=outputFile[:-4])
beadtypes = etree.SubElement(mapping, 'Beads')
bonds = etree.SubElement(mapping, 'Bonds')

etree.SubElement(beadtypes, "Bead", index="0", beadtype='chead',
        map='0 2')
etree.SubElement(beadtypes, "Bead", index="1", beadtype='ring',
    map='0 1 15 16 17 18')
etree.SubElement(beadtypes, "Bead", index="2", beadtype='ring',
    map='4 5 6 7 8 9 10')
etree.SubElement(beadtypes, "Bead", index="3", beadtype='ring',
    map='26 27 28 29 30 31 32 33')
etree.SubElement(beadtypes, "Bead", index="4", beadtype='ring',
    map='34 47 48')
etree.SubElement(beadtypes, "Bead", index="5", beadtype='ring',
    map='39 40 41 42 43 44 45 46')
etree.SubElement(beadtypes, "Bead", index="6", beadtype='ring',
    map='19 20 21 22 23 24 25')
etree.SubElement(beadtypes, "Bead", index="7", beadtype='ctail',
    map='49 50 51 52 53 54 55 56 57 58 59 60')
etree.SubElement(beadtypes, "Bead", index="8", beadtype='cterm',
    map='61 62 63 64 65 66 67 68 69 70 71 72 73')
etree.SubElement(beadtypes, "Bead", index="9", beadtype='chme',
    map='11 12 13 14')
etree.SubElement(beadtypes, "Bead", index="10", beadtype='chme',
    map='35 36 37 38')

etree.SubElement(bonds, 'Bond', bead1='0', bead2='1')
etree.SubElement(bonds, 'Bond', bead1='1', bead2='2')
etree.SubElement(bonds, 'Bond', bead1='1', bead2='6')
etree.SubElement(bonds, 'Bond', bead1='2', bead2='3')
etree.SubElement(bonds, 'Bond', bead1='2', bead2='9')
etree.SubElement(bonds, 'Bond', bead1='3', bead2='4')
etree.SubElement(bonds, 'Bond', bead1='3', bead2='6')
etree.SubElement(bonds, 'Bond', bead1='4', bead2='5')
etree.SubElement(bonds, 'Bond', bead1='4', bead2='7')
etree.SubElement(bonds, 'Bond', bead1='4', bead2='10')
etree.SubElement(bonds, 'Bond', bead1='5', bead2='6')
etree.SubElement(bonds, 'Bond', bead1='7', bead2='8')
etree.SubElement(bonds, 'Bond', bead1='0', bead2='3')
etree.SubElement(bonds, 'Bond', bead1='0', bead2='4')
etree.SubElement(bonds, 'Bond', bead1='0', bead2='5')
etree.SubElement(bonds, 'Bond', bead1='0', bead2='7')
etree.SubElement(bonds, 'Bond', bead1='0', bead2='8')
etree.SubElement(bonds, 'Bond', bead1='0', bead2='9')
etree.SubElement(bonds, 'Bond', bead1='0', bead2='10')
etree.SubElement(bonds, 'Bond', bead1='1', bead2='4')
etree.SubElement(bonds, 'Bond', bead1='1', bead2='7')
etree.SubElement(bonds, 'Bond', bead1='1', bead2='8')
etree.SubElement(bonds, 'Bond', bead1='1', bead2='10')
etree.SubElement(bonds, 'Bond', bead1='2', bead2='5')
etree.SubElement(bonds, 'Bond', bead1='2', bead2='7')
etree.SubElement(bonds, 'Bond', bead1='2', bead2='8')
etree.SubElement(bonds, 'Bond', bead1='2', bead2='10')
etree.SubElement(bonds, 'Bond', bead1='3', bead2='8')
etree.SubElement(bonds, 'Bond', bead1='4', bead2='9')
etree.SubElement(bonds, 'Bond', bead1='5', bead2='8')
etree.SubElement(bonds, 'Bond', bead1='5', bead2='9')
etree.SubElement(bonds, 'Bond', bead1='6', bead2='7')
etree.SubElement(bonds, 'Bond', bead1='6', bead2='8')
etree.SubElement(bonds, 'Bond', bead1='6', bead2='9')
etree.SubElement(bonds, 'Bond', bead1='6', bead2='10')
etree.SubElement(bonds, 'Bond', bead1='7', bead2='9')
etree.SubElement(bonds, 'Bond', bead1='8', bead2='9')
etree.SubElement(bonds, 'Bond', bead1='8', bead2='10')
etree.SubElement(bonds, 'Bond', bead1='9', bead2='10')





bonds.set('n_bonds', str(len(bonds.getchildren())))
beadtypes.set('n_beads', str(len(beadtypes.getchildren())))
tree =  etree.ElementTree(root)
tree.write(outputFile, pretty_print=True)
