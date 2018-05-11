from lxml import etree

# This is a hard-coded script to desing the DSPC.xml file

outputFile = "DSPC.xml"
root = etree.Element('Mapping')
mapping = etree.SubElement(root,"Molecule", name=outputFile[:-4])
beadtypes = etree.SubElement(mapping, 'Beads')
bonds = etree.SubElement(mapping, 'Bonds')

# Add beads
etree.SubElement(beadtypes, "Bead", index="0", beadtype='PCN', 
                map="0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15")
etree.SubElement(beadtypes, "Bead", index="1", beadtype='PCP',
                map="16 17 18 19 20 21 22 23 24 25 26")
etree.SubElement(beadtypes, "Bead", index="2", beadtype="E1",
                map="27 28 29 30 31")
etree.SubElement(beadtypes, "Bead", index="3", beadtype="E1",
                map="35 36 37 38 39 40")
etree.SubElement(beadtypes, "Bead", index="4", beadtype="C3",
                map="32 33 34 44 45 46 47 48 49")
etree.SubElement(beadtypes, "Bead", index="5", beadtype="C3",
                map="50 51 52 53 54 55 56 57 58")
etree.SubElement(beadtypes, "Bead", index="6", beadtype="C3",
                map="59 60 61 62 63 64 65 66 67")
etree.SubElement(beadtypes, "Bead", index="7", beadtype="C3",
                map="68 69 70 71 72 73 74 75 76")
etree.SubElement(beadtypes, "Bead", index="8", beadtype="C3",
                map="77 78 79 80 81 82 83 84 85")
etree.SubElement(beadtypes, "Bead", index="9", beadtype="C2",
                map="86 87 88 89 90 91 92")
etree.SubElement(beadtypes, "Bead", index="10", beadtype="C3",
                map="41 42 43 93 94 95 96 97 98")
etree.SubElement(beadtypes, "Bead", index="11", beadtype="C3",
                map="99 100 101 102 103 104 105 106 107")
etree.SubElement(beadtypes, "Bead", index="12", beadtype="C3",
                map="108 109 110 111 112 113 114 115 116")
etree.SubElement(beadtypes, "Bead", index="13", beadtype="C3",
                map="117 118 119 120 121 122 123 124 125")
etree.SubElement(beadtypes, "Bead", index="14", beadtype="C3",
                map="126 127 128 129 130 131 132 133 134")
etree.SubElement(beadtypes, "Bead", index="15", beadtype="C2",
                map="135 136 137 138 139 140 141")


# Add bonds
etree.SubElement(bonds, "Bond", bead1="0", bead2="1")
etree.SubElement(bonds, "Bond", bead1="1", bead2="2")
etree.SubElement(bonds, "Bond", bead1="2", bead2="3")
etree.SubElement(bonds, "Bond", bead1="2", bead2="4")
etree.SubElement(bonds, "Bond", bead1="4", bead2="5")
etree.SubElement(bonds, "Bond", bead1="5", bead2="6")
etree.SubElement(bonds, "Bond", bead1="6", bead2="7")
etree.SubElement(bonds, "Bond", bead1="7", bead2="8")
etree.SubElement(bonds, "Bond", bead1="8", bead2="9")
etree.SubElement(bonds, "Bond", bead1="3", bead2="10")
etree.SubElement(bonds, "Bond", bead1="10", bead2="11")
etree.SubElement(bonds, "Bond", bead1="11", bead2="12")
etree.SubElement(bonds, "Bond", bead1="12", bead2="13")
etree.SubElement(bonds, "Bond", bead1="13", bead2="14")
etree.SubElement(bonds, "Bond", bead1="14", bead2="15")


bonds.set('n_bonds', str(len(bonds.getchildren())))
beadtypes.set('n_beads', str(len(beadtypes.getchildren())))
tree =  etree.ElementTree(root)
tree.write(outputFile, pretty_print=True)
