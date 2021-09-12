import xml.etree.ElementTree as ET
import csv

def ipr2pf(infile):
	with open('ipr2pf.csv', mode = 'w') as outfile:
		fieldnames = ['interpro', 'short_name', 'type', 'xrefdb', 'xref']
		outfile = csv.DictWriter(outfile, fieldnames = fieldnames)
		outfile.writeheader()
	
		tree = ET.parse(infile)
		root = tree.getroot()

		for interpro in root.findall('interpro/member_list/db_xref[@db="PFAM"]/../..'):
			for db_xref in interpro.findall('member_list/db_xref'):
				outfile.writerow({'interpro': interpro.get('id'), 'short_name': interpro.get('short_name'), 'type': interpro.get('type'), 'xrefdb': db_xref.get('db'), 'xref': db_xref.get('dbkey')})
	
ipr2pf('interpro.xml')
