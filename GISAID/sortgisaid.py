# Sort metadata.tsv into reverse date order, and remove entries with unknown date, to make it quick to process recent entries
# Memory requirement: approx 2.9 * (input file size) for typical GISAID file as of 2022-04-07
#
# Example usage:
# tar xf metadata_tsv_2022_07_07.tar.xz metadata.tsv -O | python3 sortgisaid.py > metadata_sorted.tsv

import csv,sys

reader=csv.reader(sys.stdin,delimiter='\t')
writer=csv.writer(sys.stdout,delimiter='\t')

headings=next(reader)
writer.writerow(headings)
col=headings.index('Collection date')
rest=[]
for row in reader:
  if len(row[col])==10: rest.append(row)
rest.sort(key=lambda row: row[col], reverse=True)
for row in rest: writer.writerow(row)
