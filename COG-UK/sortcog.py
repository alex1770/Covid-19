# Sort cog_metadata.csv into reverse date order to make it quick to process entries more recent than X

import csv,sys

reader=csv.reader(sys.stdin)
writer=csv.writer(sys.stdout)

headings=next(reader)
writer.writerow(headings)
col=headings.index('sample_date')
rest=list(reader)
rest.sort(key=lambda row: row[col], reverse=True)
for row in rest: writer.writerow(row)

