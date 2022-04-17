import csv,re,textwrap,sys
inputFile = sys.argv[1]
file = open(inputFile)
csvreader = csv.reader(file)
fastaname = "." + file.name.split(".")[1]
header = next(csvreader)
rows = []
for row in csvreader:
    rows.append(row)
file.close()
rowcount = 0 
output = open(fastaname + ".fasta", 'w')
print('Writing')
for row in rows:
    output.write(">" + row[0] + '\n')
    output.write(row[1].lower() + '\n\n')
print('Finished')
output.close()
