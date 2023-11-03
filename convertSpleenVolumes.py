import os
import pandas as pd

docname = "Abdominal_features_ukb672614.docx"

f = open("Abdominal_features_ukb672614.docx")
x = f.read()
lines = x.split("\n")
data = []
for elem in lines[1:]:
	data.append(elem.split("\t"))

cols = lines[0].split("\t")
print(cols)
df = pd.DataFrame(data, columns =cols) 
input(df.head())
df.to_csv('abdominal_features.csv', index = None)

#for line in lines:
#	input(line)

#import docx
import io
import pandas as pd

#content = docx.Document(docname).paragraphs[0].text
# or if all paragraphs
# content = '\n'.join([p.text for p in docx.Document('data.docx').paragraphs

#df = pd.read_csv(docname)
#print(df.head())
