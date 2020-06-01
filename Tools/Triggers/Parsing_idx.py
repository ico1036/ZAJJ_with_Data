import sys
import subprocess


# Usage: python Parsing_idx.py my_Ele_trigger_list.txt trigger_idx_2016HnB.txt



input_filename = sys.argv[1]
target_file = sys.argv[2]
arr_idx = []


with open(input_filename,'r') as fhand:
	print("## Use following triggers: ")
	for i,line in enumerate(fhand):
		line = line.rstrip()
		print("#{0} {1}".format(i+1, line))
		arr_idx.append(line)



print(" ")
print("Start Parsing..\n")
for line in arr_idx:
	args = 'grep' + ' ' +  '"' + line + '"' + ' ' + target_file + ' ' + '| ' + ' ' + 'awk' + ' ' + "'{print $4, $3}'"
	subprocess.call(args,shell=True)
