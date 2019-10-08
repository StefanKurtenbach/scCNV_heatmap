# scCNV_heatmap
version = "1.1.1"
# Stefan Kurtenbach
# Stefan.Kurtenbach@med.miami.edu
# make upper and lower ploidies accessable from terminal
# make keepSVG accessable from terminal

import os
import copy
import argparse
import csv
import cairosvg

def write_to_file(row):
    with open(filename, "a") as f:
        f.write(row)
        f.write("\n")
def draw_rect(x_coord, y_coord, color, width, hight1, opacity, stroke_thickness): #upper left corner is anchor
    return '''<rect x="''' + str(x_coord) + '''" opacity="''' + str(opacity) + '''" y="''' + str(y_coord) + '''" fill="''' + color + '''" width="''' + str(width) + '''" height="''' + str(hight1) + '''" style="stroke: #000000; stroke-width:''' + str(stroke_thickness) + ''';stroke:#000000"/>'''


colors = ["#4D80FF", "#6CA8FF", "#E2E2E2", "#FFDBA1", "#FFC373", "#FFAC33", "#FF9700", "#FF8500", "#FF5500", "#FF1E00"]

keep_output_SVG = True
desired_width = float(15000)
desired_hight = float(20000)

parser = argparse.ArgumentParser(description='CNV_args')
parser.add_argument('-o','--output_filename', help='output name', required=False, type=str)
parser.add_argument('-c','--cnv_call', help='node unmerged cnv call file', required=False, type=str)
parser.add_argument('-p','--per_cell_summary_metrics', help='per cell summary metrics CNV file', required=False, type=str)
parser.add_argument('-e','--exclude_noisy_cells', help='set to "no" if noisy cells should be included', required=False, type=str, default="yes")
args = vars(parser.parse_args())

print("")
print("scCNV heatmap version " + str(version) + " initiated")

filename = args['output_filename']
filename += ".svg"

exclude_cells = args['exclude_noisy_cells']
print(exclude_cells)
if exclude_cells == "yes":
    print("Noisy cells will be excluded (default)")
elif exclude_cells == "no":
    print("Noisy cells will NOT be excluded")

CNV_data = args['cnv_call']

cell_summary_metrics_file = args['per_cell_summary_metrics']
cells = []
with open(cell_summary_metrics_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    for x, row in enumerate(csv_reader):
        if x > 0:
            if exclude_cells == "yes":
                if str(row[17]) == "0":
                    cells.append([row[1], row[14]])
            else:
                if float(row[14]) > 1.1:
                    cells.append([row[1], row[14]])

cells.sort(key=lambda x: x[1])

cell_ID = []
ploidy = []
for x, i in enumerate(cells):
    cell_ID.append(int(i[0]))
    ploidy.append(float(i[1]))
cells_by_ploidy = [cell_ID, ploidy]

nr_of_cells_to_include = len(cell_ID) # can be altered to speed up analysis for tests
print("Number of cells detected: " + str(nr_of_cells_to_include))

chromosome_sizes = [["13", 115169878]]

#chromosome_sizes = [["1", 249250621],["2", 243199373],["3", 198022430],["4", 191154276],["5", 180915260],["6", 171115067],["7", 159138663], ["8", 146364022],["9", 141213431],["10", 135534747],["11", 135006516],["12", 133851895],["13", 115169878],["14", 107349540],["15", 102531392],["16", 90354753],["17", 81195210],["18", 78077248], ["19", 59128983],["20", 63025520], ["21", 48129895],["22", 51304566]] #, ["X", 155270560], ["Y", 59373566]]

chromosomal_length = 0
for i in chromosome_sizes:
    chromosomal_length += i[1]

ploidy_lower_threshold = 0 # make this 0 to include everything
chromosome_bar_thickness_percentage = 2
ploidy_bar_width_percentage = 4
space_ploidy_bar_percentage = 2
chromosome_bar_thickness = desired_hight*(chromosome_bar_thickness_percentage/100)
space_ploidy_bar = desired_width *(space_ploidy_bar_percentage/100)
ploidy_bar_width = desired_width * (ploidy_bar_width_percentage/100)
ploidy_bar_hight = desired_hight - chromosome_bar_thickness

real_hight = desired_hight - chromosome_bar_thickness
real_width = desired_width - ploidy_bar_width - space_ploidy_bar
length_factor = float(real_width)/chromosomal_length
stroke = 8


if os.path.exists(filename):
    os.remove(filename)

if os.path.exists(filename + ".png"):
    os.remove(filename + ".png")

write_to_file('''<svg viewBox="0 0 ''' + str(desired_width) + " " + str(desired_hight) + '''" xmlns="http://www.w3.org/2000/svg">''')
write_to_file(draw_rect(0, 0, "#FFFFFF", desired_width, desired_hight, 100, 0)) #background


# Draws chromosome sizes

color1 = "#E2E2E2"
color2 = "#5E5E5E"
start = 0
for x, i in enumerate(chromosome_sizes):  # alternate chromosome color
    if (x % 2) == 0:
        color = color1
    else:
        color = color2
    write_to_file(draw_rect(start, 0, color, i[1]*length_factor, chromosome_bar_thickness, 100, 0))
    start += i[1] * length_factor


scaling_factor_y_pos = (real_hight/(len(cells_by_ploidy[0]))) # scaling factor to determine y pos. e.g. cell 100 is not at 100, but at 120 if scalingf. is 1.2

with open(CNV_data) as f:
    for x, line in enumerate(f):
        if x > 2: # get rid of headers
            line_split = line.split("\t")
            cellID_in_data = int(line_split[3])
            if cellID_in_data in cells_by_ploidy[0]:
                ploidy_of_cellID = int(line_split[4])
                if ploidy_of_cellID > 9:
                    color = colors[-1]
                else:
                    color = colors[ploidy_of_cellID]
                start_pos = 0    
                for i in chromosome_sizes:  # this allows groups to contain random cells
                    if i[0] != line_split[0]:
                        start_pos += i[1]
                    else:
                        start_pos += int(line_split[1])
                        start_pos *= length_factor
                        break
                y_pos = cells_by_ploidy[0].index(cellID_in_data) * scaling_factor_y_pos
                write_to_file(draw_rect(start_pos, y_pos + chromosome_bar_thickness, color, (int(line_split[2])-int(line_split[1]))*length_factor, scaling_factor_y_pos, 100, 0))

### Ploidy bar
ploidy_position_percentages = []
ploidy_counter = 0.1 # define interval
percentage_change = 0
ploidies = copy.copy(cells_by_ploidy[1])

while percentage_change < 100:
    if percentage_change == 100:
        break
    nr_of_cells_with_ploidy = 0
    for i in ploidies:
        if i > ploidy_counter - 0.1:
            if i < ploidy_counter:
                nr_of_cells_with_ploidy += 1
    percentage_change += (100*nr_of_cells_with_ploidy)/nr_of_cells_to_include
    if percentage_change > 99:
        percentage_change = 100
    ploidy_position_percentages.append(percentage_change)
    ploidy_counter += 0.1

print(ploidy_position_percentages)


start = chromosome_bar_thickness
summed_height = 0

for x, i in enumerate(ploidy_position_percentages): # x is 0.1 ploidy intervals starting with 0.1
    height = (i*(ploidy_bar_hight/100)) - summed_height
    color = colors[round((x+1)/10)]
    write_to_file(draw_rect(str(real_width + space_ploidy_bar), str(start), color, str(ploidy_bar_width), str(height), 100, stroke))
    summed_height += height
    start += height
write_to_file("</svg>")

# Write output
if os.path.exists(filename + ".png"):
    os.remove(filename + ".png")

cairosvg.svg2png(url=filename, write_to=filename + ".png", parent_width=desired_width, parent_height=desired_hight)

if keep_output_SVG == False:
    if os.path.exists(filename):
        os.remove(filename)
