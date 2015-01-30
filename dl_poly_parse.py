#!/usr/bin/env python3

# dl_poly_parse
# If ran as script, takes a DL_POLY OUTPUT file and returns the physical properties as a parsed
# file of simple columns, for easy readability by plotting software.
#
# To do:
#  * give option to output as csv
#  * give option to return properties as horizontally or vertically sorted
#  * allow importing as library to get single properties as lists
#  * functions to get rolling and total averages

OUTPUT = "OUTPUT"
PARSED = "parsed.txt"
BREAK = " ------------------------------------------------------------------------------------------------------------------------\n"

def get_lines():
    "Returns lines from file truncated to tabulated output data only."
    with open(OUTPUT, "r") as _file:
        lines = _file.readlines()
    return lines[lines.index(BREAK):]

def get_headers(lines):
    "Returns output data headers, solving the problem of one of them including a space."
    headers = lines[2].split() + lines[3].split() + lines[4].split()
    headers.remove("(s)")
    headers[headers.index("cpu")] = "cpu(s)"
    return headers


def get_property(lines, headers, prop, avg=False):
    """
    Returns all values of a named property as a list.
    If avg=True, returns list of rolling averages instead.
    """
    prop_index = headers.index(prop)
    prop_list = []

    if avg:         # THIS IS SEVERELY BROKEN!
        offset = 4
    else:
        offset = 0

    for i, l in enumerate(lines):
        if l == BREAK and len(lines[i+1]) == 118:
            values = lines[i+1+offset].split() + lines[i+2+offset].split() + lines[i+3+offset].split()
            try:
                prop_list.append(float(values[prop_index]))
            except ValueError:
                prop_list.append(values[prop_index])

    return prop_list[:-1] # discard final item as it is the total average

def get_final_avg(lines, headers, prop):
    "Returns the final average of a named property."

    prop_index = headers.index(prop)

    # total averages are given after the second-to-last BREAK
    avg_line = [i for i, line in enumerate(lines) if line == BREAK][-2]
    values = lines[avg_line+1].split() + lines[avg_line+2].split() + lines[avg_line+3].split()

    try:
        return float(values[prop_index])
    except ValueError:
        return values[prop_index]

def sort_by_column(unsorted):
    """
    Returns list reading down each column of 3 in OUTPUT rather than across each row
    this puts certain values usefully adjacent to each other e.g. time, step, cpu
    but separates others e.g. alpha, beta, gamma.
    """
    sort = []
    for i in range(0,len(unsorted)):
        triple = unsorted[i::10]
        for j in triple:
            sort.append(j)
    return sort[:30]

def get_all_props(lines, avg=False):
    """
    Returns all properties as a list of lists, with same ordering as headers.
    If avg=True, returns lists of rolling averages instead.
    """
    properties = []

    if avg:
        offset = 4
    else:
        offset = 0

    for i, l in enumerate(lines):
        if l == BREAK and len(lines[i+1]) == 118: # data always found in lines of 10 after BREAK

            values = lines[i+1+offset].split() + lines[i+2+offset].split() + lines[i+3+offset].split()
            if properties == []:    # fill with lists of initial values if empty
                properties = [[float(v)] for v in values]
                properties[0][0] = int(properties[0][0])
            else:                   # append otherwise
                for j, p in enumerate(properties):
                    if j == 0: # the first value, step, is the only int
                        p.append(int(values[j]))
                    # final averages give ******** for some values, which throw ValueError
                    # these are appended to the list as strings as they will later be discarded
                    else:
                        try:
                            p.append(float(values[j]))
                        except ValueError:
                            p.append(values[j])

    return [p[:-1] for p in properties] # discard final items as they are total averages
 
def main():
    lines = get_lines()
    headers = get_headers(lines)
    properties = get_all_props(lines)
    tot_num = len(properties[0])

    sorted_headers = sort_by_column(headers)
    sorted_props = sort_by_column(properties)

    parsed = ""
    for h in sorted_headers:
        parsed += "%-12s" % (h)
    for i in range(tot_num):
        parsed += "\n"
        for p in sorted_props:
            parsed += "%-11s " % (p[i])


    with open(PARSED, "w") as _file:
        _file.write(parsed)

if __name__ == '__main__':
    main()
