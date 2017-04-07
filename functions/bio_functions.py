

def gcContent(seq):
    " return GC content as a float "
    c = 0
    for x in seq:
        if x in ["G","C"]:
            c+= 1
    return (float(c)/len(seq))