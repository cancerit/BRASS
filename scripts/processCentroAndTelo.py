import sys
import gzip
from collections import defaultdict

if __name__ == "__main__":

    d_telo= defaultdict(list)
    d_centro = {}
    with gzip.open(sys.argv[1], "rt") as ifi:
        for line in ifi:
            if not line.startswith("#"):
                tokens = line.strip().split()
                if "telo" in tokens[7]:
                    d_telo[tokens[1]].append(str(tokens[3]))
                else:
                    d_centro[tokens[1]] = (tokens[2], tokens[3])

    print("chr\tptel\tcen_start\tcen_end\tqtel\tcomment")
    for i in d_telo:
        d_telo[i] = sorted(d_telo[i])
        x = "\t".join([str(k) for k in [i, d_telo[i][0], d_centro[i][0], d_centro[i][1], d_telo[i][1], "."]])
        print(x)
