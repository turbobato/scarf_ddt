with open("avg_ddt_result_round_5bits.txt") as ddt : 
    ddt_tab = []
    i = 0
    for line in ddt :
        ddt_tab.append(list(map(float, line.strip().split(" "))))
        i+=1
    diffs = []
    # print(len(ddt),len(ddt[0]))
    for i in range(len(ddt_tab)):
        for j in range(len(ddt_tab)):
            if ddt_tab[i][j] != 0.0 :
                diffs.append((i,j))
    print(list(map((lambda xy : (hex(xy[0]),hex(xy[1]))), diffs)))