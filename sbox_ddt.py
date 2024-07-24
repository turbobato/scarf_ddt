S = [0,2,4,12,8,14,24,21,16,19,28,5,17,20,11,23,1,6,7,26,25,18,10,27,3,13,9,29,22,30,15,31]

valid_diffs = set()

for alpha in range(1<<5):
    for x in range(1<<5):
        beta = S[alpha^x]^S[x]
        if beta == alpha :
            valid_diffs.add(alpha)

for x in valid_diffs : 
    print(hex(x))
