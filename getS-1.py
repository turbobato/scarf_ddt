S = [0,2,4,12,8,14,24,21,16,19,28,5,17,20,11,23,1,6,7,26,25,18,10,27,3,13,9,29,22,30,15,31]

Sinv = [0]*(1<<5)

for i in range(1<<5):
    Sinv[S[i]] = i

print(Sinv)