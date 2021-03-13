
list1 = []
ifn1 = r'test.target'
infile1 = open(ifn1,'r')
for eachline1 in infile1.readlines():
    line1 = eachline1
    print(line1)
    list1.append(line1.replace('\n',''))
# print(list1)


list2 =[]
ifn2 =r'products prediction form transoformer'
infile2 = open(ifn2,'r')
for eachline2 in infile2.readlines():
    line2 = eachline2
    list2.append(line2.replace('\n',''))
#print(list2)

list3 = []
ifn3 = r'decode_this.txt'
infile3 = open(ifn3,'r')
for eachline3 in infile3.readlines():
    line3 = eachline3
    print(line3)
    list3.append(line3.replace('\n',''))

same_number = 0
f1 = open('right_reactant.txt','w')
f2 = open('right_product.txt','w')
# n can be 1,2,3...
for i in range(len(list1)):
    flag = False
    for l in range(i*n,i*n+n):
        if list2[l] == list1[i]:
            line5 = list2[l].rstrip()
            f1.write(list3[i] + '\n')
            f2.write(list5[i] + '\n')
            flag = True
            same_number += 1
        if flag:
            break

    if flag == False:
        f1.write('\n')
        f2.write( '\n')

print(same_number/len(list1))





